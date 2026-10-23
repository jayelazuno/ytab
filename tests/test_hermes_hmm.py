from itertools import product

import numpy as np
from scipy.optimize import approx_fprime

from ytab.essentiality import hermes_hmm as hh


def brute_force(logB, A, pi):
    T, K = logB.shape; total = 0.; gamma = np.zeros((T, K)); xi = np.zeros((K, K))
    for path in product(range(K), repeat=T):
        p = pi[path[0]]*np.exp(logB[0, path[0]])
        for t in range(1, T):p *= A[path[t-1], path[t]]*np.exp(logB[t, path[t]])
        total += p
        for t in range(T):gamma[t, path[t]] += p
        for t in range(1, T):xi[path[t-1], path[t]] += p
    return gamma/total, xi/total, np.log(total)


def test_forward_backward_matches_enumeration():
    rng = np.random.default_rng(1)
    for T in (1, 2, 5, 7):
        logB = rng.normal(-3, 2, (T, 4)); logB[min(2, T-1)] = 0.   # one missing window
        A = rng.dirichlet(np.ones(4), 4); pi = rng.dirichlet(np.ones(4))
        g, x, ll = hh.forward_backward(logB, A, pi); g0, x0, ll0 = brute_force(logB, A, pi)
        assert np.allclose(g, g0, atol=1e-10) and np.isclose(ll, ll0, atol=1e-9)
        assert np.allclose(x, x0, atol=1e-9)


def test_forward_backward_long_sequence_is_stable():
    rng = np.random.default_rng(2); T = 5000
    logB = rng.normal(-50, 20, (T, 4)); A = np.full((4, 4), .01); np.fill_diagonal(A, .97); pi = np.full(4, .25)
    g, x, ll = hh.forward_backward(logB, A, pi)
    assert np.all(np.isfinite(g)) and np.isfinite(ll) and np.allclose(g.sum(1), 1) and np.isclose(x.sum(), T-1)


def test_objective_gradient():
    rng = np.random.default_rng(3); T, L, P = 300, 2, 2
    X = rng.normal(size=(T, P)); Y = rng.poisson(2., (T, L)).astype(float); G = rng.dirichlet(np.ones(4), T)
    theta = hh.pack(np.log([.05, .4, 1., 2.5]), np.array([.3, .5]), np.array([.2, -.1]), np.log(.4))
    hist = hh.count_histograms(Y); f = lambda t:hh._objective(t, Y, X, G, hist)[0]
    assert np.allclose(hh._objective(theta, Y, X, G, hist)[1], approx_fprime(theta, f, 1e-6), rtol=1e-4, atol=1e-3)


def test_unique_kmers_detect_repeats_on_both_strands():
    rng = np.random.default_rng(4); unit = "".join(rng.choice(list("ACGT"), 40))
    comp = unit[::-1].translate(str.maketrans("ACGT", "TGCA"))
    seq = "".join(rng.choice(list("ACGT"), 200))+unit+"".join(rng.choice(list("ACGT"), 200))+comp+"".join(rng.choice(list("ACGT"), 200))
    u = hh.unique_kmer_starts({"c":hh.encode(seq)}, 20)["c"]
    assert not u[200:221].any() and not u[440:461].any()        # k-mers inside the repeat and its reverse complement
    assert u[:180].all() and not u[len(seq)-19:].any()          # unique flank; no full k-mer at the end


def simulate(rng, T=4000, L=2):
    A = np.array([[.97, .01, .02, 0], [.01, .97, .02, 0], [.01, .01, .97, .01], [0, 0, .03, .97]]); A /= A.sum(1, keepdims=True)
    states = np.zeros(T, dtype=int); states[0] = 2
    for t in range(1, T):states[t] = rng.choice(4, p=A[states[t-1]])
    X = rng.normal(size=(T, 1)); r = np.array([.03, .3, 1., 3.]); psi = .5
    mu = np.exp(np.log(2.)+.4*X[:, 0])[:, None]*r[states][:, None]*np.ones((1, L))
    Y = rng.poisson(mu*rng.gamma(1/psi, psi, (T, 1))).astype(float)
    return states, X, Y, r


def test_baum_welch_recovers_states_and_parameters():
    rng = np.random.default_rng(5); states, X, Y, r = simulate(rng)
    out = hh.fit([dict(Y=Y, X=X, keep=np.ones(len(Y), dtype=bool))], hh.Settings(max_iterations=200))
    assert np.all(np.diff(out["trace"]) > -1e-6*abs(out["trace"][0]))            # EM is monotone
    assert np.mean(out["gamma"].argmax(1) == states) > .85
    assert abs(np.log(out["r"][1]/r[1])) < .3 and abs(out["b"][0]-.4) < .1


def test_gene_windows_trim_three_prime_end():
    s = hh.Settings(window=10, cterm=10.)
    plus = hh.gene_windows(dict(blocks=[(1, 100)], strand="+"), 10, s); minus = hh.gene_windows(dict(blocks=[(1, 100)], strand="-"), 10, s)
    assert [w for w, _ in plus] == list(range(9)) and [w for w, _ in minus] == list(range(1, 10))


def test_negative_multinomial_matches_nb_and_integration():
    from scipy.integrate import quad
    from scipy.stats import gamma, poisson
    y, mu, psi = np.array([[3.]]), np.array([[2.5]]), .7
    assert np.isclose(hh.nm_logpmf(y, mu, psi)[0], hh.nb_logpmf(3., 2.5, .7))
    y2, mu2 = np.array([[2., 5.]]), np.array([[1.5, 4.]])
    f = lambda u:poisson.pmf(2, 1.5*u)*poisson.pmf(5, 4*u)*gamma.pdf(u, 1/psi, scale=psi)
    assert np.isclose(np.exp(hh.nm_logpmf(y2, mu2, psi)[0]), quad(f, 0, np.inf)[0], rtol=1e-7)


def test_run_probability_matches_enumeration():
    rng = np.random.default_rng(6); T, lo, hi, m = 8, 2, 7, 2
    logB = rng.normal(-2, 1.5, (T, 4)); A = rng.dirichlet(np.ones(4)*2, 4); pi = rng.dirichlet(np.ones(4))
    kept = np.array([1, 1, 1, 0, 1, 1, 1, 1], dtype=bool)
    gamma, _, _, B, beta = hh.forward_backward(logB, A, pi, full=True)
    got = hh.run_probability(gamma[lo], A, B[lo:hi], beta[lo:hi], kept[lo:hi], m)
    total = hit = 0.
    for path in product(range(4), repeat=T):
        p = pi[path[0]]*np.exp(logB[0, path[0]])
        for t in range(1, T):p *= A[path[t-1], path[t]]*np.exp(logB[t, path[t]])
        total += p; run = best = 0
        for t in range(lo, hi):
            if path[t] == 0:run += kept[t]
            else:run = 0
            best = max(best, run)
        hit += p*(best >= m)
    assert np.isclose(got, hit/total, atol=1e-10)


def test_learned_target_model_recovers_planted_motif_on_both_strands():
    rng = np.random.default_rng(3); L = 200_000
    codes = rng.integers(0, 4, L).astype(np.uint8)
    # Insertions only at oriented T(+1)...A(+6), i.e. nTnnnnAn read outward; half on each strand.
    pos, strand = [], []
    for i in rng.choice(np.arange(50, L-50), 4000, replace=False):
        if rng.random() < .5:
            codes[i+1], codes[i+6] = 3, 0; pos.append(i+1); strand.append("W")
        else:   # C site at p: oriented base j is the complement of genome p-j
            codes[i-1], codes[i-6] = 0, 3; pos.append(i+1); strand.append("C")
    model = hh.learn_target_model({"c": codes}, {"c": (np.array(pos), np.array(strand))})
    top = model["offsets"][np.argsort(model["information"])[-2:]]
    assert sorted(top.tolist()) == [1, 6]
    k1, k6 = list(model["offsets"]).index(1), list(model["offsets"]).index(6)
    assert "ACGT"[int(np.argmax(model["logodds"][k1, :4]))] == "T" and "ACGT"[int(np.argmax(model["logodds"][k6, :4]))] == "A"
    prop = hh.target_propensity(codes, model)
    sites = np.array(pos)-1
    assert prop[sites].mean() > 3*np.median(prop)


def test_target_propensity_is_strand_symmetric():
    rng = np.random.default_rng(5); codes = rng.integers(0, 4, 5000).astype(np.uint8)
    model = dict(offsets=np.arange(-3, 5), logodds=np.column_stack([rng.normal(0, 1, (8, 4)), np.zeros(8)]))
    rc = hh.COMPLEMENT[codes[::-1]]
    # the propensity of a position equals that of its mirror position on the reverse-complement genome
    assert np.allclose(hh.target_propensity(codes, model)[10:-10], hh.target_propensity(rc, model)[::-1][10:-10])
