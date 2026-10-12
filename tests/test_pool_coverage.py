"""Numerical and audited-data validation of inferred pool-group sampling."""
import itertools
import tempfile
import unittest
from pathlib import Path
import numpy as np
from ytab.screen_coverage.pool_coverage import (pool_identity,pool_distributions,load_pool_catalogue,
    expected_recovery,integer_allocation,optimize_allocation,monte_carlo_pools,monte_carlo_pool_alleles,
    equal_allocation_target,recovery_denominators,complexity_allocation,
    balanced_saturation_allocation,pool_gene_saturation,pool_distinct_allele_expectations)
from ytab.screen_coverage.statistics import poisson_coverage,colonies_for_target_coverage


class PoolTests(unittest.TestCase):
    def test_explicit_recovery_denominators(self):
        row=recovery_denominators(dict(represented=4379,lower_genes=4000,upper_genes=4379),4379,5220)
        self.assertEqual(row["represented_union_percent"],100)
        self.assertAlmostEqual(row["represented_reference_percent"],83.8888888889)
        self.assertEqual(row["represented_denominator"],4379)
        self.assertEqual(row["reference_denominator"],5220)
        self.assertNotEqual(row["lower_genes_union_percent"],row["lower_genes_reference_percent"])

    def test_equal_target_not_immediate_gene_optimizer(self):
        p=np.array([[.2,.001],[.15,.005],[0,.02],[0,.015]])
        target=equal_allocation_target(p,.9)
        allocation=target["allocation"];k=target["k"]
        self.assertLessEqual(np.ptp(allocation),1)
        self.assertGreaterEqual(expected_recovery(p,allocation).mean(),.9)
        self.assertLess(expected_recovery(p,integer_allocation(k-1,[1,1])).mean(),.9)
        self.assertIsNone(equal_allocation_target(np.array([[.1,0],[0,0]]),.9))

    def test_default_runner_complexity_and_manual(self):
        # Small synthetic catalogue exercises the orchestration, not just helpers.
        from unittest.mock import patch
        from ytab.screen_coverage.pool_runner import calculate_pool_screen
        p=np.array([[.8,.05],[0,.05],[.2,.9]])
        data=dict(pool_ids=["A","B"],summary=dict(represented=2,annotated=3),
            pools=[dict(pool_id="A",sites=2,coding_sites=1,noncoding_sites=1,coding_fraction=.5,label="A",unique_genes=1,unique_coding_alleles=1,ge2=0,ge3=0),dict(pool_id="B",sites=3,coding_sites=2,noncoding_sites=1,coding_fraction=2/3,label="B",unique_genes=1,unique_coding_alleles=1,ge2=0,ge3=0)],
            overlaps=[],observations=[dict(source="a"),dict(source="b")],sources=[],
            gene_probabilities={"equal":p[:2],"abundance":p[:2]},
            allele_genes=np.array([0,1,-1]),probabilities={"equal":p,"abundance":p},
            gene_ids=["g1","g2"],gene_pool_presence=p[:2]>0,observation_reads=p,
            gene_lengths=np.array([1000,2000]),gene_coordinate_counts=np.array([[1,0],[0,1]]))
        with patch("ytab.screen_coverage.pool_runner.load_pool_catalogue",return_value=data):
            result=calculate_pool_screen(dict(k=11,iterations=10,simulate=False))
            self.assertEqual(result["settings"]["primary_strategy"],"complexity")
            rows=[r for r in result["practical"] if r.get("model") and r["k"]==11 and r["strategy"]=="complexity"]
            self.assertEqual(len(rows),2)
            self.assertTrue(all(r["allocation"]==[4,7] for r in rows))
            preview=calculate_pool_screen(dict(k=11,iterations=10,simulate=False,planner_only=True,include_discovery=False))
            self.assertTrue(all(r["strategy"] in {"complexity","balanced","equal"} for r in preview["baselines"]))
            self.assertEqual(preview["discovery"],[])
            equal_preview=calculate_pool_screen(dict(k=11,iterations=10,simulate=False,planner_only=True,
                include_discovery=False,model_mode="equal"))
            self.assertEqual({r["model"] for r in equal_preview["baselines"]},{"equal"})
            self.assertTrue(all(r["model"]=="equal" for r in equal_preview["pool_saturation"]))
            equal_targets=[r for r in equal_preview["targets"] if r["method"]=="Hermes Empirical Monte Carlo — Equal Insertion"]
            self.assertTrue(equal_targets)
            self.assertTrue(all(r["k"] is not None for r in equal_targets))
            result=calculate_pool_screen(dict(k=999,iterations=10,simulate=False,allocation_mode="manual",manual_allocation={"A":3,"B":8}))
            self.assertEqual(result["settings"]["k"],11)
            rows=[r for r in result["practical"] if r.get("model") and r["k"]==11 and r["strategy"]=="manual"]
            self.assertEqual([(r["strategy"],r["allocation"]) for r in rows],[("manual",[3,8])]*2)

    def test_pool_specific_support(self):
        # Drawing only B cannot recover A's coordinate, even if globally present.
        result=monte_carlo_pool_alleles([0,1,-1],[[.5,0],[0,.5],[.5,.5]],[[0,100]],2,100,8)
        self.assertEqual(result[0][0]["represented"],1)
        self.assertEqual(result[0][1]["represented"],0)

    def test_catalogue_builder_preserves_existing_audit(self):
        from ytab.screen_coverage.build_parent_catalogue import audit
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"parent_union_audit.json";path.write_text("preserve")
            with self.assertRaises(ValueError):audit(tmp,tmp,["future_project"])
            self.assertEqual(path.read_text(),"preserve")

    def test_numbered_independent_pool_identity_ignores_replicate_label(self):
        row=dict(species="glabrata",background="replicate-A",pool="1")
        self.assertNotEqual(pool_identity(row,"library-A"),pool_identity(row,"library-B"))
        self.assertEqual(pool_identity(row,"library-A"),pool_identity(dict(row,background="replicate-B"),"library-A"))
        self.assertEqual(pool_identity(row,"library-A"),pool_identity(dict(row,sample="another sequencing observation"),"library-A"))
        self.assertEqual(pool_identity(dict(species="glabrata",pool="1"),"lib"),pool_identity(dict(species="glabrata",pool="1",background="other"),"lib"))

    def test_complexity_normalized_allocation_arbitrary_j(self):
        for counts in ([8],[8,12],[8,12,20,30,40,50],[83391,112214,72500,70958]):
            allocation=complexity_allocation(25001,counts)
            self.assertEqual(int(allocation.sum()),25001)
            np.testing.assert_allclose(allocation/25001,np.asarray(counts)/sum(counts),atol=1/25001+1e-10)
        counts=np.array([83391,112214,72500,70958]);allocation=complexity_allocation(20000,counts)
        self.assertEqual(int(np.argmax(allocation)),1)
        self.assertLess(np.ptp(allocation/counts),1/70000)

    def test_balanced_saturation_keeps_all_pools_and_balances(self):
        p=np.array([[.6,.1,.05],[.1,.4,.1],[0,.1,.35],[0,0,.2]])
        present=p>0
        allocation=balanced_saturation_allocation(300,p,present)
        self.assertEqual(int(allocation.sum()),300)
        self.assertTrue(np.all(allocation>0))
        saturation=pool_gene_saturation(p,allocation,present)
        self.assertLess(np.ptp(saturation),.03)

    def test_balanced_allocator_supports_arbitrary_pool_count(self):
        rng=np.random.default_rng(73)
        for pools in (1,2,3,7):
            p=rng.uniform(.001,.025,size=(30,pools))
            # Leave probability mass for non-target/intergenic insertion picks.
            p*=np.minimum(1.,.8/p.sum(axis=0))[None,:]
            present=p>0
            allocation=balanced_saturation_allocation(max(200,pools),p,present)
            self.assertEqual(int(allocation.sum()),max(200,pools))
            self.assertTrue(np.all(allocation>=1))

    def test_pool_saturation_and_distinct_allele_marginals(self):
        p=np.array([[.2],[.3],[.5]])
        presence=np.array([[1],[1]],dtype=bool)
        self.assertAlmostEqual(pool_gene_saturation(p[:2],1,presence)[0],(.2+.3)/2)
        expected=pool_distinct_allele_expectations([0,0,-1],[.2,.3,.5],[1,2,10],1)
        self.assertTrue(np.all(expected[:,0]>=expected[:,1]))
        self.assertTrue(np.all(expected[:,1]>=expected[:,2]))
        self.assertTrue(np.all(np.diff(expected[:,0])>=0))

    def test_baym_style_site_to_gene_occupancy_formula(self):
        # Three equally likely insertion-bearing mutant coordinates map to
        # genes A, A and B. Repeated picks count a gene only on its first hit.
        gene_probabilities=np.array([[2/3],[1/3]])
        recovered=expected_recovery(gene_probabilities,[3])
        self.assertAlmostEqual(recovered[0],1-(1/3)**3)
        self.assertAlmostEqual(recovered[1],1-(2/3)**3)
        self.assertAlmostEqual(recovered.sum(),monte_carlo_pools(
            gene_probabilities,[[3]],iterations=50000,seed=2026)[0]["represented"],delta=.015)

    def test_within_pool_probabilities_use_full_catalogue_and_normalize_observations(self):
        reads=np.array([[10,100,0],[0,100,8],[10,0,2]],dtype=float)
        presence,p=pool_distributions(reads,["A","A","B"])
        np.testing.assert_array_equal(presence,[[1,0],[1,1],[1,1]])
        np.testing.assert_allclose(p["equal"][:,0],[1/3,1/3,1/3])
        np.testing.assert_allclose(p["abundance"][:,0],[.5,.25,.25])
        np.testing.assert_allclose(p["equal"][:,1],[0,.5,.5])
        np.testing.assert_allclose(p["abundance"][:,1],[0,.8,.2])
        for values in p.values():np.testing.assert_allclose(values.sum(axis=0),1)
        _,rep=pool_distributions(np.c_[reads,reads[:,2]*1000],["A","A","B","B"])
        for model in p:np.testing.assert_allclose(p[model],rep[model])

    def test_noncoding_insertion_positive_picks_consume_budget(self):
        # One coding coordinate and one noncoding coordinate are equally likely.
        # K is the total insertion-positive pick count, not coding-conditioned K.
        p=np.array([[.5],[.5]])
        result=monte_carlo_pool_alleles([0,-1],p,[[1],[20]],1,20000,42)
        self.assertAlmostEqual(result[0][0]["mean"],.5,delta=.015)
        self.assertAlmostEqual(result[1][0]["mean"],1.,delta=1e-12)
        coding_conditioned=monte_carlo_pool_alleles([0],[[1.]],[[1]],1,20000,42)
        self.assertGreaterEqual(coding_conditioned[0][0]["mean"],result[0][0]["mean"])
        self.assertLessEqual(result[0][0]["mean"],coding_conditioned[0][0]["mean"])

    def test_arbitrary_pool_count_and_constraints(self):
        for n in (1,2,3,4,7):
            p=np.full((10,n),.025)
            fit=optimize_allocation(p,101,minimum=np.ones(n,dtype=int))
            self.assertEqual(sum(fit["allocation"]),101)
            self.assertTrue(np.all(fit["allocation"]>=1))
        a=integer_allocation(20,[1,2,3],minimum=[2,3,4],maximum=[3,8,20])
        self.assertEqual(sum(a),20);self.assertLessEqual(a[0],3)
        self.assertLessEqual(np.ptp(integer_allocation(101,np.ones(7))),1)
        with self.assertRaises(ValueError):integer_allocation(2,[1,1],minimum=[2,2])

    def test_optimizer_expected_objective_and_changing_budget(self):
        p=np.array([[.2,.001],[.15,.005],[0,.02],[0,.015]])
        ratios=[]
        for k in (5,20,100,1000):
            fit=optimize_allocation(p,k,[100,30],[20,12]);ratios.append(fit["allocation"][0]/k)
            for baseline in fit["baselines"].values():
                self.assertGreaterEqual(fit["expected"]+1e-9,expected_recovery(p,baseline).sum())
            if k<=20:
                brute=max(expected_recovery(p,[x,k-x]).sum() for x in range(k+1))
                self.assertAlmostEqual(fit["expected"],brute,places=7)
                self.assertGreaterEqual(fit["expected"]+fit["relaxation_gap"]+1e-8,brute)
        self.assertGreater(np.ptp(ratios),.1)

    def test_manual_and_equal_simulation_reproducible(self):
        p=np.array([[.2,.01],[.1,.1],[0,.2]])
        allocations=[[0,0],[10,10],[3,17],[10,20],[20,40]]
        a=monte_carlo_pools(p,allocations,3000,6)
        self.assertEqual(a,monte_carlo_pools(p,allocations,3000,6))
        for allocation,row in zip(allocations,a):
            expected=expected_recovery(p,allocation).sum()
            self.assertLess(abs(row["represented"]-expected),max(.02,6*row["sd_genes"]/np.sqrt(3000)))
        self.assertLessEqual(a[1]["mean"],a[3]["mean"])
        self.assertLessEqual(a[3]["mean"],a[4]["mean"])

    def test_coordinate_identity_shared_across_pools(self):
        # Both physical pools contain the SAME coordinate: two colonies != two alleles.
        result=monte_carlo_pool_alleles([0],[[1,1]],[[1,1],[4,7]],1,100,2)
        for rows in result:self.assertEqual([r["represented"] for r in rows],[1,0,0])
        # Different coordinates in the same gene count as independent alleles.
        result=monte_carlo_pool_alleles([0,0],[[1,0],[0,1]],[[1,1]],1,100,2)
        self.assertEqual([r["represented"] for r in result[0]],[1,1,0])
        # Noncoding coordinates are sampled too, but have no gene assignment.
        result=monte_carlo_pool_alleles([0,-1],[[.5],[.5]],[[1]],1,5000,2)
        self.assertAlmostEqual(result[0][0]["mean"],.5,delta=.04)

    def test_allele_order_ceiling_reproducibility(self):
        gids=[0,0,0,1,2,-1];p=np.array([[.1,.1],[.1,0],[0,.1],[.1,.1],[0,0],[.7,.7]])
        a=monte_carlo_pool_alleles(gids,p,[[2,3],[20,30],[100,100]],3,1000,7)
        self.assertEqual(a,monte_carlo_pool_alleles(gids,p,[[2,3],[20,30],[100,100]],3,1000,7))
        for rows in a:
            self.assertGreaterEqual(rows[0]["mean"],rows[1]["mean"])
            self.assertGreaterEqual(rows[1]["mean"],rows[2]["mean"])
            self.assertLessEqual(rows[0]["mean"],2/3)
        self.assertTrue(all(a[i][0]["mean"]<=a[i+1][0]["mean"] for i in (0,1)))

    def test_poisson_accepted_N(self):
        self.assertEqual([colonies_for_target_coverage(4379,t) for t in (.9,.95,.99,.999)],[10084,13119,20167,30250])
        self.assertAlmostEqual(poisson_coverage(3000,1000),1-np.exp(-3))

    def test_accepted_audit_and_pool_inventory(self):
        repo=Path(__file__).resolve().parents[1]
        audit=repo/"docs/codex/parent_union_audit"
        if not audit.exists():self.skipTest("Accepted project audit not installed")
        data=load_pool_catalogue(audit,repo)
        s=data["summary"]
        self.assertEqual([s[k] for k in ("annotated","sites","coding_sites","represented","ge2","ge3","missing")],
                         [5220,169099,49034,4379,3789,3343,841])
        self.assertEqual(len(data["pools"]),4)
        self.assertEqual([p["unique_genes"] for p in data["pools"]],[36,173,97,104])
        for pool in data["pools"]:
            self.assertEqual(pool["sites"],pool["coding_sites"]+pool["noncoding_sites"])
            self.assertAlmostEqual(pool["coding_fraction"],pool["coding_sites"]/pool["sites"])
        self.assertTrue(all(p["observation_count"]==2 for p in data["pools"]))
        self.assertEqual(len(set(data["coordinates"])),169099)
        self.assertTrue(all("treated" not in o["sample"] for o in data["observations"]))
        for p in data["probabilities"].values():np.testing.assert_allclose(p.sum(axis=0),1)
        noncoding=data["allele_genes"]<0
        self.assertTrue(np.all(data["probabilities"]["equal"][noncoding[:,None]&data["presence"]]>0))
        self.assertTrue(np.all(data["probabilities"]["abundance"][noncoding[:,None]&data["presence"]]>0))
        self.assertTrue(np.all(data["gene_probabilities"]["equal"].sum(axis=0)<1))
        site_counts=[p["sites"] for p in data["pools"]]
        alloc=complexity_allocation(25000,site_counts)
        np.testing.assert_allclose(alloc/25000,np.asarray(site_counts)/sum(site_counts),atol=1/25000+1e-10)
        self.assertEqual(int(np.argmax(alloc)),1)
        self.assertGreater(np.sum(data["presence"].sum(axis=1)>1),0)
        # Classification-independent representation is constructed exclusively from audit CDS mapping.
        self.assertEqual(sum(np.any(data["gene_probabilities"]["equal"]>0,axis=1)),4379)
        one=load_pool_catalogue(audit,repo,selected_pools=[data["pool_ids"][0]])
        self.assertEqual(one["summary"]["represented"],4379)
        self.assertEqual(one["summary"]["accessible_genes"],3677)
        self.assertAlmostEqual(one["summary"]["selected_pool_ceiling"],3677/4379)


if __name__=="__main__":unittest.main()
