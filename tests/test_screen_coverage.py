"""Run: python -m unittest discover -s tests -p test_screen_coverage.py"""
import math
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np

from ytab.screen_coverage.statistics import (
    poisson_probability, poisson_coverage, colonies_for_target_coverage,
    colonies_for_redundancy_target, redundancy_probability, empirical_gene_probabilities,
    weighted_poisson_coverage, exact_empirical_coverage, empirical_target,
    coverage_ceiling, monte_carlo_uniform, monte_carlo_hermes, empirical_rarefaction, monte_carlo_alleles)
from ytab.screen_coverage.library import read_coding_genes, load_library, is_parent_sample
from ytab.screen_coverage.audited_library import load_audited_library


class CoverageTests(unittest.TestCase):
    def test_distinct_alleles_not_repeated_colonies(self):
        result=monte_carlo_alleles([0],[1.],1,[1,2,20],iterations=100,seed=1)
        for row in result:self.assertEqual(row["mean"],1. if row["r"]==1 else 0.)
        result=monte_carlo_alleles([0,0],[.5,.5],1,[1,2,10],iterations=4000,seed=2)
        self.assertEqual(result,monte_carlo_alleles([0,0],[.5,.5],1,[1,2,10],iterations=4000,seed=2))
        self.assertLess(abs(next(r["mean"] for r in result if r["r"]==2 and r["k"]==2)-.5),.025)
        for row in result:
            if row["r"]==3:self.assertEqual(row["mean"],0.)
        for rr in [1,2,3]:
            self.assertTrue(np.all(np.diff([r["mean"] for r in result if r["r"]==rr])>=0))
        empty=monte_carlo_alleles([0],[0],2,[100],iterations=10)
        self.assertTrue(all(r["mean"]==0 for r in empty))

    def test_accepted_audit_validation_gates(self):
        audit=Path(__file__).resolve().parents[1]/"docs/codex/parent_union_audit"
        if not audit.exists():self.skipTest("Accepted project audit not available in this checkout")
        d=load_audited_library(audit)
        self.assertTrue(d["audit_matches"])
        self.assertEqual(d["universes"],dict(all=5220,observed=4387,informative=4379,alleles2=3789,alleles3=3343))
        self.assertEqual(d["unique_sites"],169099)
        self.assertEqual(d["informative_sites"],49034)
        self.assertEqual(d["representation"]["missing"],841)
        self.assertEqual(d["representation"]["shared_sites"],50704)
        self.assertEqual(d["classifier"],"Not used")
        self.assertEqual(len(d["_allele_genes"]),49034)
        self.assertTrue(all(is_parent_sample(r) for r in d["pools"]))
        genes=d["genes"]
        self.assertAlmostEqual(sum(g["p_sites"] for g in genes),49034/169099)
        self.assertEqual(sum(g["p_h2o2"]>0 for g in genes),4184)
        self.assertEqual(sum(g["p_zn"]>0 for g in genes),3661)
        for b in ("sites","h2o2","zn"):
            np.testing.assert_allclose(np.bincount(d["_allele_genes"],weights=d["_allele_p"][b],minlength=len(genes)),[g[f"p_{b}"] for g in genes])
        self.assertFalse(np.allclose([g["p_h2o2"] for g in genes],[g["p_zn"] for g in genes]))
        conditional=load_audited_library(audit,conditioning="conditional")
        for b in ("sites","h2o2","zn"):
            self.assertAlmostEqual(sum(g[f"p_{b}"] for g in conditional["genes"]),1.)
        self.assertEqual(d["universes"],conditional["universes"])
        with self.assertRaises(ValueError):load_audited_library(audit,selected_sources=["Zn_toxicity_screen::parent-treated"])

    def test_parent_roles(self):
        self.assertFalse(is_parent_sample(dict(sample="yH298-parent-pool1-treated",control_or_treated="treated")))
        self.assertFalse(is_parent_sample(dict(sample="unknown")))
        self.assertTrue(is_parent_sample(dict(sample="yH298-parent-pool1",condition="parent")))
        self.assertTrue(is_parent_sample(dict(sample="yH298-parent-pool1-mock",control_or_treated="control")))

    def test_known_poisson_values(self):
        for depth, expected in [(1,.63212056),(2,.86466472),(3,.95021293),(4,.98168436),(4.60517,.99),(6.90776,.999)]:
            self.assertAlmostEqual(poisson_coverage(1000*depth,1000),expected,places=6)
        self.assertAlmostEqual(poisson_probability(0,1),math.exp(-1))
        for target in [.9,.95,.99,.999]:
            k=colonies_for_target_coverage(1000,target)
            self.assertGreaterEqual(poisson_coverage(k,1000),target)
            self.assertLess(poisson_coverage(k-1,1000),target)

    def test_redundancy_inverse(self):
        for r in [1,2,3,5,30,100]:
            k=colonies_for_redundancy_target(1000,.95,r)
            self.assertGreaterEqual(poisson_coverage(k,1000,r),.95)
            self.assertLess(poisson_coverage(k-1,1000,r),.95)

    def test_weighted_exact_invariants(self):
        uniform=np.full(1000,.001)
        for k in [0,1000,3000,10000]:
            self.assertAlmostEqual(weighted_poisson_coverage(k,uniform),poisson_coverage(k,1000))
            self.assertLess(abs(exact_empirical_coverage(k,uniform)-poisson_coverage(k,1000)),.0003)
        p=np.array([.3,.02,.001,0])
        for fn in [weighted_poisson_coverage,exact_empirical_coverage]:
            values=[fn(k,p) for k in [0,1,10,100,10000]]
            self.assertTrue(np.all(np.diff(values)>=0))
            for k in [0,1,10,100,10000]:
                self.assertLessEqual(fn(k,p,3),fn(k,p,2))
                self.assertLessEqual(fn(k,p,2),fn(k,p,1))
                self.assertLessEqual(fn(k,p,1),coverage_ceiling(p))
            self.assertEqual(fn(10**8,[0]),0)
        self.assertIsNone(empirical_target(p,.9)["colonies"])
        self.assertIsNone(empirical_target(p,.75)["colonies"])
        for exact in [True,False]:
            for r in [1,2,3,8]:
                k=empirical_target(p,.7,r,exact)["colonies"]
                fn=exact_empirical_coverage if exact else weighted_poisson_coverage
                self.assertGreaterEqual(fn(k,p,r),.7)
                self.assertLess(fn(k-1,p,r),.7)

    def test_probability_mass_retains_non_target_draws(self):
        np.testing.assert_allclose(empirical_gene_probabilities([2,3,0],10),[.2,.3,0])
        np.testing.assert_allclose(empirical_gene_probabilities([2,3,0]),[.4,.6,0])
        for p in [[.8,.8],[-.1,.3],[np.nan],[]]:
            with self.assertRaises(ValueError):coverage_ceiling(p)
        with self.assertRaises(ValueError):exact_empirical_coverage(1.5,[.2])

    def test_multinomial_convergence_and_reproducibility(self):
        a=monte_carlo_uniform(1000,[1000,3000],iterations=1000,seed=4)
        self.assertEqual(a,monte_carlo_uniform(1000,[1000,3000],iterations=1000,seed=4))
        for row in a:
            expected=exact_empirical_coverage(row["k"],np.full(1000,.001),row["r"])
            self.assertLess(abs(row["mean"]-expected),max(.001,5*row["sd"]/math.sqrt(1000)))
        b=monte_carlo_hermes([.1,.001,0],[0,10,1000],iterations=1000)
        for row in b:self.assertLessEqual(row["upper"],2/3)
        self.assertLess(abs(next(x["mean"] for x in b if x["k"]==10 and x["r"]==1)-exact_empirical_coverage(10,[.1,.001,0])),.02)

    def test_rarefaction_finite_endpoint(self):
        rows=empirical_rarefaction([1,2,4,0],[0,3,10],iterations=1000,other_sites=3)
        for row in rows:
            self.assertLess(abs(row["mean"]-row["exact_mean"]),.025)
            if row["k"]==10:
                self.assertEqual(row["sd"],0)
                self.assertEqual(row["mean"],{1:.75,2:.5,3:.25}[row["r"]])
        with self.assertRaises(ValueError):empirical_rarefaction([1,2],[4])

    def test_cds_union_not_gene_span(self):
        with tempfile.TemporaryDirectory() as directory:
            f=Path(directory)/"genes.gff"
            f.write_text("chr1\tx\tgene\t1\t100\t.\t-\t.\tlocus_tag=g\n"
                "chr1\tx\tCDS\t1\t10\t.\t-\t0\tlocus_tag=g;gene=ABC\n"
                "chr1\tx\tCDS\t90\t100\t.\t-\t0\tlocus_tag=g\n"
                "chr1\tx\tCDS\t90\t100\t.\t-\t0\tlocus_tag=g\n")
            g=read_coding_genes(f)["g"]
            self.assertEqual(g["length"],21)
            self.assertEqual(g["blocks"],[(1,10),(90,100)])

    def test_library_mapping_and_physical_mixture(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/"src/ytab").mkdir(parents=True)
            project=root/"output/projects/example"
            (project/"config").mkdir(parents=True)
            (project/"classifier/T1").mkdir(parents=True)
            (root/"genes.gff").write_text(
                "chr1\tx\tCDS\t1\t5\t.\t+\t0\tlocus_tag=a\n"
                "chr1\tx\tCDS\t10\t15\t.\t+\t0\tlocus_tag=a\n"
                "chr1\tx\tCDS\t30\t40\t.\t-\t0\tlocus_tag=b\n"
                "chr1\tx\tCDS\t40\t50\t.\t+\t0\tlocus_tag=c\n")
            cfg=project/"config/project.yaml"
            cfg.write_text(json.dumps(dict(output_project_dir="output/projects/example",sample_sheet="samples.csv",reference=dict(gff="genes.gff"))))
            (root/"samples.csv").write_text("sample,include,hit_file,condition\np1,true,p1.txt,parent\np2,true,p2.txt,mock\np3,true,p3.txt,parent\nparent-treated,true,p1.txt,treated\n")
            (root/"p1.txt").write_text("Chromosome\tHit position\tHit count\nchr1\t2\t2\nchr1\t2\t3\nchr1\t7\t5\nchr1\t40\t10\n")
            (root/"p2.txt").write_text("Chromosome\tHit position\tHit count\nchr1\t35\t80\nchr1\t99\t20\n")
            classifier=project/"classifier/T1/essentiality_predictions.T1.csv"
            classifier.write_text("Standard name,RF - G4 - ess. for FPR 0.100\na,No\nb,No\nc,Yes\n")
            d=load_library(cfg,["p1","p2"],classifier,weights=[1,3],universe="nonessential",basis="reads",catalogue_mode="mixture")
            self.assertEqual(d["universes"],dict(all=3,nonessential=2,observed=2,informative=2))
            self.assertEqual(d["unique_sites"],5)
            self.assertEqual(d["informative_sites"],2)
            self.assertEqual(d["pools"][0]["ambiguous_sites"],1)
            self.assertAlmostEqual(d["genes"][0]["p"],.25*5/20)
            self.assertAlmostEqual(d["genes"][1]["p"],.75*80/100)
            self.assertEqual(d["genes"][0]["informative_sites"],1)
            conditional=load_library(cfg,["p1","p2"],classifier,weights=[1,3],universe="nonessential",conditioning="conditional")
            self.assertAlmostEqual(sum(g["p"] for g in conditional["genes"]),1)
            all_genes=load_library(cfg,["p1"],universe="all")
            self.assertEqual(all_genes["unknown_genes"],3)
            classified_all=load_library(cfg,["p1"],classifier,universe="all")
            self.assertEqual(classified_all["nonessential_without_informative"],1)
            self.assertEqual(classified_all["informative_genes"],1)
            with self.assertRaises(ValueError):load_library(cfg,["p1"],universe="nonessential")
            with self.assertRaisesRegex(ValueError,"parent/control"):
                load_library(cfg,["parent-treated"],universe="all")
            other=root/"output/projects/second"
            (other/"config").mkdir(parents=True)
            cfg2=other/"config/project.yaml"
            cfg2.write_text(json.dumps(dict(output_project_dir="output/projects/second",sample_sheet="samples.csv",reference=dict(gff="genes.gff"))))
            (root/"p3.txt").write_text("Chromosome\tHit position\tHit count\nchr1\t2\t10\nchr1\t3\t10\nchr1\t99\t80\n")
            sources=[dict(project=str(cfg),samples=["p1","p2"]),dict(project=str(cfg2),samples=["p3"])]
            union=load_library(cfg,[],classifier,parent_sources=sources,universe="nonessential")
            self.assertEqual(union["unique_sites"],6)
            self.assertEqual(union["repeated_site_observations"],2)
            self.assertEqual(union["selected_informative_sites"],3)
            self.assertAlmostEqual(union["genes"][0]["p"],2/6)
            self.assertAlmostEqual(union["genes"][1]["p"],1/6)
            self.assertAlmostEqual(union["genes"][0]["p_reads"],(5/20+0+20/100)/3)
            with self.assertRaisesRegex(ValueError,"only once"):
                load_library(cfg,[],classifier,parent_sources=sources+sources)
            (root/"other.gff").write_text((root/"genes.gff").read_text()+"# different reference\n")
            cfg2.write_text(json.dumps(dict(output_project_dir="output/projects/second",sample_sheet="samples.csv",reference=dict(gff="other.gff"))))
            with self.assertRaisesRegex(ValueError,"matching reference"):
                load_library(cfg,[],classifier,parent_sources=sources)


if __name__=="__main__":unittest.main()
