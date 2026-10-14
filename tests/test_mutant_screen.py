"""Sanger parsing, archive safety, coordinate semantics and result linking."""
import csv
import gzip
import io
import json
import tempfile
import tarfile
import unittest
import zipfile
import pysam
from pathlib import Path
from unittest.mock import Mock, patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ytab.mutant_screen.identify import (load_reads, input_members, junction_position,
    annotated_candidates, link_results, clip_reads, adapters_from_reference, junction_review, junction_status,
    request_fingerprint, cached_result, junction_proximal_candidates, assign_flank_candidates,
    annotate_gene_structure, gene_region)

class MutantScreenTests(unittest.TestCase):
    def test_explicit_cds_exon_intron_and_utr_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"genes.gff"
            p.write_text("##gff-version 3\n"
                "chr\ttest\tgene\t100\t400\t.\t+\t.\tID=gene-g;locus_tag=g\n"
                "chr\ttest\tmRNA\t100\t400\t.\t+\t.\tID=tx;Parent=gene-g\n"
                "chr\ttest\texon\t100\t180\t.\t+\t.\tID=e1;Parent=tx\n"
                "chr\ttest\texon\t320\t400\t.\t+\t.\tID=e2;Parent=tx\n"
                "chr\ttest\tCDS\t120\t180\t.\t+\t0\tID=c;Parent=tx\n"
                "chr\ttest\tCDS\t320\t380\t.\t+\t0\tID=c;Parent=tx\n"
                "chr\ttest\tfive_prime_UTR\t100\t119\t.\t+\t.\tParent=tx\n")
            genes,meta=annotate_gene_structure([dict(gene_id="g",chromosome="chr",start=100,end=400,exons=[[100,400]])],p)
            g=genes[0]
            self.assertEqual(g["cds_blocks"],[(120,180),(320,380)])
            self.assertEqual(gene_region(g,120),"CDS");self.assertEqual(gene_region(g,180),"CDS")
            self.assertEqual(gene_region(g,181),"Intron");self.assertEqual(gene_region(g,319),"Intron")
            self.assertEqual(gene_region(g,110),"UTR");self.assertEqual(gene_region(g,390),"Noncoding exon / UTR")
            self.assertEqual(gene_region(g,401),"Intergenic")
            self.assertEqual(meta["coordinates"],"1-based, inclusive")

    def test_intronic_alignment_not_counted_as_cds_overlap(self):
        g=dict(gene_id="g",display_name="G",chromosome="chr",start=100,end=500,
               cds_blocks=[(100,180),(420,500)],transcript_exons=[(100,180),(420,500)])
        a=self.proximal_alignment();a.reference_start=199;a.query_sequence="A"*100;a.cigartuples=[(0,100)]
        rows,info=assign_flank_candidates(a,[g],[],None)
        self.assertEqual(rows,[]);self.assertEqual(info["overlapping_cds_count"],0)
        self.assertEqual(info["rule"],"noncoding_gene_overlap")
        self.assertEqual(info["gene_span_ids"],["g"])

    def test_confirmed_intronic_coordinate_not_a_cds_hit(self):
        g=dict(gene_id="g",display_name="G",chromosome="chr",start=100,end=500,
               cds_blocks=[(100,180),(420,500)],transcript_exons=[(100,180),(420,500)])
        rows,info=assign_flank_candidates(self.proximal_alignment(),[g],[],250)
        self.assertEqual(rows,[])
        self.assertEqual(info["boundary_gene_regions"],[dict(gene_id="g",region="Intron")])
    def test_single_cds_flank_is_retained_without_proximal_filter(self):
        genes=[dict(gene_id="only",display_name="Only",chromosome="chr",start=250,end=450)]
        rows,info=assign_flank_candidates(self.proximal_alignment(),genes,[],None)
        self.assertEqual([r["gene_id"] for r in rows],["only"])
        self.assertEqual(info["rule"],"single_cds_flank")
        self.assertIn("unconfirmed",rows[0]["context"])

    def test_multi_cds_resolution_uses_junction_facing_orientation(self):
        genes=[dict(gene_id="left",display_name="Left",chromosome="chr",start=90,end=220),
               dict(gene_id="right",display_name="Right",chromosome="chr",start=250,end=450)]
        for reverse,expected in ((False,"left"),(True,"right")):
            rows,info=assign_flank_candidates(self.proximal_alignment(reverse),genes,[],None)
            self.assertEqual([r["gene_id"] for r in rows],[expected])
            self.assertEqual(info["overlapping_cds_count"],2)
            self.assertEqual(info["rule"],"multiple_cds_junction_proximal")

    def test_confirmed_intergenic_junction_does_not_assign_later_cds(self):
        genes=[dict(gene_id="later",display_name="Later",chromosome="chr",start=250,end=450)]
        rows,info=assign_flank_candidates(self.proximal_alignment(),genes,[],100)
        self.assertEqual(rows,[]);self.assertEqual(info["rule"],"confirmed_junction")

    def proximal_alignment(self,reverse=False,cigar=None):
        header=pysam.AlignmentHeader.from_dict({"SQ":[{"SN":"chr","LN":10000}]})
        a=pysam.AlignedSegment(header);a.query_name="read";a.query_sequence="A"*300
        a.reference_id=0;a.reference_start=99;a.flag=16 if reverse else 0
        a.cigartuples=cigar or [(0,300)]
        return a

    def test_proximal_assignment_excludes_downstream_gene(self):
        genes=[dict(gene_id="near",display_name="Near",chromosome="chr",start=100,end=210),
               dict(gene_id="far",display_name="Far",chromosome="chr",start=250,end=450)]
        rows,info=junction_proximal_candidates(self.proximal_alignment(),genes,[],None)
        self.assertEqual([r["gene_id"] for r in rows],["near"])
        self.assertEqual(info["aligned_bases"],100);self.assertFalse(info["confirmed_junction"])
        self.assertIn("Provisional",rows[0]["context"])
        rows,info=junction_proximal_candidates(self.proximal_alignment(True),genes,[],None)
        self.assertEqual([r["gene_id"] for r in rows],["far"])
        self.assertEqual(info["putative_adjacent_base"],399)

    def test_proximal_window_counts_soft_clipping_and_requires_orientation(self):
        genes=[dict(gene_id="g",display_name="G",chromosome="chr",start=100,end=1000)]
        rows,info=junction_proximal_candidates(self.proximal_alignment(cigar=[(4,70),(0,230)]),genes,[],None)
        self.assertEqual(rows,[]);self.assertEqual(info["aligned_bases"],30)
        rows,info=junction_proximal_candidates(self.proximal_alignment(),genes,[],None,default_end="require")
        self.assertEqual(rows,[]);self.assertIn("undetected",info["reason"])

    def test_intergenic_edge_does_not_assign_gene_later_in_window(self):
        genes=[dict(gene_id="downstream",display_name="G",chromosome="chr",start=120,end=450)]
        rows,info=junction_proximal_candidates(self.proximal_alignment(),genes,[],None)
        self.assertEqual(rows,[])

    def test_confirmed_junction_near_gene_boundary_remains_assigned(self):
        genes=[dict(gene_id="edge",display_name="G",chromosome="chr",start=90,end=102)]
        evidence=[dict(adapter="hermes_front_right",start=0,end=33)]
        rows,info=junction_proximal_candidates(self.proximal_alignment(),genes,evidence,100)
        self.assertEqual([r["gene_id"] for r in rows],["edge"])
        self.assertEqual(rows[0]["junction_proximal_aligned_bases"],3)
        rows,info=junction_proximal_candidates(self.proximal_alignment(),genes,evidence,None)
        self.assertEqual(rows,[])

    def test_paired_seq_abi_uses_quality_and_one_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);text=root/"clone.seq";abi=root/"clone.ab1"
            text.write_text(">clone\nAACCGGTT\n");abi.write_bytes(b"ABI fixture")
            record=SeqRecord(Seq("AACCGGTT"),id="clone-vendor")
            record.letter_annotations["phred_quality"]=[1,1,40,40,40,40,1,1]
            with patch("ytab.mutant_screen.identify.SeqIO.read",return_value=record):
                paired=load_reads([text,abi]);reversed_inputs=load_reads([abi,text]);abi_only=load_reads([abi])
            self.assertEqual(len(paired),1);self.assertEqual(paired,reversed_inputs)
            self.assertEqual(paired[0]["mutant"],"clone")
            self.assertEqual(paired[0]["sequence"],"CCGG")
            self.assertEqual(paired[0]["pre_qc_sequence"],"AACCGGTT")
            self.assertTrue(paired[0]["sequence_export_agrees"])
            self.assertEqual(len(paired[0]["input_sources"]),2)
            self.assertTrue(abi_only[0]["quality_available"])
            self.assertFalse(load_reads([text])[0]["quality_available"])
            text.write_text(">clone\nAATCGGTT\n")
            with patch("ytab.mutant_screen.identify.SeqIO.read",return_value=record):conflict=load_reads([text,abi])[0]
            self.assertFalse(conflict["sequence_export_agrees"])
            self.assertEqual(conflict["sequence"],"CCGG")
            self.assertIn("review required",conflict["input_merge"])

    def test_exact_duplicate_inputs_not_duplicate_mutants(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"clone.seq";p.write_text(">clone\nACGT\n")
            self.assertEqual(len(load_reads([p,p])),1)

    def test_fingerprint_and_reuse_require_same_content_settings_and_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/"clone.seq";p.write_text(">clone\nACGT\n")
            reads=load_reads([p]);args=(reads,"project","ref","hermes",[],{"quality":20},[],{})
            key=request_fingerprint(*args)
            self.assertNotEqual(key,request_fingerprint(reads,"project","ref","hermes",[],{"quality":0},[],{}))
            p.write_text(">clone\nACGA\n")
            self.assertNotEqual(key,request_fingerprint(load_reads([p]),*args[1:]))
            run=root/"run";run.mkdir()
            data=dict(status="success",output_dir=str(run),request_fingerprint=key)
            (run/"results.json").write_text(json.dumps(data))
            self.assertIsNone(cached_result(root/"new",key))
            (run/"mutant_identifications.csv").write_text("mutant\nclone\n")
            (run/"genomic_flanks.fasta").write_text(">clone\nACGT\n")
            self.assertTrue(cached_result(root/"new",key)["reused"])
            self.assertFalse((root/"new").exists())

    def test_vendor_seq_and_gzip_preserve_ids_and_no_invented_quality(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"example.seq.gz"
            with gzip.open(p,"wt") as f:f.write(">clone-A 1000 56 478 0.05\nNNNACGTACGTNN\n")
            r=load_reads([p])[0]
            self.assertEqual((r["mutant"],r["sequence"]),("clone-A","ACGTACGT"))
            self.assertFalse(r["quality_available"])
            self.assertEqual(r["end_trimmed_bases"],5)
            self.assertEqual(len(r["source_member_sha256"]),64)

    def test_fastq_end_quality(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"q.fastq";p.write_text("@clone\nACGTACGT\n+\n!!IIII!!\n")
            r=load_reads([p],quality=20)[0]
            self.assertEqual(r["sequence"],"GTAC");self.assertTrue(r["quality_available"])

    def test_tar_zip_and_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seq=b">clone\nACGTACGT\n"
            with tarfile.open(root/"samples.tar.gz","w:gz") as a:
                info=tarfile.TarInfo("folder/a.seq");info.size=len(seq);a.addfile(info,io.BytesIO(seq))
            with zipfile.ZipFile(root/"samples.zip","w") as a:a.writestr("folder/b.fasta",seq)
            self.assertEqual(len(load_reads([root])),2)

    def test_archive_traversal_and_links_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"bad.tar";seq=b">clone\nACGT\n"
            with tarfile.open(p,"w") as a:
                info=tarfile.TarInfo("../outside.seq");info.size=len(seq);a.addfile(info,io.BytesIO(seq))
            with self.assertRaisesRegex(ValueError,"Unsafe"):list(input_members([p]))
            z=Path(tmp)/"bad.zip"
            with zipfile.ZipFile(z,"w") as a:a.writestr("/outside.seq",seq)
            with self.assertRaisesRegex(ValueError,"Unsafe"):list(input_members([z]))

    def test_junction_coordinate_orientations_and_local_clipping(self):
        front=[dict(adapter="hermes_front_right",start=0,end=33)]
        back=[dict(adapter="hermes_back_right",start=100,end=133)]
        a=Mock(reference_start=99,reference_end=200,cigartuples=[(0,101)],is_reverse=False)
        self.assertEqual(junction_position(a,front),100)
        self.assertEqual(junction_position(a,back),200)
        a.is_reverse=True
        self.assertEqual(junction_position(a,front),200)
        self.assertEqual(junction_position(a,back),100)
        a.cigartuples=[(0,98),(4,3)]
        self.assertIsNone(junction_position(a,front))
        self.assertIsNone(junction_position(a,[dict(adapter="hermes_front_right_short",start=5,end=17)]))
        self.assertIsNone(junction_position(a,[]))

    def test_gene_and_intergenic_context_never_claims_nearest_disruption(self):
        genes=[dict(gene_id="g",display_name="Gene",chromosome="chr",start=100,end=200,exons=[(100,140),(160,200)])]
        self.assertEqual(annotated_candidates(genes,"chr",110,130,110)[0]["context"],"Coding exon")
        self.assertEqual(annotated_candidates(genes,"chr",140,170,150)[0]["context"],"Noncoding gene span")
        self.assertIn("nearby",annotated_candidates(genes,"chr",300,400,None)[0]["context"])
        self.assertIn("unconfirmed",annotated_candidates(genes,"chr",120,180,None)[0]["context"])

    def test_short_hermes_overlap_is_a_review_flag_not_alignment_rejection(self):
        evidence=[dict(adapter="hermes_front_right_short",start=5,end=17)]
        a=Mock(reference_start=99,reference_end=200,cigartuples=[(0,101)],is_reverse=False)
        overlap,reason=junction_review(a,evidence)
        self.assertEqual(overlap,12);self.assertIn("<20",reason);self.assertIn("alignment retained",reason)
        self.assertIsNone(junction_position(a,evidence))
        self.assertEqual(junction_position(a,evidence,min_overlap=12),100)

    def test_sufficient_hermes_overlap_is_supported_without_inventing_boundary(self):
        evidence=[dict(adapter="hermes_front_right",start=0,end=33)]
        a=Mock(reference_start=99,reference_end=200,cigartuples=[(0,101)],is_reverse=False)
        self.assertEqual(junction_status(junction_position(a,evidence),evidence),"Hermes junction supported")
        self.assertIn("adjacent genomic base aligns",junction_review(a,evidence)[1])
        a.cigartuples=[(4,3),(0,101)]
        self.assertIsNone(junction_position(a,evidence))
        self.assertEqual(junction_status(None,evidence),"Hermes overlap supported; exact junction position unresolved")
        self.assertNotIn("unconfirmed",junction_review(a,evidence)[1])
        self.assertIn("unconfirmed",junction_status(None,[dict(adapter="hermes_front_right_short",start=0,end=12)]))

    def test_screen_links_use_native_id_and_keep_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/"treated_vs_parent/run/treated_vs_parent_results.csv";p.parent.mkdir(parents=True)
            p.write_text("feature_id,final_call,mean_log2FC\ng,depleted,-2\nother,enriched,1\n")
            rows=link_results(root,[dict(gene_id="g")],p)
            self.assertEqual(len(rows),1);self.assertEqual(rows[0]["mean_log2FC"],"-2")
            self.assertEqual(rows[0]["analysis_or_sample"],"run")

    def test_empty_sequence_can_be_reported_without_running_clipper(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"n.seq";p.write_text("NNNN")
            h=Path(tmp)/"h.fa";h.write_text(">h\n"+"ACGT"*100)
            reads=load_reads([p]);commands=[]
            clip_reads(reads,h,Path(tmp),commands)
            self.assertEqual(commands,[]);self.assertEqual(reads[0]["clean_sequence"],"")

if __name__=="__main__":unittest.main()
