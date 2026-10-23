"""Designated-mate Hermes cleanup, exact endpoints, and real aligner integration."""
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

import dnaio
import pysam
import pytest
import yaml

from ytab.mapping.hermes_preprocess import HermesCleaner, prepare_fastqs, reverse_complement
from ytab.pipeline.mapfastq_runner import build_mapfastq_command, mapping_fingerprint, run_mapfastq_sample, validate_mapping_inputs

ROOT=Path(__file__).resolve().parents[1]
GENOME=''.join(random.Random(117).choices('ACGT',k=3000))
FLANK=GENOME[300:375]
ADAPTER='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'


def record(seq,name='read'):
    return dnaio.SequenceRecord(name,seq,'I'*len(seq))


def write_fastq(path,seqs):
    with dnaio.open(str(path),mode='w') as out:
        for i,seq in enumerate(seqs):
            out.write(record(seq,f'q{i}'))


@pytest.fixture
def cleaner():
    return HermesCleaner()


def test_catalogue_sequencing_primer_ends_exactly_at_boundary(cleaner):
    c=cleaner.catalogue
    primer=c['primers']['Hermes-sequencing']
    assert cleaner.hermes.endswith(primer[-56:])
    assert c['sequencing_design']['Hermes-sequencing']['remaining_hermes_bases']==0
    assert cleaner.hermes[2693:2716] == c['primers']['HR-2-PCR'][-23:]
    assert len(cleaner.hermes[2716:])==33


def test_genomic_start_is_not_trimmed_or_discarded(cleaner):
    result,stats=cleaner.clean(record(FLANK))
    assert result.sequence==FLANK
    assert not stats


@pytest.mark.parametrize('terminal', ['right','left','primer','hr2','hr1'])
def test_terminal_or_retained_primer_is_removed(cleaner,terminal):
    tail={'right':cleaner.hermes[-33:],'left':reverse_complement(cleaner.hermes[:56]),
          'primer':cleaner.catalogue['primers']['Hermes-sequencing'],
          'hr2':cleaner.catalogue['primers']['HR-2-PCR']+cleaner.hermes[-33:],
          'hr1':cleaner.catalogue['primers']['HR-1']+cleaner.hermes[-60:]}[terminal]
    result,stats=cleaner.clean(record(tail+FLANK))
    assert result.sequence==FLANK
    assert stats['trimmed_5p_bases']==len(tail)
    assert stats['hermes_5p_reads']==1


def test_mismatch_in_hermes_tail_is_tolerated(cleaner):
    tail=list(cleaner.hermes[-33:]);tail[8]='C' if tail[8]!='C' else 'A'
    result,stats=cleaner.clean(record(''.join(tail)+FLANK))
    assert result.sequence==FLANK and stats['hermes_5p_reads']==1


def test_adapter_readthrough_only_is_clipped(cleaner):
    result,stats=cleaner.clean(record(FLANK+ADAPTER))
    assert result.sequence==FLANK
    assert stats['trimmed_3p_bases']==len(ADAPTER)


def test_complete_5prime_adapter_can_be_removed(cleaner):
    adapter=cleaner.catalogue['adaptors']['Illumina-P5']
    result,stats=cleaner.clean(record(adapter+FLANK))
    assert result.sequence==FLANK and stats['adapter_5p_reads']==1


def test_no_genomic_quality_trimming_or_short_index_clipping(cleaner):
    seq='CGTGATT'+FLANK
    r=dnaio.SequenceRecord('r',seq,'!'*5+'I'*(len(seq)-5))
    result,stats=cleaner.clean(r)
    assert result.sequence==seq and result.qualities==r.qualities
    assert not stats


def test_explicit_junction_end_normalizes_sequence_and_quality(cleaner):
    seq=cleaner.hermes[-33:]+FLANK
    r=record(reverse_complement(seq))
    result,stats=cleaner.clean(r,read_end='end')
    assert result.sequence==FLANK and stats['reverse_complemented']==1


def test_too_short_and_empty_samples_fail_clearly(tmp_path,cleaner):
    r,stats=cleaner.clean(record(cleaner.hermes[-33:]+'ACGT'))
    assert r is None and stats['too_short']==1
    source=tmp_path/'empty.fq';source.touch()
    with pytest.raises(ValueError,match='No reads'):
        prepare_fastqs([source],[tmp_path/'out.fq'])


def config(tmp_path,prefix='/tmp/no-index'):
    return dict(_repo_root=str(ROOT),threads=2,output_project_dir=str(tmp_path/'project'),
                output_export_dir=str(tmp_path/'exports'),reference={'bowtie2_index_prefix':str(prefix)})


def sample(tmp_path):
    return dict(sample='test',layout='paired',fastq_1=str(tmp_path/'r1.fq'),fastq_2=str(tmp_path/'r2.fq'))


def test_default_r1_and_explicit_r2_without_mate_agreement(tmp_path):
    s=sample(tmp_path);cfg=config(tmp_path)
    cmd=build_mapfastq_command(s,cfg)
    assert '--r1' not in cmd and '--r2' not in cmd
    assert cmd[cmd.index('--input-file-name')+1]==s['fastq_1']
    assert '--hermes-preprocess' in cmd and '--discard-untrimmed' not in cmd
    s['junction_mate']='R2';cmd=build_mapfastq_command(s,cfg)
    assert cmd[cmd.index('--input-file-name')+1]==s['fastq_2']
    assert cmd[cmd.index('--source-junction-mate')+1]=='R2'


def test_explicit_agreement_and_project_defaults(tmp_path):
    cfg=config(tmp_path);cfg['mapping']={'junction_mate':'R2','require_mate_alignment':True}
    s=sample(tmp_path);cmd=build_mapfastq_command(s,cfg)
    assert '--require-mate-alignment' in cmd
    assert cmd[cmd.index('--junction-mate')+1]=='R2'
    s.update(junction_mate='R1',require_mate_alignment=False)
    cmd=build_mapfastq_command(s,cfg)
    assert '--r1' not in cmd and cmd[cmd.index('--input-file-name')+1]==s['fastq_1']
    s['junction_mate']='both'
    with pytest.raises(ValueError,match='one junction'):
        build_mapfastq_command(s,cfg)


def test_empty_sample_sheet_options_do_not_override_project_defaults(tmp_path):
    cfg=config(tmp_path);cfg['mapping']={'junction_mate':'R2','hermes_cleanup':False}
    s=sample(tmp_path);s.update(junction_mate='',hermes_cleanup='')
    cmd=build_mapfastq_command(s,cfg)
    assert '--hermes-preprocess' not in cmd
    assert cmd[cmd.index('--input-file-name')+1]==s['fastq_2']


def test_fingerprint_changes_when_input_or_settings_change(tmp_path):
    s=sample(tmp_path);cfg=config(tmp_path)
    write_fastq(Path(s['fastq_1']),[FLANK]);write_fastq(Path(s['fastq_2']),[FLANK])
    first=mapping_fingerprint(build_mapfastq_command(s,cfg))
    write_fastq(Path(s['fastq_1']),[FLANK,FLANK])
    assert mapping_fingerprint(build_mapfastq_command(s,cfg))!=first


def test_catalogue_content_changes_invalidate_cache(tmp_path):
    s=sample(tmp_path);cfg=config(tmp_path)
    cleaner=HermesCleaner()
    catalogue=tmp_path/'primers.json'
    data=cleaner.catalogue;data['reference']=str(cleaner.reference)
    catalogue.write_text(json.dumps(data));s['primer_catalogue']=str(catalogue)
    command=build_mapfastq_command(s,cfg);first=mapping_fingerprint(command)
    data['adaptors']['custom']='GATTACAGATTACAGATTACAGATTACA'
    catalogue.write_text(json.dumps(data))
    assert mapping_fingerprint(command)!=first
    s['junction_mate']='R2'
    assert mapping_fingerprint(build_mapfastq_command(s,cfg))!=first


def test_pair_names_and_lengths_are_checked(tmp_path):
    r1,r2=tmp_path/'r1.fq',tmp_path/'r2.fq'
    write_fastq(r1,[FLANK,FLANK]);write_fastq(r2,[FLANK])
    with pytest.raises(ValueError,match='unequal'):
        prepare_fastqs([r1,r2],[tmp_path/'one.fq',tmp_path/'two.fq'])


@pytest.fixture(scope='module')
def index(tmp_path_factory):
    if not all(shutil.which(t) for t in ['bowtie2','bowtie2-build','samtools']):
        pytest.skip('Real alignment integration requires Bowtie2 and samtools')
    root=tmp_path_factory.mktemp('hermes-index');fasta=root/'ref.fa'
    fasta.write_text('>chrSynthetic\n'+GENOME+'\n')
    subprocess.run(['bowtie2-build',str(fasta),str(root/'ref')],check=True,capture_output=True)
    return root/'ref'


def run_alignment(tmp_path,index,r1,r2=None,**settings):
    s=sample(tmp_path);cfg=config(tmp_path,index)
    write_fastq(Path(s['fastq_1']),r1)
    if r2 is not None:write_fastq(Path(s['fastq_2']),r2)
    else:s['layout']='single';s['fastq_2']=''
    s.update(settings)
    result=run_mapfastq_sample(s,cfg)
    assert result['status']=='success',Path(result['log_file']).read_text()
    with pysam.AlignmentFile(Path(result['output_dir'])/'test.sorted.bam','rb') as bam:
        alignments=list(bam)
    return result,alignments,s,cfg


def test_real_cleaning_preserves_forward_and_reverse_insertion_endpoints(tmp_path,index,cleaner):
    forward=cleaner.hermes[-33:]+FLANK+ADAPTER
    reverse=reverse_complement(GENOME[600:675])
    result,rows,s,cfg=run_alignment(tmp_path,index,[forward,reverse])
    assert len(rows)==2
    assert {(a.reference_end if a.is_reverse else a.reference_start+1) for a in rows}=={301,675}
    report=result['hermes_preprocessing']
    assert report['mates'][0]['hermes_5p_reads']==1
    assert report['mates'][0]['adapter_3p_reads']==1
    assert report['retained_reads_or_pairs']==2
    assert run_mapfastq_sample(s,cfg)['status']=='skipped'


def test_real_r2_selection_does_not_require_r1_to_map(tmp_path,index):
    result,rows,_,_=run_alignment(tmp_path,index,['N'*75],[FLANK],junction_mate='R2')
    assert len(rows)==1 and rows[0].reference_start==300
    assert result['junction_mate']=='R2'


def test_unselected_mate_is_not_a_required_input(tmp_path,index):
    s=sample(tmp_path);cfg=config(tmp_path,index)
    write_fastq(Path(s['fastq_1']),[FLANK]);cfg['samples']=[s]
    assert not validate_mapping_inputs(cfg)
    s['require_mate_alignment']=True
    assert any('R2 missing' in error for error in validate_mapping_inputs(cfg))
    s.update(require_mate_alignment=False,junction_mate='R2',fastq_1='/absent/r1.fq')
    write_fastq(Path(s['fastq_2']),[FLANK])
    assert not validate_mapping_inputs(cfg)


def test_real_default_keeps_junction_when_other_mate_unmapped(tmp_path,index):
    _,rows,_,_=run_alignment(tmp_path,index,[FLANK],['N'*75])
    assert len(rows)==1


def test_real_concordance_rejects_bad_partner_and_never_double_counts(tmp_path,index):
    mate=reverse_complement(GENOME[500:575])
    _,rows,_,_=run_alignment(tmp_path,index,[FLANK,FLANK],[mate,'N'*75],require_mate_alignment=True)
    assert len(rows)==1 and rows[0].is_read1 and rows[0].is_proper_pair
    assert rows[0].reference_start==300


def test_real_concordance_r2_only_emits_r2_endpoint(tmp_path,index):
    mate=reverse_complement(GENOME[500:575])
    _,rows,_,_=run_alignment(tmp_path,index,[FLANK],[mate],junction_mate='R2',require_mate_alignment=True)
    assert len(rows)==1 and rows[0].is_read2 and rows[0].is_proper_pair
    assert rows[0].reference_end==575


def test_preview_cannot_turn_old_bam_into_valid_cache(tmp_path,index):
    _,_,s,cfg=run_alignment(tmp_path,index,[FLANK])
    s['min_read_length']=30
    preview=run_mapfastq_sample(s,cfg,dry_run=True)
    assert preview['status']=='dry-run'
    assert run_mapfastq_sample(s,cfg)['status']=='success'


def test_real_unified_cli_project_overrides_use_r2_without_r1(tmp_path,index):
    s=sample(tmp_path);cfg=config(tmp_path,index)
    write_fastq(Path(s['fastq_2']),[FLANK])
    cfg['samples']=[s];cfg['repo_root']=str(ROOT)
    project=tmp_path/'project.yaml';project.write_text(yaml.safe_dump(cfg))
    result=subprocess.run([sys.executable,str(ROOT/'scripts/local/ytab.py'),'mapfastq',
        '--project-config',str(project),'--junction-mate','R2','--no-require-mate-alignment','--threads','2'],
        capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    with pysam.AlignmentFile(tmp_path/'project/mapfastq/test/test.sorted.bam','rb') as bam:
        rows=list(bam)
    assert len(rows)==1 and rows[0].reference_start==300


def test_transposon_designs_set_trimming_and_hermes_default_is_unchanged(tmp_path):
    cfg=config(tmp_path);s=sample(tmp_path)
    default=build_mapfastq_command(dict(s),cfg)
    assert build_mapfastq_command(dict(s,transposon='hermes',library_design='splinkerette_pcr'),cfg)==default
    cmd=build_mapfastq_command(dict(s,transposon='ac_ds',library_design='ds_primer_tail'),cfg)
    assert '--hermes-preprocess' not in cmd and '--discard-untrimmed' in cmd
    assert cmd[cmd.index('--trim-5p-seq')+1]=='GTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA'
    cmd=build_mapfastq_command(dict(s,transposon='ac_ds',library_design='satay_restriction'),cfg)
    assert '--trim-5p-seq' not in cmd and cmd.count('--trim-3p-seq')==3
    cfg['mapping']={'transposon':'hermes','library_design':'terminal_10bp'}
    cmd=build_mapfastq_command(dict(s),cfg)
    assert cmd[cmd.index('--trim-5p-seq')+1]=='^AAGTTCTCTG'
    # an explicit setting overrides the design's
    cmd=build_mapfastq_command(dict(s,discard_untrimmed=False),cfg)
    assert '--discard-untrimmed' not in cmd


def test_unknown_transposon_or_design_fails_clearly(tmp_path):
    cfg=config(tmp_path);s=sample(tmp_path)
    with pytest.raises(ValueError,match='Unknown transposon'):
        build_mapfastq_command(dict(s,transposon='himar1'),cfg)
    with pytest.raises(ValueError,match='Unknown library design'):
        build_mapfastq_command(dict(s,transposon='ac_ds',library_design='nope'),cfg)


def test_mapfastq_trims_every_3prime_sequence_in_one_pass(tmp_path):
    if not shutil.which('cutadapt'):
        pytest.skip('cutadapt not installed')
    sys.path.insert(0,str(ROOT/'src/ytab/mapping'))
    import MapFastq
    source=tmp_path/'in.fq';out=tmp_path/'out.fq'
    reads=['ACGTACGTACGTACGTGATCGTATCGGTTTTCGATTACCGTATTTATCC','TTTTGGGGCCCCAAAACATGTGTGATTTTACCGAACAAAAATACCGG']
    source.write_text(''.join(f'@r{i}\n{r}\n+\n{"I"*len(r)}\n' for i,r in enumerate(reads)))
    MapFastq.TrimAdapter3_SE(str(source),['GTATCGGTTTTCGATTACCGTATTTATCCCGTTCG','TGTGATTTTACCGAACAAAAATACCGGTTCCCGTCCG'],str(out),shutil.which('cutadapt'))
    seqs=out.read_text().splitlines()[1::4]
    assert seqs==['ACGTACGTACGTACGTGATC','TTTTGGGGCCCCAAAACATG']
