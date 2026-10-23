"""Argon planning, qsub contracts and failure gates; no cluster needed."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('argon_pipeline', ROOT/'pipeline/hpc/argon_pipeline.py')
argon = importlib.util.module_from_spec(spec)
spec.loader.exec_module(argon)


@pytest.fixture
def project(tmp_path, monkeypatch):
    root = tmp_path/'checkout with spaces'
    cfg = root/'output/projects/demo/config/project.yaml'
    cfg.parent.mkdir(parents=True)
    sheet = cfg.with_name('sample_sheet.csv')
    sheet.write_text('sample,include\na,True\n')
    cfg.write_text(yaml.safe_dump(dict(project_id='demo', repo_root=str(root), species='glabrata',
        threads=8, output_project_dir='output/projects/demo', output_export_dir='output/exports/demo',
        sample_sheet=str(sheet), reference={})))
    scripts = root/'scripts/local'
    scripts.mkdir(parents=True)
    stub = ('import os,sys,json\nfrom pathlib import Path\n'
            'Path(os.environ["CAPTURE"]).write_text(json.dumps(sys.argv))\n'
            'raise SystemExit(int(os.environ.get("TEST_EXIT","0")))\n')
    for name in [*argon._command.__globals__['CLI'].values(), 'ytab.py']:
        (scripts/name).write_text(stub)
    monkeypatch.setattr(argon, 'ROOT', root)
    monkeypatch.setenv('CAPTURE', str(tmp_path/'capture.json'))
    monkeypatch.setenv('NSLOTS', '4')
    return root, cfg


def submit(project, monkeypatch, stages='mapfastq,create_hit_file,summary', fail_at=None):
    root, cfg = project
    calls = []

    def qsub(cmd, **kw):
        calls.append((cmd, kw))
        if len(calls) == fail_at:
            return SimpleNamespace(returncode=1, stdout='', stderr='queue unavailable')
        return SimpleNamespace(returncode=0, stdout=f'{1000+len(calls)}\n', stderr='')

    monkeypatch.setattr(argon.subprocess, 'run', qsub)
    rc = argon.main(['submit', '--project-config', str(cfg), '--stages', stages])
    active = cfg.parents[1]/'manifests/argon/active.json'
    plan = Path(json.loads(active.read_text())['plan']) if active.exists() else None
    return rc, plan, calls


def worker(plan, stage):
    return argon.main(['execute', '--plan', str(plan), '--stage', stage])


def test_profiles_keep_existing_scientific_stages():
    for profile in argon.PROFILES:
        assert argon.choose_stages(profile) == [*argon.PROFILES[profile], 'report', 'export']
    assert argon.choose_stages('all', 'summary,mapfastq') == ['mapfastq', 'summary']
    for names in ['unknown', 'summary,summary', ',']:
        with pytest.raises(ValueError):
            argon.choose_stages('core', names)


def test_master_qsub_chain_and_unique_logs(project, monkeypatch):
    rc, plan, calls = submit(project, monkeypatch)
    assert rc == 0
    data = json.loads(plan.read_text())
    assert [j['job_id'] for j in data['jobs']] == ['1001', '1002', '1003']
    for i, (cmd, kw) in enumerate(calls):
        assert cmd[:2] == ['qsub', '-terse']
        assert cmd[cmd.index('-pe')+1:cmd.index('-pe')+3] == ['smp', '8']
        assert cmd[-4:] == ['--plan', str(plan), '--stage', data['stages'][i]]
        assert kw['cwd'] == project[0] and not kw.get('shell')
        assert '-V' not in cmd
        if i:
            assert cmd[cmd.index('-hold_jid')+1] == str(1000+i)
        else:
            assert '-hold_jid' not in cmd
    assert len({c[0][c[0].index('-o')+1] for c in calls}) == 3


def test_workers_use_allocated_slots_and_shared_scripts(project, monkeypatch):
    rc, plan, calls = submit(project, monkeypatch)
    for stage in ['mapfastq', 'create_hit_file', 'summary']:
        assert worker(plan, stage) == 0
        record = json.loads((plan.parent/f'{stage}.json').read_text())
        assert record['status'] == 'success'
        cmd = record['command']
        assert cmd[cmd.index('--threads')+1] == '4'
        assert '--project-config' in cmd
        assert cmd[0] == sys.executable
    assert (project[1].parents[1]/'mapfastq').exists() is False


def test_failure_blocks_all_downstream_and_preserves_exit(project, monkeypatch):
    rc, plan, _ = submit(project, monkeypatch)
    monkeypatch.setenv('TEST_EXIT', '7')
    assert worker(plan, 'mapfastq') == 7
    captured = Path(os.environ['CAPTURE']).read_text()
    assert worker(plan, 'create_hit_file') == 1
    assert worker(plan, 'summary') == 1
    assert Path(os.environ['CAPTURE']).read_text() == captured
    assert json.loads((plan.parent/'summary.json').read_text())['status'] == 'blocked'
    assert worker(plan, 'mapfastq') == 1  # Receipts cannot be overwritten.


def test_missing_receipt_after_qdel_blocks_next_stage(project, monkeypatch):
    rc, plan, _ = submit(project, monkeypatch)
    assert worker(plan, 'create_hit_file') == 1
    assert not Path(os.environ['CAPTURE']).exists()
    assert 'missing' in json.loads((plan.parent/'create_hit_file.json').read_text())['reason']


def test_config_changes_after_submission_fail_closed(project, monkeypatch):
    rc, plan, _ = submit(project, monkeypatch)
    project[1].with_name('sample_sheet.csv').write_text('sample,include\nb,True\n')
    assert worker(plan, 'mapfastq') == 1
    assert not Path(os.environ['CAPTURE']).exists()
    assert json.loads((plan.parent/'mapfastq.json').read_text())['status'] == 'failed'


def test_new_submission_does_not_reuse_receipts_or_overlap(project, monkeypatch):
    rc, plan, calls = submit(project, monkeypatch)
    before = len(calls)
    assert argon.main(['submit', '--project-config', str(project[1]), '--stages', 'summary']) == 1
    assert len(calls) == before
    assert argon.main(['submit', '--project-config', str(project[1]), '--stages', 'summary', '--replace-active']) == 0
    active = Path(json.loads((plan.parent.parent/'active.json').read_text())['plan'])
    assert active != plan and not (active.parent/'summary.json').exists()


def test_partial_submission_failure_records_already_submitted_jobs(project, monkeypatch):
    rc, plan, calls = submit(project, monkeypatch, fail_at=2)
    assert rc == 1 and len(calls) == 2
    data = json.loads(plan.read_text())
    assert data['submission_status'] == 'partial_failure'
    assert data['jobs'][0]['job_id'] == '1001'


def test_preview_no_writes_or_qsub(project, monkeypatch, tmp_path):
    fastqs = tmp_path/'fastqs/sub/run'
    fastqs.mkdir(parents=True)
    (fastqs/'sample_1.fastq.gz').write_bytes(b'fixture')
    monkeypatch.setattr(argon.subprocess, 'run', lambda *a, **k: pytest.fail('qsub called in preview'))
    assert argon.main(['submit', '--project-id', 'new', '--fastq-dir', str(fastqs.parents[1]), '--dry-run']) == 0
    assert not (project[0]/'output/projects/new').exists()
    assert argon.main(['submit', '--project-config', str(project[1]), '--profile', 'all', '--dry-run']) == 0
    assert not (project[1].parents[1]/'manifests').exists()


def test_nested_inputs_duplicate_rejection_and_symlink_initialization(project, tmp_path, monkeypatch):
    from ytab.pipeline import project_config
    source = tmp_path/'reads'
    (source/'parent/run').mkdir(parents=True)
    read = source/'parent/run/SRR1_1.fastq.gz';read.write_bytes(b'fixture')
    (source/'parent/run/SRR1_2.fastq.gz').write_bytes(b'fixture')
    monkeypatch.setattr(project_config, 'resolve_reference', lambda *a: SimpleNamespace(to_dict=lambda: {}))
    monkeypatch.setattr(argon, 'validate_project_config', lambda config: [])
    cfg = project[0]/'output/projects/new/config/project.yaml'
    cfg.parent.mkdir(parents=True)
    args = SimpleNamespace(fastq_dir=source, project_id='new', species='glabrata', threads=8, sample_sheet=None)
    argon.initialize(args, cfg, argon.catalog(source))
    view = cfg.parent/'argon_fastq_inputs'
    assert (view/read.name).is_symlink() and (view/read.name).resolve() == read
    loaded = argon.load_project_for_mapping(cfg)
    assert loaded['samples'][0]['sample'] == 'SRR1'
    assert loaded['samples'][0]['layout'] == 'paired'
    (source/'other').mkdir()
    (source/'other'/read.name).write_bytes(b'duplicate')
    with pytest.raises(ValueError, match='Duplicate FASTQ'):
        argon.catalog(source)


def test_metadata_csv_merges_roles_not_file_paths(project, tmp_path, monkeypatch):
    from ytab.pipeline import project_config
    source = tmp_path/'reads';source.mkdir()
    read = source/'SRR1.fastq.gz';read.write_bytes(b'fixture')
    sheet = tmp_path/'metadata.csv'
    sheet.write_text('sample,condition,pool,background,fastq_1\nSRR1,parent,1,BG2,/wrong.fastq\n')
    monkeypatch.setattr(project_config, 'resolve_reference', lambda *a: SimpleNamespace(to_dict=lambda: {}))
    monkeypatch.setattr(argon, 'validate_project_config', lambda config: [])
    cfg = project[0]/'output/projects/new/config/project.yaml'
    args = SimpleNamespace(fastq_dir=source, project_id='new', species='glabrata', threads=8, sample_sheet=sheet)
    argon.initialize(args, cfg, argon.catalog(source))
    row = argon.load_project_for_mapping(cfg)['samples'][0]
    assert row['condition'] == 'parent' and row['pool'] == '1'
    assert row['fastq_1'] == str(read) and row['include'] == 'True'


def test_command_parity_for_all_stages_and_existing_classifier(project):
    for stage in argon.STAGES:
        cmd, cwd = argon.stage_command(stage, project[1], 8, {})
        assert cwd == project[0] and Path(cmd[1]).parent == project[0]/'scripts/local'
        assert '--force' not in cmd
        if stage == 'classifier':
            assert Path(cmd[1]).name == 'ytab_run_classifier.py'
            assert cmd[cmd.index('--target')+1] == 'recommended'
        if stage == 'treated_vs_parent':
            assert cmd[cmd.index('--input-mode')+1] == 'raw-summary'
            assert '--target' not in cmd
        if stage == 'export':
            assert '--include-hit-files' in cmd and '--include-bams' not in cmd
    cmd, _ = argon.stage_command('mapfastq', project[1], 4,
        dict(junction_mate='R2', require_mate_alignment=False, force=True))
    assert '--no-require-mate-alignment' in cmd and '--force' in cmd
    assert cmd[cmd.index('--junction-mate')+1] == 'R2'


def test_cancelled_worker_stops_child_and_records_receipt(project, monkeypatch):
    rc, plan, _ = submit(project, monkeypatch)
    child = SimpleNamespace(pid=43210, poll=lambda: None)
    counts = []

    def wait(**kwargs):
        counts.append(kwargs)
        if not kwargs:
            raise KeyboardInterrupt
        return -15

    child.wait = wait
    killed = []
    monkeypatch.setattr(argon.subprocess, 'Popen', lambda *a, **k: child)
    monkeypatch.setattr(argon.os, 'killpg', lambda *a: killed.append(a))
    assert worker(plan, 'mapfastq') == 130
    assert killed and json.loads((plan.parent/'mapfastq.json').read_text())['status'] == 'cancelled'


def test_standalone_stage_and_invalid_slots(project, monkeypatch):
    monkeypatch.setenv('FORCE', '1')
    assert argon.main(['execute', '--project-config', str(project[1]), '--stage', 'summary']) == 0
    cmd = json.loads(Path(os.environ['CAPTURE']).read_text())
    assert '--force' in cmd
    monkeypatch.setenv('NSLOTS', '0')
    assert argon.main(['execute', '--project-config', str(project[1]), '--stage', 'summary']) == 1


def test_shell_wrappers_syntax():
    for filename in ['argon_environment.sh', 'submit_pipeline.sh', 'ytab_stage.txt']:
        subprocess.run(['bash', '-n', str(ROOT/'pipeline/hpc'/filename)], check=True)


def test_real_shell_launcher_and_batch_worker_with_fake_qsub(project, tmp_path):
    """Exercise the actual shell/Python boundary, rather than mocking dispatch."""
    import shutil
    root, cfg = project
    hpc = root/'pipeline/hpc';hpc.mkdir(parents=True)
    for filename in ['argon_environment.sh', 'submit_pipeline.sh', 'ytab_stage.txt', 'argon_pipeline.py']:
        shutil.copyfile(ROOT/'pipeline/hpc'/filename, hpc/filename)
    conda = tmp_path/'conda.sh';conda.write_text('conda() { return 0; }\n')
    fake_qsub = tmp_path/'qsub'
    fake_qsub.write_text('#!' + sys.executable + '\nprint("7654321")\n')
    fake_qsub.chmod(0o755)
    env = dict(os.environ, YTAB=str(root), PYTHON_BIN=sys.executable,
               YTAB_ENV_MODE='conda', YTAB_CONDA_SH=str(conda),
               PYTHONPATH=str(ROOT/'src'), SCRATCH_ROOT=str(tmp_path/'scratch'))
    result = subprocess.run(['bash', str(hpc/'submit_pipeline.sh'), '--project-config', str(cfg),
                             '--stages', 'summary', '--qsub-bin', str(fake_qsub)],
                            env=env, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    plan = Path(json.loads((cfg.parents[1]/'manifests/argon/active.json').read_text())['plan'])
    result = subprocess.run(['bash', str(hpc/'ytab_stage.txt'), '--plan', str(plan), '--stage', 'summary'],
                            env=env, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    assert json.loads((plan.parent/'summary.json').read_text())['status'] == 'success'


def test_submission_lock_excludes_second_master(project):
    args = SimpleNamespace(dry_run=False, project_config=project[1])
    with argon.submission_lock(args):
        with pytest.raises(ValueError, match='Another master'):
            with argon.submission_lock(args):
                pytest.fail('second master acquired the lock')


def test_per_stage_slot_override(project, monkeypatch):
    rc, plan, calls = submit(project, monkeypatch, stages='summary')
    assert worker(plan, 'summary') == 0
    assert argon.main(['submit', '--project-config', str(project[1]), '--stages', 'summary',
                       '--stage-threads', 'summary=1']) == 0
    cmd = calls[-1][0]
    assert cmd[cmd.index('-pe')+2] == '1'


def test_export_opt_ins_include_actual_hits_tracks_and_logs(project):
    from ytab.pipeline.export_bundle import plan_export
    root, cfg = project
    base = cfg.parents[1]
    names = ['create_hit_file/a/a_hits.txt', 'logs/summary/a.log',
             'create_hit_file/a/a.insertions.wig', 'mapfastq/a/a.sorted.bam',
             'mapfastq/a/a.fastq.gz']
    for name in names:
        p = base/name;p.parent.mkdir(parents=True, exist_ok=True);p.write_bytes(b'fixture')
    context, compact = plan_export(cfg)
    assert not set(names).intersection(str(p.relative_to(base)) for p, arc in compact)
    context, included = plan_export(cfg, include_logs=True, include_browser_tracks=True, include_hit_files=True)
    paths = {str(p.relative_to(base)) for p, arc in included}
    assert set(names[:3]) <= paths
    assert not set(names[3:]).intersection(paths)
