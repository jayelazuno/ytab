#!/usr/bin/env python3
"""Personal Argon/SGE submission and success-gated execution of existing YTAB CLIs."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
from datetime import datetime, timezone
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from ytab.pipeline.orchestrator import ORDER, PROFILES, _command
from ytab.pipeline.mapfastq_runner import load_project_for_mapping
from ytab.pipeline.project_config import create_project_config, validate_project_config, update_fastq_base
from ytab.pipeline.sample_discovery import FASTQ_RE

STAGES = [*ORDER, 'report', 'export']
SAFE_ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f'.{uuid4().hex}.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def catalog(directory):
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f'FASTQ directory does not exist: {directory}')
    files = sorted(p for p in directory.rglob('*') if p.is_file() and FASTQ_RE.search(p.name))
    if not files:
        raise ValueError(f'No FASTQs found: {directory}')
    by_name = {}
    for path in files:
        if path.name in by_name:
            raise ValueError(f'Duplicate FASTQ basename; prepare an unambiguous input directory: {path.name}')
        by_name[path.name] = path
    return files


def choose_stages(profile, selected=None):
    names = [s.strip() for s in selected.split(',')] if selected else [*PROFILES[profile], 'report', 'export']
    if not names or len(set(names)) != len(names) or set(names) - set(STAGES):
        raise ValueError('Stages must be unique names from: ' + ', '.join(STAGES))
    # A conservative serial chain is intentional: no concurrent writes to project status.
    return [s for s in STAGES if s in names]


def stage_command(stage, config_path, threads, options):
    if stage not in STAGES:
        raise ValueError(f'Unknown stage: {stage}')
    if threads < (2 if stage in {'mapfastq', 'reference'} else 1):
        raise ValueError(f'Insufficient slots for {stage}: {threads}')
    config = load_project_for_mapping(Path(config_path))
    root = Path(config['_repo_root'])
    context = dict(root=root, project_config=Path(config_path), config=config)
    if stage in ORDER:
        cmd = _command(stage, context, options.get('target', 'recommended'),
                       options.get('analysis_id', 'treated_vs_parent_raw_cpm'),
                       threads, options.get('force', False), False, False)
    else:
        cmd = [sys.executable, str(root/'scripts/local/ytab.py'), stage,
               '--project-config', str(config_path)]
        if stage == 'export':
            cmd += ['--include-hit-files', '--include-browser-tracks', '--include-logs']
        if options.get('force'):
            cmd.append('--force')
    if stage == 'mapfastq':
        if options.get('junction_mate'):
            cmd += ['--junction-mate', options['junction_mate']]
        if options.get('require_mate_alignment') is not None:
            cmd.append('--require-mate-alignment' if options['require_mate_alignment'] else '--no-require-mate-alignment')
    return cmd, root


def snapshot(config_path):
    config = load_project_for_mapping(config_path)
    paths = [Path(config_path)]
    if config.get('sample_sheet'):
        p = Path(config['sample_sheet'])
        paths.append(p if p.is_absolute() else Path(config['_repo_root'])/p)
    return {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def execute(args):
    stage = args.stage or os.environ.get('STAGE', '')
    receipt = None
    record = dict(stage=stage, status='running', started_at=now(), job_id=os.environ.get('JOB_ID'))
    if args.plan:
        plan_path = args.plan.resolve()
        plan = json.loads(plan_path.read_text())
        if stage not in plan['stages']:
            raise ValueError('Stage is not part of this submission plan.')
        receipt = plan_path.parent/f'{stage}.json'
        if receipt.exists():
            raise ValueError('This stage already has a receipt; submit a new run instead of reusing a plan.')
        # SGE hold_jid releases after failure/qdel too. Require this run's success receipt.
        prior = plan['stages'][:plan['stages'].index(stage)]
        for dependency in prior:
            path = plan_path.parent/f'{dependency}.json'
            state = json.loads(path.read_text()).get('status') if path.is_file() else 'missing'
            if state != 'success':
                record.update(status='blocked', finished_at=now(), reason=f'{dependency}: {state}')
                write_json(receipt, record)
                print('BLOCKED: ' + record['reason'], flush=True)
                return 1
        config_path = Path(plan['project_config'])
        options = plan['options']
    else:
        value = args.project_config or os.environ.get('PROJECT_CONFIG')
        if not value:
            project = os.environ.get('PROJECT_ID', '')
            if not SAFE_ID.fullmatch(project):
                raise ValueError('Supply PROJECT_CONFIG or a valid PROJECT_ID.')
            value = ROOT/'output/projects'/project/'config/project.yaml'
        config_path = Path(value).expanduser().resolve()
        force = os.environ.get('FORCE', '0')
        if force not in {'0', '1'}:
            raise ValueError('FORCE must be 0 or 1.')
        options = dict(force=force == '1', target=os.environ.get('TARGET', 'recommended'),
                       analysis_id=os.environ.get('ANALYSIS_ID', 'treated_vs_parent_raw_cpm'),
                       junction_mate=os.environ.get('JUNCTION_MATE') or None)
        agreement = os.environ.get('REQUIRE_MATE_ALIGNMENT', '')
        if agreement not in {'', '0', '1'}:
            raise ValueError('REQUIRE_MATE_ALIGNMENT must be 0 or 1.')
        options['require_mate_alignment'] = None if not agreement else agreement == '1'
    child = None
    old_handlers = {}

    def cancelled(signum, frame):
        raise KeyboardInterrupt

    try:
        if receipt:
            write_json(receipt, record)
            for filename, expected in plan['inputs'].items():
                if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != expected:
                    raise ValueError(f'Configuration/sample sheet changed after submission: {filename}')
        threads = int(os.environ.get('NSLOTS', '0'))
        if options.get('junction_mate') not in {None, 'R1', 'R2'}:
            raise ValueError('JUNCTION_MATE must be R1 or R2.')
        cmd, root = stage_command(stage, config_path, threads, options)
        record['command'] = cmd
        print('Command: ' + shlex.join(cmd), flush=True)
        for sig in (signal.SIGTERM, signal.SIGINT):
            old_handlers[sig] = signal.signal(sig, cancelled)
        child = subprocess.Popen(cmd, cwd=root, start_new_session=True)
        rc = child.wait()
        record.update(status='success' if rc == 0 else 'failed', returncode=rc)
        return rc if rc >= 0 else 128 - rc
    except KeyboardInterrupt:
        if child and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        record.update(status='cancelled', returncode=130)
        return 130
    except Exception as exc:
        record.update(status='failed', reason=str(exc), returncode=1)
        print(f'ERROR: {exc}', file=sys.stderr, flush=True)
        return 1
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        record['finished_at'] = now()
        if receipt:
            write_json(receipt, record)
        print(f"Finished {stage}: {record['status']}", flush=True)


def initialize(args, config_path, files):
    # The shared initializer discovers flat inputs. A symlink-only view supports
    # nested SRA downloads without copying reads or altering the local pipeline.
    nested = any(p.parent != args.fastq_dir.expanduser().resolve() for p in files)
    source = args.fastq_dir.expanduser().resolve()
    if nested:
        source = config_path.parent/'argon_fastq_inputs'
        source.mkdir(parents=True, exist_ok=False)
        for path in files:
            (source/path.name).symlink_to(path.resolve())
    config = create_project_config(args.project_id, source, args.species, ROOT, args.threads)
    if any('Multiple FASTQs' in str(r.get('warnings', '')) for r in config['samples']):
        raise ValueError('Ambiguous sample/read filenames; fix the sample sheet before submitting.')
    if args.sample_sheet:
        import pandas as pd
        sheet = pd.read_csv(args.sample_sheet, dtype=str, keep_default_na=False)
        discovered = {r['sample']: r for r in config['samples']}
        if 'sample' not in sheet or sheet['sample'].duplicated().any() or set(sheet['sample']) != set(discovered):
            raise ValueError('Sample sheet must contain each discovered sample ID exactly once.')
        for i, row in sheet.iterrows():
            for key in ('fastq_1', 'fastq_2', 'layout'):
                # Sequencing files/layout come from discovery; user sheet adds roles/metadata.
                sheet.loc[i, key] = discovered[row['sample']][key]
        if 'include' not in sheet:
            sheet['include'] = True
        sheet.to_csv(config_path.with_name('sample_sheet.csv'), index=False)
        config['samples'] = sheet.to_dict(orient='records')
    warnings = validate_project_config(config)
    for warning in warnings:
        print(warning)
    if any(w.startswith('ERROR:') for w in warnings):
        raise ValueError('Project validation failed. Fix its configuration before submitting.')


@contextmanager
def submission_lock(args):
    """Prevent simultaneous masters from submitting overlapping project chains."""
    if args.dry_run:
        yield
        return
    if args.project_config:
        config = load_project_for_mapping(args.project_config)
        project = Path(config['output_project_dir'])
        if not project.is_absolute():
            project = Path(config['_repo_root'])/project
    else:
        if not args.project_id or not SAFE_ID.fullmatch(args.project_id) or not args.fastq_dir:
            raise ValueError('New project requires a safe --project-id and --fastq-dir.')
        catalog(args.fastq_dir)  # Fail on invalid inputs before creating lock directories.
        project = ROOT/'output/projects'/args.project_id
    directory = project/'manifests/argon'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory/'.submission.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('Another master is submitting this project.') from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def submit(args):
    if args.project_config and args.project_id:
        raise ValueError('Use either --project-config or --project-id, not both.')
    if args.threads < 2:
        raise ValueError('--threads must be at least 2; use --stage-threads for single-slot stages.')
    stages = choose_stages(args.profile, args.stages)
    slots = {s: args.threads for s in stages}
    for override in args.stage_threads:
        name, value = override.split('=', 1)
        if name not in slots:
            raise ValueError(f'Thread override is not in the selected plan: {name}')
        slots[name] = int(value)
    for name, value in slots.items():
        if value < (2 if name in {'mapfastq', 'reference'} else 1):
            raise ValueError(f'Insufficient slots: {name}={value}')
    if args.project_config:
        if args.sample_sheet:
            raise ValueError('Existing project: edit its sample sheet rather than supplying --sample-sheet.')
        config_path = args.project_config.expanduser().resolve()
        config = load_project_for_mapping(config_path)
        project_id = config['project_id']
        if config.get('start_stage') == 'create_hit_file' and not args.stages:
            stages = [s for s in stages if s not in {'reference', 'mapfastq', 'create_hit_file'}]
            print('Imported-hit project: mapping and hit creation are already supplied; starting downstream.')
    else:
        if not args.project_id or not SAFE_ID.fullmatch(args.project_id) or not args.fastq_dir:
            raise ValueError('New project requires a safe --project-id and --fastq-dir.')
        project_id = args.project_id
        config_path = ROOT/'output/projects'/project_id/'config/project.yaml'
        if config_path.exists():
            raise ValueError(f'Project exists; use --project-config {config_path} to preserve its settings.')
        files = catalog(args.fastq_dir)
        print(f'FASTQs found: {len(files)} (including nested directories); sample IDs follow filenames.')
        if not args.dry_run:
            initialize(args, config_path, files)
        config = dict(output_project_dir=str(config_path.parents[1]), _repo_root=str(ROOT))
    if not SAFE_ID.fullmatch(project_id):
        raise ValueError('Project ID must be a safe filename identifier.')
    project = Path(config['output_project_dir'])
    if not project.is_absolute():
        project = Path(config['_repo_root'])/project
    options = dict(force=args.force, target=args.target, analysis_id=args.analysis_id,
                   junction_mate=args.junction_mate, require_mate_alignment=args.require_mate_alignment)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid4().hex[:8]
    run_dir = project/'manifests/argon'/stamp
    log_dir = project/'logs/argon'/stamp
    plan_path = run_dir/'plan.json'
    plan = dict(version=1, created_at=now(), project_config=str(config_path),
                profile=args.profile, stages=stages, options=options, slots=slots, jobs=[])
    env = dict(os.environ, YTAB=str(ROOT))
    # Copy the environment that was explicitly activated by the launcher, not base.
    env.setdefault('YTAB_ENV_MODE', 'conda')
    if args.dry_run:
        if args.project_config and args.fastq_dir:
            print(f'Would update only the FASTQ base: {args.fastq_dir.expanduser().resolve()}')
        for stage in stages:
            print(f"{stage}: {slots[stage]} slots; waits for previous stage's successful receipt")
        print(f'Would write submission plan: {plan_path}\nDRY RUN: no initialization, files or qsub jobs created.')
        return 0
    active = project/'manifests/argon/active.json'
    if active.is_file() and not args.replace_active:
        old = Path(json.loads(active.read_text())['plan'])
        old_plan = json.loads(old.read_text())
        if any(not (old.parent/f'{s}.json').is_file() or
               json.loads((old.parent/f'{s}.json').read_text()).get('status') == 'running'
               for s in old_plan['stages']):
            raise ValueError('A previous run may still be queued/running. Cancel/check it before using --replace-active.')
    run_dir.mkdir(parents=True)
    log_dir.mkdir(parents=True)
    if args.project_config and args.fastq_dir:
        update_fastq_base(config_path, args.fastq_dir)
    plan['inputs'] = snapshot(config_path)
    write_json(plan_path, plan)
    write_json(active, dict(plan=str(plan_path)))
    previous = None
    try:
        for stage in stages:
            cmd = [args.qsub_bin, '-terse', '-cwd', '-q', args.queue, '-pe', args.parallel_env,
                   str(slots[stage]), '-N', f'YTAB_{stage}', '-o', str(log_dir/f'{stage}.log'), '-j', 'y']
            if args.email:
                cmd += ['-M', args.email, '-m', 'e']
            # Export only wrapper configuration, not the entire login environment.
            keys = ['YTAB', 'YTAB_ENV_MODE', 'YTAB_CONDA_SH', 'YTAB_CONDA_ENV',
                    'YTAB_MODULES', 'PYTHON_BIN', 'SCRATCH_ROOT']
            exported = [k for k in keys if env.get(k)]
            cmd += ['-v', ','.join(exported)]
            if previous:
                cmd += ['-hold_jid', previous]
            cmd += [str(ROOT/'pipeline/hpc/ytab_stage.txt'), '--plan', str(plan_path), '--stage', stage]
            print(shlex.join(cmd), flush=True)
            result = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr.strip() or result.stdout.strip() or 'qsub failed')
            lines = result.stdout.strip().splitlines()
            job_id = lines[-1] if lines else ''
            if not re.fullmatch(r'\d+', job_id):
                raise RuntimeError(f'Unexpected qsub -terse response: {result.stdout!r}; check qstat before resubmitting.')
            plan['jobs'].append(dict(stage=stage, job_id=job_id, command=cmd))
            write_json(plan_path, plan)
            previous = job_id
    except BaseException:
        plan['submission_status'] = 'partial_failure'
        write_json(plan_path, plan)
        print('Submission interrupted. Already submitted jobs are NOT automatically deleted. Check qstat; IDs are in ' + str(plan_path), file=sys.stderr)
        raise
    plan['submission_status'] = 'submitted'
    write_json(plan_path, plan)
    print(f'Submitted {len(stages)} jobs. Plan: {plan_path}\nLogs: {log_dir}')
    print('Monitor: qstat -u "$USER". Per-stage JSON receipts record success, failure or blocked status.')
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('submit', help='Initialize/reuse a project and submit a success-gated SGE chain')
    p.add_argument('--project-id')
    p.add_argument('--project-config', type=Path)
    p.add_argument('--fastq-dir', type=Path)
    p.add_argument('--species', default='glabrata')
    p.add_argument('--sample-sheet', type=Path, help='Optional metadata CSV; must match discovered sample IDs')
    p.add_argument('--profile', choices=PROFILES, default='core')
    p.add_argument('--stages', help='Explicit comma-separated stage subset, in standard pipeline order')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--stage-threads', action='append', default=[], metavar='STAGE=N')
    p.add_argument('--queue', default='UI')
    p.add_argument('--parallel-env', default='smp')
    p.add_argument('--email', default='jayelazuno@uiowa.edu')
    p.add_argument('--qsub-bin', default='qsub')
    p.add_argument('--target', default='recommended')
    p.add_argument('--analysis-id', default='treated_vs_parent_raw_cpm')
    p.add_argument('--junction-mate', choices=['R1', 'R2'])
    p.add_argument('--require-mate-alignment', action=argparse.BooleanOptionalAction, default=None)
    p.add_argument('--force', action='store_true', help='Forward force to existing scientific runners; never rewrite project config')
    p.add_argument('--dry-run', action='store_true', help='Read-only preview; no project creation or qsub')
    p.add_argument('--replace-active', action='store_true', help='Allow a new chain after checking/cancelling a stale active run')
    p = sub.add_parser('execute', help='Batch worker or standalone stage (normally invoked by qsub)')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--plan', type=Path)
    group.add_argument('--project-config', type=Path)
    p.add_argument('--stage', choices=STAGES)
    args = parser.parse_args(argv)
    try:
        if args.action == 'submit':
            with submission_lock(args):
                return submit(args)
        return execute(args)
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
