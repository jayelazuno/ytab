"""Test the personal SGE wrapper without submitting jobs or running real analyses."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
JOB=ROOT/'pipeline/hpc/ytab_mapFastq.txt'
HIT_JOB=ROOT/'pipeline/hpc/ytab_CreateHitFile.txt'


@pytest.fixture
def environment(tmp_path):
    repo=tmp_path/'checkout with spaces'
    (repo/'scripts/local').mkdir(parents=True)
    (repo/'scripts/local/ytab.py').write_text(
        'import json,os,sys\nfrom pathlib import Path\n'
        'Path(os.environ["CAPTURE"]).write_text(json.dumps(dict(args=sys.argv[1:],cwd=os.getcwd(),tmpdir=os.environ["TMPDIR"])))\n'
        'raise SystemExit(int(os.environ.get("FAKE_EXIT","0")))\n')
    project=repo/'output/projects/example/config/project.yaml'
    project.parent.mkdir(parents=True);project.write_text('project_id: example\n')
    conda=tmp_path/'conda.sh';conda.write_text('conda() { return 0; }\n')
    bin_dir=tmp_path/'bin';bin_dir.mkdir()
    for name in ['bowtie2','samtools']:
        tool=bin_dir/name;tool.write_text('#!/bin/bash\necho "test executable"\n');tool.chmod(0o755)
    env=dict(os.environ,YTAB=str(repo),YTAB_CONDA_SH=str(conda),YTAB_CONDA_ENV='test',
             YTAB_ENV_MODE='conda',PYTHON_BIN=sys.executable,PROJECT_ID='example',PROJECT_CONFIG='',
             NSLOTS='8',SCRATCH_ROOT=str(tmp_path/'scratch'),CAPTURE=str(tmp_path/'capture.json'),
             PATH=str(bin_dir)+os.pathsep+os.environ['PATH'],JUNCTION_MATE='',
             REQUIRE_MATE_ALIGNMENT='',FORCE='0',DRY_RUN='0',KEEP_GOING='1',SAMPLES='')
    return env,tmp_path


def execute(env,cwd,job=JOB):
    return subprocess.run(['bash',str(job)],cwd=cwd,env=env,capture_output=True,text=True)


def test_bash_syntax_and_sge_header():
    assert JOB.read_text().startswith('#!/bin/bash\n')
    subprocess.run(['bash','-n',str(JOB)],check=True)
    assert '#$ -pe smp 8' in JOB.read_text()


def test_same_cli_routes_and_allocated_threads(environment):
    env,cwd=environment;result=execute(env,cwd)
    assert result.returncode==0,result.stdout+result.stderr
    captured=json.loads(Path(env['CAPTURE']).read_text());args=captured['args']
    assert args[0]=='mapfastq'
    assert args[args.index('--threads')+1]=='8'
    assert args[args.index('--project-config')+1]==str(Path(env['YTAB'])/'output/projects/example/config/project.yaml')
    assert '--junction-mate' not in args and '--require-mate-alignment' not in args
    assert '--output' not in args and '--keep-going' in args
    assert captured['cwd']==env['YTAB']
    assert Path(captured['tmpdir']).parent==Path(env['SCRATCH_ROOT'])
    assert not (Path(env['YTAB'])/'output/projects/example/mapfastq').exists()


def test_explicit_settings_and_config_path(environment):
    env,cwd=environment
    cfg=cwd/'other project.yaml';cfg.write_text('project_id: other\n')
    env.update(PROJECT_ID='',PROJECT_CONFIG=cfg.name,JUNCTION_MATE='R2',
               REQUIRE_MATE_ALIGNMENT='0',FORCE='1',DRY_RUN='1',SAMPLES='a,b',KEEP_GOING='0')
    result=execute(env,cwd)
    assert result.returncode==0,result.stdout+result.stderr
    args=json.loads(Path(env['CAPTURE']).read_text())['args']
    assert args[args.index('--project-config')+1]==str(cfg)
    assert args[args.index('--junction-mate')+1]=='R2'
    assert args[args.index('--samples')+1]=='a,b'
    assert '--no-require-mate-alignment' in args
    assert '--force' in args and '--dry-run' in args and '--keep-going' not in args


def test_explicit_concordance_and_failure_propagation(environment):
    env,cwd=environment;env.update(REQUIRE_MATE_ALIGNMENT='1',FAKE_EXIT='7')
    result=execute(env,cwd)
    assert result.returncode==7
    assert '--require-mate-alignment' in json.loads(Path(env['CAPTURE']).read_text())['args']
    assert 'SUCCESS:' not in result.stdout


@pytest.mark.parametrize('overrides',[
    {'PROJECT_ID':''},{'PROJECT_ID':'../bad'}, {'NSLOTS':'1'},
    {'NSLOTS':'invalid'}, {'JUNCTION_MATE':'both'}, {'FORCE':'yes'},
    {'REQUIRE_MATE_ALIGNMENT':'maybe'}, {'YTAB_ENV_MODE':'unknown'},
])
def test_bad_configuration_stops_before_cli(environment,overrides):
    env,cwd=environment;env.update(overrides)
    result=execute(env,cwd)
    assert result.returncode==2,result.stdout+result.stderr
    assert not Path(env['CAPTURE']).exists()


def test_module_mode_uses_declared_modules(environment):
    env,cwd=environment
    init=cwd/'module-init.sh';init.write_text('module() { echo "test module: $*"; }\n')
    env.update(YTAB_ENV_MODE='modules',YTAB_MODULES='python/test bowtie2/test samtools/test',BASH_ENV=str(init))
    result=execute(env,cwd)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'test module: load python/test bowtie2/test samtools/test' in result.stdout
    assert Path(env['CAPTURE']).is_file()


def test_hit_job_header_and_syntax():
    assert HIT_JOB.read_text().startswith('#!/bin/bash\n')
    assert '#$ -N YTAB_CreateHitFile' in HIT_JOB.read_text()
    subprocess.run(['bash','-n',str(HIT_JOB)],check=True)


def test_hit_job_uses_shared_cli_and_project_routes(environment):
    env,cwd=environment
    result=execute(env,cwd,HIT_JOB)
    assert result.returncode==0,result.stdout+result.stderr
    captured=json.loads(Path(env['CAPTURE']).read_text())
    assert captured['args']==['hit-files','--project-config',
        str(Path(env['YTAB'])/'output/projects/example/config/project.yaml'),
        '--threads','8','--keep-going']
    assert captured['cwd']==env['YTAB']
    assert Path(captured['tmpdir']).parent==Path(env['SCRATCH_ROOT'])
    assert not (Path(env['YTAB'])/'output/projects/example/create_hit_file').exists()


def test_hit_job_explicit_settings_and_single_slot(environment):
    env,cwd=environment
    cfg=cwd/'other project.yaml';cfg.write_text('project_id: other\n')
    env.update(PROJECT_ID='',PROJECT_CONFIG=cfg.name,NSLOTS='1',FORCE='1',
               DRY_RUN='1',SAMPLES='a,b',KEEP_GOING='0')
    result=execute(env,cwd,HIT_JOB)
    assert result.returncode==0,result.stdout+result.stderr
    assert json.loads(Path(env['CAPTURE']).read_text())['args']==[
        'hit-files','--project-config',str(cfg),'--threads','1',
        '--samples','a,b','--force','--dry-run']


def test_hit_job_failure_propagates(environment):
    env,cwd=environment;env['FAKE_EXIT']='7'
    result=execute(env,cwd,HIT_JOB)
    assert result.returncode==7
    assert 'SUCCESS:' not in result.stdout


@pytest.mark.parametrize('overrides',[
    {'PROJECT_ID':''},{'PROJECT_ID':'../bad'},{'NSLOTS':'0'},
    {'NSLOTS':'invalid'},{'FORCE':'yes'},{'DRY_RUN':'maybe'},
    {'KEEP_GOING':'yes'},{'YTAB_ENV_MODE':'unknown'},
])
def test_hit_job_bad_configuration(environment,overrides):
    env,cwd=environment;env.update(overrides)
    result=execute(env,cwd,HIT_JOB)
    assert result.returncode==2,result.stdout+result.stderr
    assert not Path(env['CAPTURE']).exists()


def test_hit_job_module_mode(environment):
    env,cwd=environment
    init=cwd/'module-init.sh';init.write_text('module() { echo "test module: $*"; }\n')
    env.update(YTAB_ENV_MODE='modules',YTAB_MODULES='python/test',BASH_ENV=str(init))
    result=execute(env,cwd,HIT_JOB)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'test module: load python/test' in result.stdout
