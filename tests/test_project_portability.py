"""Portable downstream routes and single-setting FASTQ relocation."""
from pathlib import Path
import shutil
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from ytab.pipeline import project_config
from ytab.pipeline.mapfastq_runner import load_project_for_mapping
from ytab.pipeline.create_hit_file_runner import find_sample_bam
from ytab.pipeline.summary_table_runner import find_sample_hits_file


@pytest.fixture
def initialized(tmp_path, monkeypatch):
    root = tmp_path/'old checkout';root.mkdir()
    fastqs = tmp_path/'old fastqs';fastqs.mkdir()
    for mate in (1, 2):
        (fastqs/f'parent_pool1_{mate}.fastq.gz').write_bytes(b'fixture')
    monkeypatch.setattr(project_config, 'resolve_reference',
        lambda *a: SimpleNamespace(to_dict=lambda: dict(fasta='resources/ref.fna')))
    project_config.create_project_config('demo', fastqs, 'glabrata', root, 2)
    cfg = root/'output/projects/demo/config/project.yaml'
    return root, fastqs, cfg


def test_new_configs_keep_only_fastq_base_machine_specific(initialized):
    root, fastqs, cfg = initialized
    stored = yaml.safe_load(cfg.read_text())
    assert stored['repo_root'] == '../../../..'
    assert stored['fastq_paths_relative_to'] == 'fastq_dir'
    sheet = pd.read_csv(cfg.with_name('sample_sheet.csv'), keep_default_na=False)
    assert sheet.iloc[0]['fastq_1'] == 'parent_pool1_1.fastq.gz'
    loaded = load_project_for_mapping(cfg)
    assert Path(loaded['samples'][0]['fastq_1']) == fastqs/'parent_pool1_1.fastq.gz'
    assert Path(loaded['_repo_root']) == root


def test_move_checkout_and_change_only_fastq_base(initialized, tmp_path, monkeypatch):
    root, fastqs, cfg = initialized
    moved = tmp_path/'new checkout';shutil.move(root, moved)
    moved_reads = tmp_path/'new reads';shutil.copytree(fastqs, moved_reads)
    cfg = moved/'output/projects/demo/config/project.yaml'
    stored = yaml.safe_load(cfg.read_text());stored['fastq_dir'] = str(moved_reads)
    cfg.write_text(yaml.safe_dump(stored))
    monkeypatch.chdir(tmp_path)
    loaded = load_project_for_mapping(cfg)
    sample = loaded['samples'][0]
    assert Path(sample['fastq_1']) == moved_reads/'parent_pool1_1.fastq.gz'
    assert find_sample_bam(sample, loaded) == moved/'output/projects/demo/mapfastq/parent_pool1/parent_pool1.sorted.bam'
    assert find_sample_hits_file(sample, loaded) == moved/'output/projects/demo/create_hit_file/parent_pool1/parent_pool1_hits.txt'
    assert Path(project_config.load_project_config(cfg)['_repo_root']) == moved


def test_legacy_rebase_preserves_metadata_and_downstream_settings(initialized, tmp_path):
    root, fastqs, cfg = initialized
    stored = yaml.safe_load(cfg.read_text())
    stored.pop('fastq_paths_relative_to')
    stored['repo_root'] = str(root)
    stored['mapping'] = dict(junction_mate='R2')
    for row in stored['samples']:
        for key in ('fastq_1', 'fastq_2'):
            row[key] = str(fastqs/row[key]) if row[key] else ''
    cfg.write_text(yaml.safe_dump(stored))
    sheet_path = cfg.with_name('sample_sheet.csv')
    sheet = pd.read_csv(sheet_path, dtype=str, keep_default_na=False)
    sheet['condition'] = 'parent';sheet['pool'] = '1';sheet['notes'] = 'reviewed metadata'
    for key in ('fastq_1', 'fastq_2'):
        sheet[key] = sheet[key].map(lambda p: str(fastqs/p) if p else '')
    sheet.to_csv(sheet_path, index=False)
    relocated = tmp_path/'hpc reads';relocated.mkdir()
    project_config.update_fastq_base(cfg, relocated)
    after = yaml.safe_load(cfg.read_text())
    assert after['mapping'] == stored['mapping']
    assert after['output_project_dir'] == stored['output_project_dir']
    assert after['reference'] == stored['reference']
    loaded = load_project_for_mapping(cfg)
    row = loaded['samples'][0]
    assert Path(row['fastq_1']) == relocated/'parent_pool1_1.fastq.gz'
    assert row['condition'] == 'parent' and row['notes'] == 'reviewed metadata'
    assert row['pool'] == '1'


def test_rebase_refuses_unrelated_paths_without_writes(initialized, tmp_path):
    root, fastqs, cfg = initialized
    sheet = cfg.with_name('sample_sheet.csv')
    frame = pd.read_csv(sheet, dtype=str, keep_default_na=False)
    frame.loc[0, 'fastq_1'] = '/unrelated/library.fastq.gz'
    frame.to_csv(sheet, index=False)
    before = cfg.read_bytes(), sheet.read_bytes()
    with pytest.raises(ValueError, match='outside'):
        project_config.update_fastq_base(cfg, fastqs)
    assert before == (cfg.read_bytes(), sheet.read_bytes())


def test_imported_hits_do_not_become_fastq_projects(initialized):
    root, fastqs, cfg = initialized
    stored = yaml.safe_load(cfg.read_text());stored['start_stage'] = 'create_hit_file'
    cfg.write_text(yaml.safe_dump(stored));before = cfg.read_bytes()
    with pytest.raises(ValueError, match='Imported-hit'):
        project_config.update_fastq_base(cfg, fastqs)
    assert cfg.read_bytes() == before
