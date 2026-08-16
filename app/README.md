# Local Shiny application

Launch from any working directory with:

```bash
./scripts/local/ytab_launch_app.sh --project-config output/projects/<PROJECT_ID>/config/project.yaml
```

The launcher checks the environment, selects another port if the requested one is occupied, and
binds locally to `127.0.0.1` unless another host is explicitly supplied. Pipeline profile runs use a
single background process with live output and cancellation; partial outputs remain available for
cache-aware resumption.
# Local Shiny application

Launch with `./scripts/local/ytab_launch_app.sh`. The app begins on the YTAB landing page, where users can initialize a project from a local FASTQ directory or open any project discovered under `output/projects/*/config/project.yaml`. A launcher-provided `--project-config` is preselected but never starts a stage automatically.

Projects without a completed raw SummaryTable enter Preprocessing; analysis-ready projects enter Overview. FASTQ data is referenced in place rather than uploaded. On low-memory systems, use two threads and consider completing the core profile in the terminal before launching the app.

## Gene & Domain Insertion Explorer

The app includes a **Gene & Domain Insertion Explorer** tab for analysis-ready
projects with existing hit files. Search by gene name, systematic name, locus
tag, gene ID, or any alias available in the project reference annotation. After
resolving a gene, choose raw CreateHitFile tracks, a combined parent hit file
when available, or all available tracks to plot insertion positions across the
gene plus the selected flank.

Domain annotations are drawn only when real domain data are available in the
YTAB reference resources. Genes without domain annotation still plot
successfully and show a no-domain message. Generated PNG figures, insertion
tables, and manifests are saved under:

```text
output/projects/<PROJECT_ID>/gene_domain_explorer/
```
