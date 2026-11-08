#!/usr/bin/env python3
"""Build 1:1 ortholog groups (C. glabrata, S. pombe, C. albicans, each linked through S. cerevisiae) for cross-species library comparison.

Sources (inputs, not modified):
- C. glabrata - S. cerevisiae: 1:1 orthogroups from resources/species/glabrata/essentiality_labels/labels.tsv
  (orthology == "1:1"; built from resources/species/glabrata/orthologs/C_glabrata_BG2_S_cerevisiae_orthologs.txt);
- S. pombe - S. cerevisiae: PomBase curated orthologs (release 2025-10-01,
  resources/species/pombe/orthologs/pombe_cerevisiae_orthologs.txt). A pombe gene is 1:1 when it has exactly
  one S. cerevisiae gene with no fusion qualifier (no '|', '+', '(N)', '(C)', 'NONE'), and that S. cerevisiae
  gene is listed for no other pombe gene;
- C. albicans - S. cerevisiae: CGD orthologs (resources/species/albicans/orthologs/C_albicans_SC5314_S_cerevisiae_orthologs.txt,
  genome A22, created 2017-01-18); haplotype-A gene IDs (C*_A, as in the A22 GFF). 1:1 when each ID occurs once.
Outputs (group, species, gene; each S. cerevisiae gene is a group, listing every species with a 1:1 gene):
- resources/species/orthologs/glabrata_cerevisiae_pombe_1to1.tsv (three species; kept unchanged for earlier runs);
- resources/species/orthologs/glabrata_cerevisiae_pombe_albicans_1to1.tsv (four species).
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LABELS = ROOT/"resources/species/glabrata/essentiality_labels/labels.tsv"
POMBE = ROOT/"resources/species/pombe/orthologs/pombe_cerevisiae_orthologs.txt"
ALBICANS = ROOT/"resources/species/albicans/orthologs/C_albicans_SC5314_S_cerevisiae_orthologs.txt"
OUT = ROOT/"resources/species/orthologs"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    glabrata = {}
    with LABELS.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["orthology"] == "1:1" and r["scer_id"]:glabrata[r["scer_id"]] = r["cglab_id"]
    pombe, seen = {}, {}
    for line in POMBE.read_text().splitlines():
        if not line.strip() or line.startswith("#"):continue
        sp, sc = line.split("\t")[:2]
        if any(t in sc for t in ("|", "+", "(", "NONE")) or not sc.strip():continue
        seen[sc] = seen.get(sc, 0)+1; pombe[sc] = sp
    pombe = {sc:sp for sc, sp in pombe.items() if seen[sc] == 1}
    pairs = [l.split("\t") for l in ALBICANS.read_text().splitlines() if l.strip() and not l.startswith("#")]
    pairs = [(p[0], p[3]) for p in pairs if len(p) > 3 and p[0] and p[3]]
    count_ca, count_sc = {}, {}
    for ca, sc in pairs:count_ca[ca] = count_ca.get(ca, 0)+1; count_sc[sc] = count_sc.get(sc, 0)+1
    albicans = {sc:ca for ca, sc in pairs if count_ca[ca] == 1 and count_sc[sc] == 1}
    OUT.mkdir(parents=True, exist_ok=True)
    def write(name, species):
        rows = []
        for sc in sorted(set().union(*species.values())):
            rows.append(dict(group=sc, species="cerevisiae", gene=sc))
            rows += [dict(group=sc, species=sp, gene=table[sc]) for sp, table in species.items() if sc in table]
        path = OUT/name
        with path.open("w", newline="") as handle:
            w = csv.DictWriter(handle, ["group", "species", "gene"], delimiter="\t"); w.writeheader(); w.writerows(rows)
        return path
    sources = {str(LABELS.relative_to(ROOT)):sha(LABELS), str(POMBE.relative_to(ROOT)):sha(POMBE), str(ALBICANS.relative_to(ROOT)):sha(ALBICANS)}
    notes = dict(pombe_source="PomBase curated orthologs, release 2025-10-01 (https://www.pombase.org/monthly_releases/2025/pombase-2025-10-01/curated_orthologs/)",
                 albicans_source="CGD C_albicans_SC5314_S_cerevisiae_orthologs.txt (genome A22-s07-m01-r18, created 2017-01-18); haplotype A gene IDs",
                 script=str(Path(__file__).relative_to(ROOT)))
    three = write("glabrata_cerevisiae_pombe_1to1.tsv", dict(glabrata=glabrata, pombe=pombe))
    summary3 = dict(glabrata_cerevisiae_1to1=len(glabrata), pombe_cerevisiae_1to1=len(pombe), all_three=len(set(glabrata) & set(pombe)),
                    sources={k:v for k, v in sources.items() if "albicans" not in k}, output_sha256=sha(three)) | {k:v for k, v in notes.items() if k != "albicans_source"}
    (OUT/"glabrata_cerevisiae_pombe_1to1.json").write_text(json.dumps(summary3, indent=2)+"\n")
    four = write("glabrata_cerevisiae_pombe_albicans_1to1.tsv", dict(glabrata=glabrata, pombe=pombe, albicans=albicans))
    g, p, a = set(glabrata), set(pombe), set(albicans)
    summary4 = dict(glabrata_cerevisiae_1to1=len(g), pombe_cerevisiae_1to1=len(p), albicans_cerevisiae_1to1=len(a),
                    glabrata_albicans=len(g & a), glabrata_pombe=len(g & p), pombe_albicans=len(p & a), all_four=len(g & p & a),
                    sources=sources, output_sha256=sha(four)) | notes
    (OUT/"glabrata_cerevisiae_pombe_albicans_1to1.json").write_text(json.dumps(summary4, indent=2)+"\n"); print(json.dumps(summary4, indent=1))


if __name__ == "__main__":main()
