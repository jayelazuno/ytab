#!/usr/bin/env python3
"""Draw/query gene-domain insertion figures from existing YTAB project outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from ytab.domain_explorer.domain_insertion_explorer import (  # noqa: E402
    build_gene_lookup,
    list_insertion_tracks,
    query_gene,
    run_gene_domain_explorer,
)


def _split_samples(value: str | None) -> list[str] | None:
    if value is None or value.strip() == "" or value.strip().lower() == "all":
        return None
    return [part.strip() for part in value.split(",") if part.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True, type=Path)
    parser.add_argument("--gene", help="Gene identifier/name to resolve and plot")
    parser.add_argument("--query", help="Query only; print candidate genes")
    parser.add_argument("--list-tracks", action="store_true", help="List available insertion tracks")
    parser.add_argument("--track-source", default="raw", choices=("raw", "combined", "combined_parent", "all"))
    parser.add_argument(
        "--track-preset",
        default="all",
        choices=("all", "parents", "treated", "matched_pairs", "pool1_pair", "pool2_pair", "pool3_pair", "pool4_pair", "custom"),
        help="Preset track selection; use custom with --samples for explicit tracks",
    )
    parser.add_argument("--samples", default="all", help="Comma-separated sample/track names, or all")
    parser.add_argument("--flank-bp", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--width-px", type=int, default=1800)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--label-mode", choices=("full", "compact"), default="full")
    parser.add_argument("--hide-site-counts", action="store_true")
    parser.add_argument("--hide-domains", action="store_true")
    parser.add_argument("--hide-direction", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    return parser.parse_args()


def _print_json(data: object) -> None:
    print(json.dumps(data, indent=2, sort_keys=True))


def main() -> int:
    args = parse_args()
    try:
        if args.list_tracks:
            tracks = list_insertion_tracks(args.project_config, args.track_source)
            if args.json:
                _print_json({"status": "success", "tracks": tracks})
            else:
                print(f"Available {args.track_source} insertion tracks: {len(tracks)}")
                for track in tracks:
                    print(f"  {track['sample']}\t{track['track_source']}\t{track['source_file']}")
            return 0

        if args.query:
            lookup = build_gene_lookup(args.project_config)
            candidates = query_gene(args.project_config, args.query)
            status = "no_match" if not candidates else ("resolved" if len(candidates) == 1 else "ambiguous")
            data = {
                "status": status,
                "query": args.query,
                "project_id": lookup["project_id"],
                "species": lookup["species"],
                "annotation_source": lookup["annotation_source"],
                "candidates": candidates,
            }
            if args.json:
                _print_json(data)
            else:
                print(f"Query: {args.query}")
                print(f"Status: {status}")
                if not candidates:
                    print("No matching gene found.")
                for row in candidates:
                    print(
                        "\t".join(
                            str(row.get(key, ""))
                            for key in ("gene_id", "display_name", "chromosome", "start", "end", "strand", "product")
                        )
                    )
            return 0

        if not args.gene:
            print("ERROR: provide --gene, --query, or --list-tracks", file=sys.stderr)
            return 2

        result = run_gene_domain_explorer(
            args.project_config,
            args.gene,
            samples=_split_samples(args.samples),
            track_source=args.track_source,
            track_preset=args.track_preset,
            flank_bp=args.flank_bp,
            show_domains=not args.hide_domains,
            show_direction=not args.hide_direction,
            width_px=args.width_px,
            dpi=args.dpi,
            label_mode=args.label_mode,
            show_site_counts=not args.hide_site_counts,
            output=args.output,
            force=args.force,
        )
    except Exception as exc:
        if args.json:
            _print_json({"status": "failure", "error": str(exc)})
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if result.get("status") in {"no_match", "ambiguous"}:
        if args.json:
            _print_json(result)
        else:
            print(f"Status: {result['status']}")
            if result["status"] == "no_match":
                print("No matching gene found.")
            else:
                print("Candidate genes:")
                for row in result.get("candidates", []):
                    print(f"  {row.get('gene_id')}  {row.get('display_name')}  {row.get('chromosome')}:{row.get('start')}-{row.get('end')}")
        return 0

    manifest = result["manifest"]
    if args.json:
        _print_json({"status": manifest.get("status", "success"), "gene": result["gene"], "manifest": manifest})
    else:
        print(f"Matched gene: {manifest['resolved_gene_id']} ({manifest['resolved_gene_name']})")
        print(f"Output PNG: {manifest['figure_path']}")
        print(f"Output insertion table: {manifest['table_path']}")
        print(f"Output manifest: {manifest['manifest_path']}")
        print(manifest.get("domain_message") or "")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
