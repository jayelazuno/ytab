"""Standalone exports of the existing analytical Poisson benchmark."""
import argparse
import csv
import json
import sys
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Analytical equal-gene Poisson benchmark: expected recovery, not confidence that every gene is recovered.")
    parser.add_argument("--genes", type=int, required=True, help="Number of genes in the selected benchmark universe N")
    parser.add_argument("--target", type=float, default=95, help="Expected recovery percentage, strictly between 0 and 100 (default: 95)")
    parser.add_argument("--output", type=Path, default=Path("screen_coverage"))
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--json", action="store_true", help="Print machine-readable result summary")
    args = parser.parse_args(argv)
    if args.genes < 1 or not 0 < args.target < 100:
        parser.error("--genes must be positive and --target must be strictly between 0 and 100")
    from ytab.screen_coverage.statistics import colonies_for_target_coverage, poisson_coverage
    target = args.target / 100
    k = colonies_for_target_coverage(args.genes, target)
    coverage = poisson_coverage(k, args.genes)
    result = dict(model="Theoretical Poisson / Equal-Gene Benchmark", genes=args.genes,
                  target_percent=args.target, required_picks=k, expected_recovery_percent=100*coverage,
                  expected_genes=args.genes*coverage, pick_units="idealized equal-probability gene-hit opportunities",
                  interpretation="Expected recovery under independent equal-gene sampling; not an experimental recommendation or all-genes confidence.")
    output = args.output.resolve()
    try:
        output.mkdir(parents=True, exist_ok=True)
        metadata = output/"coverage.json"
        artifacts = [metadata, output/"summary.csv", output/"curve.csv"]
        if not args.no_plot:
            artifacts.append(output/"coverage.png")
        if metadata.exists():
            same = json.loads(metadata.read_text()) == result
            if not same and not args.force:
                raise ValueError("Different benchmark already exists at --output; choose another directory or use --force")
            reused = same and all(p.is_file() for p in artifacts)
        else:
            if any(p.exists() for p in artifacts) and not args.force:
                raise ValueError("Existing exports without matching metadata; choose another directory or use --force")
            reused = False
        if not reused or args.force:
            # Preserve the scientific implementation: all values use the shared accepted helpers.
            points = sorted({0, k, *(round(k*1.25*i/200) for i in range(201))})
            values = [100*poisson_coverage(x, args.genes) for x in points]
            with (output/"summary.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(result))
                writer.writeheader(); writer.writerow(result)
            with (output/"curve.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["idealized_gene_hit_picks", "expected_recovery_percent"])
                writer.writerows(zip(points, values))
            if not args.no_plot:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
                ax.plot(points, values, color="#F8766D")
                ax.axhline(args.target, color="#D55E00", linestyle="--")
                ax.axvline(k, color="#0072B2", linestyle=":")
                ax.scatter([k], [100*coverage], color="#0072B2")
                ax.set(xlabel="Idealized gene-hit picks K", ylabel="Expected genes recovered (%)",
                       title="Theoretical Poisson / Equal-Gene Benchmark", ylim=(0, 102))
                fig.savefig(output/"coverage.png", dpi=150)
                plt.close(fig)
            metadata.write_text(json.dumps(result, indent=2)+"\n")
        summary = dict(**result, output_dir=str(output), reused=reused and not args.force)
        print(json.dumps(summary) if args.json else f"N={args.genes:,}; expected recovery target={args.target:g}%; K={k:,} idealized gene-hit picks\nOutputs: {output}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
