
"""
coverage of my aligner against bowtie2, from the samtools depth files
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_custom_alignments(filepath):
    """Load custom aligner TSV results, return list of (pos, read_len)."""
    alignments = []
    with open(filepath) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split('\t')
            pos = int(parts[2])
            rlen = int(parts[5])
            if pos >= 0:
                alignments.append((pos, rlen))
    return alignments

def load_samtools_depth(filepath, genome_len):
    """Load samtools depth output into coverage array."""
    coverage = np.zeros(genome_len, dtype=np.int32)
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            pos = int(parts[1]) - 1  # samtools is 1-based
            depth = int(parts[2])
            if 0 <= pos < genome_len:
                coverage[pos] = depth
    return coverage

def compute_coverage(alignments, genome_len):
    coverage = np.zeros(genome_len, dtype=np.int32)
    for pos, rlen in alignments:
        end = min(pos + rlen, genome_len)
        coverage[pos:end] += 1
    return coverage

def plot_comparison(cov_custom, cov_bowtie, regions, output_file):
    """Plot  coverage comparison."""
    n = len(regions)
    fig, axes = plt.subplots(n, 2, figsize=(16, 4 * n), sharey='row')
    if n == 1:
        axes = axes.reshape(1, 2)

    for i, (name, start, end) in enumerate(regions):
        positions = np.arange(start, end)


        axes[i, 0].fill_between(positions, cov_custom[start:end], alpha=0.6, color='#2196F3')
        axes[i, 0].plot(positions, cov_custom[start:end], linewidth=0.5, color='#1565C0')
        axes[i, 0].set_title(f'{name} - Custom Aligner', fontsize=10)
        axes[i, 0].set_ylabel('Coverage')
        axes[i, 0].grid(True, alpha=0.3)
        axes[i, 0].set_xlim(start, end)

        # Bowtie2
        axes[i, 1].fill_between(positions, cov_bowtie[start:end], alpha=0.6, color='#4CAF50')
        axes[i, 1].plot(positions, cov_bowtie[start:end], linewidth=0.5, color='#2E7D32')
        axes[i, 1].set_title(f'{name} - Bowtie2', fontsize=10)
        axes[i, 1].grid(True, alpha=0.3)
        axes[i, 1].set_xlim(start, end)

    axes[-1, 0].set_xlabel('Genomic Position (chr19)')
    axes[-1, 1].set_xlabel('Genomic Position (chr19)')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved comparison plot: {output_file}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_results.py <custom_results.tsv> <bowtie2_depth.tsv> [genome_len]")
        sys.exit(1)

    custom_file = sys.argv[1]
    bowtie_depth_file = sys.argv[2]
    genome_len = int(sys.argv[3]) if len(sys.argv) > 3 else 59128983

    # Load custom aligner results
    custom_aligns = load_custom_alignments(custom_file)
    print(f"Custom aligner: {len(custom_aligns):,} mapped reads")
    cov_custom = compute_coverage(custom_aligns, genome_len)

    # Load bowtie2 depth
    print("Loading bowtie2 depth...")
    cov_bowtie = load_samtools_depth(bowtie_depth_file, genome_len)


    regions = [
        ("Region 1 - chr19:27,731,500-27,735,500", 27731500, 27735500),
        ("Region 2 - chr19:42,770,500-42,774,500", 42770500, 42774500),
        ("Region 3 - chr19:39,338,500-39,342,500", 39338500, 39342500),
        ("Region 4 - chr19:42,745,000-42,749,000", 42745000, 42749000),
        ("Region 5 - chr19:54,692,500-54,696,500", 54692500, 54696500),
    ]

    # Try loading all the regions from coverage_summary.tsv
    try:
        with open("coverage_summary.tsv") as f:
            f.readline()
            regions = []
            for line in f:
                parts = line.strip().split('\t')
                rank = int(parts[0])
                start = int(parts[1])
                end = int(parts[2])
                pad = 500
                regions.append((f"Enriched Region {rank}", max(0, start - pad), min(genome_len, end + pad)))
                if len(regions) >= 5:
                    break
        print(f"Using {len(regions)} regions from coverage_summary.tsv")
    except FileNotFoundError:
        print("Using default example regions")

    plot_comparison(cov_custom, cov_bowtie, regions, "comparison_coverage.png")

    # Print summary comparison
    total_custom = cov_custom.sum()
    total_bowtie = cov_bowtie.sum()
    covered_custom = (cov_custom > 0).sum()
    covered_bowtie = (cov_bowtie > 0).sum()
    
    print(f"\n{'='*50}")
    print(f"COMPARISON SUMMARY")
    print(f"{'='*50}")
    print(f"{'Metric':<30} {'Custom':>12} {'Bowtie2':>12}")
    print(f"{'-'*54}")
    print(f"{'Total coverage (bp×reads)':<30} {total_custom:>12,} {total_bowtie:>12,}")
    print(f"{'Bases covered (≥1 read)':<30} {covered_custom:>12,} {covered_bowtie:>12,}")
    print(f"{'Max depth':<30} {cov_custom.max():>12,} {cov_bowtie.max():>12,}")
    print(f"{'Mean depth (covered bases)':<30} {cov_custom[cov_custom>0].mean():>12.1f} {cov_bowtie[cov_bowtie>0].mean():>12.1f}")

if __name__ == '__main__':
    main()
