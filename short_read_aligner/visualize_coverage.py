
"""
coverage of my alignments and plots of the top enriched regions
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

def load_alignments(filepath):
    """Load alignment results TSV file."""
    alignments = []
    with open(filepath) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            name, fl, pos, strand, mm, rlen = parts
            pos = int(pos)
            rlen = int(rlen)
            mm = int(mm)
            if pos >= 0:  # mapped
                alignments.append((name, pos, strand, rlen))
    return alignments

def compute_coverage(alignments, genome_len):
    """Compute per-base coverage array."""
    print("Computing coverage...")
    coverage = np.zeros(genome_len, dtype=np.int32)
    for name, pos, strand, rlen in alignments:
        end = min(pos + rlen, genome_len)
        coverage[pos:end] += 1
    return coverage

def find_enriched_regions(coverage, window_size=2000, top_n=20):
    """Find top enriched windows by summing coverage."""
    print(f"Finding enriched regions (window={window_size})...")
    best = []
    step = window_size // 2  # 50% overlap
    for start in range(0, len(coverage) - window_size, step):
        end = start + window_size
        total = coverage[start:end].sum()
        best.append((total, start, end))
    best.sort(reverse=True)
    
    # Merge overlapping top regions
    selected = []
    for total, start, end in best:
        overlap = False
        for _, s2, e2 in selected:
            if start < e2 and end > s2:
                overlap = True
                break
        if not overlap:
            selected.append((total, start, end))
        if len(selected) >= top_n:
            break
    return selected

def plot_region(coverage, start, end, region_name, output_file):
    """Plot coverage for a specific genomic region."""
    fig, ax = plt.subplots(figsize=(12, 4))
    positions = np.arange(start, end)
    cov = coverage[start:end]
    
    ax.fill_between(positions, cov, alpha=0.6, color='#2196F3')
    ax.plot(positions, cov, linewidth=0.5, color='#1565C0')
    ax.set_xlabel('Genomic Position (chr19)', fontsize=11)
    ax.set_ylabel('Coverage (reads)', fontsize=11)
    ax.set_title(f'{region_name}\nchr19:{start:,}-{end:,} ({end-start:,} bp)', fontsize=12)
    ax.set_xlim(start, end)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_coverage.py <alignment_results.tsv> [genome_length]")
        sys.exit(1)

    alignment_file = sys.argv[1]
    genome_len = int(sys.argv[2]) if len(sys.argv) > 2 else 59128983  # chr19 hg19 length

    alignments = load_alignments(alignment_file)
    print(f"Loaded {len(alignments):,} mapped reads")

    coverage = compute_coverage(alignments, genome_len)

    # Find top enriched regions
    enriched = find_enriched_regions(coverage, window_size=3000, top_n=10)

    print(f"\nTop enriched regions:")
    for i, (total, start, end) in enumerate(enriched[:5]):
        max_cov = coverage[start:end].max()
        print(f"  {i+1}. chr19:{start:,}-{end:,}  total_cov={total:,}  max_depth={max_cov}")

    # Plot top 5 regions
    for i, (total, start, end) in enumerate(enriched[:5]):
        # Expand slightly for context
        pad = 500
        plot_start = max(0, start - pad)
        plot_end = min(genome_len, end + pad)
        plot_region(coverage, plot_start, plot_end,
                   f"Region {i+1} (Total coverage: {total:,})",
                   f"coverage_region_{i+1}.png")

    # Save coverage summary
    with open("coverage_summary.tsv", 'w') as f:
        f.write("rank\tstart\tend\tsize\ttotal_coverage\tmax_depth\n")
        for i, (total, start, end) in enumerate(enriched[:10]):
            max_cov = coverage[start:end].max()
            f.write(f"{i+1}\t{start}\t{end}\t{end-start}\t{total}\t{max_cov}\n")

    print("\nDone. Coverage plots saved.")

if __name__ == '__main__':
    main()
