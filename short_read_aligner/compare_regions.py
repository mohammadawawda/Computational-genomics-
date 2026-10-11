import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

regions = [
    ("Region 1", 27731500, 27735500),
    ("Region 2", 42770500, 42774500),
    ("Region 3", 39338500, 39342500),
    ("Region 4", 42745000, 42749000),
    ("Region 5", 54692500, 54696500),
]

def load_custom_region(filepath, start, end):
    cov = np.zeros(end - start, dtype=np.int32)
    with open(filepath) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            pos = int(parts[2])
            rlen = int(parts[5])
            if pos < 0:
                continue
            read_end = pos + rlen
            if read_end < start or pos >= end:
                continue
            s = max(pos, start) - start
            e = min(read_end, end) - start
            cov[s:e] += 1
    return cov

def load_depth_file(filepath, start, end):
    cov = np.zeros(end - start, dtype=np.int32)
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            pos = int(parts[1]) - 1
            depth = int(parts[2])
            if start <= pos < end:
                cov[pos - start] = depth
    return cov

custom_file = sys.argv[1]

for i, (name, start, end) in enumerate(regions):
    depth_file = f"depth_region{i+1}.tsv"
    
    print(f"Processing {name} (chr19:{start}-{end})...")
    cov_custom = load_custom_region(custom_file, start, end)
    cov_bowtie = load_depth_file(depth_file, start, end)
    
    positions = np.arange(start, end)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 4), sharey=True)
    
    axes[0].fill_between(positions, cov_custom, alpha=0.6, color='#2196F3')
    axes[0].plot(positions, cov_custom, linewidth=0.5, color='#1565C0')
    axes[0].set_title(f'{name} - Custom Aligner')
    axes[0].set_ylabel('Coverage')
    axes[0].set_xlabel('Position (chr19)')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].fill_between(positions, cov_bowtie, alpha=0.6, color='#4CAF50')
    axes[1].plot(positions, cov_bowtie, linewidth=0.5, color='#2E7D32')
    axes[1].set_title(f'{name} - Bowtie2')
    axes[1].set_xlabel('Position (chr19)')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"comparison_region_{i+1}.png", dpi=150)
    plt.close()
    print(f"  Saved comparison_region_{i+1}.png")
    
    print(f"  Custom: max={cov_custom.max()}, mean={cov_custom[cov_custom>0].mean():.1f}")
    print(f"  Bowtie: max={cov_bowtie.max()}, mean={cov_bowtie[cov_bowtie>0].mean():.1f}")
