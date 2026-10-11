# Short Read Alignment

My own seed based aligner for ATAC-seq short reads on human chromosome 19 (hg19), and a comparison of its results with bowtie2.

## Files

| File | What it does |
|------|-------------|
| `aligner.py` | my custom seed-and-extend aligner |
| `visualize_coverage.py` | computes coverage and plots the top enriched regions |
| `compare_regions.py` | side-by-side coverage comparison between my aligner and bowtie2 |
| `compare_results.py` | coverage plots of my aligner vs bowtie2 from the samtools depth files |
| `report.pdf` | The full report with explanations, pseudocode, plots, and analysis |

## Data

The data files (not included):
- `hg19-chr19.fa.gz` - reference genome (chr19, hg19)
- `ATAC.chr19.R01.fastq.gz` - read 1 (40bp, 677,996 reads)
- `ATAC.chr19.R02.fastq.gz` - read 2 (36bp, 677,996 reads)

## How to run

### Custom aligner

```
pip install numpy matplotlib
python aligner.py hg19-chr19.fa.gz ATAC.chr19.R01.fastq.gz ATAC.chr19.R02.fastq.gz alignment_results.tsv
```

This takes about 10-30 minutes and needs ~4-6 GB of RAM. It outputs `alignment_results.tsv` with the mapping results.

To measure time and memory on Linux:
```
/usr/bin/time -v python aligner.py hg19-chr19.fa.gz ATAC.chr19.R01.fastq.gz ATAC.chr19.R02.fastq.gz
```

### Coverage plots

```
python visualize_coverage.py alignment_results.tsv 59128983
```

This produces `coverage_region_1.png` through `coverage_region_5.png` and `coverage_summary.tsv`.

### Bowtie2 and the comparison

You need bowtie2 and samtools installed. On Ubuntu/WSL:
```
sudo apt install bowtie2 samtools
```

Then decompress the genome and run the pipeline:
```
gzip -dk hg19-chr19.fa.gz
bowtie2-build hg19-chr19.fa chr19_bt2_index
bowtie2 -x chr19_bt2_index -1 ATAC.chr19.R01.fastq.gz -2 ATAC.chr19.R02.fastq.gz --maxins 2000 --no-mixed --no-discordant --threads 4 -S bowtie2_output.sam
samtools view -bS bowtie2_output.sam -o bowtie2_output.bam
samtools sort bowtie2_output.bam -o bowtie2_output.sorted.bam
samtools index bowtie2_output.sorted.bam
samtools flagstat bowtie2_output.sorted.bam
```

Extract depth for the 5 enriched regions:
```
samtools depth -a -r chr19:27731500-27735500 bowtie2_output.sorted.bam > depth_region1.tsv
samtools depth -a -r chr19:42770500-42774500 bowtie2_output.sorted.bam > depth_region2.tsv
samtools depth -a -r chr19:39338500-39342500 bowtie2_output.sorted.bam > depth_region3.tsv
samtools depth -a -r chr19:42745000-42749000 bowtie2_output.sorted.bam > depth_region4.tsv
samtools depth -a -r chr19:54692500-54696500 bowtie2_output.sorted.bam > depth_region5.tsv
```

Then run the comparison:
```
python compare_regions.py alignment_results.tsv
```

This produces `comparison_region_1.png` through `comparison_region_5.png`.

## Results

- Custom aligner mapping rate: **98.10%** (R01: 98.58%, R02: 97.61%)
- Both aligners produce very similar coverage profiles across all 5 enriched regions
- The peak shapes and positions match closely, confirming the seed-and-extend approach works well

## Requirements

- Python 3.x
- numpy
- matplotlib
- bowtie2 (for the comparison, Linux/WSL only)
- samtools (for the comparison, Linux/WSL only)
