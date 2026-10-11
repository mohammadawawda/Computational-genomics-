
"""
seed and extend aligner
"""

import gzip
import sys
import time
import os
from collections import defaultdict


K = 14                  # kmer length
MAX_MISMATCHES = 2      # max mismatches allowed in full alignment
MAX_SEED_HITS = 500     # skip repetitive seeds
NUM_SEEDS = 3           # number of non overlapping seeds per read



COMP = str.maketrans('ACGTN', 'TGCAN')

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMP)[::-1]

def load_fasta(filepath):
    """Load a FASTA file  and return the sequence in uppercase"""
    print(f"Loading genome from {filepath}...")
    opener = gzip.open if filepath.endswith('.gz') else open
    seq_parts = []
    with opener(filepath, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq_parts.append(line.strip().upper())
    genome = ''.join(seq_parts)
    print(f"  Genome length: {len(genome):,} bp")
    return genome

def load_fastq(filepath):
    """Load a FASTQ file and return list of (name, sequence, quality)"""
    print(f"Loading reads from {filepath}...")
    opener = gzip.open if filepath.endswith('.gz') else open
    reads = []
    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip().upper()
            f.readline()  # + line
            qual = f.readline().strip()
            name = header[1:].split()[0]  # remove @ and keep first token
            reads.append((name, seq, qual))
    print(f"  Loaded {len(reads):,} reads")
    return reads


def build_kmer_index(genome, k):
    """
    Build a hash table mapping each k-mer in the genome to a list of positions while skipping k-mers containing 'N'
    """
    print(f"Building k-mer index (k={k})...")
    t0 = time.time()
    index = defaultdict(list)
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        if 'N' in kmer:
            continue
        index[kmer].append(i)
    elapsed = time.time() - t0
    print(f"  Index built in {elapsed:.1f}s, {len(index):,} unique k-mers")
    return index



def count_mismatches(read, genome, pos, max_mm):
    """Count mismatches between read and genome at given position"""
    mm = 0
    for i in range(len(read)):
        if read[i] != genome[pos + i]:
            mm += 1
            if mm > max_mm:
                return mm
    return mm

def get_seed_offsets(read_len, k, num_seeds):
    """Compute seed starting positions by dividing the read into equal segments so no two seeds overlap."""
    if num_seeds == 1:
        return [0]
    offsets = []
    step = (read_len - k) // (num_seeds - 1)
    for i in range(num_seeds):
        offset = min(i * step, read_len - k)
        offsets.append(offset)
    return sorted(set(offsets))

def align_read(read_seq, genome, index, k, max_mm=MAX_MISMATCHES):
    """
    Try to find the best matching position for a read on the genome by looking up
     seeds and extending. Returns the position, strand, and mismatch count or None if no valid match is found.
    """
    genome_len = len(genome)
    read_len = len(read_seq)
    best_pos = None
    best_mm = max_mm + 1
    best_strand = None

    seed_offsets = get_seed_offsets(read_len, k, NUM_SEEDS)

    for strand in ('+', '-'):
        query = read_seq if strand == '+' else reverse_complement(read_seq)

        for offset in seed_offsets:
            seed = query[offset:offset + k]
            if seed not in index:
                continue
            hits = index[seed]
            # Skip repetitive seeds
            if len(hits) > MAX_SEED_HITS:
                continue

            for genome_pos in hits:
                read_start = genome_pos - offset
                if read_start < 0 or read_start + read_len > genome_len:
                    continue
                mm = count_mismatches(query, genome, read_start, best_mm - 1)
                if mm < best_mm:
                    best_mm = mm
                    best_pos = read_start
                    best_strand = strand
                    if mm == 0:
                        return (best_pos, best_strand, best_mm)

    if best_pos is not None and best_mm <= max_mm:
        return (best_pos, best_strand, best_mm)
    return None

# ============== MAIN ==============

def main():
    if len(sys.argv) < 4:
        print("Usage: python aligner.py <genome.fa.gz> <reads_R01.fastq.gz> <reads_R02.fastq.gz> [output.tsv]")
        sys.exit(1)

    genome_file = sys.argv[1]
    r01_file = sys.argv[2]
    r02_file = sys.argv[3]
    output_file = sys.argv[4] if len(sys.argv) > 4 else "alignment_results.tsv"

    total_start = time.time()

    # Load genome
    genome = load_fasta(genome_file)

    # Build index
    index = build_kmer_index(genome, K)

    # Load reads
    reads_r01 = load_fastq(r01_file)
    reads_r02 = load_fastq(r02_file)

    # Align all reads
    print("\nAligning reads...")
    align_start = time.time()
    mapped = 0
    total = 0
    results = []

    for file_label, reads in [("R01", reads_r01), ("R02", reads_r02)]:
        file_mapped = 0
        for i, (name, seq, qual) in enumerate(reads):
            total += 1
            result = align_read(seq, genome, index, K)
            if result is not None:
                pos, strand, mm = result
                mapped += 1
                file_mapped += 1
                results.append((name, file_label, pos, strand, mm, len(seq)))
            else:
                results.append((name, file_label, -1, '.', -1, len(seq)))

            if (i + 1) % 50000 == 0:
                print(f"  [{file_label}] Processed {i+1:,} / {len(reads):,} reads...")

        print(f"  [{file_label}] Mapped: {file_mapped:,} / {len(reads):,} ({100*file_mapped/len(reads):.2f}%)")

    align_time = time.time() - align_start
    total_time = time.time() - total_start

    # Summary
    print(f"\n{'='*50}")
    print(f"ALIGNMENT SUMMARY")
    print(f"{'='*50}")
    print(f"Total reads:     {total:,}")
    print(f"Mapped reads:    {mapped:,}")
    print(f"Unmapped reads:  {total - mapped:,}")
    print(f"Mapping rate:    {100*mapped/total:.2f}%")
    print(f"Alignment time:  {align_time:.1f}s")
    print(f"Total time:      {total_time:.1f}s")
    print(f"k-mer size:      {K}")
    print(f"Max mismatches:  {MAX_MISMATCHES}")
    print(f"Seeds per read:  {NUM_SEEDS}")

    # Save results
    print(f"\nSaving results to {output_file}...")
    with open(output_file, 'w') as f:
        f.write("read_name\tfile\tposition\tstrand\tmismatches\tread_length\n")
        for name, fl, pos, strand, mm, rlen in results:
            f.write(f"{name}\t{fl}\t{pos}\t{strand}\t{mm}\t{rlen}\n")
    print("Done.")

if __name__ == '__main__':
    main()
