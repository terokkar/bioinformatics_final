from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- CONFIG ---
fastq_file = "second_trim_2.fastq"
max_reads = 5000

# --- Initialize ---
quality_scores = []
base_counts = {"A": [], "T": [], "G": [], "C": [], "N": []}

# --- Parse FASTQ and collect data ---
for i, record in enumerate(SeqIO.parse(fastq_file, "fastq")):
    if i >= max_reads:
        break
    quals = record.letter_annotations["phred_quality"]
    seq = record.seq

    quality_scores.append(quals)

    for pos, base in enumerate(seq):
        if base not in base_counts:
            base_counts[base] = []
        while len(base_counts[base]) <= pos:
            base_counts[base].append(0)
        base_counts[base][pos] += 1

# --- Plot 1: Quality Score Boxplot ---
max_len = max(len(q) for q in quality_scores)
quality_padded = [q + [np.nan] * (max_len - len(q)) for q in quality_scores]
quality_array = np.array(quality_padded)

plt.figure(figsize=(12, 6))
sns.boxplot(data=quality_array, showfliers=False)
plt.xlabel("Base Position")
plt.ylabel("Phred Quality Score")
plt.title("Quality Scores Across Read Positions (Pre-Filtering)")
plt.tight_layout()
plt.show()

# --- Plot 2: Base Composition Line Plot ---
plt.figure(figsize=(12, 6))
for base, counts in base_counts.items():
    plt.plot(range(len(counts)), counts, label=base)
plt.xlabel("Base Position")
plt.ylabel("Base Count")
plt.title("Base Composition by Position (Pre-Filtering)")
plt.legend()
plt.tight_layout()
plt.show()
