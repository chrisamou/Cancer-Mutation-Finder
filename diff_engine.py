import csv
import matplotlib.pyplot as plt
from Bio import SeqIO

print("Initializing Advanced Cancer Mutation Diff Engine...\n")

# 1. Load the DNA
healthy_dna = SeqIO.read("healthy.fasta", "fasta").seq
tumor_dna = SeqIO.read("tumor.fasta", "fasta").seq

mutations = []

# 2. Scan for Mutations
for i, (h_base, t_base) in enumerate(zip(healthy_dna, tumor_dna)):
    if h_base != t_base:
        mutations.append({
            "Position": i + 1,
            "Healthy_Allele": h_base,
            "Tumor_Allele": t_base
        })

print(f"Scan complete. Found {len(mutations)} mutation(s).\n")

# 3. Option 1: Generate the Clinical CSV Report
csv_filename = "clinical_report.csv"
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["Position", "Healthy_Allele", "Tumor_Allele"])
    writer.writeheader()
    writer.writerows(mutations)

print(f"-> [SUCCESS] Clinical CSV report saved to '{csv_filename}'")

# 4. Option 2: Generate the Visual Heatmap
plt.figure(figsize=(10, 2))

import csv
import matplotlib.pyplot as plt
from Bio import SeqIO

print("Initializing Advanced Cancer Mutation Diff Engine...\n")

# 1. Load the DNA
healthy_dna = SeqIO.read("healthy.fasta", "fasta").seq
tumor_dna = SeqIO.read("tumor.fasta", "fasta").seq

mutations = []

# 2. Scan for Mutations
for i, (h_base, t_base) in enumerate(zip(healthy_dna, tumor_dna)):
    if h_base != t_base:
        mutations.append({
            "Position": i + 1,
            "Healthy_Allele": h_base,
            "Tumor_Allele": t_base
        })

print(f"Scan complete. Found {len(mutations)} mutation(s).\n")

# 3. Option 1: Generate the Clinical CSV Report
csv_filename = "clinical_report.csv"
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["Position", "Healthy_Allele", "Tumor_Allele"])
    writer.writeheader()
    writer.writerows(mutations)

print(f"-> [SUCCESS] Clinical CSV report saved to '{csv_filename}'")

# 4. Option 2: Generate the Visual Heatmap
plt.figure(figsize=(10, 2))

# Draw the healthy DNA backbone as a gray line
plt.plot([1, len(healthy_dna)], [0, 0], color='lightgray', linewidth=6, zorder=1, label="Healthy Genome")

# Drop a red dot wherever a mutation was found
if mutations:
    x_positions = [m["Position"] for m in mutations]
    y_positions = [0] * len(mutations)
    plt.scatter(x_positions, y_positions, color='red', s=150, zorder=2, label="Tumor Mutation")

plt.title("Patient 01: Genomic Mutation Heatmap")
plt.xlabel("Genomic Position (Base Pairs)")
plt.yticks([]) # Hide the Y axis since it's just a 1D map
plt.legend(loc="upper right")
plt.tight_layout()

image_filename = "mutation_heatmap.png"
plt.savefig(image_filename)
print(f"-> [SUCCESS] Visual heatmap saved to '{image_filename}'\n")
