import csv
import matplotlib.pyplot as plt
from Bio import SeqIO

print("Initializing Precision Oncology Diff Engine...\n")

# 1. The Pharmacogenomics Database (Drug Matcher)
# We map the exact genomic position of the mutation to the clinical treatment
DRUG_DB = {
    31: {"gene": "BRAF", "variant": "V600E", "drug": "Vemurafenib", "type": "BRAF Inhibitor"}
}

# 2. Load the specific gene targets
healthy_dna = SeqIO.read("braf_healthy.fasta", "fasta").seq
tumor_dna = SeqIO.read("braf_tumor.fasta", "fasta").seq

mutations = []

print("Scanning Patient 02 (BRAF Amplicon) for somatic mutations...\n")

# 3. Scan the sequences side-by-side
for i, (h_base, t_base) in enumerate(zip(healthy_dna, tumor_dna)):
    pos = i + 1
    if h_base != t_base:
        # Default report if we don't know the drug
        mutation_info = {
            "Position": pos,
            "Healthy_Allele": h_base,
            "Tumor_Allele": t_base,
            "Target_Gene": "Unknown Variant",
            "Recommended_Therapy": "None - Consult Oncologist"
        }
        
        # 4. THE DRUG MATCHER LOGIC
        # If the mutation position is in our database, update the clinical report!
        if pos in DRUG_DB:
            match = DRUG_DB[pos]
            mutation_info["Target_Gene"] = f"{match['gene']} {match['variant']}"
            mutation_info["Recommended_Therapy"] = f"{match['drug']} ({match['type']})"
            
            print(f"🚨 CRITICAL MATCH: {match['gene']} {match['variant']} mutation detected!")
            print(f"💊 RECOMMENDED THERAPY: {match['drug']} ({match['type']})\n")
            
        mutations.append(mutation_info)

# 5. Export the Advanced Clinical Report
csv_filename = "oncology_report.csv"
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["Position", "Healthy_Allele", "Tumor_Allele", "Target_Gene", "Recommended_Therapy"])
    writer.writeheader()
    writer.writerows(mutations)

print(f"-> [SUCCESS] Clinical CSV report saved to '{csv_filename}'")
