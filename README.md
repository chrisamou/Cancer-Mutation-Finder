# 🧬 SeqSense: Precision Genomics Engine

A bioinformatics web application built with Python and Streamlit that dynamically compares healthy and tumor DNA sequences to identify actionable somatic mutations. 

## 🚀 Features
* **Dynamic Sequence Alignment:** Upload any two FASTA files (Healthy vs. Tumor) for real-time base-by-base comparison.
* **Point Mutation Detection:** Accurately identifies indels and single nucleotide polymorphisms (SNPs).
* **Pharmacogenomics Database:** Cross-references detected genomic coordinates with a clinical database to recommend targeted therapies (e.g., matching the PIK3CA H1047R mutation with Alpelisib).
* **Clinical Reporting:** Generates a human-readable diagnostic table and a visual genomic heatmap using Matplotlib.

## 🛠️ Tech Stack
* **Language:** Python 3.14
* **Frontend:** Streamlit
* **Bioinformatics:** BioPython
* **Data Visualization:** Matplotlib

## 💻 How to Run Locally

1. Clone this repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/Cancer-Mutation-Finder.git](https://github.com/YOUR_USERNAME/Cancer-Mutation-Finder.git)
