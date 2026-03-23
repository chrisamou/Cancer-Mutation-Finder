# 🧬 SeqSense: Precision Genomics Engine (V2.0)

**SeqSense** is an enterprise-grade bioinformatics pipeline and interactive web application designed to identify somatic mutations, translate biological impacts, and generate clinical diagnostic reports in real-time. 

Built for modern clinical research, SeqSense handles everything from single-patient diagnostics to massive clinical trial cohorts by comparing wild-type reference genomes against patient-derived sequences.

🌐 **[Live Application Server](https://seqsense.streamlit.app)** *(Note: Ensure this is your actual URL!)*

---

## ✨ V2.0 Enterprise Features

* **📡 Ensembl API Auto-Fetch:** Eliminates the need for manual healthy reference uploads. Simply input a target gene (e.g., *PIK3CA*, *BRAF*), and the engine autonomously queries the European Bioinformatics Institute (Ensembl) to download the official wild-type CDS transcript in real-time.
* **📄 Clinical PDF Dossiers:** Dynamically generates formatted, hospital-ready PDF diagnostic reports (via `fpdf`) alongside standard CSV data exports.
* **📊 Interactive Visualizations:** Upgraded from static charts to interactive `Plotly` dashboards, featuring hoverable mutation locus maps and cohort demographic pie charts.
* **🧬 Dynamic Sequence Alignment:** Base-by-base comparison to pinpoint exact genomic coordinates of somatic mutations.
* **🔬 Amino Acid Translation:** Extracts affected codons and utilizes BioPython to translate DNA into proteins, revealing structural impacts (e.g., `S → R`).
* **💊 Targeted Therapy Database:** Cross-references detected variants with an internal database to recommend FDA-approved inhibitors, supplemented by real-time locus data from the NIH MyGene API.
* **📦 Big Data Cohort Processing:** Ingests compressed `.zip` archives containing dozens of patient sequences, runs a batch pipeline with a dynamic progress bar, and generates a master trial dataframe.

---

## 🛠️ Technology Stack

* **Frontend & Cloud Hosting:** Streamlit, Streamlit Community Cloud
* **Bioinformatics Engine:** BioPython (`SeqIO`, `Seq`)
* **Data & Analytics:** Pandas, Plotly (`express`, `graph_objects`)
* **API Routing:** Python `requests` (NIH MyGene, Ensembl REST API)
* **Document Generation:** FPDF

---

## 🚀 Running the Engine Locally

To run SeqSense on your local machine, follow these steps:

**1. Clone the repository:**
```bash
git clone [https://github.com/your-username/Cancer-Mutation-Finder.git](https://github.com/your-username/Cancer-Mutation-Finder.git)
cd Cancer-Mutation-Finder
