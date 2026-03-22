import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO
import io

# 1. Page Configuration
st.set_page_config(page_title="Universal Oncology Engine", layout="wide")
st.title("🧬 Universal Precision Oncology Engine")
st.write("Upload any two FASTA files to compare sequences and identify clinical variants.")

# 2. SIDEBAR: File Uploaders
st.sidebar.header("📁 Upload Genomic Data")
healthy_file = st.sidebar.file_uploader("Upload Healthy Sequence (FASTA)", type=["fasta", "fa"])
tumor_file = st.sidebar.file_uploader("Upload Tumor Sequence (FASTA)", type=["fasta", "fa"])

# 3. Pharmacogenomics Database (Updated with PIK3CA)
DRUG_DB = {
    31: {"gene": "BRAF", "variant": "V600E", "drug": "Vemurafenib", "type": "BRAF Inhibitor"},
    1047: {"gene": "PIK3CA", "variant": "H1047R", "drug": "Alpelisib", "type": "PI3K Inhibitor"}
}

# 4. Main Application Logic
if healthy_file and tumor_file:
    if st.button("Run Diagnostic Scan"):
        with st.spinner("Analyzing custom sequences..."):
            
            # Convert uploaded bytes into BioPython readable format
            h_str = io.StringIO(healthy_file.getvalue().decode("utf-8"))
            t_str = io.StringIO(tumor_file.getvalue().decode("utf-8"))
            
            healthy_dna = SeqIO.read(h_str, "fasta").seq
            tumor_dna = SeqIO.read(t_str, "fasta").seq
            
            mutations = []
            
            # Compare sequences side-by-side
            min_len = min(len(healthy_dna), len(tumor_dna))
            for i in range(min_len):
                h_base = healthy_dna[i]
                t_base = tumor_dna[i]
                pos = i + 1
                if h_base != t_base:
                    # Check the database for a match
                    match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown", "drug": "Research Required", "type": "N/A"})
                    
                    mutations.append({
                        "Position": pos,
                        "Healthy": h_base,
                        "Tumor": t_base,
                        "Gene": f"{match['gene']} {match['variant']}",
                        "Therapy": f"{match['drug']} ({match['type']})"
                    })
            
            # 5. Display Results
            if mutations:
                st.error(f"🚨 {len(mutations)} Somatic Mutation(s) Detected!")
                
                st.subheader("📑 Clinical Diagnostic Report")
                # Bulletproof Markdown Table
                md_table = "| Position | Healthy | Tumor | Target Gene | Recommended Therapy |\n|---|---|---|---|---|\n"
                for m in mutations:
                    md_table += f"| {m['Position']} | {m['Healthy']} | {m['Tumor']} | {m['Gene']} | {m['Therapy']} |\n"
                st.markdown(md_table)
                
                st.subheader("📊 Genomic Mutation Heatmap")
                fig, ax = plt.subplots(figsize=(10, 2))
                ax.plot([1, min_len], [0, 0], color='lightgray', linewidth=6, zorder=1)
                
                x_positions = [m["Position"] for m in mutations]
                ax.scatter(x_positions, [0]*len(mutations), color='red', s=150, zorder=2)
                
                ax.set_yticks([])
                ax.set_xlabel("Genomic Position (Base Pairs)")
                st.pyplot(fig)
            else:
                st.success("✅ No mutations detected. Tissue is genetically identical.")
else:
    st.info("👈 Please upload both a Healthy and a Tumor FASTA file in the sidebar to begin.")
