import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO
import io

# 1. Page Configuration
st.set_page_config(page_title="Universal Oncology Engine", layout="wide")
st.title("🧬 Universal Precision Oncology Engine")
st.write("Upload Healthy vs. Tumor FASTA files to identify clinical variants.")

# 2. Pharmacogenomics Database (The Clinical Map)
DRUG_DB = {
    31: {"gene": "BRAF", "variant": "V600E", "drug": "Vemurafenib", "type": "BRAF Inhibitor"},
    1047: {"gene": "PIK3CA", "variant": "H1047R", "drug": "Alpelisib", "type": "PI3K Inhibitor"}
}

# 3. SIDEBAR: File Uploaders
st.sidebar.header("📁 Upload Genomic Data")
healthy_file = st.sidebar.file_uploader("Upload Healthy Sequence", type=["fasta", "fa"])
tumor_file = st.sidebar.file_uploader("Upload Tumor Sequence", type=["fasta", "fa"])

# 4. Main Application Logic
if healthy_file and tumor_file:
    if st.button("Run Diagnostic Scan"):
        with st.spinner("Analyzing sequences..."):
            
            # Read and Clean Data: Removing any hidden formatting issues
            h_raw = io.StringIO(healthy_file.getvalue().decode("utf-8"))
            t_raw = io.StringIO(tumor_file.getvalue().decode("utf-8"))
            
            # Using str().strip().upper() ensures the comparison is 100% accurate
            healthy_dna = str(SeqIO.read(h_raw, "fasta").seq).strip().upper()
            tumor_dna = str(SeqIO.read(t_raw, "fasta").seq).strip().upper()
            
            mutations = []
            min_len = min(len(healthy_dna), len(tumor_dna))
            
            # Base-by-Base Comparison
            for i in range(min_len):
                if healthy_dna[i] != tumor_dna[i]:
                    pos = i + 1
                    # Look up position in our Clinical Database
                    match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown", "drug": "Research Required", "type": "N/A"})
                    
                    mutations.append({
                        "Position": pos,
                        "Healthy": healthy_dna[i],
                        "Tumor": tumor_dna[i],
                        "Gene": f"{match['gene']} {match['variant']}",
                        "Therapy": f"{match['drug']} ({match['type']})"
                    })
            
            # 5. Display Results
            if mutations:
                st.error(f"🚨 {len(mutations)} Somatic Mutation(s) Detected!")
                
                st.subheader("📑 Clinical Diagnostic Report")
                md_table = "| Position | Healthy | Tumor | Target Gene | Recommended Therapy |\n|---|---|---|---|---|\n"
                for m in mutations:
                    md_table += f"| {m['Position']} | {m['Healthy']} | {m['Tumor']} | {m['Gene']} | {m['Therapy']} |\n"
                st.markdown(md_table)
                
                # Visual Heatmap
                st.subheader("📊 Genomic Mutation Heatmap")
                fig, ax = plt.subplots(figsize=(10, 2))
                ax.plot([1, min_len], [0, 0], color='lightgray', linewidth=6, zorder=1)
                
                x_positions = [m["Position"] for m in mutations]
                ax.scatter(x_positions, [0]*len(mutations), color='red', s=150, zorder=2)
                
                ax.set_yticks([])
                ax.set_xlabel("Genomic Position (Base Pairs)")
                st.pyplot(fig)
            else:
                st.success("✅ No mutations detected. Sequences are identical.")
                st.info(f"Analyzed {min_len} base pairs.")
else:
    st.info("👈 Please upload both FASTA files in the sidebar to begin.")
