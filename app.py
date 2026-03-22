import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO
import io

st.set_page_config(page_title="Universal Oncology Engine", layout="wide")
st.title("🧬 Universal Precision Oncology Engine")

st.sidebar.header("📁 Upload Genomic Data")
healthy_file = st.sidebar.file_uploader("Upload Healthy Sequence (FASTA)", type=["fasta", "fa"])
tumor_file = st.sidebar.file_uploader("Upload Tumor Sequence (FASTA)", type=["fasta", "fa"])

DRUG_DB = {
    31: {"gene": "BRAF", "variant": "V600E", "drug": "Vemurafenib", "type": "BRAF Inhibitor"},
    1047: {"gene": "PIK3CA", "variant": "H1047R", "drug": "Alpelisib", "type": "PI3K Inhibitor"}
}

if healthy_file and tumor_file:
    if st.button("Run Diagnostic Scan"):
        with st.spinner("Analyzing custom sequences..."):
            h_str = io.StringIO(healthy_file.getvalue().decode("utf-8").strip())
            t_str = io.StringIO(tumor_file.getvalue().decode("utf-8").strip())
            
            # Extract and CLEAN the sequences
            healthy_dna = str(SeqIO.read(h_str, "fasta").seq).strip().upper()
            tumor_dna = str(SeqIO.read(t_str, "fasta").seq).strip().upper()
            
            mutations = []
            min_len = min(len(healthy_dna), len(tumor_dna))
            
            for i in range(min_len):
                if healthy_dna[i] != tumor_dna[i]:
                    pos = i + 1
                    match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown", "drug": "Research Required", "type": "N/A"})
                    mutations.append({
                        "Position": pos,
                        "Healthy": healthy_dna[i],
                        "Tumor": tumor_dna[i],
                        "Gene": f"{match['gene']} {match['variant']}",
                        "Therapy": f"{match['drug']} ({match['type']})"
                    })
            
            if mutations:
                st.error(f"🚨 {len(mutations)} Somatic Mutation(s) Detected!")
                md_table = "| Position | Healthy | Tumor | Target Gene | Recommended Therapy |\n|---|---|---|---|---|\n"
                for m in mutations:
                    md_table += f"| {m['Position']} | {m['Healthy']} | {m['Tumor']} | {m['Gene']} | {m['Therapy']} |\n"
                st.markdown(md_table)
                
                fig, ax = plt.subplots(figsize=(10, 2))
                ax.plot([1, min_len], [0, 0], color='lightgray', linewidth=6)
                x_positions = [m["Position"] for m in mutations]
                ax.scatter(x_positions, [0]*len(mutations), color='red', s=150)
                ax.set_yticks([])
                ax.set_xlabel("Base Pair Position")
                st.pyplot(fig)
            else:
                # Debug info to see what the computer is actually seeing
                st.warning("No differences found. Checking sequence lengths...")
                st.write(f"Healthy Length: {len(healthy_dna)}")
                st.write(f"Tumor Length: {len(tumor_dna)}")
