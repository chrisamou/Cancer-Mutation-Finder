import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO

# 1. Set up the Webpage style
st.set_page_config(page_title="Oncology Engine", layout="wide")
st.title("🧬 Precision Oncology Diff Engine")
st.write("Upload patient DNA sequences to identify somatic mutations and recommended targeted therapies.")

# 2. The Pharmacogenomics Database
DRUG_DB = {
    31: {"gene": "BRAF", "variant": "V600E", "drug": "Vemurafenib", "type": "BRAF Inhibitor"}
}

# 3. Create a big "Run" button on the website
if st.button("Run Genomic Scan"):
    
    with st.spinner("Aligning sequences and scanning for mutations..."):
        # Load the files
        healthy_dna = SeqIO.read("braf_healthy.fasta", "fasta").seq
        tumor_dna = SeqIO.read("braf_tumor.fasta", "fasta").seq
        
        mutations = []
        
        # Scan the sequences side-by-side
        for i, (h_base, t_base) in enumerate(zip(healthy_dna, tumor_dna)):
            pos = i + 1
            if h_base != t_base:
                # Check our Drug Database
                match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown", "drug": "Consult Oncologist", "type": "N/A"})
                
                mutations.append({
                    "Position": pos,
                    "Healthy Allele": h_base,
                    "Tumor Allele": t_base,
                    "Target Gene": f"{match['gene']} {match['variant']}",
                    "Recommended Therapy": f"{match['drug']} ({match['type']})"
                })
        
        # 4. Display the Results on the Webpage!
        if mutations:
            st.error(f"🚨 CRITICAL ALERT: {len(mutations)} Somatic Mutation(s) Detected!")
            
            # THE BULLETPROOF FIX: Draw a Markdown Table instead of a Pandas/PyArrow table
            st.subheader("📑 Clinical Diagnostic Report")
            
            # We literally just write the table as a text string!
            md_table = "| Position | Healthy | Tumor | Target Gene | Recommended Therapy |\n"
            md_table += "|---|---|---|---|---|\n"
            for m in mutations:
                md_table += f"| {m['Position']} | {m['Healthy Allele']} | {m['Tumor Allele']} | {m['Target Gene']} | {m['Recommended Therapy']} |\n"
            
            st.markdown(md_table)
            
            # Draw the Visual Heatmap
            st.subheader("📊 Genomic Mutation Heatmap")
            fig, ax = plt.subplots(figsize=(10, 2))
            ax.plot([1, len(healthy_dna)], [0, 0], color='lightgray', linewidth=6, zorder=1)
            
            x_positions = [m["Position"] for m in mutations]
            ax.scatter(x_positions, [0]*len(mutations), color='red', s=150, zorder=2)
            
            ax.set_yticks([])
            ax.set_xlabel("Genomic Position (Base Pairs)")
            st.pyplot(fig) # Tells Streamlit to draw the Matplotlib graph!
            
        else:
            st.success("✅ No mutations detected. Tissue is genetically identical.")
