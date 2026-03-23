import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO
import io
import requests
import pandas as pd

# 1. Page Configuration & Branding
st.set_page_config(page_title="SeqSense Engine", layout="wide", page_icon="🧬")
st.title("🧬 SeqSense: Precision Genomics Engine")
st.write("Upload Healthy vs. Tumor FASTA files to identify clinical variants and fetch live biological data.")

# 2. Local Database (The Clinical Map)
DRUG_DB = {
    31: {"gene": "BRAF", "variant": "V600E", "drug": "Vemurafenib", "type": "BRAF Inhibitor"},
    1047: {"gene": "PIK3CA", "variant": "H1047R", "drug": "Alpelisib", "type": "PI3K Inhibitor"}
}

# 3. LIVE API CONNECTION
@st.cache_data
def fetch_live_gene_data(gene_symbol):
    url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&fields=symbol,name,map_location"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            if 'hits' in data and len(data['hits']) > 0:
                hit = data['hits'][0]
                return {
                    "chromosome": hit.get("map_location", "Unknown"),
                    "description": hit.get("name", "Description pending")
                }
    except requests.exceptions.RequestException:
        pass
    return {"chromosome": "API Offline", "description": "Could not connect to live database."}

# 4. SIDEBAR: File Uploaders
st.sidebar.header("📁 Upload Genomic Data")
healthy_file = st.sidebar.file_uploader("Upload Healthy Sequence", type=["fasta", "fa"])
tumor_file = st.sidebar.file_uploader("Upload Tumor Sequence", type=["fasta", "fa"])

# 5. Main Application Logic
if healthy_file and tumor_file:
    if st.button("Run Diagnostic Scan"):
        with st.spinner("Aligning sequences and querying NIH databases..."):
            
            h_raw = io.StringIO(healthy_file.getvalue().decode("utf-8"))
            t_raw = io.StringIO(tumor_file.getvalue().decode("utf-8"))
            
            healthy_dna = str(SeqIO.read(h_raw, "fasta").seq).strip().upper()
            tumor_dna = str(SeqIO.read(t_raw, "fasta").seq).strip().upper()
            
            mutations = []
            min_len = min(len(healthy_dna), len(tumor_dna))
            
            for i in range(min_len):
                if healthy_dna[i] != tumor_dna[i]:
                    pos = i + 1
                    match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown", "drug": "Research Required", "type": "N/A"})
                    
                    live_info = {"chromosome": "N/A", "description": "N/A"}
                    if match["gene"] != "Unknown":
                        live_info = fetch_live_gene_data(match["gene"])
                    
                    mutations.append({
                        "Position": pos,
                        "Mutation": f"{healthy_dna[i]} -> {tumor_dna[i]}",
                        "Gene": f"{match['gene']} {match['variant']}",
                        "Therapy": f"{match['drug']} ({match['type']})",
                        "Chromosome": live_info["chromosome"],
                        "Description": live_info["description"]
                    })
            
            # 6. Display Dynamic Results & Export
            if mutations:
                st.error(f"🚨 {len(mutations)} Somatic Mutation(s) Detected!")
                
                st.subheader("📑 Clinical Diagnostic Report (Live Data)")
                md_table = "| Position | Change | Target Gene | 📍 Locus | 📖 Protein Name | Recommended Therapy |\n|---|---|---|---|---|---|\n"
                for m in mutations:
                    md_table += f"| {m['Position']} | **{m['Mutation']}** | {m['Gene']} | {m['Chromosome']} | _{m['Description']}_ | **{m['Therapy']}** |\n"
                st.markdown(md_table)
                
                # ⚡ NEW FEATURE: CSV EXPORT
                st.subheader("💾 Export Clinical Data")
                df = pd.DataFrame(mutations)
                csv = df.to_csv(index=False).encode('utf-8')
                
                st.download_button(
                    label="📥 Download Report as CSV",
                    data=csv,
                    file_name="seqsense_diagnostic_report.csv",
                    mime="text/csv",
                )
                
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
else:
    st.info("👈 Please upload both FASTA files in the sidebar to begin.")
