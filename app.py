import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import io
import requests
import pandas as pd
import zipfile
import plotly.express as px
import plotly.graph_objects as go

# 1. Page Configuration & Branding
st.set_page_config(page_title="SeqSense Engine", layout="wide", page_icon="🧬")
st.title("🧬 SeqSense: Precision Genomics Engine")
st.write("Identify clinical variants and fetch live biological data from single patients or massive clinical cohorts.")

# 2. Local Database
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
                return {"chromosome": hit.get("map_location", "Unknown"), "description": hit.get("name", "Pending")}
    except requests.exceptions.RequestException:
        pass
    return {"chromosome": "API Offline", "description": "Offline"}

# 4. SIDEBAR: The Mode Switcher
st.sidebar.header("⚙️ Analysis Mode")
mode = st.sidebar.radio("Select Scale:", ["Single Patient", "Batch Cohort (ZIP)"])

st.sidebar.divider()

# ---------------------------------------------------------
# MODE 1: SINGLE PATIENT
# ---------------------------------------------------------
if mode == "Single Patient":
    st.sidebar.header("📁 Upload Patient Data")
    healthy_file = st.sidebar.file_uploader("Upload Reference (Healthy)", type=["fasta", "fa"])
    tumor_file = st.sidebar.file_uploader("Upload Patient Sequence", type=["fasta", "fa"])

    if healthy_file and tumor_file:
        if st.button("Run Diagnostic Scan"):
            with st.spinner("Aligning sequence and translating proteins..."):
                h_raw = io.StringIO(healthy_file.getvalue().decode("utf-8"))
                t_raw = io.StringIO(tumor_file.getvalue().decode("utf-8"))
                
                healthy_dna = str(SeqIO.read(h_raw, "fasta").seq).strip().upper()
                tumor_dna = str(SeqIO.read(t_raw, "fasta").seq).strip().upper()
                
                mutations = []
                min_len = min(len(healthy_dna), len(tumor_dna))
                
                for i in range(min_len):
                    if healthy_dna[i] != tumor_dna[i]:
                        pos = i + 1
                        codon_start = (i // 3) * 3
                        
                        h_aa = str(Seq(healthy_dna[codon_start:codon_start+3]).translate()) if len(healthy_dna[codon_start:codon_start+3]) == 3 else "?"
                        t_aa = str(Seq(tumor_dna[codon_start:codon_start+3]).translate()) if len(tumor_dna[codon_start:codon_start+3]) == 3 else "?"
                        
                        match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown", "drug": "Research Required", "type": "N/A"})
                        live_info = fetch_live_gene_data(match["gene"]) if match["gene"] != "Unknown" else {"chromosome": "N/A", "description": "N/A"}
                        
                        mutations.append({
                            "Position": pos,
                            "DNA Change": f"{healthy_dna[i]} → {tumor_dna[i]}",
                            "Protein Impact": f"{h_aa} → {t_aa}",
                            "Gene": f"{match['gene']} {match['variant']}",
                            "Therapy": f"{match['drug']} ({match['type']})",
                            "Locus": live_info["chromosome"]
                        })
                
                if mutations:
                    st.error(f"🚨 {len(mutations)} Somatic Mutation(s) Detected!")
                    st.dataframe(pd.DataFrame(mutations), use_container_width=True)
                    
                    csv = pd.DataFrame(mutations).to_csv(index=False).encode('utf-8')
                    st.download_button("📥 Download Patient Report", data=csv, file_name="patient_report.csv", mime="text/csv")
                    
                    # ⚡ NEW: INTERACTIVE PLOTLY HEATMAP
                    st.subheader("🔍 Interactive Mutation Viewer")
                    fig = go.Figure()
                    # Draw the DNA strand
                    fig.add_trace(go.Scatter(x=[1, min_len], y=[0, 0], mode='lines', line=dict(color='#E5E7EB', width=12), hoverinfo='skip'))
                    # Plot the mutations as interactive dots
                    fig.add_trace(go.Scatter(
                        x=[m["Position"] for m in mutations], y=[0]*len(mutations),
                        mode='markers', marker=dict(color='#EF4444', size=16, line=dict(color='darkred', width=2)),
                        text=[f"<b>{m['Gene']}</b><br>Pos: {m['Position']}<br>DNA: {m['DNA Change']}<br>Protein: {m['Protein Impact']}" for m in mutations],
                        hoverinfo='text'
                    ))
                    fig.update_layout(height=250, yaxis=dict(showticklabels=False, showgrid=False, zeroline=False), xaxis_title="Genomic Position (Base Pairs)", showlegend=False, plot_bgcolor='rgba(0,0,0,0)')
                    st.plotly_chart(fig, use_container_width=True)
                    
                else:
                    st.success("✅ No mutations detected.")

# ---------------------------------------------------------
# MODE 2: BATCH COHORT (The Big Data Upgrade)
# ---------------------------------------------------------
elif mode == "Batch Cohort (ZIP)":
    st.sidebar.header("📁 Upload Clinical Trial Data")
    healthy_file = st.sidebar.file_uploader("Upload Reference (Healthy)", type=["fasta", "fa"])
    zip_file = st.sidebar.file_uploader("Upload Patient Cohort (.zip)", type=["zip"])

    if healthy_file and zip_file:
        if st.button("Run Cohort Batch Scan"):
            h_raw = io.StringIO(healthy_file.getvalue().decode("utf-8"))
            healthy_dna = str(SeqIO.read(h_raw, "fasta").seq).strip().upper()
            
            master_mutations = []
            
            with zipfile.ZipFile(zip_file, 'r') as z:
                fasta_files = [f for f in z.namelist() if f.endswith(('.fasta', '.fa'))]
                total_files = len(fasta_files)
                
                st.write(f"🔍 Found {total_files} patient sequences. Scanning...")
                
                # ⚡ NEW: DYNAMIC PROGRESS BAR
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                for idx, file_name in enumerate(fasta_files):
                    status_text.text(f"Processing: {file_name} ({idx+1}/{total_files})")
                    
                    with z.open(file_name) as f:
                        t_raw = io.StringIO(f.read().decode("utf-8"))
                        tumor_dna = str(SeqIO.read(t_raw, "fasta").seq).strip().upper()
                        
                        min_len = min(len(healthy_dna), len(tumor_dna))
                        for i in range(min_len):
                            if healthy_dna[i] != tumor_dna[i]:
                                pos = i + 1
                                codon_start = (i // 3) * 3
                                h_aa = str(Seq(healthy_dna[codon_start:codon_start+3]).translate()) if len(healthy_dna[codon_start:codon_start+3]) == 3 else "?"
                                t_aa = str(Seq(tumor_dna[codon_start:codon_start+3]).translate()) if len(tumor_dna[codon_start:codon_start+3]) == 3 else "?"
                                
                                match = DRUG_DB.get(pos, {"gene": "Unknown", "variant": "Unknown"})
                                
                                master_mutations.append({
                                    "Patient ID": file_name,
                                    "Position": pos,
                                    "DNA": f"{healthy_dna[i]} → {tumor_dna[i]}",
                                    "Protein": f"{h_aa} → {t_aa}",
                                    "Target Gene": match['gene']
                                })
                    
                    # Update the progress bar
                    progress_bar.progress((idx + 1) / total_files)
            
            status_text.text("✅ Cohort scan complete!")
            
            if master_mutations:
                df_master = pd.DataFrame(master_mutations)
                
                # ⚡ NEW: COHORT DEMOGRAPHICS PIE CHART
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    st.subheader("📑 Master Cohort Table")
                    st.dataframe(df_master, use_container_width=True)
                    csv = df_master.to_csv(index=False).encode('utf-8')
                    st.download_button("📥 Download Cohort CSV", data=csv, file_name="cohort_report.csv", mime="text/csv")
                
                with col2:
                    st.subheader("📊 Demographics")
                    gene_counts = df_master['Target Gene'].value_counts().reset_index()
                    gene_counts.columns = ['Gene', 'Count']
                    pie_fig = px.pie(gene_counts, names='Gene', values='Count', hole=0.4, color_discrete_sequence=px.colors.qualitative.Pastel)
                    pie_fig.update_layout(margin=dict(t=0, b=0, l=0, r=0))
                    st.plotly_chart(pie_fig, use_container_width=True)
            else:
                st.success("✅ No mutations detected in the entire cohort.")
