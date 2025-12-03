import streamlit as st
import streamlit.components.v1 as components
import base64
from util import hvar_to_output, mvar_to_output, score_ortho_vars

# ---- TITLE ----

def get_image_base64(path):
    with open(path, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

mgi_base64 = get_image_base64('logos/mgi_logo.png')

title_html = f"""
<div style="display: flex; align-items: center; justify-content: center; gap: 30px;">
    <img src="data:image/png;base64,{mgi_base64}" alt="MGI Logo"
         style="width: 100px; height: auto;">
    <h2 style="margin: 0; font-size: 30px; font-family: Verdana, sans-serif; line-height: 1; text-align: center; font-weight: 550;">
        Human Variant Mapping: Mouse Model Discovery
    </h2>
</div>"""

components.html(title_html, height=100)



gene = st.text_input('Gene Symbol', key='gene_input')
chrom = st.text_input('Chromosome', key='chrom_input')
start = int(st.text_input('Genomic Start', value=1, key='start_input'))
stop = int(st.text_input('Genomic Stop', value=1, key='stop_input'))
ref = st.text_input('Reference Allele(s)', key='ref_input')
alt = st.text_input('Alternate Allele(s)', key='alt_input')

if st.button('Submit Query'):
    human_gene_df, human_prt_df = hvar_to_output(gene, chrom, start, stop, ref, alt)
    mouse_gene_df, mouse_prt_df = mvar_to_output(gene)

    var_scores_df = score_ortho_vars(human_prt_df, mouse_prt_df)

    st.subheader('Human Gene Information')
    st.dataframe(human_gene_df)

    st.subheader('Human Protein Information')
    st.dataframe(human_prt_df)

    st.subheader('Mouse Gene Information')
    st.dataframe(mouse_gene_df)

    st.subheader('Mouse Transcript Information')
    st.dataframe(mouse_prt_df)

    st.subheader('Variant Orthology Scores')
    st.dataframe(var_scores_df)

