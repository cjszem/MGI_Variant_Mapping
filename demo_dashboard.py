import streamlit as st
import streamlit.components.v1 as components
import base64
from util import hvar_to_output, mvar_to_output, score_output, process_batch_query, batch_hvar_to_output


# ---- CONFIGURATION ----
st.set_page_config(page_title="Variant Mapping Query", layout="wide", )


# ---- TITLE ----
# Extract the abse64 string of MGI logo for html display
def get_image_base64(path):
    with open(path, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

mgi_base64 = get_image_base64('logos/mgi_logo.png')

# Title and Logo
title_html = f"""
<div style="display: flex; align-items: center; justify-content: center; gap: 20px; flex-wrap: wrap;">
    <h2 style="margin: 0; font-size: 28px; font-family: Verdana, sans-serif; font-weight: 550; text-align: center;">
        Human Variant Mapping:<br>
        Mouse Model Discovery
    </h2>
    <img src="data:image/png;base64,{mgi_base64}" alt="MGI Logo"
         style="width: 80px; height: auto;">
</div>
"""
components.html(title_html, height=100)

single_tab, batch_tab = st.tabs(['Single Variant', 'Batch Variants'])

with single_tab:
    # ---- INPUT FIELDS ----
    # Use columns to prevent input fields from stretching across whole screen
    col1, col2, col3 = st.columns([1, 3, 1])  # Middle column is wider for inputs

    with col2:
        # Input fields
        gene = st.text_input('Gene Symbol', key='gene_input')
        chrom = st.text_input('Chromosome', key='chrom_input')
        start = int(st.text_input('Genomic Start', value=1, key='start_input'))
        stop = int(st.text_input('Genomic Stop', value=1, key='stop_input'))
        ref = st.text_input('Reference Allele(s)', key='ref_input')
        alt = st.text_input('Alternate Allele(s)', key='alt_input')

        # Assembly Warning
        warning_html = """
        <div style="display: flex; align-items: center; justify-content: center; max-width: 800px; margin: auto;">
            <p style="margin: 0; font-size: 11px; font-family: Verdana, sans-serif; text-align: center;">
                Please note that all genomic coordinates should be based on the GRCh38 assembly. All mouse results will be based on the GRCm39 assembly.
                If you are using a different assembly, please convert your coordinates before continuing.
            </p>
        </div>
        """
        components.html(warning_html, height=50)

        # Submit button
        submit_single = st.button('Submit Query')


    # ---- OUTPUT ----
    if submit_single:

        human_gene_df, human_prt_df = hvar_to_output(gene, chrom, start, stop, ref, alt)
        mouse_gene_df, mouse_prt_df = mvar_to_output(gene)

        var_scores_df = score_output(human_gene_df, human_prt_df, mouse_gene_df, mouse_prt_df)

        # HUMAN
        st.subheader('Human Gene Information')
        st.dataframe(human_gene_df, hide_index=True)

        st.subheader('Human Protein Information')
        st.dataframe(human_prt_df, hide_index=True)

        # MOUSE
        st.subheader('Mouse Gene Information')
        st.dataframe(mouse_gene_df, hide_index=True)

        st.subheader('Mouse Protein Information')
        mouse_prt_df['AlleleSymbol'] = mouse_prt_df['AlleleSymbol'].str.replace(r'<sup>(.*?)</sup>', r'[\1]', regex=True) # Make allele symbol readable - superscript in brackets
        st.dataframe(mouse_prt_df, hide_index=True)

        # SCORES
        st.subheader('Allele Equivalence Scores')
        var_scores_df['AlleleSymbol'] = var_scores_df['AlleleSymbol'].str.replace(r'<sup>(.*?)</sup>', r'[\1]', regex=True) # Make allele symbol readable - superscript in brackets
        st.dataframe(var_scores_df, hide_index=True)

with batch_tab:
    batch_input = st.text_area('Input all variants here, one per line, in the format: chromosome:start-end:ref/alt', key='input',
                               height=300, placeholder='Example:\n2:157774114-157774114:C/T\n1:39468726-39468726:T/G')

    # Submit button
    submit_batch = st.button('Submit Query', key='batch')

    if submit_batch:
        variants = process_batch_query(batch_input)
        hum_gene_df, hum_prt_df = batch_hvar_to_output(variants)

        st.dataframe(hum_gene_df, hide_index=True)
        st.dataframe(hum_prt_df, hide_index=True)