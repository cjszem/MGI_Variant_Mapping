import streamlit as st
import streamlit.components.v1 as components
import base64
from app.main_util import batch_hvar, batch_mvar, batch_score
from app.processing_util import process_batch_query


# ---- CONFIGURATION ----
st.set_page_config(page_title='Variant Mapping Query', layout='wide')


# ---- TITLE ----
# Extract the abse64 string of MGI logo for html display
def get_image_base64(path):
    with open(path, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode()

mgi_base64 = get_image_base64('static/logos/mgi_logo.png')

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

human_tab,  = st.tabs(['Human Variant Search'])

with human_tab:
    batch_input = st.text_area('Input all variants here, one per line, in the format: chromosome:start-end:ref/alt', key='input',
                               height=300, placeholder='Example:\n2:157774114-157774114:C/T\n1:39468726-39468726:T/G')

    # Submit button
    submit_batch = st.button('Submit Query', key='batch')

    if submit_batch:
        variants = process_batch_query(batch_input)
        hum_gene_df, hum_prt_df, input_gene_df = batch_hvar(variants)
        mouse_gene_df, mouse_prt_df, phenotype_df, gene_input_df = batch_mvar(input_gene_df)
        score_df = batch_score(hum_prt_df, mouse_prt_df, gene_input_df)


        st.dataframe(hum_gene_df, hide_index=True)
        st.dataframe(hum_prt_df, hide_index=True)
        st.dataframe(mouse_gene_df, hide_index=True)
        st.dataframe(mouse_prt_df, hide_index=True)
        st.dataframe(phenotype_df, hide_index=True)
        st.dataframe(score_df, hide_index=True)