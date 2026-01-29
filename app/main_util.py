import yaml
import logging
import pandas as pd
from app.log_util import log_batch_query
from app.processing_util import prepare_input
from app.vep_util import prepare_vep_output, run_vep
from app.phenotype_util import allele_phenotype_match
from app.gene_util import fetch_gene_info, fetch_mus_alleles, fetch_homologous_gene, process_mus_alleles
from app.disease_util import assign_clinvar, fetch_mus_doid, clndisdb_to_mondo, map_doids_to_mondo

import time

# Load config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


# Load homology
homology_df = pd.read_csv(config['paths']['mgi_homology'],
                          sep='\t', names=['MGI_ID', 'MusGeneSymbol', 'MusEntrezGeneID', 
                                           'MusHGNC_ID', 'HumGeneSymbol', 'HumEntrezGeneID'])
homology_dict = dict(zip(homology_df['HumGeneSymbol'], homology_df['MusGeneSymbol']))



# ---- Batch ----
def batch_hvar(variants, assembly='GRCh38'):
    '''
    Annotate a batch of human variant using local gene metadata, Ensembl VEP, and ClinVar.

    Parameters

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
        input_gene_map
    '''
    s = time.time()
    # Check assembly
    if assembly != 'GRCh38':
        logging.info(f'Assembly Error: {assembly} is not supported. Only GRCh38 is supported.')
        raise ValueError('Assembly must be GRCh38')

    # Log the query
    log_batch_query(variants)
    
    # Prepare HGVS notation
    variants, submission_map = prepare_input(variants)


    # Query VEP
    vep_df = run_vep(variants, 'homo_sapiens')


    # Assign input field
    vep_df['Input'] = vep_df['Submission'].map(submission_map)

    # Clean VEP output into protein DataFrame
    protein_df = prepare_vep_output(vep_df, 'human')


    # Create a human gene to input mapping DataFrame
    input_gene_df = (protein_df[['Input', 'Gene Symbol']].drop_duplicates()
                     .rename(columns={'Gene Symbol': 'Hum Gene'}).copy())

    # Build gene table
    genes = protein_df['Gene Symbol'].unique()
    gene_df = fetch_gene_info(genes, species='human')


    # Disease Associations
    disease_df = assign_clinvar(variants)
    disease_df = clndisdb_to_mondo(disease_df)
    protein_df = protein_df.merge(disease_df, on='Input', how='left')

    gene_df = gene_df.where(pd.notnull(gene_df), None)
    protein_df = protein_df.where(pd.notnull(protein_df), None)

    gene_df.to_csv('app/results/hum_gene_df.csv', index=False)
    protein_df.to_csv('app/results/hum_protein_df.csv', index=False)

    # Print resulting tables
    # print(gene_df)
    # print('----------')
    # print(protein_df)
    # print('----------')
    # print(input_gene_df)
    # print('----------')
    e = time.time()

    print(f'Human Time: {e-s}')

    return gene_df, protein_df, input_gene_df


def batch_mvar(input_mapping_df, assembly='GRCm39'):
    '''
    Annotate a mouse variant using local gene metadata, Ensembl VEP, and MouseMine.

    This function performs the following steps:
        1. Fetch mouse gene metadata.
        2. Identify all variants for the gene in MGI local database.
        3. Call Ensembl VEP for each variant → transcript-level annotations.
        4. Retrieve domain names from InterPro.
        5. Query MouseMine for ontology annotations.

    Parameters
        hum_gene_df: df. Human gene DF 
        assembly : string, optional. Genome assembly. Only "GRCm39" is supported.

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
    '''
    s = time.time()
    # Check assembly
    if assembly != 'GRCm39':
        raise ValueError('Assembly must be GRCm39')
    

    # Extract Mus gene symbol
    gene_input_df = fetch_homologous_gene(input_mapping_df)

    # Build orthologous gene table
    genes = gene_input_df['Mus Gene'].unique()
    mus_gene_df = fetch_gene_info(genes, species='mouse')

    gene_input_df.to_csv('app/testing_results/gene_input_df.csv', index=False)


    # Extract gene ids
    MGI_gene_ids = mus_gene_df['Accession'].unique()

    # Fetch mouse alleles
    mouse_allele_df = fetch_mus_alleles(MGI_gene_ids)
    mouse_prt_df = process_mus_alleles(mouse_allele_df)
    

    # Query MGD for ontology associations
    doid_map = fetch_mus_doid(mouse_prt_df['AlleleID'].unique())

    # Drop alleles without disease associations
    mouse_prt_df = mouse_prt_df[mouse_prt_df['AlleleID'].isin(doid_map.keys())]


    # Perpare HGVS for VEP
    mouse_prt_df, _ = prepare_input(mouse_prt_df)
    allele_map = mouse_prt_df[['Submission', 'Gene Symbol', 'AlleleID', 'AlleleSymbol']]


    # Query VEP for each variant
    variant_vep = run_vep(mouse_prt_df, 'mus_musculus')

    # Clean VEP output
    mouse_prt_df = prepare_vep_output(variant_vep, 'mouse')

    # Drop transcripts not in homologous gene
    mouse_prt_df = mouse_prt_df[mouse_prt_df['Gene Symbol'].isin(genes)]

    # Map Alleles back to VEP output
    allele_map = allele_map.rename(columns={'Gene Symbol': 'Gene Symbol Allele'})
    mouse_prt_df = mouse_prt_df.merge(allele_map, on='Submission', how='left')

    # Only keep Allele Gene Symbol
    mouse_prt_df.drop('Gene Symbol', inplace=True, axis=1)
    mouse_prt_df.rename(columns={'Gene Symbol Allele': 'Gene Symbol'}, inplace=True)


    # Compute MONDO mapping
    mondo_id_map, mondo_term_map = map_doids_to_mondo(doid_map)

    # Insert into the variant df
    mouse_prt_df['MONDO'] = mouse_prt_df['AlleleID'].map(mondo_id_map)
    mouse_prt_df['Associated Diseases'] = mouse_prt_df['AlleleID'].map(mondo_term_map)


    # Rearrange columns
    mouse_prt_df = mouse_prt_df[['Gene Symbol', 'AlleleID', 'AlleleSymbol', 'Transcript ID', 'Biotype', 
                                 'Exon Rank', 'Pfam Domain ID', 'Pfam Domain Name', 'Molecular Consequence', 
                                 'Codon Switch', 'Amino Acids', 'refAA', 'varAA', 'Associated Diseases', 'MONDO']]
    

    # Extract phenotypes
    phenotype_df = allele_phenotype_match(mouse_allele_df)

    # Replace all null with None
    mus_gene_df = mus_gene_df.where(pd.notnull(mus_gene_df), None)
    mouse_prt_df = mouse_prt_df.where(pd.notnull(mouse_prt_df), None)
    phenotype_df = phenotype_df.where(pd.notnull(phenotype_df), None)

    # Save debugging tables to CSVs
    mus_gene_df.to_csv('app/results/mouse_gene_df.csv', index=False)
    mouse_prt_df.to_csv('app/results/mouse_protein_df.csv', index=False)

    # Print resulting tables
    # print(mus_gene_df)
    # print('----------')
    # print(mouse_prt_df)
    # print('----------')

    e = time.time()
    print(f'Mouse Time: {e-s}')

    return mus_gene_df, mouse_prt_df, phenotype_df, gene_input_df


def batch_score(hum_prt, mouse_prt, gene_inputs):
    '''
    '''
    s = time.time()
    hum_prt['MONDO_set'] = hum_prt['MONDO'].apply(set)
    mouse_prt['MONDO_set'] = mouse_prt['MONDO'].apply(set)

    expanded_df = gene_inputs.merge(hum_prt, on='Input', suffixes=('', '_human'))
    expanded_df = expanded_df.merge(mouse_prt, left_on='Mus Gene', right_on='Gene Symbol', suffixes=('_human', '_mouse'))

    score_df = pd.DataFrame()
    score_df['Input'] = expanded_df['Input']
    score_df['Allele ID'] = expanded_df['AlleleID']
    score_df['Allele Symbol'] = expanded_df['AlleleSymbol']
    score_df['Transcript ID'] = expanded_df['Transcript ID_mouse']
    score_df['Biotype Match'] = expanded_df['Biotype_human'] == expanded_df['Biotype_mouse']
    score_df['Consequence Match'] = expanded_df['Molecular Consequence_human'] == expanded_df['Molecular Consequence_mouse']
    score_df['AA Match'] = (expanded_df['refAA_human'] == expanded_df['refAA_mouse']) & (expanded_df['varAA_human'] == expanded_df['varAA_mouse'])
    score_df['AA Position Match'] = expanded_df['Amino Acids_human'] == expanded_df['Amino Acids_mouse']
    score_df['Exon Match'] = expanded_df['Exon Rank_human'] == expanded_df['Exon Rank_mouse']
    score_df['Domain Match'] = expanded_df['Pfam Domain ID_human'] == expanded_df['Pfam Domain ID_mouse']
    score_df['Disease Match'] = expanded_df.apply(lambda r: bool(r['MONDO_set_human'] & r['MONDO_set_mouse']), axis=1)

    # columns which can be attributed to score
    match_cols = ['Biotype Match', 'Consequence Match', 'AA Match', 'AA Position Match', 'Exon Match', 'Domain Match', 'Disease Match']

    # Calculate precentage of hits
    score_df['total_score'] = score_df[match_cols].sum(axis=1) / score_df[match_cols].notna().sum(axis=1) * 100


    score_df = score_df.where(pd.notnull(score_df), None)

    e = time.time()
    print(f'Score Time: {e-s}')

    return score_df