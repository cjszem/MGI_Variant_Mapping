import re
import yaml
import logging
import pandas as pd
from log_util import log_query, log_batch_query, log_results
from vep_util import get_vep_data, fetch_vep_data, prepare_vep_output, get_domain_name
from gene_util import get_gene_info, fetch_gene_info, fetch_mus_alleles
from disease_util import get_clinvar_ids, get_disease_associations, fetch_assign_clinvar, search_clinvar, \
                         fetch_mus_disease_map, mouse_disease_associations, get_doid_name, fetch_doid_name


# Load config
with open('config.yaml') as f:
    config = yaml.safe_load(f)


# Load homology
homology_df = pd.read_csv(config['paths']['mgi_homology'],
                          sep='\t', names=['MGI_ID', 'MusGeneSymbol', 'MusEntrezGeneID', 
                                           'MusHGNC_ID', 'HumGeneSymbol', 'HumEntrezGeneID'])
homology_dict = dict(zip(homology_df['HumGeneSymbol'], homology_df['MusGeneSymbol']))


# ---- Input ----
def process_batch_query(input):
    '''
    Processes a string resulting from the batch input of variants with formatting:
    chromosome:start-end:ref/allele

    Parameters:
        input: string. Multiline string containg all variants to process.

    Returns: DataFrame. Containing Chromosome, Start, Ref, Alt, for each variant.
    '''
    # Regex pattern to match each line: gene:chromosome:start-end:ref/allele
    pattern = r'^(\w+):(\d+)-(\d+):([ACGTacgt]+)/([ACGTacgt]+)$'

    variants = []
    # Go through each line
    for line in input.splitlines():
        if line.strip() == '':
            continue
        
        # Extract inputted variant fields
        match = re.match(pattern, line)
        if match:
            variants.append({'Chromosome': match.group(1),
                             'Start': int(match.group(2)),
                             'Stop': int(match.group(3)),
                             'Ref': match.group(4),
                             'Alt': match.group(5)})
        
    # Create pandas dataframe
    variants = pd.DataFrame(variants)
    variants.to_csv('./testing_results/batch_variants.csv', index=False)

    return variants


def prepare_hgvs(variants):
    '''
    Converts a DataFrame of variants into a list of HGHVS notations for VEP querying
    
    Parameters:
        variants: DataFrame. containing Chromosome, Start, Ref, Alt.

    Returns:
        list of strings. HGVS notations (chrom:g.startRef>Alt).
    '''
    variants['HGVS'] = (variants['Chromosome'].astype(str) + ':g.' + variants['Start'].astype(str) + variants['Ref'] + '>' + variants['Alt'])

    return variants


# ---- Single ----
def hvar_to_output(gene, chrom, start, end, ref, alt, assembly='GRCh38'):
    '''
    Annotate a human variant using local gene metadata, Ensembl VEP, and ClinVar.

    This function performs the following steps:
        1. Fetch human gene metadata for the specified gene.
        2. Query Ensembl VEP REST API for transcript-level consequences.
        3. Process amino acid changes (ref/alt AA, protein position).
        4. Retrieve Pfam/InterPro domain names for protein annotations.
        5. Build two structured DataFrames:
                - gene_df: basic gene + variant information.
                - protein_df: transcript-level VEP annotations.
        6. Fetch ClinVar associations (protein changes → disease names).
        7. Map ClinVar diseases into the per-transcript annotations.

    Parameters
        gene: string. Gene symbol to extract info for.
        chrom: str. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        ref: string. Reference allele.
        alt: string. Alternate allele.
        assembly: string, optional. Genome assembly. Only "GRCh38" is supported.

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
    '''
    # Log the query
    log_query(gene, chrom, start, end, ref, alt)

    # Check assembly
    if assembly != 'GRCh38':
        logging.info(f'Assembly Error: {assembly} is not supported. Only GRCh38 is supported.')
        raise ValueError('Assembly must be GRCh38')

    # Build gene table
    gene_info = get_gene_info(gene)
    gene_df = pd.DataFrame([{'Organism': 'Human', 
                             'Gene': gene_info.get('name'), 
                             'Description': gene_info.get('description'), 
                             'HGNC ID': gene_info.get('HGNC_id'), 
                             'Ensembl ID': gene_info.get('ENSMBL_id'), 
                             'Biotype': gene_info.get('biotype'), 
                             'Strand': gene_info.get('strand'), 
                             'Chr': chrom, 
                             'Start': start, 
                             'End': end, 
                             'Ref': ref, 
                             'Alt': alt}])


    # Use VEP to get pathogenicity and protein info
    vep_df = get_vep_data(chrom, start, end, alt)
    vep_df = vep_df[vep_df['biotype'] == 'protein_coding']


    # Split amino_acids to REFAA and VARAA if in X/Y format
    aa_split = vep_df['amino_acids'].astype(str).str.split('/', expand=True)
    vep_df['refAA'] = aa_split[0]
    vep_df['varAA'] = aa_split[1]


    # Add protein location of amino_acids
    vep_df['protein_start'] = vep_df['protein_start'].astype(pd.Int64Dtype())
    vep_df['amino_acids'] = vep_df.apply(lambda row: f"{row['refAA']}{row['protein_start']}{row['varAA']}", axis=1)


    # Merge consequence terms into a single string
    vep_df['consequence_terms'] = vep_df['consequence_terms'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)


    # Fetch domain names
    for domain in vep_df['Domain'].unique():
        if domain == None:
            dom_name = None
        else:
            dom_name = get_domain_name(domain)
        vep_df.loc[vep_df['Domain'] == domain, 'domain_name'] = dom_name

    
    # Create protein DataFrame
    keep = ['transcript_id', 'biotype', 'exon', 'Domain', 'domain_name', 'polyphen_prediction', 
                'polyphen_score', 'consequence_terms', 'codons', 'amino_acids', 'refAA', 'varAA']
    protein_df = vep_df[keep]

    # Update column names
    protein_df.rename(columns={'transcript_id': 'Transcript ID', 'biotype': 'Biotype', 
                       'exon': 'Exon Rank', 'Domain': 'Pfam Domain ID', 
                       'domain_name': 'Pfam Domain Name', 'polyphen_prediction': 'Polyphen Prediction', 
                       'polyphen_score': 'Polyphen Score', 'consequence_terms': 'Molecular Consequence', 
                       'codons': 'Codon Switch', 'amino_acids': 'Amino Acids', 
                       'refAA': 'refAA', 'varAA': 'varAA'}, inplace=True)

    # Disease Associations
    clinvar_ids = get_clinvar_ids(gene, start)

    # Build lookup dictionary of protein change and consequence to diseases
    disease_dict = {}
    for clinvar in clinvar_ids: # for each ClinVar ID
        prt_changes, consequence, diseases = get_disease_associations(clinvar) # extract the protein information and diseases
        for prt_change in prt_changes:
            disease_dict[f'{prt_change}_{consequence}'] = diseases # map protein changes and consequence to diseases

    # Attach diseases to each transcript row
    associated_diseases = []
    for index, prt_row in protein_df.iterrows():
        key = f'{prt_row['Amino Acids']}_{prt_row['Molecular Consequence']}'
        d = disease_dict.get(key, None)
        if d is None: associated_diseases.append(None)
        else: associated_diseases.append('; '.join(d))

    protein_df['Associated Diseases'] = associated_diseases

    # Print resulting tables
    print(gene_df)
    print('----------')
    print(protein_df)
    print('----------')
    
    return gene_df, protein_df


def mvar_to_output(gene, assembly='GRCm39'):
    '''
    Annotate a mouse variant using local gene metadata, Ensembl VEP, and MouseMine.

    This function performs the following steps:
        1. Fetch mouse gene metadata.
        2. Identify all variants for the gene in MGI local database.
        3. Call Ensembl VEP for each variant → transcript-level annotations.
        4. Retrieve domain names from InterPro.
        5. Query MouseMine for ontology annotations.

    Parameters
        gene: string. Gene symbol to extract info for.
        assembly : string, optional. Genome assembly. Only "GRCm39" is supported.

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
    '''
    # Check assembly
    if assembly != 'GRCm39':
        raise ValueError('Assembly must be GRCm39')
    
    # Extract Mus gene symbol
    mus_gene = homology_dict.get(gene, None)
    
    # Catch for is there is no orthologous Mus gene found
    if mus_gene is None:
        logging.warning(f'No mouse ortholog found for human gene {gene}')

        # Create empty results DatFrames
        mouse_gene_df = pd.DataFrame(columns = ['Organism', 'Gene', 'Description', 'HGNC', 'Ensembl ID', 'Biotype', 'Strand'])
        mouse_prt_df = pd.DataFrame(columns = ['AlleleID', 'AlleleSymbol', 'Transcript ID', 'Biotype', 'Exon Rank', 'Pfam Domain ID', 'Pfam Domain Name', 
                                               'Molecular Consequence', 'Codon Switch', 'Amino Acids', 'refAA', 'varAA', 'ontologyName', 'ontologyID'])

        # Return empty results
        return mouse_gene_df, mouse_prt_df


    # Build orthologous gene table
    mouse_gene_info = get_gene_info(mus_gene, species='mouse')
    mouse_gene_df = pd.DataFrame([{'Organism': 'Mouse', 
                                   'Gene': mouse_gene_info.get('name'), 
                                   'Description': mouse_gene_info.get('description'), 
                                   'MGI ID': mouse_gene_info.get('MGI_id'), 
                                   'Ensembl ID': mouse_gene_info.get('ENSMBL_id'), 
                                   'Biotype': mouse_gene_info.get('biotype'), 
                                   'Strand': mouse_gene_info.get('strand')}])

    # Extract gene ids
    MGI_gene_ids = mouse_gene_df['MGI ID'].unique()

    # Fetch mouse alleles
    mouse_allele_df = fetch_mus_alleles(MGI_gene_ids)


    # Perpare HGVS for VEP
    mouse_allele_df = prepare_hgvs(mouse_allele_df)
    allele_map = mouse_allele_df[['HGVS', 'AlleleID', 'AlleleSymbol']]

    print(allele_map)

    # Query VEP for each variant
    variant_vep = fetch_vep_data(mouse_allele_df['HGVS'].tolist(), 'mouse')

    # Clean VEP output
    mouse_prt_df = prepare_vep_output(variant_vep, 'mouse')

    # Map Alleles back to VEP output
    mouse_prt_df = mouse_prt_df.merge(allele_map, on='HGVS', how='left')


    # Query MGD for ontology associations
    mouse_disease_map = fetch_mus_disease_map(mouse_prt_df['AlleleID'].unique())

    # Drop alleles without disease associations
    mouse_prt_df = mouse_prt_df[mouse_prt_df['AlleleID'].isin(mouse_disease_map.keys())]
    
    # Map diseases to mouse alleles
    mouse_prt_df['Disease Association'] = mouse_prt_df['AlleleID'].map(mouse_disease_map)
    mouse_prt_df['Disease Association'] = mouse_prt_df['Disease Association'].apply(lambda x: ','.join(x) if isinstance(x, set) else x)


    # Initialize domain name column
    mouse_prt_df['Pfam Domain Name'] = None

    # Fetch domain names
    for domain in mouse_prt_df['Pfam Domain ID'].unique():
        if pd.isna(domain):
            dom_name = None
        else:
            dom_name = get_domain_name(domain)
        mouse_prt_df.loc[mouse_prt_df['Pfam Domain ID'] == domain, 'Pfam Domain Name'] = dom_name


    # Rearrange columns
    mouse_prt_df = mouse_prt_df[['AlleleID', 'AlleleSymbol', 'Transcript ID', 'Biotype', 
                                 'Exon Rank', 'Pfam Domain ID', 'Pfam Domain Name', 'Molecular Consequence', 
                                 'Codon Switch', 'Amino Acids', 'refAA', 'varAA', 'Disease Association']]
    

    # Save debugging tables to CSVs
    mouse_allele_df.to_csv('./testing_results/mouse_allele_df.csv', index=False)

    # Print resulting tables
    print(mouse_gene_df)
    print('----------')
    print(mouse_allele_df)
    print('----------')
    print(mouse_prt_df)
    print('----------')

    return mouse_gene_df, mouse_prt_df


def score_output(human_gene_df, human_prt_df, mouse_gene_df, mouse_prt_df):
    '''
    Score orthologous mouse variant similarity against human variant.

    Parameter:
        human_prt_df: pandas.DataFrame. A single-row DataFrame describing the human variant.
        mouse_prt_df : pandas.DataFrame. A multi-row DataFrame of all orthologous mouse alleles to compare.

    Returns:
        pandas.DataFrame. Containing boolean match columns and total score for each mouse allele
    '''
    # If Mus protein DatFrame is empty, return empty scoring DataFrame
    if mouse_prt_df.empty:
        # Create empty DataFrame
        score_df = pd.DataFrame(columns=['AlleleID', 'AlleleSymbol', 'Transcript ID', 'biotype_match', 'consequence_match', 
                                     'AA_match', 'exon_match', 'domain_match', 'disease_match', 'total_score'])
        
        # Log all results
        log_results(human_gene_df, human_prt_df, mouse_gene_df, mouse_prt_df, score_df)

        return score_df

    # Extract human variant features
    biotype = human_prt_df['Biotype'][0]
    exon_rank = human_prt_df['Exon Rank'][0]
    domain = human_prt_df['Pfam Domain ID'][0]
    consequence = human_prt_df['Molecular Consequence'][0]
    refAA = human_prt_df['refAA'][0]
    altAA = human_prt_df['varAA'][0]
    amino_acids = human_prt_df['Amino Acids'][0]
    diseases = human_prt_df['Associated Diseases'][0]

    # Convert human diseases to DOID names
    diseases_doid = set()
    # Go through each correlated disease
    for disease in diseases.split('; '):
        doid = fetch_doid_name(disease)
        if doid:
            diseases_doid.add(doid) # append disease name if found
        else:
            diseases_doid.add(disease) # append original disease name if no DOID found
    
    # Create boolean score dataframe
    score_df = pd.DataFrame()
    score_df['AlleleID'] = mouse_prt_df['AlleleID']
    score_df['AlleleSymbol'] = mouse_prt_df['AlleleSymbol']
    score_df['Transcript ID'] = mouse_prt_df['Transcript ID']
    score_df['biotype_match'] = mouse_prt_df['Biotype'] == biotype
    score_df['consequence_match'] = mouse_prt_df['Molecular Consequence'] == consequence
    score_df['AA_match'] = (mouse_prt_df['refAA'] == refAA) & (mouse_prt_df['varAA'] == altAA)
    score_df['AA_Position_match'] = mouse_prt_df['Amino Acids'] == amino_acids
    score_df['exon_match'] = mouse_prt_df['Exon Rank'] == exon_rank
    score_df['domain_match'] = mouse_prt_df['Pfam Domain ID'].apply(lambda x: (x == domain) if (x is not None and domain is not None) else None)
    score_df['disease_match'] = mouse_prt_df['Disease Association'].isin(set(diseases_doid))

    # columns which can be attributed to score
    match_cols = ['biotype_match', 'consequence_match', 'AA_match', 'AA_Position_match', 'exon_match', 'domain_match', 'disease_match']

    # Calculate precentage of hits
    score_df['total_score'] = score_df[match_cols].sum(axis=1) / score_df[match_cols].notna().sum(axis=1) * 100

    # Sort results by score
    score_df = score_df.sort_values(by='total_score', ascending=False)

    print(score_df)
    print('----------')

    # Log and save results
    log_results(human_gene_df, human_prt_df, mouse_gene_df, mouse_prt_df, score_df)

    return score_df


# ---- Batch ----
def batch_hvar_to_output(variants, assembly='GRCh38'):
    '''
    Annotate a batch of human variant using local gene metadata, Ensembl VEP, and ClinVar.

    Parameters

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
    '''
    # Check assembly
    if assembly != 'GRCh38':
        logging.info(f'Assembly Error: {assembly} is not supported. Only GRCh38 is supported.')
        raise ValueError('Assembly must be GRCh38')

    # Log the query
    log_batch_query(variants)
    
    # Prepare HGVS notation
    variants = prepare_hgvs(variants)

    # Query VEP
    vep_df = fetch_vep_data(variants['HGVS'].tolist(), 'human')

    # Build gene table
    genes = vep_df['gene_symbol'].unique()
    gene_df = fetch_gene_info(genes, species='human')

    # Clean VEP output into protein DataFrame
    protein_df = prepare_vep_output(vep_df, 'human')

    # Disease Associations
    disease_df = fetch_assign_clinvar(variants)
    protein_df = protein_df.merge(disease_df, on='HGVS', how='left')
    protein_df.rename({'diseases': 'Associated Diseases'}, axis=1, inplace=True)

    # Print resulting tables
    print(gene_df)
    print('----------')
    print(protein_df)
    print('----------')
    
    return gene_df, protein_df

