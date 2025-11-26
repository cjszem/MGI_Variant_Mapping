'''
Cole Szeman
The Jackson Laboratory
2025
Utility functions for Mouse and Human Variant Annotation
'''

import re
import requests
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET

# Setup API Information
Entrez.email = "cole.szeman@jax.org"
ensembl_base_url = 'https://rest.ensembl.org'
interpro_base_url = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/'
mousemine_base_url = 'https://www.mousemine.org/mousemine/service/query/results'
ols_base_url = 'https://www.ebi.ac.uk/ols4/api/'

# Load MGI Allele Data
mouse_vars = pd.read_csv('VARIANT-ALLELE-MGI.tsv', sep='\t')


# ---- Helper functions ----
def get_gene_info(gene, species='human'):
    '''
    Query Ensembl REST API and extract info for a given gene.

    Parameters:
        gene: string. Gene symbol to extract info for.
        species: string. Species to extract gene from. Defaults to 'human'.

    Returns:
        dictionary. Containins HGNC ID, Esnembl ID, name, strand, biotype, and description.
    '''
    # Construct URL
    url = f'{ensembl_base_url}/lookup/symbol/{species}/{gene}?content-type=application/json'

    # Make request to Ensemble REST API
    request = requests.get(url, headers={'Content-Type': 'application/json'})
    request.raise_for_status()
    data = request.json()

    # Extract relevant data
    ens_id = data.get('id')
    name = data.get('display_name')
    biotype = data.get('biotype')
    strand = data.get('strand')
    match = re.match(r'^(.*?) \[Source:(.+) Symbol;Acc:(.+)\]$', data.get('description')) # Split string to access description and HGNC ID
    description = match.group(1).strip()
    hgnc_id = match.group(3).strip()

    # Create gene dictionary
    gene_info = {'HGNC': hgnc_id, 'id': ens_id, 'name': name, 'strand': strand, 'biotype': biotype, 'description': description}

    return gene_info

def get_domain_name(pfscan_domain_id):
    '''
    Query InterPro Pfam database to get domain name.

    Parameters:
        pfscan_domain_id: string. PFScan domain ID to extract full name for.

    Returns:
        string. The corresponding long name for given domain ID.
    '''
    # Construct URL
    url = interpro_base_url + pfscan_domain_id
    headers = {'Accept': 'application/json'}

    try:
        # Make request to InterPro REST API
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        json = r.json()

        # Extract domain long name
        metadata = json['metadata']
        name = metadata['name'].get('name') or metadata['name'].get('value')

        # Return domain long name
        return name
    
    # Return None if the request fails
    except Exception:
        return None

def get_vep_data(chromosome, start, end, alt, species='human'):
    '''
    Use Ensembl VEP REST endpoint to get consequence data.

    Parameters:
        chrom: int. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        alt: string. Alternate allele
        species: string. Species to extract gene from. Defaults to 'human'.

    Returns:
        string. The corresponding long name for given domain ID.
    '''
    # Build region string | chromosome:start-stop/alt
    region = f'{chromosome}:{start}-{end}/{alt}'

    # Construct URL
    url = f'{ensembl_base_url}/vep/{species}/region/{region}?numbers=1&domains=1'
    
    # Only choose one reference transcript for human variants
    # Example url: https://rest.ensembl.org/vep/human/region/2:157774114-157774114/T?numbers=1&domains=1&pick_order=mane_select,length&pick=1
    if species == 'human': 
        url += '&pick_order=mane_select,length&pick=1'

    try:
        # Make request to VEP
        request = requests.get(url, headers={'Content-Type': 'application/json', 'Accept': 'application/json'}, timeout=15)
        request.raise_for_status()
        data = request.json()

        # If no data is found, return empty DataFrame
        if not data:
            return pd.DataFrame()
        
        variant = data[0] # The first element contains all variant information
        consequences = variant.get('transcript_consequences', []) # Create a list of vep results

        # If no VEP results, return empty DataFrame
        if len(consequences) == 0:
            print('VEP yielded no results')
            return pd.DataFrame()
        
        # Normalize JSON result into a DataFrame
        vep_df = pd.json_normalize(consequences)

        # Keep desired columns that exist in DataFrame
        keep = ['transcript_id', 'polyphen_prediction', 'polyphen_score', 'amino_acids',
                'protein_start', 'protein_end', 'consequence_terms', 'exon', 'domains', 
                'codons', 'impact', 'biotype']
        vep_df = vep_df[[col for col in keep if col in vep_df.columns]]
        
        # Extract pfam names from dictionary of domain names and databases
        if 'domains' in vep_df.columns: 
            for domain in vep_df['domains'][0]:
                if isinstance(domain, dict) and domain.get('db') == 'Pfam':
                    pfam = domain['name']

            vep_df['domains'] = pfam

        # Return VEP DataFrame
        return vep_df
    
    # Handle request failure
    except Exception as e:
        print('VEP request failed:', e)
        return pd.DataFrame()

def get_clinvar_ids(gene, start):
    '''
    Use Entrez to fetch ClinVar IDs for a given gene and location.

    Parameters:
        gene: string. Gene symbol to extract info for.
        start: int. Genomic location of variant start.

    Returns:
        list of strings. The ClinVar IDs of variants at given location.
    '''
    # Make request to ClinVar
    handle = Entrez.esearch(db='clinvar', term=f'{gene}[Gene Name] AND {start}[Base Position]', retmax=500) 
    record = Entrez.read(handle)
    handle.close()

    # Extract ClinVar IDs
    clinvar_ids = record['IdList']

    return clinvar_ids

def get_disease_associations(clinvar_id):
    '''
    Use Entrez to fetch disease associations for given ClinVar ID.

    Parameters:
        clinvar_id: string. Numerical ClinVar ID.

    Returns:
        tuple
            prt_changes. list of strings. All reported protein changes associated with the variant 
            consequence. string. The first molecular consequence term extracted from the record. Returns None if none found.
            diseases. list of strings. A list of disease/condition names associated with the variant.
    '''
    # Fetch the VCV record from ClinVar
    handle = Entrez.efetch(db='clinvar', rettype='vcv', id=clinvar_id, from_esearch=True)
    xml_text = handle.read()
    handle.close()

    # Parse XML into an element tree
    clinvar_root = ET.fromstring(xml_text)

    # Extract all protein changes
    prt_changes = []
    for prt_section in clinvar_root.findall(".//ProteinChange"):
        prt = prt_section.text
        prt_changes.append(prt)

    # Extract the first molecular consequence
    for mol_cons in clinvar_root.findall(".//MolecularConsequence"):
        try:
            consequence = mol_cons.get("Type").replace(" ", "_")
            break # Only take the first available consequence
        except: 
            continue

    # Extract associated diseases
    diseases = []
    for rcv in clinvar_root.findall(".//RCVAccession"):

        # If only want to retain associations of certain review status
        # accepted_review_status = ['practice guideline', 'reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts']
        # review = rcv.find(".//ReviewStatus").text
        # if review not in accepted_review_status:
        #     continue

        # Get disease
        condition = rcv.find(".//ClassifiedCondition") 
        disease = condition.text

        # Get submission count
        germline_description = rcv.find(".//Description") 
        submissions = germline_description.get("SubmissionCount")

        diseases.append(disease)
    
    return prt_changes, consequence, diseases

def mouse_disease_associations(gene):
    '''
    Query MouseMine for allele and ontology associations for a given mouse gene symbol.

    Parameters:
        gene: string. Gene symbol to extract info for.
    
    Returns:
        pandas.DataFrame. Containins allele and ontology annotation records associated with gene. 
            - Gene.alleles.symbol
            - Gene.alleles.name
            - Gene.alleles.primaryIdentifier
            - Gene.ontologyAnnotations.ontologyTerm.name
            - Gene.ontologyAnnotations.ontologyTerm.identifier
            - Gene.ontologyAnnotations.subject.symbol
            - Gene.ontologyAnnotations.qualifier
    '''
    
    # Construct query XML
    query_xml = f'''
    <query name="HMM Search" model="genomic" 
        view="Gene.alleles.symbol 
            Gene.alleles.name 
            Gene.alleles.primaryIdentifier 
            Gene.ontologyAnnotations.ontologyTerm.name 
            Gene.ontologyAnnotations.ontologyTerm.identifier" 
        longDescription="" 
        sortOrder="Gene.alleles.symbol asc" 
        constraintLogic="A">
        <constraint path="Gene.symbol" code="A" op="=" value="{gene}"/>
    </query>
    '''

    # Make request to MouseMine
    response = requests.get(
        mousemine_base_url,
        params={'format': 'json',
                'query': query_xml})
    data = response.json()

    # Create results DataFrame
    ontology_df = pd.DataFrame(data['results'], columns=data['columnHeaders'])

    # Rename columns
    ontology_df.rename(columns={'Gene > Alleles > Symbol':'symbol', 
                                'Gene > Alleles > Name':'name', 
                                'Gene > Alleles > Primary Identifier':'primaryIdentifier', 
                                'Gene > Ontology Annotations > Term Name':'ontologyName', 
                                'Gene > Ontology Annotations > Ontology Term . Identifier':'ontologyID'}, inplace=True)

    # Filter DataFrame to only retain DOID ontology terms
    ontology_df = ontology_df[ontology_df['ontologyID'].str.startswith('DOID:')]

    return ontology_df

def get_doid_name(disease):
    '''
    Convert disease name to DOID using OLS API.

    Parameters:
        disease: string. Disease name to convert.

    Returns:
        dictionary. Containing label and ID of disease. Returns None if no match found.
    '''
    # Construct URL
    url = ols_base_url + f'search?q={disease}&ontology=doid&exact=true'

    # Make request to OLS
    request = requests.get(url)
    request.raise_for_status()

    # Extract data
    data = request.json()
    docs = data.get("response", {}).get("docs", [])

    # If no results return None
    if not docs:
        return None
    
    # Extract best match
    doc = docs[0]

    return {'label': doc.get('label'),
            'DOID': doc.get('obo_id')}
    

# ---- Main functions ----
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
        chrom: int. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        ref: string. Reference allele.
        alt: string. Alternate allele.
        assembly: string, optional. Genome assembly. Only "GRCh38" is supported.

    Returns
        gene_df: pandas.DataFrame. Gene-level metadata and input variant summary.
        protein_df: pandas.DataFrame. Transcript-level information.
    '''
    # Check assembly
    if assembly != 'GRCh38':
        raise ValueError('Assembly must be GRCh38')


    # Build gene table
    gene_info = get_gene_info(gene)
    gene_df = pd.DataFrame([{'Organism': 'Human', 
                             'Gene': gene_info.get('name'), 
                             'Description': gene_info.get('description'), 
                             'HGNC': gene_info.get('HGNC'), 
                             'Ensembl ID': gene_info.get('id'), 
                             'Biotype': gene_info.get('biotype'), 
                             'Strand': gene_info.get('strand'), 
                             'Chr': chrom, 
                             'Start': start, 
                             'End': end, 
                             'Ref': ref, 
                             'Alt': alt}])


    # Use VEP to get pathogenicity and protein info
    vep_df = get_vep_data(chrom, start, end, alt)
    vep_df.to_csv('vep_test.csv', index=False)
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
    for domain in vep_df['domains'].unique():
        dom_name = get_domain_name(domain)
        vep_df.loc[vep_df['domains'] == domain, 'domain_name'] = dom_name

    
    # Create protein DataFrame
    protein_df = pd.DataFrame()
    old_cols = ['transcript_id', 'biotype', 'exon', 'domains', 'domain_name', 'polyphen_prediction', 
                'polyphen_score', 'consequence_terms', 'codons', 'amino_acids', 'refAA', 'varAA'] # Column names in vep_df
    new_cols = ['Transcript ID', 'Biotype', 'Exon Rank', 'Pfam Domain ID', 'Pfam Domain Name', 'Polyphen Prediction', 
                'Polyphen Score', 'Molecular Consequence', 'Codon Switch', 'Amino Acids', 'refAA', 'varAA'] # Corresponding column names for protein_df
    protein_df[new_cols] = vep_df[old_cols] # Update column names
    

    # Disease Associations
    clinvar_ids = get_clinvar_ids(gene, start)

    # Build lookup dictionary of protein cahnge and consequence to diseases
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

    # Save resulting tables to CSVs
    gene_df.to_csv('hum_gene_df.csv', index=False)
    protein_df.to_csv('hum_protein_df.csv', index=False)
    
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

    # Build orthologous gene table
    mouse_gene_info = get_gene_info(gene, species='mouse')
    mouse_gene_df = pd.DataFrame([{'Organism': 'Mouse', 
                                   'Gene': mouse_gene_info.get('name'), 
                                   'Description': mouse_gene_info.get('description'), 
                                   'HGNC': mouse_gene_info.get('HGNC'), 
                                   'Ensembl ID': mouse_gene_info.get('id'), 
                                   'Biotype': mouse_gene_info.get('biotype'), 
                                   'Strand': mouse_gene_info.get('strand')}])


    # Only retain variants with logged genomic infromation
    mouse_select = mouse_vars[(mouse_vars['AlleleAssociatedGeneId'] == mouse_gene_info.get('HGNC')) &
                              (mouse_vars['StartPosition'].notnull())].copy()

    # List of alleles associated without genomic information
    mouse_allele_extra = mouse_vars[(mouse_vars['AlleleAssociatedGeneId'] == mouse_gene_info.get('HGNC')) &
                                    (mouse_vars['StartPosition'].isnull())]['AlleleId'].tolist()


    # Create allele DataFrame
    mouse_allele_df = pd.DataFrame()
    old_cols = ['AlleleAssociatedGeneSymbol', 'AlleleId', 'Chromosome', 'StartPosition', 
                'EndPosition', 'SequenceOfReference', 'SequenceOfVariant'] # Column names in mouse_select
    new_cols = ['Gene', 'AlleleId', 'Chromosome', 'Start', 'End', 'Ref', 'Alt'] # Corresponding column names for mouse_allele_df
    mouse_allele_df[new_cols] = mouse_select[old_cols] # Update column names


    # Extract primary Molecular Consequence term
    mouse_allele_df['Molecular Consequence'] = mouse_select['MostSevereConsequenceName'].str.split(',').str[0]


    # Convert coordinates to integers
    mouse_allele_df['Start'] = mouse_allele_df['Start'].astype(int)
    mouse_allele_df['End'] = mouse_allele_df['End'].astype(int)


    # Query VEP for each variant
    mouse_vep_df = pd.DataFrame()
    for i, row in mouse_allele_df.iterrows():
        variant_vep = get_vep_data(row['Chromosome'], row['Start'], row['End'], row['Alt'], species='mouse')
        variant_vep['AlleleId'] = row['AlleleId']
        mouse_vep_df = pd.concat([mouse_vep_df, variant_vep], ignore_index=True)
    
    # Keep only protein coding transcript consequences
    mouse_prt_df = mouse_vep_df[mouse_vep_df['biotype'] == 'protein_coding'].copy()
    

    # Split amino_acids to REFAA and VARAA if in X/Y format
    aa_split = mouse_prt_df['amino_acids'].astype(str).str.split('/', expand=True)
    mouse_prt_df['refAA'] = aa_split[0]
    mouse_prt_df['varAA'] = aa_split[1]

    # Add protein location of amino_acids
    mouse_prt_df['protein_start'] = mouse_prt_df['protein_start'].astype(pd.Int64Dtype())
    mouse_prt_df['amino_acids'] = mouse_prt_df.apply(lambda row: f"{row['refAA']}{row['protein_start']}{row['varAA']}", axis=1)

    # Fix consequence formating - originally returned as a list
    mouse_prt_df['consequence_terms'] = mouse_prt_df['consequence_terms'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)


    # Fetch domain names
    for domain in mouse_prt_df['domains'].unique():
        dom_name = get_domain_name(domain)
        mouse_prt_df.loc[mouse_prt_df['domains'] == domain, 'Pfam Domain Name'] = dom_name


    # Rearrange and rename columns
    mouse_prt_df = mouse_prt_df[['AlleleId', 'transcript_id', 'biotype', 
                                 'exon', 'domains', 'Pfam Domain Name', 
                                 'consequence_terms', 'codons', 'amino_acids', 
                                 'refAA', 'varAA']]
    
    mouse_prt_df.rename(columns={'transcript_id': 'Transcript ID', 'biotype': 'Biotype', 
                                 'exon': 'Exon Rank', 'domains': 'Pfam Domain ID', 
                                 'consequence_terms': 'Molecular Consequence', 
                                 'codons': 'Codon Switch', 'amino_acids': 'Amino Acids'})


    # Query MouseMine for ontology associations
    mouse_ontology_df = mouse_disease_associations(gene)
    mouse_ontology_df.to_csv('mouse_onntology_df.csv', index=False)

    # Merge ontology associations into prt DataFrame
    mouse_prt_df = mouse_prt_df.merge(mouse_ontology_df[['primaryIdentifier', 'ontologyName', 'ontologyID']], 
                                      left_on='AlleleId', right_on='primaryIdentifier')


    # Print resulting tables
    print(mouse_gene_df)
    print('----------')
    print(mouse_allele_df)
    print('----------')
    print(mouse_prt_df)
    print('----------')
    print(mouse_allele_extra)

    # Save resulting tables to CSVs
    mouse_gene_df.to_csv('mouse_gene_df.csv', index=False)
    mouse_prt_df.to_csv('mouse_protein_df.csv', index=False)
    mouse_allele_df.to_csv('mouse_allele_df.csv', index=False)
    
    return mouse_gene_df, mouse_prt_df

def score_ortho_vars(human_prt_df, mouse_prt_df):
    '''
    Score orthologous mouse variant similarity against human variant.

    Parameter:
        human_prt_df: pandas.DataFrame. A single-row DataFrame describing the human variant.
        mouse_prt_df : pandas.DataFrame. A multi-row DataFrame of all orthologous mouse alleles to compare.

    Returns:
        pandas.DataFrame. Containing boolean match columns and total score for each mouse allele
    '''

    # Extract human variant features
    biotype = human_prt_df['Biotype'][0]
    exon_rank = human_prt_df['Exon Rank'][0]
    domain = human_prt_df['Pfam Domain ID'][0]
    consequence = human_prt_df['Molecular Consequence'][0]
    amino_acids = human_prt_df['Amino Acids'][0]
    diseases = human_prt_df['Associated Diseases'][0]

    # Convert human diseases to DOID names
    diseases_doid = []
    # Go through each correlated disease
    for disease in diseases.split('; '):
        doid = get_doid_name(disease)
        if doid:
            diseases_doid.append(doid['label']) # append disease name if found

    # Score each mouse variant against human
    score_rows = []
    for idx, mouse_prt in mouse_prt_df.iterrows():
        match_biotype = mouse_prt['biotype'] == biotype
        match_exon = mouse_prt['exon'] == exon_rank
        match_domain = mouse_prt['domains'] == domain
        match_consequence = mouse_prt['consequence_terms'] == consequence
        match_AA = mouse_prt['amino_acids'] == amino_acids
        match_disease = mouse_prt['ontologyName'] in diseases_doid

        # Calculate total score (precentage of passed checks)
        total_score = sum([match_biotype,
                           match_exon,
                           match_domain,
                           match_consequence,
                           match_AA,
                           match_disease]) / 6 * 100

        # Append score row
        score_rows.append({'biotype_match': match_biotype,
                           'exon_match': match_exon,
                           'domain_match': match_domain,
                           'consequence_match': match_consequence,
                           'AA_match': match_AA,
                           'disease_match': match_disease,
                           'total_score': total_score})

    # Create and display score DataFrame
    score_df = pd.DataFrame(score_rows)
    score_df.to_csv('var_scores_df.csv', index=False)

    print(score_df)

    return score_df

if __name__ == '__main__':

    gene = 'ACVR1'
    start = 157774114
    end = 157774114
    chrom = 2
    ref = 'C'
    alt = 'T'
    assembly = 'GRCh38'