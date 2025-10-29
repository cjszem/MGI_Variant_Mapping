import re
import requests
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET


Entrez.email = "cole.szeman@jax.org"
ensembl_base_url = 'https://rest.ensembl.org'
interpro_base_url = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/'
mouse_vars = pd.read_csv('VARIANT-ALLELE-MGI.tsv', sep='\t')


# ---- Helper functions ----
def get_gene_info(symbol, species='human'):
    '''
    Query Ensembl REST API and extract info for a given gene.
    '''
    url = f'{ensembl_base_url}/lookup/symbol/{species}/{symbol}?content-type=application/json'
    request = requests.get(url, headers={'Content-Type': 'application/json'})
    request.raise_for_status()
    data = request.json()

    ens_id = data.get('id')
    name = data.get('display_name')
    biotype = data.get('biotype')
    strand = data.get('strand')
    match = re.match(r'^(.*?) \[Source:(.+) Symbol;Acc:(.+)\]$', data.get('description'))
    description = match.group(1).strip()
    hgnc_id = match.group(3).strip()

    return {'HGNC': hgnc_id, 'id': ens_id, 'name': name, 'strand': strand, 'biotype': biotype, 'description': description}

def get_domain_name(pfscan_domain_id):
    '''
    Query InterPro Pfam database to get domain name.
    '''
    if not pfscan_domain_id:
        return None
    url = interpro_base_url + pfscan_domain_id
    headers = {'Accept': 'application/json'}

    try:
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        json = r.json()

        # Extract long name
        metadata = json['metadata']
        name = metadata['name'].get('name') or metadata['name'].get('value')

        return name
    
    except Exception:
        return None

def get_vep_data(chrom, start, end, alt, species='human'):
    '''
    Use Ensembl VEP REST endpoint to get consequence data.
    '''
    region = f'{chrom}:{start}:{end}/{alt}' # Build region string - chrom:pos/alt
    url = f'{ensembl_base_url}/vep/{species}/region/{region}?numbers=1&domains=1'
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    try:
        request = requests.get(url, headers=headers, timeout=15)
        request.raise_for_status()
        data = request.json()

        if not data: return pd.DataFrame()
        
        variant = data[0] # The first element contains variant information
        consequences = variant.get('transcript_consequences', [])

        if len(consequences) == 0: return pd.DataFrame()
        vep_df = pd.json_normalize(consequences)
        
        keep = ['transcript_id', 'polyphen_prediction', 'polyphen_score', 'amino_acids',
                'protein_start', 'protein_end', 'consequence_terms', 'exon', 'domains', 
                'codons', 'impact', 'biotype']
        vep_df = vep_df[[col for col in keep if col in vep_df.columns]]

        # Only retain pfam domain id
        def extract_pfam(domain_list):
            if not isinstance(domain_list, list): return '' # Skip if there are no domains
            pfam_names = [domain['name'] for domain in domain_list if isinstance(domain, dict) and domain.get('db') == 'Pfam'] # Extract pfam names into list
            return ';'.join(pfam_names)

        if 'domains' in vep_df.columns: vep_df['domains'] = vep_df['domains'].apply(extract_pfam) # Extract pfam names for all domains
        else: vep_df['domains'] = ''

        return vep_df
    
    except Exception as e:
        print('VEP request failed:', e)
        return pd.DataFrame()

def get_clinVar_ids(gene, start):
    '''
    Use Entrez to fetch ClinVar IDs for a given gene and location.
    '''
    handle = Entrez.esearch(db='clinvar', term=f'{gene}[Gene Name] AND {start}[Base Position]', retmax=500) 
    record = Entrez.read(handle)
    handle.close()
    clinvar_ids = record["IdList"]

    return clinvar_ids

def get_disease_associations(clinvar_id):
    handle = Entrez.efetch(db='clinvar', rettype='vcv', id=clinvar_id, from_esearch=True)
    xml_text = handle.read()
    handle.close()
    clinvar_root = ET.fromstring(xml_text)

    prt_changes = []
    for prt_section in clinvar_root.findall(".//ProteinChange"):
        prt = prt_section.text
        prt_changes.append(prt)

    for mol_cons in clinvar_root.findall(".//MolecularConsequence"):
        try:
            consequence = mol_cons.get("Type").replace(" ", "_")
            break
        except: continue

    # accepted_review_status = ['practice guideline', 'reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts']
    diseases = []
    for rcv in clinvar_root.findall(".//RCVAccession"):
        review = rcv.find(".//ReviewStatus").text
        # if review not in accepted_review_status:
        #     continue
        condition = rcv.find(".//ClassifiedCondition")
        disease = condition.text
        germline_description = rcv.find(".//Description")
        submissions = germline_description.get("SubmissionCount")
        diseases.append(f'{disease} ({submissions})')
    
    return prt_changes, consequence, diseases


# ---- Main function ----
def hvar_to_output(gene, chrom, start, end, ref, alt, assembly='GRCh38'):
    '''
      Find transcripts overlapping variant using pyensembl (local GTF cache)
      Call VEP REST to obtain per-transcript consequences
      Attempt fetches InterPro domain names
      Return two DataFrames: gene_df and prt_df
    '''
    if assembly != 'GRCh38':
        raise ValueError('Assembly must be GRCh38')

    # Build gene table
    gene_info = get_gene_info(gene)
    gene_df = pd.DataFrame([{'Organism': 'Human', 'Gene': gene_info.get('name'), 'Description': gene_info.get('description'), 'HGNC': gene_info.get('HGNC'), 'Ensembl ID': gene_info.get('id'), 
                             'Biotype': gene_info.get('biotype'), 'Strand': gene_info.get('strand'), 'Chr': chrom, 'Start': start, 'End': end, 'Ref': ref, 'Alt': alt}])


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

    vep_df['consequence_terms'] = vep_df['consequence_terms'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)

    # Fetch domain names
    for domain in vep_df['domains'].unique():
        dom_name = get_domain_name(domain)
        vep_df.loc[vep_df['domains'] == domain, 'domain_name'] = dom_name

    protein_df = pd.DataFrame()
    protein_df[['Transcript ID', 'Biotype', 'Exon Rank', 'Pfam Domain ID', 
                'Pfam Domain Name', 'Polyphen Prediction', 'Polyphen Score', 
                'Molecular Consequence', 'Codon Switch', 'Amino Acids', 'refAA', 'varAA']] = vep_df[['transcript_id', 'biotype', 'exon', 
                                                                                                     'domains', 'domain_name', 'polyphen_prediction', 'polyphen_score', 'consequence_terms', 'codons', 
                                                                                                     'amino_acids', 'refAA', 'varAA']]
    

    # Disease Associations
    disease_dict = {}
    clinvar_ids = get_clinVar_ids(gene, start)
    for clinvar in clinvar_ids:
        prt_changes, consequence, diseases = get_disease_associations(clinvar)
        for prt in prt_changes:
            disease_dict[f'{prt}_{consequence}'] = diseases

    associated_diseases = []
    for index, prt_row in protein_df.iterrows():
        key = f'{prt_row['Amino Acids']}_{prt_row['Molecular Consequence']}'
        d = disease_dict.get(key, None)
        if d is None: associated_diseases.append(None)
        else: associated_diseases.append('; '.join(d))

    protein_df['Associated Diseases (submissions)'] = associated_diseases

    print(gene_df)
    print('----------')
    print(protein_df)

    gene_df.to_csv('hum_gene_df.csv', index=False)
    protein_df.to_csv('hum_protein_df.csv', index=False)
    
    return gene_df, protein_df

def mvar_to_output(gene, assembly='GRCm39'):
    '''
      Find transcripts overlapping variant using pyensembl (local GTF cache)
      Call VEP REST to obtain per-transcript consequences
      Attempt fetches InterPro domain names
      Return two DataFrames: gene_df and prt_df
    '''
    if assembly != 'GRCm39':
        raise ValueError('Assembly must be GRCm39')

    # Build orthologous gene table
    mouse_gene_info = get_gene_info(gene, species='mouse')
    mouse_gene_df = pd.DataFrame([{'Organism': 'Mouse', 'Gene': mouse_gene_info.get('name'), 'Description': mouse_gene_info.get('description'), 'HGNC': mouse_gene_info.get('HGNC'), 
                                   'Ensembl ID': mouse_gene_info.get('id'), 'Biotype': mouse_gene_info.get('biotype'), 'Strand': mouse_gene_info.get('strand')}])

    mouse_select = mouse_vars[(mouse_vars['AlleleAssociatedGeneId'] == mouse_gene_info.get('HGNC')) &
                              (mouse_vars['StartPosition'].notnull())].copy()

    mouse_allele_extra = mouse_vars.loc[(mouse_vars['AlleleAssociatedGeneId'] == gene) &
                                        (mouse_vars['StartPosition'].isnull()), 'AlleleId']

    mouse_allele_df = pd.DataFrame()
    mouse_allele_df[['Gene', 'AlleleId', 'Chromosome', 'Start', 'End', 'Ref', 'Alt']] = mouse_select[['AlleleAssociatedGeneSymbol', 'AlleleId', 'Chromosome',
                                                                                                      'StartPosition', 'EndPosition', 'SequenceOfReference', 'SequenceOfVariant']].copy()

    mouse_allele_df['Molecular Consequence'] = mouse_select['MostSevereConsequenceName'].str.split(',').str[0]

    mouse_allele_df['Start'] = mouse_allele_df['Start'].astype(int)
    mouse_allele_df['End'] = mouse_allele_df['End'].astype(int)

    mouse_vep_df = pd.DataFrame()
    for i, row in mouse_allele_df.iterrows():
        variant_vep = get_vep_data(row['Chromosome'], row['Start'], row['End'], row['Alt'], species='mouse')
        variant_vep['AlleleId'] = row['AlleleId']
        mouse_vep_df = pd.concat([mouse_vep_df, variant_vep], ignore_index=True)
    
    mouse_prt_df = mouse_vep_df[mouse_vep_df['biotype'] == 'protein_coding'].copy()
    
    # Split amino_acids to REFAA and VARAA if in X/Y format
    aa_split = mouse_prt_df['amino_acids'].astype(str).str.split('/', expand=True)
    mouse_prt_df['refAA'] = aa_split[0]
    mouse_prt_df['varAA'] = aa_split[1]

    # Add protein location of amino_acids
    mouse_prt_df['protein_start'] = mouse_prt_df['protein_start'].astype(pd.Int64Dtype())
    mouse_prt_df['amino_acids'] = mouse_prt_df.apply(lambda row: f"{row['refAA']}{row['protein_start']}{row['varAA']}", axis=1)

    # Fix consequence formating -- returned as a list
    mouse_prt_df['consequence_terms'] = mouse_prt_df['consequence_terms'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)

    # Fetch domain names
    for domain in mouse_prt_df['domains'].unique():
        dom_name = get_domain_name(domain)
        mouse_prt_df.loc[mouse_prt_df['domains'] == domain, 'domain_name'] = dom_name

    print(mouse_gene_df)
    print('----------')
    print(mouse_allele_df)
    print('----------')
    print(mouse_prt_df)
    print('----------')
    print(mouse_allele_extra)

    mouse_gene_df.to_csv('mouse_gene_df.csv', index=False)
    mouse_prt_df.to_csv('mouse_protein_df.csv', index=False)
    mouse_allele_df.to_csv('mouse_allele_df.csv', index=False)
    
    return mouse_gene_df, mouse_prt_df


if __name__ == '__main__':

    gene = 'ACVR1'
    start = 157774114
    end = 157774114
    chrom = 2
    ref = 'C'
    alt = 'T'
    assembly = 'GRCh38'