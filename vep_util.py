import yaml
import logging
import requests
import numpy as np
import pandas as pd
import subprocess
import json

# Load Config
with open('config.yaml') as f:
    config = yaml.safe_load(f)

# Load APIs
ensembl_base_url = config['api']['ensembl_base_url']
interpro_base_url = config['api']['interpro_base_url']


def get_domain_name(pfam_domain_id):
    '''
    Query InterPro Pfam database to get domain name.

    Parameters:
        pfam_domain_id: string. Pfam domain ID to extract full name for.

    Returns:
        string. The corresponding long name for given domain ID.
    '''
    # Construct URL
    url = interpro_base_url + pfam_domain_id
    headers = {'Accept': 'application/json'}

    # Log InterPro request
    logging.info(f'InterPro Request: {url}')

    try:
        # Make request to InterPro REST API
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        json = r.json()

        # Extract domain long name
        metadata = json['metadata']
        long_name = metadata['name'].get('name') or metadata['name'].get('value')

        # Return domain long name
        return long_name
    
    # Return None if the request fails
    except Exception:
        logging.error(f'InterPro Request Failed: {url}')
        return None


def get_vep_data(chromosome, start, end, alt, species='human'):
    '''
    Use Ensembl VEP REST endpoint to get consequence data.

    Parameters:
        chrom: str. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        alt: string. Alternate allele
        species: string. Species to extract gene from. Defaults to 'human'.

    Returns:
        DataFrame.
    '''
    # Build region string | chromosome:start-stop/alt
    region = f'{chromosome}:{start}-{end}/{alt}'

    # Construct URL
    url = f'{ensembl_base_url}/vep/{species}/region/{region}?numbers=1&domains=1'
    
    # Only choose one reference transcript for human variants
    # Example url: https://rest.ensembl.org/vep/human/region/2:157774114-157774114/T?numbers=1&domains=1&pick_order=mane_select,length&pick=1
    if species == 'human': 
        url += '&pick_order=mane_select,length&pick=1'

    # Log VEP request
    logging.info(f'Ensembl VEP Request: {url}')

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
            logging.warning(f'VEP yielded no results')
            return pd.DataFrame()
        
        # Normalize JSON result into a DataFrame
        vep_df = pd.json_normalize(consequences)

        # Keep desired columns that exist in DataFrame
        keep = ['transcript_id', 'polyphen_prediction', 'polyphen_score', 'amino_acids',
                'protein_start', 'protein_end', 'consequence_terms', 'exon', 'domains', 
                'codons', 'impact', 'biotype']
        vep_df = vep_df[[col for col in keep if col in vep_df.columns]]
        
        # Extract pfam names from dictionary of domain names and databases
        pfam = None

        if 'domains' in vep_df.columns: 
            for domain in vep_df['domains'][0]:
                if isinstance(domain, dict) and domain.get('db') == 'Pfam':
                    pfam = domain.get('name')
                    break
            
            if pfam is None:
                logging.warning('VEP domain failure: no pfam database in response')

        else:
            logging.warning('VEP domain failure: no "domains" information in response')

        vep_df['domains'] = pfam
        vep_df.rename(columns={'domains': 'Domain'}, inplace=True)

        # Handle empty or null amino_acid values - frequently asssociated with downstream_gene_variants and upstream_gene_variants
        vep_df['amino_acids'] = vep_df.apply(lambda row: row['amino_acids'] if pd.notnull(row['amino_acids']) else None, axis=1)

        # Return VEP DataFrame
        return vep_df
    
    # Handle request failure
    except Exception as e:
        logging.error('VEP request failed:', e)
        return pd.DataFrame()


def extract_pfam_name(domain_list):
    """Extract first Pfam domain name from a VEP domains list."""
    if not isinstance(domain_list, list):
        return None
    for d in domain_list:
        if isinstance(d, dict) and d.get('db') == 'Pfam':
            return d.get('name')
    return None


def fetch_vep_data(variants, species):
    '''
    Use Ensembl VEP REST endpoint to get consequence data.

    Parameters:
        variants: list of strings. List of variant HGVS notations to search.
        species: string. Species to extract gene from. Defaults to 'human'.

    Returns:
        DataFrame. Containing gene_symbol, transcript_id, polyphen_prediction, polyphen_score, 
        amino_acids, protein_start, consequence_terms, exon, codons, impact, biotype, Domain, HGVS.
    '''
    # Construct URL extension
    ext = f'/vep/{species}/hgvs?numbers=1&domains=1'

    # Only select MANE transcript if human
    if species == 'human':
        ext += '&pick_order=mane_select,length&pick=1'

    headers={'Content-Type': 'application/json', 'Accept': 'application/json'}
    hgvs_dict = {'hgvs_notations': list(set(variants))}

    logging.info(f'Ensembl VEP Request: {ensembl_base_url + ext}, for {hgvs_dict}')

    try:
        # Make POST request
        r = requests.post(ensembl_base_url + ext, headers=headers, json=hgvs_dict, timeout=30)
        r.raise_for_status()
        data = r.json()

        if not data:
            logging.warning("No VEP data returned")
            return pd.DataFrame()

        # Process each variant
        all_dfs = []
        for variant in data:
            consequences = variant.get('transcript_consequences', [])
            if not consequences:
                logging.warning(f"No transcript consequences for {variant.get('input')}")
                continue

            # Extract desired columns
            df = pd.json_normalize(consequences)
            keep = ['gene_symbol', 'transcript_id', 'polyphen_prediction', 'polyphen_score',
                    'amino_acids', 'protein_start', 'consequence_terms',
                    'exon', 'domains', 'codons', 'impact', 'biotype']
            df = df[[col for col in keep if col in df.columns]]

            # Extract Pfam domain names
            if 'domains' in df.columns:
                df['Domain'] = df['domains'].apply(extract_pfam_name)
                df.drop('domains', axis=1, inplace=True)
            else:
                logging.warning('VEP domain failure: no "domains" information in response')
                df['Domain'] = None

            # Add variant identifier for traceability
            df['HGVS'] = variant.get('input')

            all_dfs.append(df)

        # Combine all variant DataFrames
        vep_data = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()

        vep_data = vep_data.replace({np.nan: None})

        vep_data.to_csv('./testing_results/vep_output.csv', index=False)

        return vep_data
    
    # Handle request failure
    except Exception as e:
        logging.error('VEP request failed:', e)
        return pd.DataFrame()


def prepare_vep_output(vep_df, species):
    '''
    Cleans VEP output for desired protein information and format.

    Inputs:
        vep_df: DataFrame. Vep output from fetch_vep_data function.
    
    Returns:
        DataFrame. Containing Transcript ID, Gene Symbol, HGVS, Biotype, Exon Rank, Pfam Domain ID, Pfam Domain Name, Polyphen Prediction
        Polyphen Score, Molecular Consequence, Codon Switch, Amino Acids, refAA, varAA.
    '''
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
    keep = ['Input', 'Submission', 'gene_symbol', 'transcript_id', 'HGVS', 'biotype', 'exon', 'Domain', 'domain_name', 'polyphen_prediction', 
            'polyphen_score', 'consequence_terms', 'codons', 'amino_acids', 'refAA', 'varAA']
    cols_present = [col for col in keep if col in vep_df.columns]
    protein_df = vep_df[cols_present].copy()

    # Update column names
    protein_df.rename(columns={'transcript_id': 'Transcript ID', 'gene_symbol': 'Gene Symbol', 'biotype': 'Biotype', 
                       'exon': 'Exon Rank', 'Domain': 'Pfam Domain ID', 'domain_name': 'Pfam Domain Name', 
                       'polyphen_prediction': 'Polyphen Prediction', 'polyphen_score': 'Polyphen Score', 
                       'consequence_terms': 'Molecular Consequence', 'codons': 'Codon Switch', 'amino_acids': 'Amino Acids'}, 
                       inplace=True,  errors='ignore')
    
    return protein_df


def docker_input(variants):
    '''
    Creates a VCFv4.0 of unique variants in given DataFrame.

    Parameters:
        variants: DataFrame. Contains Chromosome, Start, Ref, Alt
    
    Creates file processing/vep/input.vcf
    '''
    unique_variants = variants.drop_duplicates(subset='Submission')

    with open('processing/vep/input.vcf', 'w') as f:
        # VCF meta-information
        f.write('##fileformat=VCFv4.0\n')

        # Header line
        f.write('#CHROM\tPOS\tID\tREF\tALT\n')

        # Write variants
        for _, row in unique_variants.iterrows():
            f.write(f'{row['Chromosome']}\t'
                    f'{row['Start']}\t'
                    f'.\t'
                    f'{row['Ref']}\t'
                    f'{row['Alt']}\n')

def docker_vep(species, input_file='input.vcf', output_file='output.json'):
    '''
    Runs local VEP on input file for species.

    Parameters:
        species: str. Either 'homo_sapiens' or 'mus_musculus'.
        input_file: str. Input file in processing dictionary. Defaults to 'input.vcf'
        output_file: str. Output file in processing dictionary. Defaults to 'output.json'
    '''
    cmd = ['docker', 'run', '--rm',
           '-v', '/Users/szemac/Desktop/MGI_Variant_Mapping/processing/vep:/processing',
           '-v', '/Users/szemac/.vep:/opt/vep/.vep',
           'ensemblorg/ensembl-vep',
           'vep',
           '--input_file', f'/processing/{input_file}',
           '--output_file', f'/processing/{output_file}',
           '--json',
           '--cache',
           '--offline',
           '--force_overwrite',
           '--species', species,
           '--domains',
           '--biotype',
           '--symbol',
           '--numbers']
    
    if species == 'homo_sapiens':
        cmd += ['--polyphen', 'b',
                '--pick',
                '--pick_order', 'mane_select,length']

    subprocess.run(cmd, check=True)

def parse_vep_json(input_file='processing/vep/output.json'):
    '''
    Parse VEP json output (Docker output) into a DataFrame.

    Parameters:
        inpute_file. str. path to json for parsing.
    
    Returns:
        DataFrame. VEP Output.
    '''

    all_dfs = []

    with open(input_file) as f:
        for line in f:  # NDJSON: one variant per line
            if not line.strip():
                continue

            variant = json.loads(line)

            consequences = variant.get('transcript_consequences', [])
            if not consequences:
                logging.warning(f"No transcript consequences for {variant.get('input')}")
                continue

            # Flatten transcript consequences
            df = pd.json_normalize(consequences)

            # Keep only columns you care about
            keep = [
                'gene_symbol', 'transcript_id',
                'polyphen_prediction', 'polyphen_score',
                'amino_acids', 'protein_start',
                'consequence_terms', 'exon',
                'domains', 'codons',
                'impact', 'biotype',
                'hgvsp', 'hgvsc'
            ]
            df = df[[c for c in keep if c in df.columns]]

            # Extract Pfam domain
            if 'domains' in df.columns:
                df['Domain'] = df['domains'].apply(extract_pfam_name)
                df.drop(columns=['domains'], inplace=True)
            else:
                df['Domain'] = None

            # Attach submitted input (optional)
            df['Submission'] = variant.get('input')

            all_dfs.append(df)

    vep_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
    vep_df = vep_df.replace({np.nan: None})

    return vep_df

def run_vep(variants, species):

    docker_input(variants)

    docker_vep(species)
    
    vep_df = parse_vep_json()

    return vep_df

