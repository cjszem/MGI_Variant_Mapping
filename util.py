'''
Cole Szeman
The Jackson Laboratory
2025
Utility functions for Mouse and Human Variant Annotation
'''

import re
import json
import logging
import requests
import pandas as pd
from cyvcf2 import VCF
from Bio import Entrez
import xml.etree.ElementTree as ET

# Setup API Information
Entrez.email = "cole.szeman@jax.org"
ensembl_base_url = 'https://rest.ensembl.org'
interpro_base_url = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/'
mousemine_base_url = 'https://www.mousemine.org/mousemine/service/query/results'
ols_base_url = 'https://www.ebi.ac.uk/ols4/api/'

# Setup logging
logging.basicConfig(filename='mapping.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Load ClinVar VCF data
clinvar_vcf = VCF('data/NCBI/clinvar.grch38.vcf.gz')

# Load MGI Allele data
mus_alleles_df = pd.read_csv('./data/MGI/VARIANT-ALLELE-MGI.tsv', sep='\t')

# Load MGI Gene Homology data
homology_df = pd.read_csv('./data/MGI/HOM_ProteinCoding.rpt', 
                          names=['MGI_ID', 'MusGeneSymbol', 'MusEntrezGeneID', 'MusHGNC_ID', 'HumGeneSymbol', 'HumEntrezGeneID'], 
                          sep='\t')
homology_dict = dict(zip(homology_df['HumGeneSymbol'], homology_df['MusGeneSymbol']))

# Load MGI Disease data
disease_cols = ['DOterm', 'DOID', 'NOTmodel', 'AllelePairs', 'StrainBackground', 'AlleleSymbol', 'AlleleID',
                'NumReferences', 'AlleleRepositoryID', 'AlleleRRID', 'MarkerSymbol', 'MarkerMGIid', 'GeneRepositoryID', '']
MGI_disease_df = pd.read_csv('./data/MGI/MGI_DiseaseMouseModel.rpt', names=disease_cols,
                             comment="#", header=None, sep='\t')

# Load DOID mappings
disease_mapping = json.load(open('data/DOID/doid_map.json'))

# Load Ensembl gene data
mouse_gene_pqt = pd.read_parquet('data/Ensembl/mouse_gene_data.parquet')
human_gene_pqt = pd.read_parquet('data/Ensembl/human_gene_data.parquet')

# ---- Helper functions ----
def log_query(gene, chrom, start, end, ref, alt):
    '''
    Logs the query parameters in mapping.log file for debugging and record keeping purposes.
    '''
    logging.info(f'Query Submitted: Gene={gene}, Chromosome={chrom}, Start={start}, End={end}, Ref={ref}, Alt={alt}')

def log_batch_query(variants):
    '''
    Logs the query in mapping.log file for debugging and record keeping purposes.
    '''
    logging.info(f'Query Submitted: processing {len(variants)} variants')

def process_batch_query(input):
    '''
    Processes a string resulting from the batch input of variants with formatting:
    chromosome:start-end:ref/allele

    Parameters:
        input: string. Multiline string containg all variants to process.

    Returns: DataFrame. Containing chrom, start, stop, ref, alt, for each variant.
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
            variants.append({'chrom': match.group(1),
                             'start': int(match.group(2)),
                             'stop': int(match.group(3)),
                             'ref': match.group(4),
                             'alt': match.group(5)})
        
    # Create pandas dataframe
    variants = pd.DataFrame(variants)
    variants.to_csv('./testing_results/batch_variants.csv', index=False)

    return variants

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

    # Log Ensembl request
    logging.info(f'Ensemble Gene Info Request: {url}')

    try:
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
        
        # Different MGI IDs for mouse and HGNC IDs for human
        if species == 'human':
            hgnc_id = match.group(3).strip()

            # Create gene dictionary
            gene_info = {'HGNC_id': hgnc_id, 'ENSMBL_id': ens_id, 'name': name, 'strand': strand, 'biotype': biotype, 'description': description}
        
        if species == 'mouse':
            hgnc_id = match.group(3).strip()

            # Create gene dictionary
            gene_info = {'MGI_id': hgnc_id, 'ENSMBL_id': ens_id, 'name': name, 'strand': strand, 'biotype': biotype, 'description': description}
        
    except Exception as e:
        logging.error(f'Ensembl Gene Info Request Failed: {e}')
        ValueError(f'Ensembl Gene Info Request Failed: {e}')

    return gene_info

def fetch_gene_info(genes, species='human'):
    '''
    Fetch gene info from saved parquets in /data/Ensembl/ for a lsit of gene symbols.

    Parameters:
        genes: list of strings. Gene symbols to extract info for.
        species: string. 'human' or 'mouse'. Defaults to 'human'.

    Returns:
        DataFrame. Containins HGNC ID, Ensembl ID, name, strand, biotype, and description.
    '''
    # Assign correct parquet based on species
    if species == 'human': pqt = human_gene_pqt
    elif species == 'mouse': pqt = mouse_gene_pqt

    # Subset parquet to requested genes
    subset = pqt[pqt['Name'].isin(genes)]

    # Preserve order of input genes
    subset = subset.set_index('Name').reindex(genes).reset_index()

    # Extract description and accession
    subset[['Description', 'Accession']] = subset['description'].str.extract(r'^(.*?) \[Source:.+ Symbol.+Acc:(.+)\]$')

    # Clean DataFrame
    subset.drop('description', axis=1, inplace=True)
    subset.rename(columns={'Name': 'Gene', 'biotype': 'Biotype', 'seqid': 'Chromosome', 'start': 'Start', 
                           'end': 'End', 'strand': 'Strand', 'gene_id': 'Ensembl_ID'}, inplace=True)
    subset = subset[['Gene', 'Description', 'Biotype', 'Chromosome', 'Start', 'End', 'Strand', 'Ensembl_ID', 'Accession']]

    return subset

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

def prepare_hgvs(variants):
    '''
    Converts a DataFrame of variants into a list of HGHVS notations for VEP querying
    
    Parameters:
        variants: DataFrame. containing chrom, start, ref, alt.

    Returns:
        list of strings. HGVS notations (chrom:g.startRef>Alt).
    '''
    hgvs_list = (variants['chrom'].astype(str) + ':g.' + variants['start'].astype(str) + variants['ref'] + '>' + variants['alt']).tolist()

    return hgvs_list

def fetch_vep_data(variants, species):
    '''
    Use Ensembl VEP REST endpoint to get consequence data.

    Parameters:
        chrom: str. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        alt: string. Alternate allele
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
    hgvs_dict = {'hgvs_notations': variants}

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
            pfam = None
            if 'domains' in df.columns and len(df['domains']) > 0:
                for domain in df['domains'][0]:
                    if isinstance(domain, dict) and domain.get('db') == 'Pfam':
                        pfam = domain.get('name')
                        break

                if pfam is None: logging.warning('VEP domain failure: no pfam database in response')

            else: logging.warning('VEP domain failure: no "domains" information in response')
            
            # Add Pfam domain and drop domain dictionary
            df['Domain'] = pfam
            df.drop('domains', axis=1, inplace=True)

            # Add variant identifier for traceability
            df['HGVS'] = variant.get('input')

            all_dfs.append(df)

        # Combine all variant DataFrames
        vep_data = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()

        vep_data.to_csv('./testing_results/vep_output.csv', index=False)

        return vep_data
    
    # Handle request failure
    except Exception as e:
        logging.error('VEP request failed:', e)
        return pd.DataFrame()

def prepare_vep_output(vep_df):
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
    keep = ['gene_symbol', 'transcript_id', 'HGVS', 'biotype', 'exon', 'Domain', 'domain_name', 'polyphen_prediction', 
            'polyphen_score', 'consequence_terms', 'codons', 'amino_acids', 'refAA', 'varAA']
    protein_df = vep_df[keep]

    # Update column names
    protein_df.rename(columns={'transcript_id': 'Transcript ID', 'gene_symbol': 'Gene Symbol', 'biotype': 'Biotype', 
                       'exon': 'Exon Rank', 'Domain': 'Pfam Domain ID', 'domain_name': 'Pfam Domain Name', 
                       'polyphen_prediction': 'Polyphen Prediction', 'polyphen_score': 'Polyphen Score', 
                       'consequence_terms': 'Molecular Consequence', 'codons': 'Codon Switch', 'amino_acids': 'Amino Acids'}, inplace=True)
    
    return protein_df

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

def get_clinvar_ids(gene, start):
    '''
    Use Entrez to fetch ClinVar IDs for a given gene and location.

    Parameters:
        gene: string. Gene symbol to extract info for.
        start: int. Genomic location of variant start.

    Returns:
        list of strings. The ClinVar IDs of variants at given location.
    '''
    # Log Request to ClinVar
    logging.info(f'ClinVar ID Request: [Gene Name]={gene}, [Base Position]={start}')

    try:
        # Make request to ClinVar
        handle = Entrez.esearch(db='clinvar', term=f'{gene}[Gene Name] AND {start}[Base Position]', retmax=500) 
        record = Entrez.read(handle)
        handle.close()

        # Extract ClinVar IDs
        clinvar_ids = record['IdList']
    
    # Handle Request Failure
    except Exception as e:
        logging.error(f'ClinVar Request Failed: {e}')
        ValueError(f'ClinVar Request Failed: {e}')

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
    # Log Request to ClinVar
    logging.info(f'ClinVar ID Request: [ClinVar ID]={clinvar_id}')
    
    try:
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
    
    # Handle Request Failure
    except Exception as e:
        logging.error(f'ClinVar Disease Request Failed: {e}')
        ValueError(f'ClinVar Disease Request Failed: {e}')
    
    return prt_changes, consequence, diseases

def fetch_assign_clinvar(variants):
    '''
    Fetches ClinVar disease associations and assigns them to a list of variants.

    Inputs:
        Variants: DataFrame. Containing chrom, start, stop, ref, alt.

    Returns:
        DataFrame. Containing HGVS, diseases for each variant.
    '''
    hgvs = prepare_hgvs(variants)
    results = []

    for x, row in variants.iterrows():
        chrom, start, stop, ref, alt = row[['chrom', 'start', 'stop', 'ref', 'alt']]
        diseases = search_clinvar(chrom, start, stop, ref, alt)
        results.append({'HGVS': hgvs[x],
                        'diseases': diseases})

    return pd.DataFrame(results)

def search_clinvar(chrom, start, end, ref, alt):
    '''
    Searches ClinVar VCF for disease associations with a given variant.

    Parameters:
        chrom: str. Chromosome of variant.
        start: int. Genomic location of variant start.
        end: int. Genomic location of variant end.
        ref: string. Reference allele.
        alt: string. Alternate allele.
    
    Returns:
        list. containing: chromosome, start, end, ref, alt, diseases
    '''
    # Craft region string
    region = f'{chrom}:{start}-{end}'
    diseases = None

    try:
        # Filter by region and search for variants
        for record in clinvar_vcf(region):
            if record.REF == ref and alt in record.ALT:

                # Extract disease names
                disease_info = record.INFO.get('CLNDN')

                diseases = disease_info.split('|') if disease_info else None

                diseases = [d.replace('_', ' ').lower() for d in diseases]
        
    except:
        logging.warning(f'ClinVar search failed for region {region} with ref {ref} and alt {alt}')
            
    return diseases

def fetch_mus_alleles(MGI_gene_ids):
    '''
    Fetches all Mouse alleles with human models for all genes from list of MGI_gene_ids.
    
    Parameters:
        MGI_gene_ids: list of str. A list of all MGI_gene_ids to retur strings for.
    
    Returns:
        DataFrame. Containing Gene, AlleleID, AlleleSymbol, Chromosome, Start, End, Ref, Alt, Molecular Consequence.
    '''
    # Only retain variants with logged genomic information
    mouse_select = mus_alleles_df[(mus_alleles_df['AlleleAssociatedGeneId'].isin(MGI_gene_ids)) &
                                  (mus_alleles_df['StartPosition'].notnull())]

    # Create allele DataFrame
    mouse_allele_df = pd.DataFrame()
    old_cols = ['AlleleAssociatedGeneSymbol', 'AlleleId', 'AlleleSymbol', 'Chromosome', 'StartPosition', 
                'EndPosition', 'SequenceOfReference', 'SequenceOfVariant'] # Column names in mouse_select
    new_cols = ['Gene', 'AlleleID', 'AlleleSymbol', 'Chromosome', 'Start', 'End', 'Ref', 'Alt'] # Corresponding column names for mouse_allele_df
    mouse_allele_df[new_cols] = mouse_select[old_cols] # Update column names


    # Extract primary Molecular Consequence term
    mouse_allele_df['Molecular Consequence'] = mouse_select['MostSevereConsequenceName'].str.split(',').str[0]

    # Convert coordinates to integers
    mouse_allele_df['Start'] = mouse_allele_df['Start'].astype(int)
    mouse_allele_df['End'] = mouse_allele_df['End'].astype(int)

    return mouse_allele_df

def fetch_mus_disease_map(MGI_allele_ids):
    '''
    Fetches all Mouse disease associations with human models for all alleles from list of MGI_gene_ids.
    
    Parameters:
        MGI_gene_ids: list of str. A list of all MGI_gene_ids to return strings for.
    
    Returns:
        dictionary. Mapping of AlleleID to set of DO terms.
    '''

    filtered_disease_df = MGI_disease_df[MGI_disease_df['AlleleID'].isin(MGI_allele_ids)]

    mgi_disease_dict  = (filtered_disease_df.groupby('AlleleID')['DOterm'].apply(set).to_dict())
        
    return mgi_disease_dict

def mouse_disease_associations(gene):
    '''
    Query MouseMine for allele and ontology associations for a given mouse gene symbol.

    Parameters:
        gene: string. Gene symbol to extract info for.
    
    Returns:
        DataFrame. Containins symbol, name, primaryIdentifier, ontologyName, ontologyID.
    '''

    # log Request to MouseMine
    logging.info(f'MouseMine Request: [Gene Symbol]={gene}')
    
    # Construct query XML
    query_xml = f'''
    <query name="HMM Search" model="genomic" 
        view="Gene.alleles.symbol 
            Gene.alleles.name 
            Gene.alleles.primaryIdentifier 
            Gene.alleles.ontologyAnnotations.ontologyTerm.name 
            Gene.alleles.ontologyAnnotations.ontologyTerm.identifier" 
        longDescription="" 
        sortOrder="Gene.alleles.symbol asc" 
        constraintLogic="A">
        <constraint path="Gene.symbol" code="A" op="=" value="{gene}"/>
    </query>
    '''

    try:
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
                                    'Gene > Alleles > Ontology Annotations > Term Name':'ontologyName', 
                                    'Gene > Alleles > Ontology Annotations > Ontology Term . Identifier':'ontologyID'}, inplace=True)

        # Filter DataFrame to only retain DOID ontology terms
        ontology_df = ontology_df[ontology_df['ontologyID'].str.startswith('DOID:')]

    # Handle Request Failure
    except Exception as e:
        logging.error(f'MouseMine Request Failed: {e}')
        ValueError(f'MouseMine Request Failed: {e}')

    return ontology_df

def get_doid_name(disease):
    '''
    Convert disease name to DOID using OLS API.

    Parameters:
        disease: string. Disease name to convert.

    Returns:
        dictionary. Mapping label to ID of disease. Returns None if no match found.
    '''
    # Construct URL
    url = ols_base_url + f'search?q={disease}&ontology=doid&exact=true'

    # Log Request to OLS
    logging.info(f'OLS Request: {url}')

    try:
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

    except Exception as e:
        logging.error(f'OLS Request Failed: {e}')
        ValueError(f'OLS Request Failed: {e}')

    return {'label': doc.get('label'),
            'DOID': doc.get('obo_id')}

def fetch_doid_name(disease):
    '''
    Search disease name synonyms for DO name.

    Parameters:
        disease: string. Disease name to convert.

    Returns:
        string. DO disease name. Returns None if None found.
    '''
    doid_name = disease_mapping.get(disease.lower())

    return doid_name

def log_results(hum_gene_df, hum_protein_df, mouse_gene_df, mouse_prt_df, var_scores_df):
    '''
    Logs the resulting DataFrames to CSV files for record keeping and debugging purposes.
    '''
    # Save human tables to CSVs
    hum_gene_df.to_csv('./results/hum_gene_df.csv', index=False)
    logging.info('Saved hum_gene_df to ./results/hum_gene_df.csv')

    hum_protein_df.to_csv('./results/hum_protein_df.csv', index=False)
    logging.info('Saved hum_protein_df to ./results/hum_protein_df.csv')

    # Save mouse tables to CSVs
    mouse_gene_df.to_csv('./results/mouse_gene_df.csv', index=False)
    logging.info('Saved mouse_gene_df to ./results/mouse_gene_df.csv')

    mouse_prt_df.to_csv('./results/mouse_protein_df.csv', index=False)
    logging.info('Saved mouse_protein_df to ./results/mouse_protein_df.csv')

    # Save score table to CSV
    var_scores_df.to_csv('./results/var_scores_df.csv', index=False)
    logging.info('Saved var_scores_df to ./results/var_scores_df.csv')

# ---- Draft Main Functions ----
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


    # Query VEP for each variant
    mouse_vep_df = pd.DataFrame()
    for i, row in mouse_allele_df.iterrows():
        variant_vep = get_vep_data(row['Chromosome'], row['Start'], row['End'], row['Alt'], species='mouse')
        variant_vep['AlleleID'] = row['AlleleID']
        variant_vep['AlleleSymbol'] = row['AlleleSymbol']
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


    # Query MouseMine for ontology associations
    mouse_disease_map = fetch_mus_disease_map(mouse_prt_df['AlleleID'].unique())

    # Drop alleles without disease associations
    mouse_prt_df = mouse_prt_df[mouse_prt_df['AlleleID'].isin(mouse_disease_map.keys())]\
    
    # Map diseases to mouse alleles
    mouse_prt_df['Disease Association'] = mouse_prt_df['AlleleID'].map(mouse_disease_map)
    mouse_prt_df['Disease Association'] = mouse_prt_df['Disease Association'].apply(lambda x: ','.join(x) if isinstance(x, set) else x)


    # Fetch domain names
    for domain in mouse_prt_df['Domain'].unique():
        if domain == None:
            dom_name = None
        else:
            dom_name = get_domain_name(domain)
        mouse_prt_df.loc[mouse_prt_df['Domain'] == domain, 'Pfam Domain Name'] = dom_name


    # Rearrange and rename columns
    mouse_prt_df = mouse_prt_df[['AlleleID', 'AlleleSymbol', 'transcript_id', 'biotype', 
                                 'exon', 'Domain', 'Pfam Domain Name', 'consequence_terms', 
                                 'codons', 'amino_acids', 'refAA', 'varAA', 'Disease Association']]
    
    mouse_prt_df.rename(columns={'transcript_id': 'Transcript ID', 'biotype': 'Biotype', 'exon': 'Exon Rank', 
                                 'Domain': 'Pfam Domain ID', 'consequence_terms': 'Molecular Consequence', 
                                 'codons': 'Codon Switch', 'amino_acids': 'Amino Acids'}, inplace=True)
    

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

# ---- Batch Main Functions ----
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
    hgvs_list = prepare_hgvs(variants)

    # Query VEP
    vep_df = fetch_vep_data(hgvs_list, species='human')

    # Build gene table
    genes = vep_df['gene_symbol'].unique()
    gene_df = fetch_gene_info(genes, species='human')

    # Clean VEP output into protein DataFrame
    protein_df = prepare_vep_output(vep_df)

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

