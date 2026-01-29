import yaml
import json
import logging
import pandas as pd
from Bio import Entrez
from cyvcf2 import VCF
import xml.etree.ElementTree as ET


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


clinvar_vcf = VCF(config['paths']['clinvar_vcf'])


MGI_disease_df = pd.read_csv(config['paths']['mgi_disease'],
                             sep='\t', comment='#', header=None,
                             names=['DOterm', 'DOID', 'NOTmodel', 'AllelePairs', 'StrainBackground',
                                    'AlleleSymbol', 'AlleleID', 'NumReferences', 'AlleleRepositoryID',
                                    'AlleleRRID', 'MarkerSymbol', 'MarkerMGIid', 'GeneRepositoryID', ''])


mondo_xref_map = json.load(open('app/data/MONDO/mondo_xref_map.json'))
mondo_term_map = json.load(open('app/data/MONDO/mondo_term_map.json'))


disease_mapping = json.load(open(config['paths']['doid_map']))


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


def assign_clinvar(variants):
    '''
    Fetch ClinVar disease associations for a list of variants.

    Parameters:
        variants: DataFrame. Expects Input, Chromosome, STart, Stop, Ref, Alt.
    
    Returns:
        DataFrame. Contains Input and CLINDISDB.
    '''
    results = []

    # Precompute HGVS list
    hgvs_list = variants['Input'].tolist()

    # Cache to avoid redundant queries
    cache = {}

    for idx, row in variants.iterrows():
        chrom, start, stop, ref, alt = row.Chromosome, row.Start, row.Stop, row.Ref, row.Alt
        key = (chrom, start, stop, ref, alt)

        if key not in cache:
            cache[key] = search_clinvar(chrom, start, stop, ref, alt)

        results.append({'Input': hgvs_list[idx], 'CLNDISDB': cache[key]})

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
                disease_info = record.INFO.get('CLNDISDB')

                if disease_info: 
                    diseases = [d.split(',') for d in disease_info.split('|')]


    except Exception as e:
        logging.warning(f'ClinVar search failed for {region}: {e}')
            
    return diseases


def fetch_mus_doid(allele_ids):
    '''
    Fetches all Mouse disease associations with human models for all alleles from list of MGI_gene_ids.
    
    Parameters:
        MGI_gene_ids: list of str. A list of all MGI_gene_ids to return strings for.
    
    Returns:
        dictionary. Mapping of AlleleID to set of DOIDs.
    '''

    filtered_disease_df = MGI_disease_df[MGI_disease_df['AlleleID'].isin(allele_ids)]

    doid_allele_dict  = (filtered_disease_df.groupby('AlleleID')['DOID'].apply(set).to_dict())
        
    return doid_allele_dict

def map_doids_to_mondo(doid_map):
    '''
    Convert {AlleleID: {DOID}} to a long table, merge MONDO info, then collapse back into per-allele lists.

    Parameters:
        doid_map: dict. Mapping of AlleleID to DOID.

    Returns:
        tuple of dicts. Mapping of AlleleID to MONDO ID and mapping of AlleleID to Disease Association terms.
    '''
    # Normalize to long df
    rows = [{'AlleleID': allele, 'DOID': doid}
            for allele, doids in doid_map.items()
            for doid in doids]
    
    df = pd.DataFrame(rows)

    # Map MONDO ID and MONDO term
    df['MONDO'] = df['DOID'].map(mondo_xref_map)
    df['Disease Association'] = df['MONDO'].map(mondo_term_map)

    # Collapse back to allele level
    mondo_by_allele = df.groupby('AlleleID')['MONDO'].apply(lambda s: sorted(set(s.dropna())))
    term_by_allele  = df.groupby('AlleleID')['Disease Association'].apply(lambda s: sorted(set(s.dropna())))

    return mondo_by_allele.to_dict(), term_by_allele.to_dict()


def clndisdb_to_mondo(clinvar_results):
    '''
    Extracts MONDO IDs and terms from ClinVar CLNDISDB entries.

    Parameters:
        clinvar_results: DataFrame. containing 'CLNDISDB' column with disease entries.
    
    Returns:
        DataFrame. with added 'MONDO' and 'Associated Diseases' columns.
    '''
    mondos_col = []
    terms_col = []

    # Iterate rows
    for _, disease_list in enumerate(clinvar_results['CLNDISDB'].values):

        mondos = []
        terms = []

        # Iterate diseases
        for disease_tokens in disease_list:

            # Extract MONDO and term
            mondo, term = resolve_disease_tokens(disease_tokens)

            if mondo:
                mondos.append(mondo)
                terms.append(term)

        mondos_col.append(mondos)
        terms_col.append(terms)

    # Assign results to HGVS DataFrame
    clinvar_results['Associated Diseases'] = terms_col
    clinvar_results['MONDO'] = mondos_col

    clinvar_results.drop('CLNDISDB', inplace=True, axis=1)

    return clinvar_results

def resolve_disease_tokens(tokens):
    '''
    Given a list of tokens from a CLNDISDB disease entry, extracts MONDO term.

    Parameters:
        tokens: list of str. Tokens from a CLNDISDB entry (disease IDs).

    Returns:
        tuple: MONDO ID (str) and MONDO term (str). If none found, (None, None).
    '''
    for token in tokens:

        # Direct MONDO entry
        if token.startswith("MONDO:"):
            _, mondo = token.split(":", 1)
            return mondo, mondo_term_map.get(mondo)

        # Normalize MedGen to UMLS
        if token.startswith("MedGen:"):
            token = "UMLS:" + token.split(":")[1]

        # Xref lookup
        if token in mondo_xref_map:
            mondo = mondo_xref_map[token]
            return mondo, mondo_term_map.get(mondo)

    # No valid MONDO mapping
    return None, None

