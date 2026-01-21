import yaml
import json
import logging
import requests
import pandas as pd
from Bio import Entrez
from cyvcf2 import VCF
import xml.etree.ElementTree as ET


# Load Config
with open('config.yaml') as f:
    config = yaml.safe_load(f)


clinvar_vcf = VCF(config['paths']['clinvar_vcf'])


MGI_disease_df = pd.read_csv(config['paths']['mgi_disease'],
                             sep='\t', comment='#', header=None,
                             names=['DOterm', 'DOID', 'NOTmodel', 'AllelePairs', 'StrainBackground',
                                    'AlleleSymbol', 'AlleleID', 'NumReferences', 'AlleleRepositoryID',
                                    'AlleleRRID', 'MarkerSymbol', 'MarkerMGIid', 'GeneRepositoryID', ''])


mousemine_base_url = config['api']['mousemine_base_url']
ols_base_url = config['api']['ols_base_url']


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


def fetch_assign_clinvar(variants):
    '''
    Fetches ClinVar disease associations and assigns them to a list of variants.

    Inputs:
        Variants: DataFrame. Containing chrom, start, stop, ref, alt.

    Returns:
        DataFrame. Containing HGVS, diseases for each variant.
    '''
    hgvs = variants['HGVS'].tolist()
    results = []

    for x, row in variants.iterrows():
        chrom, start, stop, ref, alt = row[['Chromosome', 'Start', 'Stop', 'Ref', 'Alt']]
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

