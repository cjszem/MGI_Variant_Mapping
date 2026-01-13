import re
import sys
import json
import logging
from xmlrpc import server
import requests
import pandas as pd
from cyvcf2 import VCF
from Bio import Entrez
import xml.etree.ElementTree as ET

clinvar_vcf = VCF('data/NCBI/clinvar.grch38.vcf.gz')


def batch_clinvar(variants):
    '''

    '''
    diseases = []

    for x, row in variants.iterrows():
        disease = search_clinvar(row['chrom'], row['start'], row['stop'], row['ref'], row['alt'])
        diseases.append(disease)
    
    return diseases

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
    print(region)
    diseases = None

    try:
        # Filter by region and search for variants
        for record in clinvar_vcf(region):
            if record.REF == ref and alt in record.ALT:

                # Extract disease names
                disease_info = record.INFO.get('CLNDN')

                diseases = disease_info.split('|') if disease_info else None
        
    except:
        logging.warning(f'ClinVar search failed for region {region} with ref {ref} and alt {alt}')
            
    return diseases

# def alt_clinvar():
#     records = []
#     for rec in clinvar_vcf:
#         records.append({
#             'chrom': rec.CHROM,
#             'start': rec.POS,
#             'ref': rec.REF,
#             'alt': rec.ALT,  # list of ALT alleles
#             'disease': rec.INFO.get('CLNDN')
#         })

#     vcf_df = pd.DataFrame(records)

#     merged = variants.merge(vcf_df, on=['chrom', 'start', 'ref', 'alt'], how='left')
#     print(merged)

if __name__ == '__main__':

    variants = pd.DataFrame({'chrom': ['2', '1', '11'], 
                             'start': [157774114, 39468726, 17394295], 
                             'stop': [157774114, 39468726, 17394295], 
                             'ref': ['C', 'T', 'C'], 
                             'alt': ['T', 'G', 'T']})
    


    r = search_clinvar('2', 157774114, 157774114, 'C', 'T')
    print(r)

    rese = batch_clinvar(variants)
    print(rese)