import re
import pandas as pd


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
            variants.append({'Input': match.group(0),
                             'Chromosome': match.group(1),
                             'Start': int(match.group(2)),
                             'Stop': int(match.group(3)),
                             'Ref': match.group(4),
                             'Alt': match.group(5)})
        
    # Create pandas dataframe
    variants = pd.DataFrame(variants)

    return variants


def prepare_input(variants):
    '''
    Converts a DataFrame of variants into a list of HGHVS notations for VEP querying
    
    Parameters:
        variants: DataFrame. containing Chromosome, Start, Ref, Alt.

    Returns:
        tuple. list of strings and dictionary. HGVS notations (chrom:g.startRef>Alt) and mapping of HGVS to input.
    '''
    variants['Submission'] = (variants['Chromosome'].astype(str) + '\t' + variants['Start'].astype(str) + '\t.\t' + variants['Ref'] + '\t' + variants['Alt'])

    try: submission_map = dict(zip(variants['Submission'], variants['Input']))

    except: submission_map = None

    return variants, submission_map

