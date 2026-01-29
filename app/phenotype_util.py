import yaml
import pandas as pd


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


# MGI_phenotype_df = pd.read_csv(config['paths']['mgi_phenotype'], sep='\t')
# MGI_phenotype_df = MGI_phenotype_df.set_index('Allele ID')

def allele_phenotype_match(allele_df):
    MGI_phenotype_df = pd.read_csv(config['paths']['mgi_phenotype'], sep='\t', usecols=['Allele ID', 'Phenotypes']) 
    MGI_phenotype_df = MGI_phenotype_df.set_index('Allele ID')
    
    phenotype_df = allele_df.join(MGI_phenotype_df, on='AlleleID', how='inner')

    phenotype_df.drop(['Chromosome', 'Start', 'End', 'Ref', 'Alt', 'Molecular Consequence'], inplace=True, axis=1)

    phenotype_df['Phenotypes'] = (phenotype_df['Phenotypes'].replace(',', '-').fillna('').str.split('|'))

    return phenotype_df