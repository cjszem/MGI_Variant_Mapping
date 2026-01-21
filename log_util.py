import yaml
import logging


# Load Config
with open('config.yaml') as f:
    config = yaml.safe_load(f)


# Logging setup
logging.basicConfig(filename=config['logging']['file'],
                    level=config['logging']['level'],
                    format=config['logging']['format'])


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

