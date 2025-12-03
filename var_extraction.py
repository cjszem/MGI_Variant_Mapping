from util import hvar_to_output, mvar_to_output, score_ortho_vars
import pandas as pd

gene = 'ACVR1'
start = 157774114
end = 157774114
chrom = '2'
ref = 'C'
alt = 'T'

# gene = 'MACF1'
# chrom = 1
# start = 39468726
# end = 39468726
# ref = 'T'
# alt = 'G'
# assembly = 'GRCh38'

human_gene_df, human_prt_df = hvar_to_output(gene, chrom, start, end, ref, alt)
mouse_gene_df, mouse_prt_df = mvar_to_output(gene)

# human_prt_df = pd.read_csv('hum_protein_df.csv')
# mouse_prt_df = pd.read_csv('mouse_protein_df.csv')

var_scores_df = score_ortho_vars(human_prt_df, mouse_prt_df)