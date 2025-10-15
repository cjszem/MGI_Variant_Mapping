from human_util import hvar_to_output

gene = 'ACVR1'
loc = 157774114
chrom = 2
ref = 'C'
alt = 'T'

# gene = 'MACF1'
# chrom = 1
# loc = 39468726
# ref = 'T'
# alt = 'G'
# assembly = 'GRCh38'

hvar_to_output(gene, chrom, loc, ref, alt)

