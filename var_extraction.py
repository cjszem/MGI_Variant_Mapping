from util import hvar_to_output, mvar_to_output

gene = 'ACVR1'
start = 157774114
end = 157774114
chrom = 2
ref = 'C'
alt = 'T'

# gene = 'MACF1'
# chrom = 1
# start = 39468726
# end = 39468726
# ref = 'T'
# alt = 'G'
# assembly = 'GRCh38'

hvar_to_output(gene, chrom, start, end, ref, alt)
mvar_to_output(gene)