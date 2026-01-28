import json
import obonet
import pandas as pd


def parse_doid_obo(doid_path):
    '''
    Parses the DOID OBO file to create a mapping of disease synonyms to their primary names.
    Creates a JSON file 'doid_map.json' in the 'data/DOID/' directory.
    '''
    graph = obonet.read_obo(doid_path)

    synonym_to_name = {}
    for doid, data in graph.nodes(data=True):
        name = data.get('name').lower()

        for syn in data.get('synonym', []):
            # Synonym strings look like: '"--------" EXACT []'
            # Extract the quoted text:
            parts = syn.split('"')

            if len(parts) >= 3:
                syn_text = parts[1].strip().lower()
                synonym_to_name[syn_text] = name

    json.dump(synonym_to_name, open('data/DOID/doid_map.json', 'w'), indent=4)


def parse_mondo_obo(mondo_path):
    graph = obonet.read_obo(mondo_path)

    mondo_term_map = {}
    mondo_xref_map = {}
    for term, data in graph.nodes(data=True):
        name = data.get("name")

        id = term

        mondo_term_map[id] = name

        for x_id in data.get('xref', []):
            x_ref = x_id.split(' ')[0]
            mondo_xref_map[x_ref] = id

    json.dump(mondo_term_map, open('data/MONDO/mondo_term_map.json', 'w'), indent=4)
    json.dump(mondo_xref_map, open('data/MONDO/mondo_xref_map.json', 'w'), indent=4)


def parse_ensembl_gff3(gff3_path, species):
    '''
    Parses Ensembl GFF3 file to extract gene information and saves it as a parquet file.
    '''
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    gene_df = pd.read_csv(gff3_path, sep='\t', comment='#', compression='gzip', header=None, names=cols, low_memory=False)
    gene_df = gene_df[(gene_df['type'] == 'gene')]

    # Extract fields using regex
    gene_df['Name'] = gene_df['attributes'].str.extract(r'Name=([^;]+)')
    gene_df['description'] = gene_df['attributes'].str.extract(r'description=([^;]+)')
    gene_df['biotype'] = gene_df['attributes'].str.extract(r'biotype=([^;]+)')
    gene_df['gene_id'] = gene_df['attributes'].str.extract(r'ID=gene:([^;]+)')

    # Keep only needed columns
    gene_df = gene_df[['Name', 'description', 'biotype', 'seqid', 'start', 'end', 'strand', 'gene_id']]

    gene_df.to_parquet(f'data/Ensembl/{species}_gene_data.parquet')

