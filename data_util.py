import obonet, json

def parse_doid_obo():
    '''
    Parses the DOID OBO file to create a mapping of disease synonyms to their primary names.
    Creates a JSON file 'doid_map.json' in the 'data/DOID/' directory.
    '''
    graph = obonet.read_obo('data/doid/doid.obo')

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


if __name__ == '__main__':
    parse_doid_obo()