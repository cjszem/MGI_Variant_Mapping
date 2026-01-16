# Data Sources for Variant-to-Allele Mapping
## MGI Data

1. data/MGI/VARIANT-ALLELE-MGI.tsv

    Description: MGI transcribed allele entries with genomic positions.

    Usage: Provides allele-level genomic coordinates for mapping variants to mouse alleles.

2. data/MGI/HOM_ProteinCoding.rpt

    Description: MGI gene homology report.

    Format:
    ```
    MGI_ID    MusGeneSymbol    MusEntrezGeneID    MusHGNC_ID    HumGeneSymbol    HumEntrezGeneID
    ```
    

    Source: https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt

    Usage: Links mouse and human genes for cross-species mapping.

3. data/MGI/MGI_DiseaseMouseModel.rpt

    Description: MGI disease models.

    Usage: Associates mouse models with disease phenotypes.

4. data/MGI/MGI_PhenotypicAllele.rpt

    Description: MGI alleles and phenotype information.
    
    Usage: Provides allele type and phenotype annotations.


## Disease Ontology (DOID)

1. data/DOID/doid.obo

    Description: DOID ontology file containing disease terms and hierarchy.

    Source: https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo

2. data/DOID/doid_map.json

    Description: DOID synonyms mapped to disease names.

    Created by: ```data_util.parse_doid_obo(doid_path)```

    Usage: Enables synonym resolution for disease terms.


## Ensembl Gene Data

1. data/Ensembl/human_gene_data.parquet

    Description: Human gene data extracted from Ensembl GFF3.

    Created by: ```parse_ensembl_gff3('data/Ensembl/Homo_sapiens.GRCh38.115.gff3.gz', species)```


2. data/Ensembl/mouse_gene_data.parquet

    Description: Mouse gene data extracted from Ensembl GFF3.

    Created by: ```parse_ensembl_gff3('data/Ensembl/Mus_musculus.GRCm39.115.gff3.gz', species)```



## ClinVar

1. data/NCBI/clinvar.grch38.vcf.gz

    Description: ClinVar entries for GRCh38 variants.

    Usage: Provides clinically relevant human variants for mapping.

2. clinvar.grch38.vcf.gz.tbi

    Description: Tabix file for quick indexing of data/NCBI/clinvar.grch38.vcf.gz