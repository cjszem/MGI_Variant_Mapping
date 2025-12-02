# MGI Variant Mapping
The Mouse-Human Variant Mapping project is to develop a software pipeline to identify equivalent genome sequence variants in protein-coding genes in mouse and human genomes.

## Abstract
This thesis will be completed in partnership with Dr. Carol Bult at JAX Lab, aiming to create a software pipeline to identify equivalent genome sequence variants in protein-coding genes in mouse and human genomes. This software will report mouse variants equivalent to the inputted human variant, ranking them by whether the variant occurs in an orthologous gene, the same gene exon, the same protein domain, results in the same molecular consequence and causes the same amino acid change. Each discovered mouse variant will receive a similarity score to the inputted human variant. This software will then be evaluated by its performance in a clinical setting utilizing real genomic variants from clinically performed genome wide association studies and its effectiveness in a clinical practice will be judged.

## Background
Genome sequencing has become increasingly prevalent in clinical practice; however, efficiently connecting sequence variants identified in patients with clinical significance remains a significant challenge. The Jackson Lab Informatics team has developed the Mouse Genome Database (MGD, http://www.informatics.jax.org), which contains information acquired from the scientific literature about genotype-to-phenotype associations in the laboratory mouse (1). These data can be used to assert clinical relevance of human variants; however, it is currently very labor intensive to perform this on a large collection of genomic variants. The software developed for this project will bridge the gap between human sequence variants and mouse models with known phenotypes or human disease relevance. This software will allow for the immediate discovery of clinically relevant mouse models, which can be used to determine the best course of action and assess possible health consequences of variants discovered in patients, with a final goal of providing more successful health outcomes and personalized treatment for patients.

## Testing
The steps for testing the proejct in is current state are as follows:

1. Download the repository and ceate a working directory.

2. Create the conda environment using the provided yml file. All remaining steps occur within this environment.
```
conda env create -f environment.yml
```

3. Insert test variables into var_extraction.py.
```
gene = 'ACVR1'
start = 157774114
end = 157774114
chrom = 2
ref = 'C'
alt = 'T'
```

4. Run var_extraction, either in terminal or using IDE 'Run' button.
```
python var_extraction.py
```

5. Examine output (both printed and generated CSV files available).

## Output

There are currently 7 printed outputs for examination:
1. Human gene table: Provides information on the given variant on a gene level.
```
  Organism   Gene                Description      HGNC       Ensembl ID         Biotype  Strand  Chr      Start        End Ref Alt
0    Human  ACVR1  activin A receptor type 1  HGNC:171  ENSG00000115170  protein_coding      -1    2  157774114  157774114   C   T
```

2. Human protein table: Provides information on the given variant on a transcript level.
```
     Transcript ID         Biotype Exon Rank Pfam Domain ID  ... Amino Acids refAA  varAA                                Associated Diseases
0  ENST00000434821  protein_coding      6/11        PF08515  ...       R206H     R      H  Progressive myositis ossificans; not provided;...
```
3. Mouse gene table: Provides information on the given orthologous variant on a gene level.
```
  Organism   Gene                 Description       HGNC          Ensembl ID         Biotype  Strand
0    Mouse  Acvr1  activin A receptor, type 1  MGI:87911  ENSMUSG00000026836  protein_coding      -1
```
4. Mouse allele table: Provides all extracted mouse alleles with genomic data.
```
         Gene     AlleleId Chromosome     Start       End Ref Alt Molecular Consequence
55198   Acvr1  MGI:6140231          2  58364210  58364211  GC  AT      missense_variant
64079   Acvr1  MGI:5763014          2  58364211  58364211   C   T      missense_variant
67727   Acvr1  MGI:5471642          2  58364211  58364211   C   T      missense_variant
101505  Acvr1  MGI:6414911          2  58352976  58352976   C   A      missense_variant
112037  Acvr1  MGI:6414907          2  58352976  58352976   C   A      missense_variant
```
5. Mouse protein table: Provides information on the given orthologous variant on a transcript level.
```
       AlleleId       transcript_id         biotype   exon  domains  ... refAA varAA primaryIdentifier                           ontologyName  ontologyID
0   MGI:6140231  ENSMUST00000056376  protein_coding   6/11  PF08515  ...     R     H       MGI:6140231  fibrodysplasia ossificans progressiva  DOID:13374
1   MGI:6140231  ENSMUST00000090935  protein_coding   8/13  PF08515  ...     R     H       MGI:6140231  fibrodysplasia ossificans progressiva  DOID:13374
2   MGI:6140231  ENSMUST00000112599  protein_coding   5/10  PF08515  ...     R     H       MGI:6140231  fibrodysplasia ossificans progressiva  DOID:13374
3   MGI:6140231  ENSMUST00000112601  protein_coding   7/12  PF08515  ...     R     H       MGI:6140231  fibrodysplasia ossificans progressiva  DOID:13374
4   MGI:6140231  ENSMUST00000126407  protein_coding    NaN  PF08515  ...   nan  None       MGI:6140231  fibrodysplasia ossificans progressiva  DOID:13374
5   MGI:5763014  ENSMUST00000056376  protein_coding   6/11  PF08515  ...     R     H       MGI:5763014  fibrodysplasia ossificans progressiva  DOID:13374
...
```
6. Unused mouse alleles: List of all possible orthologous mouse alleles without documented genomic information.
```
['MGI:7486439', 'MGI:7309408', 'MGI:2651395', 'MGI:1857711', 'MGI:1931734', 'MGI:3696672', 'MGI:7309409', 'MGI:7376166', 'MGI:4399476', 'MGI:5548914']
```
7. Variant scoring table: Boolean match values for each mouse transcript
```
    biotype_match  exon_match  domain_match  consequence_match  AA_match  disease_match  total_score
0            True        True          True               True      True           True   100.000000
1            True       False          True               True      True           True    83.333333
2            True       False          True               True      True           True    83.333333
3            True       False          True               True      True           True    83.333333
4            True       False          True              False     False           True    50.000000
5            True        True          True               True      True           True   100.000000
...
```

All outputs are also stored in the corresponding CSV files.
