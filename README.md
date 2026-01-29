# MGI Variant Mapping
The Mouse-Human Variant Mapping project is to develop a software pipeline to identify equivalent genome sequence variants in protein-coding genes in mouse and human genomes.

## Abstract
This thesis will be completed in partnership with Dr. Carol Bult at JAX Lab, aiming to create a software pipeline to identify equivalent genome sequence variants in protein-coding genes in mouse and human genomes. This software will report mouse variants equivalent to the inputted human variant, ranking them by whether the variant occurs in an orthologous gene, the same gene exon, the same protein domain, results in the same molecular consequence and causes the same amino acid change. Each discovered mouse variant will receive a similarity score to the inputted human variant. This software will then be evaluated by its performance in a clinical setting utilizing real genomic variants from clinically performed genome wide association studies and its effectiveness in a clinical practice will be judged.

## Background
Genome sequencing has become increasingly prevalent in clinical practice; however, efficiently connecting sequence variants identified in patients with clinical significance remains a significant challenge. The Jackson Lab Informatics team has developed the Mouse Genome Database (MGD, http://www.informatics.jax.org), which contains information acquired from the scientific literature about genotype-to-phenotype associations in the laboratory mouse (1). These data can be used to assert clinical relevance of human variants; however, it is currently very labor intensive to perform this on a large collection of genomic variants. The software developed for this project will bridge the gap between human sequence variants and mouse models with known phenotypes or human disease relevance. This software will allow for the immediate discovery of clinically relevant mouse models, which can be used to determine the best course of action and assess possible health consequences of variants discovered in patients, with a final goal of providing more successful health outcomes and personalized treatment for patients.

## Running
To run the server locally follow the steps below:
1. Download the repository and ceate a working directory.

2. Download a dockerized version of ensembl-vep and homo_sapiens and mus_musculous caches.

2. Create the conda environment using the provided yml file. All remaining steps occur within this environment.
```
conda env create -f environment.yml
```

3. Run the server using uvicorn.
```
uvicorn server:app --reload
```

4. Insert a list of variants into the search field and submit the 

5. Examine output.

## Output

1. Human gene table: Provides information on the given variant on a gene level.
```
Gene Symbol Description Biotype Chromosome  Start	End	Strand	Ensembl_ID	Accession
ACVR1 activin A receptor type 1	protein_coding	2	157736249	157876347	-	ENSG00000115170	HGNC:171
```

2. Human variant table: Provides information on the given variant on a transcript level.
```
Gene Symbol	Transcript ID	Biotype	Exon Rank	Pfam Domain ID	Pfam Domain Name	Polyphen Prediction	Polyphen Score	Molecular Consequence	Codon Switch	Amino Acids	Associated Diseases	MONDO
ACVR1	ENST00000434821	protein_coding	6/11	PF08515	Transforming growth factor beta type I GS-motif	probably_damaging	0.999	missense_variant	cGc/cAc	R206H	fibrodysplasia ossificans progressiva MONDO:0007606
```
3. Mouse gene table: Provides information on the given orthologous variant on a gene level.
```
Gene Symbol	Description	Biotype	Chromosome	Start	End	Strand	Ensembl_ID	Accession
Acvr1	activin A receptor%2C type 1	protein_coding	2	58278656	58457169	-	ENSMUSG00000026836	MGI:87911
```
4. Mouse model table: Provides information on the given orthologous variant on a transcript level.
```
Gene Symbol	AlleleID	AlleleSymbol	Transcript ID	Biotype	Exon Rank	Pfam Domain ID	Pfam Domain Name	Molecular Consequence	Codon Switch	Amino Acids	Associated Diseases	MONDO
Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000056376	protein_coding	6/11	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000090935	protein_coding	8/13	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000112599	protein_coding	5/10	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
Acvr1	MGI:6140231	Acvr1tm1Glh	ENSMUST00000112601	protein_coding	7/12	PF08515	Transforming growth factor beta type I GS-motif	missense_variant	cGC/cAT	R206H	fibrodysplasia ossificans progressiva	MONDO:0007606
...
```
5. Model equivalence table: Boolean match values for each mouse transcript
```
Allele ID	Allele Symbol	Transcript ID	Biotype Match	Consequence Match	AA Match	AA Position Match	Exon Match	Domain Match	Disease Match	total_score
MGI:6140231	Acvr1tm1Glh	ENSMUST00000056376	true	true	true	true	true	true	true	100
MGI:6140231	Acvr1tm1Glh	ENSMUST00000090935	true	true	true	true	false	true	true	85
MGI:6140231	Acvr1tm1Glh	ENSMUST00000112599	true	true	true	true	false	true	true	85
MGI:6140231	Acvr1tm1Glh	ENSMUST00000112601	true	true	true	true	false	true	true	85
...
```
6. Mouse phenotype table: Provides detailed phenotype information on all alleles in gene.
```
Gene Symbol	AlleleID	AlleleSymbol	Allele Symbol	Phenotypes
Acvr1	MGI:6140231	Acvr1tm1Glh	Acvr1	MP:0005390,skeleton phenotype
Acvr1	MGI:1857711	Acvr1tm1Enl	Acvr1	MP:0001672,abnormal embryo development,MP:0001675,abnormal ectoderm development,MP:0001683,absent mesoderm,MP:0001695,abnormal gastrulation,MP:0001698,decreased embryo size,MP:0001710,absent amniotic folds,MP:0002230,abnormal primitive streak formation,MP:0003087,absent allantois,MP:0005030,absent amnion,MP:0009593,absent chorion,MP:0011098,embryonic lethality during organogenesis, complete penetrance,MP:0011186,abnormal visceral endoderm morphology,MP:0011190,thick embryonic epiblast
...
```
All outputs are also stored in the corresponding CSV files.
