1. urgent stuff


PREGUNTA ESTO A DAVID NEXT TIME

- pheno_data
    - samples without pheno data
        - also about the 41 samples the bottom with no phenotyipic daata
        - I guess these samples should be removed from all analyses?
            - remove after quality control? maybe it can be useul for PCAs and pop structure, see tutorials
    - ask david that from sample 1161 to 1376, age is integer, not float, in contrast with almost all the rest samples. This is ok?
    - in some phenotypes, some some samples have value of 0 and others have no value. I guess zero should be NA, right?
        - body mass week 1 and 8
        - VO2 max week 1
    - sample 1194, the value for week 8 beep includes a letter: 11.1O. I guess I can safely change that "o" letter by zero.
    - I guess that the sheet "DNA with only wk1" includes genotyped samples with only data for the first week, not week 8. So I should only use the sheet "All DNA samples" and discard the 42 samples at the bottomn with NA for all columns except the AGRF code.
    - some rows are coloured, there is something special about these samples it could be relevant for the analysis?
    - there is data about ancestry? I am, african, european... I will do analyses to detect genetic outliers but it would be good if we have this data.


- Things to check when QC:
    - batch data
        - There are two batches (ILGSA24-17303 and ILGSA24-17873), being the data separated for these. 
         
        - In ILGSA24-17303.zip, we have the final reports for each 216 samples, along with the sample and snp maps.
            - In the initial_stuff folder there is a zip called "ILGSA24-17303.zip" that I may downloaded from the initial location where this data was stored in summer 2022. There are Plink files, but I am not sure this is the correct data and I cannot find the final_report files.
            - In 17873, we have the IDAT files with probs intensity from the microarrays used to genotype (first zips), the final reports (CAGRF20093767.zip) and a inputs for plink. But all of this only for 1248 individuals, not the whole cohort.
            - CAGRF20093767.zip includes the final reports of 1248 individuals, along with the sample and snp maps.
                warning [CAGRF20093767.zip]:  32332459452 extra bytes at beginning or within zipfile (attempting to process anyway)
                error [CAGRF20093767.zip]:  start of central directory not found; zipfile corrupt. (please check that you have transferred or created the zipfile in the appropriate BINARY mode and that you have compiled UnZip properly)
                - CHECK WARNING
                    - solution for the error: https://askubuntu.com/questions/54904/unzip-error-end-of-central-directory-signature-not-found
    - Check the genome build. I have manually checked some SNPs and they have the position of hg38.
    - YOU ARE NOT ADDED CENTIMORGANS, SO BE CAREFUL WITH WHAT PLINK DO ABOUT LD.
    - Remove zero chromosomes
        - also MT, X, Y and XY?
    - think how to combine the batches
        - https://github.com/satijalab/seurat/issues/3631
        - https://groups.google.com/g/plink2-users/c/SPRgIOBENgI/m/T-4S3TNcAAAJ
    - Check strand, because the foward strand of illumina not always matches that of dbSNP
        - StrandScript is a software to check and solve that
    - Check the illumina pdf report for the first batch (zip file in initial stuff), because they say there are three problematic samples.
            - these three sample also have unknown sex according to illumina data!!
        - check sex between illumina and real sex data (pheno_data)
            - remove uknown sex individuals?
        - we have the three problematic cases with sex unknown and two more with the opposite sex between illumina and phenotype data
    - ASK DAVID about what sex to use. I guess we should stick to the reported in the phenotype data
    - there also 2 males and two snps in X that have two genotypes!
        -       FID        ID   snp_name
        - 0  combat  7691CPSO  rs5917469
        - 1  combat  7684BSAO  rs5955017
    - check genotype calls that are I or D
        - these are indels?
            - 1:207754848-GATAA-G
                - rs5780395 is its name and it is an indel
                - https://www.ncbi.nlm.nih.gov/snp/?term=rs5780395
                - https://useast.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=1:207581500-207581510;source=dbSNP;v=rs5780395;vdb=variation;vf=3120483
        - how plink deals with this?
            - plink has flag  --snps-only to keep snps
                - https://www.biostars.org/p/378475/
    - quitar mitocontrodial and other strange chromosomes?


    - IMPUTATION?
        - using high cove 1000 KGP1?
            - One of the major applications of the phase 3 1kGP call set has been its widespread use as a reference panel for variant imputation in sparse, array-based genotyping data with a goal of improving the statistical power of downstream genome-wide association studies (GWAS) and facilitating fine-mapping of causal variants. As part of this publication, we release an improved reference imputation panel based on the high-coverage WGS consisting of SNV, INDEL, and SV calls across the 3,202 1kGP samples, including full trios.
            - https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue

- PRS
    - We need a base cohort and a target cohort
    - See explanations in "03a_association_analyses.py"

- Machine learning bishop, tneemls miestra sufociente para training eval? Si no cogemos snps previamentr asociados y con eso hacemos el polygenic score en caso de que no tengamos muestra para hacer nuestra propia seleccion y luego testarla.

- repasar en general los pasos? como la perdida de 300K sNPs due to maf filtering?

2. Tutorials
- Deep learning-based polygenic risk analysis for Alzheimer’s disease prediction
    - https://www.nature.com/articles/s43856-023-00269-x
- cool tutorial merging data with 1KGP using plink!!
    - https://www.biostars.org/p/335605/
- tutorial christopher chang
    - https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3
- general tutorial of QC with plink
    - https://genomicsbootcamp.github.io/book/
- Omics Data Preprocessing for Machine Learning: A Case Study in Childhood Obesity
    - paper augusto, shap values.
- Identifying and mitigating batch effects in whole genome sequencing data
- Quality Control Procedures for Genome-Wide Association Studies
    - there is a version of 2011 and new one in 2022 with github repo
    - https://github.com/RitchieLab/GWAS-QC#basic-overview
- Integrative polygenic risk score improves the prediction accuracy of complex traits and diseases 
- PGSbuilder: An end-to-end platform for human genome association analysis and polygenic risk score predictions
    - use the paper as s guide of good practice PRS?
- Interpreting population and family-based genome-wide association studies in the presence of confounding
- A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis
- Genetic prediction of complex traits with polygenic scores: a statistical review
- Addressing the challenges of polygenic scores in human genetic research
- https://www.varianteffect.org/
- functional tools
    - Regulatory dissection of the severe COVID-19 risk locus introgressed by Neanderthals
- PRSet: Pathway-based polygenic risk score analyses and software
- Testing for differences in polygenic scores in the presence of confounding
- Reply to: Genome-wide association studies of polygenic risk score-derived phenotypes may lead to inflated false positive rates
- Integrative polygenic risk score improves the prediction accuracy of complex traits and diseases
- The PPARGC1A Gly482Ser polymorphism is associated with elite long-distance running performance
    - 2023, Robert M. Erskine
- to annotate genes
    - Gene-metabolite annotation with shortest reactional distance enhances metabolite genome-wide association studies result