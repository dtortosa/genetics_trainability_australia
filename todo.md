1. urgent stuff

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
        - we have the three problematic cases with sex unknown and two more with the opposite sex
    - ASK DAVID about what sex to use. I guess we should stick to the reported in the phenotype data
    - there also 2 males and two snps in X that have two genotypes!
        -       FID        ID   snp_name
        - 0  combat  7691CPSO  rs5917469
        - 1  combat  7684BSAO  rs5955017
    - check genotype calls that are I or D
        - how plink deals with this?
    - ask David about the samples without phenotype in the excel file?
        - the last 42
    - Ask david about the sample included in first bath but without phenotype data


2. Tutorials

- general tutorial of QC with plink
    - https://genomicsbootcamp.github.io/book/
- Quality Control Procedures for Genome-Wide Association Studies
    - there is a version of 2011 and new one in 2022 with github repo
    - https://github.com/RitchieLab/GWAS-QC#basic-overview
- A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis
- Omics Data Preprocessing for Machine Learning: A Case Study in Childhood Obesity
- Genetic prediction of complex traits with polygenic scores: a statistical review
- Addressing the challenges of polygenic scores in human genetic research