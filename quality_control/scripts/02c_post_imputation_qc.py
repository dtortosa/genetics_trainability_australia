#!/usr/bin/env python3.9
# coding: utf-8
    #to run this script: chmod +x script.py; ./script.py
    #!/bin/sh does not work with my terminal en msi of David.
    #if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
    #you can save the output and the errors
        #./script.py > script.out #only output
        #./script.py 2> error.out #only error
        #./script.py > script.out 2> error.out #both in different files
        #./script.py > script.out 2>&1 #both in the same file
        #https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/



####################################
######## POST IMPUTATION QC ########
####################################

#This script is based on the following tutorials:
    #1. Quality Control Procedures for Genome-Wide Association Studies
        #https://github.com/RitchieLab/GWAS-QC
        #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
    #2. Tutorial: a guide to performing polygenic risk score analyses
        #https://www.nature.com/articles/s41596-020-0353-1
        #https://choishingwan.github.io/PRS-Tutorial/
    #3. A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
    #4. Genome-wide association studies
        #https://www.nature.com/articles/s43586-021-00056-9
    #5. A guide to genome-wide association analysis and post-analytic interrogation
        #https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
    #6. Data Management and Summary Statistics with PLINK
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
    #7. Genomics Boot Camp
        #https://genomicsbootcamp.github.io/book/
    #8. Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
    #9. Omics Data Preprocessing for Machine Learning: A Case Study in Childhood Obesity
        #https://www.mdpi.com/2073-4425/14/2/248






#######################################
# region INITIAL ANOTATIONS AND STEPS #
#######################################


########################################
# define function to print text nicely #
########################################

#text="checking function to print nicely"
#header=1
def print_text(text, header=2):
    if header==1:
        print("\n#######################################\n#######################################")
        print(text)
        print("#######################################\n#######################################")
    elif header==2:
        print("\n###### " + text + " ######")
    elif header==3:
        print("\n## " + text + " ##")
print_text("checking function to print nicely", header=1)
print_text("checking function to print nicely", header=2)



########################################
# define function to run bash commands #
########################################

#create a wrapper for subprocess.run in order to define a set of arguments and avoid typing them each time. We will ensure that we are using bash always and not sh.
from subprocess import run, PIPE
#command="ls"
def run_bash(command, return_value=False):

    #run the command
    complete_process = run(
        command, 
        shell=True,
        executable="/bin/bash", 
        stdout=PIPE,
        stderr=PIPE, 
        text=True)
    #we have to use popen in order to ensure we use bash, os.system does not allow that
        #shell=True to execute the command through the shell. This means that the command is passed as a string, and shell-specific features, such as wildcard expansion and variable substitution, can be used.
            #THIS IS DANGEROUS IF UNTRUSTED DATA
        #executable="/bin/bash" to ensure the use of bash instead of sh
        #stdout=PIPE to capture the output into an python object. You can also capture the error doing stderr=PIPE. stdout and stderr are the standard output and error
            #you could also use capture_output=True to capture both stdout and stderr
        #text=True will return the stdout and stderr as string, otherwise as bytes
            #https://www.datacamp.com/tutorial/python-subprocess
            #https://docs.python.org/3/library/subprocess.html#subprocess.run

    #this generated a CompletedProcess instance where you can get
        #args: The command and arguments that were run.
        #returncode: The return code of the subprocess.
        #stdout: The standard output of the subprocess, as a bytes object.
        #stderr: The standard error of the subprocess, as a bytes object.

    #if the return code is 0, i.e., success 
        #https://askubuntu.com/a/892605
    if complete_process.returncode==0:

        #if stderr is empty
        if complete_process.stderr=="":

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
        else:

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #print the standard error without stopping
            print(complete_process.stderr)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
    else:
        #print the standard error and stop
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

#test it
print_text("check behaviour run_bash", header=1)
print_text("see working directory", header=2)
run_bash("pwd")
print_text("list files/folders there", header=2)
run_bash("ls")

# endregion






##########################################
# region STEPS FOLLOWED IN TOPMED SERVER #
##########################################

#Genotype imputation is the statistical infer- ence of unobserved genotype, which enables scientists to reconstruct the missing data in each genome and accurately evaluate the ev- idence for association at genetic markers that were not genotyped. Imputation has become an essential component of GWAS because it increases power, facilitates meta-analysis, and aids in the interpretation of signals which can directly influence the results of an analysis (Das et al., 2016; Li, Willer, Sanna, & Abecasis, 2009; Marchini & Howie, 2010; Verma et al., 2014). Genotype imputation is achieved by comparing short stretches of an individual genome against stretches of previously characterized reference genomes. It is usually performed on SNPs, which are the most com mon type of genetic variation (Verma et al., 2014).
    #From Richies tutorial

#If you check the slides from the TOPMed server, you can see how they are basically matching haplotypes of our data with haplotypes in the reference panel to fill the gaps. In other words, if you have A....A....G in your sample, and in the reference panel these A - A - G are included in the haplotype ATATG, you can fill the gaps in your data saying that you ATATG. This is of course base in the linkage disqueilibrium between SNPs found in a given ancestry. Because of that it is important to use a reference panel that correctlty matches our data.
    #https://raw.githubusercontent.com/genepi/imputationserver-ashg/main/slides/MIS_Workshop_2023.pdf

#############################################
## options selected and input requeriments ##
#############################################

#Run/Genotype Imputation (Minimac4) 2.0.0-beta3
    #Minimac4 is a lower memory and more computationally efficient implementation of the genotype imputation algorithms in minimac/mininac2/minimac3.
    #https://github.com/statgen/Minimac4

#TOPMed Imputation Server provides a free genotype imputation service using Minimac4. You can upload phased or unphased GWAS genotypes and receive PHASED and IMPUTED genomes in return. This server offers imputation from the TOPMed reference panel. For all uploaded datasets an extensive QC is performed.

#Citation:
    #Das S, Forer L, Schönherr S, Sidore C, Locke AE, Kwong A, Vrieze S, Chew EY, Levy S, McGue M, Schlessinger D, Stambolian D, Loh PR, Iacono WG, Swaroop A, Scott LJ, Cucca F, Kronenberg F, Boehnke M, Abecasis GR, Fuchsberger C. Next-generation genotype imputation service and methods. Nature Genetics 48, 1284–1287 (2016).

#Reference panel: TOPMed (Version R3 on GRC38). 
    #The TOPMed Imputation Server offers genotype imputation for the TOPMed reference panel, which is the largest and most accurate panel available amongst the two imputation servers.
    #Number of Samples: 133,597
    #Sites (chr1-22): 445,600,184
    #Chromosomes: 1-22, X
    #Imputation Server: https://imputation.biodatacatalyst.nhlbi.nih.gov
    #Website: https://www.nhlbiwgs.org/

#Input files: TOPMed Imputation Server accepts VCF files compressed with bgzip. Please make sure the following requirements are met:
    #Create a separate vcf.gz file for each chromosome.
        #YES
    #Variants must be sorted by genomic position.
        #YES
    #GRCh37 or GRCh38 coordinates are required.
        #HG38
    #If your input data is GRCh37/hg19, please ensure chromosomes are encoded without prefix (e.g. 20). If your input data is GRCh38/hg38, please ensure chromosomes are encoded with prefix 'chr' (e.g. chr20).
        #chr20 as we are using hg38
    #VCF files need to be version 4.2 (or lower). This is specified in the VCF file header section.
        #VCF4.2, plink´s default
    #Must contain GT field in the FORMAT column. All other FORMAT fields will be ignored. (if you are seeing problems with very large uploads, it may help to remove other FORMAT fields)
        #YES
    #Due to server resource requirements, there is a maximum of 25k samples per chromosome per job (and a minimum of 20 samples). Please see the FAQ for details.
        #1203 samples.

#Build
    #Please select the build of your data. Currently, the options hg19 and hg38 are supported. The TOPMed Imputation Server automatically updates the genome positions of your data (liftOver). The TOPMed reference panel is based on hg38 coordinates.
    #HG38 in our case.

#rsq filter
    #For each variant, how confident can we be that the imputation dosages are sufficiently “accurate” for association analyses? Measure of confidence in imputed dosages: “Rsq” column [range 0-1]
    #To minimize the file size, the Imputation Server includes a r2 filter option, excluding all imputed SNPs with a r2-value (= imputation quality) smaller than the specified value.
    #It seems this is applied after imputation. So I guess this calculate the correlation between the genotypes of the imputed SNP and the SNP in the reference, which would be the true genotype. According to Ritiche "the Rsq filter must be set, which computes the estimated value of squared correlations between imputed genotypes and true genotypes"
        #https://genome.sph.umich.edu/wiki/MaCH_FAQ
    #Detailed definition:
        #This is the estimated value of the squared correlation between imputed genotypes and true, unobserved genotypes. Since true genotypes are not available, this calculation is based on the idea that poorly imputed genotype counts will shrink towards their expectations based on population allele frequencies alone; specifically 2p where p is the frequency of the allele being imputed.
        #https://genome.sph.umich.edu/wiki/Minimac3_Info_File
    #We are using 0.3: This removes >70% of poorly-imputed SNPs at the cost of <0.5% well-imputed SNPs. This is the recommendation of TOPMed and of Ritchies tutorial.
        #https://genome.sph.umich.edu/wiki/MaCH_FAQ
    #Also the TOPMed slides explain the use of Rsq as a filter after imputation. Removing cases under 0.3 makes it enough for common variants without losing a lot of variants.
        #Minimal Rsq value for common variants
            #≥ 0.30:
        #Minimal Rsq value for low frequency/rare variants
            #≥0.50: Before performing GWAS, remove variants that do not meet these thresholds
        #https://raw.githubusercontent.com/genepi/imputationserver-ashg/main/slides/MIS_Workshop_2023.pdf

#checking the strand:
    #First and foremost, high quality imputation requires that allele calls be on the same physical strand of DNA for both study and reference human genome data. The variability between study sites can yield differences in genotyping platform and calling algorithm. Thus, several algorithms and platforms exist to check strand and perform strand flip (e.g., BEAGLE, SHAPEIT2). Strand check is performed based on three criteria: a) observed alleles, b) minor allele frequencies (MAF), and c) LD pattern across 100-SNP lengths of the genome (Verma et al., 2014). SNPs are discarded from the data accordingly if inconsistent MAF and LD patterns cannot be resolved by flipping the strand. In preparation for phasing, data are subsetted by chromosome, and strand flip is executed to ensure that the SNPs correctly align with the reference panel “+” strand.
        #We have used the Ritchie approach of their github to solve this problem.

#Phasing
    #If your uploaded data is unphased, Eagle v2.4 will be used for phasing. In case your uploaded VCF file already contains phased genotypes, please select the "No phasing" option.
    #The Eagle algorithm estimates haplotype phase using the HRC reference panel. This method is also suitable for single sample imputation. After phasing or imputation you will receive phased genotypes in your VCF files.
        #For haplotype phasing Eagle2 is used. Eagle2 attains high accuracy across a broad range of cohort sizes by efficiently leveraging information from large external reference panels (such as the Haplotype Reference Consortium; HRC).
    #Historically, whole-genome sequencing generated a single consensus sequence without distinguishing between variants on homologous chromosomes. Phased sequencing, or genome phasing, addresses this limitation by identifying alleles on maternal and paternal chromosomes.
    #info from Ritchie
        #After checking for strand, co-localized alleles on the same copy of a chromosome need to be identified in a process known as haplotype phasing. Prior to imputation, pre-phasing can be performed where the haplotype phase is estimated for all of the alleles. The phased data improves imputation speed and can be used for any future data imputation as reference panel improve over time. Although pre-phasing can speed up the imputation process, this step can also introduce error because of any haplotype uncertainty (Howie, Fuchsberger, Stephens, Marchini, & Abecasis, 2012). On the Michigan and NHLBI Trans-Omics for Precision Medicine (TOPMed) Imputation Servers, the algorithm Eagle v2.4 is used for phasing to estimate haplotype phase using the Haplotype Reference Consortium (HRC) reference panel and returns phased genotypes in VCF format (Das et al., 2016).
    #Therefore, we not have phased data. We selected Eagle v2.4, phased output because our data should not be phased.

#Population
    #Please select whether to compare allele frequencies between your data and the reference panel. In case your samples are mixed from different populations, please select Skip to skip the allele frequency check. For mixed populations, no QC-Report will be created.
    #In our case, we have an homogeneus population after cleaning the dataset, so we did the allele frequency comparison with the TOPMed Panel.

#Mode
    #Please select if you want to run Quality Control & Imputation, Quality Control & Phasing Only or Quality Control Only.
    #We select quality control and imputation, everything (this also includes phasing).

#AES 256 encryption
    #All Imputation Server results are returned as an encrypted .zip file by default. If you select this option, we will use stronger AES 256 encryption instead of the default encryption method. However, note that AES encryption does not work with standard unzip programs. If this option is selected, we recommend using 7-zip to open your results.
    #We used this option to enhance security.



####################
## steps followed ##
####################

#https://topmedimpute.readthedocs.io/en/latest/

#Quality Control steps performed by the Imputation Server
    #Create chunks with a size of 10 Mb
    #For each 10Mb chunk we perform the following checks:
        #On Sample level:
            #For chr1-22, a chunk is excluded if one sample has a call rate < 50 %. Only complete chunks are excluded, not samples (see "On Chunk level" above)
        #On Chunk level:
            #Determine amount of valid variants: A variant is valid if it is included in the reference panel. At least 3 variants must be included.
                #This is the absolute number of variants
            #Determine amount of variants found in the reference panel: At least 50% of the variants must be included in the reference panel.
                #While this is in terms of percentage. You could have 3 variants of the chunk present in the reference panel, but having 4 additional variants in that chunk that are not in the panel, so more than 50% of the variants would not be in the panel, i.e., overlap<50%.
            #Determine sample call rate: At least 50 % of the variants must be called for each sample.
            #Chunk exclusion: if (#variants < 3 || overlap < 50% || sampleCallRate < 50%)
        #On Variant level:
            #Check alleles: Only A,C,G,T are allowed
            #Calculate alternative allele frequency (AF): Mark all with a AF > 0.5.
            #Calculate SNP call rate
            #Calculate chi square for each variant (reference panel vs. study data)
                #Difference of allele frequencies between reference panel and our data
            #Determine allele switches: Compare ref and alt of reference panel with study data (A/T and C/G variants are ignored).
            #Determine strand flips: After eliminating possible allele switches, flip and compare ref/alt from reference panel with study data.
            #Determine allele switches in combination with strand flips: Combine the two rules from above.
            #Variant exclusion: Variants are excluded in case of: [a] invalid alleles occur (!(A,C,G,T)), [b] duplicates (DUP filter or (pos - 1 == pos)), [c] indels, [d] monomorphic sites, [e] allele mismatch between reference panel and study, [f] SNP call rate < 90%.
        #Perform a liftOver step, if build of input data and reference panel does not match (b37 vs b38).
        #Our case
            #1 Chunk(s) excluded: reference overlap < 50.0%
            #4 Chunk(s) excluded: < 20 SNPs
            #Remaining chunk(s): 288

#Phasing
    #Execute for each chunk eagle2.4 is used (we use an overlap of 5 Mb). For example, chr20:1-20000000 and reference population EUR:
        #./eagle
            #--vcfRef HRC.r1-1.GRCh37.chr20.shapeit3.mac5.aa.genotypes.bcf
            #--vcfTarget chunk_20_0000000001_0020000000.vcf.gz  
            #--geneticMapFile genetic_map_chr20_combined_b37.txt
                #this is an old example, but they are now using hg38
            #--outPrefix chunk_20_0000000001_0020000000.phased
            #--bpStart 1 
            #--bpEnd 25000000
            #--allowRefAltSwap
            #--vcfOutFormat z
    #Please note: Target-only sites for unphased data are not included in the final output.

#Imputation
    #For each chunk, run minimac in order to impute the phased data (we use a window of 500 kb):
        #./Minimac4 
            #--refHaps HRC.r1-1.GRCh38.chr1.shapeit3.mac5.aa.genotypes.m3vcf.gz
            #--haps chunk_1_0000000001_0020000000.phased.vcf
            #--start 1
            #--end 20000000
            #--window 500000
            #--prefix chunk_1_0000000001_0020000000
            #--cpus 4
            #--chr 20
            #--noPhoneHome
            #--format GT,DS,GP
            #--allTypedSites
            #--meta
            #--minRatio 0.00001
            #--referenceEstimates
            #--map B38_MAP_FILE.map

#Compression and Encryption
    #Merge all chunks of one chromosome into one single vcf.gz
    #Encrypt data with one-time password
    #Unzip command is not working
        #Please check the following points:
            #(1) When selecting AES256 encryption, please use 7z to unzip your files (Debian: sudo apt-get install p7zip-full). For the default encryption option, all common .zip decompression programs should work. 
            #(2) If your password includes special characters (e.g. \), please put single or double quotes around the password when extracting it from the command line (e.g. 7z x -p"PASSWORD" chr_22.zip).

#Security measures
    #https://topmedimpute.readthedocs.io/en/latest/data-sensitivity/

#FAQ:
    #https://topmedimpute.readthedocs.io/en/latest/faq/

# endregion





###############################################
# region INITIAL RESULTS OF IMPUTATION - TEXT #
###############################################

######################
## input validation ##
######################

#22 valid VCF file(s) found
#Samples: 1203
#Chromosomes: 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9
#SNPs: 246628
#Chunks: 293
#Datatype: unphased
#Build: hg38
#Reference Panel: apps@topmed-r3@1.0.0 (hg38)
#Population: all
#Phasing: eagle
#Mode: imputation
#Rsq filter: 0.3

#Summary: these are the options selected in Ritchies tutorial
    #https://github.com/RitchieLab/GWAS-QC?tab=readme-ov-file#step-14----imputation-using-topmed-imputation-server


###################
## QC Statistics ##
###################

#Alternative allele frequency > 0.5 sites: 68,904
#Reference Overlap: 98.44 %
#Match: 242,778
#Allele switch: 0
#Strand flip: 0
#Strand flip and allele switch: 0
#A/T, C/G genotypes: 0
#Filter flag set: 0
#Invalid alleles: 0
#Multiallelic sites: 0
#Duplicated sites: 0
#NonSNP sites: 0
#Monomorphic sites: 0
#Allele mismatch: 13
#SNPs call rate < 90%: 0

#Excluded sites in total: 13
#Remaining sites in total: 242,778
#See snps-excluded.txt for details
#Typed only sites: 3,837
#See typed-only.txt for details

#Warning: 4 Chunk(s) excluded: < 20 SNPs (see chunks-excluded.txt for details).
#Warning: 1 Chunk(s) excluded: reference overlap < 50.0% (see chunks-excluded.txt for details).
#Remaining chunk(s): 288

#QC html report:
    #./quality_control/data/genetic_data/quality_control/20_imputation_results/01_qc_reports/qcreport.html
    #The plot show a correlation of 0.924 between the allele frequencies of the reference panel and our data. There are not strange data points going in opposite direction, everything is clean like Ritchie´s plot. Although they have a 0.972 correlation but, still, our correlation is pretty high. Also the highest density is in the upper right side of the plot, like in our case.
    #We have 4326 SNPs with differences in allele frequency between the reference panel and our data. This is a half of what we got in our initial imputation tests. Ritchie nor TOPMed says anything about removing these SNPs, so we are just going to leave this as it is.
    #The important thing is the overall high correlation between the reference panel and our data. Also, we will remove SNPs based on imputation quality (Estimated Imputation Accuracy (R-square)) so if these SNP generated bad imputed genomes, these will be removed.

#Summary:
    #input 246628 but matched 242778, so we lose a total of 3850
    #from there, we lose 13 variants due to allele mismatches
    #therefore, typed only sites are 3837
        #From the total number of SNPS, the type SNPs are SNP not present in the reference panel
        #typed only: https://www.biostars.org/p/446894/
    #this means that the 13 mismatches and the typed sites should not be included in the final count of 242778
    #In other words: 242778+13+3837=246628
    #lose 5 chunks due to low overlap or less than 20 SNPs
    #I have checked this is the case in the corresponding files in 
        #./quality_control/data/genetic_data/quality_control/20_imputation_results/01_qc_reports/
    #The correlation between the reference panel and our data is high.

# endregion






############################
# region MERGE CHROMOSOMES #
############################

print_text("MERGE CHROMOSOMES", header=1)
print_text("create folder", header=2)
run_bash(" \
    mkdir \
        -p \
        ./data/genetic_data/quality_control/21_post_imputation_qc/00_plink_filesets \
")


print_text("decompress the VCF files of each chromosome and convert to plink fileset", header=2)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    for file in ./20_imputation_results/04_uncompressed_vcf_files/chr*.dose.vcf.gz; do \
        gunzip \
            --keep \
            --stdout \
            ${file} >  ./00_plink_filesets/$(basename \"${file%.gz}\")\
    done \
")




print_text("make mergelist.txt file with the list of chromosomes to merge", header=2)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc; \
    seq 2 22 | sed 's/^/chr/' > mergelist.txt \
")
    #seq 2 22 to get the numbers from 2 to 22
    #sed is a stream editor used for text manipulation.
        #s stands for substitution.
        #^ is a regular expression that matches the beginning of a line.
            #I understand this will match any row because you are looking for rows starting with "" and not pattern is used, so you are basically selecting the start of each row
            #to substitute only rows starting with 1, you would do  sed 's/^1/chr1/'
        #chr is the string to be inserted at the beginning of each line.

print_text("check we have the correct versions of plink and plink2. We are using plink1.9 for mergning", header=2)
plink_version=run_bash("plink --version", return_value=True).strip()
plink2_version=run_bash("plink2 --version", return_value=True).strip()
if((plink_version!="PLINK v1.90b7 64-bit (16 Jan 2023)") | (plink2_version!="PLINK v2.00a4.2LM 64-bit Intel (31 May 2023)")):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH PLINK VERSIONS")
else:
    print("PLINK VERSIONS ARE OK")

print_text("Merge chromosome files into 1 bim/bed/bam file THIS TAKES A LONG TIME", header=2)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./20_imputation_results/04_uncompressed_vcf_files/chr1 \
        --merge-list ./21_post_imputation_qc/mergelist.txt \
        --make-bed \
        --out ./21_post_imputation_qc/merged \
")


###por aquii
### FALTA UN Step, tenemos que descomprimir los VCF de dentro y convert a bim

##once merged, REMOVED DECOMPRESSED VCF FILES




#merge chromosomes
    #follow ritchie
    #https://github.com/RitchieLab/GWAS-QC?tab=readme-ov-file#part-4----post-imputation-qc

#inputs in imputation results
    #./quality_control/data/genetic_data/quality_control/20_imputation_results/04_uncompressed_vcf_files/



