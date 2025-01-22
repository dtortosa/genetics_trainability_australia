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

##################
# import modules #
##################

import pandas as pd



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






###################
# region APPLY QC #
###################

print_text("APPLY QC", header=1)
print_text("prepare folders", header=2)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    mkdir -p ./01_first_qc_step; \
    mkdir -p ./02_second_qc_step; \
    mkdir -p ./03_third_qc_step; \
")


print_text("FIRST QC STEP", header=2)
print_text("filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/01_first_qc_step/; \
    plink2 \
        --bfile ../00_plink_filesets/merged_file_sets/merged \
        --geno 0.01 \
        --make-bed \
        --out merged_1_geno; \
    ls ./merged_1_geno.* \
")
    #Want to QC on geno 0.01 (geno filters out all variants with missing call rates exceeding the provided value to be removed)
        #we are using the same value used in Ritchie´s tutorial
    #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.

print_text("check we have not lost any sample", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    n_snps_before_filtering=$(awk 'END{print NR}' ./00_plink_filesets/merged_file_sets/merged.bim); \
    n_snps_after_filtering=$(awk 'END{print NR}' ./01_first_qc_step/merged_1_geno.bim); \
    if [[ $n_snps_after_filtering -ne $n_snps_before_filtering ]]; then \
        echo 'ERROR: FALSE! WE HAVE A PROBLEM WITH THE FIRST QC FILTER'; \
    else \
        echo 'FIRST QC FILTER WORKED FINE'; \
    fi \
")


print_text("SECOND QC STEP", header=2)
print_text("filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    plink2 \
        --bfile ./01_first_qc_step/merged_1_geno \
        --mind 0.01 \
        --maf 0.05 \
        --hwe 0.05 \
        --make-bed \
        --out ./02_second_qc_step/merged_2_geno; \
    ls ./02_second_qc_step/merged_2_geno.* \
")
    #Want to QC on
        #mind 0.01 (filters out al lthe samples with missing call rates exceeding the provided value to be removed)
        #maf 0.05 (filters out all variants with minor allele frequency below the provided threshold)
        #hwe: filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below p·10-nk, where n is the sample size, and k is 0 if unspecified.

print_text("check we have the same number of samples and at least 5M SNPs", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    n_snps_after_filtering=$(awk 'END{print NR}' ./02_second_qc_step/merged_2_geno.bim); \
    n_samples_before_filtering=$(awk 'END{print NR}' ./01_first_qc_step/merged_1_geno.fam); \
    n_samples_after_filtering=$(awk 'END{print NR}' ./02_second_qc_step/merged_2_geno.fam); \
    if [[ $n_snps_after_filtering -lt 5000000  || $n_samples_after_filtering -ne $n_samples_before_filtering ]]; then \
        echo 'ERROR: FALSE! WE HAVE A PROBLEM WITH THE SECOND QC FILTER'; \
    else \
        echo 'SECOND QC FILTER WORKED FINE'; \
    fi \
")


print_text("THIRD QC STEP", header=2)
print_text("filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    plink2 \
        --bfile ./02_second_qc_step/merged_2_geno \
        --mind 0.01 \
        --make-bed \
        --out ./03_third_qc_step/merged_3_geno; \
    ls ./03_third_qc_step/merged_3_geno.* \
")

print_text("check we do not lose any sample", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    n_samples_before_filtering=$(awk 'END{print NR}' ./02_second_qc_step/merged_2_geno.fam); \
    n_samples_after_filtering=$(awk 'END{print NR}' ./03_third_qc_step/merged_3_geno.fam); \
    if [[ $n_samples_after_filtering -ne $n_samples_before_filtering ]]; then \
        echo 'ERROR: FALSE! WE HAVE A PROBLEM WITH THE THIRD QC FILTER'; \
    else \
        echo 'THIRD QC FILTER WORKED FINE'; \
    fi \
")
    #In the previous steps no sample was removed, so the MAF and HWE of any SNP should be affected, rembember that removing a SNP becuase it has low MAF does not influence the MAF of other SNP.
    #However, as we have removed millions of SNPs in the pervious step, it is possible that the proportion of missing could increase in some sample, so we repeat the --mind filter again to check that no sample is removed with the reduced set of SNPs.


print_text("check that we do not have ANY snp with an RSQ equal or lower than 0.3", header=3)
run_bash(" \
    for chrom_numb in {1..22}; do \
        echo \"Starting chromosome ${chrom_numb}\"; \
        awk \
            'BEGIN{FS=\";\"; start=\"no\"}{ \
                if (start==\"yes\"){ \
                    split($5, a, \"=\"); \
                    if (a[2] <= 0.3){ \
                        exit 1; \
                    } \
                }; \
                if (/^#CHROM/){ \
                    start=\"yes\"; \
                } \
            }' \
            <(gunzip --keep --stdout ./data/genetic_data/quality_control/20_imputation_results/04_uncompressed_vcf_files/chr${chrom_numb}.info.gz); \
        echo \"Chromosome ${chrom_numb} is OK\"; \
    done; \
")
    #For each chromosome
        #decompres the INFO file and send it to stdout on the fly, being used as input for AWK
        #split by ";" so we can separate the INFO fields which are separated by ";". Also create a variable called start being equal to "no".
        #do stuff only if start is "yes", and that would be only the case after we pass the row with the names of the columns of the VCF file, i.e., we skip the header.
        #in each row, split the 5th fields (R2=...) using "=" and save the two elements, R2 and the number in "a", if a2, i.e., the value of R2 is equal or lower than 0.3, then stop the execution.




#CHECK RITCHIE PAPER TO SEE IF MISSING STEP
    #DECIDE WHETHER TO USE 0.05 OR 0.01 FOR MAF... IT IS PERCENTAGE?


#CHECK THE CHROSOMOE CODES! WE HAVE USED PLINK2 AND THAT CHANGES IT



# endregion