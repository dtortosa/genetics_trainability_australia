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
    plink \
        --bfile ../00_plink_filesets/merged_file_sets/merged \
        --geno 0.01 \
        --make-bed \
        --out merged_1_geno; \
    ls ./merged_1_geno.* \
")
    #Want to QC on geno 0.01 (geno filters out all variants with missing call rates exceeding the provided value to be removed)
        #we are using the same value used in Ritchie´s tutorial
    #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.

print_text("check we have not lost any SNP due to missingness", header=3)
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
    plink \
        --bfile ./01_first_qc_step/merged_1_geno \
        --maf 0.05 \
        --hwe 1e-6 midp \
        --mind 0.01 \
        --make-bed \
        --out ./02_second_qc_step/merged_2_geno; \
    ls ./02_second_qc_step/merged_2_geno.* \
")
    #Want to QC on
        #maf 0.05 (filters out all variants with minor allele frequency below the provided threshold)
            #0.05 is the same MAF filter used pre-imputation, which makes sense given our sample size, not small but also not very big. Check line 1835 in 02a_pre_imputation_qc.py for futher details about the reasoning behind the decision. That took into account Ritchie´s explanations for MAF filtering after imputation.
            #If you check you will see that the question is between 0.05 and 0.01, we selected the most stringent option pre-imputation, i.e., only retaining SNPs with a MAF above 0.05, removing all SNPs between 0.05 and 0.01.
            #If you check Figure 8 in Ritchie´s paper, you will see that the power is very low for SNPs with a MAF under 0.05 even if they have a strong odds-ratio (1.3 or more) and the sample size is big (5000 cases and 5000 controls).
                #http://csg.sph.umich.edu/abecasis/cats/gas_power_calculator/index.html
            #It does not make sense to retain these SNPs given their low power and the increase in computational and multiple-test burden. So we are going for 0.05.
        #hwe:
            #filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below p·10-nk, where n is the sample size, and k is 0 if unspecified.
            #--hwe's 'midp' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium. The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data. We recommend its use.
            #I think there is an error in Ritchie´s github here, because they say they are using --hardy to test for HWE deviations but then the code show -hwe 0.05, which is a very stringent filter.
            #In the Ritchie tutorial, they initially say that one of the steps after imputation is "Flag/Remove Variants Out of Hardy-Weinberg Equilibrium (HWE)" but then in the corresponding section, they say that the deviations should not removed as it has been consistently shown that the number of SNPs deviating from HWE at any given significance threshold than would be expected by chance. They could be indeed SNPs under selection and, hence, implicated in traits. 
            #Maares and O´really tutorials suggest to remove cases with P lower than 1e-6, while plinkt suggests a more lower threshold (1e-25) with the idea to retain anything legit but remove clear deviations as these are likely genotyping errors. The point is that THEY ARE REFERING TO PRE-IMPUTATION! However, Augusto applied the 10−6 filter AFTER IMPUTATION.
            #Given that there is likely an error in Ritchie´s Github due to inconsistence between Github and the paper, and the fact that Maares and O´reilly's tutorials suggest to use 10-6 and Augusto used that threshold after imputation, we are going to use that. It is also an approach that does not remove many SNPs, so we are removing likely genotyping/problematic errors while living potentially relevant SNPs.
        #mind 0.01 (filters out al the samples with missing call rates exceeding the provided value to be removed)
            #We are using the same value used pre-imputation. Check line 2468 in 02a_pre_imputation_qc.py for futher details about the reasoning behind the decision
            #In any case, no sample is removed in this step, so we are good.

print_text("check we have the same number of samples and at least 5M SNPs", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    n_snps_before_filtering=$(awk 'END{print NR}' ./01_first_qc_step/merged_1_geno.bim); \
    n_snps_after_filtering=$(awk 'END{print NR}' ./02_second_qc_step/merged_2_geno.bim); \
    n_samples_before_filtering=$(awk 'END{print NR}' ./01_first_qc_step/merged_1_geno.fam); \
    n_samples_after_filtering=$(awk 'END{print NR}' ./02_second_qc_step/merged_2_geno.fam); \
    echo \"We had ${n_snps_before_filtering} SNPs before filtering and ${n_snps_after_filtering} after filtering, i.e., we have lost $((${n_snps_before_filtering}-${n_snps_after_filtering})) \"; \
    if [[ $n_snps_after_filtering -lt 5000000  || $n_samples_after_filtering -ne $n_samples_before_filtering ]]; then \
        echo 'ERROR: FALSE! WE HAVE A PROBLEM WITH THE SECOND QC FILTER'; \
    else \
        echo 'SECOND QC FILTER WORKED FINE'; \
    fi \
")

print_text("check how many SNPs we retain due to setting the HWE threshold at 1-e6 instead of 0.05 and how many are under 1e-6", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/01_first_qc_step/; \
    plink \
        --bfile ./merged_1_geno \
        --hardy \
        --out ./merged_1_geno; \
    snps_between_thresholds=$( \
        awk \
            'BEGIN{FS=\" \"}{ \
                if(NR>1){ \
                    if($9 >= 1e-6 && $9 <=0.05){ \
                        count++; \
                    } \
                } \
            }END{print count}' \
            ./merged_1_geno.hwe; \
    ); \
    snps_under_06=$( \
        awk \
            'BEGIN{FS=\" \"}{ \
                if(NR>1){ \
                    if($9 < 1e-6){ \
                        count++; \
                    } \
                } \
            }END{print count}' \
            ./merged_1_geno.hwe; \
    ); \
    echo \"We have ${snps_between_thresholds} SNPs with a HWE p-value between 1e-6 and 0.05\"; \
    echo \"We have ${snps_under_06} SNPs with a HWE p-value under 1e-6\"; \
    if [[ ${snps_between_thresholds} -lt 300000 || ${snps_under_06} -gt 200 ]]; then \
        echo 'ERROR! FALSE! WE HAVE A PROBLEM WITH THE SECOND QC FILTER, WE HAVE NOT RECOVERED AS MANY SNPS AS EXPECTED BY DECREASING THE HWE THRESHOLD OR THE NUMBER OF SNPS BELOW HWE 1E-6 IS TOO HIGH'; \
    else \
        echo 'OK!'; \
    fi \
")
    #create a hwe file with plink where the 9th column contains the HWE p-value
    #start procesing after the first row to skip the header
    #if the P-value columns is between 0.01 and 0.05, add 1 to the count
    #save the count as a variable and check it is not lower than 300K
    #do the same with the SNPs under 1-6 and check is not larger than 200

print_text("check we do not lose too many SNPs for increasing MAF filter from 0.01 to 0.05", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/01_first_qc_step/; \
    plink \
        --bfile ./merged_1_geno \
        --freq \
        --out ./merged_1_geno; \
    snps_between_01_05=$( \
        awk \
            'BEGIN{FS=\" \"}{ \
                if(NR>1){ \
                    if($5 >= 0.01 && $5 <= 0.05){ \
                        count++; \
                    } \
                } \
            }END{print count}' \
            ./merged_1_geno.frq; \
    ); \
    echo \"We have ${snps_between_01_05} SNPs with a MAF between 0.01 and 0.05\"; \
    if [[ $snps_between_01_05 -gt 3000000 ]]; then \
        echo 'ERROR! FALSE! WE HAVE A PROBLEM WITH THE SECOND QC FILTER, WE LOST MORE THAN 3M OF SNPS DUE TO INCREASE MAF FROM 0.01 TO 0.05'; \
    else \
        echo 'OK! WE LOSE A REASONABLE AMOUNT OF SNPS DUE TO INCREASE MAF FROM 0.01 TO 0.05'; \
    fi \
")
    #create a freq file with plink where the 5th column contains the MAF
    #start procesing after the first row to skip the header
    #if the MAF column is between 0.01 and 0.05, add 1 to the count
    #save the count as a variable and check it is not larger than 3M


print_text("THIRD QC STEP", header=2)
print_text("filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/; \
    plink \
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
    if [[ ${n_samples_after_filtering} -ne ${n_samples_before_filtering} ]]; then \
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
    #Note that Augusto applied an R2 filter under 0.9 instead of 0.7. We have used the most stringent filter option in TOPMed and the recommendation from Ritchie´s tutorial.
        #AUGUSTO POST-IMPTUATION: Once our genomic data were imputed, a second quality control analysis was performed with PLINK 1.9 software [13]. The second quality control exclusion criteria were: low imputation quality (𝑅2<0.9); variants that did not meet the Hardy–Weinberg equilibrium (HWE-P>10−6); and low minor allele frequency (MAF<0.01) [14].

# endregion




print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/quality_control
#chmod +x ./scripts/02d_post_imputation_qc.py
#singularity exec ./singularity_containers/02c_post_input_qc.sif ./scripts/02d_post_imputation_qc.py > 02d_post_imputation_qc.out 2>&1
