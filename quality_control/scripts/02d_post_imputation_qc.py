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
    mkdir -p ./04_batch_effects; \
")


print_text("FIRST QC STEP", header=2)
print_text("filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/01_first_qc_step/; \
    plink \
        --bfile ../00_plink_filesets/merged_file_sets/merged \
        --geno 0.01 \
        --make-bed \
        --out ./merged_1_geno; \
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
            #Given that there is likely an error in Ritchie´s Github due to inconsistence between Github (P=0.05) and the paper, the fact that deviation from HWE can mean impact on the gentoype and the fact that Maares and O´reilly's tutorials suggest to use 1e-6 and Augusto used that threshold after imputation, we are going to use that. We apply the less stringent filter of 1e-6. It is also an approach that does not remove many SNPs, so we are removing likely genotyping/problematic errors while living potentially relevant SNPs.
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
    fi; \
    gzip ./merged_1_geno.hwe; \
")
    #create a hwe file with plink where the 9th column contains the HWE p-value
    #start procesing after the first row to skip the header
    #if the P-value columns is between 0.01 and 0.05, add 1 to the count
    #save the count as a variable and check it is not lower than 300K
    #do the same with the SNPs under 1-6 and check is not larger than 200
    #compress the HWE file.

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
    fi; \
    gzip ./merged_1_geno.frq; \
")
    #create a freq file with plink where the 5th column contains the MAF
    #start procesing after the first row to skip the header
    #if the MAF column is between 0.01 and 0.05, add 1 to the count
    #save the count as a variable and check it is not larger than 3M
    #compress the freq file


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

print_text("See the final number of samples and variants", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/03_third_qc_step/; \
    final_n_samples=$(awk 'END{print NR}' ./merged_3_geno.fam); \
    final_n_snps=$(awk 'END{print NR}' ./merged_3_geno.bim); \
    echo \"We have ${final_n_samples} samples and ${final_n_snps} SNPs after the third QC step\"; \
")

# endregion



########################
# region BATCH EFFECTS #
########################

print_text("BATCH EFFECTS", header=1)
#From Ritchie´s tutorial: Thousands of DNA samples are typically genotyped in a GWAS, which requires partitioning samples into small batches for genotyping in the laboratory (e.g., the set of samples on a 96-well plate). Genotype imputation strategies now provide an opportunity to impute genotypes from multiple platforms to a reference genotyping platform. Therefore, data from different sources can be merged after imputation. Notably, the imputations are not perfect, and the quality and genotyping coverage of input data affect the quality of imputed datasets (Verma et al., 2014). Therefore, it is crucial to evaluate batch effects for the merged datasets. The precise size and composition of the sample batch depends on the array and lab process used. Systematic differences among the compositions of individuals in a batch (i.e., the case to control ratio or ancestry of individuals on plates), and the within-plate accuracy and efficiency, can result in batch effects—apparent associations confounded by batch. The problem is in essence the same problem observed with population stratification—namely, if there is an imbalance of cases and controls on a plate, and there are nonrandom (unknown) biases or inaccuracies in genotyping that differ from plate to plate, spurious associations will result. Ideally, no batch effect will be present be- cause individuals with different phenotypes, sex, ancestry, and other confounders should be plated randomly or imputed, and because modern high-throughput genotyping technology is much more accurate, efficient, and con- sistent than earlier generations of genotyping assays. There are several approaches for exam- ining a dataset for potential batch effects. 
    #not sure if this applies to us, but we are going to do a quick analysis.


print_text("check averages", header=2)
#One simple approach is to calculate the average minor allele frequency and average genotyping call rate across all SNPs for each batch. Gross differences in either of these on any batch can easily be identified.

print_text("first MAF", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/04_batch_effects/; \
    plink \
        --bfile ../03_third_qc_step/merged_3_geno \
        --family \
        --freq \
        --out ./merged_3_geno_batches; \
    head ./merged_3_geno_batches.frq.strat; \
    averages=$( \
        awk \
            'BEGIN{FS=\" \"; count1=0; sum1=0; count2=0; sum2=0}{ \
                if($3==\"combat_ILGSA24-17303\"){ \
                    sum1 += $6; \
                    count1++; \
                } \
                if($3==\"combat_ILGSA24-17873\"){ \
                    sum2 += $6; \
                    count2++; \
                } \
            }END{ \
                if(count1 > 0 && count2 > 0){ \
                    average1=sum1/count1; \
                    average2=sum2/count2; \
                    print average1, average2; \
                }else { \
                    exit 1 \
                } \
            }' \
            ./merged_3_geno_batches.frq.strat \
    ); \
    read average1 average2 <<< ${averages}; \
    difference=$(echo \"${average1} - ${average2}\" | bc -l); \
    echo \"MAF-Avg1: ${average1}; MAF-Avg2:${average2}; MAF-Diff:${difference}\"; \
    diff_check=$(echo \"${difference} > 0.001\" | bc -l); \
    if [[ ${diff_check} -eq 1 ]]; then \
        echo \"ERROR! FALSE! MAF AVERAGES ARE VERY DIFFERENT\"; \
    fi; \
    gzip --force ./merged_3_geno_batches.frq.strat; \
")
    #calculate a freq file in a stratified within each batch, i.e., family
    #calculate the average MAF for each batch using awk
        #BEGIN by creating two empty counts and two empty sums
        #for each row extract the MAF and add it to the sum of the corresponding batch, also count the row for the corresponding count
        #at the END, divided sum by count to get the average MAF of each batch
    #split the variable with the two averages
        #read: it is used to read input from a file descriptor or standard input and assign it to one or more variables.
        #average1 and average2 are the variables that will store the values read from the input. In this case, average1 will store the first value, and average2 will store the second value.
        #"<<<": This is a here-string operator. It allows you to pass a string directly to the read command as if it were input from a file or standard input. The string is treated as a single line of input.
    #calculate the difference
        #print the math calculation and input it to bc in math mode ("-l")
    #check if the difference is higher than 0.001 using bc again
    #if the output of bc is 1, means the difference is higher than 0.001, so we have a problem.

print_text("then call rate", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/21_post_imputation_qc/04_batch_effects/; \
    plink \
        --bfile ../03_third_qc_step/merged_3_geno \
        --family \
        --missing \
        --out ./merged_3_geno_batches; \
    head ./merged_3_geno_batches.lmiss; \
    averages=$( \
        awk \
            'BEGIN{FS=\" \"; count1=0; sum1=0; count2=0; sum2=0}{ \
                if($3==\"combat_ILGSA24-17303\"){ \
                    sum1 += $7; \
                    count1++; \
                } \
                if($3==\"combat_ILGSA24-17873\"){ \
                    sum2 += $7; \
                    count2++; \
                } \
            }END{ \
                if(count1 > 0 && count2 > 0){ \
                    average1=sum1/count1; \
                    average2=sum2/count2; \
                    print average1, average2; \
                }else { \
                    exit 1 \
                } \
            }' \
            ./merged_3_geno_batches.lmiss \
    ); \
    read average1 average2 <<< ${averages}; \
    difference=$(echo \"${average1} - ${average2}\" | bc -l); \
    echo \"MISSING-Avg1: ${average1}; MISSING-Avg2:${average2}; MISSING-Diff:${difference}\"; \
    diff_check=$(echo \"${difference} > 0.001\" | bc -l); \
    if [[ ${diff_check} -eq 1 ]]; then \
        echo \"ERROR! FALSE! MISSING RATE AVERAGES ARE VERY DIFFERENT\"; \
    fi; \
    gzip --force ./merged_3_geno_batches.lmiss; \
    rm ./merged_3_geno_batches.imiss; \
")
    #the same than the previous awk script but using an lmiss file and the F_MISS column
    #lmiss file
        #CHR: Chromosome code
        #SNP: Variant identifier
        #CLST: Cluster identifier. Only present with --within/--family.
        #N_MISS: Number of missing genotype call(s), not counting obligatory missings or het. haploids
        #N_CLST: Cluster size (does not include nonmales on chrY). Only present with --within/--family.
        #N_GENO: Number of potentially valid call(s)
        #F_MISS: Missing call rate


print_text("check about the fact that we have NO missing genotypes after imputation", header=2)
#I noticed that we do not have any missing genotype when I was doing the batch checks
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./19_topmed_prep/01_second_step/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated \
        --missing \
        --out ./21_post_imputation_qc/04_batch_effects/missing_before_imputation; \
    rm ./21_post_imputation_qc/04_batch_effects/missing_before_imputation.imiss; \
")
    #create a missing report using the last combine plink fileset before splitting across chromosomes and sending to the imputation server
average_missing_before_imput=pd.read_csv("./data/genetic_data/quality_control/21_post_imputation_qc/04_batch_effects/missing_before_imputation.lmiss", sep="\s+").loc[:,"F_MISS"].mean()
    #sep="\s+" because some columns are separated by 1 space while others are separated by several spaces.
if (average_missing_before_imput>0.00039891746854450993):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE SNP MISSING RATE BEFORE IMPUTATION")
else:
    print("MISSING RATE BEFORE IMPUTATION WAS: " + str(average_missing_before_imput))
    #We have NO missing genotypes. This makes sense because the average missing rate across SNPs in the plink fileset used as input for imputation already had a very low missing rate (0.000398917 being 1 the max possible). So it makes sense that imputation has been able to fill the few remaining gaps in the SNPs used as input. Indeed, Augusto et al., did NOT apply a missingness filter after imptuation in their paper. So all this makes sense.
        #Once our genomic data were imputed, a second quality control analysis was performed with PLINK 1.9 software [13]. The second quality control exclusion criteria were: low imputation quality (𝑅2<0.9); variants that did not meet the Hardy–Weinberg equilibrium (𝐻𝑊𝐸−𝑃>10−6); and low minor allele frequency (𝑀𝐴𝐹<0.01) [14].
        #This behavior is expected. Genetic imputation not only predicts the genotypes of SNPs that were not originally genotyped but also fills in missing values for SNPs that were genotyped but had missing data points. This process leverages the linkage disequilibrium (LD) patterns and reference panels to infer the most likely genotypes for the missing data[1](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-5-404)[2](https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-5-S7-P61).
        #By using imputation, you can achieve a more complete dataset, which is particularly useful for downstream analyses like genome-wide association studies (GWAS). The imputation process improves the accuracy and power of these analyses by providing a more comprehensive set of genotypic data[2](https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-5-S7-P61)[3](https://www.nature.com/articles/5201988.pdf).


print_text("check for association between the batch assignation and SNP frquencies", header=2)
#Another method involves coding case/control status by batch followed by running the GWAS analysis testing each batch against all other batches. For example, the status of all samples on batch 1 will be coded as case, while the status of every other sample is to be coded control. A GWAS analysis is performed (e.g., using the --assoc option in PLINK), and both the average p-value and the number of results significant at a given threshold (e.g., p <1 × 10-4 ) can be recorded. SNPs with low minor allele frequency (i.e., <5%) should be removed before this analysis is performed to improve the stability of test statistics. This procedure should be repeated for each batch in the study. If any single batch has many more or many fewer significant re- sults or has an average p-value <0.5 (under the null, the average p-value will be 0.5 over many tests), then this batch should be further inves- tigated for genotyping, imputation, or compo- sition problems. If batch effects are present, methods like those employed for population stratification (e.g., genomic control) may be used to mitigate the confounding effects.

print_text("calculate the association between snp frequencies and batches ", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/;\
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{ \
            pheno=0; \
            if($1 == \"combat_ILGSA24-17303\"){ \
                pheno=1; \
            } else if ($1 == \"combat_ILGSA24-17873\") { \
                pheno=2; \
            } \
            print $1, $2, pheno; \
        }' \
        ./21_post_imputation_qc/03_third_qc_step/merged_3_geno.fam > ./21_post_imputation_qc/04_batch_effects/case_control.tsv; \
    plink \
        --bfile ./21_post_imputation_qc/03_third_qc_step/merged_3_geno \
        --pheno ./21_post_imputation_qc/04_batch_effects/case_control.tsv \
        --assoc \
        --out ./21_post_imputation_qc/04_batch_effects/merged_3_geno_batch_assoc; \
")
    #from the last fam file, get the family and within family ID, also create a case/control variable. Initially the variable is 0 (no pheno), and then a value is assigned: 1 to samples of the first batch (control) and 2 to samples of second batch (cases). This willl be the input for --pheno
        #Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
    #--pheno causes phenotype values to be read from the 3rd column of the specified space- or tab-delimited file, instead of the .fam or .ped file. The first and second columns of that file must contain family and within-family IDs, respectively.
    #--assoc:
        #Given a case/control phenotype, --assoc writes the results of a 1df chi-square allelic test to plink.assoc (or .assoc.fisher with 'fisher'/'fisher-midp').

print_text("then calculate average p-value and number of signifncant SNPs ", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/;\
    p_stats=$( \
        awk \
            'BEGIN{FS=\" \"; p_count=0; p_sum=0; p_signi_count=0}{ \
                if(NR>1){ \
                    p_sum += $9; \
                    p_count++; \
                    if($9 < 1e-4){ \
                        p_signi_count++; \
                    } \
                } \
            }END{ \
                if(p_count>0){ \
                    p_average=p_sum/p_count \
                } else { \
                    exit 1; \
                } \
                print p_average, p_signi_count; \
            }' \
            ./21_post_imputation_qc/04_batch_effects/merged_3_geno_batch_assoc.assoc; \
    ); \
    read p_average p_signi_count <<< ${p_stats}; \
    echo \"The average p-value is ${p_average}, while the number of SNPs with P<1e-4 is ${p_signi_count}\"; \
    if [[ $(echo \"${p_average} < 0.05\" | bc -l) -eq 1 || ${p_signi_count} -gt 1000 ]]; then \
        echo \"ERROR! FALSE! WE HAVE SIGIFNICANT DIFFERENCES IN THE ALLELE FREQUENCIES OF THE SNPS BETWEEN BATCHES\"; \
    else \
        echo \"OK! WE DO NOT HAVE SIGIFNICANT DIFFERENCES IN THE ALLELE FREQUENCIES OF THE SNPS BETWEEN BATCHES\"; \
    fi; \
")
    #for each row of the assoc results
        #get the p-value and sum it to the previous values, add 1 to the count to calculate the average later
        #also add 1 to the count if the p-value is lower than 1e-4
        #at the END, calculate teh average p-value dividing the sum of p-values by the total number of p-values (count)
        #print the average and the number of signifcant p-values
    #split the two results into two different variables
    #return error if p_average is lower than 0.05 or the number of very significant SNPs is above 1000

# endregion


print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/quality_control
#chmod +x ./scripts/02d_post_imputation_qc.py
#singularity exec ./singularity_containers/02c_post_input_qc.sif ./scripts/02d_post_imputation_qc.py > ./02d_post_imputation_qc.out 2>&1

