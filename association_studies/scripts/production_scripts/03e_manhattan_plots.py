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



#############################################
######## PLOT TO VISUALIZE PRS RESULTS ######
#############################################

#This script will be used to create plot that help to visualize the results from the PRSs.

#This and previous scripts are based on the following tutorials:
    #1. Quality Control Procedures for Genome-Wide Association Studies
        #https://github.com/RitchieLab/GWAS-QC
        #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
    #2. Genome-wide association studies
        #https://www.nature.com/articles/s43586-021-00056-9
    #3. Data Management and Summary Statistics with PLINK
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
    #4. Genomics Boot Camp
        #https://genomicsbootcamp.github.io/book/
    #5. Tutorial: a guide to performing polygenic risk score analyses
        #https://www.nature.com/articles/s41596-020-0353-1
        #https://choishingwan.github.io/PRS-Tutorial/
    #6. Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
    #7. Omics Data Preprocessing for Machine Learning: A Case Study in Childhood Obesity
        #https://www.mdpi.com/2073-4425/14/2/248
    #8. LDAK
        #https://dougspeed.com/
        #https://www.nature.com/articles/s41467-021-24485-y




#######################################
# region INITIAL ANOTATIONS AND STEPS #
#######################################


###########
# imports #
###########

import pandas as pd
import numpy as np
from itertools import product
import re


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
    elif header==4:
        print("\n# " + text + " #")
print_text("checking function to print nicely: header 1", header=1)
print_text("checking function to print nicely: header 2", header=2)
print_text("checking function to print nicely: header 3", header=3)
print_text("checking function to print nicely: header 4", header=4)



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






###########################################
# region RUN THE PRSs ON THE FULL DATASET #
###########################################

#Both elastic and linear across thresholds

#We are going to use small dataset only. We have already seen no meaninful differences between small and large when running the PRSs in 100 training/evaluation sets

#response_variable="beep_change"
def prs_calc(response_variable):

    print_text(f"For phenotype {response_variable}, and the small dataset of covariates", header=1)
    print_text("initial preparations", header=2)
    print_text("create folders for results", header=3)
    run_bash(" \
        mkdir -p \
            ./results/final_results/analysis_full_data/" + response_variable + "; \
    ")

    print_text("load the phenotype data", header=3)
    pheno_subset_transform = pd.read_csv( \
        "./data/full_set_transform/small_set_predictors/" + response_variable + "/" + response_variable + "_full_set_transform.tsv", \
        header=0, \
        sep="\t" \
    )
    print(pheno_subset_transform)

    print_text("specify the covariates", header=3)
    selected_covariates = pheno_subset_transform.columns[~pheno_subset_transform.columns.isin(["family_id", "AGRF code", response_variable])]
    print(selected_covariates)
    #check we have correct covariates
    total_list_covariates = ["Age", "sex_code", "Week 1 Body Mass", "Week 1 Beep test", "Week 1 Distance (m)", "Week 1 Pred VO2max"] + [f"PCA{i}" for i in range(1,21)]
    if(sum([1 for cov in selected_covariates if cov not in total_list_covariates])!=0):
        raise ValueError("ERROR: FALSE! WE HAVE COVARIATES THAT ARE NOT IN THE TOTAL LIST")

    print_text("decompress bim and bed files with all sample generated after NA cleaning", header=3)
    run_bash(" \
        cd ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/; \
        gunzip \
            --keep \
            --force \
            ./" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing.bed.gz \
            ./" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing.bim.gz; \
    ")

    print_text("create dict to change names of covariates that are problematic for LDAK", header=3)
    dict_change_names={
        "family_id": "FID",
        "AGRF code": "IID",
        "Age": "age",
        "Week 1 Body Mass": "week_1_weight",
        "Week 1 Beep test": "week_1_beep",
        "Week 1 Distance (m)": "week_1_distance",
        "Week 1 Pred VO2max": "week_1_vo2"
    }

    print_text("create dict to change names of responses", header=3)
    dict_change_responses={
        "Change in Body Mass": "weight_change",
        "Change in Beep Test": "beep_change",
        "Change in Distance (m)": "distance_change",
        "Change in VO2max": "vo2_change"
    }


    print_text("prepare LDAK inputs", header=2)
    #save the response variable after changing the names
    pheno_subset_transform[["family_id", "AGRF code", response_variable]].rename(columns=dict_change_names).to_csv( \
        "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_response.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )
        #In LDAK, phenotype files should be in PLINK format. The first two columns should provide sample IDs, with subsequent columns providing values for each phenotype. The files can contain a header row, but then the first two elements must be named either FID & IID or ID1 & ID2.
        #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype). We do not have NAs anyways.
        #Also, we are going to have just 1 response variable, so no problem with filling NA cases with the average. We have just 1 response and it is free of NAs.

    #save the covariates that are factors
    if ("sex_code" in selected_covariates):

        #save sex_code
        pheno_subset_transform[["family_id", "AGRF code", "sex_code"]].rename(columns=dict_change_names).to_csv( \
            "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_factors.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )
            #In LDAK, phenotype files should be in PLINK format. The first two columns should provide sample IDs, with subsequent columns providing values for each phenotype. The files can contain a header row, but then the first two elements must be named either FID & IID or ID1 & ID2.
            #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype). We do not have NAs anyways.

    #save the covariates that are continuous
    selected_covariates_cont = [cov for cov in selected_covariates if cov != "sex_code"]
    pheno_subset_transform[["family_id", "AGRF code"] + selected_covariates_cont].rename(columns=dict_change_names).to_csv( \
        "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_cont.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )
        #In LDAK, phenotype files should be in PLINK format. The first two columns should provide sample IDs, with subsequent columns providing values for each phenotype. The files can contain a header row, but then the first two elements must be named either FID & IID or ID1 & ID2.
        #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype). We do not have NAs anyways.

    print_text("calculate elastic PRS", header=2)
    #Remember we are not using anymore training/evaluation sets, we are using the full dataset so the --elastic is already calculating the PRS in the input dataset, no need to use --calc-scores to predict the phenotype in new samples based on the PRS
    print_text("run the PRS with --elastic", header=3)
    run_bash(" \
        ldak6.1.linux \
          --elastic ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_elastic \
          --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
          --pheno ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_response.tsv \
          --covar ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_cont.tsv \
          --factors ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_factors.tsv \
          --LOCO NO \
    ")
        #Dr. Speed: Yes, your steps sound right - perform QC, then use the PLINK files with --elastic (adding --LOCO NO, because you only care about prediction) and then --calc-scores
            #We DO NOT CARE ABOUT HERIABILITY BECAUSE WE CANNOT CALCULATE IT PROPERLY!! JUST PREDICTION.
        #Here we explain how to construct an elastic net prediction model. These instructions assume you are analysing individual-level data.
            #So I understand we are combining here ridge and lasso, which is the elastic net.
        #construct an elastic net prediction model by analysing individual-level data (if instead you are analysing summary statistics, you should use MegaPRS).
            #--bfile/--speed <datastem> or --bgen <datafile> - to specify the genetic data files (see File Formats). If your genetic data are in a different format, you should first Make Data.
                #LDAK accepts genetic data in many formats. When analyzing SNP data, it is easiest to use bed format (binary PLINK) which is designed to store hard-coded SNP genotypes (all values must be 0, 1, 2 or NA).
                #LDAK considers the last two columns of the bim file are the A1 and A2 alleles, respectively. During, QC, we checked for strand issues and flip alleles if needed. Then, we select REF and ALT from the VCF files generated from the imputation and converted them to plink format using plink.
                #Note that plink swaps REF by ALT if the REF is less frequent, but this should not be a problem because we are not going to use the alleles, just PRS. Anyways, if REF is "A" and ALT is "T", as long as we are using the foward strand, then we can compare with other studies, because "T" is indeed "T", so we could see if that allele associate with a given phenotype in the same way than other studies. And we have ensured this during the quality control also removing problematic palindromic cases...
                    #If needed, you can use "--ref-allele" and input the REF allele list from the source VCF files to correct this.
            #--pheno <phenofile> - to specify phenotypes (in PLINK format; see above). Samples without a phenotype will be excluded. If <phenofile> contains more than one phenotype, specify which should be used with --mpheno <integer>, or use --mpheno ALL to analyze all phenotypes.
                #We are using just one phenotype, so we do not need to specify --mpheno
            #You can use --covar <covarfile> or --factors <factorfile> to provide quantitative or categorical covariates (in PLINK format; see above); the phenotype will be regressed on these prior to estimating effect sizes.
                #Note that if a categorical covariate has U unique values, LDAK will (internally) replace it with U-1 indicator variables (LDAK will give an error if the total number of indicator variables is greater than half the sample size).
                #For example, sex code is 1 and 2, so in LDAK would be 0 and 1.
            #--LOCO NO - to tell LDAK to focus on creating the genome-wide prediction model (instead of creating leave-one-chromosome-out models for use with LDAK-KVIK).
                #This was the recommendation of Dr. Speed, and we are not using the LDAK-KVIK model anyways.
            #By default, LDAK will estimate the heritability and the power parameter alpha; to instead specify their values use --power <float> and --her <float> (note that if you use --her, you must also use --power).
            #By default, LDAK will use 90%/10% cross-validation to determine suitable prior distribution parameters. You can change the fraction of test samples uing --cv-proportion <float>,  specify the test samples using --cv-samples <cvsampsfile>, or turn off cross-validation, using --skip-cv YES (LDAK will then output multiple models, each trained using 100% of samples).
                #So it uses CV internally, which is great. We do CV to tune the hyperparameters of the model, and then we use the best model to predict the test set.
                #This makes each run to have slighlty different results
                    #When using elastic, the cross-validation samples are randomly picked, and so it is possible different parameters are picked each time (there might also be other small differences between versions, such as maybe I changed the tolerance or something)
            #By default, LDAK will assign all predictors (i.e., SNPs) weighting one (equivalent to using --ignore-weights YES). If you prefer to provide your own weightings, use --weights <weightsfile> or --ind-hers <indhersfile> (note that if using --ind-hers, you can not use --her or --power).
                #we are just going to use the default weights. To my understanding, LDAK will give different weights to the SNPs depending on their MAF, so we do not need to worry about that.
                #Doug: Yes, that is right (no need to use weights, they will automatically be set to one, and ldak automatically estimates alpha, the MAF scaling parameger)
            #You can use --keep <keepfile> and/or --remove <removefile> to restrict to a subset of samples, and --extract <extractfile> and/or --exclude <excludefile> to restrict to a subset of predictors (for more details, see Data Filtering).
            #output:
                #The estimated prediction model is saved in <outfile>.effects. Usually, this file has five columns, providing the predictor name (SNP), its A1 and A2 alleles, the average number of A1 alleles, then its estimated effect (relative to the A1 allele). If you used --skip-cv YES, there will be effect sizes for each of the different prior parameters. This file is ready to be used for Calculating Scores (i.e., to predict the phenotypes of new samples).
                #Regarding phenotypes in .PRS file:
                    #Yes, these are the phenotypes after regressing on covariates (e.g., they will have mean zero)
                    #I understand this is the previous step where LDAK remove the impact of covariates to then analyze the impact of the SNPs. So, in the .prs file, the phenotype is already adjusted for the covariates.
            #https://dougspeed.com/elastic-net/

    print_text("do the calculation of the PRS in the same set of samples", header=3)
    print_text("first using the phenotype as input", header=3)
    run_bash(" \
        ldak6.1.linux \
            --calc-scores ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_elastic_prs_calc_with_pheno \
            --scorefile ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_elastic.effects \
            --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
            --pheno ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_response.tsv \
            --power 0 \
    ")
        #calculate linear combinations of predictor values (the linear projection of genetic data onto predictor effect sizes). These are most commonly used for creating polygenic risk scores and testing the performance of prediction models. Please note that if you include a phenotype when calculating scores, then the resulting profile file can be used with Jackknife (in order to compute additional measures of accuracy and corresponding estimates of precision).
        #--scorefile <scorefile> - to provide the predictor effect sizes. The score file should have at least five columns. The first 4 columns provide the name of each predictor, its two alleles (test allele followed by reference allele), then its centre (the mean number of test alleles); the remaining columns provide sets of predictor effect sizes. The file should have a header row, whose first element must be "Predictor" or "SNP". Note that if centre is "NA" for a predictor, then LDAK will centre values based on the mean number of test alleles in the genetic data.
        #–bfile/–gen/–sp/–speed <prefix> or --bgen <datafile> - to specify genetic data files (see File Formats)
        #-pheno <phenofile> - to specify phenotypes (in PLINK format; see above).
            #Typically, the sets of effect sizes in the score file will correspond to different prediction models, and we are interested in measuring their accuracy. We have phenotypes for samples in the genetic data (i.e., you have individual-level validation data), in which case you should provide these using --pheno <phenofile>. LDAK will then calculate the correlation between scores and phenotypic values for the samples in the genetic data.
        #Sometimes you will have two sets of prediction models, for example, one computed using training samples and one computed using all samples. You can then use --scorefile <scorefile> to provide the first set of prediction models, --pheno <phenofile> or --summary <sumsfile> to provide phenotypic values or summary statistics, and --final-effects <finaleffectsfile> to provide the second set of prediction models. LDAK will save to <outfile>.effects.best the prediction model from the second set that corresponds to the most accurate model from the first set (for example, if Model 5 from the first set has highest correlation, LDAK will save Model 5 from the second set).
        #--power <float> - to specify how predictors are scaled (see below). Usually, the score file contains raw effects, so you should use --power 0.
            #When power equals one, predictors are divided by their expected standard deviation (i.e., are standardised). Therefore, you should use --power=-1 when the score file contains standardised effect sizes. By contrast, when power equals zero, predictors are no longer scaled. Therefore, you should use --power=0 when the score file contains raw effect sizes.
            #Our effect sizes are NOT standardized!!!
                #Doug: I think the latest version does not require power, but it should be zero, because you will use scorefile to provide raw effect sizes.
        #I guess it is recommended to use the default options for --calc-scores and NOT use "--hwe-stand NO"
            #Doug: Go with default, but normally makes no difference
        #output:
            #The profile scores will be saved in <outfile>.profile, the (estimated) correlation between scores and phenotypic values will be saved in <outfile>.cors.
        #https://dougspeed.com/profile-scores/

    print_text("then without phenotype, so we can check the PRS calculation of input reposne and covariates and hence we do not need to add covariates", header=3)
    run_bash(" \
        ldak6.1.linux \
            --calc-scores ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_elastic_prs_calc_without_pheno \
            --scorefile ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_elastic.effects \
            --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
            --power 0 \
    ")

    print_text("check the PRS is the same", header=3)
    run_bash(" \
        cd ./results/final_results/analysis_full_data/" + response_variable + "; \
        n_differ=$( \
            awk \
                'BEGIN{FS=\"\t\"}{ \
                    if(NR==FNR){ \
                        a[$5]; \
                        next \
                    } \
                    if(!($5 in a)){ \
                        print $0 \
                    } \
                }' \
                " + response_variable + "_small_set_predictors_set_elastic_prs_calc_with_pheno.profile \
                " + response_variable + "_small_set_predictors_set_elastic_prs_calc_without_pheno.profile | \
            wc -l \
        ); \
        if [[ $n_differ -ne 0 ]]; then \
            echo \"ERROR: FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLES FOR ELASTIC\"; \
        fi; \
    ")
        #do operations in the first file, i.e., when NR=FNR
            #NR: The total number of records (lines) processed so far across all input files.
            #FNR: The number of records (lines) processed in the current file.
            #when they are the same it means we are processing the first file, and when they are different it means we are processing the second file.
        #the operation to do is:
            #a[$5]: The fifth field (PRS for not thresholding) of the current line is used as an index in the array a. This means that the script is storing the values of the fifth field from the first file in the array a.
            #the "next" statement is used to skip the remaining commands in the current iteration of the loop and move to the next line of the input file. This makes sense because if NR=FNR, then we are processing the first file and we do not need to run the code dedicated to the second file.
        #the code after the "next" statement is executed only when we are processing the second file (i.e., when NR!=FNR). In this case, it checks if the value of the fifth field ($5) is not present in the array a. If it is not present, it prints the entire line.
        #numbre of problematic lines is stored in the variable n_differ, and if it is not equal to 0, then we have a problem.


    print_text("calculate a PRS using the classical approach", header=2)
    #Dr. Speed: You may wish to do a classical PRS just for interest / comparison, in which case, you can run --linear (on the training samples), then the file with suffix .score gives classical PRS .corresponding to 7 p-value thresholds (that you can then provide to --calc-scores)
    print_text("calculate linear association between SNPs and the phenotype", header=3)
    run_bash(" \
        ldak6.1.linux \
            --linear ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_linear_raw \
            --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
            --pheno ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_response.tsv \
            --covar ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_cont.tsv \
            --factors ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_factors.tsv \
            --permute NO \
    ")
        #perform one-predictor-at-a-time analysis, using linear regression (classical model, i.e., not mixed model).
            #--bfile/--gen/--sp/--speed <datastem> or --bgen <datafile> - to specify the genetic data files (see File Formats).
            #--pheno <phenofile> - to specify phenotypes (in PLINK format). Samples without a phenotype will be excluded. If <phenofile> contains more than one phenotype, specify which should be used with --mpheno <integer>, or use --mpheno ALL to analyze all phenotypes.
            #You can use --covar <covarfile> or --factors <factorfile> to provide quantitative or categorical covariates (in PLINK format) as fixed effect in the regression.
            #If you add --spa-test YES, LDAK will recompute p-values for the most associated predictors using a SaddlePoint Approximation.
            #It remains possible to  add  --grm <kinfile>, in which case LDAK will perform mixed-model linear regression. However, we instead recommend using LDAK-KVIK, which is usually faster and more powerful.
            #To perform a within-family analysis, add --families YES (for more details of this analysis, see Howe et al.). Note that LDAK will infer families based on the 1st column of the fam file (the FID).
            #To perform a trio analysis, add --trios YES. Note that LDAK will infer trios based on the 3rd and 4th column of the fam file (the PID and MID).
            #To perform weighted linear regression, use --sample-weights <sampleweightfile>. The file <sampleweightfile> should have three columns, where each row provides two sample IDs followed by a positive float. Note that by default, LDAK will use the sandwich estimator of the effect size variance (see this page for an explanation); to instead revert to the standard estimator of variance, add --sandwich NO.
            #If you add --permute YES - the phenotypic values will be shuffled. This is useful if wishing to perform permutation analysis to see the distribution of p-values or test statistics when there is no true signal.
                #SO WHEN YOU SET THIS AS YES, YOU ARE GETTING THE P-VALUE WHEN THERE IS NO SIGNAL
                #Doug confirmed this.
            #output:
                #When performing a standard analysis, LDAK produces five output files: <outfile>.assoc contains the main results; <outfile>.summaries contains summary statistics (in the format required for use with SumHer, and MegaPRS); <outfile>.pvalues contains p-values (useful if you wish to Thin Predictors); <outfile>.coeff contains estimates of the fixed effects; <outfile>.score contains simple prediction models corresponding to six different p-value thresholds.
            #https://dougspeed.com/single-predictor-analysis/

    print_text("define p-values cut-offs for thresholding", header=3)
    #Doug: Clumping normally helps, so I would use --thin-tops (e.g., window-kb 1000, window-prune .2)
    print_text("load the p-values", header=4)
    training_linear_raw_pvalues = pd.read_csv( \
        "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_linear_raw.pvalues", \
        sep="\t", \
        header=0 \
    )

    print_text("get the minimum p-value", header=4)
    min_p_value = training_linear_raw_pvalues["P"].min()

    print_text("make a list of thresholds that are above the minimum p-value", header=4)
    list_thresholds_raw = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 5e-8]
        #if no SNP is under 0.001, then it does not make sense to filter by 0.001 as no SNPs would be left. Indeed this generates an error in LDAK.
    list_thresholds = [i for i in list_thresholds_raw if i>min_p_value]

    print_text("obtain scores in a reduced set of SNPs after thresholding and cumpling", header=3)
    #threshold=list_thresholds[0]
    for threshold in list_thresholds:

        print_text("perform thresholding considering the threshold " + str(threshold) + " and then perform clumping", header=4)
        run_bash(" \
            mkdir \
                -p \
                ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "; \
            ldak6.1.linux \
                --thin-tops ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_small_set_predictors_set_linear_clump_thresholding_" + str(threshold) + "_predictors \
                --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
                --window-prune 0.2 \
                --window-kb 1000 \
                --pvalues ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_linear_raw.pvalues \
                --cutoff " + str(threshold) + " \
        ")
            #Doug: Clumping normally helps, so I would use --thin-tops (e.g., window-kb 1000, window-prune .2)
            #--bfile/--gen/--sp/--speed <datastem> or --bgen <datafile> - to specify the genetic data files (see File Formats).
                #"--thin-tops" is first thresholding by considering only those SNPs below the threshold indicated in "--cutoff". From the set of SNPs below the threshold, it then selects the most significant SNP in each genomic window (window size specified by --window-kb) that are correlated (correlation threshold set in --window-prune) and discards all other SNPs in that block.
            #--window-prune <float> - to specify the correlation squared threshold.
            #--window-cm <float>, --window-kb <float> or --window-length <integer> - to specify the window size (how far to search for correlated predictors, where the units are centiMorgans, kilobase or number of predictors, respectively). Note that --window-length ALL will tell LDAK to consider all predictors on the same chromosome.
            #--pvalues <pvalues> - to provide p-values for each predictor (the file <pvalues> should have two columns, that provide predictor names then p-values). When LDAK finds two highly correlated predictors, it will discard the one with the highest p-value.
            #--cutoff <float> - to provide the p-value threshold (LDAK will only consider predictors with p-values below this threshold). In other words, predictors with p-values below the cutoff will be treated as top predictors.
                #Note that this is thresholding, just remove SNPs above a threshold, and then clumping will be applied on the resulting list of SNPs.
                #See for example this examples from LDAK page:
                    #./ldak.out --thin-tops clump --bfile human --window-prune 0.05 --window-cm 1 --pvalues quant.pvalues --cutoff 5e-8
                    #In total, there were six predictors with p-value less than 5e-8; after clumping, four of these remain (listed in the file clump.in).
                    #Therefore, you got first 6 predictors after thresholding, and then after clumping you lose 2 more, so --cutoff is doing thresholding while --window-prune/--window-kb are doing clumping.
            #output:
                #The lists of retained and discarded predictors are saved in the files <output>.in and <output>.out, respectively.
            #https://dougspeed.com/clumping/

        print_text("run again linear on the reduced set of SNPs", header=4)
        run_bash(" \
            ldak6.1.linux \
                --linear ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_linear_clump_thresholding_" + str(threshold) + " \
                --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
                --pheno ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_response.tsv \
                --covar ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_cont.tsv \
                --factors ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_covars_factors.tsv \
                --extract ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_small_set_predictors_set_linear_clump_thresholding_" + str(threshold) + "_predictors.in \
                --permute NO \
        ")
            #You can use --keep <keepfile> and/or --remove <removefile> to restrict to a subset of samples, and --extract <extractfile> and/or --exclude <excludefile> to restrict to a subset of predictors (for more details, see Data Filtering).
                #I am using as input just a file with one column, the SNP names. There is no problem of having similar IDs across chromosomes because our IDs include the chromosome name, position and alelles.
                #https://dougspeed.com/data-filtering/

        print_text("do the calculation of the PRS in the same set of samples", header=3)
        print_text("first using the phenotype as input", header=3)
        run_bash(" \
            ldak6.1.linux \
                --calc-scores ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_linear_clump_thresholding_" + str(threshold) + "_prs_calc_with_pheno \
                --scorefile ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".score \
                --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
                --pheno ./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_transform_subset_response.tsv \
                --power 0 \
        ")
        print_text("then without phenotype, so we can check the PRS calculation of input reposne and covariates and hence we do not need to add covariates", header=3)
        run_bash(" \
            ldak6.1.linux \
                --calc-scores ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_linear_clump_thresholding_" + str(threshold) + "_prs_calc_without_pheno \
                --scorefile ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".score \
                --bfile ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
                --power 0 \
        ")

        print_text("check the PRS is the same", header=4)
        run_bash(" \
            cd ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/; \
            n_differ=$( \
                awk \
                    'BEGIN{FS=\"\t\"}{ \
                        if(NR==FNR){ \
                            a[$5]; \
                            next \
                        } \
                        if(!($5 in a)){ \
                            print $0 \
                        } \
                    }' \
                    " + response_variable + "_linear_clump_thresholding_" + str(threshold) + "_prs_calc_with_pheno.profile \
                    " + response_variable + "_linear_clump_thresholding_" + str(threshold) + "_prs_calc_without_pheno.profile | \
                wc -l \
            ); \
            if [[ $n_differ -ne 0 ]]; then \
                echo \"ERROR: FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLES FOR LINEAR\"; \
            fi; \
        ")
            #do operations in the first file, i.e., when NR=FNR
                #NR: The total number of records (lines) processed so far across all input files.
                #FNR: The number of records (lines) processed in the current file.
                #when they are the same it means we are processing the first file, and when they are different it means we are processing the second file.
            #the operation to do is:
                #a[$5]: The fifth field (PRS for not thresholding) of the current line is used as an index in the array a. This means that the script is storing the values of the fifth field from the first file in the array a.
                #the "next" statement is used to skip the remaining commands in the current iteration of the loop and move to the next line of the input file. This makes sense because if NR=FNR, then we are processing the first file and we do not need to run the code dedicated to the second file.
            #the code after the "next" statement is executed only when we are processing the second file (i.e., when NR!=FNR). In this case, it checks if the value of the fifth field ($5) is not present in the array a. If it is not present, it prints the entire line.
            #numbre of problematic lines is stored in the variable n_differ, and if it is not equal to 0, then we have a problem.


    print_text("check we have used the correct samples in both analyses", header=2)
    print_text("load the FAM files used for training and test", header=3)
    fam_file = pd.read_csv( \
        "./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing.fam", \
        sep=" ", \
        header=None \
    )

    print_text("samples in files generated by elastic net", header=3)
    combined_elastic = pd.read_csv( \
        "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_elastic.combined", \
        sep=" ", \
        header=0 \
    )
    if ( \
        (not fam_file[1].rename("IID").equals(combined_elastic["IID"])) \
    ):
        raise ValueError("ERROR: FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLES FOR ELASTIC NET")

    print_text("samples in files generated by linear", header=3)
    #threshold=list_thresholds[0]
    for threshold in list_thresholds:
        combined_linear = pd.read_csv( \
            "./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".combined", \
            sep=" ", \
            header=0 \
        )
        if ( \
            (not fam_file[1].rename("IID").equals(combined_linear["IID"])) \
        ):
            raise ValueError("ERROR: FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLES FOR ELASTIC NET")


    print_text("compress the results", header=2)
    print_text("remove bed and bim plink files initialles used as input (compressed files already created)", header=3)
    run_bash(" \
        cd ./data/plink_filesets/small_set_predictors/" + response_variable + "_filesets/; \
        rm \
            ./" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing.bed \
            ./" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing.bim; \
    ")
    
    print_text("elastic outputs", header=3)
    #not needed

    print_text("linear outputs", header=3)
    print_text("first raw outputs before clumping", header=4)
    run_bash(" \
        cd ./results/final_results/analysis_full_data/" + response_variable + "/; \
        gzip \
            --force \
                " + response_variable + "_small_set_predictors_set_linear_raw.assoc \
                " + response_variable + "_small_set_predictors_set_linear_raw.pvalues \
                " + response_variable + "_small_set_predictors_set_linear_raw.score \
                " + response_variable + "_small_set_predictors_set_linear_raw.summaries; \
    ")

    print_text("then outputs after clumping", header=4)
    #threshold=list_thresholds[0]
    for threshold in list_thresholds:
        run_bash(" \
            cd ./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(threshold) + "/; \
            gzip \
                --force \
                ./" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".assoc \
                ./" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".pvalues \
                ./" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".score \
                ./" + response_variable + "_linear_clump_thresholding_" + str(threshold) + ".summaries; \
        ")


    print_text("plot PRS against phenotype", header=2)
    print_text("load the phenotype data", header=3)
    original_cleaned_data = pd.read_csv( \
        "./data/pheno_data/pheno_data_cleaned.tsv", \
        sep="\t", \
        header=0 \
    )


    original_cleaned_data_response = original_cleaned_data.loc[original_cleaned_data["ID"].isin(pheno_subset_transform["family_id"] + ":" + pheno_subset_transform["AGRF code"]), ["ID", response_variable]]

    original_cleaned_data_response[["ID1", "ID2"]] = original_cleaned_data_response["ID"].str.split(":", expand=True)

    # Create a figure with 8 subplots arranged in a grid (e.g., 2 rows x 4 columns)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))  # Adjust figsize as needed

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    #model_type=list_thresholds[0]
    for model_type in ["elastic"]+list_thresholds:

        selected_ax = axes[list_thresholds.index(model_type)]

        if model_type == "elastic":
            prs_file = pd.read_csv( \
                "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_elastic_prs_calc_without_pheno.profile", \
                sep="\t", \
                header=0 \
            )[["ID1", "ID2", "Profile_1"]]
        else:
            prs_file = pd.read_csv( \
                "./results/final_results/analysis_full_data/" + response_variable + "/clump_thresholding_" + str(model_type) + "/" + response_variable + "_linear_clump_thresholding_" + str(model_type) + "_prs_calc_without_pheno.profile", \
                sep="\t", \
                header=0 \
            )[["ID1", "ID2", "Profile_1"]]

        pheno_prs = prs_file.merge(original_cleaned_data_response, on=["ID1", "ID2"])


        pheno_prs["ID"].equals(pheno_prs["ID1"] + ":" + pheno_prs["ID2"])



        # Calculate 20 quantiles based on "response_variable"
        pheno_prs["quantile"] = pd.qcut(pheno_prs[response_variable], q=20, labels=False, duplicates="drop")

        # Initialize a list to store results
        results = []

        # Iterate over each quantile
        #quantile=pheno_linear_prs["quantile"].unique()[0]
        for quantile in pheno_prs["quantile"].unique():
            # Select samples within the current quantile
            quantile_data = pheno_prs[pheno_prs["quantile"] == quantile]
            
            # Calculate mean and standard deviation for the "prs" column
            lower_ci_pheno_linear = quantile_data[response_variable].quantile(0.025)
            median_pheno_linear = quantile_data[response_variable].quantile(0.5)
            higher_ci_pheno_linear = quantile_data[response_variable].quantile(0.975)


            # Append the results
            results.append({"quantile": quantile, "lower_ci_prs": lower_ci_pheno_linear, "median_prs": median_pheno_linear, "higher_ci_prs": higher_ci_pheno_linear})

        # Convert results to a DataFrame
        quantile_stats = pd.DataFrame(results)

        # Print the results
        print(quantile_stats)


        import matplotlib.pyplot as plt


        

        # Plot the mean as dots
        selected_ax.scatter(quantile_stats["quantile"], quantile_stats["median_prs"], color="black", label="Median PRS")

        # Plot the error bars (mean ± sd)
        selected_ax.errorbar(
            quantile_stats["quantile"], 
            quantile_stats["median_prs"], 
            yerr=quantile_stats["higher_ci_prs"]-quantile_stats["lower_ci_prs"],
            fmt="o", 
            color="black", 
            capsize=5, 
            label="Median ± 95CI"
        )

        # Add labels and title
        selected_ax.set_xlabel("Quantile")
        selected_ax.set_ylabel([k for k, v in dict_change_responses.items() if v==response_variable][0])
        selected_ax.set_title("Thresholding (P<" + str(model_type) + ") + Clumping")



    # Adjust layout to prevent overlapping
    plt.tight_layout()

    # Show the plot
    plt.savefig("./results/final_results/analysis_full_data/" + response_variable + "/prs_quantiles.png", dpi=300)
    plt.close()





#parallelize per phenotype


# endregion




#manhattan plot showing low number of snps betqeen nominal and bonferroni
    #low significance
#correlation between PRS and phenotype in full dataset
    #high correlation for elastic and linear P=1, problem with stringent thresholds, low signifiance
#CV across 100 iterations, lack of correlation specially for stringent thresholds, this and hers estimate suggests we do not have enough power
    #compare hers with bouchard...
    #maybe hers noisy because of low sample size but still thresholds results say teh same, underpower


    #For instance, if the GWAS data are relatively underpowered, then the optimal threshold is more likely to be P = 1 (all SNPs) even if a small fraction of SNPs are causal (see ref. 5 for details).
        #OUR CASE, argument in favour of limited power in our case, as P=1 is the best threshold. we get very few significant SNPs, and the best threshold is 1, i.e., just clumping.
        #hers is very high also, but cor is low, same line but not sure we can trust hers estimate withe this very low sample size
        #https://www.nature.com/articles/s41596-020-0353-1


"""
Notes about Genomic control and confounding effects from Dr. Speed's slides

- Confounding due to population structure in GWAS
    Compared with observational epidemiology, GWAS have few opportunities for confounding bias. The main problem is population structure, which refers to 
        - mating patterns within a pop -> subpops (more relatedness within than between)
        - allele frequency differences across subpops 
        - environmental exposures may also vary across subpops. 

    Phenotypes can also vary across subpops, because 
        - the causal alleles vary in frequency and/or
        - they vary with environmental factors correlated with pop structure, and/or 
        - ascertainment bias: recruitment of phenotypic groups differs across subpops 

    These can lead to significant associations not with CVs but with SNPs whose allele frequencies correlate with trait across subpops

- Adjustments for pop structure confounding: genomic control (GC) + PCA
    GC was an early approach to adjusting for confounding, based on the idea that pop structure can lead to many significant SNPs genome-wide 
        - all association test stats (with χ21 null distribution) are divided by the ratio of empirical to null medians (called a genomic inflation factor, GIF) provided GIF > 1 
        - assumes sparsity: true causals are rare, so there are few non-null test stats, so median test statistic is close to null value. 
    
    However, the omnigenic nature of many complex traits means that the assumption is false and GC is overly conservative. 
    
    The first few eigenvectors (or principal components) of XX T often reflect pop structure   
        - Included as covariates in GWAS regression models, they can absorb pop structure effects on the trait. 
    
    Now Mixed Model Association Analysis (MMAA) is the preferred approach to adjusting for pop structure effects in tests of association (see slide 40 from Dr. Speed's slides; teaching_slides.pdf)

- Given we have performed a fine analysis of population structure, removing many PCA outliers (we only have 1 ancestry) and then calculating PCAs on very clean data and considering these PCAs in the models, we have covered this point, so no need for genomic control.
"""

"""
Quantile plot of transformed phenotype against PRS quantiles

Quantile plots corresponding to the effect of a PRS on a normally distributed target trait should reflect the S-shape of the probit function (Fig. 5a). This is because the trait values are more spread out between quantiles at the tails of a normal distribution. Thus, plotting quantiles of PRS versus (absolute) effect on trait shows increasingly larger jumps up/down the y-axis from the median to the extreme upper/lower quantiles. When unequal strata are plotted, with the smallest strata at the tails, then this effect appears stronger. When the target outcome is disease status and prevalence or OR are plotted on the y-axis, then the shape is expected to be different: here, the shape is asymmetrical, showing a marked inflection at the upper end (Fig. 5b), since cases are enriched at the upper end only. Thus, inflections of risk at the tails of the PRS distribution82,83 should be interpreted according to these statistical expectations and not as interesting in themselves.

https://www.nature.com/articles/s41596-020-0353-1
"""


def manhattan_plot(covariate_dataset, response_variable):
        

    dict_change_names={
        "family_id": "FID",
        "AGRF code": "IID",
        "Age": "age",
        "Week 1 Body Mass": "week_1_weight",
        "Week 1 Beep Test": "week_1_beep",
        "Week 1 Distance (m)": "week_1_distance",
        "Week 1 Pred VO2max": "week_1_vo2"
    }

    print_text("specify the covariates", header=3)
    selected_covariates = pheno_subset.columns[~pheno_subset.columns.isin(["family_id", "AGRF code", response_variable])]
    print(selected_covariates)
    #check we have correct covariates
    total_list_covariates = ["Age", "sex_code", "Week 1 Body Mass", "Week 1 Beep Test", "Week 1 Distance (m)", "Week 1 Pred VO2max"] + [f"PCA{i}" for i in range(1,21)]
    if(sum([1 for cov in selected_covariates if cov not in total_list_covariates])!=0):
        raise ValueError("ERROR: FALSE! WE HAVE COVARIATES THAT ARE NOT IN THE TOTAL LIST")


    pheno_subset_transform = pd.read_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )

    #save the response variable after changing the names
    pheno_subset_transform[["family_id", "AGRF code", response_variable]].rename(columns=dict_change_names).to_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_response.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )

    #save the covariates that are factors
    if ("sex_code" in selected_covariates):

        #save sex_code
        pheno_subset_transform[["family_id", "AGRF code", "sex_code"]].rename(columns=dict_change_names).to_csv( \
            "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_factors.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )

    #save the covariates that are continuous
    selected_covariates_cont = [cov for cov in selected_covariates if cov != "sex_code"]
    pheno_subset_transform[["family_id", "AGRF code"] + selected_covariates_cont].rename(columns=dict_change_names).to_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_cont.tsv", \
        sep="\t", \
        header=True, \
        index=False, \
        na_rep="NA" \
    )

    #create a file with the samples of the selected set
    pheno_subset_transform[["family_id", "AGRF code"]].to_csv( \
        "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_samples_in.tsv", \
        sep="\t", \
        header=False, \
        index=False, \
        na_rep="NA" \
    )

    #use that file to select the corresponding sample from the plink fileset
    run_bash(" \
        plink \
            --bfile ./data/plink_filesets/" + covariate_dataset + "/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
            --keep ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_samples_in.tsv \
            --make-bed \
            --out ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_plink_fileset \
    ")
        #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.


    run_bash(" \
        ldak6.1.linux \
            --linear ./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/" + response_variable + "_full_dataset_linear \
            --bfile ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_plink_fileset \
            --pheno ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_response.tsv \
            --covar ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_cont.tsv \
            --factors ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform_covars_factors.tsv \
            --permute NO \
    ")

    #--linear
        #https://dougspeed.com/single-predictor-analysis/


    #clumping?

    print("load assoc results to pandas")
    assoc_results = pd.read_csv( \
        "./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/" + response_variable + "_full_dataset_linear.assoc", \
        sep="\t", \
        header=0, \
        low_memory=False)
    print(assoc_results)
        #CHECK COLUMNS
        #compare wald P with the .pvalue file
    

    print("check we have the correct dtypes for dash_bio.ManhattanPlot (see code below)")
    print(assoc_results["Chromosome"].dtype == "int64")
    print(assoc_results["Basepair"].dtype == "int64")
    print(assoc_results["Wald_P"].dtype == "float64")
    print(assoc_results["Predictor"].dtype == "O")

    dict_titles = { \
        "distance_change": "Change in distance (m)", \
        "beep_test": "Change in beep test", \
        "vo2max": "Change in predicted VO2max", \
        "weight": "Change in body mass", \
    }

    # import libraries
    from scipy.stats import uniform
    from scipy.stats import randint
    import matplotlib.pyplot as plt


    # -log_10(pvalue)
    assoc_results['minuslog10pvalue'] = -np.log10(assoc_results["Wald_P"])
    assoc_results["Chromosome"] = assoc_results["Chromosome"].astype('category')
    assoc_results = assoc_results.sort_values('Chromosome')

    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    assoc_results['ind'] = range(len(assoc_results))
    df_grouped = assoc_results.groupby(('Chromosome'))

    # manhattan plot
    fig = plt.figure(figsize=(22, 8)) # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    genome_wide_threshold = -np.log10(0.05/assoc_results.shape[0])
    ax.axhline(y=genome_wide_threshold, color='red', linestyle='--', linewidth=1.5, label=f'Genome-wide threshold ({genome_wide_threshold:.2f})')
    
    # set axis limits
    ax.set_xlim([0, len(assoc_results)])
    if (assoc_results['minuslog10pvalue'].max() > genome_wide_threshold):
        ax.set_ylim([0, assoc_results['minuslog10pvalue'].max() + 2])
    else:
        ax.set_ylim([0, genome_wide_threshold + 2])

    # x axis label
    ax.set_xlabel('Chromosome')

    # Save the plot as a static image
    plt.savefig("./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/" + response_variable + ".png", dpi=300)

    #code from:
        #https://python-graph-gallery.com/manhattan-plot-with-matplotlib/



    #check manhatan plots
        #there is a strange gap in VO2 max for one of the first chromosomes in the prelim results



