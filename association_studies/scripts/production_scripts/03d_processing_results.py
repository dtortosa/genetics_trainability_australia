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
######## PROCESS ASSOCIATION RESULTS ########
#############################################

#This script will process the results of the association studies, summarizing the predictive power of the PRS across CV iterations and generating Manhattan plots.

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






###########################
# region Summary strategy #
###########################



# endregion


#in the results of --linear
    #I think you only have look at the first row of cors because we are already applying a threshold with --thin-tops, if you used a threshold during clumping of 0.0001, no snps above 0.0001 are going to be present, thus it does not make sense to see the different profiles generated by --calc-score thorugh 0.01, 0.001... we are already iterating over that list of thresholds, see previous code ofr further details

#we need to check inflation? maybe just bonferroni and show no snp is signifncat, so PRSs are useful

#check overlap across iterations for the samples with the highest PRS values in the training or test set?


##WHEN INTERPRETING THE RESULTS OF THE PRS, LOOK SLIDES FROM DOUG
    #https://dougspeed.com/short/
###WHEN DONE, YOU CAN CHECK THE SECTION OF INTERPRETATION FROM ORRELLY
    #Interpretation and presentation of results


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


    ##INFLATION FACTOR?

    #check manhatan plots
        #there is a strange gap in VO2 max for one of the first chromosomes in the prelim results



#FOR BAT ANALYSES
    #this would an additional step in this project that would be outside of the paper
    #take the 1000kb gene windows for all coding genes, liftover to hg38. If the USCS tool accepts genomic ranges, just use them as input, if not, split in two datasets the start and the end of the gene windows
    #for each phenotype (VO2, beep....), calculate the average (better than median because want influence of outliers within gene like in iHS, if a SNPs is veery important in a gene that should influence the info about the whole gene) p-value for the association of SNPs inside each gene
    #then, calculate 1000 random sets of genes, within each set, calculate the median association of all genes inside the set and compare with the BAT set to obtain a distribution and empirical p-value (is association lower in BAT? LOOF BAT PAPER SCRIPTS FOR THIS). Here I want median because i do not want a gene outliser change things, I want the overall impact of BAT genes in general, not just a few genes.
