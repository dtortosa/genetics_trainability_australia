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

#This script will process the results of the association studies, summarizing the predictive power of the PRS across CV iterations.

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






###########################
# region Summary strategy #
###########################

#We are going to create the following table
    #40 rows:
        #2 dataset types x 4 phenotype for elastic
        #2 dataset types x 4 phenotype x 4 thresholds for linear
    #12 columns
        #3 columns for 95CI x 4 evaluation metrics
    #First will a go distance_change with small set and linear, then the same phentoype and dataset with elastic net, then the same phenotype with large and linear/elastic and so on. In this way, we can easily compare the classic PRS approach with the newer one for each phentoype and also small vs large predictors set. The next phenotypes will be beep_test and VO2 max, finishing with weight change.

# endregion






#################################
# region Check the output files #
#################################
print_text("Check the output files", header=1)
print_text("calculate all combinations of predictor set type and phenotypes", header=2)
print_text("define list with the dataset types and phenotypes", header=3)
dataset_types = ["small_set_predictors", "large_set_predictors"]
response_variables = ["distance_change", "beep_change", "vo2_change", "weight_change"]

print_text("Calculate all combinations", header=3)
combinations_pheno_dataset = list(product(response_variables, dataset_types))
    #Product is a function from the itertools module that computes the Cartesian product of the two lists.
    #It generates all possible pairs of elements, where the first element comes from dataset_type and the second element comes from response_variables
    #"product" generates a generator, which means it produces items one at a time as you iterate over it, rather than storing all the items in memory at once. Then, you have to call list() to do the operation to be done and save the results in a list, i.e., a list of tuples with each pair of dataset type and phenotype.
if(len(combinations_pheno_dataset)!=len(dataset_types)*len(response_variables)):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE COMBINATIONS")

print_text("iterate across combinations", header=3)
#response_variable="distance_change"; dataset_type="small_set_predictors"
for response_variable, dataset_type in combinations_pheno_dataset:
    
    print_text("dataset type: " + dataset_type + ", phenotype: " + response_variable, header=4)

    print("Define the file path")
    file_path="./03b_prs_calculation_100_iter_" + dataset_type + "_" + response_variable + ".out"

    print("Open and read the file into a list of lines")
    with open(file_path, "r") as file:
        lines = file.readlines()

    print("Check for lines with a problem")
    #line=lines[0]
    #line=lines[100]
    for line in lines:

        #Searches for the words in the line, ignoring case sensitivity.
        search_word = re.search(pattern=r"(error|false|fail)", string=line, flags=re.IGNORECASE)
            #pattern: A regular expression
                #The r in r"(error|false|fail)" indicates that the string is a raw string literal in Python. A raw string treats backslashes (\) as literal characters and does not interpret them as escape characters. In regular expressions, backslashes (\) are commonly used for special sequences (e.g., \b for word boundaries, \d for digits). Without the r prefix, Python would interpret the backslashes as escape characters, which could lead to errors or unintended behavior.
                #Matches any of the words but not as whole words
                #The \b in the pattern (at the beginning and at the end) ensures that "error" is matched as a whole word (not part of another word like "errors").
                #We do not want that. We explicitly explicitly want to match substrings of other words (e.g., "errors", "failing", "falsely"), so avoid b
                #The re.IGNORECASE flag makes the search case-insensitive, meaning it will match "Error", "ERROR", "error", etc.
        
        # Check for "error", if no match, then "search_word" will be None
        if search_word is not None:
            
            #Converts the line to lowercase and checks if any of the excluded phrases ("mean absolute error" or "mean squared error") are present in the line.
            #phrase=exclude_phrases[0]
            if not any(phrase in line.lower() for phrase in ["mean absolute error", "mean squared error"]):

                #if the error case does not include any of these phrases, then we have an actual error
                raise ValueError(f"ERROR! FALSE! WE FOUND AN ERROR FOR {dataset_type} and {response_variable}: {line.strip()}")

# endregion






#############################
# region generate the table #
#############################
print_text("Generate the table", header=1)
print_text("calculate all combinations of phenotypes, predictor dataset type and models", header=2)
print_text("define list with model types (see above for the rest of lists)", header=3)
response_variables = ["distance_change", "beep_change", "vo2_change", "weight_change"]
dataset_types = ["small_set_predictors", "large_set_predictors"]

print_text("Calculate all combinations by model type to avoid problems with thresholds", header=3)
combinations_pheno_dataset_iter = list(product(response_variables, dataset_types))
if(len(combinations_pheno_dataset_iter)!=len(response_variables)*len(dataset_types)):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE COMBINATIONS")

print_text("define function to extract the evaluation metrics", header=2)
#response_variable="beep_change"; dataset_type="large_set_predictors"; model_type="elastic"; iteration=1
#response_variable="beep_change"; dataset_type="large_set_predictors"; model_type="linear"; iteration=1
def extract_metrics(response_variable, dataset_type, model_type, iteration):

    if model_type=="linear":

        threshold_folders = run_bash(" \
            ls ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/test_set/linear/ \
        ", return_value=True).split("\n")[0:-1]

        threshold_list = sorted([float(i.split("_")[-1]) for i in threshold_folders], reverse=True)




        evaluation_metrics_raw_list = list()

        #threshold=threshold_list[0]
        for threshold in threshold_list:

            evaluation_metrics_raw = run_bash(" \
                awk \
                    'BEGIN{FS=\" \"}{ \
                        if ($1==1) { \
                            print $3; \
                        } \
                    }' \
                    ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/test_set/" + model_type + "/clump_thresholding_" + str(threshold) + "/" + response_variable + "_prs_calculation_classical_jacknife_eval.jack \
            ", return_value=True).split("\n")[:-1]

            evaluation_metrics_raw_list.append(tuple(evaluation_metrics_raw))


        evaluation_metrics = tuple([float(item) for sublist in evaluation_metrics_raw_list for item in sublist])

        evaluation_metrics

        column_names = [f"{metric}_{threshold}" for threshold in threshold_list for metric in ["correlation", "squared_correlation", "mean_squared_error", "mean_absolute_error"]]
        evaluation_metrics = [tuple(column_names), evaluation_metrics]
        
        
    elif model_type=="elastic":

        evaluation_metrics_raw = run_bash(" \
            awk \
                'BEGIN{FS=\" \"}{ \
                    if ($1==1) { \
                        print $3; \
                    } \
                }' \
                ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/test_set/" + model_type + "/" + response_variable + "_prs_jacknife_eval.jack \
        ", return_value=True).split("\n")[:-1]

        evaluation_metrics = tuple([float(metric) for metric in evaluation_metrics_raw])

        column_names = ["correlation", "squared_correlation", "mean_squared_error", "mean_absolute_error"]
        evaluation_metrics = [tuple(column_names), evaluation_metrics]


    return evaluation_metrics


# Initialize an empty DataFrame for the reshaped data
eval_metrics = pd.DataFrame()

#response_variable="beep_change"; dataset_type="large_set_predictors"
for response_variable, dataset_type in combinations_pheno_dataset_iter:

    ##get elastic
    elastic_metrics_list = list()
    for iteration in range(1, 2):
        elastic_metrics_list.append(extract_metrics(response_variable=response_variable, dataset_type=dataset_type, model_type="elastic", iteration=iteration))

    # Extract column names from the first tuple of the first sublist
    column_names = elastic_metrics_list[0][0]

    # Extract all float tuples (second element of each sublist)
    rows = [sublist[1] for sublist in elastic_metrics_list]

    # Create the DataFrame
    elastic_metrics_df = pd.DataFrame(rows, columns=column_names)

    # Extract all float tuples (second element of each sublist)
    rows = [sublist[1] for sublist in elastic_metrics_list]

    # Create the DataFrame
    elastic_metrics_df = pd.DataFrame(rows, columns=column_names)


    # Calculate percentiles for each column
    percentiles = [0.025, 0.5, 0.975]
    #p=0.025
    result_elastic = pd.DataFrame(
        {f"{p}th_percentile": elastic_metrics_df.quantile(p) for p in percentiles}
    )

    result_elastic["metric"] = result_elastic.index

    # Flatten the rows into a single row
    row = {}
    #row_data=result_elastic.iloc[0,:]
    for _, row_data in result_elastic.iterrows():
        metric = row_data["metric"]
        row["phenotype"] = response_variable
        row["dataset_type"] = dataset_type
        row["model_type"] = "elastic"
        row["threshold"] = np.nan
        row[f"{metric}_0.025th_percentile"] = row_data["0.025th_percentile"]
        row[f"{metric}_0.5th_percentile"] = row_data["0.5th_percentile"]
        row[f"{metric}_0.975th_percentile"] = row_data["0.975th_percentile"]
    
    eval_metrics = pd.concat([eval_metrics, pd.DataFrame([row])])


    ##get linear
    linear_metrics_list = list()
    for iteration in range(1, 2):
        linear_metrics_list.append(extract_metrics(response_variable=response_variable, dataset_type=dataset_type, model_type="linear", iteration=iteration))

    # Extract column names from the first tuple of the first sublist
    column_names = linear_metrics_list[0][0]

    #CHECK WE AHVE RTHE SAME COLUMN NAMES, I.E., THE SAME THRESHOLDS IN ALL CASES

    # Extract all float tuples (second element of each sublist)
    rows = [sublist[1] for sublist in linear_metrics_list]

    # Create the DataFrame
    linear_metrics_df = pd.DataFrame(rows, columns=column_names)


    # Calculate percentiles for each column
    percentiles = [0.025, 0.5, 0.975]
    #p=0.025
    result = pd.DataFrame(
        {f"{p}th_percentile": linear_metrics_df.quantile(p) for p in percentiles}
    )


    # Extract thresholds and metrics from the index
    result["threshold"] = result.index.str.split("_").str[-1]
    result["metric"] = ["_".join(i[0:-1]) for i in result.index.str.split("_")]

    


    # Iterate over unique thresholds
    #threshold=result["threshold"].unique()[0]
    for threshold in result["threshold"].unique():
        
        # Filter rows for the current threshold
        subset = result[result["threshold"] == threshold]
        
        # Flatten the rows into a single row
        row = {}
        #row_data=subset.iloc[0,:]
        for _, row_data in subset.iterrows():
            metric = row_data["metric"]
            row["phenotype"] = response_variable
            row["dataset_type"] = dataset_type
            row["model_type"] = "linear"
            row["threshold"] = threshold
            row[f"{metric}_0.025th_percentile"] = row_data["0.025th_percentile"]
            row[f"{metric}_0.5th_percentile"] = row_data["0.5th_percentile"]
            row[f"{metric}_0.975th_percentile"] = row_data["0.975th_percentile"]
        
        
        # Append the row to the reshaped DataFrame
        eval_metrics = pd.concat([eval_metrics, pd.DataFrame([row])])




# Print the reshaped DataFrame
print(eval_metrics.iloc[:, 0:6])

eval_metrics.to_csv("./results/final_results/evaluation_metrics.tsv", sep="\t", header=True, index=False, na_rep="NA")


#in the results of --linear
    #I think you only have look at the first row of cors because we are already applying a threshold with --thin-tops, if you used a threshold during clumping of 0.0001, no snps above 0.0001 are going to be present, thus it does not make sense to see the different profiles generated by --calc-score thorugh 0.01, 0.001... we are already iterating over that list of thresholds, see previous code ofr further details

#we need to check inflation? maybe just bonferroni and show no snp is signifncat, so PRSs are useful

#check overlap across iterations for the samples with the highest PRS values in the training or test set?


##WHEN INTERPRETING THE RESULTS OF THE PRS, LOOK SLIDES FROM DOUG
    #https://dougspeed.com/short/
###WHEN DONE, YOU CAN CHECK THE SECTION OF INTERPRETATION FROM ORRELLY
    #Interpretation and presentation of results


