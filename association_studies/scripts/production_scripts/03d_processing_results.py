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
    #16 columns
        #3 columns for 95CI x 4 evaluation metrics + 4 columns indicating phenotype, dataset type, model type and threshold
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

    print("check the last line indicates the end of the process")
    if(lines[-3]!="FINISH\n"):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE LAST LINE OF THE FILE. IT SHOULD BE 'FINISH', but it is not")

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


print_text("summarize the evaluation metrics and store them in a table", header=2)
print_text("define function to extract the values of the metrics", header=3)
#response_variable="beep_change"; dataset_type="large_set_predictors"; model_type="elastic"; iteration=1
#response_variable="beep_change"; dataset_type="large_set_predictors"; model_type="linear"; iteration=1
def extract_metrics(response_variable, dataset_type, model_type, iteration):

    #if the model type is linear, we need to get the thresholds from the folders
    if model_type=="linear":

        #get the names of folders with thresholds results
        threshold_folders = run_bash(" \
            ls ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/test_set/linear/ \
        ", return_value=True).split("\n")[0:-1]
            #split the output by new line and avoid the last element which is an empty space

        #get the value of each threshold
        #threshold=threshold_folders[0]
        threshold_list = sorted([int(threshold.split("_")[-1]) if threshold=="clump_thresholding_1" else float(threshold.split("_")[-1]) for threshold in threshold_folders], reverse=True)
            #for each folder name, split the name by "_" and get the last element, which is the threshold value, then convert the string to int if the threshols is 1 (we do not want to have 1.0, which is not the actual name of the folder with results but 1), or to float if the threshodl is not 1, then sort the list in descending order

        #iterate across thresholds to get the evaluation metrics
        #empty list
        evaluation_metrics_raw_list = list()
        #threshold=threshold_list[0]
        for threshold in threshold_list:

            #get the evaluation metrics for each threshold
            evaluation_metrics_raw = run_bash(" \
                awk \
                    'BEGIN{FS=\" \"}{ \
                        if ($1==1) { \
                            print $3; \
                        } \
                    }' \
                    ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/test_set/linear/clump_thresholding_" + str(threshold) + "/" + response_variable + "_prs_calculation_classical_jacknife_eval.jack \
            ", return_value=True).split("\n")[:-1]
                #get only the rows with the first column equal to 1, which are the evaluation metrics for the first PRS. From these rows, get the third column, which is the estimate of the evaluation metric
                    #I understand that we have 7 scores corresponding with ALL (P<1), P<0.1, P<0.01, P<0.001, P<0.0001, P<0.00001, P<0.000001. So, given we have already filtered in previous steps with --thin-tops, we can use the first row considering all SNPs that were used as input (which in our case are already filtered)
                    #Doug has confirmed that the 7 scores you get in jackknife correspond to these 7 thresholds.
                #split the output by new line and avoid the last element which is an empty space
            
            #convert the list to a tuple and append to the list with the results across thresholds
            evaluation_metrics_raw_list.append(tuple(evaluation_metrics_raw))

        #convert all metrics to float and flatten the list to get all of them in a single tuple
        #sublist=evaluation_metrics_raw_list[0]; item=sublist[0]
        evaluation_metrics = tuple([float(item) for sublist in evaluation_metrics_raw_list for item in sublist])

        #get the names of the columns
        column_names = [f"{metric}_{threshold}" for threshold in threshold_list for metric in ["correlation", "squared_correlation", "mean_squared_error", "mean_absolute_error"]]
        
        #add the threshold to the column names
        evaluation_metrics = [tuple(column_names), evaluation_metrics]
        
    #else, if the model is elastic  
    elif model_type=="elastic":

        #get the heritability metrics from elastic
        heritability_elastic_raw = run_bash(" \
            awk \
                'BEGIN{FS=\" \"}{ \
                    if ($1==\"RHE\") { \
                        print $1, $4; \
                    } else if ($1==\"MCREML\") { \
                        print $1, $4; \
                    } \
                }' \
                ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/training_set/elastic/" + response_variable + "_training_elastic.hers \
        ", return_value=True) 

        #get a list with each heritability metric as a different element, i.e., split the string 
        heritability_elastic = re.split(r"[ \n]+", heritability_elastic_raw)[:-1]
            #split the resulting string by space and new line and avoid the last element which is an empty space

        #get the evaluation metrics for the elastic model
        evaluation_metrics_raw = run_bash(" \
            awk \
                'BEGIN{FS=\" \"}{ \
                    if ($1==1) { \
                        print $3; \
                    } \
                }' \
                ./results/final_results/" + dataset_type + "/" + response_variable + "/train_test_iter_" + str(iteration) + "/test_set/elastic/" + response_variable + "_prs_jacknife_eval.jack \
        ", return_value=True).split("\n")[:-1]

        #convert the list of evaluation metrics to a tuple of floats
        #metric=evaluation_metrics_raw[0]
        evaluation_metrics_raw_2 = tuple([float(metric) for metric in evaluation_metrics_raw])

        #make a single tuple with the evaluation metrics and the heritability metrics
        #i=heritability_elastic[0]
        #i=heritability_elastic[1]
        evaluation_metrics = evaluation_metrics_raw_2 + tuple(float(i) for i in heritability_elastic if i not in ["RHE", "MCREML"])
            #add to the evaluation metrics the heritability metrics, but only the values, not the names

        #get the names of the columns
        #i=heritability_elastic[0]
        #i=heritability_elastic[1]
        column_names = ["correlation", "squared_correlation", "mean_squared_error", "mean_absolute_error"] + [i for i in heritability_elastic if i in ["RHE", "MCREML"]]
            #we first have the evaluation metrics and then the heritability metrics, so we can just add the names of the heritability metrics to the end of the list

        #make a list with two tuples, the column names and the values
        evaluation_metrics = [tuple(column_names), evaluation_metrics]

    #return the evaluation metrics
    return evaluation_metrics

print_text("run the function across all combinations", header=3)
print_text("Initialize an empty DataFrame for the reshaped data", header=4)
eval_metrics = pd.DataFrame()

print_text("loop the function", header=4)
#response_variable="beep_change"; dataset_type="large_set_predictors"
#[i for i in combinations_pheno_dataset_iter if i[1]=="small_set_predictors"]
for response_variable, dataset_type in combinations_pheno_dataset_iter:

    print("get elastic metrics for " + dataset_type + " and " + response_variable)
    #empty list to save results
    elastic_metrics_list = list()

    #loop across iterations to extract the metrics and save in the empty list
    #iteration=1
    for iteration in range(1, 101):
        elastic_metrics_list.append(extract_metrics(response_variable=response_variable, dataset_type=dataset_type, model_type="elastic", iteration=iteration))

    #check we have the same column names in all iterations
    #i=elastic_metrics_list[0]
    if(sum([1 for i in elastic_metrics_list if i[0]!=elastic_metrics_list[0][0]])!=0):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE COLUMN NAMES FOR ELASTIC METRICS")
            #for each iteration, check whether the names (first element) is not the same than the names in the first iteration, if not, then add 1. 
            #if the sum is not 0, then we have a problem.

    #check we have the correct number of columns
    #i=elastic_metrics_list[0]
    if(sum([1 for i in elastic_metrics_list if (len(i[0])!=6) | (len(i[1])!=6)])!=0):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF COLUMNS FOR ELASTIC METRICS")
    
    #extract column names from the first tuple of the first sublist
    column_names = elastic_metrics_list[0][0]

    #extract all float tuples with the metrics (second element of each sublist)
    #sublist=elastic_metrics_list[0]
    elastic_tuples_metrics = [sublist[1] for sublist in elastic_metrics_list]

    #create the DataFrame
    elastic_metrics_df = pd.DataFrame(elastic_tuples_metrics, columns=column_names)

    #calculate percentiles for each column
    percentiles = [2.5, 50, 97.5]
    #p=2.5
    result_elastic = pd.DataFrame(
        {f"{p}th_percentile": elastic_metrics_df.quantile(p/100) for p in percentiles}
    )
        #with a comprehension, calculate the three percentiles for each column using "elastic_metrics_df.quantile()". 
        #each percentile is store as a key/value pair in the dictionary, where the key is the name of the column with the percentile and the value is the result of the quantile calculation for that column.
        #the dict is converted to a DF where each key (percentile) is a column

    #flatten the rows into a single row, i.e., all metrics of a given phenotype and dataset type are in a single row
    elastic_row = {}
        #this empty dict will store the values of each metric and percentile as a different key
    
    #add the phenotype, model.... to the empty dict
    elastic_row["phenotype"] = response_variable
    elastic_row["dataset_type"] = dataset_type
    elastic_row["model_type"] = "elastic"
    elastic_row["threshold"] = np.nan
    
    #iterate over the rows of the data.frame
    #row_index=result_elastic.index[0]; row_data=result_elastic.iloc[0,:]
    for row_index, row_data in result_elastic.iterrows():

        #get the metric name
        metric = row_index

        #add the percentiles
        elastic_row[f"{metric}_2.5th_percentile"] = row_data["2.5th_percentile"]
        elastic_row[f"{metric}_50th_percentile"] = row_data["50th_percentile"]
        elastic_row[f"{metric}_97.5th_percentile"] = row_data["97.5th_percentile"]

    #convert the dict to a DF (each key as a column) and the concatenate to the eval metrics
    eval_metrics = pd.concat([eval_metrics, pd.DataFrame([elastic_row])])

    print("get linear metrics for " + dataset_type + " and " + response_variable)
    #empty list to save results
    linear_metrics_list = list()

    #loop across iterations to extract the metrics and save in the empty list
    #iteration=1
    for iteration in range(1, 101):
        linear_metrics_list.append(extract_metrics(response_variable=response_variable, dataset_type=dataset_type, model_type="linear", iteration=iteration))

    #check we have the same column names in all iterations
    #we cannot do this for linear because we can have different thresholds for different iterations if there smaller p-values in some iterations that are below more stringent thresholds
    #i=linear_metrics_list[0]
    #if(sum([1 for i in linear_metrics_list if i[0]!=linear_metrics_list[0][0]])!=0):
        #raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE COLUMN NAMES FOR LINEAR METRICS")
            #for each iteration, check whether the names (first element) is not the same than the names in the first iteration, if not, then add 1. 
            #if the sum is not 0, then we have a problem.

    #check we have the correct number of columns
    #we cannot do this check for the same reason as above
    #i=linear_metrics_list[0]
    #if(sum([1 for i in linear_metrics_list if (len(i[0])!=6) | (len(i[1])!=6)])!=0):
        #raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF COLUMNS FOR LINEAR METRICS")

    #get the max number of metric names we have
    #sublist=linear_metrics_list[0]
    max_n_names = max([len(sublist[0]) for sublist in linear_metrics_list])
        #for each sublist, which is a list of two tuples, the names of the metrics and the actual values, calculate the length of the first element, i.e., the names of the metrics
        #get the maximum length of the names of the metrics across all iterations
        #we are doing this because we can have different thresholds for different iterations if there smaller p-values in some iterations that are below more stringent thresholds and hence we would have metrics for more thresholds in one iterations than in others

    #get the column names for the cases with the max number of thresholds
    #sublist=linear_metrics_list[0]
    column_names = [sublist[0] for sublist in linear_metrics_list if len(sublist[0])==max_n_names]
        #select the sublist with the maximum number of names and get the names of the metrics

    #extract all float tuples (second element of each sublist)
    #sublist=linear_metrics_list[0]
    linear_tuples_metrics = [sublist[1] for sublist in linear_metrics_list]

    #create the DataFrame
    linear_metrics_df = pd.DataFrame(linear_tuples_metrics)
        #we have NA for the iterations where the number of thresholds is slower
        #the NAs will be always at the end? if you have data at 0.0001 for an iteration, you should have data at 0.001 in that iteration because signicaint snps with the more stringent threshold should be also significant with the less stringent one. So, the NAs should be at the end of the list.

    #add the column names to the DataFrame
    linear_metrics_df.columns = column_names[0]
        #we use column_names[0] because we may have several iterations with the max number of names and hence we would have a list of tuples with the same sames, so we just take the first one

    #count the number of NAs in each column and select those columns with more than 1 NAs
    linear_metrics_df = linear_metrics_df.loc[:, ~(linear_metrics_df.isna().sum() > 1)]
        #If a column has more than 1 NAs, it means that this threshold has not been used for at least all iterations and the corresponding median would not be based on 100 iterations, so we discart that threshold.
        #If we have 1 or less NAs, it means that we have at least 99 iterations with that threshold and we can use the median of those iterations to calculate the percentiles.

    #calculate percentiles for each column
    percentiles = [2.5, 50, 97.5]
    #p=2.5
    results_linear = pd.DataFrame(
        {f"{p}th_percentile": linear_metrics_df.quantile(p/100) for p in percentiles}
    )
        #with a comprehension, calculate the three percentiles for each column using "elastic_metrics_df.quantile()". 
        #each percentile is store as a key/value pair in the dictionary, where the key is the name of the column with the percentile and the value is the result of the quantile calculation for that column.
        #the dict is converted to a DF where each key (percentile) is a column

    #split the index to get thresholds and metrics
    split_index = results_linear.index.str.split("_")

    #extract metrics and the thresholds from the index
    #i=split_index[1]
    results_linear["metric"] = ["_".join(i[0:-1]) for i in split_index]
        #for each splitted index, get all elements execpt the last one (which is the threshold) and join them with "_"
    results_linear["threshold"] = [i[-1] for i in split_index]
        #just the last element

    #iterate over unique thresholds
    #threshold=results_linear["threshold"].unique()[0]
    for threshold in results_linear["threshold"].unique():
        
        #filter rows for the current threshold
        subset_threshold = results_linear[results_linear["threshold"] == threshold]
        
        #flatten the rows into a single row, i.e., all metrics of a given phenotype and dataset type are in a single row
        #this empty dict will store the values of each metric and percentile as a different key
        linear_row = {}

        #add the phenotype, dataset set, model type and threshold to the empty dict
        linear_row["phenotype"] = response_variable
        linear_row["dataset_type"] = dataset_type
        linear_row["model_type"] = "linear"
        linear_row["threshold"] = threshold

        #iterate across rows
            #"_," is used as a throwaway variable. It indicates that the value it holds is not going to be used in the loop.
        #row_data=subset_threshold.iloc[0,:]
        for _, row_data in subset_threshold.iterrows():
            
            #get the metric name
            metric = row_data["metric"]

            #save the percentiles
            linear_row[f"{metric}_2.5th_percentile"] = row_data["2.5th_percentile"]
            linear_row[f"{metric}_50th_percentile"] = row_data["50th_percentile"]
            linear_row[f"{metric}_97.5th_percentile"] = row_data["97.5th_percentile"]

        #add the columns for percentiles of heritability metrics
        linear_row[f"RHE_2.5th_percentile"] = np.nan
        linear_row[f"RHE_50th_percentile"] = np.nan
        linear_row[f"RHE_97.5th_percentile"] = np.nan
        linear_row[f"MCREML_2.5th_percentile"] = np.nan
        linear_row[f"MCREML_50th_percentile"] = np.nan
        linear_row[f"MCREML_97.5th_percentile"] = np.nan
            #these should be NA because we do not have heritability metrics for linear models
                #It seems the heritability can be calculated from linear models in ldak using the Human Default Model (./ldak.out --sum-hers snpher1 --summary quant.summaries --tagfile HumDef.tagging), but we are not doing it now, just we will use the estimations from elastic models.
            #we still need them to be able to concatenate the DF we are going to generate from this dict with the DF we generated for elastic models

        #convert the dict to a DF (each key as a column) and the concatenate to the eval metrics
        eval_metrics = pd.concat([eval_metrics, pd.DataFrame([linear_row])])

print_text("print the results for the correlation between PRSs and phenotypes", header=4)
print(eval_metrics.iloc[:, 0:7])

print_text("print the results for heritability metrics", header=4)
print(eval_metrics.iloc[:, -6:])

print_text("check number rows and columns", header=4)
#IF:
    #the number of rows IS NOT the number of phenotypes times number of dataset types times 7 (1 for elastic and 6 for linear as we have 6 thresholds)
        #our thresholds are 1, 0.1, 0.01, 0.001, 0.0001 and 0.00001. The last threshold (0.000001) is not used never becuase there are no SNPs below that threhold for any phenotype
    #OR
    #the number of columns IS NOT 4 evaluation metrics and 2 heritability metrics times 3 percentiles plus 4 columns indicating phenotype, dataset type, model type and threshold
if ( \
    (eval_metrics.shape[0]!=len(response_variables)*len(dataset_types)*7) | \
    (eval_metrics.shape[1]!=6*3+4) \
):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE FINAL TABLE")

print_text("save results", header=4)
eval_metrics.to_csv( \
    "./results/final_results/evaluation_metrics.tsv", \
    sep="\t", \
    header=True, \
    index=False, \
    na_rep="NA" \
)






##########################
# region Overlap top PRS #
##########################
print_text("Overlap top PRS", header=1)
print("""
It was suggested we can check the overlap across iterations for the samples with the highest PRS values in the training or test set.

We could do this by checking the top 1% of samples in the training set and the test set for each iteration and then calculate wath percentage of sample is the same across ALL iterations.
      
The problem here is that we have several phenotypes, datasets, models, thresholds.... so we should do it for each combination of these variables. Considering this and the fact that we already have several evaluation metrics that specifically target the predictive power of the PRSs, I think it is not worth to do this analysis.
""")






#########################
# region Interpretation #
#########################

#Relevant links
    #slides from Dr. Speed
        #https://dougspeed.com/short/
    #O´reilly paper
        #https://www.nature.com/articles/s41596-020-0353-1

#Clinical utility of PRSs
    #“Polygenic risk scores: from research tools to clinical instruments”, 
        #Lewis and Vassos, Genome Medicine, 2020 (figure below) 
    #“Could Polygenic Risk Scores Be Useful in Psychiatry? A Review”, 
        #Murray et al., JAMA Psychiatry, 2020 (figure on next page)

"""
- I have calculated the correlation between the PRSs and phenotype residuals in 100 different test sets (i.e., 100 different training/test partitions) and then calculated the 95CI across the 100 values. For most phenotypes, --elastic shows a higher correlation compared to --linear, but the thing is that most correlations are pretty low. The median across iterations is between 0.03 and 0.08 and, in all cases expect one phenotype, the 95CI overlaps with zero. Even for the best phenotype, the percentile 2.5th is almost at zero (0.0027). Note that I ran this using the new PCAs based on the clean data, and also checked that adding more covariates into the model (i.e., add all PCAs) does not have an influence. Would you say this level of correlation preclude to consider these scores informative for the phenotypes under study?
    - Dr. Speed: yes, these seem weak scores (but good you were careful and used 100 replicates). You can get an idea friom the screen output (or .hers file) whether this is purely due to your low sample size, or also reflects that the traits have small heritability

- Using --elastic, I got a median heritability around 0.7-0.8 according to both RHE and MCREML. So I guess this means that the heritability of the phenotype is high but we do not have enough power due to low sample size, right? Or do you think we cannot conclude even that because the low sample size could make the h2 estimates too noisy?
"""

# endregion






print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/association_studies/
#chmod +x ./scripts/production_scripts/03d_processing_results.py
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03d_processing_results.py > ./03d_processing_results.out 2>&1
#grep -Ei 'error|false|fail' ./03d_processing_results.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.

