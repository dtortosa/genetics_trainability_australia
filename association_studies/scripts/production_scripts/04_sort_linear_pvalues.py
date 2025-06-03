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



########################################################
######## GET THE FIRST 10K SNPS ORDER BY PVALUE ########
########################################################

#David needs the first 10K SNPs according linear p-value association 






#######################################
# region INITIAL ANOTATIONS AND STEPS #
#######################################


###########
# imports #
###########

import pandas as pd
import numpy as np
import argparse
from multiprocessing import Pool
from functools import partial



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
    elif ("warning" in complete_process.stderr):
        
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






#################################
######## DEFINE FUNCTION ########
#################################

#define the function
#response_variable="distance_change"
def sorting_pvalues(response_variable):

    print_text("starting " + response_variable, header=1)
    print_text("load assoc results to pandas", header=2)
    assoc_results = pd.read_csv( \
        "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_linear_raw.assoc.gz", \
        sep="\t", \
        header=0, \
        low_memory=False)
    print(assoc_results)

    print_text("sort in increaseing order by Wald_P", header=2)
    assoc_results_sorted = assoc_results.sort_values(by="Wald_P", ascending=True)

    print_text("check we correlectly sorted", header=2)
    if (assoc_results_sorted["Wald_P"].equals(assoc_results["Wald_P"].sort_values(ascending=True))) & (assoc_results_sorted["Wald_P"].min() == assoc_results["Wald_P"].min()) :
        print(assoc_results_sorted)
    else:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM SORTING THE DATA")

    print_text("get the first 10K SNPs", header=2)
    assoc_results_sorted_10k = assoc_results_sorted.head(10000)
    
    print_text("check the number of rows", header=2)
    if assoc_results_sorted_10k.shape[0]!=10000:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM SUBSETING THE DATA")

    print_text("save the results", header=2)
    assoc_results_sorted_10k.to_csv( \
        "./results/final_results/analysis_full_data/" + response_variable + "/" + response_variable + "_small_set_predictors_set_linear_raw_sorted_10k.assoc.gz", \
        sep="\t", \
        index=False, \
        compression="gzip" \
    )






#################################
######## RUN THE FUCNTION #######
#################################

#phenotype="distance_change"
for phenotype in ["distance_change", "beep_change", "vo2_change", "weight_change"]:
    sorting_pvalues(phenotype)






print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/association_studies/
#chmod +x ./scripts/production_scripts/04_sort_linear_pvalues.py
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/04_sort_linear_pvalues.py > ./04_sort_linear_pvalues.out 2>&1
#grep -Ei 'error|false|fail' ./04_sort_linear_pvalues.out
