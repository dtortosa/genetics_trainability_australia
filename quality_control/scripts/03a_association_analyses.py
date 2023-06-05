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



######################################
######## ASSOCIATION ANALYSIS ########
######################################

#This script will perform association analyses. 

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

### APPROACH
#VO2 max trainibility
    #Make a Polygenic predictor score of VO2 max trainibility
        #Obtain effect sizes and p-values from previous genome-wide studies done on VO2 max trainibility. This is important to avoid using the same data to derive the score and validate, if not we could overestimate its predictive power of the score.
            #Genes to predict VO2max trainability: a systematic review
                #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
            #Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
                #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7#Sec21
            #Genotype-Phenotype Models Predicting VO2max Response to High-Intensity Interval Training in Physically Inactive Chinese
        #if we have data from different gwas we can obtain meta-p-values or something like that
            #Deep learning-based polygenic risk analysis for Alzheimer’s disease prediction
        #use this information to create algorithms that predict trainability on our study.
            #We could use simple techniques like P-value prunning, or use machine learning from linear models to neural networks
        #we train the algorithm and evaluate it with our data, so we can have a measure of predictive power for trainability
        #In summary, we would take advantage of your experiment and machine learning to use previous genetic knowledge in the prediction of trainability.
        #problems
            #the GWAS studies performed by Bouchard et al do not report p-values or effect sizes. 
            #I have found an intervention study, in which David is indeed a co-author, that says the data could be shared if requested. I am talking about the Predict-HIIT study. We could use their data to develop different polygenic predictive scores, not only using their approach, but also others (e.g., machine learning), and validate the different scores in our cohort. I know that the protocols and the populations used are different, but we do not have the same protocol applied to other cohort with also genetic data...
                #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
        #solutions
            #obtain effects sizes of variants from GWAS studies of VO2 max which I think there are more. It is not ideal to use a base line instead of change, but I have seen other studies doing similar things with BMI and weight loss. 
                #Polygenic risk score for predicting weight loss after bariatric surgery
                    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171810/
            #just use VO2 max directly
    #Discover new variants associated with trainability
        #I would not got for this option as we do not have a validation cohort. 
        #we would need genome-wide data from an intervention where VO2 was measured.
        #PPS more interesting because we are directly predicting trainability.
    #Make a Polygenic predictor score of VO2 max
        #It is less interesting, but it would be easier because we could use the previous GWAS studies to obtain effect sizes, then train in our dataset and even validate in a third dataset, as the UK biobank has VO2 max data along with genetics.








###IMPORTANT, THINK WHAT chromosomes include, X, Y, PAR and mito should be analyzed?
    #I guess in independent associations is ok, an SNP in X is not affecting chromosome 1,
        #but if you use any tool that consider all snps, then they can influence
    #but what about polygenic scores? there all snps are combined!

###GWAS VALIDATION
#the gold standard of validation is to use an independent and diverse cohort, but in cases where it is difficult to find it (like in rare diseases or our case, i.e., a experimental manipulation with not so common traits), internal validation is used
    #Evaluation of a two-step iterative resampling procedure for internal validation of genome-wide association studies
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4859941/
    #In "Genome-wide association studies", nature review methods, they say that internal validation is a possibility, but no the gold standard
        #https://www.nature.com/articles/s43586-021-00056-9

##POLY SCORE CALC
#you can do polygenic score with your own data? CHECK this before making a deciison
    #I guess we cannot do GWAS, validation of GWAS and PRS with the same datasset, too many splits for a small dataset... 
        #"Accuracy measurements can be inflated if the discovery GWAS and the target cohort share individuals"
            #Genome-wide association studies
                #https://www.nature.com/articles/s43586-021-00056-9
    #we could do PRS if we find summary statistics for VO2 max for our snps in previous GWAS, so they can be used as input for PRS, deep learning (like in "Deep learning-based polygenic risk analysis for Alzheimer’s disease prediction" )
    #Tutorial: a guide to performing polygenic risk score analyses, nsture, o reilly
        #https://www.nature.com/articles/s41596-020-0353-1
