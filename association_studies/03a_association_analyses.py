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
    #5. Tutorial: a guide to performing polygenic risk score analyses
        #https://www.nature.com/articles/s41596-020-0353-1
        #https://choishingwan.github.io/PRS-Tutorial/
    #6. Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
    #7. Omics Data Preprocessing for Machine Learning: A Case Study in Childhood Obesity
        #https://www.mdpi.com/2073-4425/14/2/248



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



#################################
# approach to calculate the PRS #
#################################

#beside our interesting alzerhemie paper with neural nets
    #other option
        #SNPRS:Stacked Neural network for predicting Polygenic Risk Score

#VO2 max trainibility
    #Make a Polygenic predictor score (PPS) of VO2 max trainibility
        #In order to calculate a PPS, we need effect sizes and p-values for the association between SNPs and the trait obtained from a different study than the one used to train the scores, so we avoid overestimating its predictive power.
            #"Accuracy measurements can be inflated if the discovery GWAS and the target cohort share individuals"
                #Genome-wide association studies
                    #https://www.nature.com/articles/s43586-021-00056-9
            #if the score has a small effect size but it is associated and robust to confounders, you can at least say that there is a connection between genetics and the trait
                #"However, if the results are shown to be robust to confounding (see Population genetic structure and the generalizability of PRSs), then the effect size is not important if the aim is only to establish whether an association exists, which may provide etiological insight."
                    #Interpretation of PRS-trait associations
                        #https://www.nature.com/articles/s41596-020-0353-1
        #I have checked previous GWAS about VO2 max trainability and they all have much lower sample size than this study and hence, less power to detect true associations than us. Therefore, I think makes sense to use this study to associate genes and trainability rather than use p-values from previous GWAS to calculate the scores.
        #The problem is that we then need at least a second cohort to train the score and see how they predict. 
            #In the PRS tutorial, they say that you need at least two cohort, the base to obtain GWAS summary statistics and then the target to do cross-validation with the PRS. 
            #Of course, a test set would be then necessary for truly real predictive power, but in the absence of that, we could just use CV.
                #Overfitting in PRS-trait association testing
                    #https://www.nature.com/articles/s41596-020-0353-1
            #maybe Improve-HIIT could be used as a test set?
        #David, I have seen you are co-author of one of the previous GWAS on V02 max trainability, the HIIT-Predict study. Do you think there would have been any possibility to talk with them about using their data for training polygenic scores purposes? They say in the paper data is avaiable upon reasonable request.
            #we should check that samples from both cohorts are not related, as both come from Australia
            #maybe Improve-HIIT could be used as a test set?
            #Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
                #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7#Sec21
        #maybe we could even use our PPS developed in the cohort to test whether it interacts with the intervention type in HIIT?
            #In other words, try also models with the interaction PPSxTrainingType
            #maybe high PPS values associate with higher increase in VO2 max only in high volumen HIIT, while not in constant training. This would be a small window to training personalization based on genetics.
            #I am just speculating in this point because I am not sure if there is enough sample size in Predict-HIIT to do test interaction and also split data for training-evaluation.
        #In summary, we would take advantage of your experiment and machine learning to develop a new PPS of trainability that would be then tested.
    #Alternatives
        #just ask for p-values/betas from previous GWAS on VO2 max change like the HERITAGE of Predict-HIIT studies. 
            #This is not genotypes, just the summary statistics.
            #the problem is that we would be using studies with smaller sample size to discover variants. Indeed, the Predict-HIIT study did not have any SNP with a significant p-value after correction and the paper say that they are probably underpowered. Their power calculations indicates that "a cohort of 2960 samples would have 80% power to detect a quantitative trait with a true heritability of 30%.". Therefore, it would be better to use our study to discover, given its larger sample size.
            #if we have data from different gwas we can obtain meta-p-values or something like that
                #Deep learning-based polygenic risk analysis for Alzheimer’s disease prediction
            #previos GWAS about trainability
                #some GWAS of VO2 max trainability
                    #Genes to predict VO2max trainability: a systematic review
                        #https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4192-6
                    #Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
                        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
                    #Genotype-Phenotype Models Predicting VO2max Response to High-Intensity Interval Training in Physically Inactive Chinese
        #Obtain p-values/betas from GWAS of baseline VO2 max in the UK Biobank. 
            #This biobank includes estimates of VO2 max based on a submaximal cycle ramp test and it has a large sample size, which is good for power in the base cohort.
            #there is already a study that have done a GWAS on this data, so we could just ask them for the summary statistics. In the worst case scenario, we could do the GWAS on our own.
            #then we can test the polygenic scores in our cohort.
                #The genetic case for cardiorespiratory fitness as a clinical vital sign and the routine prescription of physical activity in healthcare
                    #https://www.medrxiv.org/content/10.1101/2020.12.08.20243337v2
            #problems
                #we would use a correlate of VO2 change, not the trait itself. this can be done, but it is not optimal
                    #See "Predicting different traits and exploiting multiple PRSs" 
                        #https://www.nature.com/articles/s41596-020-0353-1
                #we would not take advantage of the greater sample size we have in this study compared to previous interventions in order to derive the polygenic score.
        #Just do GWAS of trainability with our cohort
            #I would not got for this option as we do not have a validation cohort, only the discovery cohort.
            #In some cases like studies of rare diseases, they use interval validation of the GWAS, i.e., using the same cohort but applying resampling methods. 
                #This is an option but it would be a clear limitation of the study.
                #the gold standard in this point is to use an independent and diverse cohort for validation
            #info
                #Evaluation of a two-step iterative resampling procedure for internal validation of genome-wide association studies
                    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4859941/
                #In "Genome-wide association studies", nature review methods, they say that internal validation is a possibility, but no the gold standard
                    #https://www.nature.com/articles/s43586-021-00056-9
        #Make a Polygenic predictor score of VO2 max
            #It is less interesting, but it would be easier because we could use the previous GWAS studies to obtain effect sizes, then train in our dataset and even validate in a third dataset, as the UK biobank has VO2 max data along with genetics.


###you have 1057 men!! sex imbalance

###IMPORTANT, THINK WHAT chromosomes include, X, Y, PAR and mito should be analyzed?
    #I guess in independent associations is ok, an SNP in X is not affecting chromosome 1,
        #but if you use any tool that consider all snps, then they can influence
    #but what about polygenic scores? there all snps are combined!

#### add genotyping batch as a factor in the model???
    #done here:
    #https://academic.oup.com/cardiovascres/article/118/Supplement_1/cvac066.013/6605381