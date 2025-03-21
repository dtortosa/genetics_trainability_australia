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
######## CALCULATE FINAL PCAS ########
######################################

#we are going to re-calculate PCAs
    #we are going to do it with the final dataset after QC just before imputation as recommended by Dr. Speed
    #we have also going to include thinning to avoid very high loadings in one specific region of chromosome 8 for PCA 5, 4 and 3.

#I asked for advice to Dr. Speed regarding the selection of covariates in the context of his PRS tool. He basically said the same as Jonatan: we could apply a sensible criteria (e.g., select covariates significantly associated) and then check what happens doing the opposite, i.e., increasing the number of covariates to check our results are robust. In this context, we ended up discussing the weights (influence) of specific SNPs on the different PCAs as I mentioned that, for some PCAs, there was an abnormal accumulation of SNPs with very high weights in a narrow region of chromosome 8. What I originally did in the QC was just to discard these PCAs as covariates for future modeling, but Dr. Speed recommended recalculating the PCAs removing SNPs in high linkage-disequilibrium regions as these cause this kind of problem. This opened a can of worms as it could imply I should have to repeat all the QC steps, including imputation!

#Dr. Speed recommended me to not do that as long as the PCAs have not been used for anything else besides outlier removal. Even if the number of outliers slightly differ when using the new PCAs, that should be ok as imputation results should be robust to that. I have checked what happens if, in the QC, we would have removed samples using the "new PCAs" and the difference is that 8 additional samples would remain. So our current approach does not add problematic samples, just remove maybe a few more. Note that "saving" these samples would require repeating everything with the whole dataset.

#In summary, I do not think it makes sense to repeat all the QC steps and imputation due to this matter. As suggested by Dr. Speed, we can just add an additional step before the PRS calculation where we recalculate the PCAs with the clean set of SNPs and samples after QC, obtaining more robust and clean PCAs for modeling. The initial approach I followed had as consequence the removal of a few additional samples, which should not be a problem given our sample size (1008 final samples). 


#POR AQUIII


    #Is there any specific recommendation about the use of covariates for LDAK besides the need of splitting factors vs continuous covariates? I have the usual covariates (age, sex, body mass, phenotype baseline, PCAs...) and I am using only those covariates that are individually associated with the phenotype of interest to limit model complexity.
        #I examined this question about 15 years ago, to try and find a more justified solution (instead of simply "include 5 PCAs, age and sex"). I remember I experimented with only including PCs that were significantly associated with the phenotype, or only PCs with significant eigenvalues (a software called EigenSoft proposed a test for significance) . Howveer, my conclusion was that it is generally not too important (albeit I had larger sample sizes).
        #Then about 10 years ago, I switched to using heritability analysis to decide number of PCs - see the left and right analyses on this page https://dougspeed.com/quality-control/ - basically, include enough PCs that the sum of estimates from left and right halves of genome equal estimate from whole genome.
        #Nowdays, I normally include a lot (eg 10 or 20 PCs, maybe 10 internal PCs and 10 from 1000 GP), and if I wanted to check, I would repeat the analysis with fewer. Now again, this is larger sample sizes - the impact is larger for smaller sample sizes. But perhaps it remains a good strategy - I suggest you first perform the analysis using the best choice you can make (whether that is using 5 or 10, because another paper did similar, or choosing the significant ones only - provided you have a sensible criterium). Then repeat using a different choice as a sensitivity analysis (so if you included a lot first time, repeat with fewer). Provided the results are not very different (e.g., you have a SNP that is highly significant in one, but not in the other), then your choices should be good
    #Thank you for your detailed answer! Yes, that makes a lot of sense. I was indeed considering the p-values from eigensoft, specifically those axes with a p-value<0.01 (first 5 in my case). Then I checked the distribution of SNPs loadings across the genome and found a suspicious increase of positive and negative loadings in a narrow region of chromosome 8 that was very marked for PC5, but was also visible for PC4 and PC3. So I decided to use only PC1 and PC2 for subsequent analyses. Then, as I noted before, I included these two PCAs if they were associated with the phenotype. Maybe this is a very stringent approach, but we are using clear criteria and it could make sense for a small dataset. I will follow your advice and check what happens including more axes (e.g., all significant axes according to eigensoft, 10, etc...)
        #Yes, very good to check the loadings (so be suspicous if you find a very dense patch - similar to what you observed on chr 8). Similarly, check the PCs and be suspicoius if you find groupings (e.g, divide into 2 or 3) that do not seem to correlate with ancestry. Both of these normally indicate problems due to LD - that the loadings are maiinly concentrated on a region of highLD, so instead of picking up genome-wide variation, they just record the genotypes for that small region
        #If this occurs, it is best to repeat with more thinning - and if human, I recommend excluding long range ld regions https://dougspeed.com/high-ld-regions/
    #Amazing all these insights! Yes, I initially checked the existence of clusters and thankfully no clear cluster was present, but I completely missed the option to remove the SNPs with the high loadings instead of the PCAs showing the problem. After increasing thinning (plink2 --indep-pairwise from 500kb 0.2 to 600kb 0.1) the peak completely disappears and the rest of regions/PCAs look perfect in this regard. At this point, after performing all QC steps and imputation, do you think it would necessary to re-run all my pipeline from this step? or it would be reasonable to just recalculate the PCAs with the increased thinning (to be used as covariates) and remove the high-LD from the plink fileset just before running LDAK?
        #Thanks, that is good to hear
        #In terms of your pipeline, my quick answer would be that you only need to repeat the steps that used the PCs - from looking above, you did the following: 1.2K samples and 3M of genotyped+imputed SNPs after QC - then make PRS using --elastic
        #So here I dont think the imputation depends on the PCs, so it sounds fine to leave, especially as I imagine that is the slowest step (well, maybe you used PCs to decide which individuals to exclude - but imputation results should be robust to changing the individuals slightly).
        #Then for --elastic, you can repeat with the new PCs. Note that it is fine to include high ld regions for --elastic (ie, remove them when computing the PCS, to avoid the uneven loadings, but you can retain them when making the PRS)
    #Exactly, I am only using the PCAs as covariates and to exclude outliers (iterative removal with smartpca). I have checked what happens doing the removal with higher thinning and the difference is that 8 samples are retained, while with less thinning they are removed. I understand this should not have a great impact on imputation considering the total size of my cohort. So yes, I think I will calculate the PCAs avoiding high ld regions and then use them as covariates for --elastic (retaining in that case high-ld regions). Would you recommend to calculate the PCAs with the fileset generated after imputation and QC? or with the fileset I originally used for the PCAs during pre-imputation QC?
        #I think it should not matter - normally, I compute PCs using only high quality (eg directly genotyped ) snps, rather than also include imputed
        #But I think general advice is to use only the final list of samples (and not the original list of samples). A simple case is if you had a few ethnic outliers, that were excluded during qc, then these would dominate PCs computed using all samples, so it would probably be best to use only the samples after qc
    #Ok, so I can calculate the PCAs using the data generated at the end of QC (when all sample and SNP filters have been applied) but before the imputation. Thank you so much for your insights!

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



######################
# folder preparation #
######################
run_bash(" \
    mkdir -p ./data/final_pcas; \
    ls -l \
")

# endregion






#####################################################################
# region  #
#####################################################################


#take the last plink setfile before imputation
#do LD-subset increasing thinning
#also remove the areas indicated by Dr.Speed in his comment


# endregion