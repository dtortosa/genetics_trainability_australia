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



################################################
######## PRE-IMPUTATION QUALITY CONTROL ########
################################################

######################################################
# region INITIAL ANOTATIONS AND STEPS ################
######################################################


#This script will perform pre-imputation QC. 

#This and next scripts are based on the following tutorials:
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



###########################################
# conservative approach for performing QC #
###########################################

#we should be conservatives in the QC analysis because small errors can become inflated when aggregated across SNPs in PRS calculation.
    #"QC of base and target data" in:
        #https://www.nature.com/articles/s41596-020-0353-1

#It makes sense because when doing polygenic score, you are analyzing results of previous analyses, maybe several GWAS combined (meta-p-values..) so you can add many errors...

#this applies if our cohort is going to be used as base, thus it does not matter if we are going to use it as discovery because it can aggregate errors.

#if our cohort is going to be used as target, then it is final, the last step where we want to check how useful our PRS is with cross-validation, thus we still need to be conservative.



##########################################
# rationale of merging before imputation #
##########################################

#It is clear that we can do merging after imputation
    #From reading the recent GWAS protocol of the Ritchie lab. I understand that the usual approach is to do the imputation of different sources (e.g., batches) separately and then do the merging. In the section "Batches effect", they say that "Genotype imputation strategies now provide an opportunity to impute genotypes from multiple platforms to a reference genotyping platform. THEREFORE, DATA FROM DIFFERENT SOURCES CAN BE MERGED AFTER IMPUTATION."
        #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view?usp=sharing
    
    #In addition, in other paper they say 
        #"This process (imputation) is particularly important when combining or performing meta-analysis on data generated using multiple different genotyping platforms"
        #"After imputation and merging of the datasets, quality control procedures were implemented to create high quality, analysis-ready data set for genome-wide association studies."
        #https://www.frontiersin.org/articles/10.3389/fgene.2014.00370/full

#but it seems to be done for data from different studies
    #I understand that this is specially relevant when you have data from different sources, it is very important to do imputation because many SNPs are NOT going to be shared across sources, you are going to lose them unless you impute them to have the same SNPs across all panels.

    #In our particular case, this should not be very important because the same SNPs are genotyped in both batches, so it is very unlikely that all samples of a batch have missing for a given SNP so that SNP is not present in one batch but it is present in the other one. 

#We have a problem because our sample size is not very big, so if we split even more between batches and filter by missingness, MAF... it is more likely to lose more SNPs. Indeed, a filter of MAF<0.05, which is recommended for sample sizes < 1000, leads to the lost of half of SNPs in both batches!! We should try at least to combine the batches to gain more minor homozygotes for these low-MAF SNPs.

#It does not make sense to lose SNPs just because the merging was done after. If the batch effects are carefully checked, we should not have a problem. Both batches should be a random sample of samples, so they would be equally influenced by ancestry, etc...

#We can merge, then apply the filters required before PCA and then run the PCA to check for batch effects. We can also include the batch as a factor in the models as done in this paper:
    #New genetic determinants of VO2max-level identified by GWAS: The HUNT Study
        #https://academic.oup.com/cardiovascres/article/118/Supplement_1/cvac066.013/6605381



#############################################################
# final decision about the 1100JHJM, 1200JPJM and 7800AGSO #
#############################################################

#the three duplicated IDs
    #For 1100JHJM and 1200JPJM, David's postdoc agrees we should remove them as we have the same ID for two different rows in the excel file.
        #"This is correct. In session 844, there are two individuals with the same code (1100JHJM, both male), and in session 845, there are two individuals with the same code (1200JPJM, both male). I agree that we will need to remove these from our analysis."
    #For 7800AGSO, David's postdoc asked if I could know which is the sample from the excel file that doesn’t appear to have a DNA sample, but I cannot know that.
        #"There is only one individual (male) with code 7800 from session 878, so there appears to have been a labelling error while preparing the DNA samples. Are you able to tell us which is the sample from the excel file that doesn’t appear to have a DNA sample? That might help us to work this out. If not, I agree that we will need to remove this sample also."
            #We have genetic data for 1464 samples, while in the excel we have 1463 samples. So the problem is the other way around, all samples in excels have genetic data, but 1 sample in illumina does not have phenotypic data.
                #There are 41 of these samples with genetic but not phenotypic data.
            #7800AGSO in the excel could be 7800AGSO_1 or 7800AGSO_2, but I cannot know which of these is.
                #I do not have phenotypic information in illumina data, only the ID and then genetic data. The only information shared between the two datasets is the ID.

#AGRF confirmed that they included these samples (1100JHJM, 1200JPJM and 7800AGSO) in the reproducibility report but they are duplicates, not replicates. See 01b_illumina_report_to_plink.py for further details.
    #"As you have highlighted, these three pairs were not highlighted in the submission as technical replicates, but the pairs were assigned the same sample name for the sample manifest and tubes provided"

#Therefore, these three pairs of samples should be removed. They have been removed by 01b_illumina_report_to_plink.py.



########################
# notes about 2397LDJA #
########################

#I have an Illumina report for one sample (2397LDJA; ILGSA24-17303) that is not present in the excel file. In contrast, the sample ID 2399LDJA is present in the excel file but there is no illumina report for this ID. Both IDs are the same except for the 4th digit. Maybe they are the same sample, but it would probably be a good idea to remove both IDs.
    #"Yes! I had exactly the same note, I think it is a mislabelling of the last digit of the number (the labelling was very hard to read on some of the blood samples). So, I think 2397LDJA; ILGSA24-17303 is 2399LDJA in the excel file"



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



#######################################
# Passing arguments of python program #
#######################################

#define input arguments to be passed when running this script in bash
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--n_cores", type=int, default=4, help="Number of cores/threads requested. Integer always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
n_cores = args.n_cores




############################
# starting with the script #
############################
print_text("starting with the pre-imputation QC using " + str(n_cores) + " cores ", header=1)

# endregion




##############################################
# region MERGING BATCHES #####################
##############################################
print_text("merging batches", header=2)
print_text("check we have 216 and 1242 samples for batch 1 and 2 (batch 2 lost 6 samples that were duplicated), respectively looking at the merged fam file of each batch and the list of samples to be merged", header=3)
run_bash(" \
    n_samples_first_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/ILGSA24-17303/03_merged_data/ILGSA24-17303_merged_data.fam.gz | \
        awk 'END{print NR}'); \
    n_samples_second_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/ILGSA24-17873/03_merged_data/ILGSA24-17873_merged_data.fam.gz | \
        awk 'END{print NR}'); \
    if [[ $n_samples_first_batch -eq 216 && $n_samples_second_batch -eq 1242 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #decompress the fam file of each merged file, i.e., for each batch, then count the number of lines and save the output
        #get the number of the row at the end, i.e., the total number of rows
            #use awk instead wc -l because "wc -l" counts "new line" character, and it is possible to have a line without that character. If at the end of the file, you do not have an empty line, your actual last line will not have "new line" character so it will not be counted by wc -l.
                #https://unix.stackexchange.com/a/362280
                #https://stackoverflow.com/a/28038682/12772630
        #if the number of lines (samples) in the first and the second batch is 216 and 1242, respectively, then OK
            #remember that we lose 6 duplicated samples in the second batch.
            #note that "==" is only for strings, you have to use "-eq" for numbers
                #https://stackoverflow.com/questions/20449543/shell-equality-operators-eq
run_bash(" \
    n_samples_first_batch=$( \
        awk 'END{print NR}' ./data/genetic_data/plink_bed_files/ILGSA24-17303/02_data_to_merge/list_to_merge.txt); \
    n_samples_second_batch=$( \
        awk 'END{print NR}' ./data/genetic_data/plink_bed_files/ILGSA24-17873/02_data_to_merge/list_to_merge.txt); \
    if [[ $n_samples_first_batch -eq 216 && $n_samples_second_batch -eq 1242 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #count the number of rows in the file with the list of samples to be merged in both batches. We use NR (number of records) of awk instead of wc -l because the file does not end with an empty line, so wc does not count the last one, and we get 1 row less
            #https://stackoverflow.com/questions/32481877/what-are-nr-and-fnr-and-what-does-nr-fnr-imply
            #https://stackoverflow.com/questions/12616039/wc-command-of-mac-showing-one-less-result



print_text("create new folders to store files", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/; \
    mkdir \
        --parents \
        merged_batches/data_to_merge; \
    mkdir \
        --parents \
        merged_batches/merged_plink_files; \
    ls -l ./merged_batches")
        #--parents to make folders recursively and avoid errors if the folder already exists



print_text("copy the plink files of both batches", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files; \
    cp ./ILGSA24-17303/03_merged_data/ILGSA24-17303_merged_data* ./merged_batches/data_to_merge/; \
    cp ./ILGSA24-17873/03_merged_data/ILGSA24-17873_merged_data* ./merged_batches/data_to_merge/; \
    ls -l ./merged_batches/data_to_merge")



print_text("create a .txt with the name of the plink inputs for each batch", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/; \
    echo ILGSA24-17303_merged_data > list_files_to_merge.txt; \
    echo ILGSA24-17873_merged_data >> list_files_to_merge.txt; \
    ls -l; \
    cat list_files_to_merge.txt")
        #for --merge-list
            #If a line contains only one name, it is assumed to be the prefix for a binary fileset
            #https://www.cog-genomics.org/plink/1.9/data#merge_list
        #">" creates a new file where the output of the command it is saved
        #">>" appends the output of the command to an existing file
            #https://unix.stackexchange.com/questions/159513/what-are-the-shells-control-and-redirection-operators
            #https://unix.stackexchange.com/questions/77277/how-to-append-multiple-lines-to-a-file



print_text("decompress the plink inputs", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/; \
    gunzip \
        --keep \
        --force \
        ./ILGSA24-17303_merged_data*.gz; \
    gunzip \
        --keep \
        --force \
        ./ILGSA24-17873_merged_data*.gz; \
    ls -l")
        #-keep to keep the original compressed file and --force to overwrite if the decompressed file exists



print_text("merge the two batches using two different approaches in plink: --merge-list and --bmerge", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge; \
    plink \
        --merge-list ./list_files_to_merge.txt \
        --out ../merged_plink_files/merged_batches")
        #we can use --merge-list to merge batches listed in a .txt file.
            #https://www.cog-genomics.org/plink/1.9/data#merge_list
        #you can first call a fileset (used as reference) doing "plink --bfile NAME_REFERENCE" and then add "--merge-list LIST_TO_MERGE" where LIST_TO_MERGE has the names of the filesets to be merged with the reference fileset.
            #https://www.biostars.org/p/9467886/
        #alternatively, this can be used without a reference, i.e., all files are included in LIST_TO_MERGE. In that case, the newly created fileset is then treated as the reference by most other PLINK operations.
        #we are using the second option as we did for merging plink files of each sample within each batch.
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge; \
    plink \
        --bfile ILGSA24-17303_merged_data\
        --bmerge ILGSA24-17873_merged_data \
        --out ../merged_plink_files/merged_batches_2")
        #we can also call first one of the batches with "--bfile" and use it as reference without using --merge-list, then add the second batch with --bmerge. We can directly use the prefix of the fileset, because the three files (bed, bim and fam) should be named the same way except for the extension.
        #we get a warning about SNPs with the same position, but this is ok because later we will remove duplicates by position and allele code.



print_text("check that merging the two batches with --merge-list and --bmerge gives the same", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    bed_status=$(cmp --silent merged_batches.bed merged_batches_2.bed; echo $?); \
    bim_status=$(cmp --silent merged_batches.bim merged_batches_2.bim; echo $?); \
    fam_status=$(cmp --silent merged_batches.fam merged_batches_2.fam; echo $?); \
    if [[ $bed_status -eq 0 && $bim_status -eq 0 && $fam_status -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #check that files (bed, bim, fam) obtained with --merged-list (merged_batches) and --bmerge (merged_batches_2) are the same.
        #check byte by byte whether the two files are the same
            #cmp takes two files and compare them until 1 byte is different
            #we make it silent and get the final status
            #remember that "$?" gives the return value of the last run command.
                #For example, 
                    #ls somefile
                    #echo $?
                    #If somefile exists (regardless whether it is a file or directory), you will get the return value thrown by the ls command, which should be 0 (default "success" return value). If it doesn't exist, you should get a number other then 0. The exact number depends on the program.
                #https://stackoverflow.com/a/6834572/12772630
            #the return value of cmp will be 0 if the two files are identical, if not, then we have differences between the files
                #https://stackoverflow.com/a/53529649/12772630
        #if the status of the three comparisons is zero, then perfect.



print_text("remove the second fileset created only for the check", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    rm merged_batches_2*; \
    n_files=$(ls merged_batches* | wc -l); \
    if [[ $n_files -eq 7 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")



print_text("check that the number of unique sample IDs in the merged fam file is equal to the total number of samples, i.e., no ID is duplicated", header=3)
run_bash("\
    n_uniq_ids=$(\
        awk \
            'BEGIN{ \
                FS=\"\t\"; \
                OFS=\"\t\"};\
            { \
                print $2\
            }' \
            ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.fam | \
        uniq | \
        wc -l); \
    n_samples_batch_1=$(\
        awk \
            'BEGIN{ \
                FS=\"\t\"; \
                OFS=\"\t\"};\
            { \
                if($1==\"combat_ILGSA24-17303\"){ \
                    print $0 \
                } \
            }' \
            ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.fam | \
        wc -l); \
    n_samples_batch_2=$(\
        awk \
            'BEGIN{ \
                FS=\"\t\"; \
                OFS=\"\t\"};\
            { \
                if($1==\"combat_ILGSA24-17873\"){ \
                    print $0 \
                } \
            }' \
            ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.fam | \
        wc -l); \
    total_number_samples=$(($n_samples_batch_1+$n_samples_batch_2));\
    if [[ $total_number_samples -eq 1458 && $n_uniq_ids -eq $total_number_samples && n_samples_batch_1 -eq 216 && n_samples_batch_2 -eq 1242 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #load the fam file with awk using tabs for the input and also for the output, then select the second column, IDs, get uniq cases and count them.
        #also count the number of rows for which the first column, family ID, is the name of the first and the second batch, count cases.
        #sum the number of samples of both batches using shell arithmetic expansion ("$(())")
            #https://phoenixnap.com/kb/bash-math
        #Each bath has 216 and 1242 samples, respectively, totalling 1458 samples. Therefore, the number of of unique IDs should be also 1458.
            #Remember that we lose 6 duplicated samples in the second batch.



## more info about merging ##
    #The new fileset plink.bed + .bim + .fam is automatically created in the process. (Corner case exception: if "--recode lgen" is part of the same run, the prefix is plink-merge instead.) Thus, it is no longer necessary to combine --merge with --make-bed if you aren't simultaneously applying some filters to the data.
        #https://www.cog-genomics.org/plink/1.9/data#merge
    #The order of sample IDs in the new fileset can now be controlled with --indiv-sort. 
        #the default is 'natural'/'n'
            #"Natural sort" of family and within-family IDs, similar to the logic used in macOS and Windows file selection dialogs; e.g. 'id2' < 'ID3' < 'id10'. This is the PLINK 1.9 default when merging datasets.
            #this is what we want, because we have actually a pattern in the IDs, first numbers and then letters. According the manual, if the IDs mix letters and numbers randomly, then other modes are better, but this is not our case.
                #Therefore, first samples of the first batch and then samples of the second batch, as the number ID of the second batch is larger (sort between families, which is the batch in our case). Within batch, we also want natural order by numberic and alphabetic.
            #https://www.cog-genomics.org/plink/1.9/data#indiv_sort
    #the default merging mode is (default) Ignore missing calls, otherwise set mismatches to missing.
        #I guess the default mode will set as missing those SNPs that are present in one batch but not in other.
        #In an example (https://www.biostars.org/p/496434/), this guy select shared SNPs between the two batches and then do the merge with --merge-list, but in our case we have the same number of SNPs in both batches (checked in the next lines), so we should not have problems with this.
    #If binary merging fails because at least one variant would have more than two alleles, a list of offending variant(s) will be written to plink.missnp
        #we do not have a plink.missnp file, so we should be good with this.



print_text("check we do NOT have a .missnp file, where variants with more than two alleles are indicated. We should not have this file", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    count_miss_file=$(ls -l *.missnp | wc -l); \
    if [[ $count_miss_file -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #list files with extension ".missnp" and count them. This should be zero.



print_text("check that the bim files (SNP maps) of both batches and the merged file have the same SNPs, although the minor can be different", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    awk -F ' ' '{print $1, $2, $3, $4}' ./merged_plink_files/merged_batches.bim > merged_columns.txt; \
    awk -F ' ' '{print $1, $2, $3, $4}' ./data_to_merge/ILGSA24-17303_merged_data.bim > batch_1_columns.txt; \
    awk -F ' ' '{print $1, $2, $3, $4}' ./data_to_merge/ILGSA24-17873_merged_data.bim > batch_2_columns.txt; \
    bim_check_1=$(cmp --silent merged_columns.txt batch_1_columns.txt; echo $?); \
    bim_check_2=$(cmp --silent merged_columns.txt batch_2_columns.txt; echo $?); \
    bim_check_3=$(cmp --silent batch_1_columns.txt batch_2_columns.txt; echo $?); \
    if [[ $bim_check_1 -eq 0 && $bim_check_2 -eq 0 && $bim_check_3 -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi; \
    rm merged_columns.txt; rm batch_1_columns.txt; rm batch_2_columns.txt")
        #select the first 4 columns of the three bim files using awk (-F is the delimiter, space in this case). Save the result as three different .txt files.
            #the 4 first columns are the chromosome code, Variant identifier, Position in morgans or centimorgans (safe to use dummy value of '0'), Base-pair coordinate (1-based).
            #these columns should be the same in both batches, the difference can be in the last two columns with the name of the minor and major alleles, because it would be possible that one allele is not present in any of the samples of a batch, but in the other batch it has both alleles (see below).
        #compare the merged file with both batches and then the two batches between them using cmp for that.
            #cmp takes two files and compare them until 1 byte is different
            #we make it silent and get the final status
            #remember that "$?" gives the return value of the last run command.
                #For example, 
                    #ls somefile
                    #echo $?
                    #If somefile exists (regardless whether it is a file or directory), you will get the return value thrown by the ls command, which should be 0 (default "success" return value). If it doesn't exist, you should get a number other then 0. The exact number depends on the program.
                #https://stackoverflow.com/a/6834572/12772630
            #the return value of cmp will be 0 if the two files are identical, if not, then we have differences between the files
                #https://stackoverflow.com/a/53529649/12772630
        #if the status of the three comparisons is zero, then perfect.
        #finally remove the files generated for this check



## note about the differences in minor/major alleles between batches ##
    #the four first columns are identical between batches, but the minor/major alleles are different.
    #the first difference between the two batches is a SNP in row 140. It has only C in first bath (monomorphic), but it has also T in the second. Accordingly the merged batch has both alleles, because this included first and second batch
        #first batch
            #0       GSA-10:47138541 0       0       0       C
        #second batch
            #0       GSA-10:47138541 0       0       T       C
        #merged file
            #0       GSA-10:47138541 0       0       T       C
    #opposite case with a snp in row 165, which is the first line of difference between merged and the second batch. It is monomorphic in second batch, but it has two alleles in first batch, so when merging, you get two alleles, not 1.
        #first batch
            #0       GSA-rs145797772 0       0       A       C
        #second batch
            #0       GSA-rs145797772 0       0       0       C
        #merged file
            #0       GSA-rs145797772 0       0       A       C



print_text("see if the first SNPs that have different minor/major allele between batches are indeed different in terms of frequency and alleles according to plink --freq", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    plink \
        --bfile ./data_to_merge/ILGSA24-17303_merged_data \
        --freq; \
    sed -n -e 1p -e 141p plink.frq; \
    sed -n 166p plink.frq;")
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    plink \
        --bfile ./data_to_merge/ILGSA24-17873_merged_data \
        --freq; \
    sed -n -e 1p -e 141p plink.frq; \
    sed -n 166p plink.frq;")
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    plink \
        --bfile ./merged_plink_files/merged_batches\
        --freq; \
    sed -n -e 1p -e 141p plink.frq; \
    sed -n 166p plink.frq;")
        #Calculate minor allele frequency of all SNPs in each batch using plink --freq
            #By itself, --freq writes a minor allele frequency report to plink.frq. 
                #https://www.cog-genomics.org/plink/1.9/basic_stats
        #We use sed -n to show specific lines. -e lets you to select several lines. 
            #note that line 140 in the bim file is 141 in this file because we have header in freq file, but no header is present in the bim file.
            #https://unix.stackexchange.com/questions/288521/with-the-linux-cat-command-how-do-i-show-only-certain-lines-by-number
        #the merged file has updated information about the minor/major alleles considering all samples, because of this, these two snps are biallelic in the merged file when in one of the batches is monomorphic, and the MAF changes when merging also, because now you have more samples.



print_text("remove freq files", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    rm plink.frq; \
    rm plink.hh; \
    rm plink.log; \
    rm plink.nosex; \
    ls -l")



print_text("compress the merged files", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/; \
    gzip \
        --force \
        ./merged_batches.bed; \
    gzip \
        --force \
        ./merged_batches.bim; \
    gzip \
        --force \
        ./merged_batches.fam; \
    ls -l")
        #--force: force overwrite of output file and compress links


print_text("remove the decompressed input files", header=2)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/;\
    rm ./ILGSA24-17303_merged_data.bed; \
    rm ./ILGSA24-17303_merged_data.bim; \
    rm ./ILGSA24-17303_merged_data.fam; \
    rm ./ILGSA24-17873_merged_data.bed; \
    rm ./ILGSA24-17873_merged_data.bim; \
    rm ./ILGSA24-17873_merged_data.fam; \
    ls -l")



print_text("check again we have the exact sum of samples from both batches", header=3)
run_bash("\
    n_samples_merged=$(\
        gunzip -c ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.fam.gz | \
        wc -l);\
    n_samples_first_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/ILGSA24-17303_merged_data.fam.gz | \
        wc -l); \
    n_samples_second_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/ILGSA24-17873_merged_data.fam.gz | \
        wc -l); \
    sum_samples=$((n_samples_first_batch+$n_samples_second_batch)); \
    if [[ $n_samples_merged -eq $sum_samples ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")



print_text("check again after merging we still have 654027 snps, and this is the number of variants than in each of the batches", header=3)
run_bash("\
    n_snps_merged=$(\
        gunzip -c ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.bim.gz | \
        wc -l);\
    n_snps_first_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/ILGSA24-17303_merged_data.bim.gz | \
        wc -l); \
    n_snps_second_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/ILGSA24-17873_merged_data.bim.gz | \
        wc -l); \
    if [[ $n_snps_merged -eq 654027 && $n_snps_merged -eq $n_snps_first_batch && $n_snps_merged -eq $n_snps_second_batch ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")




print_text("do some checks on the plink files of the batch", header=2)
print_text("perform missing report with the merged data", header=3)
print_text("calculate missing reports of the plink file sets generated in the previous step", header=4)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    gunzip \
        --keep \
        --force \
        ./merged_batches*;\
    plink \
        --bfile ./merged_batches \
        --missing; \
    ls -l")
        #gunzip all the files of the file set
        #use them as input to calculate the missing report
            #--missing produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.
                #https://www.cog-genomics.org/plink/1.9/basic_stats#missing


print_text("see first rows of the missing report for samples", header=4)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    awk \
        'BEGIN{FS=\" \"}{if(NR<10){print $0}}' \
        plink.imiss")
        #see lines below for details about missing report format
            #https://www.cog-genomics.org/plink/1.9/formats#imiss



print_text("check about the number of genotypes in the DNA report file", header=3)
#I have detected something odd in the DNA report of BOTH batches. For ALL samples, if you sum the number of genotypes with GS_score < 0.15, i.e., no calls, plus the number of called genotypes, we do not get 654027, but 650181. Therefore, we have 654027-650181=3846 missing variants.
#Note that both DNA reports indicate that the total number of loci is 654027, but then the sums are not 654027.
#This is caused by the SNPs where both the allele1 and allele2 are 0, see below.
print_text("load the bim file", header=4)
bim_file_before_any_filter = pd.read_csv( \
    "./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.bim", \
    sep="\t", \
    header=None, \
    low_memory=False)


print_text("check that the number of variants with 0 for the allele1 and allele2 is 3846", header=4)
n_zero_alleles = sum((bim_file_before_any_filter.loc[:,4]=="0") & (bim_file_before_any_filter.loc[:,5]=="0"))
if (n_zero_alleles == 3846):
    print("the number of SNPs with allele1 and allele1 as zero is 3846, which is exactly the number of SNPs lacking in the sum of SNPs called+non-called of the DNA report for each individual. This explains that this sum does not total to 654027. These cases with allele1 and allele2 as zero cannot have genotypes (there are no alleles), and without genotypes you cannot calculate the GS_score, which is used to determine whether a SNP is called or not.")
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE dna REPORT")
#Conclusion about missing vs DNA report: 
#I am not going to do more comparisons. The DNA report is not counting cases with allele1 and allele2 as zero, while the missing report it is. The missing report count every snp except those with hetero problems (heterozigous in non-PAR regions of sex chromosomes for males) or those that are obligatory (e.g., Y snps in women). This is what we want, as we will deal with hetero problems later.


print_text("check what is going on with the samples in the missing file whose potential number of total SNPs is lower", header=3)
print_text("see the distinct values for the number of potential SNPs across the whole sample missing report", header=4)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    awk \
        'BEGIN{FS=\" \"}{if(!a[$5]++){print $5}}' \
        plink.imiss")
            #this is the way to get uniq cases from a colum in awk
                #if the value in $5 has not been previously included in "a" before get True (without the "!" you would get false). Then, print the value of $5 in those cases, i.e., the unique value of potential number of SNPs
                #https://stackoverflow.com/questions/1915636/is-there-a-way-to-uniq-by-column


print_text("process the sample missing file with awk in order to read it with pandas", header=4)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6}' \
        plink.imiss > plink_awk_processed.imiss")
            #I have problems with reading the missing report directly in pandas due to the delimiter. I have to do this to have all 6 columns correctly detected.


print_text("load the missing sample report and the fam file into python", header=4)
sample_missing_report_before_filtering = pd.read_csv( \
    "./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/plink_awk_processed.imiss", \
    sep="\t", \
    header=0, \
    low_memory=False)
merged_fam_file = pd.read_csv( \
    "./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.fam", \
    sep="\t", \
    header=None, \
    low_memory=False)

print_text("check we have the correct number of samples in each batch", header=4)
print(sum(merged_fam_file[0] == "combat_ILGSA24-17303") == 216)
print(sum(merged_fam_file[0] == "combat_ILGSA24-17873") == 1242)


print_text("calculate the number of variants excluding those in Y chromosome", header=4)
n_variants_excluding_y = bim_file_before_any_filter.shape[0]-sum(bim_file_before_any_filter[0]==24)
print(n_variants_excluding_y)
    #in plink
        #Unless you specify otherwise, PLINK interprets chromosome codes as if you were working with human data: 1–22 (or “chr1”–“chr22”) refer to the respective autosomes, 23 refers to the X chromosome, 24 refers to the Y chromosome, 25 refers to pseudoautosomal regions, and 26 refers to mitochondria.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #https://www.cog-genomics.org/plink/1.9/formats#bim


print_text("check that the number of unique cases for n_snps in sample missing is only 2, being one the total number of snps and the second the number of snps minus Y snps. Also check that the samples in the missing report with less than 654027 are non-males, because they should not have SNPs in Y chromosome (chromosome 24)", header=4)
unique_n_variants_missing_rep = sample_missing_report_before_filtering["N_GENO"].unique()
check_1 = len(unique_n_variants_missing_rep) == 2
    #only two different number of total variants in the missing report
check_2 = bim_file_before_any_filter.shape[0] in unique_n_variants_missing_rep
    #the total number of snps in the whole snp map is one of these two numbers
check_3 = unique_n_variants_missing_rep[unique_n_variants_missing_rep !=bim_file_before_any_filter.shape[0]] == n_variants_excluding_y
    #the other number is the number of non-Y variants
check_4 = \
    sum( \
        sample_missing_report_before_filtering \
            .loc[ \
                sample_missing_report_before_filtering["N_GENO"] == n_variants_excluding_y, \
                "IID"] \
            .isin( \
                merged_fam_file \
                    .loc[ \
                        merged_fam_file[4].isin([0,2]), \
                        1])) == \
    sum(sample_missing_report_before_filtering["N_GENO"] == n_variants_excluding_y)
            #get the IDs of all samples in the missing report whose number of variants is 654027 - number of Y SNPs. 
            #these IDs should be all included in those IDs of non-males (0 and 2 which is unknown and female) in the fam file
            #the cases with True for this should be the same than ALL samples in the report that have the number of non-Y variants as the potential number of genotypes 
if(check_1 & check_2 & check_3 & check_4):
    print("We do have the expected number of potential genotypes for samples in the missing report according to the number of non-Y variants in non-male samples. Y SNPs in females are obligatory missing, so they are not counted.")
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE DIFFERENT NUMBER OF VARIANTS BETWEEN MEN AND WOMEN")



print_text("check differences in call rate and MAF between batches. According to Ritchie's tutorial, one simple approach to detect batch effects is to calculate the average minor allele frequency and average genotyping call rate across all SNPs for each batch. Gross differences in either of these on any batch can easily be identified", header=3)
print_text("call rate", header=4)
batch_1_call_rate = 1 - sample_missing_report_before_filtering \
    .loc[ \
        sample_missing_report_before_filtering["FID"]=="combat_ILGSA24-17303", \
        "F_MISS"] \
    .mean()
batch_2_call_rate = 1 - sample_missing_report_before_filtering \
    .loc[ \
        sample_missing_report_before_filtering["FID"]=="combat_ILGSA24-17873", \
        "F_MISS"] \
    .mean()
        #calculate the average frequency of missing (number of non-obligatory missing / number potential genotypes) per sample in each batch and then subtract from 1 to get the call rate.
        #The results are going to be over 1, so 0.98, 0.99....
print(f"average call rate for batch 1 and 2 is {round(batch_1_call_rate, 3)} and {round(batch_2_call_rate, 3)}, respectively")
import numpy as np
diff_call = np.abs(batch_1_call_rate-batch_2_call_rate)
if (diff_call <= 0.01):
    print(f"The difference in call rate between batches is {round(diff_call, 3)}%")
else:
    raise ValueError("ERROR: FALSE! WE HAVE MORE THAN 1% OF DIFFERENCE IN CALL RATE BETWEEN BATCHES")
    

print_text("MAF", header=4)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    plink \
        --bfile ./merged_batches \
        --family \
        --freq; \
    ls -l")
        #calculate MAF per SNP but stratifying by batch
            #By itself, --freq writes a minor allele frequency report to plink.frq. 
            #you can use --freq with --within/--family to write a cluster-stratified frequency report to plink.frq.strat
                #https://www.cog-genomics.org/plink/1.9/basic_stats
                #https://www.cog-genomics.org/plink/1.9/input#within
#process the file with awk in order to correctly specify the delimiter and avoid problems with pandas
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6,$7,$8}' \
        plink.frq.strat > plink_awk_processed.frq.strat; \
    head plink_awk_processed.frq.strat")
#read the stratified frequency file
freqency_report_before_filtering = pd.read_csv( \
    "./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/plink_awk_processed.frq.strat", \
    sep="\t", \
    header=0, \
    low_memory=False)
        #A text file with a header line, and then C lines per variant (where C is the number of clusters; batches in our case) with the following 8-9 lines:
            #CHR Chromosome code
            #SNP Variant identifier
            #CLST    Cluster identifier
            #A1  Allele 1 (usually minor)
            #A2  Allele 2 (usually major)
            #MAF Allele 1 frequency in cluster
            #MAC Allele 1 count in cluster
            #NCHROBS Number of allele observations in cluster
                #https://www.cog-genomics.org/plink/1.9/formats#frq_strat
#calculate the average MAF within each batch 
avg_maf_batch_1 = freqency_report_before_filtering \
    .loc[ \
        freqency_report_before_filtering["CLST"] == "combat_ILGSA24-17303", \
        "MAF"] \
    .mean()
avg_maf_batch_2 = freqency_report_before_filtering \
    .loc[ \
        freqency_report_before_filtering["CLST"] == "combat_ILGSA24-17873", \
        "MAF"] \
    .mean()
print(f"average MAF for batch 1 and 2 is {round(avg_maf_batch_1, 3)} and {round(avg_maf_batch_2, 3)}, respectively")
diff_maf = np.abs(avg_maf_batch_1-avg_maf_batch_2)
if (diff_maf <= 0.01):
    print(f"The difference in call rate between batches is {round(diff_maf, 3)}%")
else:
    raise ValueError("ERROR: FALSE! WE HAVE MORE THAN 1% OF DIFFERENCE IN MAF BETWEEN BATCHES")



print_text("check we have the correct number of samples looking at the merged fam file", header=3)
run_bash(" \
    n_samples_batch=$( \
        awk \
            -F '\t' \
            'END{print NR}' \
            ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.fam); \
    if [[ $n_samples_batch -eq 1458 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #decompress the fam file of each merged file, i.e., for each batch, then count the number of lines (load to awk indicating sep of fields and count number of rows at the end) and save the output
        #get the number of the row at the end, i.e., the total number of rows
            #use awk instead wc -l because "wc -l" counts "new line" character, and it is possible to have a line without that character. If at the end of the file, you do not have an empty line, your actual last line will not have "new line" character so it will not be counted by wc -l.
                #https://unix.stackexchange.com/a/362280
                #https://stackoverflow.com/a/28038682/12772630
        #if the number of lines (samples) in the first and the second batch is 216 and 1242, respectively, then OK
            #remember that we lose 6 duplicated samples in the second batch.
            #note that "==" is only for strings, you have to use "-eq" for numbers
                #https://stackoverflow.com/questions/20449543/shell-equality-operators-eq



print_text("check that we do NOT have the duplicate samples (1100JHJM, 1200JPJM and 7800AGSO) in the second batch", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    n_no_dups=$( \
        awk \
            -F '\t' \
            '{ \
                if($2!=\"1100JHJM_1\" && $2!=\"1200JPJM_1\" && $2!=\"7800AGSO_1\" && $2!=\"1100JHJM_2\" && $2!=\"1200JPJM_2\" && $2!=\"7800AGSO_2\"){\
                    count ++\
                } \
            } \
            END {print count}' \
            merged_batches.fam); \
    total_n=$( \
        awk \
            -F '\t' \
            'END{print NR}' \
            merged_batches.fam); \
    if [[ $n_no_dups -eq $total_n ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #first create a count in awk adding 1 for each row in the fam file that do not contain the ID of the duplicated samples (both suffixes), then print the count at the end.
            #https://stackoverflow.com/questions/14739057/using-awk-with-column-value-conditions
            #https://stackoverflow.com/questions/15839723/awk-or-statement
            #https://stackoverflow.com/questions/12809909/efficient-way-to-count-the-amount-lines-obeying-some-condition
        #second calculate the number of rows of the whole fam file
        #if both numbers are the same, then number of rows without the duplicated samples is exactly the same than the total number of rows, i.e., we do NOT have duplicated samples.

# endregion




##############################################
# region REMOVE DUPLICATED SNPS ##############
##############################################
print_text("Remove SNPs duplicated. We are going to remove first duplicates because this should not be affected by population structure", header=2)
#Remember that plink considers duplicates by POSITION and ALLELES CODES. By default, this ignores A1/A2 allele assignments, since PLINK 1 normally does not preserve them. Therefore, two variants with identical positions and reversed allele assignments are considered duplicates. In our case, I am not sure what information uses Illumina to set the Allele 1 and 2 in the forward strand, but likely it is not REF/ANCESTRAL and plink 1 does not preserve A1/A2 allele assignments anyway. Therefore, we should not use this information and just consider as duplicates two SNPs with the same position and allele codes irrespectively of the allele1/allele2 assignment.
#if two SNPs have the positions and allele codes, this is going to be irrespectively from the ancestry, this is in the SNP map of the whole batch, so all samples share the same information for these two SNPs in the map, except the genotype values, of course, and what allele is minor/major, but as I said, we are only using the position of the SNP the and the names of the alleles, not their assignment, so we should be fine.
#Also, we do not want to use duplicated/error data for calculating PCAs and doing other checks, so it is better to just remove it now.
    #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars



print_text("check we have 654027 SNPs, saving the number of SNPs into a variable", header=3)
n_snps = run_bash(" \
    awk \
        -F '\t' \
        'END{print NR}' \
        ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.bim", \
    return_value=True).strip()
n_snps = int(n_snps)
print(n_snps == 654027)



print_text("check SNPs IDs are not duplicated, as plink 1.9 does not look for duplicates using the ID", header=3)
run_bash(" \
    n_snp_ids_non_dup=$( \
        awk \
            -F '\t' \
            '!a[$2]++{count++} \
            END{print count}' \
            ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.bim); \
        if [[ $n_snp_ids_non_dup -eq " + str(n_snps) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #detect duplicates of ID
                #a[$2] is the name of the array that holds $2 (SNP ID) as keys.
                #uses the current value of $2 as key to the array "a", taking the value stored there. If this particular key was never referenced before, a[$2] evaluates to the empty string.
                #In "!a[$2]", the ! negates the value from before. If it was empty or zero (i.e., the key was never referenced before; false), we now have a true result. If it was non-zero (i.e., the key was referenced before; true in the original test), we have a false result. If the whole expression evaluated to true, meaning that a[$2] was not set to begin with, the whole line is printed as the default action. Therefore, it is only printing a value of $2 if was never references before, i.e., non-duplicate.
                #In other words:
                    #a[$2]: look at the value of key $2, in associative array a. If it does not exist, automatically create it with an empty string.
                    #!a[$2]++: negate the value of expression. If a[$2]++ returned 0 (a false value), the whole expression evaluates to true, and makes awk perform the default action print $2. Otherwise, if the whole expression evaluates to false, no further action is taken.
            #a[$1]++: 
                #increment the value of a[$2], return the old value as value of expression. The ++ operator returns a numeric value, so if a[$2] was empty to begin with, 0 is returned and a[$2] incremented to 1.
            #save into count
            #the total number should be equal to the number of SNPs, meaning we do not have any duplicated SNP.
            #https://unix.stackexchange.com/a/159697
            #https://stackoverflow.com/a/32085039/12772630



print_text("check again SNPs IDs are not duplicated using pandas", header=3)
sum( \
    pd.read_csv( \
        "./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.bim", \
        sep="\t", \
        header=None, \
        low_memory=False) \
    .iloc[:,1] \
    .duplicated( \
        keep=False)) == 0
        #load the bim file
        #get all rows of the second column (IDs)
        #check for duplicates
            #keep=False:
                #Mark all duplicates as ``True``, not only the second occurrence.
        #the sum of Trues should be zero



print_text("list duplicated SNPs", header=3)
run_bash(" \
    cd ./data/genetic_data/; \
    mkdir \
        --parents \
        ./quality_control/00_snp_duplicates/; \
    plink \
        --bfile ./plink_bed_files/merged_batches/merged_plink_files/merged_batches \
        --list-duplicate-vars suppress-first ids-only\
        --out  ./quality_control/00_snp_duplicates/list_duplicates; \
    ls -l ./quality_control/00_snp_duplicates/")
        #list-duplicate-vars to list duplicates by POSITION and ALLELES CODES. By default, this ignores A1/A2 allele assignments, since PLINK 1 normally does not preserve them. Therefore, two variants with identical positions and reversed allele assignments are considered duplicates. To avoid this, you have to use "use the 'require-same-ref' modifier, along with --keep-allele-order/--a2-allele". In our case, I am not sure what information uses Illumina to set the Allele 1 and 2 in the forward strand, but likely it is not REF/ANCESTRAL and plink 1 does not preserve A1/A2 allele assignments anyway. Therefore, we should not use this information and just consider as duplicates two SNPs with the same position and allele codes irrespectively of the allele1/allele2 assignment.
            #suppress-first prevents the first variant in each group from being reported (since, if you're removing duplicates, you probably want to keep one member of each group).
            #ids-only modifier removes the header and the position/allele columns, so the generated list can be used as input for --exclude
                #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
        #in plink 2 you can also look for duplicates in ID, but we already check SNP names duplicates with spark and it is seems is ok



print_text("load the files with duplicates into pandas", header=3)
duplicate_cases = pd.read_csv(
    "./data/genetic_data/quality_control/00_snp_duplicates/list_duplicates.dupvar", 
    sep="\t", 
    header=None,
    low_memory=False)
print(duplicate_cases)



print_text("the file generated has in each row ALL snps that have the same position and alleles, except the first one, thus, the number of rows it is not the total number of duplicates. Therefore, we need also to know the number of snps in each row.", header=3)
print("For example: GSA-rs142775857 rs184047537 rs200286606")
#row="GSA-rs142775857 rs184047537 rs200286606"
n_duplicates_plink = sum([len(row.split(" ")) if " " in row else 1 for row in duplicate_cases[0]])
    #select each row in the first and only column of duplicate cases (only one column due to "ids-only" flag in plink)
    #split the row by space and count the number of pieces, i.e., snps, if there are spaces in the row
    #if not, just count 1
    #then sum
    #https://stackoverflow.com/questions/4406389/if-else-in-a-list-comprehension

print_text("see number of plink duplicates", header=3)
print(n_duplicates_plink)
print(f"The number of duplicates is {(n_duplicates_plink/n_snps)*100} percent respect to the total number of SNPs")



print_text("calculate the number of duplicates with pandas and check is the same than plink duplicates. As we want to consider as duplicates SNPs with the same alleles, irrespectively of the assigment, we need to create a new column where we get both alleles in the same order independently if allele 1 is A or T. Therefore A1=A and A2=T will be the same than A1=T and A2=A, i.e., AT.", header=3)
bim_file_dup_check = pd.read_csv( \
        "./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/merged_batches.bim", \
        sep="\t", \
        header=None, \
        low_memory=False)
from natsort import natsorted
#row=bim_file_dup_check.iloc[0,:]
def alleles_no_assigment(row):
    return "".join( \
        natsorted( \
            [str(row[4]), str(row[5])]))
        #make a list with the two alleles as strings
        #then apply natural order, so we have the same order independently of the order between columns
        #join without space
            #https://stackoverflow.com/a/12453584/12772630
#apply the function per row so you can access the alleles of each SNP
bim_file_dup_check["new_alleles"] = bim_file_dup_check \
    .apply(alleles_no_assigment, axis=1)
#check for duplicates considering the chromosome (0), physical position (3) and new variable about alleles
    #https://www.cog-genomics.org/plink/1.9/formats#bim
bool_dups_pandas = bim_file_dup_check \
    .duplicated( \
        subset=[0, 3, "new_alleles"], \
        keep="first")
            #subset: Only consider certain columns for identifying duplicates
                #https://stackoverflow.com/a/46640992/12772630
            #keep="first": Mark duplicates as ``True`` except for the first occurrence.
                #this is what we did in plink to get the list of duplicates.
#count True cases for dups in pandas
n_duplicates_pandas = sum(bool_dups_pandas)
print(n_duplicates_plink == n_duplicates_pandas)



print_text("if there are less than 2% of duplicates, cool, else stop", header=3)
#duplicated positions should be merged or removed. In our case, we are talking about 1% of the snps, so it should not be a problem to remove them.
if (n_duplicates_plink/n_snps)*100 < 2:
    print("less than 2% of duplicates, good to go!")
else:
    raise ValueError("ERROR: FALSE! WE HAVE MORE THAN 2% OF SNPS WITH DUPLICATED POSITION")



print_text("filter these snps", header=3)
run_bash(
    "cd ./data/genetic_data/; \
    plink \
        --bfile ./plink_bed_files/merged_batches/merged_plink_files/merged_batches \
        --exclude ./quality_control/00_snp_duplicates/list_duplicates.dupvar \
        --make-bed \
        --out ./quality_control/00_snp_duplicates/merged_batches_no_snp_duplicates; \
    ls -l ./quality_control/00_snp_duplicates/")
        #--bfile: 
            #This flag causes the binary fileset plink.bed + plink.bim + plink.fam to be referenced. If a prefix is given, it replaces all instances of 'plink', i.e., it looks for bed, bim and fam files having that suffix instead of "plink".
            #https://www.cog-genomics.org/plink/1.9/input#bed
        #--exclude: Exclude ALL variants named in the file.
            #We have in the ".dupvar" file only the second and next occurrences of each position/allele duplicate, not including the first one. Therefore, we can remove these SNPs.
            #https://www.cog-genomics.org/plink/1.9/filter#snp
        #--make-bed creates a new PLINK 1 binary fileset, AFTER applying sample/variant filters and other operations
            #https://www.cog-genomics.org/plink/1.9/data#make_bed



print_text("remove the non-compressed fileset of the input folder")
run_bash(
    "cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/; \
    rm ./merged_batches.bed; \
    rm ./merged_batches.bim; \
    rm ./merged_batches.fam; \
    ls -l")


print_text("check again for duplicates", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/00_snp_duplicates/; \
    mkdir \
        --parents \
        ./check_no_duplicates/; \
    plink \
        --bfile ./merged_batches_no_snp_duplicates \
        --list-duplicate-vars suppress-first ids-only\
        --out  ./check_no_duplicates/check_duplicates; \
    ls -l ./check_no_duplicates")
        #--list-duplicate-vars
            #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
            #This is not based on variant IDs. Use PLINK 2.0's --rm-dup for ID-based deduplication.
            #By default, this ignores A1/A2 allele assignments, since PLINK 1 normally does not preserve them. If you want two variants with identical positions and reversed allele assignments to not be considered duplicates, use the 'require-same-ref' modifier, along with --keep-allele-order/--a2-allele.
            #Normally, the report has a header line, and contains positions and allele codes in the first 3-4 columns. However, if you just want an input file for --extract/--exclude, the 'ids-only' modifier removes the header and the position/allele columns, and 'suppress-first' prevents the first variant in each group from being reported (since, if you're removing duplicates, you probably want to keep one member of each group).
            #--list-duplicate-vars fails in 'ids-only' mode if any of the reported variant IDs are not unique.


print_text("Do we have zero duplicated positions after filtering?", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/00_snp_duplicates/check_no_duplicates/; \
    n_lines=$( \
        awk \
            -F '\t' \
            'END{print NR}' \
            check_duplicates.dupvar); \
    if [[ $n_lines -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")



print_text("load the bim file after removing duplicates", header=3)
bim_no_dups = pd.read_csv( \
    "./data/genetic_data/quality_control/00_snp_duplicates/merged_batches_no_snp_duplicates.bim", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(bim_no_dups)



print_text("check that the IDs in the original bim file after excluding pandas dups are the same than the remaining IDs in the new bim file obtained with plink", header=3)
print(bim_file_dup_check \
    .loc[~bool_dups_pandas,1] \
    .reset_index(drop=True) \
    .equals( \
        bim_no_dups.loc[:,1]))
        #exclude duplicated rows according to pandas and select the ID column
            #"bool_dups_pandas" contains booleans, with True for second and next occurrences of each group of duplicates
        #reset the index to have the same than in the new cleaned bim file created with plink
        #check that this is identical than the IDs in the new cleaned bim file



print_text("check that we do not have any NA in the SNP IDs in the cleaned bim file", header=3)
if sum(bim_no_dups.loc[:, 1].isna()) == 0:
    print("We do NOT have NA in SNP IDs, so good to go!")
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE NA IN SNP IDS")

# endregion




##############################################
# region FILTER SNPS BY CHROMOSOME TYPE ######
##############################################
print_text("filter SNPs by chromosome type", header=2)
print_text("create folder to do operations", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./01_chromosome_filtering/; \
    ls -l")



print_text("see the unique chromosome names in the SNP map", header=3)
print_text("select the name of the zip file with batch data based on the batch name because ILGSA24-17873 is called CAGRF20093767.zip, not ILGSA24-17873.zip", header=4)
#batch_name = "ILGSA24-17873"
for batch_name in ["ILGSA24-17873", "ILGSA24-17303"]:
    print_text(batch_name + ": starting", header=4)
    if batch_name=="ILGSA24-17873":
        zip_name = "CAGRF20093767"
    elif batch_name=="ILGSA24-17303":
        zip_name = "ILGSA24-17303"
    print("zip name is: " + zip_name)

    print_text(batch_name + ": load zip info from the zip file", header=4)
    import zipfile
    zipdata = zipfile.ZipFile("./data/genetic_data/illumina_batches/" + zip_name + ".zip")
    zipinfos = zipdata.infolist()

    print_text(batch_name + ": extract the zipinfo of the SNP_map", header=4)
    import numpy as np
    #zipinfo=zipinfos[0]
    for zipinfo in zipinfos:
        if zipinfo.filename == zip_name + "/SNP_Map.txt":
            zipinfo.filename = zipinfo.filename.split(zip_name+"/")[1]
            zipinfo_snp_map = zipinfo
                #select the zipinfo of SNP after removing the first part with the parent folder name
    print(zipinfo_snp_map)

    print_text(batch_name + ": extract the SNP map", header=4)
    zipdata.extract(zipinfo_snp_map, "./data/genetic_data/quality_control/01_chromosome_filtering/")

    print_text(batch_name + ": add batch name to the SNP map file", header=4)
    run_bash(" \
        cd ./data/genetic_data/quality_control/01_chromosome_filtering/; \
        mv SNP_Map.txt " + batch_name + "_SNP_Map.txt; \
        ls -l")

print_text("check both SNP maps identical", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/01_chromosome_filtering/; \
    head ILGSA24-17303_SNP_Map.txt; \
    head ILGSA24-17873_SNP_Map.txt; \
    map_status=$(cmp --silent ILGSA24-17303_SNP_Map.txt ILGSA24-17873_SNP_Map.txt; echo $?); \
    if [[ $map_status -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")

print_text("remove the map of one batch are retain just one to do heck as they are identical", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/01_chromosome_filtering/; \
    rm ILGSA24-17303_SNP_Map.txt; \
    mv ILGSA24-17873_SNP_Map.txt SNP_Map.txt; \
    ls -l")

print_text("extract the unique chromosome names using awk", header=4)
unique_chr_map = run_bash(" \
    cd ./data/genetic_data/quality_control/01_chromosome_filtering/; \
    tail \
        -n +2 \
        SNP_Map.txt | \
    awk \
        -F '\t' \
        '{if(!a[$3]++){print $3}}'", return_value=True)
        #remove the column names with tail (first row)
            #if you use tail with "+number_line" you can get all lines starting from the selected number of line
                #tail is faster than sed
                #https://stackoverflow.com/a/339941/12772630
        #see unique cases of column 3
            #a[$3] is the name of the array that holds $3 (chr name) as keys.
                #uses the current value of $3 as key to the array "a", taking the value stored there. If this particular key was never referenced before, a[$3] evaluates to the empty string.
            #In "!a[$3]", the ! negates the value from before. If it was empty or zero (i.e., the key was never referenced before; false), we now have a true result. If it was non-zero (i.e., the key was referenced before; true in the original test), we have a false result. If the whole expression evaluated to true, meaning that a[$3] was not set to begin with, the whole line is printed as the default action. Therefore, it is only printing a value of $3 if was never references before, i.e., non-duplicate.
            #In other words:
                #a[$3]: look at the value of key $3, in associative array a. If it does not exist, automatically create it with an empty string.
                #!a[$3]++: negate the value of expression. If a[$3]++ returned 0 (a false value), the whole expression evaluates to true, and makes awk perform the default action print $3. Otherwise, if the whole expression evaluates to false, no further action is taken.
            #a[$3]++: 
                #increment the value of a[$3], return the old value as value of expression. The ++ operator returns a numeric value, so if a[$3] was empty to begin with, 0 is returned and a[$3] incremented to 1.
            #other explanation calling "a" as "seen"
                #seen is the arbitrary name of an associative array. It is not an option of any kind. You could use a or b or most other names in its place.
                #The code !seen[$0]++ consists of a TEST and an INCREMENT.
                #If seen[$0], i.e. the value of the array element associated with the key $0, the current line of input, is zero (or empty), then the test is true.
                    #I guess this means if the value of the current line is NOT present as a key in seen, then true, because it has not been seen before
                #The value in the array corresponding to the key $0 is then incremented, which means that the test will be false all other times that the same value of $0 is found.
                    #This value is now added as a key
                #The effect is that the test is true the first time a particular line is seen in the input, and false all other times.
                #Whenever the a test with no associated action is true, the default action is triggered. The default action is the equivalent of { print } or { print $0 }, which prints the current record, which for all accounts and purposes in this example is the current unmodified line of input.
                #https://unix.stackexchange.com/a/604294

print_text("split the string generated with the awk unique chromsome names using '\n', but avoid the last element which is the last empty line, then order using natural sorting", header=4)
unique_chr_map = natsorted(unique_chr_map.split("\n")[0:-1])
print(unique_chr_map)

print_text("create a list with the expected chromosome names, i.e., unknown (0), 22 autosomals, mito, sex and pseudo-autosomal", header=4)
expected_chr = [str(i) for i in range(0, 23, 1)]
[expected_chr.append(i) for i in ["MT", "X", "XY", "Y"]]
    #first create list with unknown and autosomes (range from 0 to 23, without including 23)
    #then add mito, sex and PAR to that list

print_text("check we have the expect unique chromosome names", header=4)
print(unique_chr_map == expected_chr)

print_text("check again we have the expect unique chromosome names using pandas this time", header=4)
print( \
    natsorted( \
        pd.read_csv( \
            "./data/genetic_data/quality_control/01_chromosome_filtering/SNP_Map.txt", \
            sep="\t", \
            header=0, \
            low_memory=False) \
        ["Chromosome"] \
        .unique()) == \
    expected_chr)
        #load SNP_Map as pandas DF
        #get the chromosome column
        #select unique cases
        #apply natural sorting
        #check is the same than expected unique chromosome names


print_text("we will use an empty list to save the IDs of the SNPs with strange chromosomes, then use plink to exclude", header=3)
list_snp_ids_to_exclude = []
print(list_snp_ids_to_exclude)
    #https://www.biostars.org/p/281334/


print_text("explore chromosome 0", header=3)
#chromosome with 0
    #in the SNP_Map file of illumina
        #During manifest creation, the probe sequence is mapped to the genome build notated in the GenomeBuild column in the manifest. The single nucleotide polymorphisms (SNP) or indel coordinate is then recorded in the MapInfo column. If a marker is annotated with Chr = 0 or MapInfo = 0 in the manifest, there are two possible explanations:
            #No valid mapping for the probe.
            #More than 1 best-scoring mapping for the probe.
        #https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001423
    #in plink
        #Unless you specify otherwise, PLINK interprets chromosome codes as if you were working with human data: 1–22 (or “chr1”–“chr22”) refer to the respective autosomes, 23 refers to the X chromosome, 24 refers to the Y chromosome, 25 refers to pseudoautosomal regions, and 26 refers to mitochondria.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #https://www.cog-genomics.org/plink/1.9/formats#bim
    #I understand then that the SNP with chromosome=0 is not correctly mapped to hg38 (our genome build), so we do not know where it is located. Therefore we cannot use SNPs with zero
    #these are only 12 cases anyway

print_text("check we have 0 chromosome in the SNP_map", header=4)
print("0" in unique_chr_map)

print_text("check 0 in plink is 0 in illumina snp_map", header=4)
snp_map = pd.read_csv( \
    "./data/genetic_data/quality_control/01_chromosome_filtering/SNP_Map.txt", \
    sep="\t", \
    header=0, \
    low_memory=False)
print( \
    natsorted( \
        snp_map \
            .loc[snp_map["Chromosome"] == "0", "Name"] \
            .to_list()) == \
    natsorted( \
        bim_file_dup_check \
            .loc[bim_file_dup_check[0]==0, 1] \
            .to_list()))
                #get the IDs of SNPs that have as chromosome XY and 25 in the SNP_map and plink (original bim file without any filter), respectively.
                    #I used bim original because we are using the SNP map, which is not filtered by duplicates, so bim_no_dups would have a different (reduced) set of SNPs
                #convert to list
                #apply natural sorting
                #check identical

print_text("count SNPs for which we have 0 chromosome", header=4)
chr_0_bool = bim_no_dups[0]==0
print(f"We have {sum(chr_0_bool)} SNPs with chromosome zero")
if sum(chr_0_bool) < 50:
    print("we have less than 50 variants with chromosome zero, no problem! add them to the list to SNPs to be excluded")
    chrom_zero_cases = bim_no_dups.loc[chr_0_bool, 1].to_list()
        #select ID of these cases and convert to list
    [list_snp_ids_to_exclude.append(snp) for snp in chrom_zero_cases]
        #take each of these snps and save it in the result list
    print(list_snp_ids_to_exclude)
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE TOO MUCH VARIANTS WITH CHROMOSOME 0")


print_text("explore SNPs with XY", header=3)
#I understand these are SNPs in pseudoautosomal regions of the sex chromsomes, i.e., regions of the X and Y chromosome that are homologues and recombinate during meiosis. Therefore, they behave like autosomes in these regions.
    #in illumina
        #I have snps whose chromosome is XY
    #in plink
        #Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name.
        #Unless you specify otherwise, PLINK interprets chromosome codes as if you were working with human data: 1–22 (or “chr1”–“chr22”) refer to the respective autosomes, 23 refers to the X chromosome, 24 refers to the Y chromosome, 25 refers to pseudoautosomal regions, and 26 refers to mitochondria.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #https://www.cog-genomics.org/plink/1.9/formats#bim

print_text("check we have XY chromosome in the SNP_map", header=4)
print("XY" in unique_chr_map)

print_text("check 25 in plink is XY in illumina snp_map", header=4)
print( \
    natsorted( \
        snp_map \
            .loc[snp_map["Chromosome"] == "XY", "Name"] \
            .to_list()) == \
    natsorted( \
        bim_file_dup_check \
            .loc[bim_file_dup_check[0]==25, 1] \
            .to_list()))
                #get the IDs of SNPs that have as chromosome XY and 25 in the SNP_map and plink, respectively.
                    #I used bim original because we are using the SNP map, which is not filtered by duplicates, so bim_no_dups would have a different (reduced) set of SNPs
                #convert to list
                #apply natural sorting
                #check identical

print_text("count SNPs for which we have XY chromosome, i.e., pseudo-autosomal", header=4)
chr_XY_bool = bim_no_dups[0]==25
print(f"We have {sum(chr_XY_bool)} SNPs pseudo-autosomal regions")
if sum(chr_XY_bool) < 1000:
    print("we have less than 1000 variants in pseudo-autosomal regions, no problem!")
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE MORE THAN 1000 VARIANTS IN PSEUDO-AUTOSOMAL REGIONS")


print_text("explore SNPs in MT", header=3)
#I understand these are SNPs in mitochondrial genome
    #in illumina
        #I have snps whose chromosome is MT
    #in plink
        #Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name.
        #Unless you specify otherwise, PLINK interprets chromosome codes as if you were working with human data: 1–22 (or “chr1”–“chr22”) refer to the respective autosomes, 23 refers to the X chromosome, 24 refers to the Y chromosome, 25 refers to pseudoautosomal regions, and 26 refers to mitochondria.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #https://www.cog-genomics.org/plink/1.9/formats#bim

print_text("check we have MT chromosome in the SNP_map", header=4)
print("MT" in unique_chr_map)

print_text("check 26 in plink is MT in illumina snp_map", header=4)
print( \
    natsorted( \
        snp_map \
            .loc[snp_map["Chromosome"] == "MT", "Name"] \
            .to_list()) == \
    natsorted( \
        bim_file_dup_check \
            .loc[bim_file_dup_check[0]==26, 1] \
            .to_list()))
                #get the IDs of SNPs that have as chromosome XY and 25 in the SNP_map and plink, respectively.
                    #I used bim original because we are using the SNP map, which is not filtered by duplicates, so bim_no_dups would have a different (reduced) set of SNPs
                #convert to list
                #apply natural sorting
                #check identical

print_text("count SNPs for which we have MT chromosome", header=4)
chr_mito_bool = bim_no_dups[0]==26
print(f"We have {sum(chr_mito_bool)} SNPs in mito genome")
if sum(chr_mito_bool) < 1500:
    print("we have less than 1500 variants in mito genome, no problem!")
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE MORE THAN 1500 VARIANTS IN mito genome")


print_text("explore SNPs in X", header=3)
#I understand these are SNPs in X
    #in illumina
        #I have snps whose chromosome is X
    #in plink
        #Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name.
        #Unless you specify otherwise, PLINK interprets chromosome codes as if you were working with human data: 1–22 (or “chr1”–“chr22”) refer to the respective autosomes, 23 refers to the X chromosome, 24 refers to the Y chromosome, 25 refers to pseudoautosomal regions, and 26 refers to mitochondria.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #https://www.cog-genomics.org/plink/1.9/formats#bim

print_text("check we have X chromosome in the SNP_map", header=4)
print("X" in unique_chr_map)

print_text("check 23 in plink is X in illumina snp_map", header=4)
print( \
    natsorted( \
        snp_map \
            .loc[snp_map["Chromosome"] == "X", "Name"] \
            .to_list()) == \
    natsorted( \
        bim_file_dup_check \
            .loc[bim_file_dup_check[0]==23, 1] \
            .to_list()))
                #get the IDs of SNPs that have as chromosome XY and 25 in the SNP_map and plink, respectively.
                    #I used bim original because we are using the SNP map, which is not filtered by duplicates, so bim_no_dups would have a different (reduced) set of SNPs
                #convert to list
                #apply natural sorting
                #check identical

print_text("count SNPs in X", header=4)
chr_X_bool = bim_no_dups[0]==23
print(f"We have {sum(chr_X_bool)} SNPs in X")
if sum(chr_X_bool) > 20000:
    print("we have more than 20000 variants in X, no problem!")
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE LESS THAN 20000")


print_text("explore SNPs in Y", header=3)
#I understand these are SNPs in Y
    #in illumina
        #I have snps whose chromosome is Y
    #in plink
        #Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name.
        #Unless you specify otherwise, PLINK interprets chromosome codes as if you were working with human data: 1–22 (or “chr1”–“chr22”) refer to the respective autosomes, 23 refers to the X chromosome, 24 refers to the Y chromosome, 25 refers to pseudoautosomal regions, and 26 refers to mitochondria.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #https://www.cog-genomics.org/plink/1.9/formats#bim

print_text("check we have Y chromosome in the SNP_map", header=4)
print("Y" in unique_chr_map)

print_text("check 24 in plink is Y in illumina snp_map", header=4)
print( \
    natsorted( \
        snp_map \
            .loc[snp_map["Chromosome"] == "Y", "Name"] \
            .to_list()) == \
    natsorted( \
        bim_file_dup_check \
            .loc[bim_file_dup_check[0]==24, 1] \
            .to_list()))
                #get the IDs of SNPs that have as chromosome XY and 25 in the SNP_map and plink, respectively.
                    #I used bim original because we are using the SNP map, which is not filtered by duplicates, so bim_no_dups would have a different (reduced) set of SNPs
                #convert to list
                #apply natural sorting
                #check identical

print_text("count SNPs in Y", header=4)
chr_Y_bool = bim_no_dups[0]==24
print(f"We have {sum(chr_Y_bool)} SNPs in Y")
if (sum(chr_Y_bool) > 3000) & (sum(chr_Y_bool) < 10000):
    print("variants in Y are more than 3K and less than 10K, no problem!")
else:
    raise ValueError("ERROR: FALSE! WE VARIANTS IN Y ARE NOT BETWEEN 3 AND 10K")


print_text("remove selected SNPs", header=3)
print_text("save to a file the list of SNPs to be excluded, each SNP in anew line", header=4)
with open("./data/genetic_data/quality_control/01_chromosome_filtering/snps_to_exclude.txt", "w") as f:
    #snp=list_snp_ids_to_exclude[0]
    for snp in list_snp_ids_to_exclude:
        f.write(f"{snp}\n") #print the SNP and add new line
        #https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python-with-newlines
run_bash(" \
    cat ./data/genetic_data/quality_control/01_chromosome_filtering/snps_to_exclude.txt")

print_text("exclude these SNPs using --exclude of plink", header=4)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./00_snp_duplicates/merged_batches_no_snp_duplicates \
        --exclude ./01_chromosome_filtering/snps_to_exclude.txt \
        --make-bed \
        --out  ./01_chromosome_filtering/merged_batches_filtered_chromosomes; \
    ls -l ./01_chromosome_filtering/")
            #--bfile: 
                #This flag causes the binary fileset plink.bed + plink.bim + plink.fam to be referenced. If a prefix is given, it replaces all instances of 'plink', i.e., it looks for bed, bim and fam files having that suffix instead of "plink".
                #https://www.cog-genomics.org/plink/1.9/input#bed
            #--exclude
                #Exclude ALL variants named in the file
                #--exclude normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all listed variants from the current analysis. 
                #Note that this is slightly different from PLINK 1.07's behavior when the main input fileset contains duplicate variant IDs: PLINK 1.9 removes all matches, while PLINK 1.07 just removes one of the matching variants. If your intention is to resolve duplicates, you should now use --bmerge instead of --exclude.
                    #we do not have duplicate SNP IDs, so no problem here.
                #https://www.cog-genomics.org/plink/1.9/filter#snp
                #https://www.biostars.org/p/281334/
            #--make-bed creates a new PLINK 1 binary fileset, AFTER applying sample/variant filters and other operations
                #https://www.cog-genomics.org/plink/1.9/data#make_bed

print_text("check that the excluded SNPs are actually out of the new bim file", header=4)
bim_chrom_filter = pd.read_csv( \
    "./data/genetic_data/quality_control/01_chromosome_filtering/merged_batches_filtered_chromosomes.bim", \
    sep="\t", \
    header=None, \
    low_memory=False)
print( \
    sum(bim_chrom_filter \
        .loc[:,1] \
        .isin(values=list_snp_ids_to_exclude) \
        .to_list()) == 0)
        #get the column with IDs from the bim file (second column)
        #check if any of the IDs is included in the list of exclusions, generating a pandas series with booleans
            #values: set or list-like
                #The sequence of values to test
        #convert to list
        #the sum should be zero, becuase no True should be present

print_text("check that the number of snps removed from the bim file is the same than the number of snps in the exclusion list", header=4)
print(bim_no_dups.shape[0] - bim_chrom_filter.shape[0] == len(list_snp_ids_to_exclude))

print_text("check that snps present in the first bim file but not in the new one are exactly those of the exclusion list", header=4)
print( \
    bim_no_dups \
        .loc[ \
            ~bim_no_dups \
                .loc[:,1] \
                .isin( \
                    bim_chrom_filter \
                        .loc[:,1]), \
            1] \
        .to_list() == \
    list_snp_ids_to_exclude)
    #from the original bim file, select those rows
        #whose SNP ID is NOT present in the SNP ID column of the new bim file
    #convert to list
    #the list should be identical to list_snp_ids_to_exclude

print_text("check all the unwanted chromosomes are removed", header=4)
unique_chr_after_filter = bim_chrom_filter.loc[:,0].unique()
print(0 not in unique_chr_after_filter)
print(np.array_equal(unique_chr_after_filter, [i for i in range(1, 27, 1)]))
    #the array with unique chromosomes should be 1 (0 was removed) to 26, so range 1-27 because in python the last element of the range is not included.

print_text("see the remaining chromosomes", header=4)
print(unique_chr_after_filter)

print_text("check again we do not have duplicated snp ids", header=4)
print( \
    sum( \
        bim_chrom_filter \
        .loc[:,1] \
        .duplicated(keep=False)) == 0)
    #select the column with SNP ids
    #check for duplicates
        #keep=False:
            #Mark all duplicates as ``True``, not only the second occurrence
    #no true should be present, so the sum should be zero

# endregion




##############################################
# region GeneCall and GeneTrain scores #######
##############################################
print_text("GeneCall and GeneTrain scores", header=2)
print_text("We are not going to use GeneCall and GeneTrain scores. See script for further details", header=3)
#The GenCall score (GC) is a confidence measure assigned to each call which can be used to filter poor quality calls, SNPs or samples. Illumina generally recommend that calls with GC ≤ 0.15 represent failed genotypes. Averaged GC scores over all SNPs from a given sample, or across all samples for a given SNP can be used as sample or SNP quality metrics. A more commonly used sample quality metric is the 'no call rate'. For GenCall, genotypes with GC score less than a given threshold (0.15 in our analyses) are declared as missing. The proportion of missing values, or 'no calls' in each sample gives the no call rate; samples with higher rates are deemed less reliable than samples with lower rates. No call rates less than 1% should be expected for good quality samples which have been properly processed (Illumina Technical Support, personal communication).
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3063825/

#I have found an study of animal breeding where they test whether these two scores are correlated with the concordance, i.e., whether they are markers of genotyping errors. They associate, but there are differences between animal species.
    #https://onlinelibrary.wiley.com/doi/full/10.1111/age.13043

#I have also seen these scores as mentioned in papers about the previous steps, i.e., the steps done by AGRF to convert prob intensities to genotypes. For example, to check these scores to see if manual-re clustering is needed
    #see section "Manual re-clustering"
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/

#what AGRF did was to set no call for snps with GC score < 0.15 filter. Then, we can discard samples for which more than 1% of SNPs have GC score < 0.15, i.e., no call. Indeed, in the PDF report of the first batch, they say "Four samples are below the illumina expected 99% SNP call rate". We are going to stick to that.
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3063825/

#A tutorial of admixture filter by 0.25
    #A few technical details are that the genotypes were filtered with a cutoff of 0.25 for the Illumina GenCall score [13] (a quality score generated by the basic genotype calling software).
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #Kermani BG (2006) Artificial intelligence and global normalization methods for genotyping. U.S. Patent No. 7,035,740. Washington, DC: U.S. Patent and Trademark Office

#Remember that prob intensity is used to split samples between the three genotypes so you plot prob intesnity in try to define three clusters of samples. Sometimes the clusters of hetero and one homozygous can be very close making it difficult to separate... T
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/

#I have NOT checked in each final report whether all genotypes per sampe have gene call > 0.15.
    #However, the DNAReport of both batches says that "Low GenCall Score Cutoff" is 0.15. 
    #I have also checked in these reports that the "Call_Rate" column is very similar to divide the column "#Calls" (number of genotypes with GS score > 0.15) by the total number of SNPs.
    #this call rate is what in the PDF of the first batch they say that is too low for 4 samples, being below the 99% call rate expected for illumina.
    #Indeed, I have checked that those samples with call rate < 0.99 in the DNA report of the first batch are indeed those with call rate < 0.99 in the PDF report.
    #Therefore, I understand that AGFR has applied the filter of GC score < 0.15 to set as no call a given genotype and now we can use this to calculate the call rate per sample and per snp in order to filter by missingness.
    
# endregion




######################################################################
# region check no genetic position in map and select only snps #######
######################################################################
print_text("check no genetic position in map and select only snps", header=2)
print_text("check no genetic position in the map, it should be 0 always", header=3)
print(bim_chrom_filter.loc[:,2].unique() == 0)
    #Column 3 (named as 2 by python as is 0-based):
        #Position in morgans or centimorgans (safe to use dummy value of '0')
            #https://www.cog-genomics.org/plink/1.9/formats#bim

print_text("check no SNP has empty ID", header=3)
print(sum(bim_chrom_filter.loc[:,1].isin(["", " ", "."]))==0)
print(sum(bim_chrom_filter.loc[:,1].isna())==0)

print_text("select only SNPs", header=3)
print_text("see unique alleles: we have insertions and deletions (indels) along with 0, which means all samples are major homozigous and no minor allele is present", header=4)
print(bim_chrom_filter.loc[:, 4].unique())
print(bim_chrom_filter.loc[:, 5].unique())
    #Zero indicates missing.  In the old .ped format, a missing call was indicated by "0 0".  A .bim file should only have both allele codes set to '0' when every sample has a missing call at that site.  
    #More commonly, when every sample is homozygous major at a site, PLINK might not have any of knowing what the minor allele is supposed to be (for example, the old .ped/.map format didn't keep track of that), so it provisionally sets the minor allele code to '0' until the dataset is merged with sample(s) which aren't homozygous major there.
        #https://groups.google.com/g/plink2-users/c/JRrOneIuKRQ

print_text("make new dir", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        ./02_remove_non_snps/ \
        --parents; \
    ls -l")

print_text("remove non-SNPs", header=4)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./01_chromosome_filtering/merged_batches_filtered_chromosomes \
        --snps-only just-acgt\
        --make-bed \
        --out ./02_remove_non_snps/merged_batches_remove_non_snps; \
    ls -l ./02_remove_non_snps")
            #--snps-only excludes all variants with one or more multi-character allele codes. With 'just-acgt', variants with single-character allele codes outside of {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', <missing code>} are also excluded.
                #we are going to remove the Insertion and Deletions, leaving only ACTG.
                #missing cases are retained (see below), and we will check them after filtering by MAF, i.e., after removing SNPs with low minor frequencies.
                #https://www.cog-genomics.org/plink/1.9/filter
            #just using --snp-only does not remove INDELS. My guess is that the format of INDELS is not AAGCC, i.e., having several letter in the allele column of the bim file, but instead I or D. We have only one letter, so I guess plinl cannot detect it. When using "just-acgt", these rows are removed because they do not have single ACGT allele names in the allele columns.

print_text("check no indels are present and we only have ACTG + 0", header=4)
bim_remove_non_snps = pd.read_csv(
    "./data/genetic_data/quality_control/02_remove_non_snps/merged_batches_remove_non_snps.bim", \
    sep="\t", \
    header=None, \
    low_memory=False)
print([element in ["A", "C", "T", "G", "0"] for element in bim_remove_non_snps.loc[:,4].unique()])
print([element in ["A", "C", "T", "G", "0"] for element in bim_remove_non_snps.loc[:,5].unique()])
n_indels = bim_chrom_filter.shape[0]-bim_remove_non_snps.shape[0]
if (n_indels/bim_chrom_filter.shape[0])*100 < 5:
    print(f"We have lost {n_indels} rows after removing INDELS")
else:
    raise ValueError("ERROR: FALSE! We have more than 5% of indels!!")

# endregion




################################################################################
# region remove again duplicates but this time only considering position #######
################################################################################
print_text("remove again duplicates but this time only considering position", header=2)
#plink1.9 removes duplicates considering position AND alleles. Therefore, SNPs with the same position but different allele names are not removed. For example, rs387907306 and rs387907306.1 have the exact same position, but different alleles, AG and TG, respectively. Remember that plink1.9 does not consider the order (A1-A2) for duplicates, so AG and GA are considered duplicates. However, in this case, we have AG and TG, A is not present in the second SNP.  
#We are going to retain only the first instance of each duplicate group, removing the rest of instances.
#Some duplications could be caused by strand issues. For example, 0T and 0A. The second SNP can have the opposite strand. This will be solved later in the pipeline when we check strand. The remaining SNP will be checked for correct strand.
    #https://www.biostars.org/p/310290/
print_text("see duplicates by position", header=3)
print(bim_remove_non_snps \
    .loc[ \
        bim_remove_non_snps \
            .duplicated( \
                subset=[0,3], \
                keep=False), :])
                #subset: Only consider certain columns for identifying duplicates
                    #https://stackoverflow.com/a/46640992/12772630
                #keep=False: Mark all duplicates as ``True``


print_text("make new dir", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        ./03_remove_pos_dups/ \
        --parents; \
    ls -l")


print_text("select those SNPs with the same position", header=3)
list_snps_pos_dup = bim_remove_non_snps \
    .loc[ \
        bim_remove_non_snps \
            .duplicated( \
                subset=[0,3], \
                keep="first"), 1].to_list()
                #keep="first": Mark duplicates as ``True`` except for the first occurrence.
print(list_snps_pos_dup)


print_text("check the number of SNPs with this problem", header=3)
if (len(list_snps_pos_dup)/bim_remove_non_snps.shape[0])*100 < 1:
    print(f"We have {len(list_snps_pos_dup)} snps with duplicated position")
else:
    raise ValueError("ERROR: FALSE! WE HAVE MORE THAN 1% OF SNPS WITH THE SAME POSITION")


print_text("save as txt", header=3)
with open("./data/genetic_data/quality_control/03_remove_pos_dups/snps_pos_dup.txt", "w") as f:
    #snp=list_snp_ids_to_exclude[0]
    for snp in list_snps_pos_dup:
        f.write(f"{snp}\n") #print the SNP and add new line
        #https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python-with-newlines
run_bash("head ./data/genetic_data/quality_control/03_remove_pos_dups/snps_pos_dup.txt")


print_text("check we have the correct number of SNPs in the list", header=3)
run_bash(" \
    snps_to_remove=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./data/genetic_data/quality_control/03_remove_pos_dups/snps_pos_dup.txt); \
    if [[ $snps_to_remove -eq " + str(len(list_snps_pos_dup)) + " ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #just count the number of rows at the end of the file with awk and check this is the same than the length of the original python list with the SNPs to be removed


print_text("remove these SNPs", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./02_remove_non_snps/merged_batches_remove_non_snps \
        --exclude ./03_remove_pos_dups/snps_pos_dup.txt \
        --make-bed \
        --out ./03_remove_pos_dups/merged_batches_remove_dup_pos; \
    ls -l ./03_remove_pos_dups/")


print_text("check for pos dups after filtering", header=3)
print_text("load the bim file after filtering", header=4)
bim_remove_non_pos_dup = pd.read_csv(
    "./data/genetic_data/quality_control/03_remove_pos_dups/merged_batches_remove_dup_pos.bim", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(bim_remove_non_pos_dup)
print_text("check for position duplicates", header=4)
print(sum(bim_remove_non_pos_dup.duplicated(subset=[0,3], keep=False)) == 0)


print_text("check also for duplicates in ID", header=3)
print(sum(bim_remove_non_pos_dup.loc[:,1].duplicated(keep=False)) == 0)

# endregion




##############################
# region MAF filtering #######
##############################
print_text("MAF filtering", header=2)
#Typically, SNP-level filtering based on a large amount of missing data and lower variability is performed first. This is followed by sample-level filtering (see step 3 in the succeeding texts), and finally, SNP-level filtering based on possible genotyping errors (see step 4 in the succeeding texts) is performed. The rationale for this is that both sample-level relatedness and substructure (for which we filter in step 3) can influence the Hardy–Weinberg equilibrium (HWE) criterion (step 4) used for filtering SNPs based on genotyping errors. 
    #https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
#It also makes sense to follow this order because we want to remove first low-quality SNPs before removing samples. SNPs with low MAF or low call rate are more likely to have genotyping errors. We should remove them before removing sample based on call rate. In this way, we avoid the removal of samples because they have missing calls in low-quality SNPs.
#why this is important and what threshold use
    #SNPs with a low MAF are rare, therefore power is lacking for detecting SNP‐phenotype associations. These SNPs are also more prone to genotyping errors. The MAF threshold should depend on your sample size, larger samples can use lower MAF thresholds.
    #Ritchie's tutorial says "SNPs with frequency too low to yield reasonable statistical power (e.g., below 1%) may be removed from the analysis to lighten the computational and MULTIPLE TESTING CORRECTION BURDEN. However, in studies with very large sample sizes, it may be beneficial to change the threshold to 0.05%."
    #This makes sense as Alicia's tutorial says "...low minor allele frequency variants are typically not in linkage disequilibrium with common variants and, therefore, add a greater multiple testing burden". These are variants are independent from others, greatly increasing the number of independent tests.
    #Alicia's tutorial also says "GWAS’ mostly refers to genome-wide studies of common variants ... Generally, common variants are those with a minor allele frequency above 10%, although as population sizes grow this threshold can be as low as 1% as researchers typically adhere to a minimum minor allele count; for example, at least 100 individuals who carry at least one copy of the minor allele"
    #In the same vein, Ritchie's tutorial says "However, in studies with very large sample sizes, it may be beneficial to change the threshold to 0.05% (i.e., 0.005). More recently, in larger sample-size studies such as Pan UKBB data analyses, minor allele count (MAC) of 20 is also recommended (W. Zhou et al., 2020)."
    #The PRS tutorial recommends "minor allele frequency (MAF) >1% (MAF >5% if target sample N <1000)".
        #As you have more samples, low MAF SNPs will have a higher minor allele count in absolute terms, as we have more samples.
    #the QC tutorial recommended by the PRS tutorial says "Respectively, for large (N = 100.000) vs. moderate samples (N = 10000), 0.01 and 0.05 are commonly used as MAF threshold."
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
    #In this other tutorial, they apply 0.01 and lose aroudn 200K SNPs, so it is plausible to lose hundreds of thousands of SNPs. They also say that "Here, we remove SNPs for which the MAF is less than 1%. In some instances, particularly small sample settings, a cut point of 5% is applied."
        #https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
    #I have found cases with low sample size and MAF filter of 0.01
        #Augusto's paper says MAF<0.2 pre- and MAF<0.01 post-imputation and they only have 90 samples
        #heritage study used maf<0.01 and they had 500 samples
            #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3098655/
    #Given than all tutorials say that MAF filter depends on sample size and they only recommend 0.01 for relatively large sample sizes, I am going to use 0.05 as filter. Yes, we lose many SNPs, but this is the data we have....
#when applying this step:
    #I wasn't sure whether to remove related samples before or after the MAF filtering:
        #The VO2 max paper seems to remove related samples before and plink tutorial remove related samples before MAF filtering
        #In contrast, the Ritche tutorial consider the filter by MAF, SNP and sample call rate as a package called "Standard Quality Control", being explained before sample relatedness.
        #The most important argument I have found to do MAF before is the following:
            #The methods we are gonna use to remove related samples (KING-robust) in plink2 seems to require decent MAF. See this from plink2 help:
                #"The relationship matrix computed by --make-rel/--make-grm-list/--make-grm-bin can be used to reliably identify close relations within a single population, IF YOUR MAFS ARE DECENT"
            #Therefore, I understand that we should not have very low MAFs if we want to robustly calculate relatedness between our samples.
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./04_remove_low_maf_snps/; \
    plink \
        --bfile ./03_remove_pos_dups/merged_batches_remove_dup_pos \
        --maf 0.05 \
        --make-bed \
        --out ./04_remove_low_maf_snps/merged_batches_remove_low_maf; \
    ls -l ./04_remove_low_maf_snps")
            #create a new folder to save filesets after removing related samples
            #apply the MAF filter to the fileset filtered for logR SD
                #--maf filters out all variants with minor allele frequency below the provided threshold (default 0.01), while --max-maf imposes an upper MAF bound. Similarly, --mac and --max-mac impose lower and upper minor allele count bounds, respectively.
                    #https://www.cog-genomics.org/plink/1.9/filter


print_text("see the number of SNPs removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_snps_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./03_remove_pos_dups/merged_batches_remove_dup_pos.bim); \
    n_snps_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./04_remove_low_maf_snps/merged_batches_remove_low_maf.bim); \
    lost_snps=$(($n_snps_before - $n_snps_after)); \
    lost_snps_percent=$( \
        awk \
            -v l=$lost_snps \
            -v b=$n_snps_before \
            'BEGIN{print (l/b)*100}'); \
    if [[ $lost_snps_percent < 60 ]]; then \
        printf 'The number of SNPs lost due to MAF below 0.05 is: %s' \"$lost_snps\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 60% of SNPs due low-call rate';\
    fi")
        #First calculate the number of SNPs, i.e., rows, in the bim files before and after filtering using for that "++" to count cases.
        #Calculate the difference of snps, i.e., SNPs removed.
        #calculate the percentage respect to the initial number of SNPs
            #use the -v flag to load variables into awk, so we can use the number of SNPs removed and those previously present before the filtering. Use awk to calculate the percentage.
                #https://stackoverflow.com/a/12147154/12772630
                #https://www.unix.com/unix-for-dummies-questions-and-answers/15651-awk-v.html
            #bash cannot do operations with floats
        #if the percentage of lost SNPs is less than 1%, we are fine. So print the number of SNPs using printf, so you can combine string message with a variable, indicated as %s because it is loaded as string.
            #https://phoenixnap.com/kb/bash-printf


print_text("create freq report", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/04_remove_low_maf_snps/; \
    plink \
        --bfile merged_batches_remove_low_maf \
        --freq; \
    head plink.frq")
        #--frq produces a frequency report
            #A text file with a header line, and then one line per variant with the following six fields:
                #CHR Chromosome code
                #SNP Variant identifier
                #A1  Allele 1 (usually minor)
                #A2  Allele 2 (usually major)
                #MAF Allele 1 frequency
                #NCHROBS Number of allele observations
            #https://www.cog-genomics.org/plink/1.9/formats#frq


print_text("check that no SNP has a MAF below the selected threshold", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/04_remove_low_maf_snps/; \
    total_number_snps=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>1){count++}}END{print count}' \
            plink.frq);\
    n_snps_above_threshold=$( \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 >= 0.05)){count++}}END{print count}' \
            plink.frq); \
    if [[ $n_snps_above_threshold -eq $total_number_snps ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #count the number of rows of the frequency report, i.e., total number of SNPs
        #then count the number of SNPs in that file that have a MAF equal or higher than 0.05, which is our selected threshold.
        #the number of snps meeting this condition should be the same than the total number of SNPs, because we have already applied the filter.

# endregion




################################
# region SNP missingness #######
################################
print_text("filter by SNP missingness", header=2)
print_text("create missing report per SNP and sample", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/04_remove_low_maf_snps/; \
    plink \
        --bfile merged_batches_remove_low_maf \
        --missing; \
    ls -l")
        #--missing produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.
            #https://www.cog-genomics.org/plink/1.9/basic_stats#missing


print_text("See the first lines of the lmiss file", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/04_remove_low_maf_snps/; \
    head plink.lmiss")
    #lmiss: A text file with a header line, and K line(s) per variant with the following 5-7 fields (where K is the number of cluster(s) if --within/--family was specified, or 1 if it wasn't):
        #CHR Chromosome code
        #SNP Variant identifier
        #CLST    Cluster identifier. Only present with --within/--family.
        #N_MISS  Number of missing genotype call(s), not counting obligatory missings or het. haploids
        #N_CLST  Cluster size (does not include nonmales on chrY). Only present with --within/--family.
        #N_GENO  Number of potentially valid call(s)
        #F_MISS  Missing call rate
            #https://www.cog-genomics.org/plink/1.9/formats#lmiss
    #IMPORTANT: 
        #I understand that heterozigous cases for hayploid regions (i.e., non-PAR X-Y regions) are not considered in the missing count, we should take care of this later! The cool thing is that they do not affect when calculating missing samples, so we would not lose samples because of these problematic cases as they do not increase the number of missing per sample
        #Also, obligatory missing are are not counted, for example, Y snps are obligatory missing in females, we should not count these.
        #https://www.cog-genomics.org/plink/1.9/formats#imiss


print_text("remove SNPs with low call rate to avoid their influence when filtering samples", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./05_remove_missing_snps/;\
    plink \
        --bfile ./04_remove_low_maf_snps/merged_batches_remove_low_maf \
        --geno 0.01 \
        --make-bed \
        --out ./05_remove_missing_snps/merged_batches_remove_missing_snps; \
    ls -l ./05_remove_missing_snps/")
        #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
            #https://www.cog-genomics.org/plink/1.9/filter
        #genotyping rate >0.99 recommended for PRS tutorial, i.e., remove SNPs with more than 1% (0.01) missing. this is more stringent than used by Ritchie tutorial, but it is ok. We lose a few thousand SNPs (see below) and it is important to be conservative given we are going to use this data for PRS.
            #https://www.nature.com/articles/s41596-020-0353-1#Sec4
            #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view


print_text("see the number of SNPs removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_snps_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./04_remove_low_maf_snps/merged_batches_remove_low_maf.bim); \
    n_snps_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./05_remove_missing_snps/merged_batches_remove_missing_snps.bim); \
    lost_snps=$(($n_snps_before - $n_snps_after)); \
    lost_snps_percent=$( \
        awk \
            -v l=$lost_snps \
            -v b=$n_snps_before \
            'BEGIN{print (l/b)*100}'); \
    if [[ $lost_snps_percent < 4 ]]; then \
        printf 'The number of SNPs lost due to call rate below 0.99 is: %s' \"$lost_snps\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 4% of SNPs due low-call rate';\
    fi")


print_text("create again the missing report", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/05_remove_missing_snps/; \
    plink \
        --bfile merged_batches_remove_missing_snps \
        --missing; \
    ls -l")
        #--missing produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.
            #https://www.cog-genomics.org/plink/1.9/basic_stats#missing


print_text("check that no SNP has a missing % above the selected threshold", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/05_remove_missing_snps/; \
    total_number_snps=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>1){count++}}END{print count}' \
            plink.lmiss);\
    n_snps_below_threshold=$( \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 <= 0.01)){count++}}END{print count}' \
            plink.lmiss); \
    if [[ $n_snps_below_threshold -eq $total_number_snps ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")

# endregion




############################
# region HWE INITIAL #######
############################

#If there are far more heterozygous calls than would be expected under Hardy–Weinberg equilibrium, that is usually due to a systematic variant calling error. Any such variants should be removed from the dataset.

#Plink tutorial does this after MAF filtering
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22

#Also note this is relevant for KING. In a thread, Christopher said to someone trying to fitler related samples with KING that he/she was not performing QC correctly due to the kinship values he/she was obtaining. He recommended to use HWE
    #"This indicates that you aren't performing proper QC.  More precisely, it is important to exclude variants that are not close to Hardy-Weinberg equilibrium, especially when that's because the variant-caller is screwing up royally and calling most or all samples heterozygous; you're in a tough position with only 29 samples, but you may still be able to detect and filter enough instances of this.  The Hardy-Weinberg p-value calculator at https://www.cog-genomics.org/software/stats/ may be helpful re: setting a --hwe threshold."
        #https://groups.google.com/g/plink2-users/c/NKBcu2cC160/m/3_IftZbEAwAJ

print_text("Initial HWE filter", header=2)
print_text("create folder for this step", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./06_first_hwe_filter/; \
    ls -l")

print_text("apply HWE filter", header=3)
#According to the plink tutorial, the following PLINK 2 command line is appropriate during initial quality control:
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./05_remove_missing_snps/merged_batches_remove_missing_snps \
        --hwe 1e-25 keep-fewhet midp \
        --make-bed \
        --out ./06_first_hwe_filter/merged_batches_hwe_filtered"
)
    #--hwe filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. 
        #https://www.cog-genomics.org/plink/1.9/filter#hwe
        #We recommend setting a low threshold. Serious genotyping errors often yield extreme p-values like 1e-50 which are detected by any reasonable configuration of this test, while genuine SNP-trait associations can be expected to deviate slightly from Hardy-Weinberg equilibrium (so it's dangerous to choose a threshold that filters out too many variants).
        #Therefore, using a very low p-value would only remove SNPs extremely deviates from HWE and hence, being potential genotpying errors.
        #According to the plink tutorial: The 1e-25 threshold is extreme enough that we are unlikely to remove anything legitimate. (Unless the dataset is primarily composed of F1 hybrids from an artificial breeding program.). 
    #The “keep-fewhet” modifier causes this filter to be applied in a one-sided manner. Therefore, the variants with fewer heterozygotes than expected under HWE possibly caused by population stratification WOULD NOT be filtered out by this command. We are only filtering SNPs with higher heterozygotes than expected under HWE.
        #Later, after you have a good idea of population structure in your dataset, you may want to follow up with a round of two-sided –hwe filtering, since large (see Note 5) violations of Hardy–Weinberg equilibrium in the fewer-hets-than-expected direction within a subpopulation are also likely to be variant calling errors; with multiple subpopulations, the –write-snplist and –extract flags can help you keep just the SNPs which pass all subpopulation HWE filters.
    #--hwe's 'midp' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium. The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data. We recommend its use.
        #We are using it, but it seems it does not change anything using it or not.
    #Because of the missing data issue, you should not apply a single p-value threshold across a batch of variants with highly variable missing call rates. A warning is now given whenever observation counts vary by more than 10%.
        #we should not have this problem because we have already clean the SNPs based on missingness
    #Note the notation of non-autosomal chromosomes si changed here from 23 to X, 24 to Y and so on... This is because we are using here plink2 to filter by HWE. Adding this step with plink2 changes the chromosome notation.
    #Once we create new bed files with plink1.9, then the previous notation comeback again.



print_text("Less than 22 SNPs have been removed by the initial HWE filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    n_snps_before_hwe=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./05_remove_missing_snps/merged_batches_remove_missing_snps.bim); \
    n_snps_after_hwe=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./06_first_hwe_filter/merged_batches_hwe_filtered.bim); \
    n_snp_diff=$(($n_snps_before_hwe - $n_snps_after_hwe)); \
    if [[ $n_snp_diff -lt 22 ]]; then \
        echo \"TRUE\"; \
    else \
        echo \"FALSE! ERROR!\"; \
    fi" \
)

# endregion




###############################
# region PROBLEM PAR XY #######
###############################

print_text("remove SNPs that are considered to be pseudoautosomals but in reality are not in pseudo-autosomal regions of XY", header=2)
##XY chromosome in final reports (chromosome 25 in bim file)
#General info about XY regions
#In an Illumina Final Report, SNPs in the XY chromosome refer to genetic variations that occur in the pseudoautosomal regions (PARs) of the X and Y chromosomes[1]. 
#These regions are present on both X and Y chromosomes and are capable of undergoing recombination during meiosis, similar to autosomal chromosomes[1]. Therefore, these SNPs may show male heterozygotes[1].
#The Illumina genotyping arrays, such as the ones used in the GenomeStudio software, identify these SNPs using pre-defined oligonucleotide probes designed to hybridize specific regions of genomic DNA[1]. The identity of alleles is determined by automated clustering of samples based on the similarity of fluorescent intensity[1].
#It's important to note that the XY chromosome should be treated as an autosomal chromosome when analyzing these SNPs[1]. This is because the PARs behave more like autosomal regions rather than sex-determining regions. Because of this, SNPs in these regions have a different code (XY or 25).
#Source: Conversation with Bing, 5/1/2024
#(1) GenomeStudio Genotyping QC SOP v.1.6 - GitHub Pages. https://khp-informatics.github.io/COPILOT/GenomeStudio_genotyping_SOP.html.
#(2) How to interpret DNA strand and allele information ... - Illumina Knowledge. https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001489.
#(3) tutorials:population-diversity:snp-chips [ILRI Research Computing] - CGIAR. https://hpc.ilri.cgiar.org/tutorials/population-diversity/snp-chips.
#(4) Infinium Genotyping Data Analysis - Illumina. https://www.illumina.com/documents/products/technotes/technote_infinium_genotyping_data_analysis.pdf.

#Specific details about our data
#It seems that the SNPs marked as XY in the illumina report (and as 25 in the bim files of plink) have the physical position in the chromosome X. Sometimes, the position is the same in both X and Y, but if the are not, then the position showed in the illumina report is the X position. See examples:
#rs1883078 has code XY in the final report and position 155870488. According to ncbi (variant details page; https://www.ncbi.nlm.nih.gov/snp/rs1883078#variant_details) this SNP is in position 155870488 of chrX (hg38), while it is in position 57057008 for chrY. In the bim file, this SNP has code 25 and position 155870488, i.e., chrX.
#rs5946743 has code XY in the final report and position 733497. According to ncbi (variant details page; https://www.ncbi.nlm.nih.gov/snp/rs5946743/#variant_details) this SNP is in position 733497 for both chrX and chrY (hg38). In the bim file, this SNP has code 25 and position 733497.
#rs28416357 has code XY in the final report and position 1308288. According to ncbi (variant details page; https://www.ncbi.nlm.nih.gov/snp/rs28416357/#variant_details) this SNP is in position 1308288 for both chrX and chrY (hg38). In the bim file, this SNP has code 25 and position 1308288.
#kgp22824102 has code XY in the final report and position 825992. According to ncbi (variant details page; https://www.ncbi.nlm.nih.gov/snp/rs5946480#variant_details) this SNP is in position 825992 for both chrX and chrY (hg38). In the bim file, this SNP has code 25 and position 825992. Note, SNPs that start with “kgp” are also genetic variations, similar to those that start with “rs”. The “kgp” prefix is used by the 1000 Genomes Project. The number following “kgp” is a unique identifier for that specific SNP. Please note that not all “kgp” SNPs may have a corresponding “rs” identifier. The “rs” identifiers are assigned by the dbSNP database when a SNP is submitted to them, and not all SNPs identified by the 1000 Genomes Project may have been submitted to dbSNP.

#Comparison of the PAR region in our data and the PAR region defined for hg38 by plink
#According to plink, in GRCh38/UCSC human genome 38, the boundaries are 2781479 and 155701383 (https://www.cog-genomics.org/plink/1.9/data#split_x). It seems these are the limits of the two PAR regions in the X chromosome. 
#Accoring to hg38 data (https://www.ncbi.nlm.nih.gov/grc/human), in chrX, the first pseudo-autosomal region starts at basepair 10001 and ends at basepair 2781479, being the latter the first boundary used by plink. The second region starts at basepair 155701383 and ends at basepair 156030895, being the former the second boundary indicated by Plink. In other words, plink considers the end of the first PAR region and the start of the second PAR region.
#Therefore, we can assume these are the correct coordinates of pseudo-autosomal regions in the X chromosome. Note that I have previously check in several cases that coordinates of XY SNPs in illumina reports are chrX coordinates.



print_text("create folder for this step", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./07_par_problem_cleanup/; \
    ls -l")



#in the last BIM file after MAF filters, select ID of those SNPs with code 25 that are outside the PAR regions previously indicated.
print_text("see SNPs that are considered to be pseudoautosomals (chr=25) but in reality are not in pseudo-autosomal regions of XY", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    awk \
        'BEGIN{FS=\"\t\"}{ \
            if($1 == 25 || $1 == \"XY\"){ \
                if($4 < 10001 || $4 > 2781479 && $4 < 155701383 || $4 > 156030895){ \
                    print $2 \
                } \
            } else if($1 == 23 || $1 == \"X\") { \
                if($4 >= 10001 && $4 <= 2781479 || $4 >= 155701383 && $4 <= 156030895){ \
                    print $2 \
                } \
            } else if($1 == 24 || $1 == \"Y\") { \
                if($4 >= 10001 && $4 <= 2781479 || $4 >= 56887903 && $4 <= 57217415){ \
                    print $2 \
                } \
            } \
        }' \
        ./06_first_hwe_filter/merged_batches_hwe_filtered.bim > ./07_par_problem_cleanup/snps_par_problem.txt"
)
#load bim file after MAF filtering to awk using tabs as delimiter (checked this is the delimiter in the file)
#select those SNPs in pseudo autosomal regions (code 25) that are located: 1) before the start of the first PAR region; 2) After the first PAR region and before the second one; 3) After the second PAR region. Then print their ID. 
#else, if the SNP is considered to be in the X chromosome, but within a PAR region, also print its ID.
    #Note I am using the PAR coordinates accoring to GRC38 for chromosome X.
#else, if the SNP is considered to be in the Y chromosome, but within a PAR region, also print its ID.
    #Note I am using the PAR coordinates accoring to GRC38 for chromosome Y.
        #https://www.ncbi.nlm.nih.gov/grc/human



#check we only have 21 problematic cases
print_text("check we ONLY have 21 of these problematic cases", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/07_par_problem_cleanup; \
    n_par_problem=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./snps_par_problem.txt); \
    if [[ $n_par_problem -eq 21 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi"
)



#make a file with ID of NON problematic SNPs
print_text("create a list WITHOUT SNPs that are problematic for PAR: We discard SNPs that are considered to be pseudoautosomals but in reality are not in pseudo-autosomal regions of XY", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    awk \
        'BEGIN{FS=\"\t\"}{ \
            if($1 == 25 || $1 == \"XY\"){ \
                if($4 >= 10001 && $4 <= 2781479 || $4 >= 155701383 && $4 <= 156030895){ \
                    print $2 \
                } \
            } else if($1 == 23 || $1 == \"X\"){ \
                if($4 < 10001 || $4 > 2781479 && $4 < 155701383 || $4 > 156030895){ \
                    print $2 \
                } \
            } else if($1 == 24 || $1 == \"Y\"){ \
                if($4 < 10001 || $4 > 2781479 && $4 < 56887903 || $4 > 57217415){ \
                    print $2 \
                } \
            } else { \
                print $2 \
            } \
        }' \
        ./06_first_hwe_filter/merged_batches_hwe_filtered.bim > ./07_par_problem_cleanup/snps_par_no_problem.txt"
)
#if the chromosome is 25 (XY) and the SNPs is within the PAR limits accoridng to NCBI print the ID (second column)
#if the chromosome is 23 (X) and the SNPs is outside the PAR limits for X, print its ID
    #using hg38 coordinates of PAR for X
#if the chromosome is 24 (Y) and the SNPs is outside the PAR limits for Y, print its ID
    #using hg38 coordinates of PAR for Y
    #https://www.ncbi.nlm.nih.gov/grc/human
#if the SNP is not considered PAR, X nor Y, we can just print its ID



#retain only these SNPs without the problem
#We cannot be sure if there is something wrong with these SNPs
#They are in PAr regions or not? should they be treated as PAR o like
#sex chromosomes? so we are going to check the number is not very high
#and then remove all of them. We should be ok with the remaining PAR SNPs
#as they are considered separately from autosomals and sex chromosomes
print_text("retain only these SNPs without the problem", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./06_first_hwe_filter/merged_batches_hwe_filtered \
        --extract ./07_par_problem_cleanup/snps_par_no_problem.txt \
        --make-bed \
        --out ./07_par_problem_cleanup/merged_batches_hwe_par; \
    ls -l ")
        #--bfile: 
            #This flag causes the binary fileset plink.bed + plink.bim + plink.fam to be referenced. If a prefix is given, it replaces all instances of 'plink', i.e., it looks for bed, bim and fam files having that suffix instead of "plink".
            #https://www.cog-genomics.org/plink/1.9/input#bed
        #--extract:
            #normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces) and removes all unlisted variants from the current analysis
            #https://www.cog-genomics.org/plink/1.9/filter#snp
        #--make-bed creates a new PLINK 1 binary fileset, AFTER applying sample/variant filters and other operations
            #https://www.cog-genomics.org/plink/1.9/data#make_bed



#quick check
print_text("quick check", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    n_snps_before_par=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./06_first_hwe_filter/merged_batches_hwe_filtered.bim); \
    n_snps_after_par=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./07_par_problem_cleanup/merged_batches_hwe_par.bim); \
    n_snps_par_problem=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./07_par_problem_cleanup/snps_par_problem.txt); \
    sum_snps=$(($n_snps_after_par + $n_snps_par_problem)); \
    if [[ $n_snps_before_par -eq $sum_snps ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi"
)

# endregion




################################
# region sample missingness ####
################################
print_text("filter by sample missingness", header=2)
#The usual practice is to filter out the samples and variants with high missing-entry frequencies; these tend to be caused by mistakes in the lab, bad SNP probes, variant calling limitations, and similar issues where throwing out the entire row/column is an appropriate solution. (You still have lots of other rows and columns to work with, so at least in population genomics it usually is not worth the effort to try to salvage it.) 
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22

print_text("create folder for this step", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./08_remove_low_call_samples/; \
    ls -l")



print_text("make sample missing report after previous filters", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./07_par_problem_cleanup/merged_batches_hwe_par \
        --missing \
        --out ./08_remove_low_call_samples/merged_batches_hwe_par_miss_report;" \
)
    #--missing produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.
        #https://www.cog-genomics.org/plink/1.9/basic_stats#missing
run_bash("head ./data/genetic_data/quality_control/08_remove_low_call_samples/merged_batches_hwe_par_miss_report.imiss")
    #imiss: A text file with a header line, and one line per sample with the following six fields:
        #FID    Family ID
        #IID Within-family ID
        #MISS_PHENO  Phenotype missing? (Y/N)
        #N_MISS  Number of missing genotype call(s), not including obligatory missings or heterozygous haploids
        #N_GENO  Number of potentially valid call(s)
        #F_MISS  Missing call rate
            #https://www.cog-genomics.org/plink/1.9/formats#imiss
    #IMPORTANT: 
        #I understand that heterozigous cases for haploid regions (i.e., non-PAR X-Y regions that should not have a second allele because there is no correspondance between X and Y) are not considered in the missing count, we should take care of this later! The cool thing is that they do not affect when calculating missing samples, so we would not lose samples because of these problematic cases as they do not increase the number of missing per sample
        #Also, obligatory missing are are not counted, for example, Y snps are obligatory missing in females, we should not count these.


print_text("check that F_MISS is just the number of missing genotype divided by the number of potentially valid calls", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/08_remove_low_call_samples/; \
    n_samples_freq_file=$( \
        awk \
            'BEGIN{FS=\" \"; OFS=\" \"} \
            {if(NR>1)(count++)}\
            END{print count}'\
            merged_batches_hwe_par_miss_report.imiss \
    ); \
    n_correct_freq_miss=$( \
        awk \
            'BEGIN{ \
                FS=\" \"; OFS=\" \" \
            } \
            { \
                if((int(($4/$5)*100+0.5)/100==int($6*100+0.5)/100) && (NR>1)){ \
                    count++ \
                } \
            } \
            END { \
                print count \
            }' \
            merged_batches_hwe_par_miss_report.imiss \
    ); \
    if [[ $n_correct_freq_miss -eq $n_samples_freq_file ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #just check if freq miss is the ratio of missing calls respect to the potential number of valid calls per sample. We count the number of samples for which this is the case, N_MISS/N_GENO is equal to F_MISS and then check this is the total number of samples, i.e., all meet the condition.
        #see the SNP missing part for explanations about the script


print_text("load the report with awk and then save it controlling the delimiter. If I load it directly into pandas, I got problems separating the columns", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/08_remove_low_call_samples/; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{if(NR>0){print $1,$2,$3,$4,$5,$6}}' \
        merged_batches_hwe_par_miss_report.imiss \
    > merged_batches_hwe_par_miss_report_awk_processed.imiss; \
    head merged_batches_hwe_par_miss_report_awk_processed.imiss")
        #specify the delimiter of the output (OFS="\t") and just print all fields (1 to 6) for all rows, then save as a new file


print_text("load in pandas", header=3)
sample_missing_report = pd.read_csv( \
    "./data/genetic_data/quality_control/08_remove_low_call_samples/merged_batches_hwe_par_miss_report_awk_processed.imiss", \
    sep="\t", \
    header=0, \
    low_memory=False)
print(sample_missing_report)


print_text("check we have the correct number of samples and columns in the missing sample report", header=3)
print(sample_missing_report.shape[0]==1458) 
print(sample_missing_report.shape[1]==6) 


print_text("calculate the proportion of samples retained at different missing thresholds", header=3)
#calculate the thresholds over 1
thresholds = [x / 1000 for x in range(980, 1000, 1)]
    #we want from proportion 0.985 to 1. In order to get floats with range, we need to first use integers that include all the numbers we want to include as decimals (e.g., 985 for 0.985) and then divide by 1000 (985/1000=0.985).
        #https://stackoverflow.com/a/7267287/12772630
print("see the thresholds selected")
print(thresholds)

#calculate the proportion of retained samples
list_tuples_thresholds = [] #empty list to save results
#threshold=thresholds[0]
for threshold in thresholds:
    
    #get the number of samples above the missing threshold    
    n_retained_samples = sample_missing_report \
        .loc[(1-sample_missing_report["F_MISS"]) >= threshold,:]\
        .shape[0]
            #1 minus F_MISS gives the call rate over 1
    
    #calculate proportion of retained samples dividing retained samples by the total number of samples
    proportion_remaining = n_retained_samples/sample_missing_report.shape[0]
    
    #append to the result list a tuple with threshold and proportion
    list_tuples_thresholds.append((threshold, proportion_remaining))

#convert the list to DF
pd_tuples_thresholds = pd.DataFrame(list_tuples_thresholds, columns=["threshold", "proportion_remaining"])
print("see thresholds and proportion of retained samples")
print(pd_tuples_thresholds)



print_text("plot sample missing thresholds against the proportion of samples retained", header=3)
import matplotlib.pyplot as plt
plt.plot( \
    pd_tuples_thresholds["threshold"], \
    pd_tuples_thresholds["proportion_remaining"], \
    marker="o", \
    linestyle="dashed")
plt.xlabel("Call Rate Threshold")
plt.ylabel("Proportion Remaining")
plt.savefig( \
    fname="./data/genetic_data/quality_control/08_remove_low_call_samples/merged_batches_plot_sample_retention_thresholds.png")
plt.close()


print_text("remove samples with low call rate", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./07_par_problem_cleanup/merged_batches_hwe_par \
        --mind 0.01 \
        --make-bed \
        --out ./08_remove_low_call_samples/merged_batches_hwe_par_callrate")
            #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
                #https://www.cog-genomics.org/plink/1.9/filter
            #recommended threshold
                #the PRS tutorial recommends to retain samples with "sample missingness <0.02", i.e., samples with more than 2% missing or less than 98% call rate should be removed. 
                #The Ritchie tutorials says "A recommended threshold is 98% to 99% efficiency after first removing markers that have a low genotype call rate across samples."
                #The illumina report of the first batch says "Four samples are below the illumina expected 99% SNP call rate (values expected for typical projects, excluding tumour samples)".
                #we are going to use 99%, i.e., < 0.01 missingness. This is higher than recommended by PRS tutoria, but within the range recommended by Ritchie tutorial and the value for Illumina. In the second batch, this threshold does not increase the number of samples removed compared to 98%, while in the first one only makes 1 more sample to be removed, one with call rate=0.981. Remember that the FPD report says that this is below the expectation of Illunina. Therefore, I think we are ok removing these samples.
            #The family and sample IDs of removed individuals are indicated in a file with the extension ".irem"



print_text("see the number of samples removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_samples_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./07_par_problem_cleanup/merged_batches_hwe_par.fam); \
    n_samples_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./08_remove_low_call_samples/merged_batches_hwe_par_callrate.fam); \
    n_samples_removed=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./08_remove_low_call_samples/merged_batches_hwe_par_callrate.irem); \
    lost_samples=$(($n_samples_before - $n_samples_after)); \
    lost_samples_percent=$( \
        awk \
            -v l=$lost_samples \
            -v b=$n_samples_before \
            'BEGIN{print (l/b)*100}'); \
    if [[ $lost_samples_percent < 2 && $n_samples_removed -eq $lost_samples ]]; then \
        printf 'The number of samples lost due to call rate below 0.99 is: %s' \"$lost_samples\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 2% of samples due low-call rate';\
    fi")
        #see SNP missing code to details about awk steps



print_text("create again the missing report", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/08_remove_low_call_samples/; \
    plink \
        --bfile merged_batches_hwe_par_callrate \
        --missing \
        --out merged_batches_hwe_par_callrate_missing; \
    ls -lh")
        #--missing produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.
            #https://www.cog-genomics.org/plink/1.9/basic_stats#missing


print_text("check that no sample has a missing % above the selected threshold", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/08_remove_low_call_samples/; \
    total_number_samples=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>1){count++}}END{print count}' \
            merged_batches_hwe_par_callrate_missing.imiss);\
    n_samples_below_threshold=$( \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($6 <= 0.01)){count++}}END{print count}' \
            merged_batches_hwe_par_callrate_missing.imiss); \
    if [[ $n_samples_below_threshold -eq $total_number_samples ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #see SNP missing code to details about awk steps

# endregion




####################################################################
# region samples with LogR SD above the illumina expectation #######
####################################################################
print_text("remove samples with LogR SD above the illumina expectation", header=2)
#In the PDF summary of the first batch they say that "Three samples are above the illumina expectations of < 0.3 for LogR SD. Details can be found on page 3 of this report.". Remember that LogR is a parameter mentioned in the tutorials to detect sex inconsistences, and it is mentioned in the pdf report of batch 1, so it makes sense to use it.
#Some the samples above the expectation have also a call rate below 0.99, but not others, thus a specific filter is needed for this. Given that Illumina says that level of LogR SD is not normal, we should remove these individuals.


print_text("obtain the samples above the expectation for logR dev in each batch", header=3)
#batch_name = "ILGSA24-17873"
for batch_name in ["ILGSA24-17873", "ILGSA24-17303"]:
    
    print_text(batch_name + ": process the file with the logR deviation", header=4)
    if(batch_name == "ILGSA24-17873"):
        cnmetricts_file_name="CAGRF20093767_CNMetrics"
    elif(batch_name == "ILGSA24-17303"):
        cnmetricts_file_name="ILGSA24-17303_CNMetrics"
    print(cnmetricts_file_name)
    
    print_text(batch_name + ": get the ID of samples above the illumina expectation for logR SD", header=3)
    removed_samples_logRdev = run_bash(" \
        cd ./data/genetic_data/cn_metrics/; \
        awk \
            'BEGIN{ \
                FS=\",\"; \
                OFS=\"\t\"} \
            { \
                if((NR>2) && ($9 >= 0.3)){ \
                    print \"combat_" + batch_name + "\", $1 \
                }\
            }' \
            " + cnmetricts_file_name + ".csv > " + batch_name + "_samples_filter_out_LogRDev.tsv; \
        awk \
            'BEGIN{FS=\"\t\"}{print $2}' \
            " + batch_name + "_samples_filter_out_LogRDev.tsv", return_value=True).strip().split("\n")
            #load the CNMetrics CSV file into awk indicating that sep is ",", but the output should be tab
            #print only from the third row (the two first have complicated headers) and only those rows for which the 9th column (the one with logR SD) is above the expectation (i.e., 0.3).
                #from that rows, only print the first column, which is the ID.
                #also print the name of the family first, because plink --remove needs a file with two columns, the family ID and the within-family ID
                    #in our case, the family ID is the batch name, which is "combat_" + batch_name
                    #this can be add as another column just with print and ""
                        #https://stackoverflow.com/questions/7551991/add-a-new-column-to-the-file
            #the resulting IDs of samples above the expectation can be saved as a TSV file
            #also get the output in python
                #remove the empty spaces with strip and then split by "\n", in this way we get all IDs as different elements of a list and the space at the end is not considered.
    print(removed_samples_logRdev)


print_text("combine the IDs from both batched into one single file", header=3)
run_bash(" \
    cd ./data/genetic_data/cn_metrics/; \
    cat \
        ILGSA24-17303_samples_filter_out_LogRDev.tsv  \
        ILGSA24-17873_samples_filter_out_LogRDev.tsv > \
    merged_batches_samples_filter_out_LogRDev.tsv; \
    cat merged_batches_samples_filter_out_LogRDev.tsv")
        #https://stackoverflow.com/a/2150794/12772630


print_text("check we have the correct IDs in the merged file", header=3)
print_text("get IDs from the merged file", header=4)
total_removed_samples_logRdev = run_bash(" \
    awk \
        'BEGIN{FS=\"\t\"}{print $2}' \
        ./data/genetic_data/cn_metrics/merged_batches_samples_filter_out_LogRDev.tsv", return_value=True).strip().split("\n")

print_text("combine the list of ID samples from both batches and check it is identical than the merged list", header=4)
print( \
    run_bash(" \
        awk \
            'BEGIN{FS=\"\t\"}{print $2}' \
            ./data/genetic_data/cn_metrics/ILGSA24-17303_samples_filter_out_LogRDev.tsv", return_value=True).strip().split("\n") + \
    run_bash(" \
        awk \
            'BEGIN{FS=\"\t\"}{print $2}' \
            ./data/genetic_data/cn_metrics/ILGSA24-17873_samples_filter_out_LogRDev.tsv", return_value=True).strip().split("\n") == \
    total_removed_samples_logRdev)



print_text("remove these samples", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        09_remove_high_LogRDev_samples; \
    plink \
        --bfile ./08_remove_low_call_samples/merged_batches_hwe_par_callrate \
        --remove ../cn_metrics/merged_batches_samples_filter_out_LogRDev.tsv \
        --make-bed \
        --out ./09_remove_high_LogRDev_samples/merged_batches_hwe_par_callrate_LogRDev; \
    ls -l ./09_remove_high_LogRDev_samples/")
            #create a new folder to save plink filesets after filtering
            #then use --remove to remove samples included in a file with two columns: family id and within-family id.
                #this file is in a different parent folder so we have to use "../"
                #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.
                    #https://www.cog-genomics.org/plink/1.9/filter


print_text("check that the corresponding samples are indeed not included in the new fam file", header=3)
fam_file_after_logRdev_filtering = pd.read_csv( \
    "./data/genetic_data/quality_control/09_remove_high_LogRDev_samples/merged_batches_hwe_par_callrate_LogRDev.fam", \
    sep=" ", \
    header=None, \
    low_memory=False)
print(sum(fam_file_after_logRdev_filtering.iloc[:, 1].isin(total_removed_samples_logRdev)) == 0)
    #no sample ID (second column in fam file) should be included in the list of IDs to be removed


print_text("see the number of samples removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_samples_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./08_remove_low_call_samples/merged_batches_hwe_par_callrate.fam); \
    n_samples_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./09_remove_high_LogRDev_samples/merged_batches_hwe_par_callrate_LogRDev.fam); \
    lost_samples=$(($n_samples_before - $n_samples_after)); \
    lost_samples_percent=$( \
        awk \
            -v l=$lost_samples \
            -v b=$n_samples_before \
            'BEGIN{print (l/b)*100}'); \
    if [[ $lost_samples_percent < 2 && $lost_samples_percent > 0 ]]; then \
        printf 'The number of samples lost due to logR SD above illumina expectations is: %s; Note that we previously filtered by missingness, so we already removed there some samples with logR dev problems' \"$lost_samples\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 2% of samples due logR SD above illumina expectations OR just zero samples removed, which would be also strange, maybe you are using the wrong CNMetric file';\
    fi")
        #see SNP missing code to details about awk steps
        #just added another condition, if the number of removed samples is zero, get error because we have problematic samples for this in both batches, so maybe you are using the incorrect CNMetric file

# endregion




##############################
# region MAF SECOND ROUND ####
##############################
print_text("MAF filtering second round", header=2)
#We need decent MAFs for the next step (sample relatdness with King), so we are to check again MAF just in case the previous removal of samples has made some SNPs to be now below our MAF threshold.
#see previous MAF round for info about the selected threshold

print_text("apply the filter", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./10_remove_low_maf_snps/; \
    plink \
        --bfile ./09_remove_high_LogRDev_samples/merged_batches_hwe_par_callrate_LogRDev \
        --maf 0.05 \
        --make-bed \
        --out ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf; \
    ls -l ./10_remove_low_maf_snps")
            #create a new folder to save filesets after removing related samples
            #apply the MAF filter to the fileset filtered for logR SD
                #--maf filters out all variants with minor allele frequency below the provided threshold (default 0.01), while --max-maf imposes an upper MAF bound. Similarly, --mac and --max-mac impose lower and upper minor allele count bounds, respectively.
                    #https://www.cog-genomics.org/plink/1.9/filter


print_text("see the number of SNPs removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_snps_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./09_remove_high_LogRDev_samples/merged_batches_hwe_par_callrate_LogRDev.bim); \
    n_snps_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf.bim); \
    lost_snps=$(($n_snps_before - $n_snps_after)); \
    lost_snps_percent=$( \
        awk \
            -v l=$lost_snps \
            -v b=$n_snps_before \
            'BEGIN{print (l/b)*100}'); \
    if [[ $lost_snps_percent < 1 ]]; then \
        printf 'The number of SNPs lost due to MAF below 0.05 is: %s' \"$lost_snps\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 60% of SNPs due low-call rate';\
    fi")
        #First calculate the number of SNPs, i.e., rows, in the bim files before and after filtering using for that "++" to count cases.
        #Calculate the difference of snps, i.e., SNPs removed.
        #calculate the percentage respect to the initial number of SNPs
            #use the -v flag to load variables into awk, so we can use the number of SNPs removed and those previously present before the filtering. Use awk to calculate the percentage.
                #https://stackoverflow.com/a/12147154/12772630
                #https://www.unix.com/unix-for-dummies-questions-and-answers/15651-awk-v.html
            #bash cannot do operations with floats
        #if the percentage of lost SNPs is less than 1%, we are fine. So print the number of SNPs using printf, so you can combine string message with a variable, indicated as %s because it is loaded as string.
            #https://phoenixnap.com/kb/bash-printf


print_text("create freq report", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/10_remove_low_maf_snps/; \
    plink \
        --bfile merged_batches_hwe_par_callrate_LogRDev_maf \
        --freq; \
    head plink.frq")
        #--frq produces a frequency report
            #A text file with a header line, and then one line per variant with the following six fields:
                #CHR Chromosome code
                #SNP Variant identifier
                #A1  Allele 1 (usually minor)
                #A2  Allele 2 (usually major)
                #MAF Allele 1 frequency
                #NCHROBS Number of allele observations
            #https://www.cog-genomics.org/plink/1.9/formats#frq


print_text("check that no SNP has a MAF below the selected threshold", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/10_remove_low_maf_snps/; \
    total_number_snps=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>1){count++}}END{print count}' \
            plink.frq);\
    n_snps_above_threshold=$( \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 >= 0.05)){count++}}END{print count}' \
            plink.frq); \
    if [[ $n_snps_above_threshold -eq $total_number_snps ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #count the number of rows of the frequency report, i.e., total number of SNPs
        #then count the number of SNPs in that file that have a MAF equal or higher than 0.05, which is our selected threshold.
        #the number of snps meeting this condition should be the same than the total number of SNPs, because we have already applied the filter.

# endregion




###################################
# region sample relatedness #######
###################################
print_text("sample relatedness", header=2)
#the plink tutorial first remove related samples before filtering by MAF and calculate the PCA.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
#sample relatedness is also explained before than population substructure in Ritchie's tutorial. Indeed, they also talk about LD pruning, which is an step needed for selecting SNPs for the PCA.
    #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
#In the VO2 max paper, they also seem to remove the related individuals using pi-hat before doing the PCA
    #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7

#We have already filtered SNPs by MAF so we have enough MAF for King (see summary of steps above). We filter by relatedness and then filter again by MAF ensuring the removal of all individuals has not make some SNPs to have very low MAF. Then we can go to pop stratification,
    #The methods we are gonna use to remove related samples (KING-robust) in plink2 seems to require decent MAF. See this from plink2 help:
        #"The relationship matrix computed by --make-rel/--make-grm-list/--make-grm-bin can be used to reliably identify close relations within a single population, IF YOUR MAFS ARE DECENT"

#Note about the removal of samples
    #First, we remove samples by missing and by logRDev. From the remaining individuals, we calculate kindship, in this way we avoid to consider samples that would be removed anyways due to missingness
    #Note that having related samples should not influence missingness. It should NOT influence the amount of genotypes correctly genotype in a sample, so we can do it before removeing related samples.
    #The only problem here is MAF, because the previous removal of samples could be bad for KING. Because of this, we apply again the MAF filter before this step.
    #Yes, the new MAF filter could make some samples to have now low missingness, but this should not be a big problem for king. If a sample has low call rate, it could make difficult to check its relationthip with other samples, but that is ok, because that sample would be removed anyways, in this step (we remove the sample with lowest cal rrat ein each pair) or later with the last MAF loop.
        #Also note that we are going to use King-robust which should be robust to this and also robust to mixed populations.
            #The relationship matrix computed by --make-rel/--make-grm-list/--make-grm-bin can be used to reliably identify close relations within a single population, if your MAFs are decent. However, Manichaikul et al.'s KING-robust estimator can also be mostly trusted on mixed-population datasets (with one uncommon exception noted below), and doesn't require MAFs at all. Therefore, we have added this computation to PLINK 2, and the relationship-based pruner is now based on KING-robust.

#Rationale:
    #Many population-genomic statistics (such as the allele frequencies) and analyses are distorted when there are lots of very close relatives in the dataset; you are generally trying to make inferences about the population as a whole, rather than a few families that you oversampled. For example, PLINK 2 includes an implementation of the KING-robust [5] pairwise relatedness estimator, which can be used to prune all related pairs. This does not mean that both samples in each related pair are thrown out. Instead, –king-cutoff tries to keep as much data as possible, and as a consequence it usually keeps one sample out of each pair (see below).
    #These related samples, if treated as independent samples in the downstream analyses, having many related samples in the dataset would result in increased type I and type II errors. The options are tu use of mixed-regression models while considering in place of simple linear/logistic regression or remove one member of each pair.
    #Retaining a maximal set of unrelated individuals is computationally expensive (NP-hard), but an efficient greedy approximation is available in PLINK 2.0 using the --king-cutoff flag (as opposed to just removing one of each pair of related individuals)
        #see below about KING
    #Cryptic relatedness can interfere with the association analysis. If you have a family‐based sample (e.g., parent‐offspring), you do not need to remove related pairs but the statistical analysis should take family relatedness into account. However, for a population based sample we suggest to use a pi‐hat threshold of 0.2, which in line with the literature (Anderson et al., 2010; Guo et al., 2014).
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/

print_text("create folder for this step", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./11_remove_related_samples/; \
    ls -l")



print_text("calculate kinship coefficient", header=3)
#Using dense marker data obtained in GWAS, it is easy to pairwise kinship estimates between every individual in the study using the --genome option in PLINK 1.9. This procedure need not be performed on the entire GWAS dataset; using a linkage-disequilibrium (LD) pruned dataset consisting of independent loci yields stable estimates of relationship in the form of kinship coefficient or ˆπ (pi-hat) values.
    #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
#This is a common measure of relatedness (or duplication) between pairs of samples and it is based on identity by descent (IBD). Typically, the individual of a related pair with lower genotype call rate is removed.
    #https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
print_text("LD-pruning", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf \
        --indep-pairwise 500kb 1 0.2 \
        --out ./11_remove_related_samples/ldpruned_snplist; \
    ls -l ./11_remove_related_samples/")
        #This command produces a pruned subset of variants that are in approximate linkage equilibrium with each other, writing the IDs to plink2.prune.in (and the IDs of all excluded variants to plink2.prune.out). These files are valid input for --extract/--exclude in a future PLINK run; and, for backward compatibility, they do not affect the set of variants in the current run.
        #Since the only output of these commands is a pair of variant-ID lists, they now error out when variant IDs are not unique.
        #--indep-pairwise is the simplest approach, which only considers correlations between unphased-hardcall allele counts. We cannot use the alternative, which is --indep-pairphase and requires phased data. --indep-pairwise takes three parameters: 
            #a required window size in variant count or kilobase (if the 'kb' modifier is present) units, 
            #an optional variant count to shift the window at the end of each step (default 1, and now required to be 1 when a kilobase window is used).
                #I guess the window is shifted (moved forward) until the variant count changes in 1 unit?
            #a required r2 threshold.
        #At each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, and variants are greedily pruned from the window until no such pairs remain.
        #For example:
            #"--indep-pairwise 100kb 1 0.8"
            #removes SNPs so that no pair within 100 kilobases have squared-allele-count-correlation (r2) greater than 0.8, and saves the IDs of the remaining SNPs.
        #By default, when given a choice, the variant-pruning commands preferentially keep variants with higher nonmajor allele frequencies. However, if you provide a list of variant IDs to --indep-preferred, all variants in that list are prioritized over all variants outside it. (Allele frequencies will still be used for tiebreaking.)
            #it seems reasonable to select the SNP with the higher MAF within pairs of correlated SNPs.
        #On human data, some reasonable parameter settings are, in order of increasing strictness:
            #"--indep-pairwise 100kb 1 0.8"
            #"--indep-pairwise 200kb 1 0.5"
            #"--indep-pairwise 500kb 1 0.2"
        #As we get more strict, the number of selected SNPs is reduced. We will go for the most strict option, i.e., allowing a lower correlation between SNPs in larger windows. This means that we will have much less variants correlated across larger chunks of the genome. In this way, we will have a very clean set in terms of LD. This is important because we are going to do operations that are not LD-aware.
        #https://www.cog-genomics.org/plink/2.0/ld
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
print("Do we have at least 70K SNPs after LD pruning like in Ritchie's tutorial?")
#In the Ritchie's tutorial, they ended up with 67,000 autosomal indepent variants in order to calculate IBD and pi_hat.
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    echo 'see first lines of the .prune.in list:'; \
    head ./ldpruned_snplist.prune.in; \
    ld_snps_in=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./ldpruned_snplist.prune.in); \
    printf 'The number of included SNPs is %s\n' \"$ld_snps_in\"; \
    echo 'Is this number greater than 70K?'; \
    if [[ $ld_snps_in -gt 70000 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #you can print a variable with text using printf and %s for "string". Then you call the variable you want to add within "". You could use two variables: "$var1 $var2"
            #https://phoenixnap.com/kb/bash-printf
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf \
        --extract ./11_remove_related_samples/ldpruned_snplist.prune.in \
        --make-bed \
        --out ./11_remove_related_samples/merged_batches_hwe_par_callrate_LogRDev_maf_ld_prunned;\
    ls -l ./11_remove_related_samples/")
        #from the current fileset, select only those SNPs included in .prune.in
        #we use extract for that
            #--extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.
                #https://www.cog-genomics.org/plink/1.9/filter
        #make a new fileset

print_text("remove sex chromsomes", header=4)
#According to Marees et al. (2018), we should check for sample relatedness not only using independent SNPs (pruning), but also limiting the analysis to autosomal chromosomes only.
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
#In Ritchie's tutorial (figure 5), they show an histogram for the distribution of pairwise pi-hat. They calculated IBD after "removing sex inconsistent individuals, 95% SNP call rate, 90% sample call rate, 10% MAF, and pruning to 67,000 AUTOSOMAL variants.". Also, in Ritchie's GitHub, they remove sex chromosomes before PCA: "Exclude any SNPs that do not liftOver and non-somatic chromosomes (X, Y)"
    #Therefore they do not use SNPs in sex chromosomes!
    #Respect sex inconsistencies, we can have minorities within the sample, and as plink help says (--check-sex), imbalanced ancestries can give problems so in that case you have to do the check of sex within each ancestry group. Therefore we need to check the PCA before.
    #plink tutorial check sex after pop stratification
#They are talking specifically about autosomals, so we are going to consider ONLY autosomals, no X, Y, MT nor PAR regions. The case I was doubting more was PAR regions because these behave like autosomal chromosomes, but they are sexual chromosomes. They are just 500 in total, so we are going to remove them.

#Therefore, I think we can use this set of pruned autosomal SNPs for kinship and PCA.
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    plink \
        --bfile ./merged_batches_hwe_par_callrate_LogRDev_maf_ld_prunned \
        --autosome \
        --make-bed \
        --out ./merged_batches_hwe_par_callrate_LogRDev_maf_ld_prunned_autosomals;\
    ls -l")
        #--autosome excludes all unplaced and non-autosomal variants, while --autosome-xy does not exclude the pseudo-autosomal region of X
            #https://www.cog-genomics.org/plink/1.9/filter
print("Do we have at least 70K autosomal SNPs after LD pruning like in Ritchie's tutorial?")
#In the Ritchie's tutorial, they ended up with 67,000 autosomal variants in linkage equilibrium in order to calculate IBD and pi_hat.
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    auto_ld_snps_in=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./merged_batches_hwe_par_callrate_LogRDev_maf_ld_prunned_autosomals.bim); \
    printf 'The number of included autosomal SNPs in linkage equilibrium is %s\n' \"$auto_ld_snps_in\"; \
    echo 'Is this number greater than 70K?'; \
    if [[ $auto_ld_snps_in -gt 70000 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #you can print a variable with text using printf and %s for "string". Then you call the variable you want to add within "". You could use two variables: "$var1 $var2"
            #https://phoenixnap.com/kb/bash-printf
print("All non-autosomals SNPs have been removed from the LD pruned dataset?")
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    n_non_auto_snps=$( \
        awk \
            'BEGIN{FS=\" \"}{if($1==0 || $1==23 || $1==24 || $1==25 || $1==26 || $1== \"X\"|| $1==\"Y\" || $1==\"XY\" || $1==\"MT\"){count++}}END{print count}'\
            ./merged_batches_hwe_par_callrate_LogRDev_maf_ld_prunned_autosomals.bim); \
    if [[ $n_non_auto_snps -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
    #with awk, count rows (i.e., SNPs) having in field 1 (i.e., chromosome) number 0 (unplaced SNP), number 23 (chromosome X), number 24 (chromosome Y), number 25 (XY autosomal) or number 26 (MT).
    #that count should be zero.

print_text("calculate the kindship", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./11_remove_related_samples/merged_batches_hwe_par_callrate_LogRDev_maf_ld_prunned_autosomals \
        --genome \
        --out ./11_remove_related_samples/ibd_report;\
    ls -l ./11_remove_related_samples/")
        #These calculations ARE NOT LD-aware. It is usually a good idea to perform some form of LD-based pruning before invoking them.
        #--genome invokes an IBS/IBD computation, and then writes a report plink.genome
            #Note that there is one entry per pair of samples, so this file can be very large.
        #This estimator requires fairly accurate minor allele frequencies to work properly. Use --read-freq if you do not think your immediate dataset's empirical MAFs are representative.
        #modifiers for genome
            #The 'gz' modifier causes the output to be gzipped
            #The 'full' modifier causes the following fields to be added:
                #IBS0    Number of IBS 0 nonmissing variants
                #IBS1    Number of IBS 1 nonmissing variants
                #IBS2    Number of IBS 2 nonmissing variants
                #HOMHOM  Number of IBS 0 SNP pairs used in PPC test
                #HETHET  Number of IBS 2 het/het SNP pairs used in PPC test
            #'rel-check' removes pairs of samples with different FIDs
                #I understand that this avoid comparison between samples of different families (batches in our case). We do NOT want that because we want to calculate relatdness between all samples across both batches.
                #rel-check in make-king-table does exactly that, so I understand this is the function of this modifier. One asked for using the same modifier 'rel-check' of --genome to make-king-table
                    #I was just wondering if you would consider adding the rel-check flag from plink1 --genome to plink2 --make-king-table? When working with family data I often find I am interested in two subsets of the full relatedness matrix: all relations within FID, and all relations between FID with kinship above some threshold. The latter can of course be easily produced using --king-table-filter (and then ignoring entries with FID1==FID2). The within-FID table on the other hand, I can currently only get by filtering the full (often huge) king table, or by crafting a dummy .kin0 file listing only pairs within FID and using with --king-table-subset. Both methods are feasible, but somewhat inconvenient. Thus my feature request.
                        #https://groups.google.com/g/plink2-users/c/aCxpf3j6tXg/m/hACv878gCwAJ
                    #The 'rel-check' modifier causes only same-FID pairs to be reported
                        #https://www.cog-genomics.org/plink/2.0/distance#make_king
            #unbounded and nudge
                #The underlying P(IBD=0/1/2) estimator sometimes yields numbers outside the range [0,1]; by default, these are clipped. The 'unbounded' modifier turns off this clipping. Then, if PI_HAT^2 < P(IBD=2), nudge adjusts the final estimates to P(IBD=0) := (1-p2), P(IBD=1) := 2p(1-p), and P(IBD=2) := p2, where p is the current PI_HAT.
                #not fully understood, just use default for this, i.e., not using these modifiers
        #flags
            #--min/--max 
                #removes lines with PI_HAT values below/above the given cutoff(s).
            #--ppc-gap
                #By default, the minimum distance between informative pairs of SNPs used in the pairwise population concordance (PPC) test is 500k base pairs; you can change this with the --ppc-gap flag.
                #see cluster analyses in stratificaction section for detailes about PPC
        #https://www.cog-genomics.org/plink/1.9/ibd

print_text("explore the kindship file looking for pairs with pi_hat > 0.2", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    echo 'See first lines of the IBD report'; \
    head ibd_report.genome; \
    echo '\n See now only those rows for which pi_hat > 0.2'; \
    awk \
        'BEGIN{FS=\" \"}{if((NR>1) && ($10 > 0.2)){print $0}}' \
        ibd_report.genome")
        #fields of plink.genome:
            #FID1    Family ID for first sample
            #IID1    Individual ID for first sample
            #FID2    Family ID for second sample
            #IID2    Individual ID for second sample
            #RT  Relationship type inferred from .fam/.ped file
                #FS = full sibling
                #HS = half-sibling
                #PO = parent-offspring
                #OT = other
                #https://groups.google.com/g/plink2-users/c/ckrUIamlxeE/m/RC0r8CXiBQAJ
            #EZ: IBD sharing expected value, based on just .fam/.ped relationship
            #Z0  P(IBD=0)
                #Individuals sharing close to zero alleles IBD at every locus (Z0~1) are unrelated.
            #Z1  P(IBD=1)
                #Individuals sharing one allele IBD at every locus (Z1~1) are parent-child pairs. 
                #On average, siblings share zero, one, and two alleles IBD at 25%, 50%, and 25% of the genome, respectively. Therefore, Z0=0.25, Z1=0.5 and Z2=0.25.
            #Z2  P(IBD=2)
                #(IBD). Individuals sharing two alleles IBD at nearly every locus (Z2~1) are monozygotic twins, or the pair is a single sample processed twice.
            #PI_HAT  Proportion IBD, i.e., P(IBD=2) + 0.5*P(IBD=1)
                #This counts the whole probability of having 2 shared alleles IBD and then half of the probability of having 1 shared allele IBD.
                #I guess this is the whole probability of have a shared allele at all.
                #particular cases
                    #Identical twins, and duplicates, are 100% identical by descent (Pihat 1.0)
                        #P(IBD=1)=0 and P(IBD=2)=1, thus PI_HAT=1+0.5*0=1
                        #Because they share two alelles in every loci, i.e., Z1=0 and Z2=1.
                    #First-degree relatives are 50% IBD (Pihat 0.5)
                        #Parent-child:
                            #P(IBD=1)=1 and P(IBD=2)=0, thus PI_HAT=0+0.5*1=0.5
                            #Becuase they share 1 allele in almost every locus, i.e., Z1=1 and Z2=0
                        #Siblings:
                            #P(IBD=1)=0.5 and P(IBD=2)=0.25, thus PI_HAT=0.25+0.5*0.5=0.5
                            #Because they share 1 allele in half of loci and 2 alleles in 25% of loci, i.e., Z1=0.5 and Z2=0.25
                    #Second-degree relatives are 25% IBD (Pihat 0.25)
                        #grandparents, grandchildren, aunts, uncles, nieces, nephews, and half-siblings
                    #Third-degree relatives are 12.5% equal IBD (Pihat 0.125).
                        #great-grandparents, great-grandchildren, first cousins, great-uncles, great-aunts, and the great-nieces and great-nephews
                        #https://www.biostars.org/p/58663/
                        #https://www.biostars.org/p/75335/
                    #A PI_HAT value of 0.2 is usually used as threshold so we remove second degree relatives and above. This is according to Marees et al., review. Also the trainability paper used that threshold.
                        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
                        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
            #PHE Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)
            #DST IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
            #PPC IBS binomial test
            #RATIO   HETHET : IBS0 SNP ratio (expected value 2)

print_text("convert delimiter of ibd report", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \
        ibd_report.genome > ibd_report_awk_processed.genome.tsv; \
    gzip \
        --force \
        ./ibd_report_awk_processed.genome.tsv; \
    ls -l")

print_text("load ibd report as pandas DF", header=4)
ibd_report = pd.read_csv( \
    "./data/genetic_data/quality_control/11_remove_related_samples/ibd_report_awk_processed.genome.tsv.gz", \
    sep="\t", 
    header=0, 
    low_memory=False)
print(ibd_report)

print_text("plot Z0 vs Z1 and Z2", header=4)
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2)
axs[0].scatter( \
    ibd_report["Z0"], \
    ibd_report["Z1"], \
    s=1)
axs[0].set_ylabel("Z1")
axs[1].scatter( \
    ibd_report["Z0"], \
    ibd_report["Z2"], \
    s=1)
axs[1].set_xlabel("Z0")
axs[1].set_ylabel("Z2")
plt.savefig( \
    fname="./data/genetic_data/quality_control/11_remove_related_samples/pairs_relatdness_before_filtering.png")
plt.tight_layout()
plt.close()
    #Our results:
        #We have many samples with Z0~1, i.e., unrelated individuals.
        #We do not have cases with Z1 close to 1, so to parent-child pairs, but we have cases with Z1=0.5, thus potential siblings (they share 1 allele in almost all loci) 
        #We have one pair with Z2=1 (they share 2 alleles in almost all loci), thus we have twins or the same sample duplicated.

print("check that the cases with Z0 and Z1 equal to 0 are those with Z2=1")
if(sum(ibd_report.loc[(ibd_report["Z0"]==0) & (ibd_report["Z1"]==0), "Z2"]!=1)!=0):
    raise ValueError("ERROR! FALSE! that the cases with Z0 and Z1 equal to 0 are NOT those with Z2=1")
else:
    print("TRUE")

print_text("plot histogram of pi_hat", header=4)
ibd_report["PI_HAT"].hist(bins=100)
    #PI_HAT  Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
        #This counts the whole probability of having 2 shared alleles IBD and then half of the probability of having 1 shared allele IBD.
        #see above for calculation in specific cases.
plt.yscale('log') 
    #to gain visibility, like in Ritche tutorial
        #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
plt.savefig( \
    fname="./data/genetic_data/quality_control/11_remove_related_samples/pi_hat_hist.png")
plt.close()
    #we have several cases above pi_hat of 0.2, i.e., at least second degree relatives, possibly including first-degree relatives (PI_HAT=0.5) and twins/duplicated samples (PI_HAT=1).

print_text("We have multiple cases with PI_HAT around 0.2, 0.5 and 1, which are potential second, first-degree relatives and twins/duplicated samples, respectively", header=4)
pi_hat_above_0_2=ibd_report.loc[ibd_report["PI_HAT"]>0.2, ["FID1", "IID1", "FID2", "IID2", "PI_HAT"]]
print(pi_hat_above_0_2)

print_text("print unique related samples with pi_hat > 0.2", header=4)
samples_pi_hat_above_0_2=pd.concat( \
    [pi_hat_above_0_2["IID1"], pi_hat_above_0_2["IID2"]], \
    axis=0).unique()
print("We have %s samples with PI_HAT > 0.2" % len(samples_pi_hat_above_0_2))
print(samples_pi_hat_above_0_2)


print_text("run KING-robust with plink2", header=3)
print_text("make subfolder for KING", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/11_remove_related_samples; \
    mkdir \
        --parents \
        ./king_analyses/; \
    ls -l")

print_text("intro to KING", header=4)
#we need to use other method to ensure I am not making a mistake so we can confidently say that these samples are related.
#We are going to KING, because, according to Christopher, is much less error prone. This method is also recommended in the Ricthie´s tutorial.
    #If you need to ask, then you should ignore PI_HAT. A major advantage of the KING method is that it is much harder to misuse.
        #https://www.biostars.org/p/434832/#434898

#PLINK 2 includes an implementation of the KING-robust [5] pairwise relatedness estimator, which can be used to prune all related pairs. This does not mean that both samples in each related pair are thrown out. Instead, –king-cutoff tries to keep as much data as possible, and as a consequence it usually keeps one sample out of each pair.
    #apply the KING-robust method
        #https://www.cog-genomics.org/plink/2.0/distance#make_king

#The relationship matrix computed by --make-rel/--make-grm-list/--make-grm-bin can be used to reliably identify close relations within a single population, if your MAFs are decent. However, Manichaikul et al.'s KING-robust estimator can also be mostly trusted on mixed-population datasets (with one uncommon exception noted below), and doesn't require MAFs at all. Therefore, we have added this computation to PLINK 2, and the relationship-based pruner is now based on KING-robust.
    #The exception is that KING-robust underestimates kinship when the parents are from very different populations. You may want to have some special handling of this case; --pca can help detect it.
    #We should not have parents here. We are using the second method in plink2, which is robust to mixed pops and MAF problems.

#Note that KING kinship coefficients are scaled such that duplicate samples have kinship 0.5, not 1 (pi_hat for twins is 1). First-degree relations (parent-child, full siblings; pi_hat for first-degree is 0.5) correspond to ~0.25, second-degree relations correspond to ~0.125 (pi_hat for second-degree is 0.25). Therefore, I understand that third-degree relatives will have 0.125/2=0.0625, because pi_hat of these relatives is 0.125.
    #see above for pi_hat calculations

#According to Christopher, we can consider king-cutoff as half of pi_hat if our purpose is just to remove relatives from first to third-degree
    #As long as you're just concerned with finding/pruning close relations (1st-2nd degree, maybe third degree), you can think of KING kinship as half of PI_HAT.
    #https://groups.google.com/g/plink2-users/c/z2HRffl-6k8/m/NL0pqc-7AQAJ

#It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) to screen for monozygotic twins and duplicate samples, ~0.177 (the geometric mean of 0.25 and 0.125) to add first-degree relations, etc. Therefore, I understand that in order to filter out second-degree relatives, we would need to calculate the geometric mean of 0.125 and 0.0625, as these are the king-cutoffs for second and third-degree relatives. This would make a cutoff value of 0.088. Indeed, Christopher specifically said this in the forum: "A threshold of 0.177 corresponds to removing 1st-degree relations, and 0.088 also removes 2nd-degree relations"
    #https://www.cog-genomics.org/plink/2.0/distance#king_coefs
    #https://groups.google.com/g/plink2-users/c/938B07i8AXQ/m/dMFI8-GLAwAJ
import numpy as np
from scipy.stats import gmean
print("The geometric mean of king-cutoffs for second and third-degree relatives (0.125 and 0.0625, respectively) is %s. This could be used to filter second-degree relatives" % gmean([0.125, 0.0625]))

#Important note about LD-prunning
#The KING documentation says "Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING."
    #https://www.kingrelatedness.com/manual.shtml
#Christopher also made a reference to this.
    #https://groups.google.com/g/plink2-users/c/NKBcu2cC160/m/3_IftZbEAwAJ
#The plink tutorial filter by relatedness in the whole dataset just after filtering samples and SNPs missingness
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
#So we are going to use the dataset without LD prunning and without excluding any chromosome. Although it seems that KING only consider autosomals anyways.

print_text("run --make-king-table to get a table with kindship coefficients", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf \
        --make-king-table \
        --king-table-filter 0.088 \
        --out ./11_remove_related_samples/king_analyses/king_table; \
    ls -l ./11_remove_related_samples")
    #--make-king
        #It seems that "--make-king" writes the KING-robust kindship coefficients into a matrix for all sample pairs, while --make-king-table creates a table. The removal of samples is done with "--king-cutoff", see below.
        #We select only those sample with a kindship coefficient greater than 0.088
            #This is usually the threshold used to consider sample above as related in at least second degree (see --king-cutoff below for details)
        #--make-king writes KING-robust coefficients in matrix form to plink2.king[.zst] or plink2.king.bin, while --make-king-table writes them in table form to plink2.kin0[.zst].
            #--make-king's matrix shape ('square', 'square0', 'triangle') and encoding ('bin', 'bin4') modifiers have the same behavior as --make-rel's, except that the diagonal is excluded from triangular output. Default shape+encoding is triangle+text, but note that bin/bin4 encoding changes the default shape to square.
                #The 'square', 'square0', and 'triangle' modifiers affect the shape of the output matrix. 'square' yields a symmetric matrix; 'triangle' (normally the default) yields a lower-trianglar matrix where the first row contains only the <sample #1-sample #1> relationship, the second row has the <sample #1-sample #2> and <sample #2-sample #2> relationships in that order, etc.; and 'square0' yields a square matrix with all cells in the upper right triangle zeroed out.
                #The 'bin' modifier causes the matrix to be written to plink2.rel.bin using little-endian IEEE-754 double encoding (suitable for loading from R). When using 'bin', the default output shape is 'square' instead of 'triangle'.
                #'bin4' uses IEEE-754 single-precision encoding, and is otherwise identical to 'bin'. This saves disk space, but you'll need to specify 4-byte single-precision input for your next analysis step. The following does so in R: readBin('<filename>', what="numeric", n=<number of entries>, size=4)
            #Only autosomes are included in this computation.
                #So I understand we do NOT have to filtered them out
            #Pedigree information is currently ignored; the between-family estimator is used for all pairs.
            #For multiallelic variants, REF allele counts are used.
            #--make-king jobs with the 'square0' or 'triangle' output shapes and all --make-king-table jobs can be subdivided with --parallel.
        #In addition, with --make-king-table,
            #The 'counts' modifier causes counts rather than 0.1 frequencies to be reported in the output columns that support both.
            #The 'rel-check' modifier causes only same-FID pairs to be reported. (The between-family KING estimator is still used.)
                #I understand these make to look only into samples of the same family, but we want to do comparisons between families, i.e., between batches, so we use the default. Indeed, I have found a correlation between one sample of one batch and one of another.
            #--king-table-filter causes only kinship coefficients ≥ the given threshold to be reported.
            #--king-table-subset causes only sample-pairs mentioned in the given .kin0 file (and optionally passing a kinship-coefficient threshold) to be processed. This allows you to start with a screening step which considers all sample pairs but only a small number of variants scattered across the genome (try --maf + --bp-space), and follow up with accurate kinship-coefficient computations for just the sample pairs identified as possible relations during the screening step. (This two-step approach remains practical with millions of samples!)
                #Not necessary for us as we have a few thousands
            #--king-table-require accepts one or more space/tab-delimited text files with sample IDs, and removes sample-pairs which don't contain at least one of the listed samples from the analysis. Similarly, --king-table-require-xor removes sample-pairs which don't contain exactly one of the listed samples.
            #Refer to the file format entry for other output details and optional columns. --make-king-table now covers much of PLINK 1.x --genome's functionality.
print_text("load the talbe in python")
king_table = pd.read_csv( \
    "./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/king_table.kin0", \
    sep="\t", 
    header=0, 
    low_memory=False)
print(king_table)
    #A text file with a header line, and one line per sample pair with kinship coefficient no smaller than the --king-table-filter value. When --king-table-filter is not specified, all sample pairs are included. The following columns are present:
        #FID1: FID of first sample in current pair
        #IID1: IID of first sample in current pair
        #FID2: FID of second sample in current pair
        #IID2: IID of second sample in current pair
        #NSNP: Number of variants considered (autosomal, neither call missing)
            #As this only consider autosomals, it is a lower number than the total SNPs that passed the filters considering also sex chromosome
            #According to the king table, the number of SNPs used for each sample is around 271K.
            #The total number of SNPs that passed the filter is 289148, while only 270790 autosomal and 511 pseudo-autosomals. If we sum 270790 plus 511, we get 271301, which is around the number of SNPs considered per sample in the king table. So it seems king is considering XY regions.
                #This is not a problem as we have ensured in a previous step that we have only XY variants in recognized XY regions for hg38.
        #HETHET: Proportion/count of considered call pairs which are het-het
        #IBS0: Proportion/count of considered call pairs which are opposite homs
        #KINSHIP: KING-robust between-family kinship estimate
print("extract the unique cases ")
samples_king_value_above_0_008=pd.concat( \
    [king_table["IID1"], king_table["IID2"]], \
    axis=0).unique()
print("We have %s samples with KING-robust between-family kinship estimate > 0.088" % len(samples_king_value_above_0_008))
print(samples_king_value_above_0_008)

print_text("do we have the same problematic samples in king (kindship>0.088) and pi_hat (pi_Hat>0.2)?", header=4)
print("check we have the exact same problematic IDs")
from natsort import natsorted
import numpy as np
check_pi_hat_king = np.array_equal( \
    natsorted(samples_king_value_above_0_008), \
    natsorted(samples_pi_hat_above_0_2) \
)
if(check_pi_hat_king):
    print(check_pi_hat_king)
else:
    raise ValueError("ERROR: FALSE! PI-HAT AND ING DO NOT SELECT THE SAME RELATED SAMPLES")
print("Count number of each case with both approaches")
print( \
    "Duplicates: PI_HAT=0.5 (%s) and KINSHIP>0.354 (%s)" % ( \
        ibd_report.loc[ibd_report["PI_HAT"]>=0.9,:].shape[0], \
        king_table.loc[king_table["KINSHIP"]>=0.354,:].shape[0]))
print( \
    "First-degree: PI_HAT around 0.5 (%s) and KINDSHIP between 0.17 and 0.35 (%s)" % ( \
        ibd_report.loc[(ibd_report["PI_HAT"]>0.4) & (ibd_report["PI_HAT"]<0.6),:].shape[0], \
        king_table.loc[(king_table["KINSHIP"]>=0.177) & (king_table["KINSHIP"]<0.354),:].shape[0]))
print( \
    "Second-degree: PI_HAT around 0.25 (%s) and KINDSHIP between 0.08 and 0.17 (%s)" % ( \
        ibd_report.loc[(ibd_report["PI_HAT"]>0.2) & (ibd_report["PI_HAT"]<0.4),:].shape[0], \
        king_table.loc[(king_table["KINSHIP"]>=0.088) & (king_table["KINSHIP"]<0.177),:].shape[0]))
print("check that NONE of the problematic cases have all NA for phenotypes, this would be a problem] because we would need to prioritize cases with pheno instead of using the greedy approach implement in PLINK 2")
#get the pheno data
#IMPORTANTE: SOME VARIABLES NEED PROCESING DUE TO ERRORS, BUT WE ARE NOT GOING TO DO IT HERE BECAUSE WE ONLY NEED THE IDS WITHOUT PHENO DATA
run_bash(" \
    cp \
        ./data/pheno_data/'Combat gene DNA GWAS 23062022_v2.xlsx' \
        ./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/pheno_data.xlsx")
import pandas as pd
import numpy as np
pheno_data = pd.read_excel(
    "./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/pheno_data.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data)
run_bash("rm ./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/pheno_data.xlsx")
#get IDs of those samples with NA for all phenotypes, which we are goin 
ids_no_pheno=pheno_data.loc[ \
    (pheno_data.drop("AGRF code", axis=1).isna().all(axis=1)) & \
    (pheno_data["AGRF code"].notna()), \
    "AGRF code"]
    #Extract the ID of those rows for which all columns are NA (except AGRF code) AND do not have NA for the ID
if(len(ids_no_pheno)!=41):
    raise ValueError("ERROR! FALSE! WE DO NOT HAVE THE CORRECT NUMBER OF SAMPLES WITH NA FOR ALL PHENOTYPES")
#check the problematic IDs do not include the mislabelled sampled and the samples with all NA for phenotypes
print(sum(ids_no_pheno.isin(["2397LDJA", "2399LDJA"]))==0)
print(sum(ids_no_pheno.isin(samples_king_value_above_0_008))==0)
print(sum(ids_no_pheno.isin(samples_pi_hat_above_0_2))==0)

print_text("Remove samples with a kindship above 0.354", header=4)
#My explanation
    #I am performing QC on a dataset with approximately 1400 samples and around 600K SNPs (hg38). I am using the KING-robust method implemented in PLINK 2 to detect and filter out second-degree relatives and above. As you explained in the documentation, I am using the geometric mean of the kinship coefficient of second and third-degree relatives (i.e., 0.088). The resulting KING table is attached.
    #As you can see, I have 15 pairs with a coefficient of 0.5, suggesting these are duplicates. Is this number common? or my KING results are strange? It seems a bit high to me, as it implies a significant number of errors during processing and/or genotyping the samples…
    #Please note that before this step, I removed duplicated SNPs and filtered SNPs by MAF (0.05), missingness (0.01), HWE (1e-25), and sample missingness (0.01). Around 270K SNPs passed all the filters. The KING table remains mostly the same even if I perform the analysis without these previous filters. Also, note that filtering by the usual pi_hat threshold for second-degree relatives (0.2) yields the exact same results. For pi_hat I used LD-prunned data (95K autosomal SNPs), but not for KING.
    #In case I could have an extremely abnormal proportion of twins in my data, I have checked that each pair of samples have different phenotypes values and, specifically, different ages, which is the case. Therefore, these are not twins. 
    #Finally, I have found that 4 of these samples have a mismatch between reported-biological sex and sex based on genetic data. In these cases, it is very clear that the genetic data of one sample has been duplicated and matched with another sample having a different self-reported sex.
    #All this leads me to think that I could actually have a high number of duplicated samples. If you think this is the case, what would be your recommendation for dealing with them?
    #The removal of samples with --king-cutoff essentially removes one sample in each pair, but I am unsure whether I should remove all these samples with a coefficient of 0.5. If two samples share two alleles in almost every loci, but have different phenotypes, which phenotype does this genotype belong to?
    #My point is that, in each pair, we have two different phenotypes and 1 genome, so we cannot be sure which is the correct phenotype that matches the genotype. I guess all the 30 samples should be removed. Does this make sense?
#Answer Christopher
    #Yes, I would remove both copies of any duplicate where you’re unsure about the phenotype.
    #A ~1% rate of sample-handling errors is pretty typical.

print("from the king table, select samples above the monozygotic twins threshold")
king_duplicates = king_table.loc[king_table["KINSHIP"]>=0.354,:]
#get the fam and sample IDs for each member of each pair, set the same column name and then concat
king_duplicates_ids = pd.concat( \
    [king_duplicates[["#FID1", "IID1"]].rename(columns={"#FID1": "FID", "IID1": "IID"}), \
    king_duplicates[["FID2", "IID2"]].rename(columns={"FID2": "FID", "IID2": "IID"})], \
    axis=0)
#select those unique fam-sample IDs
king_duplicates_ids=king_duplicates_ids.loc[~king_duplicates_ids.duplicated(["FID", "IID"]), :]
#check
if((king_duplicates_ids.shape[0]!=30) | (king_duplicates.shape[0]!=15)):
    raise ValueError("ERROR! FALSE! WE DO NOT HAVE THE EXPECTED NUMBER OF DUPLICATED SAMPLES")
#save the ids to filter them out
king_duplicates_ids.to_csv("./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/id_samples_duplicates.txt",
    sep="\t",
    header=False,
    index=False)

print("remove duplicated samples from plink files")
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf \
        --remove ./11_remove_related_samples/king_analyses/id_samples_duplicates.txt \
        --make-bed \
        --out ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates"
)
    #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.
    #we cannot use --king-cutoff because we want to remove ALL samples implicated as we do not know which phenotype should associate with the duplicated genotype (see question to christopher). --king-cutoff removes one of each pair.

print("load the resulting FAM file and check we do not have the duplicated samples")
fam_file_no_dup_samples = pd.read_csv( \
    "./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates.fam", \
    sep=" ", \
    header=None, \
    low_memory=False, \
    names=["FID", "IID", "fatherID", "mother_ID", "sex_code", "phenotye_value"])
fam_file_no_dup_samples
#check
merged = fam_file_no_dup_samples.merge( \
    king_duplicates_ids, \
    on=["FID", "IID"], \
    how="left", \
    indicator=True)
    #To check whether the combination of two columns in a pandas DataFrame is present in another DataFrame, you can use the merge function with the indicator parameter set to True
    #In this code, fam_file_no_dup_samples.merge(king_duplicates_ids, on=["FID", "IID"], how="left", indicator=True) merges df1 and df2 on 'column1' and 'column2'. The indicator parameter adds a column to the resulting DataFrame that indicates whether each row is present in df1 only, df2 only, or both. The result is then assigned to the merged variable.
    #how="left":
        #In this case, how='left' means it's a left join. A left join returns all the rows from the left DataFrame and the matched rows from the right DataFrame. If there is no match, the result is NaN on the right side.
        #therefore, we are selecting all samples included in the new fam file.
    #on=["FID", "IID"] means that the merge is being performed on the 'FID' and 'IID' columns. That is, pandas is looking for rows in fam_file_no_dup_samples and king_duplicates_ids where the 'FID' and 'IID' values are the same, and combining those matching rows into a single row in the resulting DataFrame.
if(sum(merged["_merge"] == "both")!=0):
    #if a row is present in the fam file and in the list of duplicates, it will have a value of "both". We should not have any case like this.
    raise ValueError("ERROR! FALSE! WE STILL HAVE DUPLICATED SAMPLES IN OUR FILTERED DATASET")
else:
    print("True")

print("make a new king table to see the result consdering third-degree relatives and above")
#make the table in plink2
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates \
        --make-king-table \
        --king-table-filter 0.088 \
        --out ./11_remove_related_samples/king_analyses/king_table_after_dups; \
    ls -l ./11_remove_related_samples")
#load the table
king_table_after_dups = pd.read_csv( \
    "./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/king_table_after_dups.kin0", \
    sep="\t", 
    header=0, 
    low_memory=False)
print(king_table_after_dups)
#check
if( \
    (king_table_after_dups.shape[0]!=5) | \
    (sum(king_table_after_dups["KINSHIP"]>=0.354)!=0) | \
    (sum(king_table_after_dups["KINSHIP"]<=0.088)!=0)):
    raise ValueError("ERROR! FALSE! WE DO NOT HAVE EXACTLTY 5 RELATED SAMPLES AFTER APPLYING THE INITIAL FILTER OF TWINS OR WE HAVE SAMPLES WITH KINSHIP ABOVE 0.354, I.E., DUPLICATES, OR BELOW <0.088, I.E., BELOW THE THRESHOLD USED FOR THE KING TABLE")

print_text("Remove samples with first and second-degree relationship using --king-cutoff", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates \
        --king-cutoff 0.088 \
        --make-bed \
        --out ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related; \
    ls -l ./11_remove_related_samples/king_analyses")
    #code taken from:
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
    #--king-cutoff
        #relationship-based pruner. This uses KING-Robust:
        #If used in conjunction with a later calculation (I guess --make-bed suffices), --king-cutoff excludes one member of each pair of samples with kinship coefficient greater than the given threshold. Alternatively, you can invoke this on its own to write a pruned list of sample IDs to plink2.king.cutoff.in.id, and excluded IDs to plink2.king.cutoff.out.id. We are using the former option by creating a bed file.
        #According to Maares et al, it is frequent to remove second-degree relatives and above by using a pi_hat=0.2. Indeed, this is the threshold used in the trainability paper. This would correspond with a king-cutoff of 0.1 (assuming the king-cutoff is half of pi_hat, see above)
        #However, according to plink docs, we would need a king-threshold of 0.088 to remove second degree samples (see above).
        #We are going to follow the plink documentation and use the more stringent threshold.
    #PLINK tries to maximize the final sample size, but this maximum independent set problem is NP-hard, so we use a greedy algorithm which does not guarantee an optimal result as opposed to just removing one of each pair of related individuals. In practice, --king-cutoff does yield a maximum set whenever there aren't too many intertwined close relations, but if you want to try to beat it (or optimize a fancier function that takes the exact kinship-coefficient values into account), use the --make-king[-table] and --keep/--remove flags and patch your preferred algorithm in between.
        #From this, I understand that we are indeed doing better than just removing a sample from each pair, but getting a the set the maximum number of samples possibles while stying below the KING cutoff (at least if there are no too much intertwined connections). So instead just removing the member with the lowest genotyping rate, we use this.
        #so for this we can use --king-table. Calculate kindship coefficients with --make-king table, then use these coefficients as input in a algorithm that remove related samples, and the list of remaining samples being used as filter in plink with --keep/--remove.
    #--king-cutoff usually computes kinship coefficients from scratch. However, you can provide a precomputed kinship-coefficient matrix (must be --make-king binary format, triangular shape, either precision ok) as input to --king-cutoff, or a .kin0 file (can be compressed) to --king-cutoff-table; this is a time-saver when experimenting with different thresholds.
        #so if you want to compare multiple thresholds, you can save the pre-computation and save time.

print("see the samples removed in the last filter")
king_removed_samples_last_filter = pd.read_csv( \
    "./data/genetic_data/quality_control/11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.king.cutoff.out.id", \
    sep="\t", \
    header=0, \
    low_memory=False)
print(king_removed_samples_last_filter)

print_text("compare the number of removed samples with the total number of problematic samples", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    n_samples_before=$( \
        awk \
        'BEGIN{FS=\"\t\"; OFS=\"\t\"}END{print NR}' \
        ./10_remove_low_maf_snps/merged_batches_hwe_par_callrate_LogRDev_maf.fam \
    ); \
    n_samples_after=$( \
        awk \
        'BEGIN{FS=\"\t\"; OFS=\"\t\"}END{print NR}' \
        ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.fam \
    ); \
    samples_lost=$(($n_samples_before-$n_samples_after)); \
    if [[ $samples_lost -eq 35 ]]; then \
        printf 'The number of samples filtered out due to relatedness is: %s' \"$samples_lost\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 20 samples due to relatedness';\
    fi")
    #We should have 35 samples removed:
        #15 pairs of samples with kinship>0.354 from which all samples were removed, i.e., 15*2=30.
        #5 pairs of samples with kinship between 0.088 and 0.354, from which one sample per pair was removed.

#Final note: Some of the sex incosistences initially detected were indeed duplicated cases, i.e., two samples have the exact same genome, hence both have the same biological sex based on genetics, but in the phenotype one of them has a different sex:
    #combat_ILGSA24-17873	8145ORJJ	combat_ILGSA24-17873	7974SJNN	271094	0.326426	0	0.5
    #combat_ILGSA24-17873	8244DCDJ	combat_ILGSA24-17873	8244CDSJ	271170	0.327887	0	0.499992
    #combat_ILGSA24-17873	8501CGDJ	combat_ILGSA24-17873	8500NHNJ	271207	0.325301	0	0.5
    #combat_ILGSA24-17873	9187MFMM	combat_ILGSA24-17873	0197SNMM	271149	0.331242	0	0.5
    #We should not rescue these 4 cases becasue we cannot completely sure that the case with correct sex (i.e., same sex according to genetics and phenotype) is the one we have to select. We cannot be completely sure.

# endregion





#########################
# region MAF loop #######
#########################
print_text("apply again the MAF and missing filters", header=2)
#We repeat in two steps the MAF-missing snps filters and then sample filter to check that the previous removal of samples did not change allele frequencies in a way that after applying MAF + missing filters again, we lose more samples.
#In other words, 
    #imagine you have already filtered by MAF and missing (SNP and sample). Of course, after the SNP filters, you can be sure that no SNP is below the MAF threshold or above the missing threshold. 
    #But then, you remove samples. It could be the case that after removing samples, the MAF of a SNP gets below 0.05 because the few minor carriers have been removed. 
    #This makes snps with MAF close but above of 0.05 actually getting below that thershold. We remove again these SNPs.
    #This removal of SNPs can in turn influence the missing call percentage of samples because maybe a sample has a call rate of 0.99, but we have removed one SNP for which it had data (i.e., non-missing). Its call rate has decreased being now below the threshold. See the following dummy example:
        #We only have 5 SNPs, thus 10 potential genotypess
        #A sample has missing for 2 genotypes, thus its missing rate is 2/10=0.2 and its call rate is 8/10=0.8.
        #If we remove one of the SNPs, the number of potential genotypes is 8 instead of 10.
        #Imagine that the removed SNPs was indeed one of the SNPs without missing for this sample, thus the number of genotypes for this sample goes down to 6.
        #Therefore, its missing rate is now 2/8=0.25 and its call rate is 6/8=0.75.
        #We have now a higher missing rate and a lower call rate.
#We have to continue until we remove SNPs by filters and then no additional sample is removed. In that moment, we can be sure that all SNPs and samples meet the filters considered.
    #"An iterative procedure that repeats the SNP and sample-level filtering until no additional samples are removed is also common" 
        #https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6605
#the other previous steps should not be influenced by this, as we just removed SNPs duplicated or with wrong chromosome names
print_text("create a folder to save plink data after new filters", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./12_loop_maf_missing; \
    ls -l")



print_text("create missing and freq reports", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related \
        --freq \
        --missing \
        --out ./12_loop_maf_missing/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related; \
    ls -l")



print_text("check if we have SNPs below the MAF threshold after removing samples", header=3)
n_snps_below_maf_first_round_raw = run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    awk \
        'BEGIN{FS=\" \"}{if((NR>1) && ($5 < 0.05)){count++}}END{print count}' \
        ./12_loop_maf_missing/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.frq", return_value=True).strip()
n_snps_below_maf_first_round = 0 if n_snps_below_maf_first_round_raw=="" else int(n_snps_below_maf_first_round_raw)
print(f"We have {n_snps_below_maf_first_round} SNPs below the MAF threshold after the first round of filters")


print_text("check if we have SNPs below the missing threshold after removing samples", header=3)
n_snps_below_missing_first_round_raw = run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    awk \
        'BEGIN{FS=\" \"}{if((NR>1) && ($5 > 0.01)){count++}}END{print count}' \
        12_loop_maf_missing/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.lmiss", return_value=True).strip()
n_snps_below_missing_first_round = 0 if n_snps_below_missing_first_round_raw=="" else int(n_snps_below_missing_first_round_raw)
print(f"We have {n_snps_below_missing_first_round} SNPs above the missing threshold after the first round of filters")



print_text("repeat the filters if these snps exists", header=3)
if(n_snps_below_maf_first_round>0) | (n_snps_below_missing_first_round>0):
    print_text("apply MAF and missing SNP filters", header=4)
    run_bash(" \
      cd ./data/genetic_data/quality_control/; \
      plink \
          --bfile ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related \
          --maf 0.05 \
          --geno 0.01 \
          --make-bed \
          --out ./12_loop_maf_missing/loop_maf_missing_1; \
      ls -l ./12_loop_maf_missing/")  

    print_text("apply missing sample filters", header=4)
    run_bash(" \
        cd ./data/genetic_data/quality_control/; \
        plink \
            --bfile ./12_loop_maf_missing/loop_maf_missing_1 \
            --mind 0.01 \
            --make-bed \
            --out ./12_loop_maf_missing/loop_maf_missing_2; \
        ls -l ./12_loop_maf_missing/")

    print_text("create MAF and missing reports", header=4)
    run_bash(
        "cd ./data/genetic_data/quality_control/12_loop_maf_missing/; \
        plink \
            --bfile loop_maf_missing_2 \
            --freq \
            --missing \
            --out loop_maf_missing_2_reports; \
        ls -l")
    
    print_text("check reports again", header=4)
    n_snps_below_maf_second_round_raw = run_bash(" \
        cd ./data/genetic_data/quality_control/12_loop_maf_missing/; \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 < 0.05)){count++}}END{print count}' \
            loop_maf_missing_2_reports.frq", return_value=True).strip()
    n_snps_below_maf_second_round = 0 if n_snps_below_maf_second_round_raw=="" else int(n_snps_below_maf_second_round_raw)
    print(f"We have {n_snps_below_maf_second_round} SNPs below the MAF threshold after the second round of filters")    
    n_snps_below_missing_second_round_raw = run_bash(" \
        cd ./data/genetic_data/quality_control/12_loop_maf_missing/; \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 > 0.01)){count++}}END{print count}' \
            loop_maf_missing_2_reports.lmiss", return_value=True).strip()
    n_snps_below_missing_second_round = 0 if n_snps_below_missing_second_round_raw=="" else int(n_snps_below_missing_second_round_raw)
    print(f"We have {n_snps_below_missing_second_round} SNPs above the missing threshold after the second round of filters")
    n_samples_below_missing_second_round_raw = run_bash(" \
        cd ./data/genetic_data/quality_control/12_loop_maf_missing/; \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($6 > 0.01)){count++}}END{print count}' \
            loop_maf_missing_2_reports.imiss", return_value=True).strip()
    n_samples_below_missing_second_round = 0 if n_samples_below_missing_second_round_raw=="" else int(n_samples_below_missing_second_round_raw)
    print(f"We have {n_samples_below_missing_second_round} Samples above the missing threshold after the second round of filters")

    print_text("stop if we still have SNPs/samples not meeting the filters", header=4)
    if (n_snps_below_maf_second_round>0) | (n_snps_below_missing_second_round>0) | (n_samples_below_missing_second_round>0):
        raise ValueError("ERROR: FALSE! WE HAVE AN ERROR WITH THE ITERATIONS OF MAF/MISSING FILFERS!, WE STILL HAVE SNPS OR SAMPLES NOT MEETING THE CONDITIONS")
else:
    print("Step not required")



print_text("just copy plink files if no filter is required", header=3)
if(n_snps_below_maf_first_round==0) & (n_snps_below_missing_first_round==0):
    run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    cp \
        ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.bim \
        ./12_loop_maf_missing/loop_maf_missing_2.bim; \
    cp \
        ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.fam \
        ./12_loop_maf_missing/loop_maf_missing_2.fam; \
    cp \
        ./11_remove_related_samples/king_analyses/merged_batches_hwe_par_callrate_LogRDev_maf_duplicates_related.bed \
        ./12_loop_maf_missing/loop_maf_missing_2.bed; \
    ls -l ./12_loop_maf_missing/")
else:
    print("Step not required.")

# endregion




####################################
# region CHECK HETERO ##############
####################################

#After applying all filters except pop stratification and sex-mismatch, we are going to take a look to heterozigosity rate.

#It is valuable to ensure the data are free from potential inbreeding. To check for these confounders, it is necessary to carry out methods that will adjust the data based on established thresholds by first calculating the genotype call rate and heterozygosity rate. The heterozygosity rate is the ratio of the total non-missing genotypes minus the homozygous genotypes (i.e., heterozygous genotypes) compared to the total non-missing genotypes in a given sample. The --het flag in PLINK produces a file with a list of heterozygous haploid genotypes (*.hh) and a file with individual heterozygosity information (*.het). By using the *.het file, we can calculate the heterozygosity rate where the O.HOM column represents the observed homozygous genotypes while the N.NM.
    #page 6 in Ritchie´s review.
print_text("visualize heterozygosity", header=2)
print_text("create folder for this step", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./13_visual_hetero/; \
    ls -l")



print_text("calculate heterozygosity and missing reports per sample", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./12_loop_maf_missing/loop_maf_missing_2 \
        --het \
        --missing \
         --out ./13_visual_hetero/loop_maf_missing_2_hetero \
    ")



print_text("process the missing report", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    awk \
        'BEGIN{OFS=\" \"; OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6}' \
        ./13_visual_hetero/loop_maf_missing_2_hetero.imiss > ./13_visual_hetero/loop_maf_missing_2_hetero_awk_processed.imiss \
    ")
missing_data=pd.read_csv("./data/genetic_data/quality_control/13_visual_hetero/loop_maf_missing_2_hetero_awk_processed.imiss", sep="\t", low_memory=False)
print(missing_data)



print_text("process the heterozygosity report", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    awk \
        'BEGIN{OFS=\" \"; OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6}' \
        ./13_visual_hetero/loop_maf_missing_2_hetero.het > ./13_visual_hetero/loop_maf_missing_2_hetero_awk_processed.het \
    ")
hetero_data=pd.read_csv("./data/genetic_data/quality_control/13_visual_hetero/loop_maf_missing_2_hetero_awk_processed.het", sep="\t", low_memory=False)
print(hetero_data)



print_text("merge both datasets", header=3)
merged_hetero = missing_data.merge( \
    hetero_data,  \
    how="inner",  \
    on=["FID", "IID"])
    #we need to combine in the same row columns that belong to the same family ID and sample ID
    #we need only those cases with missing heterozygosity data, so we need the intersection of both dataframes
print(merged_hetero)



print_text("calculate heterozygosity rate", header=3)
#The heterozygosity rate is the ratio of the total non-missing genotypes minus the homozygous genotypes (i.e., heterozygous genotypes) compared to the total non-missing genotypes in a given sample. By using the *.het file, we can calculate the heterozygosity rate where the O.HOM column represents the observed homozygous genotypes while the N.NM column represents the non-missing genotypes
merged_hetero["hetero_rate"] = (merged_hetero["N(NM)"] - merged_hetero["O(HOM)"]) / merged_hetero["N(NM)"]
    #formula also obtained from Marees et al.



print_text("plot missingness against heterozygosity rate like in page 6 of Ritchie´s review", header=3)
#calculate heterozygosity limits
hetero_high_threshold = merged_hetero["hetero_rate"].mean()+(3*merged_hetero["hetero_rate"].std())
hetero_low_threshold = merged_hetero["hetero_rate"].mean()-(3*merged_hetero["hetero_rate"].std())
import matplotlib.pyplot as plt
#proportion of missing genotypes vs. heterozygosity rate
plt.scatter(merged_hetero["F_MISS"], merged_hetero["hetero_rate"], s=0.5)
#add horizontal lines at mean heterozygosity + and - 3 standard deviations
plt.axhline( \
    y=hetero_high_threshold, \
    color='red', \
    linestyle='dashed')
plt.axhline( \
    y=hetero_low_threshold,  \
    color='red',  \
    linestyle='dashed')
#add vertical lines at mean missingness in addition to mean +- 5% of mean missingness
plt.axvline( \
    x=merged_hetero["F_MISS"].mean(), \
    color='red', \
    linestyle='dashed')
plt.axvline( \
    x=merged_hetero["F_MISS"].mean()+(merged_hetero["F_MISS"].mean()*0.05), \
    color='blue', \
    linestyle='solid')
plt.axvline( \
    x=merged_hetero["F_MISS"].mean()-(merged_hetero["F_MISS"].mean()*0.05), \
    color='black', \
    linestyle='solid')
#add axes labels and title
plt.xlabel('F_MISS')
plt.ylabel('hetero_rate')
plt.title('F_MISS vs hetero_rate')
#save figure and close
plt.savefig( \
    fname="./data/genetic_data/quality_control/13_visual_hetero/scatter_hetero_rate_vs_missing.png")
plt.close()
    #Based on these results, it appears that the quality of the data is very clean
        #there is a low level of missingness. Low genotype call rates (3% to 7%) suggest low sample quality. Yes, we have filtered out samples with missingness above 0.01, but these were only 10. 
        #The majority of heterozygosity remains within three standard deviations of the mean (horizontal red lines indicate mean ± 3 std). Less than expected heterozygosity (mean – 3 SD) suggests possible inbreeding and greater than expected heterozygosity (mean + 3 SD) suggests possible sample contamination. However, these thresholds should take into account the expected heterozygosity rates in the ancestry group under study, as some diverse populations may exhibit different rates of heterozygosity than other populations. 
            #In our case, just 13 samples are above the mean hetero plus 3 standard deviations, while 42 are below mean hetero less 3 standard deviations. It seems we are ok because in both cases, the outliers are close to the limits (see below). At least closer compared to the plot of Ritchie´s review. In that figure you can see the limits are 0.16-0.18 but there are samples at 0.22 and 0.12.
            #We should not have relatedness as we have already filtered that using KING-robust.
            #Also note we have different ancestries, so this could cause the small deviation from the normal levels of heterozygosity.



print_text("see the number of samples above the mean heterozigosity +- 3 standard deviations", header=3)
samples_above_hetero = merged_hetero.loc[merged_hetero["hetero_rate"] > hetero_high_threshold,:]
print(samples_above_hetero)
print(samples_above_hetero.shape[0])
print("the 97.5% percentile of the heterozigosity for samples above the heterozygosity upper limit correspond with the following percentage with respect to that limit") 
print(((samples_above_hetero["hetero_rate"].quantile(q=0.975) - hetero_high_threshold)/hetero_high_threshold)*100)



print_text("see the number of samples below the mean heterozigosity +- 3 standard deviations", header=3)
samples_below_hetero = merged_hetero.loc[merged_hetero["hetero_rate"] < hetero_low_threshold,:]
print(samples_below_hetero)
print(samples_below_hetero.shape[0])
print("the 2.5% percentile of the heterozigosity for samples below the heterozygosity lower limit correspond with the following percentage with respect to that limit") 
print(((hetero_low_threshold-samples_below_hetero["hetero_rate"].quantile(q=0.025))/hetero_low_threshold)*100)

# endregion






###################################
# region POP STRATIFICACION #######
###################################
print_text("start pop stratificaction analyses to detect individuals that are outliers", header=2)
print_text("start with pca", header=3)
print_text("open folder", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./14_pop_strat/01_pca; \
    ls -l")


print_text("Filter out snps for the PCA", header=4)
#Besides filtering by MAF and selecting high quality SNPs, we need to do LD prunning. In addition, the marees review and Ritchie tutorial remove sex chromosome
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./12_loop_maf_missing/loop_maf_missing_2 \
        --indep-pairwise 500kb 1 0.2 \
        --out ./14_pop_strat/01_pca/ldpruned_pca; \
    ls -l ./14_pop_strat/01_pca/")
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./12_loop_maf_missing/loop_maf_missing_2 \
        --extract ./14_pop_strat/01_pca/ldpruned_pca.prune.in \
        --autosome \
        --make-bed \
        --out ./14_pop_strat/01_pca/loop_maf_missing_2_ldprunned_autosome_pca;\
    ls -l ./11_remove_related_samples/")
#look section of sample relatedness for details about the filters applied in this step

print_text("checks after filtering", header=4)
#In the Ritchie's tutorial, they ended up with 67,000 autosomal variants in linkage equilibrium in order to calculate IBD and pi_hat. I guess they used the same dataset for the PCA. They also say that "It is recommended that the user prune to approximately ∼100,000 SNPs for use in PCA analyses". The github of the paper shows the removal of sex chromosomes (step 8) before the PCA (step 9), while Marees et al explictily says that we should perform the PCA on autosomal SNPs. Althougl plink tutorial does not filter out sex chromosomes before, we are doing it given what the other sources are doing.
    #14_pop_strat/01_pca
#In the admixture tutorial and Admixture docs they use a more stringent LD filtering, but they say that if you want to control for population structure between related populations (i.e., within the same continent) then 10K is not enough and you should use 100K. 
    #Using the approach suggested in the Admixture manual and tutorial gives 58K SNPs.
        #--indep-pairwise 50 10 0.1
        #http://dalexander.github.io/admixture/admixture-manual.pdf
#Therefore, we need around 100K independent autosomal SNPs. This seems to be enough considering all scenarios.
#It also seems important to remove related samples before the PCA
    #"We exclude related individuals as we are not interested in population structure that is due to family relationships and methods such as PCA and ADMIXTURE can inadvertently mistake family structure for population structure"
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/; \
    auto_ld_snps_in=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./loop_maf_missing_2_ldprunned_autosome_pca.bim); \
    printf 'The number of included autosomal SNPs in linkage equilibrium is %s\n' \"$auto_ld_snps_in\"; \
    echo 'Is this number greater than 93K?'; \
    if [[ $auto_ld_snps_in -gt 93000 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #you can print a variable with text using printf and %s for "string". Then you call the variable you want to add within "". You could use two variables: "$var1 $var2"
            #https://phoenixnap.com/kb/bash-printf
print("All non-autosomals SNPs have been removed from the LD pruned dataset?")
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca; \
    n_non_auto_snps=$( \
        awk \
            'BEGIN{FS=\" \"}{if($1==0 || $1==23 || $1==24 || $1==25 || $1==26 || $1== \"X\"|| $1==\"Y\" || $1==\"XY\" || $1==\"MT\"){count++}}END{print count}'\
            ./loop_maf_missing_2_ldprunned_autosome_pca.bim); \
    if [[ $n_non_auto_snps -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")


print_text("run PCA to detect clusters of samples", header=3)
#Once you have LD-pruned and MAF-filtered your dataset, PLINK 2’s –pca command has a good shot of revealing large-scale population structure. For example,
print_text("make the run", header=4)
run_bash("\
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca; \
    plink2 \
        --pca 20 \
        --bfile ./loop_maf_missing_2_ldprunned_autosome_pca \
        --out ./pca_results")
            #--pca extracts top principal components from the variance-standardized relationship matrix computed by --make-rel/--make-grm-{bin,list}. This writes a tab-delimited table to pca_results.eigenvec, with one sample per row and one principal component per later column
                #https://www.cog-genomics.org/plink/2.0/strat#pca
                #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
            #The main plink2.eigenvec output file can be read by --covar, and can be used to correct for population stratification in --glm regressions...
                #...assuming that the top principal components in your genomic dataset actually reflect broad population structure, rather than genotyping/sequencing-batch-related error patterns, small-scale family structure or sample duplication, crazy outliers... The .eigenvec file can be easily loaded and plotted in R; this should help you find significant batch effects and outliers. --king-cutoff removes duplicate samples and close relations.
                #Since this is based on the relationship matrix, it is critical to remove very-low-MAF variants before performing this computation.
                    #DONE
                #LD pruning (using e.g. --indep-pairwise) reduces the risk of getting PCs based on just a few genomic regions, and tends to prevent deflation of --glm test statistics.
                    #DONE.
            #Technical details
                #By default, 10 PCs are extracted; you can adjust this by passing a numeric parameter.
                #This was reduced from PLINK 1.9's default of 20, since (i) the randomized algorithm would otherwise require ~4x as much memory, and (ii) in practice, 10 PCs has been effective across a wide range of studies.
                #The 'approx' modifier causes the standard deterministic computation to be replaced with the randomized algorithm originally implemented for Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ, Price AL (2016) Fast Principal-Component Analysis Reveals Convergent Evolution of ADH1B in Europe and East Asia. This can be a good idea when you have >5000 samples, and is almost required once you have >50000.
                    #NOT OUR CASE
                #The randomized algorithm always mean-imputes missing genotype calls. For comparison purposes, you can use the 'meanimpute' modifier to request this behavior for the standard computation.
                #'scols=' can be used to customize how sample IDs appear in the .eigenvec file. (maybefid, fid, maybesid, and sid column sets are supported; the default is maybefid,maybesid.)
                #The 'allele-wts' modifier requests an additional one-line-per-allele .eigenvec.allele file with PCs expressed as allele weights instead of sample scores. When it's present, 'vzs' causes the .eigenvec.allele file to be Zstd-compressed. 'vcols=' can be used to customize the .eigenvec.allele report columns; refer to the file format entry for details.
                    #WE ARE INTERESTED IN SAMPLES
                #If all your variants are biallelic, you can instead use the 'biallelic-var-wts' modifier to request the old .eigenvec.var format instead.
                #Given an allele-weight or variant-weight file, you can now use --score for PCA projection. This replaces PLINK 1.9's --pca-clusters/--pca-cluster-names projection flags
            #You may also want to look at EIGENSOFT 7, which has additional features like automatic outlier removal, LD regression, and Tracy-Widom significance testing of PCs
                #https://www.hsph.harvard.edu/alkes-price/software/


print_text("load the eigenvec file generated", header=4)
pca_results = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_results.eigenvec", \
    sep="\t", \
    header=0, \
    low_memory=False)
pca_results
        #Produced by --pca. Accompanied by an .eigenval file, which contains one eigenvalue per line.
        #The .eigenvec file 
            #The .eigenvec file is a tab-delimited file with a header line and between 1+V and 3+V columns per sample, where V is the number of requested principal components. Each row is a sample. The first columns contain the sample ID, and the rest are principal component scores in the same order as the .eigenval values (with column headers 'PC1', 'PC2', ...). 
            #With the 'allele-wts' modifier, an .eigenvec.allele file is also generated. It's a text file with a header line, followed by one line per allele with several columns (see link)
            #Alternatively, with the 'biallelic-var-wts' modifier, an old-style .eigenvec.var file is generated. It's a text file with a header line, followed by one line per variant with several columns (see link)
        #https://www.cog-genomics.org/plink/2.0/formats#eigenvec


print_text("batch effects", header=4)
#You will usually want to sanity-check the output at this point, and verify that the top principal components do not correlate too strongly with, e.g., sequencing facility or date. (A full discussion of “batch effects” and how to deal with them could take up an entire chapter; worst case, you may have to analyze your batches separately, or even redo all genotyping/sequencing from scratch. I will be optimistic here and suppose that no major problem was uncovered by PCA, but be aware that this is frequently your best chance to catch data problems that would otherwise sink your entire analysis.)

print("plot the PCA values between batches")
import seaborn as sns
#Melt the DataFrame to long format
pca_melted = pd.melt( \
    pca_results.drop("IID", axis=1),  \
    id_vars="#FID",  \
    var_name="PCA", \
    value_name="value")
    #pca_results.drop("IID", axis=1): This part of the code is using the drop function from pandas to remove the "IID" column from the pca_results DataFrame. The axis=1 parameter specifies that "IID" is a column name.
    #pd.melt(...): This is using the melt function from pandas to reshape the DataFrame from wide format to long format.
    #id_vars='#FID': This parameter specifies that the '#FID' column should be used as the identifier variable. The values in this column will stay the same during the reshaping.
    #var_name='PCA': This parameter specifies that the new column created from the reshaping that contains the column headers from the original DataFrame should be named 'PCA'.
    #value_name='value': This parameter specifies that the new column created from the reshaping that contains the values from the original DataFrame should be named 'value'.
    #So, overall, this line of code is reshaping the pca_results DataFrame (after dropping the "IID" column) from wide format to long format, with '#FID' as the identifier variable, 'PCA' as the variable column, and 'value' as the value column. Therefore, we have a row per sample and PCA value, being first the values of PC1 for batch_1

#check the melting
expected_melted_rows =  \
    (sum(pca_results["#FID"] == "combat_ILGSA24-17303")*(pca_results.shape[1]-2)) + \
    (sum(pca_results["#FID"] == "combat_ILGSA24-17873")*(pca_results.shape[1]-2))
    #in the melted DF we should have the number of samples of each batch multiplied by the number of columns miuns FID and IID (i.e., the PCs)
if(pca_melted.shape[0] != expected_melted_rows):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM PREPARING THE DATASET TO MAKE THE BOXPLOTS OF PCAs vs. batches")

#create the boxplot
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
sns.boxplot(x='PCA', y='value', hue='#FID', data=pca_melted)
plt.savefig( \
    fname="./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_vs_batches.png")
plt.close()
    #plt.figure(figsize=(10, 6)) creates a new figure with a specified size.
    #creates a boxplot using seaborn.
        #A box plot (or box-and-whisker plot) shows the distribution of quantitative data. The box shows the quartiles of the dataset while the whiskers extend to show the rest of the distribution, except for points that are determined to be "outliers" using a method that is a function of the inter-quartile range.
        #x, y: These are the names of variables in data or vector data. One of them is the categorical axis (usually x; the name of the PC), and the other is the quantitative axis (usually y; the value of the PC). If only one is given, the other axis will have the boxplot for each level in the data
        #hue: This is a name of a variable in data or a vector data (the batch in our case). It is used to group the data and produce boxes with different colors.
    #output
        #Box: The box itself represents the interquartile range (IQR). The bottom of the box indicates the first quartile (25th percentile), and the top of the box indicates the third quartile (75th percentile). Therefore, the box spans the IQR.
        #Line in the box: The line inside the box represents the median (50th percentile) of the data.
        #Whiskers: The lines extending vertically from the box (i.e., the "whiskers") indicate variability outside the upper and lower quartiles. They extend to the furthest data point within 1.5 * IQR from the box.
        #Points beyond the whiskers: Any points beyond the whiskers can be considered outliers, i.e., unusually high or low values.
    #For all PCs, the median and IQR are virtually identical between batches. In the case of outliers, most of them are in the same range between batches, except the most extreme. In other words, the ranges with the highest density of outliers tend to overlap between batches. Note that these extreme outliers could be caused by ancestry differences.

print("check the distribution of the PCs between batches as it should be similar (Mann–Whitney U test assumption)")
#define pc columns and batches
pc_columns = pca_results.drop(["#FID", "IID"], axis=1).columns
fid1 = 'combat_ILGSA24-17303'
fid2 = 'combat_ILGSA24-17873'
#run loop across PCs saving the result in a PDF
from matplotlib.backends.backend_pdf import PdfPages
#pc="PC1"
with PdfPages('./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_vs_batches_density.pdf') as pdf:
    for pc in pc_columns:

        #define the figure size
        plt.figure(figsize=(10, 6))

        #density plot of the values of selected PC for batch 1
        plt.subplot(1, 2, 1)  #1 row, 2 columns, first plot
        pca_results.loc[pca_results['#FID'] == fid1, pc].plot.kde()
        plt.title(f'Density Plot of {pc} for {fid1}')
        plt.xlabel('Values')
        plt.ylabel('Density')

        # Second subplot
        plt.subplot(1, 2, 2)  #1 row, 2 columns, second plot
        pca_results.loc[pca_results['#FID'] == fid2, pc].plot.kde()
        plt.title(f'Density Plot of {pc} for {fid2}')
        plt.xlabel('Values')
        plt.ylabel('Density')

        plt.tight_layout()  # Adjusts subplot params so that subplots fit into the figure area
        pdf.savefig()  # Saves the current figure into the pdf file
        plt.close() #to close the current figure, freeing up memory
    #Distribution tend to be the same between batches, there are only differences in PC9, but in reality the different is not so big if you look at X axis values. Note that we have large sample size, so the violation of this assumption is not so bad.

print("make Mann–Whitney U test for the difference of PC values between batches")
from scipy.stats import mannwhitneyu
#Mann–Whitney U test is a nonparametric test of the null hypothesis (H0) that, for randomly selected values X and Y from two populations, the probability of X being greater than Y is equal to the probability of Y being greater than X. 
    #https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test

#Assumptions of the Mann-Whitney U test:
    #1. **Independence of Observations**: This is an assumption of nearly all statistical tests, including the Mann-Whitney U test. It means that the observations in each group must be independent of each other. In other words, knowing the value of one observation should not provide any information about the value of any other observation.
        #PARTIAL PROBLEM: We have removed related individuals, still we can have some degree of non-independency due to population stratification
    #2. **Ordinal Measurement or Higher**: The Mann-Whitney U test requires that the data be at least ordinal. This means that the data can be logically ordered. For example, a Likert scale (e.g., strongly disagree, disagree, neutral, agree, strongly agree) is an example of ordinal data. The test can also be used with interval or ratio data, which are higher levels of measurement.
        #OK
    #3. **Similar Distribution Shapes**: The Mann-Whitney U test assumes that the distributions of both groups are similar in shape. If the shapes of the distributions are very different, the Mann-Whitney U test may not be the most appropriate test to use. Note that this assumption is less crucial when the sample sizes are large.
        #OK: According to chatGTP, it is usually considered 30 or more large sample size. In our case (with 200 and 1200 samples), we should be ok regarding this assumption.
        #Central Limit Theorem states that when independent random variables are added, their properly normalized sum tends toward a normal distribution even if the original variables themselves are not normally distributed. Therefore, we should be ok as long ass both populations, i.e., both batches, are large enough.
        #I have checked this anyways.... (see above)
    #4. **Random Sampling from Populations**: The observations should be randomly sampled from the population. This is to ensure that the samples are representative of the populations they are drawn from.

#calculate the wilcoxon U test across PCs
#pc="PC1"
wilcoxon_results = list()
for pc in pc_columns:

    #extract the values of the selected PC for each of the batches
    group1 = pca_results.loc[pca_results['#FID'] == fid1, pc]
    group2 = pca_results.loc[pca_results['#FID'] == fid2, pc]
    
    #perform the Mann-Whitney U test
    stat, p = mannwhitneyu(group1, group2)
    
    #save the results in a tuple and append to list
    wilcoxon_results.append((pc, stat, p))

#convert the results to DF
wilcoxon_results_df = pd.DataFrame(wilcoxon_results, columns=["PC", "U_statistic", "P"])
print(wilcoxon_results_df)

#check
if (wilcoxon_results_df.shape[0]!=len(pc_columns)) | (sum(wilcoxon_results_df["P"]<0.05)!=0):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE WILCOXON TEST FOR DIFFERENCES BETWEEN BATCHES")
else:
    print("The wilcoxon test for all the PCs between batches are non-significant. Therefore, we accept the null hypothesis: for randomly selected PC values X and Y from two batches, the probability of X being greater than Y is equal to the probability of Y being greater than X (see code for details)")

print_text("select PCs based on explianed variance", header=4)
print("load eigenvalues")
#In the context of Principal Component Analysis (PCA), eigenvalues represent the amount of variance explained by each principal component³. 
#PCA is a statistical procedure that uses an orthogonal transformation to convert a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables called principal components¹. This transformation is defined in such a way that the first principal component has the largest possible variance (that is, accounts for as much of the variability in the data as possible), and each succeeding component in turn has the highest variance possible under the constraint that it is orthogonal to the preceding components¹.
#The eigenvalues in the PCA output file generated by PLINK2's `--pca` function correspond to the variances along these principal components¹. The larger an eigenvalue, the more variance (information) in your data that the corresponding principal component explains³. 
#In other words, eigenvalues in this context give you an idea of how much of the total variance of your data is captured by each principal component. This can be useful when deciding how many principal components to retain for further analysis³.
#In addition, EIGENSOFT README says that "eigenvalue_k/(Sum of eigenvalues) is the proportion of variance explained by eigenvector_k."
#Please let me know if you need further clarification or have any other questions!
    #Source: Conversation with Copilot, 7/3/2024
    #(1) What are Eigenvalues & Eigenvectors in PCA? (Example) - Statistics Globe. https://statisticsglobe.com/what-are-eigenvalues-eigenvectors-pca.
    #(2) Population structure: PCA - Speciation & Population Genomics: a how-to .... https://speciationgenomics.github.io/pca/.
    #(3) plink_pca : Estimate principal components with plink2. https://rdrr.io/github/OchoaLab/genbin/man/plink_pca.html.
pca_eigenvalues = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_results.eigenval", \
    sep="\t", \
    header=None, \
    low_memory=False)
if(pca_eigenvalues.shape[0]!=(pca_results.shape[1]-2)):
    #pca_results has one column per PC and then the familiy and sample IDs as two additional columns
    raise ValueError("ERROR! FALSE! WE DO NOT HAVE THE CORRECT NUMBER OF EIGENVALUES")
else:
    #set row names as numbers 1 to 10
    pc_list = [f"P{i}" for i in range(1, pca_eigenvalues.shape[0]+1)]
        #the "end" of range is not included, so setting the end to 20 would end at 19, but we want the 20 PC axes, so we add 1.
    pca_eigenvalues.index = pc_list
    print(pca_eigenvalues)

print("calculate the proportion of explained variance")
#The proportion of total variance explained by each PC is a useful metric for understanding structure in a sample and for evaluating how many PCs one might want to include in downstream analyses. This can be computed as λi∕∑kλk, with λi being eigenvalues in decreasing order
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #The formula λi∕(∑k λk) is used in the context of **Principal Component Analysis (PCA)**, a statistical procedure that uses an orthogonal transformation to convert a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables called principal components.
    #Here's what each part of the formula means:
        #- **λi**: This represents the **i-th eigenvalue** of the covariance matrix of the dataset. Eigenvalues are a special set of scalars associated with a linear system of equations (i.e., a matrix equation) that are sometimes also known as characteristic roots, characteristic values, proper values, or latent roots.
        #- **∑k λk**: This is the **sum of all eigenvalues** of the covariance matrix. 
        #The fraction $$\frac{\lambda_i}{\sum_k \lambda_k}$$ thus represents the **proportion of the total variance** in the data that is explained by the i-th principal component. This is often used to understand the importance of each principal component in the analysis and decide how many principal components to retain for further analysis.
        #In the sentence you provided, this value is being plotted, likely to create a **scree plot**. A scree plot is a simple line segment plot that shows the fraction of total variance in the data as explained or represented by each PC. This is useful to determine the appropriate number of principal components to retain for further analysis. The plot helps to visualize the **explained variance** by each principal component and it typically decreases and eventually becomes flat with increasing components, which can help to choose the optimal number of components. 
pca_eigenvalues[1] = pca_eigenvalues[0]/pca_eigenvalues[0].sum()
#change column names
pca_eigenvalues = pca_eigenvalues.rename(columns={0:"eigenvalue", 1:"prop_variance"})
pca_eigenvalues

print("make de scree plot")
#A scree plot (Fig. 6A in Ritchie´s tutorial) can be used to evaluate and determine how many principal components would be appropriate for the covariates; the bend in the line plot, known as the “elbow”, denotes the location of the PCs that should be selected (Cattell, 1966).
#add a new column with the index (row numbers) and plot index vs the first and only column (eigenvalues)
pca_eigenvalues.reset_index().plot( \
    x='index', \
    y="prop_variance", \
    style='-o', \
    legend=None)
    #plot a line with a dot in eahc observation
    #reset_index() is called on the pca_eigenvalues DataFrame. This resets the index of the DataFrame, and the old index is added as a column named 'index'. This new DataFrame is then plotted using the plot function, with 'index' as the x-values and the first column (0) as the y-values.

#add labels
plt.title('Scatter Plot of PCA Eigenvalues')
plt.xlabel('PCs')
plt.ylabel('Proportion of variance explained')

#add xticks from the first to the last PC
plt.xticks(range(0, pca_eigenvalues.shape[0]), pc_list, fontsize=8)

#save and close
plt.savefig( \
    fname="./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_eigenvalues.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()

#The results show that the first PC explains a big fraction of the variance (46%). The next 8 (until PC9) explain between 7.5 and 2.1% of variability. After the 10, the variance explained per PC becomes relatively constant (around 2.1%).

print("see the first PC where cumulative explianed variability reaches 80% of the total")
#According to Ritchie´s review, it is recommended to use the eigenvectors that explain the greatest proportion of variance (on an average cumulative 80% variance) as covariates in any downstream analysis to correct for bias due to population stratification.

#for each eigenvalue
#pc_short="P1"; row=pca_eigenvalues.iloc[0,:]
import numpy as np
for pc_short, row in pca_eigenvalues.iterrows():

    #if the row is not the last one, i.e., the las PC
    if(pc_short!="P"+str(pca_eigenvalues.shape[0])):

        #get the position of the current row
        current_index = np.where(pca_eigenvalues.index == pc_short)[0][0]

        #take all the values from the first to the current row and sum them all
        cumulative_variability = \
            pca_eigenvalues.iloc[ \
                0:current_index+1, \
                np.where(pca_eigenvalues.columns=="prop_variance")[0][0]].sum()
            #to get the row of "index", we need need "index"+1 because the last number is not included
        
        #if the cumulative variability is above 76%
        if(cumulative_variability>=0.76):
            
            #print the name of PC and stop de loop
            print(f"The first PC where cumulative explianed variability reaches 76% of the total is {pc_short} with {cumulative_variability:.2f} of the total")
            
            #save the PC as a variable
            last_pca_to_consider = "PC"+pc_short.split("P")[1]
            break

#Summary:
#The elbow is around P9 (almost no changes after this axis), and that axis is already explaining 77% of variability, which is very close to 80%. This is below the limit usually cosnidered of 10 axis and we have enough observations to use these axes as covariates plus others confounding variables assuming at least 10 observations per covariate. We would have 9 axes plus 3 confounding variables makes 12 covariates in total. 12*10 makes 120 observations and have more than 1300.
    #Marees says that in psychiatric genetics community, up to 10 axes of MDS are accepted. It dependes of the population and your sample size, as you have more samples, you can use more axes in the glms. 
    #They say the same in this tutorial. Also they say that 10 is ok as covariates, but they only use 2 for structure detection. 
        #https://gwas-intro-bajicv.readthedocs.io/en/latest/05_pop_stratification/
    #See copiltor summary:
        #In Genome-Wide Association Studies (GWAS), Principal Component Analysis (PCA) is a standard method for estimating population structure and sample ancestry¹. The inclusion of Principal Components (PCs) as covariates in the models helps to control for population stratification¹².
        #The number of PCs to be included as covariates depends on the population structure and the sample size³⁵. Generally, the inclusion of **up to 10 components** is accepted³⁵. However, some studies have found that **five PCs** are generally sufficient to correct for stratification in simulated and real data sets⁴. Alternatively, the number of PCs may be selected through cross-validation or Tracy–Widom statistics⁴.
        #It's important to note that these are general recommendations and the optimal number of PCs to include might vary depending on the specific dataset and research question. Therefore, it's always a good idea to perform some exploratory data analysis and consider the specific characteristics of your dataset when deciding on the number of PCs to include in your models.
            #Source: Conversation with Copilot, 7/4/2024
            #(1) Controlling for stratification in (meta-)GWAS with PCA: Theory .... https://www.broadinstitute.org/talks/controlling-stratification-meta-gwas-pca-theory-applications-and-implications.
            #(2) Robust methods for population stratification in genome wide association .... https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-132.
            #(3) Population Stratifiction - Introduction to GWAS. https://gwas-intro-bajicv.readthedocs.io/en/latest/05_pop_stratification/.
            #(4) 4 Quality Control: Relatedness & Population Stratification. https://hds-sandbox.github.io/GWAS_course/Notebooks/GWAS4-QualityControlB/.
            #(5) Improving the Power of GWAS and Avoiding Confounding from Population .... https://academic.oup.com/genetics/article/197/3/1045/5935993.


print_text("PCA with smartpca", header=4)
#smartpca: general explanations
    #smartpca runs Principal Components Analysis on input genotype data and  outputs principal components (eigenvectors) and eigenvalues. We note that eigenvalue_k/(Sum of eigenvalues) is the proportion of variance explained by eigenvector_k. The method assumes that samples are unrelated. However, a small number of cryptically related individuals is usually  not a problem in practice as they will typically be discarded as outliers (we already removed related individuals).
    #The syntax of smartpca is "../bin/smartpca -p parfile"
        #The below for details about each argument

print("prepare fam file for smartpca")
#As a minor issue, smartpca ignores individuals in the .fam file if they are marked as missing in the phenotypes column. This awk command provides a new .fam file that will automatically include all individuals.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
#copy also bim and bed files
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/; \
    mkdir -p ./eigen_out; \
    awk \
        '{print $1,$2,$3,$4,$5,1}' \
        ./loop_maf_missing_2_ldprunned_autosome_pca.fam > \
    ./eigen_out/loop_maf_missing_2_ldprunned_autosome_pca.new_fam.fam; \
    cp ./loop_maf_missing_2_ldprunned_autosome_pca.bim ./eigen_out; \
    cp ./loop_maf_missing_2_ldprunned_autosome_pca.bed ./eigen_out \
")

print("decide the number of axes to output")
n_axes=20

print("smartpca parameter file")
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca; \
    PREFIX=loop_maf_missing_2_ldprunned_autosome_pca; \
    echo genotypename: ./eigen_out/$PREFIX.bed > ./eigen_out/$PREFIX.par; \
    echo snpname: ./eigen_out/$PREFIX.bim >> ./eigen_out/$PREFIX.par; \
    echo indivname: ./eigen_out/$PREFIX.new_fam.fam >> ./eigen_out/$PREFIX.par; \
    echo snpweightoutname: ./eigen_out/$PREFIX.snpeigs >> ./eigen_out/$PREFIX.par; \
    echo evecoutname: ./eigen_out/$PREFIX.eigs >> ./eigen_out/$PREFIX.par; \
    echo evaloutname: ./eigen_out/$PREFIX.eval >> ./eigen_out/$PREFIX.par; \
    echo phylipoutname: ./eigen_out/$PREFIX.fst >> ./eigen_out/$PREFIX.par; \
    echo numoutevec: " + str(n_axes) + " >> ./eigen_out/$PREFIX.par; \
    echo outliersigmathresh: 6 >> ./eigen_out/$PREFIX.par; \
    echo numoutlieriter: 11 >> ./eigen_out/$PREFIX.par; \
    echo numoutlierevec: 10 >> ./eigen_out/$PREFIX.par; \
    echo outlieroutname: ./eigen_out/$PREFIX.outliers >> ./eigen_out/$PREFIX.par; \
    echo altnormstyle: YES >> ./eigen_out/$PREFIX.par; \
    echo missingmode: NO >> ./eigen_out/$PREFIX.par; \
    echo ldregress: 0 >> ./eigen_out/$PREFIX.par; \
    echo noxdata: YES >> ./eigen_out/$PREFIX.par; \
    echo nomalexhet: YES >> ./eigen_out/$PREFIX.par; \
    echo newshrink: NO >> ./eigen_out/$PREFIX.par \
")
    #script for smartpca parameter file from ADMIXTURE TUTORIAL
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #genotypename:
        #contains genotype data for each individual at each SNP
    #snpname:
        #contains information about each SNP 
    #indivname:
        #contains information about each individual
        #accorind to the readme, the genotype and snp file can be just .bed and.bim files, respectively, from plink format, which is our format. 
        #according to the admixture tutorial, you can use the plink format, i.e., bed, bim and fam file as input. In the case of the fam file, the phenotypes have to be 1 for all samples, if not, they are not considered (see above)
            #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #snpweightoutname:
        #output file containing SNP weightings of each principal component. Note that this output file does not contain entries for monomorphic SNPs from the input .snp file. 
    #evecoutname:
        #output file of eigenvectors. See numoutevec parameter below.
    #evaloutname:
        #output file of all eigenvalues
    #phylipoutname:
        #output file containing an fst matrix which can be used as input to programs in the PHYLIP package, such as the "fitch" program for constructing phylogenetic trees.
        #we are generating this just in case.
    #numoutevec:
        #number of eigenvectors to output.  Default is 10.
        #this should affect the results, for example, the axes used for outlier removal are indicated in another argument.
    #outlier removal. According to plink tutorial (https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3)
        #It is also a good idea to throw out gross outliers at this point; any sample which is more than, say, 8 standard deviations out on any top principal component is likely to have been genotyped/sequenced improperly; you can remove such samples by creating a text file with the bad sample IDs, and then using –remove +  –make-bed:
            #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #The phrase "8 standard deviations out" refers to a data point that falls 8 standard deviations away from the mean of a dataset. In statistics, the standard deviation is a measure of the amount of variation or dispersion in a set of values. A low standard deviation means that the values tend to be close to the mean, while a high standard deviation means that the values are spread out over a wider range. When a data point is said to be "8 standard deviations out", it means it's quite far from the mean. In a normal distribution, almost all data falls within 3 standard deviations of the mean. So, a data point that is 8 standard deviations away from the mean is extremely rare and could be considered an outlier. In the context of your sentence, it suggests that any sample which falls more than 8 standard deviations away on any top principal component might have been genotyped or sequenced improperly, and thus, it might be a good idea to remove these outliers from the analysis. This is because such extreme values can significantly skew the results and interpretations of the data analysis.
        #If you do this, follow it up by repeating the PCA, since the bad samples might have distorted the principal components. Occasionally, the new principal components will reveal another bad sample, and you have to repeat these two steps, etc. We are using the iterative process of smartpca as recommended in the tutorial
        #Steps
            #The program identifies individuals whose principal component scores deviate significantly from the majority of the population. This is typically done by setting a threshold for the number of standard deviations from the mean.
            #Iterative Process: Outlier removal is often an iterative process. The program may run PCA multiple times, each time removing individuals who are identified as outliers, until no more outliers are detected.
            #Impact on Analysis: Removing outliers helps to ensure that the PCA results are not skewed by individuals who have unusual genetic backgrounds. This can provide a clearer picture of the population structure and improve the accuracy of downstream analyses.
        #outliersigmathresh: 
            #We are using the default (6).
            #number of standard deviations which an individual must exceed, along one of the top (numoutlierevec) principal components, in order for that individual to be removed as an outlier.  Default is 6.0.
        #numoutlieriter:
            #maximum number of outlier removal iterations. Default is 5.  To turn off outlier removal, set this parameter to 0.
            #I have checked that outliers are removed until iteration 7, so we just using until 8.
        #numoutlierevec:
            #number of principal components along which to remove outliers during each outlier removal iteration.  Default is 10.
            #We are using the default. The PCA axes with a significant p-value (<0.05) are the first 11, but I have checked that using 11 instead of 10 for outlier removal gives the same outliers. In addition, I have checked the eigenvectors that detect the outliers and never 10 or 11 are included, so we do not need to consider more axes. We are sticking to the default.
        #outlieroutname:
            #output logfile of outlier individuals removed. If not specified, smartpca will print this information to stdout, which is the default.
    #altnormstyle:
        #Affects very subtle details in normalization formula. Default is YES (normalization formulas of Patterson et al. 2006). To match EIGENSTRAT (normalization formulas of Price et al. 2006), set to NO.
        #I have checked with and without altnormstyle and there are almost no differences.
    #missingmode:
        #If set to YES, then instead of doing PCA on # reference alleles, do PCA on whether each data point is missing or nonmissing.  Default is NO.
    #ldregress:
        #If set to a positive integer, then LD regression is turned on, and input to PCA will be the residual of a regression involving that many previous SNPs, according to physical location.  See Patterson et al. 2006. Default is 0 (no LD regression).  If desiring LD correction, we recommend 200.
        #this is done to correct for the correlation (i.e., the linkage disequilibrium) of the SNPs. We do NOT to do this because we have already prunned our data considering LD.
    #noxdata: 
        #if set to YES, all SNPs on X chr are excluded from the data set. The smartpca default for this parameter is YES, since different variances for males vs. females on X chr may confound PCA analysis.
        #Using YES, but not required anyway because we have only autosomal SNPs.
    #nomalexhet:
        #if set to YES, any het genotypes on X chr for males are changed to missing data. The smartpca default for this parameter is YES.
        #males should be homozigous for X (except pseudo-autosomic regions)
        #Using YES, but not required anyway because we have only autosomal SNPs.
    #lsqproject
        #PCA projection is carried out by solving least squares equations rather than an orthogonal projection step. This is approriate if PCs are calculated using samples with little missing data but it is desired to project samples with much missing data onto the top PCs. In other words, most of the samples have a lot of data, but some samples have a lot of missing, in that situation, instead of filling gaps with the average, this approach does something different that solves the problem. BUT, if the sample has few missing, then this works as the default orthogonal. I have checked that setting this to NO or YES does not change the results.
        #./EIG-8.0.0/POPGEN/lsqproject.pdf
    #shrinkmode/newshrink
        #A problem with smartpca is that samples used to calculate the PC axes "stretch" the axes. So that 2 populations in fact genetically identical (2 independent samples from the same underlying population) will appear different if one is used to compute axes, and one not.  shrinkmode: YES is an attempt to solve this problem.  Details to appear later, but this has been used successfully in the Reich lab.*** warning *** shrinkmode is slow and will greatly increase the runtime. (NEW) New version:  newshrink:  YES technical variation of shrinkmode,  should be (slightly). 
        #According to chatGTP, even within a single population, if there is significant internal structure, `shrinkmode` can help ensure that the PCA axes are not unduly influenced by subgroups within your population.
        #HOWEVER, I have run with and without shrink mode and the results are EXACTLY the same, with all decimals, while the run time goes up to 40 min from 5-10. We are not using this mode.
    #Multithreading (Code added by Chris Chang).
        #smartpca now supports multithreading but NOT with fastmode: YES. By default a (hopefully) system dependent number of threads is chosen. This can be overwritten by (for example) numthreads:   10
        #it seems by default it is using all the cores.

print("run smartpca")
#we use the version 8.0.0, which is the latest at the moment of writting (sep 2024). This was released in october 2022.
    #see the contianer receipte for further details about how install in the container
    #https://github.com/DReichLab/EIG/releases
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/; \
    /bin/smartpca -p ./eigen_out/loop_maf_missing_2_ldprunned_autosome_pca.par > ./eigen_out/smart_pca_run.out")
    #smartpca also calculate the "Tracy-Widom statistics"
        #I have visually compared this table and the generated by twstats. They are the same, just a few decimals are different. Also that table has a column with eigenvalues and a p-value. Therefore, smartpca is also using twstats
            #see below for the code to ran twstats
        #The twstats program computes Tracy-Widom statistics to evaluate the statistical significance of each principal component identified by pca (Patterson et al. 2006).
        #it helps determine whether the observed eigenvalues (which correspond to the variance explained by each PC) are significantly larger than what would be expected by chance.
        #steps (from copilot):
            #Eigenvalue Calculation: After performing PCA, the eigenvalues corresponding to each principal component are calculated. These eigenvalues represent the amount of variance explained by each PC.
            #Comparison with Tracy-Widom Distribution: The largest eigenvalues are compared to the Tracy-Widom distribution. This comparison helps determine if the observed eigenvalues are significantly larger than those expected under the null hypothesis (i.e., no structure in the data, the axes are not absorbing any structure from the data).
            #Significance Testing: The program computes p-values for each eigenvalue based on the Tracy-Widom distribution. A low p-value indicates that the corresponding principal component explains a significant amount of variance, suggesting it captures meaningful structure in the data.
        #The twstats program assumes a random set of markers, and should not be used on data sets of ancestry-informative markers, as admixture-LD may violate its underlying assumptions. 
            #Ancestry-informative markers are those SNPs that significantly different between pops and are speficically used for infer ancestry of individuals.
            #This is not our case, as we have not specifically selected SNPs that differ between pops. 
            #Anyways, we are going to compare the significant PCA axes according to twstats and admixture program, so we are good here.

#in case you want to run twstats separately
'''
if False:
    run_bash(" \
        cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/; \
        /usr/bin/twstats \
            -t ../../../../../../eigensoft_versions/EIG-8.0.0/POPGEN/twtable \
            -i ./eigen_out/loop_maf_missing_2_ldprunned_autosome_pca.eval \
            -o ./eigen_out/loop_maf_missing_2_ldprunned_autosome_pca_tracy.out")
'''

#Important note about the number of SDs respect to the mean of any PC axis to remove outliers:
    #The default is 6, using this we lose around 171 samples, and we do not have clear clusters in the main axes. Also ancestry considers that we have just 1 group. 
    #Using 3 makes the PCA plots to show exactly as clouds of points, but the caveat is that we lose 500 samples. Also note that doing this still showed the problem of negatively correlated SNPs between our dataset and the reference panel in the imputation server. That was only solved after selecting the correct mode for solving flips and remove SNPs that are not possibly solved (ambiguous see below).
    #Therefore, it seems we already have a relatively homogenous dataset using the default parameters for outlier removal.

print("see significant axes")
#extract the tracy table
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/; \
    awk \
        'BEGIN{ \
            OFS=\"\t\"; \
            start_table=0; \
            end_table=0; \
        }{ \
            if($0 ~ /^## Tracy-Widom statistics: rows:/){ \
                start_table=NR \
            }; \
            if($0 ~ /kurtosis/){ \
                end_table=1 \
            }; \
            if(start_table !=0 && NR > start_table && end_table==0){ \
                print $1,$2,$3,$4,$5,$6 \
            }; \
        }' \
        smart_pca_run.out > tracy_table.tsv \
    ")
    #force the output to be tab, and create two variables as zero
    #if the has the header of tracy table
        #save the number of the row to know when the table start
    #if the row includes kurtosis
        #we are at the end of the table, so activate ending
    #if start_table is not zero, and it is greater than the first row of the tracy table and we are not yet at the end of the table 
        #print all fields separately

#load the table
tracy_table = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/tracy_table.tsv", \
    sep="\t", \
    header=0, \
    low_memory=False)
print(tracy_table)

#print significant axes
tracy_signi_axes = tracy_table.loc[tracy_table["p-value"] < 0.05,:]
print(tracy_signi_axes)
print(f"We have {tracy_signi_axes.shape[0]} significant axes")

#Results
#6 PCA axes are significant (P<0.05), meaining that they explain more genetic variance than expected by chance. Plink PCA gaves 9 relevant axes before outlier removal.

print("plot smartpca eigenvectors")
#process the eigenvector file to ensure is tab delimited
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/; \
    awk \
        'BEGIN{ \
            OFS=\"\t\" \
        }{ \
            if(NR>1){ \
                print " + ",".join(f"${i}" for i in range(1,n_axes+2)) + " \
            } \
        }'\
        ./loop_maf_missing_2_ldprunned_autosome_pca.eigs > ./loop_maf_missing_2_ldprunned_autosome_pca_tab.eigs \
")
    #force the OFS to be tabs
    #select all except the first 1, i.e., NR>1. The first one has eigenvalues of rach PCA axis. For each row
    #print the following columns
        #all columns from the first to the last one
        #the last is based on the number of PCA axes we have plus 1 because we have an ID column and 1 more because the last element of the range is not included in python
    #this effectively prints all columns as tab separated

#load the table
smartpca_eigenvectors = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/loop_maf_missing_2_ldprunned_autosome_pca_tab.eigs", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(smartpca_eigenvectors)
#check we have the correct number of columns
if(smartpca_eigenvectors.shape[1]-1!=n_axes):
    raise ValueError("ERROR! FALSE! We are not considering all the interesting axes in smartpca")

#plot the ifrst 5 axes againts each other
import matplotlib.pyplot as plt
import seaborn as sns
sns.pairplot(smartpca_eigenvectors.iloc[:,1:6])
plt.savefig("./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/eigen_vectors_pairwise.png", dpi=300)
plt.close()
    #PC1 is clearly the axis absorbing genetic structure. There are at least 2 clear groups, and maybe 3. This can be seen when plotting PC1 vs the other axes. In contrast, when plotting the other axes against each other, there are clouds with a lot of contiunity. This support the idea that PC1 is very relevant (see explained variance below).
    #We are going to use the groups based on PC1 and check with admixture whether is better to assume 1, 2 or 3 ancestry groups in our data.
        #According to the plink tutorial: "If there are obvious clusters in the first few plots, I recommend jumping ahead to Chapter 4 (on ADMIXTURE) and using it to label major subpopulations before proceeding."

print("plot the explained variance using eigenvalues")
#load the eigenvalues
smartpca_eigenvalues = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/loop_maf_missing_2_ldprunned_autosome_pca.eval", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(smartpca_eigenvalues)
    #I have checked these eigenvalues are the same than in the first row of the eigenvector table, which is named as "#eigvals"

#check eigenvalues are in decreasing order
is_decreasing = smartpca_eigenvalues.sort_values(by=0, ascending=False).equals(smartpca_eigenvalues)
if(is_decreasing == False):
    raise ValueError("ERROR! FALSE! The eigen values are not in decreasing order")

#The proportion of total variance explained by each PC is a useful metric for understanding structure in a sample and for evaluating how many PCs one might want to include in downstream analyses (see Fig. 9). This can be computed as λi∕∑kλk, with λi being eigenvalues in decreasing order
smartpca_eigenvalues[1] = smartpca_eigenvalues[0]/smartpca_eigenvalues[0].sum()
#smartpca_eigenvalues[2] = smartpca_eigenvalues.iloc[:, 1].cumsum(axis=0)
#check
if(smartpca_eigenvalues[1].sum() != 1.0):
    raise ValueError("ERROR! FALSE! We have not correctly calculated the proportion of explained variance")

#change column names
smartpca_eigenvalues = smartpca_eigenvalues.rename(columns={0:"eigenvalue", 1:"prop_variance"})
print(smartpca_eigenvalues)
    
#select the eigenvalues of the first 20 axes
smartpca_eigenvalues_20 = smartpca_eigenvalues.iloc[0:20, :]
print(smartpca_eigenvalues_20)

#add a new column with the index (row numbers) and plot index vs the first and only column (eigenvalues)
import matplotlib.pyplot as plt
smartpca_eigenvalues_20.reset_index().plot( \
    x='index', \
    y="prop_variance", \
    style='-o', \
    legend=None)
    #plot a line with a dot in eahc observation
    #reset_index() is called on the pca_eigenvalues DataFrame. This resets the index of the DataFrame, and the old index is added as a column named 'index'. This new DataFrame is then plotted using the plot function, with 'index' as the x-values and the first column (0) as the y-values.

#add labels
plt.title('Scatter Plot of PCA Eigenvalues')
plt.xlabel('PCs')
plt.ylabel('Proportion of variance explained')

#add xticks from the first to the last PC
plt.xticks( \
    ticks=smartpca_eigenvalues_20.index, \
    labels=range(1, smartpca_eigenvalues_20.shape[0]+1), \
    fontsize=8)
    #the positions are just the index of the eigenvalues
    #the labels are 1 to 20, but add 1 at the end because the last element of the range is not included

#save and close
plt.savefig( \
    fname="./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/explained_var.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()
    
#Results
    #PC1 explains 1.4% of variability, while PC2 explains 0.135%. From there, next PCAs until the number 20 explains around 0.10%. Therefore, the main reduction of explained variability occurs from PC1 to PC2 (1.27%) and a from PC2 to PC3 (0.01%). From there, the reductions are much smaller. Therefore, the main axes seems to be PC1 and PC2.
    #Note that the total explained variance is much less compared to the results of plink where variance explained by the first axis was of 40%! but it is in line what the admixture tutorial got with smartpca, so maybe this is a question of this approach and/or the removal of outliers
    #we are not calculating cumultaive percentage of explained variance because we have a small proportion explained by each axis, and the rule of Ritche of 80% would require to include hundreds of axes...

print("plot SNP weights")
#process the file to ensure it is tab delimited
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/; \
    awk \
        'BEGIN{ \
            OFS=\"\t\" \
        }{ \
            if(NR>0){ \
                print " + ",".join(f"${i}" for i in range(1,n_axes+4)) + " \
            } \
        }'\
        ./loop_maf_missing_2_ldprunned_autosome_pca.snpeigs > ./loop_maf_missing_2_ldprunned_autosome_pca_tab.snpeigs \
")
    #we to sum 1 because the last number is not included in the range
    #also 3 more because we have a column with the IDs and another with the chromosome and the phyisical position (I checked with the bim file, identical all cases)

#load the file
smartpca_snp_weights = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/loop_maf_missing_2_ldprunned_autosome_pca_tab.snpeigs", \
    sep="\t", \
    header=None, \
    low_memory=False)
smartpca_snp_weights
#check
if(smartpca_snp_weights.shape[1] != n_axes+3):
    raise ValueError("ERROR! FALSE! We have a problem with the number of axes in the SNP weights")
    #besides axes, we have 3 columns (SNP id, chromosome and position)

#sort the DF by chromosome and position
if(smartpca_snp_weights.sort_values(by=[1, 2]).equals(smartpca_snp_weights) != True):
    raise ValueError("ERROR! FALSE! Table with snp weights is not sorted by position and chromsome")

#change column names
smartpca_snp_weights.columns = ["rs_number", "chr", "pos"] + [f"pca_{i}" for i in range(1,n_axes+1)]
    #the axes are named as pca_1, pca_2... until pca_n_axes
print(smartpca_snp_weights)

#melt the DataFrame to have one row for each combination of chromosome and PCA axis
smartpca_snp_weights_melted = pd.melt( \
    smartpca_snp_weights, \
    id_vars=["rs_number", "chr", "pos"], \
    value_vars=[f"pca_{i}" for i in range(1,n_axes+1)], \
    var_name="pca_axis", \
    value_name="snp_weight")
    #This function is useful to massage a DataFrame into a format where one or more columns are identifier variables (`id_vars`), while all other columns, considered measured variables (`value_vars`), are "unpivoted" to the row axis, leaving just two non-identifier columns, 'variable' and 'value'. In other words, pd.melt transforms the DataFrame from wide format to long format, with one row for each combination of chromosome, rs_numbre, position, and PCA axis. 
        #The id_vars parameter specifies the columns to use as identifier variables. 
            #The specific combination of chromosome, SNP and position
        #the value_vars parameter specifies the columns to unpivot.
            #we are taking all the PCA_axes columns and put them all together in the same column as different rows
        #The var_name and value_name parameters specify the names of the new columns.
print(smartpca_snp_weights_melted)

#check we have all the rows
if(smartpca_snp_weights[["rs_number", "chr", "pos"]].shape[0]*n_axes!=smartpca_snp_weights_melted.shape[0]):
    raise ValueError("ERROR! FALSE! Problem when melting the DF with SNP weights, we do not have the expected number of rows")
    #94854 SNPs times 20 PCA axes makes 1897080 rows, which is the number of rows of smartpca_snp_weights_melted

#get unique PCA axes and chromosomes
pca_axes = smartpca_snp_weights_melted['pca_axis'].unique()
chromosomes = smartpca_snp_weights_melted['chr'].unique()

#create subplots
#define dimensions
import matplotlib.pyplot as plt
fig, axes = plt.subplots( \
    nrows=len(pca_axes), \
    ncols=len(chromosomes), \
    figsize=(4*len(chromosomes), 4*len(pca_axes)))
    #define the number of rows/columns of the subplot grid.

#make the plots
#for each PCA axis (i starts at 0)
#pca_axis="pca_1"; chromosome=1; i=0; j=0
for i, pca_axis in enumerate(pca_axes):
    
    #for each chromosome (j starts at 0)
    for j, chromosome in enumerate(chromosomes):
        
        #select data for this PCA axis and chromosome
        snp_weights_subset = smartpca_snp_weights_melted[(smartpca_snp_weights_melted["pca_axis"] == pca_axis) & (smartpca_snp_weights_melted["chr"] == chromosome)]
        
        #create scatter plot
        axes[i, j].scatter(x=snp_weights_subset["pos"], y=snp_weights_subset["snp_weight"])
        
        #set x-axis range extending a bit from both sides
        axes[i, j].set_xlim(snp_weights_subset["pos"].min()*0.9, snp_weights_subset["pos"].max()*1.1)
        
        #set title
        axes[i, j].set_title(f"PCA Axis: {pca_axis}, Chromosome: {chromosome}")

#adjust layout
plt.tight_layout()

#save and close
plt.savefig( \
    fname="./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/snps_weights_across_axes_chromosomes.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()

#results
    #it is useful to inspect the PC loadings to ensure that they broadly represent variation across the genome, rather than one or a small number of genomic regions. SNPs that are selected in the same direction as genome-wide structure can show high loadings, but what is particularly pathological is if the only SNPs that show high loadings are all concentrated in a single region of the genome, as might occur if the PCA is explaining local genomic structure (such as an inversion) rather than population structure.
    #Everything clean except PC5 for chromosome 8, there is a clear peak.

print("Remove outliers")
pca_outliers = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/loop_maf_missing_2_ldprunned_autosome_pca.outliers", \
    sep=" ", \
    header=None, \
    low_memory=False)
print(pca_outliers)
print(f"We have removed {pca_outliers.shape[0]} samples because they were PCA outliers")

#Results
#171 samples have been removed because they are very away from the mean of the PCA axes, specifically, more than 6 standard deviations.
#It is true that the plink tutorial says that outliers away more than 8 SDs, are "likely caused due to genotyping errors and it is recommended to remove them". However, the "predict-HIIT" study used 6 SD to remove PCA ancestry outliers. In addition, the default of smartpca is 6, not 8. Also, The predict-hiit study only conisdered the first 2 axes and losed 10 samples because of this, but smartpca considers 10 as default. Finally, most important is the fact that the removing all these outliers (6SDs and considering 10 axes) makes the data much more clear, with no genetic structure whatsoever.

#save the IDs in a file
pca_outliers[2].str.split(":", expand=True).to_csv("./data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/smartpca_outliers.txt",
    sep="\t",
    header=False,
    index=False)
    #first split the ID column into to two columns (expand=True) so we separate family and sample ID, which is the format expected by --remove of plink

#remove these outliers from the plink files of the LD subset
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca; \
    plink \
        --bfile ./loop_maf_missing_2_ldprunned_autosome_pca \
        --remove ./eigen_out/smartpca_outliers.txt \
        --make-bed \
        --out ./loop_maf_missing_2_ldprunned_autosome_pca_not_outliers"
)
    #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.

#check if the difference in samples between the original and new fam files is the number of outliers
if( \
    pd.read_csv("./data/genetic_data/quality_control/14_pop_strat/01_pca/loop_maf_missing_2_ldprunned_autosome_pca.fam", sep=" ", header=None, low_memory=False).shape[0] - \
    pd.read_csv("./data/genetic_data/quality_control/14_pop_strat/01_pca/loop_maf_missing_2_ldprunned_autosome_pca_not_outliers.fam", sep=" ", header=None, low_memory=False).shape[0] != \
    pca_outliers.shape[0]):
        raise ValueError("ERROR: FALSE! WE HAVE NOT CORRECLTY REMOVED PCA OUTLIERS ")

#remove now from the general dataset
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir ./14_pop_strat/03_ancestry_cleaned_plink;\
    plink \
        --bfile ./12_loop_maf_missing/loop_maf_missing_2 \
        --remove ./14_pop_strat/01_pca/eigen_out/smartpca_outliers.txt \
        --make-bed \
        --out ./14_pop_strat/03_ancestry_cleaned_plink/loop_maf_missing_2_pca_not_outliers"
)
#check
if( \
    pd.read_csv("./data/genetic_data/quality_control/12_loop_maf_missing/loop_maf_missing_2.fam", sep=" ", header=None, low_memory=False).shape[0] - \
    pd.read_csv("./data/genetic_data/quality_control/14_pop_strat/03_ancestry_cleaned_plink/loop_maf_missing_2_pca_not_outliers.fam", sep=" ", header=None, low_memory=False).shape[0] != \
    pca_outliers.shape[0]):
        raise ValueError("ERROR: FALSE! WE HAVE NOT CORRECLTY REMOVED PCA OUTLIERS ")

print_text("Ancestry analysis following Admixture tutorial", header=4)
    #tutorial developed by the authors of admixture
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #admixture manual
        #http://dalexander.github.io/admixture/admixture-manual.pdf
print("general info")
#ADMIXTURE is a program for estimating ancestry in a model-based manner from large autosomal SNP genotype datasets, where the individuals are unrelated (for example, the individuals in a case-control association study).

#ADMIXTURE’s input is binary PLINK (.bed), ordinary PLINK (.ped), or EIGENSTRAT (.geno) formatted files and its output is simple space-delimited files containing the parameter estimates.

#To use ADMIXTURE, you need an input file and an idea of K, your belief of the number of ancestral populations. You should also have the associated support files alongside your main input file, in the same directory. For example, if your primary input file is a .bed file, you should have the associated .bim (binary marker information file) and .fam stub file) files in the same directory.
    #IMPORTANT: They say "your belief of the number of ancestral populations". In other words, based on your previous knowledge of your sample, select reasonable numbers for K. The literaly say "if you believe that the individuals in the sample derive their ancestry from three ancestral populations then run admixture .... using k=3". We have prior information from the PCA that we have a single ancestry, so we should stay not far away from 1.

#Prune SNPs? We tend to believe this is a good idea, since our model does not explicitly take LD into consideration, and since enormous data sets take more time to analyze. It is impossible to “remove” all LD, especially in recently-admixed populations, which have a high degree of “admixture LD”. Two approaches to mitigating the effects of LD are to include markers that are separated from each other by a certain genetic distance, or to thin the markers according the observed sample correlation coefficients. The easiest way is the latter, using the --indep-pairwise option of PLINK. For example, if we start with a file rawData.bed, we could use the following commands to prune according to a correlation threshold and store the pruned dataset in prunedData.bed

#SNPs in high LD with each other contain redundant information. More worrisome is the potential for some regions of the genome to have a disproportionate influence on the results and thus distort the representation of genome-wide structure. A nice empirical example of the problem is in figure 5 of Tian et al. [30], where PC2 of the genome-wide data is shown to be reflecting the variation in a 3.8 Mb region of chromosome 8 that is known to harbor an inversion. A standard approach to address this issue is to filter out SNPs based on pairwise LD to produce a reduced set of more independent markers. Here we use plink’s commands to produce a new LD-pruned dataset with output prefix H938_Euro.LDprune. The approach considers a chromosomal window of 50 SNPs at a time, and for any pair whose genotypes have an association r2 value greater than 0.1, it removes a SNP from the pair. Then the window is shifted by 10 SNPs and the procedure is repeated
    #We have removed SNPs correlated more than 0.2 instead of 0.1, but the we are looking at 500KB windows, very big windows. We are being strict according to the plink standards

#how many markers needed? This depends on how genetically differentiated your populations are, and on what you plan to do with the estimates. It has been noted elsewhere [4] that the number of markers needed to resolve populations in this kind of analysis is inversely proportional to the genetic distance (FST ) betweeen the populations. It is also noted in that paper that more markers are needed to perform adequate GWAS correction than are needed to simply observe the population structure. As a rule of thumb, we have found that 10,000 markers suffice to perform GWAS correction for continentally separated populations (for example, African, Asian, and European pop- ulations FST > .05) while more like 100,000 markers are necessary when the populations are within a continent (Europe, for instance, FST < 0.01). 

#The ADMIXTURE software (v 1.3.0 here) comes as a pre-compiled binary executable file for either Linux or Mac operating systems. To install, simply download the package and move the executable into your standard execution path (e.g. “/usr/local/bin” on many Linux systems). Once installed, it is straightforward to run ADMIXTURE with a fixed number of source populations, commonly denoted by K. 

#There is an output file for each parameter set: Q (the ancestry fractions), and P (the allele frequencies of the inferred ancestral populations). Note that the output filenames have ‘3’ in them. This indicates the number of populations (K) that was assumed for the analysis. This filename convention makes it easy to run analyses using different values of K in the same directory.

#ADMIXTURE is a maximum-likelihood based method, so as the method runs, you will see updates to the log-likelihood as it converges on a solution for the ancestry proportions and allele frequencies that maximize the likelihood function. The algorithm will stop when the difference between successive iterations is small (the “delta” value takes a small value). A final output is an estimated FST value [11] between each of the source populations, based on the inferred allele frequencies. These estimates reflect how differentiated the source populations are, which is important for understanding whether the population structure observed in a sample is substantial or not (values closer to 0 reflect less population differentiation).


print("run admixture")
#run admixture from K=1 to K=8 using cross-validation
run_bash("mkdir -p ./data/genetic_data/quality_control/14_pop_strat/02_admixture/admixture_cv")
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/02_admixture/; \
    prefix=loop_maf_missing_2_ldprunned_autosome_pca_not_outliers; \
    Klow=1; \
    Khigh=7; \
    for ((K=$Klow; K<=$Khigh; K++)); do \
        admixture \
            --seed 854534 \
            --cv=5 \
            -j28 \
            ../01_pca/${prefix}.bed \
            $K | tee log.${prefix}.${K}.out; \
        mv ./log.${prefix}.${K}.out ./admixture_cv/log.${prefix}.K${K}.out; \
        mv ./${prefix}.${K}.Q ./admixture_cv/${prefix}.K${K}.Q; \
        mv ./${prefix}.${K}.P ./admixture_cv/${prefix}.K${K}.P; \
    done; \
")
    #--seed
        #random seed for initialization
    #--cv
        #Use ADMIXTURE’s cross-validation procedure to select the best K value. A good value of K will exhibit a low cross-validation error compared to other K values. Cross-validation is enabled by simply adding the --cv flag to the ADMIXTURE command line. In this default setting, the cross-validation procedure will perform 5-fold CV—you can get 10-fold CV, for example, using --cv=10. The cross-validation error is reported in the output.
        #ADMIXTURE includes a cross-validation procedure that allows the user to identify the value of K for which the model has best predictive accuracy, as determined by “holding out” data points.
        #If I understand correctly, they create 5 folds with all genotypes, and convert to NA one of the folds. Then predic the genotypes of that fold and calculate the difference respect to the observed. The better prediction, the better a model with K ancestries works.
    #-j
        #To split ADMIXTURE’s work among N threads, you may append the flag -jN to your ADMIXTURE command. The core algorithms will run up to N times as fast, presuming you have at least N processors
    #the last two positional arguments are
        #input file
            #a PLINK .bed file
        #K
            #the number of populations;
    #The tee command writes the output to the file log.${prefix}.${K}.out and simultaneously displays it on the terminal.

#Note about Q files:
    #The Q estimates are output as a simple matrix, so it is easy to make figures like Figure 1 in Admixture paper
    #tbl_raw = pd.read_csv("./data/genetic_data/quality_control/14_pop_strat/02_admixture/admixture_cv/loop_maf_missing_2_ldprunned_autosome_pca_not_outliers.K2.Q", sep="\s+", header=None)

#Results
    #The CV-error tend to increase from K=1 onwards (see plots below). Also the FST values are relatively low, around 0.01-0.02. The only increase to 0.04 when considering 7 groups. According to acenstry docs, continental differentiated pops have FSTs values>0.05, so it seems we are dealing with 1 single big ancestry group according to admixture.
        #continentally separated populations (for example, African, Asian, and European pop- ulations FST > .05) while more like 100,000 markers are necessary when the populations are within a continent (Europe, for instance, FST < 0.01).
    #If in a revision someone ask about the ancestry, we say that according the PCAs and the admixture analyses, we only have 1 ancestry in our sample after removing PCA outliers. The ancestry is possibly European (after comparing a few allele frequencies with ncbi freqs across pops), but we do not know for sure.

print("extract the CV errors")
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/02_admixture/; \
    Klow=1; \
    Khigh=7; \
    for ((K=$Klow; K<=$Khigh; K++)); do \
        if [ $K -eq 1 ]; then \
            > cv_errors.tsv; \
        fi; \
        awk \
            -v K=$K \
            'BEGIN{FS=\" \"; OFS=\"\t\"}{ \
                if (K==1 && NR == 1){ \
                    print \"k\", \"cv_error\"\
                }; \
                if ($0 ~ /^CV error/){ \
                    print K, $4; \
                }; \
            }' \
            ./admixture_cv/log.loop_maf_missing_2_ldprunned_autosome_pca_not_outliers.K$K.out >> cv_errors.tsv ; \
    done; \
")
    #From K=1 to K=7
        #if we are analyzing the first K value (K=1), then create a new file to store results. Using ">" makes it possible to overwrite if the file exists
        #with Awk
            #add K as a value and use spaces as input FS and tabs as output FS
            #if the analyzing the first file and its ifrst row, add the header as two strings
            #if the row starts with "CV error", then print two columns with K and the fourth field of that row, i.e., the CV error
            #add the output in the TSV file

print("plot the CV errors")
#load the data
cv_errors = pd.read_csv( \
    "./data/genetic_data/quality_control/14_pop_strat/02_admixture/cv_errors.tsv", \
    sep="\t", \
    header=0, \
    low_memory=False)

#create a scatter plot with connected dots
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.plot( \
    cv_errors["k"], \
    cv_errors["cv_error"], \
    marker="o", \
    linestyle="-")
plt.xlabel("K values")
plt.ylabel('Cross-validation error')
plt.title('Change of CV error across different K values')
plt.savefig( \
    fname="./data/genetic_data/quality_control/14_pop_strat/02_admixture/cv_errors.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()

print("evaluate the K1 results in more detail")
#evalAdmix allows to evaluate the results of an admixture analysis (i.e. the result of applying ADMIXTURE, STRUCTURE, NGSadmix and similar). It only needs the input genotype data used for the previous admixture analysis and the output of that analysis (admixture proportions and ancestral population frequencies). The genotype input data can either be called genotypes in binary plink format or genotype likelihoods in beagle format. 
    #http://popgen.dk/software/index.php/EvalAdmix

#The output is a pairwise correlation of residuals matrix between individuals. The correlation will be close to 0 in case of a good fit of the data to the admixture model. When individuals do not fit the model, individuals with similar demographic histories (i.e. usually individuals from the same population) will be positively correlated; and individuals with different histories but that are modelled as sharing one or more ancestral populations as admixture sources will have a negative correlation. Positive correlation between a pair of individuals might also be due to relatedness. 
    #I understand that if the ancestry model is not working correctly, we are leaving unexplained variance, i.e., residuals, about the relationship between samples and that makes them to correlate between them.

#run evalAdmix
run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/; \
    mkdir -p ./02_admixture/eval_admix; \
    evalAdmix \
        -plink ./01_pca/loop_maf_missing_2_ldprunned_autosome_pca_not_outliers \
        -fname ./02_admixture/admixture_cv/loop_maf_missing_2_ldprunned_autosome_pca_not_outliers.K1.P \
        -qname ./02_admixture/admixture_cv/loop_maf_missing_2_ldprunned_autosome_pca_not_outliers.K1.Q \
        -o ./02_admixture/eval_admix/evaladmixOut.K1.corres \
        -P 20 \
")
    #plink: binary plink file prefix with genotype data
    #fname: file with ancestral frequencies (space delimited, rows are sites and columns ancestral populations), so you have 1 column for K=1, 3 columns for K=3 and so on... while the number of rows is the number of SNPs used.
    #qname: file with admixture proportions (space delimited, rows are individuals and columns ancestral populations), so the number of columns is again the number of pops considered, while the number of rows is the number of samples.
    #o: prefix of output file names
    #P: Number of threads used
    #example: ./evalAdmix -plink inputPlinkPrefix -fname inputPlinkPrefix.K.P -qname inputPlinkPrefix.K.Q -o evaladmixOut.corres -P 10 

#output:
    #The analysis performed by evalAdmix produces one file, containing a tab delimited N times N symmetric correlation matrix, where column i in line j contains the correlation of residuals between individual i and j, and the diagonal values (self-correlation) are set to NA:
        #NA 0.008609 -0.006919 0.002731 0.020224
        #0.008609 NA 0.000033 0.004968 -0.008470
        #-0.006919 0.000033 NA 0.006982 0.005664
        #0.002731 0.004968 0.006982 NA 0.000521
        #0.020224 -0.008470 0.005664 0.000521 NA 

#plot the results
#I am using as reference the R code showed in the manual
    #http://popgen.dk/software/index.php/EvalAdmix
run_bash("R CMD BATCH /opt/scripts/02bb_pre_imputation_qc_eval_admix_plot.R")
    
#Results:
    #see in figure "evaladmix_plot.pdf" and compare with the figures in the evalAdmix manual
        #http://popgen.dk/software/index.php/EvalAdmix#Run_command_example
    #The plot show almost no color, like the right plot in the example. This means that there is no correlation between the residuals of the samples under the admixture model of 1 population.
    #This is completely different compared to the left plot.
    #Therefore, we can confidently say that an Admixture model assuming just 1 ancestry group fits relavitvely well our data. More evidence about the fact that our sample is already genetically homogeneous.

print_text("Summary Ancestry: Results from the PCA analyses with smartpca and the admixture analyses with Admixture strongly suggest that our sample is relatively homogeneous in genetic terms. There are not differences at the continental level, e.g., Africans vs Europeans. Sure, we can still have more subtle differences within our sample, but we are going to use the main PCAs in the GWAS analysis anyways", header=4)

# endregion






##########################
# region SEX CHECK #######
##########################
print_text("start check sex", header=2)
print_text("open folder", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./15_check_sex; \
    ls -l")

#If you have X-chromosome population-genomic data, you can employ PLINK’s –check-sex command to sanity-check the sex information in your .fam file. (The method is based on chrX heterozygosity rates.) Similarly, –impute-sex uses chrX heterozygosity rates to fill in missing sex entries when appropriate.

#If your species has pseudoautosomal regions, ensure they are not encoded as part of the X chromosome. If they are, PLINK 1.9’s –split-x or PLINK 2’s –split-par flag can be used to change the chromosome codes of the relevant variants before proceeding.
    #Yes, we have PAR regions separated as a different chromosome

print_text("LD prunning", header=3)
#we need good MAF and linkage equilibrum, so we are running in the LD subset on the full dataset after removing the ancestry outliers. We cannot use the previous LD subset because it only had autosomals, and we need sex chromosomes to do the sex check.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3
    #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
#prune LD, retainn sex
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./15_check_sex/00_ld_prunning/; \
    plink2 \
        --bfile ./14_pop_strat/03_ancestry_cleaned_plink/loop_maf_missing_2_pca_not_outliers \
        --indep-pairwise 500kb 1 0.2 \
        --out ./15_check_sex/00_ld_prunning/ldpruned; \
    ls -l ./15_check_sex/00_ld_prunning/")
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./14_pop_strat/03_ancestry_cleaned_plink/loop_maf_missing_2_pca_not_outliers \
        --extract ./15_check_sex/00_ld_prunning/ldpruned.prune.in \
        --make-bed \
        --out ./15_check_sex/00_ld_prunning/loop_maf_missing_2_pca_not_outliers_ldpruned;\
    ls -l ./15_check_sex/00_ld_prunning/")

print_text("Check we have around 100K SNPs", header=4)
    #We have used this criteria for the LD prunning during the pop structure analyses, see the begining of that section for further details
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex/00_ld_prunning; \
    sex_ld_snps_in=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./loop_maf_missing_2_pca_not_outliers_ldpruned.bim); \
    printf 'The number of included autosomal SNPs in linkage equilibrium is %s\n' \"$sex_ld_snps_in\"; \
    echo 'Is this number greater than 93K?'; \
    if [[ $sex_ld_snps_in -gt 93000 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #you can print a variable with text using printf and %s for "string". Then you call the variable you want to add within "". You could use two variables: "$var1 $var2"
            #https://phoenixnap.com/kb/bash-printf

 
print_text("Run first check sex, mostly based on F", header=3)
#Run –check-sex once without additional parameters, just to see the distribution of F (inbreeding) coefficients.
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    mkdir -p ./01_first_check_sex; \
    plink \
        --bfile ./00_ld_prunning/loop_maf_missing_2_pca_not_outliers_ldpruned \
        --check-sex \
        --out ./01_first_check_sex/f_distribution \
")
    #--check-sex normally compares sex assignments in the input dataset with those imputed from X chromosome inbreeding coefficients, and writes a report to plink.sexcheck. So you estimate sex from the genetic data.
        #It uses the F statistic as calcualted for --het/--ibc
            #It is based on the observed and expected autosomal homozygous genotype counts for each sample. The method-of-moments F coefficient is calculated as follows: (<observed hom. count> - <expected count>) / (<total observations> - <expected count>))
            #Expected counts are based on loaded (via --read-freq) or imputed MAFs; if there are very few samples in your immediate fileset, --read-freq is practically mandatory since imputed MAFs are wildly inaccurate in that case.
                #Not our case, we are using the complete dataset after filtering (see below about this matter) 
            #Also, due to the use of allele frequencies, if your dataset has a highly imbalanced ancestry distribution (e.g. >90% EUR but a few samples with ancestry primarily from other continents), you may need to process the rare-ancestry samples separately.
                #Again no problem as our sample is already homogeneous regarding ancestry (see previous section)
            #By default, the n/(n-1) multiplier in Nei's expected homozygosity formula is now omitted, since n may be unknown when using --read-freq. The 'small-sample' modifier causes the multiplier to be included, while forcing --het to use.
                #I guess if you use the MAFs from the whole dataset when you are dealing with a very small subset (again see below of why is this important), then you do not know the sample size of the MAF used is not known, as you are not using the MAF of the subset.
            #Disssect the F formula:
                #Remember that in the case of --check-sex F is calculated ONLY with genotypes of the X chromosome: "--check-sex normally compares sex assignments in the input dataset with those imputed from X chromosome inbreeding coefficients"
                #F is (<observed hom. count> - <expected count>) / (<total observations> - <expected count>)
                    #Observed Homozygosity: PLINK calculates the observed number of homozygous genotypes (e.g., AA or aa) for each sample
                    #Expected Homozygosity: It also calculates the expected number of homozygous genotypes under the assumption of random mating and no population structure
                        #Expected homozygosity is calculated with the formula "p^2 + q^2" which is derived from the Hardy-Weinberg equilibrium. As you remember, this states that in a large, randomly-mating population with no external influences (like mutation, migration, or selection), the genotype frequencies remain constant from generation to generation. Essentially, p^2 is the expected frequency of the homozygous genotype for one allele, and q^2 for the other allele.
                    #“Total observations” in the F coefficient formula refers to the total number of genotyped loci, or positions on the genome, being analyzed for each individual. It’s the sum of all observed genotypes, both homozygous and heterozygous. 
                        #In the case of --check-sex, this should the total number of genotypes in the X chromosome.
                        #Think of it as the grand tally of all genotyped data points used to calculate the homozygosity for a sample. This helps to normalize the difference between observed and expected homozygosity, providing a standardized measure of inbreeding or population structure.
                        #you calculate the difference between observed and expected homozygosity and then normalize by the total number of genotypes
                        #The idea is to measure how far the observed data deviates from what we'd expect in a perfectly random mating scenario and adjust it according to the overall amount of genotyped data points. This gives you a standardized F coefficient to assess genetic relatedness or population structure
                #As you will see below, males are expected to have a F value around 1 while females are expected to have a value around 0.
                    #Males:
                        #As males only have one X chromosome, they have only one a allele for each genotype in the X chromosome
                        #This means almost all loci on their X chromosome are effectively homozygous, leading to an F coefficient close to 1.
                        #In "(<observed hom. count> - <expected count>) / (<total observations> - <expected count>)", "expected count" is the same in numerator and denominator, thus if F=1, it means that "observed hom. count" and "total observations". In other words, the total number of genotypes in the X chromosome is identical to the total number of homozygous genotypes, all genotypes are homozygous because male only have one X chromosome.
                    #Females:
                        #Females have two X chromosomes, which means they can have both homozygous and heterozygous loci. The presence of heterozygous loci lowers the F coefficient, leading it to be closer to 0.
                        #As the number of homozygous genotypes decreases, <observed hom. count> - <expected count> becomes smaller because the count of homozygous is smaller but the exepected homozygous are the same. Remember the expaction is based on the allele frequencies which are calculated at a population level.
                        #When the number of homo genotypes is identical to the expected, we get 0 in the denominator, and hence F=0 (0 is the expectation). If homo genotypes is even smaller than the expectation, then the denominator is negative, making F value negative. It will get more negative as homo genotype are much more lower than the expectation. 
                        #This explain why females have F values around zero. As they have two X chromosomes, they can have homo but also heterozygous genotypes, meaning that they can be around the expectation, i.e,. around 0.
                    #Because all of this, we should have a tight peak around 1 for males (all homozigous genotypes with no dispersion, male=no hetero calls in the X), while the females can be around 0 with some dispersion as they can have different proportions of homo/hetero (they have two X chromosomes) but always around 0, i.e., around the expectation.
                #We need to see the empirical distribution of F in our data to detect both peaks, and then decide the best threshold to separate males and females from according to F.
        #Make sure that the X chromosome pseudo-autosomal region has been split off (with e.g. --split-x) before using this.
            #we have pseudo-autosomal as 25, and I cleaned these SNPs ensuring that SNPs with code 25 are only those inside the known PAR regions
        #By default, F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male calls. If you pass numeric parameter(s) to --check-sex (without 'y-only'), the first two control these thresholds.
            #this is useful to apply a new threshold after we have seen the actual distribution of F in our data.
        #Since this function is based on the same F coefficient as --het/--ibc, it requires reasonable MAF estimates (so it's essential to use --read-freq if there are very few samples in your immediate fileset), and it's best used on marker sets in approximate linkage equilibrium.
            #Imagine you are using just a few samples from the whole dataset and checking their sex separately. These samples are included in a separated fileset. If you do check-sex, it will calculate the F coefficient using the MAF of the SNPs in these few samples instead of considering the MAF for the whole sample. This could be misleading. Independently of the genotype call for these few individuals, the total minor allele frequency is not that of them but of the whole population. So plink docs recommend to use --read-freq to upload a file with the actual MAFs from the whole population. THIS IS NOT a problem for us as we are using the whole set of samples. 
            #We have filtered by SNPs in this steps because we need SNPs in linkage equilibrium, but removing a SNPs does not affect the MAF of other SNP.
        #As a concrete example, if you ignore the Y chromosome and don't first perform LD-based pruning, --check-sex makes flat-out incorrect calls on several 1000 Genomes phase 1 female samples: NA19332 has an F estimate just under 0.88, and there are others with F > 0.8. However, if you first run "--indep-pairphase 20000 2000 0.5", the largest female F estimate drops to about 0.66. 0.66 is, of course, still much larger than 0.2, and in most contexts it still justifies a female call. We suggest running --check-sex once without parameters, eyeballing the distribution of F estimates (there should be a clear gap between a very tight male clump at the right side of the distribution and the females everywhere else), and then rerunning with parameters corresponding to the empirical gap.
        #Due to the use of allele frequencies, if your dataset has a highly imbalanced ancestry distribution, you may need to process the rare-ancestry samples separately.
            #Our data is already homogeneus in terms of ancestry, we already remove PCA outliers.

print_text("check the F distribution", header=4)
#We suggest running --check-sex once without parameters, eyeballing the distribution of F estimates (there should be a clear gap between a very tight male clump at the right side of the distribution and the females everywhere else), and then rerunning with parameters corresponding to the empirical gap.
    #https://www.cog-genomics.org/plink/1.9/basic_stats
f_distribution = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/01_first_check_sex/f_distribution.sexcheck", \
    sep="\s+", \
    header=0, \
    low_memory=False)
        #sep
            #The separator sep="\s+" in the highlighted code is a regular expression that matches one or more whitespace characters. This means that the read_csv function will use any sequence of whitespace (spaces, tabs, etc.) as the delimiter to separate columns in the input file.
            #Using tabs or just space does not work.
        #.sexcheck: A text file with a header line, and then one line per sample with the following 6-7 fields:
            #FID	Family ID
            #IID	Within-family ID
            #PEDSEX	Sex code in input file (1 = male, 2 = female, 0 = unknown according the fam file)
            #SNPSEX	Imputed sex code (1 = male, 2 = female, 0 = unknown)
            #STATUS	'OK' if PEDSEX and SNPSEX match and are nonzero, 'PROBLEM' otherwise
            #F	Inbreeding coefficient, considering only X chromosome. Not present with 'y-only'.
                #we are using "y-count", so we still use X data, and we should have the F coefficient

print_text("explore the F distribution", header=4)
print("make a plot of the F distribution")
import matplotlib.pyplot as plt
plt.hist(f_distribution.iloc[:, 5], bins=100, edgecolor='black')
plt.title('Distribution of the F statistic')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.savefig( \
    fname="./data/genetic_data/quality_control/15_check_sex/01_first_check_sex/f_distribution.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()
    #If both genders are well-represented in the dataset, you should see a big tight clump near 1 (corresponding to the males), and a more widely dispersed set of values centered near 0 (corresponding to the females).
    #Exactly what we have, except a case with F just below 0.8, but it is only 1.
print("Sample with an F value between 0.2 and 0.8")
samples_between_f = f_distribution.loc[(f_distribution.iloc[:, 5]>0.2) & (f_distribution.iloc[:, 5]<0.8),:]
print(samples_between_f)
if(samples_between_f.shape[0]!=1):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM, WE SHOULD HAVE ONLY 1 SAMPLE WITH A F VALUE BETWEEN 0.2 AND 0.8, BUT THIS IS NOT THE CASE")
    #So we are going to continue using the default threshold values for the F statisic, i.e., 0.

#Initial exploration
    #If you change the thresholds, the F values are still going to be the same, what changes is the STATUS column, i.e,. if a sample is consider FEMALE, MALE or UNKNOWN
    #The docs indicates "--check-sex ycount [female max F] [male min F] [female max Y obs] [male min Y obs]". Therefore, after "--check-sex ycount", the first two numbers indicate the F threshold used to impute female/male calls and the next two are used to control the max number of Y genotype calls are tolerated for females and the minimum for males.
    #We are using default values for F and Y thresholds:
        #By default, F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male calls. If you pass numeric parameter(s) to --check-sex (without 'y-only'), the first two control these thresholds.
        #In 'ycount' mode, gender is still imputed from the X chromosome, but female calls are downgraded to ambiguous whenever more than 0 nonmissing Y genotypes are present, and male calls are downgraded when fewer than 0 are present. (Note that these are counts, not rates.) These thresholds are controllable with --check-sex ycount's optional 3rd and 4th numeric parameters.
        #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
    #Therefore, as we are going to still use eveything as it is, our previous analysis is like "--check-sex 0.2 0.8"

print_text("Explore the rest of results, i.e., STATUS column", header=4)
print("see samples with STATUS==PROBLEM")
import numpy as np
f'We have {np.sum(f_distribution["STATUS"]=="PROBLEM")} problematic cases'

print("count samples that have unknown sex in the fam file")
f_distribution_unknown_fam = f_distribution.loc[f_distribution["PEDSEX"]==0,:]
f'We have {f_distribution_unknown_fam.shape[0]} cases without SEX in the fam file'

print("create list with samples without pheno in the excel file ('Combat gene DNA GWAS 23062022_v2.xlsx') and check their status")    
#load the excel file
pheno_data = pd.read_excel(
    "./data/pheno_data/Combat gene DNA GWAS 23062022_v2.xlsx",
    header=0,
    sheet_name="All DNA samples")

#extract the IDs of those samples without pheno (Gender) data
samples_without_pheno = pheno_data.loc[ \
    pheno_data["AGRF code"].notna() & pheno_data["Gender"].isna(), \
    "AGRF code"].to_list()

#extract their data from the check sex file
samples_no_pheno = f_distribution.loc[f_distribution["IID"].isin(samples_without_pheno), :]
print(samples_no_pheno)

#see STATUS of these samples
samples_no_pheno_status = samples_no_pheno["STATUS"]
if( \
    (len(samples_no_pheno_status.unique())==1) & \
    (samples_no_pheno_status.unique()=="PROBLEM")[0] \
):
    print(f"ALL SAMPLES WITHOUT PHENO DATA ARE COUNTED AS PROBLEMATIC, being {samples_no_pheno_status.shape[0]} in total")
else:
    raise ValueError("ERROR! FALSE! There are samples without pheno data that are not considered problematic")

print("there is 1 sample with unknown sex in plink that it is NOT included in the list of samples that do not have phenotype in the excel, because it is indeed a sample not included in the excel at all:")
print(f_distribution_unknown_fam.loc[~f_distribution_unknown_fam["IID"].isin(samples_without_pheno)])
    #"~" is used to negate
print("this was the case that is not present in the excel but there is a sample in the excel with a very similar ID: 2399LDJA insted of 2397LDJA")
if(f_distribution_unknown_fam.loc[~f_distribution_unknown_fam["IID"].isin(samples_without_pheno),"IID"].to_numpy()[0] != "2397LDJA"):
    raise ValueError("ERROR! FALSE! PROBLEM WITH THE SAMPLE 2397LDJA")

print("Add samples without pheno in excel to be removed")
samples_remove_first_round = samples_no_pheno["IID"].to_list()
if(len(samples_remove_first_round)!=35):
    raise ValueError("ERROR! FALSE! THE NUMBER OF SAMPLES TO REMOVE IN THE FIRST SEX CLEANING ROUND IS NOT OK")

print("see now cases that are not included in the list of samples without pheno but still have a problem in sex. These are the real problems because if we have sex data and still have an error, it means there is a mismatch between self-reported sex and imputed sex using genetics genetics")
real_sex_problems = f_distribution.loc[(~f_distribution["IID"].isin(samples_without_pheno)) & (f_distribution["STATUS"]=="PROBLEM"), :]
print("Number of real problems: " + str(real_sex_problems.shape[0]))
if(real_sex_problems.shape[0]!=5):
    raise ValueError("ERROR! FALSE! THE NUMBER OF REAL SEX PROBLEMS IS NOT OK")

print("check we have the correct number of problematic cases")
#concatenate the problematic cases due to the lack of SEX in the pheno data PLUS the real problematic cases, i.e., cases with SEX reported but with a mismatch with the genetics
combined_problematic = pd.concat([f_distribution_unknown_fam["IID"], real_sex_problems["IID"]], ignore_index=True)

#count the number of duplicates, 1 means that 2 IDs are identical
n_duplicates = combined_problematic.duplicated().sum()

#the sum of problematic cases minus 1 duplicate should give the total number of problematic cases
if(combined_problematic.shape[0] - n_duplicates != f_distribution.loc[(f_distribution["STATUS"]=="PROBLEM"), :].shape[0]):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM, THE NUMBER OF PROBLEMATIC CASES RELATED TO SEX DO NOT MATCH")

#Results:
with pd.option_context("display.max_columns", None):
    print(real_sex_problems)
    #This will display all columns of the real_sex_problems DataFrame for this single print operation without changing the global pandas settings.
#5 sex-problematic cases, the rest of problems are caused because they do not have pheno data at all

#there are three cases (7699ISMO, 8244GBJJ, 8702EBMF) with very low F and no Y genotypes that are consdered male by the pheno data.
#8702EBMF was identified as male, but it has very low F, no Y genotypes and its name is one of females, so it is likely female.
#8244GBJJ is the same case but just one Y genotype. As we will discuss below, just 1 Y genotype is not enough reason to remove a sample that is ok in the rest of aspects
#We retain BOTH BUT CHANGING THE SEX TO FEMALE.
samples_change_sex = ["8244GBJJ", "8702EBMF"]
    #We will run Y genotype counts later, but we are considering that here. If, for any reason, we are wrong, the script would eventually fail because these samples are going to become females, so we would have female samples with Y genotypes, this will create a PROBLEM when running -check-sex y-count
#7699ISMO has F value close to 0, no Y genotypes, but it has a male name so this is really strange. It seems a biological female with a male name. We are going to REMOVE.
samples_remove_first_round.append("7699ISMO")

#A self-reported male (8500JADJ) and is considered unknow by plink because his F value is below 0.8. It has 64 Y genotypes, a male name, but the problem is that the F values is far below 0.99 (expected for males) and below our threshold, it is not included in any of the F peaks. This sample has lower homozygosity in the X than expected for a male, meaning that if could have two alleles for several positions in the X. We are REMOVING.
samples_remove_first_round.append("8500JADJ")

#one case (2397LDJA) do not have sex, but genetic data is very clear: F=0.99 and 64 Y genotypes, plus male name. This strongly suggesting male. CHANGE TO MALE.
    #this was the cases not present in the excel. There is a sample in the excel with a very similar ID (2399LDJA). Lo an behold, that sample in the excel is reported as a male. So, as David´s postdoc suggested, this is likely a mislabelled sample. We are going to set the sex as male and change the ID to the one in the excel file.
samples_change_sex.append("2397LDJA")

print("check we have the correct number of samples to remove")
if(len(samples_remove_first_round)+len(samples_change_sex)!=f_distribution.loc[f_distribution["STATUS"]=="PROBLEM", :].shape[0]):
    raise ValueError("ERROR! FALSE! THE NUMBER OF SAMPLES TO REMOVE IS NOT OK")
else:
    print("OK")

#Aside note about names
#I have checked the names of the 13 samples checked by David (sex_mismatches_bishop.xlsx) to see if some names were showing a very different ancestry. We should have similar ancestries given we have already filtered the samples by population structure.
#Most of the first and last names are of European origin, suggesting that our sample is possibly of European ancestry.
#There are two cases with a hispanic first name, but it seems these names are also used for non-hispanic people. Also note that the last names are English. So very unlikely these are people from the Americas or Phillipines. Note that Hispanic-americans and Filipinos usually have both first and last names from hispanic origin.
#There are two cases with Turkish names. One have a last name that is of Turkish origin and the other has a last name that is frequently used in Arabic-speaking and Muslim countries but also very frequent in Turkey probably due to historic (associated with Carlos). From what I have seen, it seems these persons are possibly of turkish origin, which is not far off from European ancestry.
#Therefore, in general, it seems that our filtering has possibly removed samples that are far away from European ancestries.

print_text("remove samples for this step", header=4)
print("load the fam file")
loop_maf_missing_2_pca_not_outliers_fam = pd.read_csv( \
    "./data/genetic_data/quality_control/14_pop_strat/03_ancestry_cleaned_plink/loop_maf_missing_2_pca_not_outliers.fam",
    sep=" ", \
    header=None, \
    low_memory=False \
)
print(loop_maf_missing_2_pca_not_outliers_fam)

print("extract family and individual IDs to be removed")
id_samples_first_removal = loop_maf_missing_2_pca_not_outliers_fam.loc[ \
    loop_maf_missing_2_pca_not_outliers_fam.iloc[:,1].isin(samples_remove_first_round), \
    [0, 1] \
]
#check
if(len(samples_remove_first_round)!=id_samples_first_removal.shape[0]):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE FIRST REMOVAL OF SAMPLES")

#save
id_samples_first_removal.to_csv( \
    "./data/genetic_data/quality_control/15_check_sex/01_first_check_sex/id_samples_first_removal.tsv",
    sep="\t",
    header=False,
    index=False \
)

print("remove")
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./14_pop_strat/03_ancestry_cleaned_plink/loop_maf_missing_2_pca_not_outliers \
        --remove ./15_check_sex/01_first_check_sex/id_samples_first_removal.tsv \
        --make-bed \
        --out ./15_check_sex/01_first_check_sex/loop_maf_missing_2_pca_not_outliers_first_sex_removal"
)
    #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.

print("check we removed the correct samples")
first_sex_clean = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/01_first_check_sex/loop_maf_missing_2_pca_not_outliers_first_sex_removal.fam", \
    sep=" ", \
    header=None, \
    low_memory=False \
)
if( \
    (first_sex_clean[1].isin(id_samples_first_removal[1]).sum()!=0) | \
    (first_sex_clean.shape[0] + id_samples_first_removal.shape[0] != f_distribution.shape[0]) \
):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE FIRST SEX REMOVAL OF SAMPLES")    


print_text("change sex and update IDs", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    mkdir -p ./02_sex_change_id_update; \
")

print_text("change sex samples", header=4)
print("load the fam file")
loop_maf_missing_2_pca_not_outliers_first_sex_removal_fam = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/01_first_check_sex/loop_maf_missing_2_pca_not_outliers_first_sex_removal.fam",
    sep=" ", \
    header=None, \
    low_memory=False \
)

print("extract family and individual IDs to be changed")
id_samples_sex_change = loop_maf_missing_2_pca_not_outliers_first_sex_removal_fam.loc[ \
    loop_maf_missing_2_pca_not_outliers_first_sex_removal_fam.iloc[:,1].isin(samples_change_sex), \
    [0, 1] \
]
#check we have selected the correct samples
if( \
    (id_samples_sex_change[1].isin(samples_change_sex).sum()!=id_samples_sex_change.shape[0]) | \
    (pd.Series(samples_change_sex).isin(id_samples_sex_change[1]).sum()!=len(samples_change_sex)) \
):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE FIRST REMOVAL OF SAMPLES")

print("add the new sex to the DF")
id_samples_sex_change[2] = pd.Series(dtype='int')
id_samples_sex_change.loc[id_samples_sex_change.iloc[:,1]=="8244GBJJ",2] = 2
id_samples_sex_change.loc[id_samples_sex_change.iloc[:,1]=="8702EBMF",2] = 2
id_samples_sex_change.loc[id_samples_sex_change.iloc[:,1]=="2397LDJA",2] = 1
    #--update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or F = female, 0 = missing) in the (n+2)th column. If no second parameter is provided, n defaults to 1. It is frequently useful to set n=3, since sex defaults to the 5th column in .ped and .fam files.
    #See above about the decisions to change sex
id_samples_sex_change[2] = id_samples_sex_change[2].astype('int')
print(id_samples_sex_change)
id_samples_sex_change.to_csv( \
    "./data/genetic_data/quality_control/15_check_sex/02_sex_change_id_update/id_samples_sex_change.tsv",
    sep="\t",
    header=False,
    index=False \
)

#update sex
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./15_check_sex/01_first_check_sex/loop_maf_missing_2_pca_not_outliers_first_sex_removal \
        --update-sex ./15_check_sex/02_sex_change_id_update/id_samples_sex_change.tsv \
        --make-bed \
        --out ./15_check_sex/02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update"
)
    #--update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or F = female, 0 = missing) in the (n+2)th column. If no second parameter is provided, n defaults to 1. It is frequently useful to set n=3, since sex defaults to the 5th column in .ped and .fam files.

print_text("update IDs", header=4)
print("load the fam file")
loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update.fam",
    sep=" ", \
    header=None, \
    low_memory=False \
)

print("extract family and individual IDs to be changed")
id_samples_id_update = loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update.loc[ \
    loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update.iloc[:,1]=="2397LDJA", \
    [0, 1] \
]

print("add the new ID")
id_samples_id_update[2] = id_samples_id_update[0].to_numpy()[0]
id_samples_id_update[3] = "2399LDJA"
print("New ID: ")
print(id_samples_id_update)

print("save")
id_samples_id_update.to_csv( \
    "./data/genetic_data/quality_control/15_check_sex/02_sex_change_id_update/id_samples_id_update.tsv",
    sep="\t",
    header=False,
    index=False \
)

print("update the ID")
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./15_check_sex/02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update \
        --update-ids ./15_check_sex/02_sex_change_id_update/id_samples_id_update.tsv \
        --make-bed \
        --out ./15_check_sex/02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update"
)
    #--update-ids expects input with the following four fields:
        #Old family ID
        #Old within-family ID
        #New family ID
        #New within-family ID

print("check")
loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_fam = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update.fam",
    sep=" ", \
    header=None, \
    low_memory=False \
)
if( \
    (loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_fam.iloc[:,1]=="2397LDJA").sum()!=0 | \
    (loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_fam.iloc[:,1]=="2399LDJA").sum()!=1
):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM UPDATING THE ID OF 2397LDJA TO 2399LDJA")


print_text("check sex using Y chromosome", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    mkdir -p ./03_second_check_sex; \
")

print_text("run check-sex only with Y", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    plink \
        --bfile ./02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update \
        --check-sex y-only 0 1 \
        --out ./03_second_check_sex/check_sex_y_only \
")
    #--check-sex has now two modes which consider Y chromosome data.
    #In 'ycount' mode
        #gender is still imputed from the X chromosome, but female calls are downgraded to ambiguous whenever more than 0 nonmissing Y genotypes are present.
            #In other words, if a female call has 1 or more call in the Y genotype (i.e., non-missing), then it is downgraded to unknown. This makes sense because if a sample reported as female has a Y genotype, then it is female or male? Probably a sample to be removed.
        #male calls are downgraded when fewer than 0 are present.
            #A male can never have less than 0 Y genotype calls, thus a male is not going to be downgraded because of this. 
            #By default, plink does not downgraded males that have zero Y genotype calls.
            #It seems that for imputing sex from scratch, it is not a good idea to discard males without Y genotypes because older males seems to suffer a decrease of the chromosome in some cells. So it would not be a very useful marker (see y-only mode section). If some actual males suffer Y losses in some cells, you could wrongly assigned them to "unknown when they are actual males".
        #Note that these thresholds are counts, not rates. These thresholds are controllable with --check-sex ycount's optional 3rd and 4th numeric parameters.
            #The docs indicates "--check-sex ycount [female max F] [male min F] [female max Y obs] [male min Y obs]". Therefore, after "--check-sex ycount", the first two numbers indicate the F threshold used to impute female/male calls and the next two are used to control the max number of Y genotype calls are tolerated for females and the minimum for males.
        #We want the default behaviour of "y-count" as this let us to consider both X and Y chromosomes in order to impute sex, while avoiding problems with the mosaic-loss of the Y chromosome. 
    #In 'y-only' mode
        #gender is imputed from nonmissing Y genotype counts, and the X chromosome is ignored. 
        #The male minimum threshold defaults to 1 instead of zero in this case. 
            #Now males with no Y genotypes are downgraded
        #IMPORTANT: This is intended to recover previously-determined gender, rather than determine it from scratch, since Y chromosome data may be scarce for older males: see e.g. Dumanski JP et al. (2014) Smoking is associated with mosaic loss of chromosome Y.
            #I guess that if you already have your dataset clean where strange cases like males with no Y genotypes being removed, you can count Y genotype calls and consider as males those with a given number of genotype Y calls.
    #We are using the "y-only" mode:
        #We have already used the X chromosome to calculate the F statistic and check sex. Now, we are interested only in see the females samples with non-missing Y genotypes. In the unlikely case that a male sample has a F value around 1 (we have alredy filtered by that) and have less than 1 Y genotype (one thing is a decrease in Y genotype in older and other is to have 0), the male sample will be set as PROBLEM, so we are going to add a check about that later.
        #Given we are not calculating F, we do not good MAF falues nor LD. We are just counting the number of Y genotypes, so we are going to use the whole sample without LD prunning and see the number of Y genotypes in the samples across all Y SNPs.
            #Since this function is based on the same F coefficient as --het/--ibc, it requires reasonable MAF estimates (so it's essential to use --read-freq if there are very few samples in your immediate fileset), and it's best used on marker sets in approximate linkage equilibrium.

print_text("explore checksex file and non-missing Y genotypes in females", header=4)
print("load the checksex file")
check_sex_y_only = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/03_second_check_sex/check_sex_y_only.sexcheck", \
    sep="\s+", \
    header=0, \
    low_memory=False \
)
print(check_sex_y_only)
    #.sexcheck: A text file with a header line, and then one line per sample with the following 6-7 fields:
        #FID	Family ID
        #IID	Within-family ID
        #PEDSEX	Sex code in input file (1 = male, 2 = female, 0 = unknown according the fam file)
        #SNPSEX	Imputed sex code (1 = male, 2 = female, 0 = unknown)
        #STATUS	'OK' if PEDSEX and SNPSEX match and are nonzero, 'PROBLEM' otherwise
            #we are using "y-count", so non-zero y-counts for females results in SNPSEX=1 and hence PROBLEM.
        #YCOUNT	Number of nonmissing genotype calls on Y chromosome. Requires 'ycount'/'y-only'.

print("extract problematic cases and check that these are females")
check_sex_y_only_problems = check_sex_y_only.loc[check_sex_y_only.iloc[:,4]=="PROBLEM"]
if( \
    (check_sex_y_only_problems.loc[ \
        (check_sex_y_only_problems.iloc[:,2]==1) \
    ].shape[0] != 0) | \
    (check_sex_y_only_problems.loc[ \
        (check_sex_y_only_problems.iloc[:,2]==2) & \
        (check_sex_y_only_problems.iloc[:,5]>0) \
    ].shape[0] != check_sex_y_only_problems.shape[0]) \
):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SAMPLES ACCORDING TO Y CHROMOSOME")
else:
    print("OK")

print("check the original problematic cases due to non-missing Y genotypes are included in this new list")
original_samples_y_problem =  pd.Series(["0295AMSM", "7692EOOO", "1390JMJM", "2197JWDM", "2282SODJ", "3400ISOM", "6796HGJS", "7300ECNO"])
    #Self-reported females, with F value around 0, and female names, but having 1 Y genotype.
if( \
    original_samples_y_problem.isin(check_sex_y_only_problems.iloc[:,1]).sum() != original_samples_y_problem.shape[0] \
):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SAMPLES ACCORDING TO Y CHROMOSOME")
else:
    print("OK")

print_text("Explore the problematic cases in the HH file, i.e., heterzigous alls in regions of the X chrosomoe that are not PAR and hence should be homozygous and haploid, i.e., heterozygous haploids", header=4)
print("load the HH file")
check_sex_y_only_hh = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/03_second_check_sex/check_sex_y_only.hh", \
    sep="\s+", \
    header=None, \
    low_memory=False \
)

print("count the number of unique SNPs and samples impicates")
print("snps")
print(check_sex_y_only_hh.iloc[:,2].value_counts().describe())
print("most of the SNPs are repeated 1 time and no more (percentile 75 is 1)") 
print("samples")
print(check_sex_y_only_hh.iloc[:,1].value_counts().describe())
print("Some samples are repeated several times, but most of them just 2 times (percentile 75 is 2")
print("Therefore, there are no SNPs/samples very badly affected")

print("check that all affected samples are males")
samples_hh_problems = loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_fam.loc[loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_fam[1].isin(check_sex_y_only_hh.iloc[:,1]), :]
    #get the FAM data of the het. haploids
if((samples_hh_problems.iloc[:,4]==1).sum()!=samples_hh_problems.shape[0]):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SAMPLES ACCORDING TO HETEROZYGOUS HAPLOIDS SOME OF THE SAMPLES ARE NOT MALES")
else:
    print("OK")

print("check that all affected SNPs are in X plus 1 case in the Y")
#load bim file
loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_bim = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update.bim",
    sep="\t", \
    header=None, \
    low_memory=False \
)
snps_hh_problems = loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_bim.loc[loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update_bim[1].isin(check_sex_y_only_hh.iloc[:,2]), :]
    #get the chromosome of all SNPs included in the HH file
if( \
    ((snps_hh_problems.iloc[:,0]==23) | (snps_hh_problems.iloc[:,0]==24)).sum()!=snps_hh_problems.shape[0] \
):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SNPS ACCORDING TO HETEROZYGOUS HAPLOIDS: SOME SNPS ARE NOT IN THE X (we should only have 1 case in the Y chromosome)")
else:
    print("OK")


print("there is also 1 Y case in the HH file")
y_male_problem = pd.concat( \
    [ \
        snps_hh_problems.loc[snps_hh_problems.iloc[:,0]==24,:].reset_index(drop=True), \
        check_sex_y_only_hh.loc[ \
            check_sex_y_only_hh.iloc[:,2].isin(snps_hh_problems.loc[snps_hh_problems.iloc[:,0]==24,1]), : \
        ].reset_index(drop=True) \
    ], \
    axis=1)
    #concatenate columns of: 
        #the ID of the problematic HH snp of chromosome Y 
        #the Id of the sample implicated in the problem
print(y_male_problem)
    #This is self-reported male (1198LKSM) with F=1 that has two different alleles for a SNP located in the non-PAR region of the Y chromosome.
    #combat_ILGSA24-17873  1198LKSM  JHU_Y.24444621
    #24  JHU_Y.24444621  0  22,298,475  C  T
    #this is the same than the rest of males that have heterzygous calls in non-PAR regions of the X chromosome, so the same logic applies.
if( \
    (y_male_problem.iloc[:,1].to_numpy()[0]!="JHU_Y.24444621") | \
    (y_male_problem.iloc[:,8].to_numpy()[0]!="JHU_Y.24444621") | \
    (y_male_problem.iloc[:,7].to_numpy()[0]!="1198LKSM") \
):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SNPS ACCORDING TO HETEROZYGOUS HAPLOIDS IN THE Y CHROMOSOME")

print("all cases (for X and Y) are outside or PAR regions?")
if( \
    snps_hh_problems[ \
        ((snps_hh_problems.iloc[:,3]>=10001) & (snps_hh_problems.iloc[:,3]<=2781479)) | \
        ((snps_hh_problems.iloc[:,3]>=155701383) & (snps_hh_problems.iloc[:,3]<=156030895)) \
    ].shape[0]!=0 \
):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SNPS ACCORDING TO HETEROZYGOUS HAPLOIDS: SOME SNPS ARE IN THE X PAR REGIONS")
else:
    print("OK")
    #Accoring to hg38 data (https://www.ncbi.nlm.nih.gov/grc/human), in chrX, the first pseudo-autosomal region starts at basepair 10001 and ends at basepair 2781479, being the latter the first boundary used by plink. The second region starts at basepair 155701383 and ends at basepair 156030895, being the former the second boundary indicated by Plink. In other words, plink considers the end of the first PAR region and the start of the second PAR region.
if( \
    ((y_male_problem.iloc[:,3][0]>=10001) & (y_male_problem.iloc[:,3][0]<=2781479)) | \
    ((y_male_problem.iloc[:,3][0]>=56887903) & (y_male_problem.iloc[:,3][0]<=57217415)) \
):
    raise ValueError("ERROR! FALSE! PROBLEM SELECTING THE PROBLEMATIC SNPS ACCORDING TO HETEROZYGOUS HAPLOIDS: SOME SNPS ARE IN THE Y PAR REGIONS")
else:
    print("OK")

print("If all the previous checka are OK, it means we have males that have heterozygous calls in non-PAR regions of the X or Y chromosome. We know these are males because we have already checked the sex with F and Y genotypes and know that SNPs in the PAR regions are clearly noted like that, thus these seems to be genotyping errors. We have to remove. The question is whether to just remove these specific genotypes or remove these SNPs for all samples altogether. We are going to remove only the specific genotypes implicated in this problem. See the script for further details about the decision")

#After removing the samples we decided to remove due to sex problems, I detected that plink was still giving me errors regarding the sex chromosomes.
    #These warnings from PLINK indicate issues with your genotype data:
        #1. **Heterozygous Haploid Genotypes**: This warning suggests that there are 650 instances where haploid genotypes (typically found on sex chromosomes) are showing heterozygous calls. This often happens with male samples on the X chromosome. You can check the `.hh` file mentioned for specific details. Using the `--split-x` command might help if these errors are in the pseudo-autosomal regions[1](https://www.biostars.org/p/9601535/).
        #2. **Nonmissing Nonmale Y Chromosome Genotypes**: This warning indicates that there are genotypes on the Y chromosome that are not missing in non-male samples. This could be due to incorrect sex assignments or genotyping errors. You might need to verify the sex information in your dataset and correct any discrepancies[2](https://www.biostars.org/p/98211/).
#As the calculation of the F statistic has to be performed in a subset of SNPs in linkage equilibrium, I did not count the number of Y genotypes for all SNPs present in the Y chromosome. When considering all SNPs in the Y chromosome, there are 35 additional female samples that have at least 1 genotype in that chromosome. This made me have second thoughts about removing female samples that are apparently females (F value around 0) but have a few Y genotypes. Even so after finding that this technical error is apparently fairly common.
    #From Christopher: If heterozygous haploid calls still remain, the most likely cause is nonmissing female genotype calls on the Y chromosome; others have reported that this is fairly common.  A quick way to check the number of these is to just load the Y chromosome with e.g. "plink --bfile semi_clean_fileset --chr 24 --freq".  If all the heterozygous haploid errors are on the Y chromosome, you can safely clobber them with --make-bed + --set-hh-missing.  (If some are on the X, --set-hh-missing *might* still be okay, but I'd need to know more about the data source and the --check-sex report to be sure.)
    #https://groups.google.com/g/plink2-users/c/4bpdLMdH2KA

#Another reason to give this a thought is the fact that we also have the opposite problem in some males. Some have a few heterozygous genotypes in the X chromosome, i.e., two different copies of the same SNP. Given that males only have one X chromosome, this is technically not possible. There are 285 males with this problem.
#Remember these (both male and female) are samples with the correct F value (0 for females and 1 for males), so the sample should be ok for the remaining genotypes. Considering this, I think that it would make sense to save all these samples and just remove the specific genotypes implicated in this.

#Christopher seems to concur with me:
    #I have multiple nonmissing female genotype calls on the Y chromosome. Assuming that these female samples have a F value around 0 (I checked that with --check-sex), the existence of a non-missing Y genotype would justify the removal of the whole sample? Or would it make sense to just remove that genotype but retain the rest of genotypes for these samples? As far as I know, this is a common problem, but I am not sure if this completely invalidates a sample.
    #Then, the .hh file include multiple male samples that have a few (between 1 to 9) heterozygous calls in non-PAR regions of the X chromosome. My question is the same, would it make sense to remove altogether these samples and their SNPs? Or could we just use --set-hh-missing to only remove these specific genotypes? Of course, this is assuming that these males have a F value above my threshold (i.e., 0.8).
    #In case it is relevant for the question in place: My .hh file includes the male samples with heterzygous haploid calls but not the nonmissing female genotype calls on the Y chromosome. However, "--set-hh-missing" remove all these genotypes, both the het. haploid and non-missing Y calls.
    #I am sure I am dealing with non-missing Y genotypes because I have checked the Y counts across all SNPs using "--check-sex y-only".
    #Chris: It’s a warning, not an error.  If you know your variant caller makes a few female chrY genotype calls, just use —set-hh-missing (which was renamed to —set-invalid-haploid-missing in recent plink 2.0 builds for the reason you describe).
        #It seems that "—set-hh-missing" has been changed to "—set-invalid-haploid-missing" because this command remove problematic haploid genotypes related to both heterozigous haploids in X males (non-PAR X regions that do not have counter part in the Y) and non-missing Y genotypes in females.

print_text("use --set-hh-missing to remove these problematic genotypes", header=4)
print("Run plink")
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    plink \
        --bfile ./02_sex_change_id_update/loop_maf_missing_2_pca_not_outliers_first_sex_removal_sex_update_id_update \
        --set-hh-missing \
        --make-bed \
        --out ./03_second_check_sex/loop_maf_missing_2_pca_not_outliers_sex_full_clean \
")
    #Normally, heterozygous haploid and nonmale Y chromosome genotype calls are logged to plink.hh (In my case nonmale Ys are not included in HH and Chris seems to be aware) and treated as missing by all analysis commands, but left undisturbed by --make-bed and --recode (since, once gender and/or chromosome code errors have been fixed, the calls are often valid). If you actually want --make-bed/--recode to erase this information, use --set-hh-missing. (The scope of this flag is a bit wider than for PLINK 1.07, since commands like --list and --recode-rlist which previously did not respect --set-hh-missing have been consolidated under --recode.)
    #Note that the most common source of heterozygous haploid errors is imported data which doesn't follow PLINK's convention for representing the X chromosome pseudo-autosomal region. This should be addressed with --split-x below, not --set-hh-missing.
        #https://www.cog-genomics.org/plink/1.9/data
    #We have already checked that the autosomal regions follow Plink´s convention and we know the cause of these problems, as previously disucssed we are going to remove only the genotypes implicated but not the rest of the data for the samples and snps implicated.

print("run again a sex check to ensure we do not have more nonmale Y and het. haploid cases")
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    plink \
        --bfile ./03_second_check_sex/loop_maf_missing_2_pca_not_outliers_sex_full_clean \
        --check-sex ycount \
        --out ./03_second_check_sex/last_check_sex \
")
    #IMPORTANT: Note that F values without LD prunning are not completely reliable, so you may have some PROBLEMs generated by this. We are not prunning SNPs because we need to count Y genotypes across all SNPs

print("check the log file that does not contain any warning anymore")
run_bash(" \
    awk \
        'BEGIN{FS=\"\t\"}{ \
            if($0 ~/^Warning/){ \
                print $0; \
                exit 1; \
            } \
        }' \
        ./data/genetic_data/quality_control/15_check_sex/03_second_check_sex/last_check_sex.log \
")
print("OK")

print("explore the last checksex file to do a few more checks")
#load the file
last_check_sex = pd.read_csv( \
    "./data/genetic_data/quality_control/15_check_sex/03_second_check_sex/last_check_sex.sexcheck", \
    sep="\s+", \
    header=0, \
    low_memory=False \
)

#calculate the checks
first_check = (not last_check_sex["STATUS"].isna().any()) & (last_check_sex["STATUS"]=="OK").all()
    #in the STATUS column
        #we do not have NAs
        #all cases are OK
second_check = last_check_sex.loc[last_check_sex["PEDSEX"]==2, "YCOUNT"].sum()==0
third_check = (last_check_sex.loc[last_check_sex["PEDSEX"]==2, "F"]<0.2).sum()==(last_check_sex["PEDSEX"]==2).sum()
    #check all females have
        #0 Y genotypes
        #F value below 0.2
fourth_check = (last_check_sex.loc[last_check_sex["PEDSEX"]==1, "YCOUNT"]>100).sum()==(last_check_sex["PEDSEX"]==1).sum()
    #check males have multiple Y genotypes
    #Most samples have more than 400 Y genotypes, there is only one smaple that has less, specifically 142. It is a self-reported male (F=1) without any other special characteristic.... Given the possibility of reduced Y genotypes in males (although this is more frequent in old people and this is a young) and the high F value (it is 1 also in the prunned dataset), we are going to maintain this sample. Also the X seems to be ok, and the Y is not going to be used after imputation, so...
        #last_check_sex.loc[(last_check_sex["PEDSEX"]==1) & (last_check_sex["YCOUNT"]<400),:]
fifth_check = (last_check_sex.loc[last_check_sex["PEDSEX"]==1, "F"]>0.8).sum()==(last_check_sex["PEDSEX"]==1).sum()
    #all males have an F value above 0.8

#apply the checks in an IF
if first_check & second_check & third_check & fourth_check & fifth_check:
    print("WE HAVE CORRECTLY CLEANED OUR DATA OF SEX ISSUES")
else:
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM CLEANING THE SEX PROBLEMS")
    #IMPORTANT: Note that F values without LD prunning are not completely reliable, so you may have some PROBLEMs generated by this. We are not prunning SNPs because we need to count Y genotypes across all SNPs

#If you want to see the excel files you sent to David, look at "/home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/quality_control/data/pheno_data/sex_mismatch_excels/"

#ANSWER DAVID AND JONATAN:
    #Jonatan: OK! go ahead from my side
    #David: I also agree with the suggested approach – it doesn’t affect the correlation and preserves as many samples as possible. My understanding after reading the Truong paper is that it is not necessary to eliminate individuals due to sex chromosome anomalies (which are expected in a proportion of the samples).
        #If David ask why you removed 2 problematic samples due to sex, you specify that these samples, these were problems in the F statistic for the X chromosome, i.e., not just a few genotypes but a problem with overall homozygosity in the X chromosome. Male with lower F than expected.

# endregion





##########################
# region HWE SECOND ROUND
###########################
print_text("start check sex", header=2)
print_text("open folder", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./16_hwe_second_round; \
")

#To remove variants calling errors, we initially applied a HWE filtering but it was one-sided. This means that SNPs with fewer than expected variantes are not removed. These desviations from the expectation would be expected in the case of population stratification, which we had before analyizing population structure. Therefore, we could remove legitime SNPs if we would apply this filter before considering population structure.
    #The “keep-fewhet” modifier causes this filter to be applied in a one-sided manner (so the fewer-hets-than-expected variants that one would expect from population stratification would not be filtered out by this command), and the 1e-25 threshold is extreme enough that we are unlikely to remove anything legitimate. (Unless the dataset is primarily composed of F1 hybrids from an artificial breeding program.) 
#After having a clear and homogeneous population, We apply a second round of HWE filtering. Within a population, violations of HWE fewer-hets-than-expected direction are also likely to be variant calling errors, so we are going to remove them.
    #After you have a good idea of population structure in your dataset, you may want to follow up with a round of two-sided –hwe filtering, since large (see Note 5) violations of Hardy–Weinberg equilibrium in the fewer-hets-than-expected direction within a subpopulation are also likely to be variant calling errors; with multiple subpopulations, the –write-snplist and –extract flags can help you keep just the SNPs which pass all subpopulation HWE filters.

#Regarding the threshold we are using now, Plink tutorial suggest using “ridiculous” high thresholds "1e-25 or 1e-50" for QC qhen having at least 1K samples, which is our case. Normal thresholds like 1e-4 should be used for specific analyses that assume HWE. Here we are just removing genotyping errors. WE USE 1e-25, LIKE IN THE PREVIOUS HWE STEP.
    #What p-value threshold does “large” correspond to, you may ask? Well, this depends on the size of your dataset and some other characteristics of your data. But a good rule of thumb is to use “ridiculous” thresholds like 1e-25 or 1e-50 for quality control when you have at least a thousand samples (or even 1e-200 if you have many more), while only using “normal” thresholds like 1e-4 or 1e-7 when you are preparing data for an analysis which actually assumes all your SNPs are in Hardy–Weinberg equilibrium.


print_text("Check we have at least 1K samples as this is relevant for the HWE threshold used", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/15_check_sex; \
    n_sample_after_sex_cleaning=$( \
        awk \
            'BEGIN{\" \"}END{print NR}'\
            ./03_second_check_sex/loop_maf_missing_2_pca_not_outliers_sex_full_clean.fam \
    ); \
    printf \"Number of samples after complete sex clean is: %s\n\" $n_sample_after_sex_cleaning; \
    if [[ $n_sample_after_sex_cleaning -lt 1000 ]]; then \
        exit 1; \
    fi; \
")
    #load the FAM file after fully clean sex issues and count the rows
    #if the number is lower than 1K, stop


print_text("run plink to filter by HWE", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./15_check_sex/03_second_check_sex/loop_maf_missing_2_pca_not_outliers_sex_full_clean \
        --hwe 1e-25 midp \
        --make-bed \
        --out ./16_hwe_second_round/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe \
")
    #--hwe filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. 
        #https://www.cog-genomics.org/plink/1.9/filter#hwe
        #We recommend setting a low threshold. Serious genotyping errors often yield extreme p-values like 1e-50 which are detected by any reasonable configuration of this test, while genuine SNP-trait associations can be expected to deviate slightly from Hardy-Weinberg equilibrium (so it's dangerous to choose a threshold that filters out too many variants).
        #Therefore, using a very low p-value would only remove SNPs extremely deviates from HWE and hence, being potential genotpying errors.
        #According to the plink tutorial: The 1e-25 threshold is extreme enough that we are unlikely to remove anything legitimate. (Unless the dataset is primarily composed of F1 hybrids from an artificial breeding program.). 
    #--hwe's 'midp' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium. The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data. We recommend its use.
        #We are using it, but it seems it does not change anything using it or not.
    #Because of the missing data issue, you should not apply a single p-value threshold across a batch of variants with highly variable missing call rates. A warning is now given whenever observation counts vary by more than 10%.
        #we should not have this problem because we have already clean the SNPs based on missingness
    #Note the notation of non-autosomal chromosomes si changed here from 23 to X, 24 to Y and so on... This is because we are using here plink2 to filter by HWE. Adding this step with plink2 changes the chromosome notation.


print_text("check the number of SNPs removed", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/16_hwe_second_round; \
    n_snps_remove=$( \
        awk \
            'BEGIN{\" \"}{ \
                if($0 ~/^--hwe midp:/){ \
                    print $3; \
                }; \
            }' \
        ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe.log \
    ); \
    printf \"The number of SNPs removed in the second HWE round is %s\" $n_snps_remove; \
    if [[ $n_snps_remove -gt 2 ]]; then \
        exit 1; \
    fi; \
")
    #Just load the log file and, in the row with the results of HWE filtering, extract the third fields (assuming space as delimiter) which is the number of SNPs removed.
    #print the number and if it is higher than 2, stop


print_text("update the chromosome names. in the previous step we used plink2 and the output files have X,Y,XY and MT as chromosome names, instead of number, which is what this script expects", header=3)
print_text("Create the update_chr.txt file", header=4)
with open('./data/genetic_data/quality_control/16_hwe_second_round/update_chr.txt', 'w') as f:
    f.write("X 23\n")
    f.write("Y 24\n")
    f.write("XY 25\n")
    f.write("MT 26\n")
    #by default, --update-chr expects the old chromosome name in the first column and the new one in the second one.

print_text("Run the PLINK command to update chromosome names and make a new bed file", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/16_hwe_second_round/; \
    plink \
        --bfile ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe \
        --update-chr ./update_chr.txt\
        --make-bed \
        --out ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr \
")
    #--update-chr, --update-cm, --update-map, and --update-name update variant chromosomes, centimorgan positions, base-pair positions, and IDs, respectively. By default, the new value is read from column 2 and the (old) variant ID from column 1, but you can adjust these positions with the second and third parameters. The optional fourth 'skip' parameter is either a nonnegative integer, in which case it indicates the number of lines to skip at the top of the file, or a single nonnumeric character, which causes each line with that leading character to be skipped. (Note that, if you want to specify '#' as the skip character, you need to surround it with single- or double-quotes in some Unix shells.)

# endregion






####################################
# region REMOVE NON-AUTOSOMAL ######
####################################
print_text("remove non-autosomal", header=2)
print_text("open folder", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./17_remove_non_autosomal; \
")

#Checking the multiple protocols I am following, I have detected that they usually remove all sex chromosomes after sex problems have been assesed. This is despite the potential increase in power if sex chromosomes are added. This is also the case for the Predict-HIIT study. It seems that the use of sex chromosomes for associations (and Polygenic Risk Scores; PRSs) requires different analyses starting from the beginning in the QC, for example, applying filters in a diferent way, and this continues in the "association stage". I honestly did not know about that and I have found out probably too late. At this stage, I do think it is worth it to go back multiple steps just to include the X chromosome (remember that Y and Mitochrondrial SNPs are not supported in TOPMed). The bulk of the data is in the autosomals and this is, for now, the norm in the GWAS studies. So I would just remove non-autosomal variants as a last step before imputation and move foward with that. Let me know if this makes sense for you!
    #https://github.com/RitchieLab/GWAS-QC-Internal?tab=readme-ov-file#step-8----exclude-data
    #https://pmc.ncbi.nlm.nih.gov/articles/pmid/32709988/
    #https://onlinelibrary.wiley.com/doi/10.1002/gepi.21782

#IMPORTANT: If you use X you would need to pass the X chromosome to the TOPMed server so they do the "Chromosome X pipeline" (see Ritchie paper) and then do additonal analyses (see above)


print_text("select autosomals", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./16_hwe_second_round/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr \
        --autosome \
        --make-bed \
        --out ./17_remove_non_autosomal/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals \
")
        #--autosome excludes all unplaced and non-autosomal variants, while --autosome-xy does not exclude the pseudo-autosomal region of X
            #https://www.cog-genomics.org/plink/1.9/filter

print_text("check we only have autosomals", header=3)
print_text("load BIM file", header=4)
loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_bim = pd.read_csv( \
    "./data/genetic_data/quality_control/17_remove_non_autosomal/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals.bim", \
    sep="\t", \
    header=None \
)
print_text("check we go from 1 to 22", header=4)
import numpy as np
test_autosomals = np.array_equal( \
    loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_bim.iloc[:,0].unique(), \
    np.array([i for i in range(1,23)]) \
)
if(test_autosomals == False):
    raise ValueError("ERROR! FALSE! WE HAVE NON-AUTOSOMAL SNPS")
else:
    print("OK")

#ANSWERS:
    #David: Makes sense to me and also should mean that any predictions are relevant for both sexes (I assume).
        #Not sure about that, but the point is that we are going to use just X to avoid problems.
        #besides controlling for sex, we may do the PRS calculation also separated by sex.

# endregion






######################################
# region LAST MAF-MISSING CHECK ######
######################################
print_text("start last maf-missing check before imputation", header=2)
print_text("open folder", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./18_last_maf_missing_check; \
")

#We repeat in two steps the MAF-missing snps filters and then sample filter to check that the previous removal of samples (due to sample relatedness, population structure and sex problems) did not change allele frequencies in a way that after applying MAF + missing filters again, we lose more samples. We removed samples, this could make some SNPs to have now MAF lower the threshold, you remove SNPs and then this in turn can make some samples to go below the missinginess threshold.

print_text("Remove SNPs that have low MAFs or high missingness", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./17_remove_non_autosomal/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals \
        --geno 0.01 \
        --maf 0.05 \
        --make-bed \
        --out ./18_last_maf_missing_check/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean \
")
    #we apply the same threshold than in the previous times that this filter was applied
    #It is ok to apply these filters together. You remove SNPs fue to one criteria and then SNPs due to the second criteria. The removal of one SNP due to low MAF is not going to influence the missingness of other SNP. The removal of SNPs can influence the missingness of the samples, because of that we apply the sample missingness filter at the end.
    #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
    #--maf filters out all variants with minor allele frequency below the provided threshold (default 0.01), while --max-maf imposes an upper MAF bound. Similarly, --mac and --max-mac impose lower and upper minor allele count bounds, respectively.
        #https://www.cog-genomics.org/plink/1.9/filter

print_text("check we did not lose so many SNPs", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/18_last_maf_missing_check; \
        awk \
        'BEGIN{\" \"}{ \
            if($0 ~ /variants removed due to missing genotype data \\(--geno\\)/) { \
                n_snps_missing_removed=$1; \
            }; \
            if($0 ~ /variants removed due to minor allele threshold\\(s\\)/) { \
                n_snps_maf_removed_maf=$1; \
            }; \
        }END{ \
            if(n_snps_missing_removed > 100 || n_snps_maf_removed_maf > 10000){ \
                exit 1; \
            } \
        }' \
        ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean.log \
")
    #extract the number of SNPs removed due to missing and maf from the log file. You have to use \\ to scape the parenthesis. You look for rows that END with the strings we are interested to extract there the first field, i.e., the number of SNPs removed.


print_text("Remove samples with high missingness", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/18_last_maf_missing_check/; \
    plink \
        --bfile ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean \
        --mind 0.01 \
        --make-bed \
        --out ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean \
")
    #we apply the same threshold than in the previous times that this filter was applied
    #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
        #https://www.cog-genomics.org/plink/1.9/filter

print_text("check we did not lose ANY sample at all", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/18_last_maf_missing_check; \
        awk \
        'BEGIN{\" \"}{ \
            if($0 ~ /people removed due to missing genotype data \\(--mind\\)/) { \
                n_samples_removed=$1; \
            }; \
        }END{ \
            if(n_samples_removed != 0){ \
                exit 1; \
            } \
        }' \
        ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean.log \
")
    #extract the number of samples removed due to missing from the log file. You have to use \\ to scape the parenthesis. You look for rows that END with the strings we are interested to extract there the first field, i.e., the number of samples removed.


print_text("check we do not have SNPs without ID", header=3)
#I have taken the code from Ritchie´s github
    #https://github.com/RitchieLab/GWAS-QC-Internal?tab=readme-ov-file#step-5----remove-snp-variants-that-do-not-have-snp-ids
print_text("create new fileset", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/18_last_maf_missing_check; \
    echo . > noSNP; \
    plink \
        --bfile loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean \
        --exclude noSNP \
        --make-bed \
        --out check_snp_ids; \
")
    #Creating the noSNP File: The command echo . > noSNP creates a file named noSNP with a single dot (.) in it. This file is not truly empty, but it effectively serves the purpose of an empty list for PLINK.
    #Then use that list to exclude SNPs. As the list is empty, no SNP is excluded. I understand that Ritchie is doing this because this will remove any SNP that does not have ID. If a SNP does not have ID, you cannot filter it by ID.

print_text("check we did not lose any variants because of this", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/18_last_maf_missing_check; \
        awk \
        'BEGIN{\" \"}{ \
            if($0 ~ /variants loaded from .bim file./) { \
                variants_loaded=$1; \
            }; \
            if($0 ~ /variants remaining./) { \
                variants_remaining=$2; \
            }; \
        }END{ \
            n_snps_removed=variants_loaded-variants_remaining; \
            if(n_snps_removed != 0){ \
                exit 1; \
            } \
        }' \
        ./check_snp_ids.log \
")
    #in the log file
    #extract the number of variants loaded and the number of variants remaining (rows ending with the relevant strings...). Calculate the difference and if it is different from 0, stop.

# endregion






################################
# region IMPUTATION PREP #######
################################
print_text("start last maf-missing check before imputation", header=2)
print_text("open folder", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./19_topmed_prep/00_first_step; \
")

#Several tools exist specifically for genotype imputation such as the Michigan and Trans-Omics for Precision Medicine (TOPMed) Imputation Servers where one uploads the phased or unphased GWAS genotypes in order to receive the imputed genomes in return. Each imputation server varies in terms of speed and accuracy. One of the most important considerations in imputation is the composition of the reference panel. Ritchie tutorial selected the TOPMed Imputation Reference panel (version r2) because it is one of the most diverse reference panels available and contains information from 97,256 deeply sequenced human genomes containing 308,107085 genetic variants distributed across the 22 autosomes and the X chromosome.
    #we have here the best correlation between datasets.

#We are following steps they follow in Ritchie´s Github to prepare files for imputation
    #https://github.com/RitchieLab/GWAS-QC-Internal?tab=readme-ov-file#step-11---sort-and-zip-files-to-create-vcf-files-for-imputation


print_text("Compare our list of SNPs with that of TOPMed", header=3)
print_text("download the tool and the data", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step; \
    rm HRC-1000G-check-bim-v4.3.0.zip; \
    rm HRC-1000G-check-bim.pl; \
    wget https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip; \
    unzip HRC-1000G-check-bim-v4.3.0.zip HRC-1000G-check-bim.pl; \
")
    #perl script to compare TOPMed with our data. This is created by the Wayner Tools group: https://www.well.ox.ac.uk/~wrayner/tools/
    #Make sure you've downloaded the following file: HRC-1000G-check-bim-v4.3.0.zip
    #Unzip and add the HRC-1000G-check-bim.pl script to your rawData/ directory
    #Changes V4.2.13 to V4.3.0
        #Updated the program to rework the reference panel loading so as to reduce runtime memory usage and time
        #It is highly recommended to use this version for TOPMed, previous versions will also work with TOPMed but will require >300GB RAM
        #This version still works with all reference panels listed below again with the reduction in memory usage
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step; \
    rm CreateTOPMed.zip; \
    rm CreateTOPMed.pl; \
    wget https://www.chg.ox.ac.uk/~wrayner/tools/CreateTOPMed.zip; \
    unzip ./CreateTOPMed.zip CreateTOPMed.pl \
")
    #tool to prepare the TOPMed for the comparison
run_bash("\
    cd ./data/genetic_data/; \
    cp \
        ./hrc_panel_data/bravo-dbsnp-all.vcf.gz \
        ./quality_control/19_topmed_prep/00_first_step/; \
")
    #The TOPMed reference panel is not available for direct download from this site, it needs to be created from the VCF of dbSNP submitted sites (currently ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz). This can be downloaded from the Bravo Website 
        #https://legacy.bravo.sph.umich.edu/freeze5/hg38/download
    #the file downlodad is not named as "ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz" because we used curl instead direct downloading by clicking from the website. See origin folder containing this file for further details.

print_text("confirm we have the exact checksum")
#calculate first the check sum and save as a file
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step; \
    rm checksum_bravo-dbsnp-all.vcf.txt; \
    checksum=$( \
        md5sum bravo-dbsnp-all.vcf.gz; \
    ); \
    echo $checksum > checksum_bravo-dbsnp-all.vcf.txt; \
")
#then do the check
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step; \
    awk \
        'BEGIN{\" \"}{ \
            if($0 ~ /bravo-dbsnp-all.vcf.gz/){ \
                checksum=$1; \
            } \
        }END{ \
            if(checksum != \"773e9e97759a4a5b4555c5d7e1e14313\"){ \
                exit 1; \
            } \
        }' \
        checksum_bravo-dbsnp-all.vcf.txt \
")
    #in the row with the name of the file, extract the first field with the checksum
    #MD5 checksum: 773e9e97759a4a5b4555c5d7e1e14313
        #https://legacy.bravo.sph.umich.edu/freeze5/hg38/download

print_text("prepare the TOPMed data", header=4)
#run the pearl script to do the prep
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step; \
    ./CreateTOPMed.pl -i bravo-dbsnp-all.vcf.gz \
")
    #once downloaded the VCF can be converted to an HRC formatted reference legend using the code here: CreateTOPMed.zip.
        #Usage: ./CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz 
    #If the shabang "./" doesn't work, you may need to run the command with "perl" instead. If you get an error in the above step, try this variation. Both were run successfully on a local computer and server using perl/5.30.0. Depending on your setup, this may take a few hours to run:  perl CreateTOPMed.pl -i bravo-dbsnp-all.vcf.gz	
    #By default this will create a file filtered for variants flagged as PASS only, if you wish to use all variants the -a flag overrides this. To override the default output file naming use -o filename.
    #the output should be: PASS.Variantsbravo-dbsnp-all.tab.gz
#unzip
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step; \
    gunzip --keep PASS.Variantsbravo-dbsnp-all.tab.gz \
")

print_text("obtain freq file from our data", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    input_file_set_name=loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean; \
    cp ./18_last_maf_missing_check/${input_file_set_name}.bim ./19_topmed_prep/00_first_step/; \
    cp ./18_last_maf_missing_check/${input_file_set_name}.fam ./19_topmed_prep/00_first_step/; \
    cp ./18_last_maf_missing_check/${input_file_set_name}.bed ./19_topmed_prep/00_first_step/; \
    cd ./19_topmed_prep/00_first_step/; \
    plink \
        --bfile ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean \
        --freq \
        --out ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean_freq; \
    head ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean_freq.frq \
")

print_text("do the actual comparison using HRC-1000G-check-bim.pl", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step/; \
    perl ./HRC-1000G-check-bim.pl \
        -b ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean.bim \
        -f ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean_freq.frq \
        -r ./PASS.Variantsbravo-dbsnp-all.tab \
        -h \
")
    #Usage with the TOPMed reference panel:
        #HRC-1000G-check-bim.pl is a program developed by Wayner Tools group to check a BIM file (from plink) against the HRC, 1000G or CAAPA reference SNP list in advance of imputation.
            #https://www.chg.ox.ac.uk/~wrayner/tools/ 
        #Requires the unzipped tab delimited HRC reference (currently v1.1 HRC.r1-1.GRCh37.wgs.mac5.sites.tab) from the Haplotype Reference Consortium Website here: http://www.haplotype-reference-consortium.org/site
        #Usage: 
            #-b <bim file> 
            #-f <Frequency file> 
            #-r <Reference panel> 
            #-h
            #-t <difference>, -n
                #two options for allele frequency thresholds (-t <difference>, -n)
                #-t 0.3 sets the allele difference threshold to 0.3, the default if not set is 0.2. Use this to change the allele frequency difference used to exclude SNPs in the final file, range 0 - 1, the larger the difference the fewer SNPs that will be excluded. I guess this removes SNPs that differ in allele frequency between input data and the reference panel.
                #-n flag to specify that you do not wish to exclude any SNPs based on allele frequency difference, if -n is used -t has no effect.
                #we are using the default as Ritchie.
            #-c 
                #Added new flag -c to specify checking individual chromosome(s) rather than assuming genome wide
                #Uisng default as Ritchie
            #-a
                #Added -a flag to disable automatic removal of palindromic SNPs with MAF > 0.4
    #Summary
        #Checks:
            #Strand, alleles, position, Ref/Alt assignments and frequency differences. In addition to the reference file v4 and above require the plink .bim and (from the plink --freq command) .frq files.
            #Implemented a check to ensure the same number of variants are present in the .bim and .frq files
        #Produces:
            #A set of plink commands to update or remove SNPs (see below and changes to V4.2.2) based on the checks as well as a file (FreqPlot) of cohort allele frequency vs reference panel allele frequency.
        #Updates:
            #Strand, position, ref/alt assignment
        #Removes:
            #A/T & G/C SNPs if MAF > 0.4
            #SNPs with differing alleles
            #SNPs with > 0.2 allele frequency difference (can be removed/changed in V4.2.2) between input data and reference panel
            #SNPs not in reference panel 
    #This tool is recommended by the TOPMed docs!
        #https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/

print_text("check we have used the correct version of the tool, panel, plink executable... using the LOG file", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step/; \
    fileset_name=loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean; \
    awk \
        -v fileset_name=${fileset_name} \
        'BEGIN{\" \"; count=0}{ \
            if($0 ~ /^Version/){ \
                if($2 != \"4.3\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Reference Panel:/){ \
                if($3 != \"HRC\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Bim filename:/){ \
                if($3 != \"./\" fileset_name \".bim\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Reference filename:/){ \
                if($3 != \"./PASS.Variantsbravo-dbsnp-all.tab\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Allele frequencies filename:/){ \
                if($4 != \"./\" fileset_name \"_freq.frq\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Plink executable to use:/){ \
                if($5 != \"plink\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Chromosome flag set:/){ \
                if($4 != \"No\"){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Allele frequency threshold:/){ \
                if($4 != 0.2){ \
                    exit 1; \
                } \
            } \
            if($0 ~ fileset_name){ \
                count++; \
            } \
        }END{ \
            if(count != 10){ \
                exit 1; \
            } else { \
                print \"OK\"; \
            } \
        }' \
        LOG-${fileset_name}-HRC.txt \
")
    #check in several rows we have the correct reference panel, thresholds, etc....
    #also count how many times we have the fileset name, that should be 10

print_text("check we do not have so many SNPs changed, lost... using the LOG file", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step/; \
    fileset_name=loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean; \
    awk \
        -v fileset_name=${fileset_name} \
        'BEGIN{\" \"}{ \
            if($0 ~ /^ Position different from HRC/){ \
                if($5 != 0){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^No Match to HRC/){ \
                if($5 > 7500){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^SNPs to change ref alt/){ \
                if($6 > 180000){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Strand to change/){ \
                if($4 > 31000){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Total removed for allele Frequency diff > 0.2/){ \
                if($9 > 2200){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Palindromic SNPs with Freq > 0.4/){ \
                if($7 > 400){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^Non Matching alleles/){ \
                if($4 > 4700){ \
                    exit 1; \
                } \
            } \
            if($0 ~ /^ID and allele mismatching/){ \
                if($5 > 4700){ \
                    exit 1; \
                } \
            } \
        }END{print \"OK\"}' \
        LOG-${fileset_name}-HRC.txt \
")
    #check we do not more than expected SNPs to be removed, changed, etc...
    #if everything is ok, do not stop and print OK at the end


print_text("run the bash script generated by TOPMED tool to solve problems found during the checks", header=3)
print_text("make dir", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./19_topmed_prep/01_second_step; \
")
#The TOPMed tool generates a plink script with the changes to be made, so I have taken that script and manually included here understanding all the steps followed.
#My understanding is that the script is always the same because there are two steps (chromosome and position update) that are not necessary in our case but they are still aplied, only with empty input files for snps to be modified.

print_text("move to the next step folder the required files generated during the checks", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    fileset_name=loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean; \
    mv ./00_first_step/${fileset_name}.bim ./01_second_step/; \
    mv ./00_first_step/${fileset_name}.bed ./01_second_step/; \
    mv ./00_first_step/${fileset_name}.fam ./01_second_step/; \
    mv ./00_first_step/Exclude-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/Chromosome-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/Position-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/Strand-Flip-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/Force-Allele1-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/FreqPlot-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/ID-${fileset_name}-HRC.txt ./01_second_step/; \
    mv ./00_first_step/Run-plink.sh ./01_second_step/; \
")

print_text("remove files not required", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/00_first_step/; \
    rm ./PASS.Variantsbravo-dbsnp-all.tab; \
    rm ./HRC-1000G-check-bim.pl; \
    rm ./CreateTOPMed.pl \
")

print_text("exclude SNPs required to be removed by TOPMed", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step/; \
    plink \
        --bfile ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean \
        --exclude ./Exclude-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
        --make-bed \
        --out ./TEMP1 \
")
    #Use as input the fileset with the name of the bim file used as input for the topmed script. Exclude those SNPs considered as problematic by the checks of the topmed tool. Generate a temporal file (TEMP1)

print_text("update chromosome if needed", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step/; \
    plink \
        --bfile ./TEMP1 \
        --update-map ./Chromosome-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
        --update-chr \
        --make-bed \
        --out ./TEMP2\
")
    #Use temporal file of the previous step as input, update the chromosomes that need to be updated.
    #--update-chr, --update-cm, --update-map, and --update-name update variant chromosomes, centimorgan positions, base-pair positions, and IDs, respectively.
    #You can combine --update-chr, --update-cm, and/or --update-map in the same run. So I understand that the TOPMed tools would create a file with chromosomes and positions, but not sure because in the next step they use update-map just with a position file
    #In our case it seems to be empty anyways, so no problem.
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step; \
    awk \
        'END{ \
            if(NR!=0){ \
                exit 1; \
            } \
        }' \
        Chromosome-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
")

print_text("update positions if needed", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step/; \
    plink \
        --bfile ./TEMP2 \
        --update-map ./Position-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
        --make-bed \
        --out ./TEMP3\
")
    #Use temporal file of the previous step as input, update the positions that need to be updated.
    #--update-chr, --update-cm, --update-map, and --update-name update variant chromosomes, centimorgan positions, base-pair positions, and IDs, respectively.
    #In our case it seems to be empty anyways.
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step; \
    awk \
        'END{ \
            if(NR!=0){ \
                exit 1; \
            } \
        }' \
        Position-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
")

print_text("solve flips", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step/; \
    plink \
        --bfile ./TEMP3 \
        --flip ./Strand-Flip-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
        --make-bed \
        --out ./TEMP4\
")
    #Flip alleles that the TOPMed tools pointed need to be flipped when comapring with the panel. Given a file containing a list of SNPs with A/C/G/T alleles, --flip swaps A↔T and C↔G. A warning will be given if any alleles are not named A, C, G, or T. To save the results instead of only applying the swap to the current run, combine this with --make-bed/--make-just-bim. If --make-bed is the only other operation in the run, you can also use --flip-subset, which only flips alleles for samples in the given ID list ('FID' family IDs in the first column, and 'IID' within-family IDs in the second column), and fails if any SNPs are not A/T or C/G.

print_text("force A2 alleles", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/01_second_step/; \
    plink \
        --bfile ./TEMP4 \
        --a2-allele ./Force-Allele1-loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-HRC.txt \
        --make-bed \
        --out ./loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated \
")
    #Force A1/A2 alleles using the list generated by the tool to match HRC: With --a1-allele, all alleles in the provided file are set to A1; --a2-allele does the reverse, i.e., all the alleles in the provided file are set to A2. If the original .bim file only has a single allele code and the --a1-allele/--a2-allele file names a second allele, a concurrent --make-bed will save both allele codes. If there are already two allele codes loaded and --a1-allele/--a2-allele names a third, a warning with the variant ID will be printed (you will usually want to resolve this with --exclude or --flip).

print_text("remove the TEMP files and create new folders for the next step", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    rm ./01_second_step/TEMP*; \
    mkdir -p ./02_third_step/00_filesets; \
    mkdir -p ./02_third_step/01_vcf_files; \
")

print_text("create files for each chromosome", header=4)
print("get the unique chromosomes")
import numpy as np
unique_chromosomes = pd.read_csv(
    "./data/genetic_data/quality_control/19_topmed_prep/01_second_step/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated.bim", \
    header=None, \
    sep="\t"
).iloc[:,0].sort_values().unique()
print(unique_chromosomes)

print("obtain filesets and VCF files for each chromosome")
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    for i in {"+str(np.min(unique_chromosomes))+".."+str(np.max(unique_chromosomes))+"}; do \
        plink \
            --bfile ./01_second_step/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated \
            --real-ref-alleles \
            --make-bed \
            --chr ${i} \
            --out ./02_third_step/00_filesets/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated-chr${i}; \
        plink \
            --bfile ./01_second_step/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated \
            --real-ref-alleles \
            --recode vcf \
            --output-chr chrM \
            --chr ${i} \
            --out ./02_third_step/01_vcf_files/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated-chr${i}; \
    done; \
")
    #From here, two steps are repeated across chromosomes
        #step 1
            #Generate a plink fileset of the interest chromsome but preserving the order of the alleles is preserved as we have changed in the last step to match that of HRC: If a binary fileset was originally loaded, --keep-allele-order forces the original A1/A2 allele encoding to be preserved; otherwise, the major allele is set to A2. --real-ref-alleles has that effect as well, and also removes 'PR' from the INFO values emitted by "--recode vcf{,-fid,-iid}".
        #step 2
            #Generate a VCF of the interest chromsome but preserving the order of the alleles is preserved as we have changed in the last step to match that of HRC: Note that --real-ref-alleles also removes 'PR' from the INFO values emitted by "--recode vcf{,-fid,-iid}", so I guess they want to remove this field. Also specify the notation for chromosome names that is accepted by TOPMed using --output-chr chrM. PLINK 1.9 and 2.0 support seven chromosome coding schemes in output files. You can select between them by providing the desired human mitochondrial code. chrM: Autosomes are 'chr' followed by a numeric code, X/Y/XY/M are preceded by 'chr', PAR1/PAR2 as usual. This is required for TOPMed.
                #We added this step following Ritchie github 
                    #https://github.com/RitchieLab/GWAS-QC?tab=readme-ov-file#step-10----calculate-frequency-files-and-compare-to-topmed-panel
            #Note that the 'vcf', 'vcf-fid', and 'vcf-iid' modifiers in --recode result in production of a VCFv4.2 file. 'vcf-fid' and 'vcf-iid' cause family IDs and within-family IDs respectively to be used for the sample IDs in the last header row. In other words, 'vcf-iid' takes the individual ID to be used as ID in the VCF file, while 'vcf-fid' does the same but for the family ID. In contrast, 'vcf' merges both IDs and puts an underscore between them (in this case, a warning will be given if an ID already contains an underscore). We prefer the latter so we keep track of the batch and sample IDs.
    #Note about the warning: "Underscore(s) present in sample IDs."
            #When using --recode vcf, sample IDs are formed by merging the FID and IID and placing an underscore between them. When the FID or IID already contains an underscore, this may make it difficult to reconstruct them from the VCF file; you may want to replace underscores with a different character in PLINK files (Unix tr is handy here).
            #this is ok. I prefer to maintain the format "combat_ILGSA...." for backwards compatibility. We will just use awk to split the FAM and IDs using two delimiters ("_" and "-")

print("check that we have the correct number of files generated")
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    filesets_count=$(ls -1 ./02_third_step/00_filesets | wc -l); \
    vcfs_count=$(ls -1 ./02_third_step/01_vcf_files | wc -l); \
    expected_fileset="+str(np.max(unique_chromosomes)*4)+"; \
    expected_vcfs="+str(np.max(unique_chromosomes)*2)+"; \
    if [[ $filesets_count -ne $expected_fileset || $vcfs_count -ne $expected_vcfs ]]; then \
        exit 1; \
    else \
        echo \"OK\"; \
    fi; \
")
    #count the number of files in each folder
        #ls -1 /path/to/folder: Lists all files in the specified folder, one per line.
        #wc -l: Counts the number of lines, which corresponds to the number of files.
    #calculate the expected number of files, per each chromosome
        #4 files in fileset folder (bed, bim, fam and log files)
        #2 files in the VCF folder (vcf and log files)
    #if any of the counts is not as expected, stop execution


print_text("Solve flips and order SNPS", header=3)
print_text("make dir", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir -p ./19_topmed_prep/03_fourth_step; \
")

print_text("download the reference genome and the check sum", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/03_fourth_step; \
    rm hg38.p13.fa; rm md5sum.txt; \
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz; \
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/md5sum.txt; \
    gunzip --keep hg38.p13.fa.gz \
")
    #I am using the fasta file for hg38.p13 according to UCSC. AGRF said we have hg38.p13!
        #Thank you for your patience whilst I confirmed the patch information.  I can confirm the manifest file version is hg38.p13.
        #https://mail.google.com/mail/u/1/#inbox/FMfcgzQXKWgNjbBcBnqPzcPbxfJNVzsq
    #We are going to use the fasta file of hg38.p13 to solve the strand issues, but I also did it with the latest (hg38.p14) and got the same results.
    #I downloadad the fasta from "broad resources". It seems that the Ritchie tutorial used Broad resources (resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta) but the file is not present in github. The closest thing I found in the bundle of Broad was "Homo_sapiens_assembly38.fasta", but when you compared (cmp wise) with the original hg38 fasta from USCS, there are differences. Broad it is suppose to use the original version, not patches, but still no match, so we are not suing this.
        #https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
        #https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
    #Also according to copilot, we should not use the masked versions present in "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"
        #For using bcftools fixref, you should use the unmasked reference FASTA file. In this case, you should use hg38.fa from the UCSC server1. The masked versions (e.g., hg38.fa.masked.gz) are not suitable for this purpose as they contain modifications that can interfere with the reference checking process.
    #Indeed, the Ritchie tutorial does not use masked fasta

print("do checksum of the file")
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/03_fourth_step; \
    checksum=$( \
        md5sum hg38.p13.fa.gz; \
    ); \
    echo $checksum > checksum_fasta_hg38_p13.txt; \
")

print("compare with the original checksum")
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/03_fourth_step; \
    awk \
        'BEGIN{\" \"}{ \
            if(FILENAME==\"md5sum.txt\"){ \
                if($0 ~ /hg38.p13.fa.gz/){ \
                    original_checksum=$1; \
                } \
            } \
            if(FILENAME==\"checksum_fasta_hg38_p13.txt\"){ \
                if($0 ~ /hg38.p13.fa.gz/){ \
                    new_checksum=$1; \
                } \
            } \
        }END{ \
            if(new_checksum != original_checksum){ \
                exit 1; \
            } \
        }' \
        md5sum.txt checksum_fasta_hg38_p13.txt \
")
    #first process md5sum.txt, and get the checksum (first field) of hg38.p13.fa.gz (if row ends with that name), which is the reference
    #first process checksum_fasta_hg38_p13.txt, and get the checksum (first field) of hg38.p13.fa.gz (if row ends with that name), which is the new checksum
    #stop if the new checksums are not identical

print_text("run a loop to solve flips using the USCS fasta as reference", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    mkdir -p ./03_fourth_step/vcf_files_flipped; \
    mkdir -p ./03_fourth_step/fixref_stats; \
    input_fileset_name=loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated; \
    for i in {"+str(np.min(unique_chromosomes))+".."+str(np.max(unique_chromosomes))+"}; do \
        bcftools +fixref \
            ./02_third_step/01_vcf_files/${input_fileset_name}-chr${i}.vcf \
            -Ov -o ./03_fourth_step/vcf_files_flipped/${input_fileset_name}_flipped_chr${i}.vcf -- \
            --discard \
            --fasta-ref ./03_fourth_step/hg38.p13.fa \
            --mode flip &> ./03_fourth_step/fixref_stats/fixref_stats_chr${i}.txt; \
    done; \
")
    #I initially decided this approach after reading in detail the docs of bcftools +fixref (see below), but then I checked the github of Ritchie, and they apply exactly the same approach!
    #-O, --output-type b|u|z|v[0-9]
        #Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion. The compression level of the compressed formats (b and z) can be set by by appending a number between 0-9.
        #We are using -Ov to output uncompressed VCF so we can use it in the next step for sorting
        #https://samtools.github.io/bcftools/bcftools.html
    #-d, --discard: Discard sites which could not be resolved
    #-f, --fasta-ref FILE.fa: Reference sequence
    #-mode flip: swap or flip REF/ALT columns and GTs for non-ambiguous SNPs and ignore the rest
        #WARNING FROM BCFTOOLS:
            #Do not use the program blindly, make an effort to understand what strand convention your data uses! Make sure the reason for mismatching REF alleles is not a different reference build!! Also do NOT use bcftools norm --check-ref s for this purpose, as it will result in nonsense genotypes!!
            #strand
                #I know the strand reference of my data, I completely sure I used the foward strand notation, not top/bot or other. You can check that in "01b_illumina_report_to_plink.py", where I used pyspark to ensure I selected the correct that from all FinalReports. 
                #This should the strand used by the imputation server as it is the one used by ncbi, see "01b_illumina_report_to_plink.py" about details.
            #I am also sure about the reference build, our data has hg38 as confirmed by AGRF. So we can use hg38 in the imptuation server.
                #Thank you for reaching out with your detailed query regarding the reference genome for project "CAGRF20093767". I can confirm that the reference genome used to generate your data is hg38
                #It’s great to hear that you achieved a high overlap (~98%) using the 1000 Genomes Project hg38 high-coverage panel. While the strand issues you encountered with the initial imputation are not uncommon when reconciling different datasets, using bcftools to address allele switches is a sound approach, and the resulting reduction in median allele switches for hg38 aligns with the reference genome used for this project.
                #Regarding the ~4K SNPs with differing allele frequencies, this may reflect population-specific differences or residual inconsistencies in strand alignment between panels. 
                #Also the page of TOPMED says the reference build is hg38
                    #https://topmedimpute.readthedocs.io/en/latest/getting-started/#build
            #Therefore, we can be sure that any strand issue between our dataset and the reference panel in the imputation server is not caused by the reference build or the strand. In other words, we do not have allele flips because I have made a mistake and I am comparing two genetic datasets that are using different strand formats or different builds.
        #we use the flip model ("-m flip") to swap or flip REF/ALT columns and GTs for non-ambiguous SNPs to foward and ignore the rest.
            #According to the manual, this is the following: Assuming the reference build is correct, just flip to fwd, discarding the rest
                #This is ok, I selected foward notation in our data, so all our SNPs should be foward. SNPs that are not in foward, are errors and we should solve fliping to foward if possible.
                #Also, the imputation server is going to check for flips, if more than 10000 obvious strand flips are detected, TOPMED stops the imputation, so we are ok.
            #We avoid ambiguous sites (A/T, C/G). According to chatGTP:
                #Strand Ambiguity: For non-palindromic SNPs, the reverse complement will change the alleles (e.g., A/C becomes T/G). However, for A/T and C/G SNPs, the reverse complement is identical, making it difficult to ascertain if the strand needs correction.
        #---discard: To remove the cases that have been deemed problematic but have not been solved due to ambiguity.
            #THIS IS THE KEY: Doing this removes all the strange SNPs whose allele frequencies were negatively correlated between our dataset and the refenrece panel of the imputation server.

print_text("check all REF allele match and we have a reduced number of unsolved problems", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    for i in {"+str(np.min(unique_chromosomes))+".."+str(np.max(unique_chromosomes))+"}; do \
        awk \
            -v chr_name=${i}\
            'BEGIN{\"\t\"}{ \
                if($1==\"NS\" &&  $2==\"total\"){ \
                    total_n_sites=$3; \
                } \
                if($1==\"NS\" &&  $2==\"ref\" &&  $3==\"match\"){ \
                    total_ref_match=$4; \
                } \
                if($1==\"NS\" &&  $2==\"unresolved\"){ \
                    n_unresolved=$3; \
                } \
            }END{ \
                if(total_n_sites != total_ref_match || n_unresolved>500){ \
                    exit 1; \
                } else { \
                    print \"Chr\" chr_name \" is OK\";\
                } \
            }' \
            ./03_fourth_step/fixref_stats/fixref_stats_chr${i}.txt; \
        if [[ $? -ne 0 ]]; then \
            exit 1; \
        fi; \
    done; \
")
    #In the file with fixref stats, extract the total number of sites, the number of alleles with matching REF and the number of unsolved SNPs.
    #All alleles have to have their REF matching and the number of unsolved should be small.
    #if this is not met, stop. Given we are in a loop, this will stop the current chromosome, but the loop will continue. We have to check the exit stauts
    #if [ $? -ne 0 ]: checks the exit status of the awk command, if not zero, it means we had an error in at least one chromosome and we stop the execution.

print_text("sort by position and compress", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    mkdir -p ./04_fifth_step/; \
")
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    for i in {"+str(np.min(unique_chromosomes))+".."+str(np.max(unique_chromosomes))+"}; do \
        bcftools sort \
            ./03_fourth_step/vcf_files_flipped/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated_flipped_chr${i}.vcf \
            -Oz -o ./04_fifth_step/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated_flipped_sorted_chr${i}.vcf.gz; \
    done; \
")
    #Sort VCF/BCF file. 
        #-o, --output FILE: output file name
        #-O, --output-type b|u|z|v: b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
    #This is the approach recommended by TOPMed docs.
        #https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/
    #The Ritche tutorial used "vcf-sort" (from VCFtools) instead of bcftools sort. I am not sure why, because "bcftools sort" just sorts the VCF file using the genomic position. 
        #The bcftools sort command is used to sort the variants in a VCF or BCF file based on their chromosomal positions, and the basic and only syntax of the bcftools sort command is the following one.
            #https://www.biocomputix.com/post/bcftools-sort
            #https://github.com/RitchieLab/GWAS-QC-Internal?tab=readme-ov-file#step-11---sort-and-zip-files-to-create-vcf-files-for-imputation

print("check that we have the correct number of files generated")
run_bash(" \
    cd ./data/genetic_data/quality_control/19_topmed_prep/; \
    vcfs_count=$(ls -1 ./04_fifth_step | wc -l); \
    expected_vcfs="+str(np.max(unique_chromosomes))+"; \
    if [[ $vcfs_count -ne $expected_vcfs ]]; then \
        exit 1; \
    else \
        echo \"OK\"; \
    fi; \
")


print_text("check all requeriments of TOPMed are already met", header=3)
#CHECK requeriments for inputs to the TOPMed server
    #Create a separate vcf.gz file for each chromosome.
        #DONE
    #Variants must be sorted by genomic position.
        #DONE
    #GRCh37 or GRCh38 coordinates are required.
        #DONE
    #If your input data is GRCh37/hg19, please ensure chromosomes are encoded without prefix (e.g. 20). If your input data is GRCh38/hg38, please ensure chromosomes are encoded with prefix 'chr' (e.g. chr20).
        #DONE (chrM mode of plink)
    #VCF files need to be version 4.2 (or lower). This is specified in the VCF file header section.
        #DONE (output format of "vcf" in plink)
    #Must contain GT field in the FORMAT column. All other FORMAT fields will be ignored. (if you are seeing problems with very large uploads, it may help to remove other FORMAT fields)
        #DONE (just checked the VCF files)
    #Due to server resource requirements, there is a maximum of 25k samples per chromosome per job (and a minimum of 20 samples). Please see the FAQ for details.
        #DONE
    #https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/

# endregion






##############################
# region IN CASE NEEDED ######
##############################

#In case you need to liftover the data, you can use the script of Ritchie, but they say that you do not needed in general because TOPMed accepts both hg38 and hg19
    #https://github.com/RitchieLab/GWAS-QC?tab=readme-ov-file#step-7----liftover-the-data

#if you need to use X, you can combine X and PAR in one VCF file ("--chr 23,25")
        #Respect to XY: For phasing and imputation, chrX is divided into three independent chunks (PAR1, non-PAR, PAR2). These chunks are then automatically merged by the Michigan Imputation Server 2 and returned as a single complete chromosome X file. Therefore, we are generating a VCF file with the X and the PAR regions together, being all their SNPs named as "chrX".
            #https://genepi.github.io/michigan-imputationserver/pipeline/

# endregion
