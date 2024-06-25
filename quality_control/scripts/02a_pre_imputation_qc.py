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




#############################
###### merging batches ######
#############################
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

#Remember that prob intensity is used to split samples between the three genotypes so you plot prob intesnity in try to define three clusters of samples. Sometimes the clusters of hetero and one homozygous can be very close making it difficult to separate... T
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/

#I have NOT checked in each final report whether all genotypes per sampe have gene call > 0.15.
    #However, the DNAReport of both batches says that "Low GenCall Score Cutoff" is 0.15. 
    #I have also checked in these reports that the "Call_Rate" column is very similar to divide the column "#Calls" (number of genotypes with GS score > 0.15) by the total number of SNPs.
    #this call rate is what in the PDF of the first batch they say that is too low for 4 samples, being below the 99% call rate expected for illumina.
    #Indeed, I have checked that those samples with call rate < 0.99 in the DNA report of the first batch are indeed those with call rate < 0.99 in the PDF report.
    #Therefore, I understand that AGFR has applied the filter of GC score < 0.15 to set as no call a given genotype and now we can use this to calculate the call rate per sample and per snp in order to filter by missingness.
    



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




print_text("filter by sample missingness", header=2)
print_text("make sample missing report after previous filters", header=3)
run_bash("head ./data/genetic_data/quality_control/05_remove_missing_snps/plink.imiss")
    #imiss: A text file with a header line, and one line per sample with the following six fields:
        #FID    Family ID
        #IID Within-family ID
        #MISS_PHENO  Phenotype missing? (Y/N)
        #N_MISS  Number of missing genotype call(s), not including obligatory missings or het. haploids
        #N_GENO  Number of potentially valid call(s)
        #F_MISS  Missing call rate
            #https://www.cog-genomics.org/plink/1.9/formats#imiss
    #IMPORTANT: 
        #I understand that heterozigous cases for hayploid regions (i.e., non-PAR X-Y regions) are not considered in the missing count, we should take care of this later! The cool thing is that they do not affect when calculating missing samples, so we would not lose samples because of these problematic cases as they do not increase the number of missing per sample
        #Also, obligatory missing are are not counted, for example, Y snps are obligatory missing in females, we should not count these.


print_text("check that F_MISS is just the number of missing genotype divided by the number of potentially valid calls", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/05_remove_missing_snps/; \
    n_samples_freq_file=$( \
        awk \
            'BEGIN{FS=\" \"; OFS=\" \"} \
            {if(NR>1)(count++)}\
            END{print count}'\
            plink.imiss \
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
            plink.imiss \
    ); \
    if [[ $n_correct_freq_miss -eq $n_samples_freq_file ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #just check if freq miss is the ratio of missing calls respect to the potential number of valid calls per sample. We count the number of samples for which this is the case, N_MISS/N_GENO is equal to F_MISS and then check this is the total number of samples, i.e., all meet the condition.
        #see above in the SNP missing part for explanations about the script


print_text("load the report with awk and then save it controlling the delimiter. If I load it directly into pandas, I got problems separating the columns", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/05_remove_missing_snps/; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{if(NR>0){print $1,$2,$3,$4,$5,$6}}' \
        plink.imiss \
    > plink_awk_processed.imiss; \
    head plink_awk_processed.imiss")
        #specify the delimiter of the output (OFS="\t") and just print all fields (1 to 6) for all rows, then save as a new file


print_text("load in pandas", header=3)
sample_missing_report = pd.read_csv( \
    "./data/genetic_data/quality_control/05_remove_missing_snps/plink_awk_processed.imiss", \
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

#create a new folder for doing operations about missing in samples
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./06_remove_low_call_samples/; \
    ls -l")


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
    fname="./data/genetic_data/quality_control/06_remove_low_call_samples/merged_batches_plot_sample_retention_thresholds.png")
plt.close()


print_text("remove samples with low call rate", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    plink \
        --bfile ./05_remove_missing_snps/merged_batches_remove_missing_snps \
        --mind 0.01 \
        --make-bed \
        --out ./06_remove_low_call_samples/merged_batches_remove_low_call_samples")
            #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
                #https://www.cog-genomics.org/plink/1.9/filter
            #recommended threshold
                #the PRS tutorial recommends to retain samples with "sample missingness <0.02", i.e., samples with more than 2% missing or less than 98% call rate should be removed. 
                #The Ritchie tutorials says "A recommended threshold is 98% to 99% efficiency after first removing markers that have a low genotype call rate across samples."
                #The illumina report of the first batch says "Four samples are below the illumina expected 99% SNP call rate (values expected for typical projects, excluding tumour samples)".
                #we are going to use 99%, i.e., < 0.01 missingness. This is higher than recommended by PRS tutoria, but within the range recommended by Ritchie tutorial and the value for Illumina. In the second batch, this threshold does not increase the number of samples removed compared to 98%, while in the first one only makes 1 more sample to be removed, one with call rate=0.981. Remember that the FPD report says that this is below the expectation of Illunina. Therefore, I think we are ok removing these samples.


print_text("see the number of samples removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_samples_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./05_remove_missing_snps/merged_batches_remove_missing_snps.fam); \
    n_samples_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./06_remove_low_call_samples/merged_batches_remove_low_call_samples.fam); \
    lost_samples=$(($n_samples_before - $n_samples_after)); \
    lost_samples_percent=$( \
        awk \
            -v l=$lost_samples \
            -v b=$n_samples_before \
            'BEGIN{print (l/b)*100}'); \
    if [[ $lost_samples_percent < 2 ]]; then \
        printf 'The number of samples lost due to call rate below 0.99 is: %s' \"$lost_samples\"; \
    else \
        echo 'ERROR: FALSE! We have lost more than 2% of samples due low-call rate';\
    fi")
        #see SNP missing code to details about awk steps


print_text("create again the missing report", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/06_remove_low_call_samples/; \
    plink \
        --bfile merged_batches_remove_low_call_samples \
        --missing; \
    ls -lh")
        #--missing produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.
            #https://www.cog-genomics.org/plink/1.9/basic_stats#missing


print_text("check that no sample has a missing % above the selected threshold", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/06_remove_low_call_samples/; \
    total_number_samples=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>1){count++}}END{print count}' \
            plink.imiss);\
    n_samples_below_threshold=$( \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($6 <= 0.01)){count++}}END{print count}' \
            plink.imiss); \
    if [[ $n_samples_below_threshold -eq $total_number_samples ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #see SNP missing code to details about awk steps




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
        07_remove_high_LogRDev_samples; \
    plink \
        --bfile ./06_remove_low_call_samples/merged_batches_remove_low_call_samples \
        --remove ../cn_metrics/merged_batches_samples_filter_out_LogRDev.tsv \
        --make-bed \
        --out ./07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples; \
    ls -l ./07_remove_high_LogRDev_samples/")
            #create a new folder to save plink filesets after filtering
            #then use --remove to remove samples included in a file with two columns: family id and within-family id.
                #this file is in a different parent folder so we have to use "../"
                #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.
                    #https://www.cog-genomics.org/plink/1.9/filter


print_text("check that the corresponding samples are indeed not included in the new fam file", header=3)
fam_file_after_logRdev_filtering = pd.read_csv( \
    "./data/genetic_data/quality_control/07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples.fam", \
    sep=" ", \
    header=0, \
    low_memory=False)
print(sum(fam_file_after_logRdev_filtering.iloc[:, 1].isin(total_removed_samples_logRdev)) == 0)
    #no sample ID (second column in fam file) should be included in the list of IDs to be removed


print_text("see the number of samples removed", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/; \
    n_samples_before=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./06_remove_low_call_samples/merged_batches_remove_low_call_samples.fam); \
    n_samples_after=$( \
        awk \
            'BEGIN{FS=\" \"}{if(NR>0){count++}}END{print count}'\
            ./07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples.fam); \
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
print_text("create missing and freq reports", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/07_remove_high_LogRDev_samples/; \
    plink \
        --bfile merged_batches_remove_high_LogRDev_samples \
        --freq \
        --missing \
        --out merged_batches_remove_high_LogRDev_samples_reports; \
    ls -l")




print_text("create a folder to save plink data after new filters", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./08_loop_maf_missing; \
    ls -l")


print_text("check if we have SNPs below the MAF threshold after removing samples", header=3)
n_snps_below_maf_first_round_raw = run_bash(" \
    cd ./data/genetic_data/quality_control/07_remove_high_LogRDev_samples/; \
    awk \
        'BEGIN{FS=\" \"}{if((NR>1) && ($5 < 0.05)){count++}}END{print count}' \
        merged_batches_remove_high_LogRDev_samples_reports.frq", return_value=True).strip()
n_snps_below_maf_first_round = 0 if n_snps_below_maf_first_round_raw=="" else int(n_snps_below_maf_first_round_raw)
print(f"We have {n_snps_below_maf_first_round} SNPs below the MAF threshold after the first round of filters")


print_text("check if we have SNPs below the missing threshold after removing samples", header=3)
n_snps_below_missing_first_round_raw = run_bash(" \
    cd ./data/genetic_data/quality_control/07_remove_high_LogRDev_samples/; \
    awk \
        'BEGIN{FS=\" \"}{if((NR>1) && ($5 > 0.01)){count++}}END{print count}' \
        merged_batches_remove_high_LogRDev_samples_reports.lmiss", return_value=True).strip()
n_snps_below_missing_first_round = 0 if n_snps_below_missing_first_round_raw=="" else int(n_snps_below_missing_first_round_raw)
print(f"We have {n_snps_below_missing_first_round} SNPs above the missing threshold after the first round of filters")


print_text("repeat the filters if these snps exists", header=3)
if(n_snps_below_maf_first_round>0) | (n_snps_below_missing_first_round>0):
    print_text("apply MAF and missing SNP filters", header=4)
    run_bash(" \
      cd ./data/genetic_data/quality_control/; \
      plink \
          --bfile ./07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples \
          --maf 0.05 \
          --geno 0.01 \
          --make-bed \
          --out ./08_loop_maf_missing/loop_maf_missing_1; \
      ls -l ./08_loop_maf_missing/")  

    print_text("apply missing sample filters", header=4)
    run_bash(" \
        cd ./data/genetic_data/quality_control/; \
        plink \
            --bfile ./08_loop_maf_missing/loop_maf_missing_1 \
            --mind 0.01 \
            --make-bed \
            --out ./08_loop_maf_missing/loop_maf_missing_2; \
        ls -l ./08_loop_maf_missing/")

    print_text("create MAF and missing reports", header=4)
    run_bash(
        "cd ./data/genetic_data/quality_control/08_loop_maf_missing/; \
        plink \
            --bfile loop_maf_missing_2 \
            --freq \
            --missing \
            --out loop_maf_missing_2_reports; \
        ls -l")
    
    print_text("check reports again", header=4)
    n_snps_below_maf_second_round_raw = run_bash(" \
        cd ./data/genetic_data/quality_control/08_loop_maf_missing/; \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 < 0.05)){count++}}END{print count}' \
            loop_maf_missing_2_reports.frq", return_value=True).strip()
    n_snps_below_maf_second_round = 0 if n_snps_below_maf_second_round_raw=="" else int(n_snps_below_maf_second_round_raw)
    print(f"We have {n_snps_below_maf_second_round} SNPs below the MAF threshold after the second round of filters")    
    n_snps_below_missing_second_round_raw = run_bash(" \
        cd ./data/genetic_data/quality_control/08_loop_maf_missing/; \
        awk \
            'BEGIN{FS=\" \"}{if((NR>1) && ($5 > 0.01)){count++}}END{print count}' \
            loop_maf_missing_2_reports.lmiss", return_value=True).strip()
    n_snps_below_missing_second_round = 0 if n_snps_below_missing_second_round_raw=="" else int(n_snps_below_missing_second_round_raw)
    print(f"We have {n_snps_below_missing_second_round} SNPs above the missing threshold after the second round of filters")
    n_samples_below_missing_second_round_raw = run_bash(" \
        cd ./data/genetic_data/quality_control/08_loop_maf_missing/; \
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
        ./07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples.bim \
        ./08_loop_maf_missing/loop_maf_missing_2.bim; \
    cp \
        ./07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples.fam \
        ./08_loop_maf_missing/loop_maf_missing_2.fam; \
    cp \
        ./07_remove_high_LogRDev_samples/merged_batches_remove_high_LogRDev_samples.bed \
        ./08_loop_maf_missing/loop_maf_missing_2.bed; \
    ls -l ./08_loop_maf_missing/")
else:
    print("Step not required.")




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



#in the last BIM file after MAF filters, select ID of those SNPs with code 25 that are outside the PAR regions previously indicated.
print_text("see SNPs that are considered to be pseudoautosomals (chr=25) but in reality are not in pseudo-autosomal regions of XY", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/08_loop_maf_missing; \
    awk \
        'BEGIN{FS=\"\t\"}{ \
            if($1 == 25){ \
                if($4 < 10001 || $4 > 2781479 && $4 < 155701383 || $4 > 156030895){ \
                    print $2 \
                } \
            } \
        }' \
        ./loop_maf_missing_2.bim > snps_par_problem.txt"
)
#load bim file after MAF filtering to awk using tabs as delimiter (checked this is the delimiter in the file)
#select those SNPs in pseudo autosomal regions (code 25) that
#are located: 1) before the start of the first PAR region;
#2) After the first PAR region and before the second one;
#3) After the second PAR region
#print these cases and count them
#at the END, print the count



#check we only have 36 problematic cases
print_text("check we ONLY have 36 of these problematic cases", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/08_loop_maf_missing; \
    n_par_problem=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./snps_par_problem.txt); \
    if [[ $n_par_problem -eq 36 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi"
)



#make a file with ID of NON problematic SNPs
print_text("create a list WITHOUT SNPs that are problematic for PAR: We discard SNPs that are considered to be pseudoautosomals but in reality are not in pseudo-autosomal regions of XY", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/08_loop_maf_missing; \
    awk \
        'BEGIN{FS=\"\t\"}{ \
            if($1 == 25){ \
                if($4 >= 10001 && $4 <= 2781479 || $4 >= 155701383 && $4 <= 156030895){ \
                    print $2 \
                } \
            } else { \
                print $2 \
            } \
        }' \
        ./loop_maf_missing_2.bim > snps_par_no_problem.txt"
)
#if the chromosome is 25 (PAR)
    #if the SNPs is within the PAR limits accoridng to NCBI print the ID (second column)
    #if not, then it is a SNP considered pseudo-autosomal but being outside of the PAR
        #region, so out.
#if the SNP is not considered PAR, then we can print the ID




#retain only these SNPs without the problem
#We cannot be sure if there is something wrong with these SNPs
#They are in PAr regions or not? should they be treated as PAR o like
#sex chromosomes? so we are going to check the number is not very high
#and then remove all of them. We should be ok with the remaining PAR SNPs
#as they are considered separately from autosomals and sex chromosomes
print_text("retain only these SNPs without the problem", header=3)
run_bash(
    "cd ./data/genetic_data/quality_control/08_loop_maf_missing; \
    plink \
        --bfile ./loop_maf_missing_2 \
        --extract ./snps_par_no_problem.txt \
        --make-bed \
        --out ./loop_maf_missing_3; \
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
    cd ./data/genetic_data/quality_control/08_loop_maf_missing; \
    n_snps_before_par=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./loop_maf_missing_2.bim); \
    n_snps_after_par=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./loop_maf_missing_3.bim); \
    n_snps_par_problem=$( \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}' \
            ./snps_par_problem.txt); \
    sum_snps=$(($n_snps_after_par + $n_snps_par_problem)); \
    if [[ $n_snps_before_par -eq $sum_snps ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi"
)




print_text("sample relatedness", header=2)
#the plink tutorial first remove related samples before filtering by MAF and calculate the PCA.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
#sample relatedness is also explained before than population substructure in Ritchie's tutorial. Indeed, they also talk about LD pruning, which is an step needed for selecting SNPs for the PCA.
    #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
#In the VO2 max paper, they also seem to remove the related individuals using pi-hat before doing the PCA
    #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
#we can filter by sample relatdness, select a subset of SNPs based on LD-pruning (for IBD), then filter by relatedness and use that subset also for PCA and heterozigosity
#Summary:
    #1. Sample relatedness: it is usually done as one of the first steps in many of the tutorials checked
    #2. Population stratification: It is strongly recommended by Plink´s author to check subgroups and then perform checks like sex-imbalances inside each group. Besides this, we will use the PCAs as covariates in the analyses.
    #3. Sex imbalances. The differences in allele frequencies between ancestry groups can influence the check for sex imbalances, so we have to do it after population stratification analyses. This should be ok, because the problem with sex imbalances could be that the phenotype of one sample is ineed the phenotype or other sample, i.e., they are swapped, but it should influence if we just do analyses with only genotypes like the PCA.

#Rationale:
    #Many population-genomic statistics (such as the allele frequencies) and analyses are distorted when there are lots of very close relatives in the dataset; you are generally trying to make inferences about the population as a whole, rather than a few families that you oversampled. For example, PLINK 2 includes an implementation of the KING-robust [5] pairwise relatedness estimator, which can be used to prune all related pairs. This does not mean that both samples in each related pair are thrown out. Instead, –king-cutoff tries to keep as much data as possible, and as a consequence it usually keeps one sample out of each pair (see below).
    #These related samples, if treated as independent samples in the downstream analyses, having many related samples in the dataset would result in increased type I and type II errors. The options are tu use of mixed-regression models while considering in place of simple linear or logistic regression.
    #Retaining a maximal set of unrelated individuals is computationally expensive (NP-hard), but an efficient greedy approximation is available in PLINK 2.0 using the --king-cutoff flag (as opposed to just removing one of each pair of related individuals)
        #see below about KING
    #Cryptic relatedness can interfere with the association analysis. If you have a family‐based sample (e.g., parent‐offspring), you do not need to remove related pairs but the statistical analysis should take family relatedness into account. However, for a population based sample we suggest to use a pi‐hat threshold of 0.2, which in line with the literature (Anderson et al., 2010; Guo et al., 2014).
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/

print_text("create folder for this step", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./09_remove_related_samples/; \
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
        --bfile ./08_loop_maf_missing/loop_maf_missing_3 \
        --indep-pairwise 500kb 1 0.2 \
        --out ./09_remove_related_samples/ldpruned_snplist; \
    ls -l ./09_remove_related_samples/")
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
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
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
        --bfile ./08_loop_maf_missing/loop_maf_missing_3 \
        --extract ./09_remove_related_samples/ldpruned_snplist.prune.in \
        --make-bed \
        --out ./09_remove_related_samples/loop_maf_missing_3_ld_pruned;\
    ls -l ./09_remove_related_samples/")
        #from the current fileset, select only those SNPs included in .prune.in
        #we use extract for that
            #--extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.
                #https://www.cog-genomics.org/plink/1.9/filter
        #make a new fileset

print_text("remove sex chromsomes", header=4)
#According to Marees et al. (2018), we should check for sample relatedness not only using independent SNPs (pruning), but also limiting the analysis to autosomal chromosomes only.
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
#In Ritchie's tutorial (figure 5), they show an histogram for the distribution of pairwise pi-hat. They calculated IBD after "removing sex inconsistent individuals, 95% SNP call rate, 90% sample call rate, 10% MAF, and pruning to 67,000 AUTOSOMAL variants.".
    #Therefore they do not use SNPs in sex chromosomes!
    #Respect sex inconsistencies, we can have minorities within the sample, and as plink help says (--check-sex), imbalanced ancestries can give problems so in that case you have to do the check of sex within each ancestry group. Therefore we need to check the PCA before.
#In Ritchie's GitHub, they remove sex chromosomes before PCA: "Exclude any SNPs that do not liftOver and non-somatic chromosomes (X, Y)"
#They are talking specifically about autosomals, so we are going to consider ONLY autosomals, no X, Y, MT nor PAR regions. The case I was doubting more was PAR regions because these behave like autosomal chromosomes, but they are sexual chromosomes. They are just 500 in total, so we are going to remove them.

#Therefore, I think we can use this set of pruned autosomal SNPs for kinship and PCA.
run_bash(" \
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
    plink \
        --bfile ./loop_maf_missing_3_ld_pruned \
        --autosome \
        --make-bed \
        --out ./loop_maf_missing_3_ld_pruned_autosomals;\
    ls -l")
        #--autosome excludes all unplaced and non-autosomal variants, while --autosome-xy does not exclude the pseudo-autosomal region of X
            #https://www.cog-genomics.org/plink/1.9/filter
print("Do we have at least 70K autosomal SNPs after LD pruning like in Ritchie's tutorial?")
#In the Ritchie's tutorial, they ended up with 67,000 autosomal variants in linkage equilibrium in order to calculate IBD and pi_hat.
run_bash(" \
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
    auto_ld_snps_in=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./loop_maf_missing_3_ld_pruned_autosomals.bim); \
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
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
    n_non_auto_snps=$( \
        awk \
            'BEGIN{FS=\" \"}{if($1==0 || $1==23 || $1==24 || $1==25 || $1==26){count++}}END{print count}'\
            ./loop_maf_missing_3_ld_pruned_autosomals.bim); \
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
        --bfile ./09_remove_related_samples/loop_maf_missing_3_ld_pruned_autosomals \
        --genome \
        --out ./09_remove_related_samples/ibd_report;\
    ls -l ./09_remove_related_samples/")
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
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
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
                #Individuals sharing one allele IBD at every locus (Z1~1) are parent-child pairs. On average, siblings share zero, one, and two alleles IBD at 25%, 50%, and 25% of the genome, respectively. Therefore, Z0=0.25, Z1=0.5 and Z2=0.25.
            #Z2  P(IBD=2)
                #(IBD). Individuals sharing two alleles IBD at nearly every locus (Z2~1) are monozygotic twins, or the pair is a single sample processed twice.
            #PI_HAT  Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
                #This counts the whole probability of having 2 shared alleles IBD and then half of the probability of having 1 shared allele IBD.
                #I guess this is the whole probability of have a shared allele at all.

    ##POR AQUII

    #particular cases
        #Identical twins, and duplicates, are 100% identical by descent (Pihat 1.0)
            #P(IBD=2)=1, i.e., Z2=1
        #First-degree relatives are 50% IBD (Pihat 0.5)
            #P(IBD=1)=1, i.e., Z1=1. Therefore, 0.5*1=0.5
        #Second-degree relatives are 25% IBD (Pihat 0.25)
        #Third-degree relatives are 12.5% equal IBD (Pihat 0.125).
            #https://www.biostars.org/p/58663/
            #https://www.biostars.org/p/75335/

            
            #PHE Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)
            #DST IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
            #PPC IBS binomial test
            #RATIO   HETHET : IBS0 SNP ratio (expected value 2)

print_text("convert delimiter of ibd report", header=4)
run_bash(" \
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \
        ibd_report.genome > ibd_report_awk_processed.genome.tsv; \
    gzip \
        --force \
        ./ibd_report_awk_processed.genome.tsv; \
    ls -l")

print_text("load ibd report as pandas DF", header=4)
ibd_report = pd.read_csv( \
    "./data/genetic_data/quality_control/09_remove_related_samples/ibd_report_awk_processed.genome.tsv.gz", \
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
    fname="./data/genetic_data/quality_control/09_remove_related_samples/pairs_relatdness_before_filtering.png")
plt.tight_layout()
plt.close()
    #Our results:
        #We have many samples with Z0~1, i.e., unrelated individuals.
        #We do not have cases with Z1 close to 1, so to parent-child pairs, but we have cases with Z1=0.5, thus potential siblings. 
        #We have one pair with Z2=1, thus we have twins or the same sample duplicated.

print_text("plot histogram of pi_hat", header=4)
ibd_report["PI_HAT"].hist(bins=100)
    #PI_HAT  Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
        #This counts the whole probability of having 2 shared alleles IBD and then half of the probability of having 1 shared allele IBD.
        #It sums Z2 plus half Z1.W
plt.yscale('log') 
    #to gain visibility, like in Ritche tutorial
        #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
plt.savefig( \
    fname="./data/genetic_data/quality_control/09_remove_related_samples/pi_hat_hist.png")
plt.close()
    #we have several cases above pi_hat of 0.2









##if we do it without the MAF filtering, we lose 3 samples due sibling-sibling/father-sibling, but after the MAF filter we lose 20!! We have to check with other approaches! because KING says that is better not remove SNPs (at least for LD prunning) so maybe we are making it difficult for the approach.




#As long as you're just concerned with finding/pruning close relations (1st-2nd degree, maybe third degree), you can think of KING kinship as half of PI_HAT.
    #https://groups.google.com/g/plink2-users/c/z2HRffl-6k8/m/ndOZIjpMBQAJ



print_text("run KING-robust with plink2", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    plink2 \
        --bfile ./09_remove_related_samples/loop_maf_missing_3_ld_pruned_autosomals \
        --king-cutoff 0.177 \
        --make-bed \
        --out ./09_remove_related_samples/remove_related_samples; \
    ls -l ./09_remove_related_samples")
    #PLINK 2 includes an implementation of the KING-robust [5] pairwise relatedness estimator, which can be used to prune all related pairs. This does not mean that both samples in each related pair are thrown out. Instead, –king-cutoff tries to keep as much data as possible, and as a consequence it usually keeps one sample out of each pair.
        #apply the KING-robust method
            #https://www.cog-genomics.org/plink/2.0/distance#make_king



###IMPORTANT: LD pruning is not recommended in KING


run_bash(" \
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
    awk \
        'BEGIN{FS=\" \"}{if(NR>1){count++}}END{print count}' \
        ./remove_related_samples.king.cutoff.out.id")

run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    cat ./09_remove_related_samples/remove_related_samples.king.cutoff.out.id")

#CHECK THIS REMOVED SAMPLES ARE THOSE WITH PI_HAT>0.2
    #0.2 WAS USED IN THE VO2 PAPER, BUT NOT SURE, CHECK FIGURE 5 RITCHIE where <0.2 are 3rd degree related or unrelated.


#maybe remove the sample with more missing from each related pair


        #check general usage of plink2 and king

            #it seems this approach tries to balance between getting a subset of unrelated samples and speed. This does not just remove one sample of each related pair

            #Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING.
                #https://www.kingrelatedness.com/manual.shtml

            #LOOK to this ancestry problem:
                #The exception is that KING-robust underestimates kinship when the parents are from very different populations. You may want to have some special handling of this case; --pca can help detect it.

            #The same samples are removed with threshold 0.354 and 4. Maybe these are duplicated samples?
                #only 3 more samples are removed with 0.177, i.e., when removing first-degree relations (parent–child and sibling–sibling)



#talk with Bishop once we have results of sample relatdness
#in the meantime work on teaching


















#LD equlibbrium es solo para un subset que se usa solo para PCA!! no para el resto de analisis!

#LD seem to be done before PCA, but also HWE? plink tutorial does but not ritche tutorial doing after imputation, also Augusot paper does it after imptuation

#the plink tutorials says to do it before and after PCA: 
    #The “keep-fewhet” modifier causes this filter to be applied in a one-sided manner (so the fewer-hets-than-expected variants that one would expect from population stratification would not be filtered out by this command)
    #if a snAfter you have a good idea of population structure in your dataset, you may want to follow up with a round of two-sided –hwe filtering, since large (see Note 5) violations of Hardy–Weinberg equilibrium in the fewer-hets-than-expected direction within a subpopulation are also likely to be variant calling errors; with multiple subpopulations, the –write-snplist and –extract flags can help you keep just the SNPs which pass all subpopulation HWE filters.p can violate HWE becuase pop structure, should be control (remove ancestry outliers) before filtering for HWE.
#the VO2 max paper does HWE before PCA



#Different ethnicities can be included in the same study, as long as the population substructure is considered to avoid false positive results
    #general tutorial gwas


#several tutorials say that the MAF (and LD) filtering should be done before the PCA, then filter by MAF again after imputation
    #not sure if LD prunning should be done in general or only for the snps of the PCA
    #see EIGENSOFT

    #ritche says to filter snps used in PCA by MAF and LD..., use eigensoft... see the vairance of the PCA axes to select those included in the models...
    
    #Once you have LD-pruned and MAF-filtered your dataset, PLINK 2’s –pca command has a good shot of revealing large-scale population structure
        #EIGENSOFT [7, 8] has some additional built-in principal component analysis options, including automated iterated outlier removal, and a top-eigenvalue-based test for significant population structure
        #If there are obvious clusters in the first few plots, I recommend jumping ahead to Chapter 4 (on ADMIXTURE) and using it to label major subpopulations before proceeding
        #plink tutorial

    #SNPs with Minor Allele Frequency (MAF)>0.05 were then used to perform principal component analysis (PCA) for ethnicity identification using SHELLFISH [45]. Ethnic and ancestry outliers (more than 6 standard deviations from the mean on either of the two first principal components (PCs)) were excluded (n=10). 
        #paper VO2 max




#after filtering for maf, you should not have any SNP with 0 in the allele column for minor. in that moment, we would have remove snps with very low minor frequencies.
#look at the bim file the unique cases
#bim_no_dups.loc[:, 4].unique()
#bim_no_dups.loc[:, 5].unique()



















##CHANGE NAME OF THE FOLDER PCA, THIS IS NO LONGER THE 5TH

#do case-control study for batch effects after PCA? or after all pre-QC steps?

#In Ritchie's GitHub, they remove sex chromosomes before PCA: "Exclude any SNPs that do not liftOver and non-somatic chromosomes (X, Y)"






##clustering
#https://www.cog-genomics.org/plink/1.9/strat#clustering


##clustering is IMCOMPLETE, there are clusterions otpions that change a lot and create different clusters
#https://www.cog-genomics.org/plink2/strat
#https://zzz.bwh.harvard.edu/plink/strat.shtml#options


##one of the options is PPC
#Pairwise Population Concordance (PPC) in PLINK is a method used to determine whether two individuals belong to the same random-mating population. It's a significance test used during the clustering process in population stratification analysis¹.
#The PPC test is applied as a restriction during the clustering process. The clustering process is based on pairwise identity-by-state (IBS) distance and uses complete linkage agglomerative clustering¹. The PPC test ensures that clusters are not merged if they contain individuals that significantly differ, implying they belong to different populations¹.
#The PPC test is invoked with a p-value threshold. For example, to only merge clusters that do not contain individuals differing at a p-value of 0.0001, you would use the command --ppc 0.0001¹.
#This approach helps ensure that the subsequent association tests are valid, as they would implicitly match every case with its nearest control, as long as the case and control do not show evidence of belonging to different populations¹.
#Please note that this explanation is a high-level overview of the PPC test in PLINK. For a more detailed understanding, you may want to refer to the official PLINK documentation or relevant genetic analysis textbooks.
#Source: Conversation with Copilot, 6/25/2024
#(1) PLINK: Whole genome data analysis toolset - Harvard University. https://zzz.bwh.harvard.edu/plink/strat.shtml.
#(2) Second-generation PLINK: rising to the challenge of larger and richer .... https://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0047-8.
#(3) How to study runs of homozygosity using PLINK? A guide for analyzing .... https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6463-x.




run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./XX_cluster; \
    ls -l")


run_bash("\
    cd ./data/genetic_data/quality_control/; \
    plink \
        --cluster \
        --family \
        --bfile ./09_remove_related_samples/loop_maf_missing_3_ld_pruned_autosomals \
        --out ./XX_cluster/cluster_after_logR_filter")
    #--cluster uses IBS values calculated via "--distance ibs"/--ibs-matrix/--genome to perform complete linkage clustering. The clustering process can be customized in a variety of ways.
        #IBD calculations are not LD-aware! It is usually a good idea to perform some form of LD-based pruning before invoking them
            #https://www.cog-genomics.org/plink/1.9/ibd
    #In the context of population stratification, the `--cluster` function in PLINK is used to identify and account for population structure in genetic data (1). This is important in genetic association studies because population structure can lead to false positive results if not properly accounted for.
    #PLINK uses a method called **complete-linkage hierarchical clustering** to assess population stratification³. This method starts by considering every individual as a separate cluster of size 1, then repeatedly merges the two closest clusters³. 
    #By default, the distance between two clusters is defined as the maximum pairwise distance between a member of the first cluster and a member of the second cluster. The 'group-avg' modifier causes average pairwise distance to be used instead, i.e., we took the distance between earch pair of samples (sample form cluster 1 and sample from cluster 2) and calculate the average across all pairs.
    #The distance between clusters is calculated based on pairwise identity-by-state (IBS) distance¹. Identity-by-state (IBS) refers to the genetic similarity between two individuals. Two individuals are said to be identical by state at a particular genetic locus if they have the same alleles at that locus, regardless of whether those alleles were inherited from a common ancestor¹.
    #By clustering individuals based on IBS distance, PLINK can identify groups of individuals that are genetically similar to each other. These groups can then be used to account for population structure in genetic association analyses¹.
    #The 'missing' modifier causes clustering to be based on identity-by-missingness instead of identity-by-state. In other words, samples are clustered based on their missing data so we create groups of samples with the same pattern of missing genotypes, instead of obtaining groups of similar genotypes.
    #"cc" and "within/family" flags
        #By default, in PLINK, each individual starts in their own cluster. However, if the `--within` or `--family` flag is present, that cluster assignment is used as the starting point instead. For example, all samples of the batch 1 start together while samples of batch 2 start together.
        #the `cc` modifier in the `--cluster` function of PLINK is used to control how clusters are merged during the clustering process.
        #When the `cc` modifier is used with `--cluster`, it prevents two all-case or two all-control clusters from being merged². This means that each cluster will contain at least one case and one control. This can be particularly useful in case-control studies, where you want to ensure that each cluster contains a mix of cases and controls. For consistency with --mcc, missing-phenotype samples are treated as controls (this is a change from PLINK 1.07).
        #This approach can help to account for population stratification in genetic association studies, by ensuring that the genetic similarities identified by the clustering process are not simply due to similarities in case or control status. In our case, this would mean that the cluster are not just caused by the batch.
    #(1) PLINK: Whole genome data analysis toolset - Harvard University. https://zzz.bwh.harvard.edu/plink/strat.shtml.
    #(2) PLINK: a tool set for whole-genome association and population-based .... https://europepmc.org/article/MED/17701901.
    #(3) Using PLINK for Genome-Wide Association Studies (GWAS) and ... - Springer. https://link.springer.com/protocol/10.1007/978-1-62703-447-0_8.
    #(4) https://www.cog-genomics.org/plink/1.9/strat#clustering

#results clustering after exploring the different flags and modifiers
    #Default mode puts all samples in the same cluster
    #Adding the "--family" (i.e., batch) flag makes no difference
    #Also adding the "cc" modifier force the separation between batches, but this seems to be an artifact. If we just use "cc" without "--family", every sample is put in a diferent cluster, i.e., we get 1446 clusters. Also note that when looking the PCA (see above), there is no clear grouping based on the batch, and the clustering without "cc" suggests that. Therefore, I do not think we have reasons to think we have batch effects.
    #group-avg does not change anything.
    #clustering by missing genotypes does not create new groups.
    #In summary, it seems we do not have clusters of samples with similar genotypes/missing genotypes.




run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./XX_mds; \
    ls -l")


run_bash("\
    cd ./data/genetic_data/quality_control/; \
    plink \
        --cluster \
        --ppc 0.0005 \
        --mds-plot 4 \
        --bfile ./09_remove_related_samples/loop_maf_missing_3_ld_pruned_autosomals \
        --out ./XX_mds/cluster_after_logR_filter")

    #ppc=0.0005 makes very difficult to decided NOT to merge two clusters, so if the clusters are still not merged, it means that there is a lot of difference. For example, ppc=0.0005 makes no difference between han chinese and japenese, as we increase to 0.05, now the threshold is stringent enough to make difficult to combine japanese and chinese. 
        #you have to study the ppc test
            #https://zzz.bwh.harvard.edu/plink/strat.shtml#options
            #https://zzz.bwh.harvard.edu/plink/strat.shtml#outlier
    #in our case, ppc=0.0005 already is able to detect differences and create clusters, so it means we have samples with more difference than chinese and japanese.


run_bash(" \
    cd ./data/genetic_data/quality_control/XX_mds/; \
    awk \
        'BEGIN{ \
            FS=\" \"; \
            OFS=\"\t\"};\
        { \
            if(NR>1){ \
                print $1, $2, $3, $4, $5, $6, $7 \
            } \
        }' \
        ./cluster_after_logR_filter.mds > cluster_after_logR_filter_tab.mds"
)




mds_results = pd.read_csv("./data/genetic_data/quality_control/XX_mds/cluster_after_logR_filter_tab.mds", sep="\t", header=None)
mds_results

import plotly.express as px
fig = px.scatter(\
    data_frame=mds_results, \
    x=3, \
    y=4, \
    color=2, #the family/batch (0) or the cluster (2)
    hover_data=[\
        0,1,2,3])
        #you can use columns of DF to add axis data, but also modify color, size, and show data per sample in a desplegable box
        #https://plotly.com/python/line-and-scatter/
        #https://plotly.com/python-api-reference/generated/plotly.express.scatter.html
#fig.show()
fig.write_html("./data/genetic_data/quality_control/XX_mds/mds_after_logR_filter_plot.html")
    #https://plotly.com/python/interactive-html-export/

#it is like the PCA



run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./XX_neighbour; \
    ls -l")


run_bash("\
    cd ./data/genetic_data/quality_control/; \
    plink \
        --neighbour 1 5 \
        --cluster \
        --ppc 0.0005 \
        --bfile ./09_remove_related_samples/loop_maf_missing_3_ld_pruned_autosomals \
        --out ./XX_neighbour/cluster_after_logR_filter")

    #https://www.cog-genomics.org/plink/1.9/strat#neighbour
    #https://www.cog-genomics.org/plink/1.9/formats#nearest
    #https://zzz.bwh.harvard.edu/plink/strat.shtml#outlier

    #check --ppc in cluster above! explore this in the future?
        #for ppc we need LD pruned? if this is using IBD, then it shoudl be LD prunned
        #we get hudnreds of cluster with --ppc of 0.0005 or 0.05




print_text("start PCA to detect individuals that are outliers", header=2)
print_text("prepare folder for PCA", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    mkdir \
        --parents \
        ./XX_pca; \
    ls -l")



print_text("run PCA to detect clusters of samples", header=3)
run_bash("\
    cd ./data/genetic_data/quality_control/; \
    plink \
        --pca \
            'tabs' \
            'header' \
        --bfile ./09_remove_related_samples/loop_maf_missing_3_ld_pruned_autosomals \
        --out ./XX_pca/pca_after_logR_filter")
        #dimensionality reduction
            #PLINK 1.9 provides two dimension reduction routines: --pca, for principal components analysis (PCA) based on the variance-standardized relationship matrix, and --mds-plot, for multidimensional scaling (MDS) based on raw Hamming distances. 
            #Top principal components are generally used as covariates in association analysis regressions to help correct for population stratification, while MDS coordinates help with visualizing genetic distances.
            #By default, --pca extracts the top 20 principal components of the variance-standardized relationship matrix; you can change the number by passing a numeric parameter. Eigenvectors are written to plink.eigenvec, and top eigenvalues are written to plink.eigenval. 
            #The 'header' modifier adds a header line to the .eigenvec file(s), and the 'tabs' modifier makes the .eigenvec file(s) tab- instead of space-delimited.
                #https://www.cog-genomics.org/plink/1.9/strat#pca



print_text("load the eigenvec file generated", header=3)
pca_after_logR_filter=pd.read_csv(\
    "./data/genetic_data/quality_control/XX_pca/pca_after_logR_filter.eigenvec", \
    sep="\t", \
    header=0, \
    low_memory=False)
pca_after_logR_filter
        #Produced by --pca. Accompanied by an .eigenval file, which contains one eigenvalue per line.
        #The .eigenvec file 
            #is, by default, a space-delimited text file with no header line and 2+V columns per sample, where V is the number of requested principal components. 
            #The --pca 'header' modifier causes a header line to be written, and the 'tabs' modifier makes this file tab-delimited.
            #The first two columns are the sample's FID/IID, and the rest are principal component weights in the same order as the .eigenval values (if the header line is present, these columns are titled 'PC1', 'PC2', ...).
        #https://www.cog-genomics.org/plink/1.9/formats#eigenvec



print_text("make interactive plot of the two first PCAs, so you can zoom in", header=3)
import plotly.express as px
fig = px.scatter(\
    data_frame=pca_after_logR_filter, \
    x="PC1", \
    y="PC2", \
    color="FID",
    hover_data=[\
        "IID"])
        #you can use columns of DF to add axis data, but also modify color, size, and show data per sample in a desplegable box
        #https://plotly.com/python/line-and-scatter/
        #https://plotly.com/python-api-reference/generated/plotly.express.scatter.html
#fig.show()
fig.write_html("./data/genetic_data/quality_control/XX_pca/pca_after_logR_filter_plot.html")
    #https://plotly.com/python/interactive-html-export/


#but how much is this? we should check in the conext of 1000 KGP....





##USE OTHER TECHINES TO CHECK FOR OUTLIERS AND POP STRATIFICATION?
    #look tutorials
    #https://www.cog-genomics.org/plink/1.9/strat


#check differences in pheno between batches?







#####SEX INCONSISTENCES ######

#sex inconsistences on ld_pruned
    #In the tutorial of plink Chirstopher checks sex on the prunned dataset
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec18
    #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex


#Plink says that, due to the use of allele frequencies we may need to check sex within ancestry groups, and he did that in the tutorial, so we are doing this after pop structure analysis
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec19


run_bash(" \
    cd ./data/genetic_data/quality_control/09_remove_related_samples; \
    plink \
        --bfile loop_maf_missing_3_ld_pruned \
        --check-sex; \
    ls -l")


run_bash(" \
    cd ./data/genetic_data/quality_control/09_remove_related_samples/; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1, $2, $3, $4, $5, $6}'\
        plink.sexcheck > plink.sexcheck_awk_processed.tsv; \
    head plink.sexcheck_awk_processed.tsv")


sex_check_report = pd.read_csv( \
    "./data/genetic_data/quality_control/09_remove_related_samples/plink.sexcheck_awk_processed.tsv", \
    sep="\t", 
    header=0, 
    low_memory=False)
print(sex_check_report)




import matplotlib.pyplot as plt

plt.hist(sex_check_report["F"], 50, density=True, alpha=0.4, label='Observed iHS')
plt.savefig( \
    fname="./data/genetic_data/quality_control/09_remove_related_samples/check_sex_f_distribution.png")
plt.close()
    
    #density?

    #from plink doc
        ##0.66 is, of course, still much larger than 0.2, and in most contexts it still justifies a female call. We suggest running --check-sex once without parameters, eyeballing the distribution of F estimates (there should be a clear gap between a very tight male clump at the right side of the distribution and the females everywhere else), and then rerunning with parameters corresponding to the empirical gap.
    #this is exactly what we have, a tight peak close to 10 and then, below 0.2 we have a wider peak and a small.
    #therefore, the default F thresholds work for us
        #By default, F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male calls. If you pass numeric parameter(s) to --check-sex (without 'y-only'), the first two control these thresholds.


sex_check_report.loc[(sex_check_report["STATUS"] == "PROBLEM") & (sex_check_report["PEDSEX"] !=0), :]

sex_check_report.loc[(sex_check_report["STATUS"] == "PROBLEM") & (sex_check_report["PEDSEX"] ==0), :]
    #42 cases here with unknown sex, but in the excel the number of samples without sex is 41
    #there are also 42 cases with sex zero in the fam file of the second batch before merging
    #the origin of the fam file is the map file!
    #The problem is 2399LDJA/2397LDJA, this mislabeled sample. 2397LDJA is present in the first batch but not in the pheno data so it is sex is "0". In contrast, 2399LDJA is present in the pheno data, having sex as M, but not in the first batch. Therefore, when we count the number of NaN sex in pheno data is only 41 (2399LDJA has sex data), while the genetic data has 42 (2397LDJA has no sex data).
    #we should change the name of this sample. in the genetic data?

#OJO MINORITIES
    #the results seems to make sense with most females below 0.2, CHECK THAT!
    #we need to check ancestry? see ritchie


#https://groups.google.com/g/plink2-users/c/4bpdLMdH2KA
    #If heterozygous haploid calls still remain, the most likely cause is nonmissing female genotype calls on the Y chromosome; others have reported that this is fairly common.  A quick way to check the number of these is to just load the Y chromosome with e.g. "plink --bfile semi_clean_fileset --chr 24 --freq".  If all the heterozygous haploid errors are on the Y chromosome, you can safely clobber them with --make-bed + --set-hh-missing.  (If some are on the X, --set-hh-missing *might* still be okay, but I'd need to know more about the data source and the --check-sex report to be sure.)










###when finished this script, you should check missing threholds, hetero and PCA plots












print_text("heterozigosity", header=2)


#Ritchie's tutorial says "while less than expected heterozygosity (mean – 3 SD) suggests possible inbreeding and greater than expected heterozygosity (mean + 3 SD) suggests possible sample contamination. How- ever, these thresholds should take into account the expected heterozygosity rates in the ancestry group under study, as some diverse populations may exhibit different rates of heterozygosity than other populations."

#the heterozygosity of each individual is influenced by the genotypes, and we have already cleaned low-quality calls BUT the problem can be with the threshold we use to remove samples with low/high heterozigosity, because this is influenced by ancestry.


##important
#check if you need pruned data and if we can use sex chromosomes or not
#we have two sets of LD pruned data, one with sex and another without sex chromosomes


#I guess we can then use PCA to remove outlier samples because of ancestry issues, also by heterogizogisty those factors that can be influenced by ancestry
#do plot hetero - missingness
#use the information of both approaches to remove samples because different ancestry (outlier PCA), contamination (high hetero) or inbreeding (low hetero)



#####see tutorials






##################################
###### check sex mismatches ######
##################################
print_text("check sex mismatches", header=1)




    #check sex uses allele frequencies, so you can have problems with different ancestries, yo have to prepare the data... so maybe is bettter to see the PCA for outliers and batch effects, filtering and then go to see sex once we have cleaner data
        #"Due to the use of allele frequencies, if your dataset has a highly imbalanced ancestry distribution, you may need to process the rare-ancestry samples separately."
        #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex

    #check-sex has to be used on LD-pruned data, so we can use the previous data
        #Since this function is based on the same F coefficient as --het/--ibc, it requires reasonable MAF estimates (so it's essential to use --read-freq if there are very few samples in your immediate fileset), and it's best used on marker sets in approximate linkage equilibrium.
    #This isn't implemented yet in plink2 since there would be little practical difference from the plink 1.9 implementation.  Use "--make-bed --chr X,Y" to export a .bed fileset with only chrX and chrY, and run plink 1.9 --check-sex on that.

    #BUT OF COURSE WE NEED DATA OF SEX CHROMOSOMES, select the prunnn ed dataset with sex chromosomes




##THIS IS IMPORTANT, CHECK THIS

#create temporary folder to save
import tempfile
temp_dir = tempfile.TemporaryDirectory()
#print(temp_dir.name)


#read only the zip and get list of files in it
import numpy as np
import pandas as pd
import zipfile

list_sample_maps = []
#batch="ILGSA24-17873"
for batch in ["ILGSA24-17303", "ILGSA24-17873"]:

    if batch == "ILGSA24-17303":
        zip_name = "ILGSA24-17303"
    else:
        zip_name = "CAGRF20093767"

    zipdata = zipfile.ZipFile("data/genetic_data/illumina_batches/" + zip_name + ".zip")
    zipinfos = zipdata.infolist()
    zipinfos_subset = zipinfos[np.where([zipinfo.filename == zip_name + "/Sample_Map.txt" for zipinfo in zipinfos])[0][0]]
    #extract the file in the temp dict
    zipdata.extract(zipinfos_subset, temp_dir.name)
    #

    selected_sample_map = pd.read_csv(temp_dir.name+ "/" + zip_name + "/Sample_Map.txt",
            delimiter="\t",
            header=0,
            low_memory=False)

    selected_sample_map["batch"] = batch

    list_sample_maps.append(selected_sample_map)

[all(sample_map.columns == list_sample_maps[0].columns) for sample_map in list_sample_maps]

sample_map_illumina = pd.concat(list_sample_maps)


sample_map_illumina.shape[0] == 1248+216



########################################################################
#DO NOT FORGET TO ASK DAVID QUESIONS ABOUT PHENO IN TODO.MD
########################################################################


#load pheno data, this include reported sex and VO2 max data. I have checked that the data is the same directly reading from excel than converting to csv
pheno_data = pd.read_excel(
    "data/pheno_data/combact gene DNA GWAS 23062022.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data)


pheno_data.loc[pheno_data["Gender"] == "M", "Gender"] = "Male"
pheno_data.loc[pheno_data["Gender"] == "F", "Gender"] = "Female"


merge_sample_pheno = sample_map_illumina[["batch", "ID", "Gender"]].merge(
    pheno_data[["AGRF code", "Gender"]],
    left_on="ID", #use the column with IDs in sample_map
    right_on="AGRF code", #use the column with IDs in pheno_data
    suffixes=["_illumina", "_pheno_excel"], #set the suffix for repeated columns
    how="outer")

#cases with NA for IDs
print(merge_sample_pheno.loc[merge_sample_pheno["ID"].isna(), :])
print(merge_sample_pheno.loc[merge_sample_pheno["AGRF code"].isna(), :])

#2397LDJA is in ILGSA24-17303 but not in the pheno data, while 2399LDJA is in the pheno data but not in any illumina report.

subset_mismatch = merge_sample_pheno.loc[
    (merge_sample_pheno["Gender_illumina"] != merge_sample_pheno["Gender_pheno_excel"]) &
    (~merge_sample_pheno["ID"].isna()) & 
    (~merge_sample_pheno["AGRF code"].isna()) & 
    (~merge_sample_pheno["Gender_pheno_excel"].isna()), :]

subset_mismatch.loc[subset_mismatch["Gender_illumina"] == "Unknown", :].shape
subset_mismatch.loc[subset_mismatch["Gender_illumina"] != "Unknown", :].shape

subset_mismatch.loc[subset_mismatch["Gender_illumina"] != "Unknown", :]

subset_mismatch.to_csv("sample_sex_mimatch.txt",
    sep="\t",
    header=True,
    index=False)

#sex (--check-sex) and hetero (--het) should be checked after PCA because accorindg to plink info, it can be problems if we have a sample with most of samples from one ancestry and then a few from another ancestry
    #https://www.cog-genomics.org/plink/1.9/basic_stats

    #R log is present in our data, so we could check prob intensity in X for a full detail sex determination, think about it.

    #with sex check done, tell David about mismatches

    #Warning: 30589 het. haploid genotypes present (see
    #./04_inspect_snp_dup/01_remove_dup/ILGSA24-17873_merged_data_no_snp_dup.hh );
    #many commands treat these as missing.
    #Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.


##CHECK ALL OF THIS

#sex-check plink, compare reported sex with genetics
#https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
#https://www.biostars.org/p/218520/

#PIENSA CASOS IN .HH FILES OF BOTH BATHCES
    #son casos que parecen tener Xx?
    #SI HICIERAN FALTA TENDRIAS QUE CORRER AMBOS BATCHS SIN ELEMINAR HH, Y EVITANDO CORTAR EL SEGUNDO CON EL ERROR

##nosex
    #List of samples with ambiguous sex codes
    #https://www.cog-genomics.org/plink/1.9/output

'''
#see the file
nosex_plink = pd.read_csv("data/plink_inputs_example/batch1_example_plink.nosex",
    names=["FID", "ID"],
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
print(nosex_plink)

#check
print("The IDs in nosex are the same than those IDs of the fam file with unknown sex?")
print(all(nosex_plink["ID"].isin(fam_file.loc[fam_file["Gender"]=="0", "ID"])))
print(all(fam_file.loc[fam_file["Gender"]=="0", "ID"].isin(nosex_plink["ID"])))

'''

##check hetero cases that should be homo

'''
##.hh
    #Produced automatically when the input data contains heterozygous calls where they shouldn't be possible (haploid chromosomes, male X/Y), or there are nonmissing calls for nonmales on the Y chromosome.
    #A text file with one line per error (sorted primarily by variant ID, secondarily by sample ID) with the following three fields:
        #Family ID
        #Within-family ID
        #Variant ID
    #https://www.cog-genomics.org/plink/1.9/formats#hh

#see the file
hh_plink = pd.read_csv("data/plink_inputs_example/batch1_example_plink.hh",
    names=["FID", "ID", "snp_name"],
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
print(hh_plink)

#check
print("heterozygous calls are caused by samples with problematic sex?")
print(unmatched_sex["AGRF code"].isin(hh_plink["ID"]))
print(hh_plink["ID"].isin(unmatched_sex["AGRF code"]))

#see the cases
#for each heterozygous and problematic case
for row_index, hh_case in hh_plink.iterrows():

    #see case
    print("###################################")
    print("Problematic case number " + str(row_index))
    print("###################################")
    print(hh_case)

    #See gender
    print("##########\nGender of the case according to sample map:\n##########")
    gender_case = fam_file.loc[fam_file["ID"] == hh_case["ID"], ["ID", "Gender"]]
    print(gender_case)

    #see genotype
    print("##########\nGenotype according to lgen file:\n##########")
    lgen_case = lgen_file.loc[(lgen_file["Sample ID"] == hh_case["ID"]) & (lgen_file["SNP Name"] == hh_case["snp_name"]), :]
    print(lgen_case)

    #see chromosome
    print("##########\nChromosome according map file:\n##########")
    map_case = map_file.loc[map_file["Name"] == hh_case["snp_name"], :]
    print(map_case)

    #check
    print("##########\nthe problematic case is male and X chromosome, but it has two genotypes, which is not possible for a male?\n##########")
    print(
        (gender_case["Gender"].values == "1") & 
        (map_case["Chromosome"].values == "X") & 
        (~lgen_case[["Allele1 - Forward", "Allele2 - Forward"]].isna().values).all())

'''



##then PCA with pheno to check correlation pheno - genetic structure?
#you can use your plotting approach in python to see pheno and genetic structure

#if you are going to use pheno data, clean the file using a different script!!!
    #change the name of the 2399LDJA for 2397LDJA in the excel file with phenotype data.
        #if we change the ID in illumina, we would have to change the FinalReport, SampleMap, sample_sheet.... and the new data we receive from illumina (IDAT files of the first batch) would need to be changed also. Therefore, it is much more complicated.
    #check if the fact you did not change dtype of week 8 test beep to float in script 1, could be a problem
        #save pheno_data after the cleaning?
    #CHECK THE QUESTIONS TO DAVID about pheno


#load pheno data, this include reported sex and VO2 max data. I have checked that the data is the same directly reading from excel than converting to csv
pheno_data = pd.read_excel(
    "data/pheno_data/combact gene DNA GWAS 23062022.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data)

#there are several phenotypes, see 01a_illumina_report_to_plink_DEV.py for details about possible errors.
#here we are only going to modify an entry that is clearly wrong and do some checks with the sample map. Also, we are going to use the sex indicated in pheno_data because there are some samples with differences between illumina estimated sex and the one showed in the pheno_data

#beep test 8 has an error, "o" letter instead "0" number in one sample
print("\n#####################\n#####################")
print("Week 8 beep test for a sample is 11.1O, i.e., letter O instead number 0")
print("#####################\n#####################")
import numpy as np
index_problematic_sample = np.where(pheno_data["Week 8 beep test"] == "11.1O")[0][0]
index_problematic_column = np.where(pheno_data.columns == "Week 8 beep test")[0][0]
print(pheno_data.iloc[index_problematic_sample, index_problematic_column])
print(pheno_data.iloc[index_problematic_sample,:])

#change 11.1O for 11.10
print("\n#####################\n#####################")
print("error solved")
print("#####################\n#####################")
pheno_data.iloc[index_problematic_sample, index_problematic_column] = 11.1
print(pheno_data.iloc[index_problematic_sample,:])

#change the type of phenotype that is not float but it should be
pheno_data["Week 8 beep test"] = pheno_data["Week 8 beep test"].astype("float64")
    #this column had a row with 11.1O, i.e., letter "O" instead of number "0", so python did not consider this column as float.

#calculate differences
pheno_data["body_mass_diff"] = pheno_data["Week 8 Body Mass"] - pheno_data["Week 1 Body Mass"]
pheno_data["beep_test_diff"] = pheno_data["Week 8 beep test"] - pheno_data["Week 1 Beep test"]
pheno_data["vo2max_diff"] = pheno_data["Week 8 Pred VO2max"] - pheno_data["Week 1 Pred VO2max"]

#merge genetic and phenotypic data
pca_pheno = pd.merge(
    right=first_pca,
    left=pheno_data,
    how="right",
    right_on="IID",
    left_on="AGRF code")
    #merge genetic data with phenotypes, selecting only those samples for which we have genetic data, now we are interested in cleaning genetic data
        #https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html

#
samples_lost_in_pheno = pca_pheno.loc[pca_pheno["AGRF code"].isna(), "IID"]
print("All IDs in illumina data are in pheno data? " + str(samples_lost_in_pheno.shape[0] == 0))
print("How many? " + str(samples_lost_in_pheno.shape[0]))
print("It is ok to have False in this check, because we already knew that some samples had illumina report but where not included in pheno_data. The output script for the first batch shows that 1 sample has genetic data but it is not present in pheno_data (2397LDJA). In the second batch, the output script shows 6 samples with genetic but no pheno_data. In the PCA, we only have one 1 missing sample, so I guess the 6 missing samples of the second batch were those duplicated, i.e., those named as ID_1 and ID_2....")



#CHECK STRAND BEFORE IMPUTATON
#PHASING BEFORE IMPUTATION
    #check ritchie paper



############run case-control study to check for batch effects once you have pre-imputation QC done
    #Another method involves coding case/control status by batch followed by running the GWAS analysis testing each batch against all other batches. For example, the status of all samples on batch 1 will be coded as case, while the status of every other sample is to be coded control. A GWAS analysis is performed (e.g., using the --assoc option in PLINK), and both the average p-value and the number of results significant at a given threshold (e.g., p <1 × 10-4) can be recorded. SNPs with low minor allele frequency (i.e., <5%) should be removed before this analysis is performed to improve the stability of test statistics. This procedure should be repeated for each batch in the study. If any single batch has many more or many fewer significant results or has an average p-value <0.5 (under the null, the average p-value will be 0.5 over many tests), then this batch should be further inves tigated for genotyping, imputation, or compo sition problems. If batch effects are present, methods like those employed for population stratification (e.g., genomic control) may be used to mitigate the confounding effects.
