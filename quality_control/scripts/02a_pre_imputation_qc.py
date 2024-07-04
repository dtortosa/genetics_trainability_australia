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



############################################
# region summary steps #####################
############################################

#0. Removal of duplicated SNPs:
#1. MAF, missingnes and initial HWE filtering
    #The most important argument I have found to do MAF before sample relatedenss is the following:
        #The methods we are gonna use to remove related samples (KING-robust) in plink2 seems to require decent MAF. See this from plink2 help:
            #"The relationship matrix computed by --make-rel/--make-grm-list/--make-grm-bin can be used to reliably identify close relations within a single population, IF YOUR MAFS ARE DECENT". However, if you continue reading, they say that "Manichaikul et al.'s KING-robust estimator can also be mostly trusted on mixed-population datasets (with one uncommon exception noted below), and doesn't require MAFs at all. Therefore, we have added this computation to PLINK 2, and the relationship-based pruner is now based on KING-robust."
                #The exception is that KING-robust underestimates kinship when the parents are from very different populations. You may want to have some special handling of this case; --pca can help detect it.
                #We should not have parents here....
        #Therefore, we do not need good MAFs, but we are doing it after MAF filtering still because we have already prepared the script for that and it should not be a problem, and Ritchie do it after maf filterint (see Figure 5). 
    #So we are going to make an initial clean of SNPs and then do sample removals. If we clean by MAF and then clean by HWE there is no change in the MAF, because we are removing SNPs, not samples.
#2. Sample call rate and relatedness:
    #We remove related samples after we have decent MAF data (for KING) and done some SNP cleaning and then we can remove samples with low genotyping rate. The previous removal of samples is not going to influenec the genotypiing rate of a sample because we removed samples, not SNPs...
    #We can do this before pop stratification becuase we are using King robust
#3. MAF loop to ensure we have OK MAF and genotyping rate
    #We repeat in two steps the MAF-missing snps filters and then sample filter to check that the previous removal of samples did not change allele frequencies in a way that after applying MAF + missing filters again, we lose more samples.
#2. Population stratification: It is strongly recommended by Plink´s author to check subgroups and then perform checks like sex-imbalances inside each group. Besides this, we will use the PCAs as covariates in the analyses.
#4. HWE
    #This is recommnded to be done also withjing ancestry groups
        #After you have a good idea of population structure in your dataset, you may want to follow up with a round of two-sided –hwe filtering, since large (see Note 5) violations of Hardy–Weinberg equilibrium in the fewer-hets-than-expected direction within a subpopulation are also likely to be variant calling errors; with multiple subpopulations, the –write-snplist and –extract flags can help you keep just the SNPs which pass all subpopulation HWE filters.
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
#5. Sex determination. The differences in allele frequencies between ancestry groups can influence the check for sex imbalances, so we have to do it after population stratification analyses. This should be ok, because the problem with sex imbalances could be that the phenotype of one sample is ineed the phenotype or other sample, i.e., they are swapped, but it should influence if we just do analyses with only genotypes like the PCA.
#6. MAF-Missing loop
    #We ahve removed SNPs (HWE) and samples (sex problems), this MAFs and genotyping rates can be changed, so we have to ensure again that MAF and geno rates are ok

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


#endregion

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
#Therefore, we need around 100K independent autosomal SNPs.
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
    #we should the number of samples of each batch multiplied by the number of columns miuns FID and IID (i.e., the PCs)
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


pca_eigenvalues["prop_variance"] = pca_eigenvalues[0]/pca_eigenvalues[0].sum()
pca_eigenvalues

#The proportion of total variance explained by each PC is a useful metric for understanding structure in a sample and for evaluating how many PCs one might want to include in downstream analyses. This can be computed as λi∕∑kλk, with λi being eigenvalues in decreasing order
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #The formula λi∕(∑k λk) is used in the context of **Principal Component Analysis (PCA)**, a statistical procedure that uses an orthogonal transformation to convert a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables called principal components.
    #Here's what each part of the formula means:
        #- **λi**: This represents the **i-th eigenvalue** of the covariance matrix of the dataset. Eigenvalues are a special set of scalars associated with a linear system of equations (i.e., a matrix equation) that are sometimes also known as characteristic roots, characteristic values, proper values, or latent roots.
        #- **∑k λk**: This is the **sum of all eigenvalues** of the covariance matrix. 
        #The fraction $$\frac{\lambda_i}{\sum_k \lambda_k}$$ thus represents the **proportion of the total variance** in the data that is explained by the i-th principal component. This is often used to understand the importance of each principal component in the analysis and decide how many principal components to retain for further analysis.
        #In the sentence you provided, this value is being plotted, likely to create a **scree plot**. A scree plot is a simple line segment plot that shows the fraction of total variance in the data as explained or represented by each PC. This is useful to determine the appropriate number of principal components to retain for further analysis. The plot helps to visualize the **explained variance** by each principal component and it typically decreases and eventually becomes flat with increasing components, which can help to choose the optimal number of components. 





print("plot the eigenvalues")
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

#The results show that the first PC explains a big fraction of the variance (46%). The next 6 (until PC7) explain between 7.5 and 2.4% of variability. After the 7, the variance explained per PC becomes relatively constant.




print("see the first PC where cumulative explianed variability reaches 80% of the total")
#According to Ritchie´s review, it is recommended to use the eigenvectors that explain the greatest proportion of variance (on an average cumulative 80% variance) as covariates in any downstream analysis to correct for bias due to population stratification.
#sum all variability explained and calculate the 80%
cumulative_variance_80=(pca_eigenvalues.sum()*0.8)[0]

#for each eigenvalue
#pc="P1"; row=pca_eigenvalues.iloc[0,:]
import numpy as np
for pc_short, row in pca_eigenvalues.iterrows():

    #if the row is not the last one, i.e., the las PC
    if(pc_short!="P"+str(pca_eigenvalues.shape[0])):

        #get the position of the current row
        current_index = np.where(pca_eigenvalues.index == pc_short)[0][0]

        #take all the values from the first to the current row and sum them all
        cumulative_variability = (pca_eigenvalues.iloc[0:current_index+1,:].sum())[0]
            #to get the row of "index", we need need "index"+1 because the lsat number is not included
        
        #if the cumulative variability is equal or lower than the 80%
        if(cumulative_variability>=cumulative_variance_80):
            
            #print the name of PC and stop de loop
            print(pc_short)
            break

#Summary:
#The elbow is at PC5, but the PCs absorbing 80% of variability are 10 first. So we should use the first 10 axes as covariates, but for doing checks and outliers we do not have to use all of them. Ritchies talks about the variance explianed in the context of using the axes as covariates.
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

print_text("outlier removal", header=4)
#It is also a good idea to throw out gross outliers at this point; any sample which is more than, say, 8 standard deviations out on any top principal component is likely to have been genotyped/sequenced improperly; you can remove such samples by creating a text file with the bad sample IDs, and then using –remove +  –make-bed:
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22

    #The phrase "8 standard deviations out" refers to a data point that falls 8 standard deviations away from the mean of a dataset. In statistics, the standard deviation is a measure of the amount of variation or dispersion in a set of values. A low standard deviation means that the values tend to be close to the mean, while a high standard deviation means that the values are spread out over a wider range. When a data point is said to be "8 standard deviations out", it means it's quite far from the mean. In a normal distribution, almost all data falls within 3 standard deviations of the mean. So, a data point that is 8 standard deviations away from the mean is extremely rare and could be considered an outlier. In the context of your sentence, it suggests that any sample which falls more than 8 standard deviations away on any top principal component might have been genotyped or sequenced improperly, and thus, it might be a good idea to remove these outliers from the analysis. This is because such extreme values can significantly skew the results and interpretations of the data analysis.


#remove outliers based on the first 10 axes
#pc="PC2"
list_df = list()
for pc in pc_columns:

    if(pc!="PC11"):
        high_limit = pca_results[pc].mean()+(pca_results[pc].std()*8)
        low_limit = pca_results[pc].mean()-(pca_results[pc].std()*8)

        pca_results_subset=pca_results.loc[(pca_results[pc]>high_limit) | (pca_results[pc]<low_limit), ["#FID", "IID"]]
        
        pca_results_subset=pca_results_subset.assign(PC=pc)

        list_df.append(pca_results_subset)
    else:
        break

concatenated_df = pd.concat(list_df)

concatenated_unique_df = concatenated_df.drop_duplicates(subset=["#FID", "IID"])

concatenated_unique_df.shape[0]
concatenated_unique_df

concatenated_unique_df[["#FID", "IID"]].to_csv("./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_outliers.txt",
    sep="\t",
    header=False,
    index=False)


run_bash(" \
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca; \
    plink \
        --bfile ./loop_maf_missing_2_ldprunned_autosome_pca \
        --remove ./pca_outliers.txt \
        --make-bed \
        --out ./loop_maf_missing_2_ldprunned_autosome_pca_not_outliers"
)
    #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.


run_bash("\
    cd ./data/genetic_data/quality_control/14_pop_strat/01_pca; \
    plink2 \
        --pca 20 \
        --bfile ./loop_maf_missing_2_ldprunned_autosome_pca_not_outliers \
        --out ./pca_results_no_outliers")

pca_results_no_outliers = pd.read_csv(\
    "./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_results_no_outliers.eigenvec", \
    sep="\t", \
    header=0, \
    low_memory=False)
pca_results_no_outliers


second_round_outliers = list()
for pc in pc_columns:

    if(pc!="PC11"):
        high_limit = pca_results_no_outliers[pc].mean()+(pca_results_no_outliers[pc].std()*8)
        low_limit = pca_results_no_outliers[pc].mean()-(pca_results_no_outliers[pc].std()*8)

        pca_results_subset=pca_results_no_outliers.loc[(pca_results_no_outliers[pc]>high_limit) | (pca_results_no_outliers[pc]<low_limit), ["#FID", "IID"]]
        
        pca_results_subset=pca_results_subset.assign(PC=pc)

        list_df.append(pca_results_subset)
    else:
        break


if(len(second_round_outliers)!=0):
    raise ValueError("ERROR! FALSE! WE STILL HAVE PCA OUTLIERS")



import seaborn as sns
sns.pairplot(pca_results[pc_columns[0:5]])
plt.savefig("./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_results_no_outliers.png", dpi=300)
plt.close()

#If there are obvious clusters in the first few plots, I recommend jumping ahead to Chapter 4 (on ADMIXTURE) and using it to label major subpopulations before proceeding.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22


#PLAN:
    #Base the stratificaciton analysis using admixture tutoriral by Jhon Novembre
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #From there, took info to improve the PCA analyses of plink (althoug he used eigensoft). For example, calculate teh explained variance as they do, we did it wrong In think. Check the removal of outliers and the selection of PCs. Think if eigensoft is needed or not.
        #from here, keep in mind the number of PC axis that are relevant because we will check if that matches the number of subpops predicted by admixture. maybe eigensoft is worth it for this because it calculate a p-value for PCs.
    #then, we follow the tutorial in admixtuure doing several steps to define the subpopulations, selecting the best number of them. Admixture label them (.Q file), so we can split the next two QC steps within each subpop.
        #think it is worth to merge our data with the 1KGDP to use as anchor. If you do that you have to check for strand flips. see admixture tutorial.





#Marees uses MDS analysis form plink, remove outliers, repeat mds and use the first axes in the glms..

#look this cool tutorial, they use anchoring and MDS slike marees 
    #https://gwas-intro-bajicv.readthedocs.io/en/latest/05_pop_stratification/

#eigensof tor ritchie?
#Occasionally, the new principal components will reveal another bad sample, and you have to repeat these two steps, etc. EIGENSOFT [7, 8] has some additional built-in principal component analysis options, including automated iterated outlier removal, and a top-eigenvalue-based test for significant population structure.)
#I have downloaded Eigensoft, the latest version at the time of writting (3/07/2024), which is the release 8.0.0:
    #https://github.com/DReichLab/EIG/archive/refs/tags/v8.0.0.tar.gz
    #https://github.com/DReichLab/EIG/releases


###follow plink tutorial
#then eigensotf or admixustre? loook ritchie and vo2 paper...
import plotly.express as px
fig = px.scatter(\
    data_frame=pca_results_no_outliers, \
    x="PC1", \
    y="PC2", \
    color="#FID",
    hover_data=[\
        "IID"])
        #you can use columns of DF to add axis data, but also modify color, size, and show data per sample in a desplegable box
        #https://plotly.com/python/line-and-scatter/
        #https://plotly.com/python-api-reference/generated/plotly.express.scatter.html
#fig.show()
fig.write_html("./data/genetic_data/quality_control/14_pop_strat/01_pca/pca_results_no_outliers.html")
    #https://plotly.com/python/interactive-html-export/







#SELECT THE AXES BASED ON EXPLAINED VARIANCE
    #you can use eigensoft


#it seesm you have to use the PCA axes in the glms even if there is not different ancestries becasue you coudl still ahve structure within an homogeneous pop, i think this is from ritchies...

#It is also a good idea to throw out gross outliers at this point; any sample which is more than, say, 8 standard deviations out on any top principal component is likely to have been genotyped/sequenced improperly; you can remove such samples by creating a text file with the bad sample IDs, and then using –remove +  –make-bed:

    #SNPs with Minor Allele Frequency (MAF)>0.05 were then used to perform principal component analysis (PCA) for ethnicity identification using SHELLFISH [45]. Ethnic and ancestry outliers (more than 6 standard deviations from the mean on either of the two first principal components (PCs)) were excluded (n=10). 
        #paper VO2 max

#plink recommend to remove outlier indicaitng genotyping errors and then check for clusters, it ehse are present, then perofmr a promer admixutre analysis to define the groups... or maybe just use eigensoft to define significant groups?

##we have to decide whetehr to maintain the groups or to remove ancestry oultiers like trainibiltiy paper
    #it depends of the number lost....
    #Different ethnicities can be included in the same study, as long as the population substructure is considered to avoid false positive results
        #general tutorial gwas
    #Christopher:
        #Population structure is typically accounted for by using top principal components (computed via plink --pca, or a similar function in another software package) as covariates.  Unfortunately, plink 1.07's haplotype association analysis does not support covariates; but the main allelic association command (--linear/--logistic, or --glm in plink 2.0) does support them.
            #https://groups.google.com/g/plink2-users/c/938B07i8AXQ/m/dMFI8-GLAwAJ



#If you do this, follow it up by repeating the PCA, since the bad samples might have distorted the principal components:

#(Occasionally, the new principal components will reveal another bad sample, and you have to repeat these two steps, etc. EIGENSOFT [7, 8] has some additional built-in principal component analysis options, including automated iterated outlier removal, and a top-eigenvalue-based test for significant population structure.)

#Anyway, once there is nothing obviously wrong with the PCA results, you can load the table in R and plot the top pairs of principal components against each other:

#If there are obvious clusters in the first few plots (i.e., first PCs plotted against each other), I recommend jumping ahead to Chapter 4 (on ADMIXTURE) and using it to label major subpopulations before proceeding.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#author-information
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4

#but how much is this? we should check in the conext of 1000 KGP....




#THINK ABOUT THIS, because maybe we have to do imptuation in each group separately?
    #michina imputation sever says we have to Choose a reference panel
    #it seems that the reference panel is something like hapmap or the 1KGP, so it shoudl ahndle different ancestries, but check
        #A simple alternative is to use a “cosmopolitan” reference set that includes all available haplotypes, each of which is assigned an equal chance of being copied a priori. This approach produces relatively accurate results in a variety of human populations and has therefore been proposed as a good fallback choice when the optimal panel composition is unclear (Guan and Stephens 2008; Huang et al. 2009a; Li et al. 2010). Another class of methods tries to maximize accuracy by weighting reference panels through cross-validation (Huang et al. 2009a) or ancestry estimation (Egyud et al. 2009; Pasaniuc et al. 2010); the Pasaniuc et al. approach differs from the others in that it uses local ancestry estimates to provide customized reference weights for each study individual. As an alternative, Jostins et al. (2011) suggested balancing accuracy and computation by using reference panels that “approximately cluster” with the study individuals on a plot of principal components (PC) that capture genetic ancestry.
        #https://academic.oup.com/g3journal/article/1/6/457/5986469



##USE OTHER TECHINES TO CHECK FOR OUTLIERS AND POP STRATIFICATION?
    #look tutorials
    #https://www.cog-genomics.org/plink/1.9/strat


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





# endregion




##########################
# region HWE SECOND ROUND
###########################


#After you have a good idea of population structure in your dataset, you may want to follow up with a round of two-sided –hwe filtering, since large (see Note 5) violations of Hardy–Weinberg equilibrium in the fewer-hets-than-expected direction within a subpopulation are also likely to be variant calling errors; with multiple subpopulations, the –write-snplist and –extract flags can help you keep just the SNPs which pass all subpopulation HWE filters.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22

# endregion




######################################################
# region heterozygosity wihting ancestry groups ######
######################################################

#Ritchie's tutorial says "in downstream analyuses .... less than expected heterozygosity (mean – 3 SD) suggests possible inbreeding and greater than expected heterozygosity (mean + 3 SD) suggests possible sample contamination. However, these thresholds should take into account the expected heterozygosity rates in the ancestry group under study, as some diverse populations may exhibit different rates of heterozygosity than other populations."

##important
#check if you have to use prunned LD data

#do plot hetero - missingness
#use the information to remove samples due to contamination (high hetero) or inbreeding (low hetero)

# endregion




##################################
# region SEX INCONSISTENCES ######
##################################


    #check sex uses allele frequencies, so you can have problems with different ancestries, yo have to prepare the data... so maybe is bettter to see the PCA for outliers and batch effects, filtering and then go to see sex once we have cleaner data
        #"Due to the use of allele frequencies, if your dataset has a highly imbalanced ancestry distribution, you may need to process the rare-ancestry samples separately."
        #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex

    #check-sex has to be used on LD-pruned data, so we can use the previous data
        #Since this function is based on the same F coefficient as --het/--ibc, it requires reasonable MAF estimates (so it's essential to use --read-freq if there are very few samples in your immediate fileset), and it's best used on marker sets in approximate linkage equilibrium.
    #This isn't implemented yet in plink2 since there would be little practical difference from the plink 1.9 implementation.  Use "--make-bed --chr X,Y" to export a .bed fileset with only chrX and chrY, and run plink 1.9 --check-sex on that.

    #BUT OF COURSE WE NEED DATA OF SEX CHROMOSOMES, select the prunnn ed dataset with sex chromosomes

#sex inconsistences on ld_pruned
    #In the tutorial of plink Chirstopher checks sex on the prunned dataset
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec18
    #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex



#By default, –check-sex/–impute-sex and several other PLINK commands assume that your dataset is representative of a single population, and all samples are members of that population; population allele frequencies and F coefficients are estimated under these assumptions.

#Thus, if you identified multiple subpopulations in the previous section, you should perform sex validation/imputation on one subpopulation at a time. PLINK’s –filter flag provides one way to do this; put sample IDs in the first two columns of subpops.txt and subpopulation IDs in the third column, then

#Plink says that, due to the use of allele frequencies we may need to check sex within ancestry groups, and he did that in the tutorial, so we are doing this after pop structure analysis
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec19


run_bash(" \
    cd ./data/genetic_data/quality_control/10_remove_low_maf_snps; \
    plink \
        --bfile ./merged_batches_hwe_par_callrate_LogRDev_maf \
        --check-sex; \
    ls -l")


run_bash(" \
    cd ./data/genetic_data/quality_control/10_remove_low_maf_snps/; \
    awk \
        'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1, $2, $3, $4, $5, $6}'\
        plink.sexcheck > plink.sexcheck_awk_processed.tsv; \
    head plink.sexcheck_awk_processed.tsv")


sex_check_report = pd.read_csv( \
    "./data/genetic_data/quality_control/10_remove_low_maf_snps/plink.sexcheck_awk_processed.tsv", \
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
    #we should change the name of this sample, probably easier to do it in the excel...

#OJO MINORITIES
    #the results seems to make sense with most females below 0.2, CHECK THAT!
    #we need to check ancestry? see ritchie


#https://groups.google.com/g/plink2-users/c/4bpdLMdH2KA
    #If heterozygous haploid calls still remain, the most likely cause is nonmissing female genotype calls on the Y chromosome; others have reported that this is fairly common.  A quick way to check the number of these is to just load the Y chromosome with e.g. "plink --bfile semi_clean_fileset --chr 24 --freq".  If all the heterozygous haploid errors are on the Y chromosome, you can safely clobber them with --make-bed + --set-hh-missing.  (If some are on the X, --set-hh-missing *might* still be okay, but I'd need to know more about the data source and the --check-sex report to be sure.)


# endregion




#####################################
# region LAST MAF-MISSING LOOP ######
#####################################

#We repeat in two steps the MAF-missing snps filters and then sample filter to check that the previous removal of samples did not change allele frequencies in a way that after applying MAF + missing filters again, we lose more samples.

#do it within each ancestry group?












#check differences in pheno between batches?
#do case-control study for batch effects AFTER all pre-QC steps?




###imputation separated between ancestry groups?


# endregion




#############################
# region CODE TO CHECK ######
#############################



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

# endregion

