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
    #2. Data Management and Summary Statistics with PLINK
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
    #3. Genomics Boot Camp
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



#########################################
# rationale of merging after imputation #
#########################################

#From reading the recent GWAS protocol of the Ritchie lab. I understand that the usual approach is to do the imputation of different sources (e.g., batches) separately and then do the merging. In the section "Batches effect", they say that "Genotype imputation strategies now provide an opportunity to impute genotypes from multiple platforms to a reference genotyping platform. THEREFORE, DATA FROM DIFFERENT SOURCES CAN BE MERGED AFTER IMPUTATION."
    #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view?usp=sharing

#In addition, in other paper they say 
    #"This process (imputation) is particularly important when combining or performing meta-analysis on data generated using multiple different genotyping platforms"
    #"After imputation and merging of the datasets, quality control procedures were implemented to create high quality, analysis-ready data set for genome-wide association studies."
    #https://www.frontiersin.org/articles/10.3389/fgene.2014.00370/full

#I understand that specially when you have data from different sources, it is very important to do imputation because many SNPs are NOT going to be shared across sources, you are going to lose them unless you impute them to have the same SNPs across all panels.

#In our particular case, this should not be very important because the same SNPs are genotyped in both batches, so it is very unlikely that all samples of a batch have missing for a given SNP so that SNP is not present in one batch but it is present in the other one. 

#I going to do imputation first just in case, to avoid any problems. Also, remember that the protocol mentions the possibility to merge after imputation the batch section, so it does not seem very strange to merge batches after imputation. Indeed, I have seen a post in biostarts doing PCAs and other stuff in different batches and then in the merged dataset.
    #https://www.biostars.org/p/438079/



#############################################################
# final decision about the 1100JHJM, 1200JPJM and 7800AGSO #
#############################################################

#the three duplicated IDs
    #For 1100JHJM and 1200JPJM, David's postdoc agrees we should remove them as we have the same ID for two different rows in the excel file.
        #"This is correct. In session 844, there are two individuals with the same code (1100JHJM, both male), and in session 845, there are two individuals with the same code (1200 JPJM, both male). I agree that we will need to remove these from our analysis."
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



#######################################
# Passing arguments of python program #
#######################################

#define input arguments to be passed when running this script in bash
#we will use a bash script to run this python program two times instead of doing a function and parallelize that function. In this way, the same python script is run for each batch and if we want to run again one batch, we just need to modify the bash script, not the python program. This makes sense because we have only two batches, and the parallelization will occur inside each batch. Importantly, we can also have separated .out files, so we can look at the checks separately.
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--batch_name", type=str, default="ILGSA24-17873", help="Name of the batch used as input. Always string.")
parser.add_argument("--n_cores", type=int, default=4, help="Number of cores/threads requested. Integer always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
batch_name = args.batch_name
n_cores = args.n_cores



###########################
# starting with the batch #
###########################
print_text("starting batch number " + batch_name + " using " + str(n_cores) + " cores ", header=1)


print_text("do some checks on the plink files of the batch", header=2)

print_text("get the total number of samples per batch:", header=3)
if batch_name=="ILGSA24-17873":
    total_samples=1242
elif batch_name=="ILGSA24-17303":    
    total_samples=216
print("we should have 216 and 1242 samples for batch 1 and 2 (batch 2 lost 6 samples that were duplicated), respectively")
print(f"For batch {batch_name}, we have {total_samples} total samples")


print_text("check we have the correct number of samples looking at the merged fam file and the list of samples considered to merge for each batch", header=3)
run_bash(" \
    n_samples_batch=$( \
        gunzip -c ./data/genetic_data/plink_bed_files/" + batch_name +  "/03_merged_data/" + batch_name + "_merged_data.fam.gz | \
        awk -F '\t' 'END{print NR}'); \
    if [[ $n_samples_batch -eq " + str(total_samples) + " ]]; then \
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
run_bash(" \
    n_samples=$( \
        awk 'END{print NR}' ./data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/list_to_merge.txt); \
    if [[ $n_samples -eq " + str(total_samples) + " ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #count the number of rows in the file with the list of samples considered in the merging of both batches. We use NR (number of records) of awk instead of wc -l because the file does not end with an empty line, so wc does not count the last one, and we get 1 row less
            #https://stackoverflow.com/questions/32481877/what-are-nr-and-fnr-and-what-does-nr-fnr-imply
            #https://stackoverflow.com/questions/12616039/wc-command-of-mac-showing-one-less-result


print_text("check we do NOT have a .missnp file, where variants with more than two alleles are indicated. We should not have this file", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/" + batch_name +  "/03_merged_data; \
    count_miss_file=$(ls -l *.missnp | wc -l); \
    if [[ $count_miss_file -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #list files with extension ".missnp" and count them. This should be zero.


print_text("check we have 654027 snps", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/" + batch_name +  "/03_merged_data; \
    n_snps=$( \
        gunzip \
            --stdout \
            ./" + batch_name + "_merged_data.bim.gz | \
        awk \
            -F '\t' \
            'END{print NR}'); \
    if [[ $n_snps -eq 654027 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")


print_text("check that the bim files (SNP maps) of both batches are identical, at least the first 4 columns, because the allele that is minor can be different", header=3)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/; \
        gunzip \
            --stdout \
            ./ILGSA24-17303/03_merged_data/ILGSA24-17303_merged_data.bim.gz | \
        awk \
            -F '\t' \
            '{print $1, $2, $3, $4}' > \
        batch_1_columns.txt; \
        gunzip \
            --stdout \
            ./ILGSA24-17873/03_merged_data/ILGSA24-17873_merged_data.bim.gz | \
        awk \
            -F '\t' \
            '{print $1, $2, $3, $4}' > \
        batch_2_columns.txt; \
    bim_check=$(cmp --silent batch_1_columns.txt batch_2_columns.txt; echo $?); \
    if [[ $bim_check -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi; \
    rm batch_1_columns.txt; \
    rm batch_2_columns.txt")
        #select the first 4 columns of the three bim files using awk (-F is the delimiter, tab in this case). Save the result as three different .txt files.
            #the 4 first columns are the chromosome code, Variant identifier, Position in morgans or centimorgans (safe to use dummy value of '0'), Base-pair coordinate (1-based).
            #these columns should be the same in both batches, the difference can be in the last two columns with the name of the minor and major alleles, because it would be possible that one allele is not present in any of the samples of a batch, but in the other batch it has both alleles (see below).
        #compare the two batches between them using cmp for that.
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

#note about the differences in minor/major alleles between batches
    #the four first columns are identical between batches, but the minor/major alleles are different.
    #the first difference between the two batches is a SNP in row 140. It has only C in first batch (monomorphic), but it has also T in the second. Accordingly the merged batch has both alleles, because this included first and second batch
        #first batch
            #0       GSA-10:47138541 0       0       0       C
        #second batch
            #0       GSA-10:47138541 0       0       T       C
    #opposite case with a snp in row 165, which is the first line of difference between merged (a merge file i did with both batches) and the second batch. It is monomorphic in second batch, but it has two alleles in first batch, so when merging, you get two alleles, not 1.
        #first batch
            #0       GSA-rs145797772 0       0       A       C
        #second batch
            #0       GSA-rs145797772 0       0       0       C
        #merged file
            #0       GSA-rs145797772 0       0       A       C


if batch_name == "ILGSA24-17873":
    print_text("check that we do NOT have the duplicate samples (1100JHJM, 1200JPJM and 7800AGSO) in the second batch", header=3)
    run_bash("\
        cd ./data/genetic_data/plink_bed_files/" + batch_name +  "/03_merged_data; \
        n_no_dups=$( \
            gunzip \
                --stdout \
                " + batch_name + "_merged_data.fam.gz | \
            awk \
                -F '\t' \
                '{ \
                    if($2!=\"1100JHJM_1\" && $2!=\"1200JPJM_1\" && $2!=\"7800AGSO_1\" && $2!=\"1100JHJM_2\" && $2!=\"1200JPJM_2\" && $2!=\"7800AGSO_2\"){\
                        count ++\
                    } \
                } \
                END {print count}'); \
        total_n=$( \
            gunzip \
                --stdout \
                " + batch_name + "_merged_data.fam.gz | \
            awk \
                -F '\t' \
                'END{print NR}'); \
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
    gunzip \
        --stdout \
        ./data/genetic_data/plink_bed_files/" + batch_name +  "/03_merged_data/" + batch_name + "_merged_data.bim.gz | \
    awk \
        -F '\t' \
        'END{print NR}'", return_value=True).strip()
n_snps = int(n_snps)
print(n_snps == 654027)


print_text("check SNPs IDs are not duplicated, as plink 1.9 does not look for duplicates using the ID", header=3)
run_bash(" \
    n_snp_ids_non_dup=$( \
        gunzip \
            --stdout \
            ./data/genetic_data/plink_bed_files/" + batch_name +  "/03_merged_data/" + batch_name +  "_merged_data.bim.gz | \
        awk \
            -F '\t' \
            '!a[$2]++{count++} \
            END{print count}'); \
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
        "./data/genetic_data/plink_bed_files/" + batch_name + "/03_merged_data/" + batch_name + "_merged_data.bim.gz", \
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


print_text("load the duplicates list obtained in '01b_illumina_report_to_plink.py' with the options '        --list-duplicate-vars suppress-first ids-only'. Therefore, for each group of SNPs with the same position and alleles, the first one is not included in the list, but the rest are and these are the one that will be removed", header=3)
duplicate_cases = pd.read_csv(
    "./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates.dupvar", 
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


print_text("calculate the number of duplicates with pandas and check is the same than plink duplicates. As we want to consider as duplicates SNPs with the same alleles, irrespectively of the assigment, we need to create a new column where we get both alleles in the same order independenlty if allele 1 is A or T. Therefore A1=A and A2=T will be the same than A1=T and A2=A, i.e., AT.", header=3)
bim_file_dup_check = pd.read_csv( \
        "./data/genetic_data/plink_bed_files/" + batch_name + "/03_merged_data/" + batch_name + "_merged_data.bim.gz", \
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
bim_file_dup_check["new_alleles"] = bim_file_dup_check\
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


print_text("open a folder to save filtered dataset", header=3)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    rm --recursive --force ./04_inspect_snp_dup/01_remove_dup/; \
    mkdir --parents ./04_inspect_snp_dup/01_remove_dup/; \
    ls ./04_inspect_snp_dup")
        #rm 
            #--recursive: remove directories and their contents recursively
            #--force: ignore nonexistent files and arguments, never prompt
        #mkdir
            #--parents: no error if existing, make parent directories as needed


print_text("decompress the plink binary fileset", header=3)
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/03_merged_data/; \
    ls -l; \
    gunzip --keep --force ./" + batch_name + "_merged_data.bed.gz; \
    gunzip --keep --force ./" + batch_name + "_merged_data.bim.gz; \
    gunzip --keep --force ./" + batch_name + "_merged_data.fam.gz; \
    ls -l")
        #-k: keep original compressed file
        #-f: force compression or decompression even if the file has multiple links or the corresponding file already exists


print_text("filter these snps", header=3)
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    plink \
        --bfile ./03_merged_data/" + batch_name + "_merged_data \
        --exclude ./04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates.dupvar \
        --make-bed \
        --out  ./04_inspect_snp_dup/01_remove_dup/" + batch_name + "_merged_data_no_snp_dup")
            #--bfile: 
                #This flag causes the binary fileset plink.bed + plink.bim + plink.fam to be referenced. If a prefix is given, it replaces all instances of 'plink', i.e., it looks for bed, bim and fam files having that suffix instead of "plink".
                #https://www.cog-genomics.org/plink/1.9/input#bed
            #--exclude: Exclude ALL variants named in the file.
                #We have in the ".dupvar" file only the second and next occurrences of each position/allele duplicate, not including the first one. Therefore, we can remove these SNPs.
                #https://www.cog-genomics.org/plink/1.9/filter#snp
            #--make-bed creates a new PLINK 1 binary fileset, AFTER applying sample/variant filters and other operations
                #https://www.cog-genomics.org/plink/1.9/data#make_bed


print_text("remove the non-compressed binary fileset of 03_merged_data")
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/03_merged_data/; \
    ls -l; \
    rm ./" + batch_name + "_merged_data.bed; \
    rm ./" + batch_name + "_merged_data.bim; \
    rm ./" + batch_name + "_merged_data.fam; \
    ls -l")


print_text("check again for duplicates", header=3)
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/; \
    plink \
        --bfile ./" + batch_name + "_merged_data_no_snp_dup \
        --list-duplicate-vars suppress-first ids-only\
        --out  ./" + batch_name + "_duplicates")
        #--list-duplicate-vars
            #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
            #This is not based on variant IDs. Use PLINK 2.0's --rm-dup for ID-based deduplication.
            #By default, this ignores A1/A2 allele assignments, since PLINK 1 normally does not preserve them. If you want two variants with identical positions and reversed allele assignments to not be considered duplicates, use the 'require-same-ref' modifier, along with --keep-allele-order/--a2-allele.
            #Normally, the report has a header line, and contains positions and allele codes in the first 3-4 columns. However, if you just want an input file for --extract/--exclude, the 'ids-only' modifier removes the header and the position/allele columns, and 'suppress-first' prevents the first variant in each group from being reported (since, if you're removing duplicates, you probably want to keep one member of each group).
            #--list-duplicate-vars fails in 'ids-only' mode if any of the reported variant IDs are not unique.


print_text("Do we have zero duplicated positions after filtering?", header=3)
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/; \
    n_lines=$( \
        awk \
            -F '\t' \
            'END{print NR}' \
            " + batch_name + "_duplicates.dupvar); \
    if [[ $n_lines -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")


print_text("remove the files we are not interested in", header=3)
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/; \
    ls -l; \
    rm " + batch_name + "_duplicates.*; \
    ls -l")


print_text("load the bim file after removing duplicates", header=3)
bim_no_dups = pd.read_csv( \
    "./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/" + batch_name + "_merged_data_no_snp_dup.bim", \
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

print_text("see the unique chromosome names in the SNP map of this batch", header=3)
print_text("select the name of the zip file with batch data based on the batch name because ILGSA24-17873 is called CAGRF20093767.zip, not ILGSA24-17873.zip", header=4)
if batch_name=="ILGSA24-17873":
    zip_name = "CAGRF20093767"
elif batch_name=="ILGSA24-17303":
    zip_name = "ILGSA24-17873"
print(zip_name)

print_text("load zip info from the zip file", header=4)
import zipfile
zipdata = zipfile.ZipFile("data/genetic_data/illumina_batches/" + zip_name + ".zip")
zipinfos = zipdata.infolist()
zipinfos

print_text("extract the zipinfo of the SNP_map", header=4)
import numpy as np
#zipinfo=zipinfos[0]
for zipinfo in zipinfos:
    if zipinfo.filename == zip_name + "/SNP_Map.txt":
        zipinfo.filename = zipinfo.filename.split(zip_name+"/")[1]
        zipinfo_snp_map = zipinfo
            #select the zipinfo of SNP after removing the first part with the parent folder name
zipinfo_snp_map

print_text("extract the SNP map", header=4)
zipdata.extract(zipinfo_snp_map, "./data/genetic_data/quality_control/" + batch_name + "/")

print_text("extract the unique chromosome names using awk", header=4)
unique_chr_map = run_bash(" \
    cd ./data/genetic_data/quality_control/" + batch_name + "/; \
    tail \
        -n +2 \
        SNP_Map.txt | \
    awk \
        -F '\t' \
        '!a[$3]++ {print $3}'", return_value=True)
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
            "./data/genetic_data/quality_control/" + batch_name + "/SNP_Map.txt", \
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

print_text("count SNPs for which we have XY chromosome, i.e., pseudo-autosomal", header=4)
chr_XY_bool = bim_no_dups[0]==25
print(f"We have {sum(chr_XY_bool)} SNPs pseudo-autosomal regions")
if sum(chr_XY_bool) < 1000:
    print("we have less than 1000 variants in pseudo-autosomal regions, no problem! add them to the list to SNPs to be excluded")
    chrom_xy_cases = bim_no_dups.loc[chr_XY_bool, 1].to_list()
        #select ID of these cases and convert to list
    [list_snp_ids_to_exclude.append(snp) for snp in chrom_xy_cases]
        #take each of these snps and save it in the result list
    print(list_snp_ids_to_exclude)
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE TOO MUCH VARIANTS WITH CHROMOSOME 0")



#check 25 in plink is XY in illumina snp_map

snp_map = pd.read_csv( \
            "./data/genetic_data/quality_control/" + batch_name + "/SNP_Map.txt", \
            sep="\t", \
            header=0, \
            low_memory=False)

sum(~snp_map.loc[snp_map["Chromosome"] == "XY", "Name"].isin(bim_file_dup_check.loc[bim_file_dup_check[0]==25, 1]))
sum(~bim_file_dup_check.loc[bim_file_dup_check[0]==25, 1].isin(snp_map.loc[snp_map["Chromosome"] == "XY", "Name"]))
    #we have the same PAR SNPs for XY and 25 (code for PAR in plink)
    #I used bim original because we are using the SNP which is not filtered by duplicates, so bim_no_dups would be different

#think if remove mito, par
#I think we can leave for now sex chromosome and mito for QC, indeed we need sex chromosomes in order to do checks about sex
#mito and PAR are only a few thousand SNPs...


#USE EXCLUDE IN PLINK WITH THE LIST TO REMOVE SNPS


#check the list of snps to remove is the sum of each indivudual list of snps to remove
#check all the unwanted chromosomes are removed


#we do not have duplicated IDs (checked above)

#tiene sentido eliminar snps with low-call rate, ancestria no deberia afectar




#do plot hetero - missingness, althouhg remove samples after PCA?
#remove snps and samples by call rate, check your three tutorials
#then removal of samples with PCA and in the next script removal of samples with heterozigosity?


#run again the first script do to removal of compressed files of 03?



#remove SNP map
"./data/genetic_data/quality_control/" + batch_name + "/SNP_Map.txt"








#after removing duplicates, check the text I have in the next lines and then go directly to the protocol, follow it, except not starting with sex and hetero but with PCA, see next comment.
    #https://github.com/RitchieLab/GWAS-QC

#Use GC score to filter?
    #it should be above 0.15 in all cases, right? if not the genotype is empty? They calculated this with GenomeStudio 
    #https://www.illumina.com/Documents/products/technotes/technote_infinium_genotyping_data_analysis.pdf

#marker quality
    #check that the number of missing SNPs is 654027-650181, because the sum of genotypes that pass quality fileters and tjhose not passing the filter in the DNAreport is 650181, not 654027, for the second batch.
        #If that is the case, add it to line "Note about Calls - No_Calls" in the first script
    #in the tutorial they say that "Note that similarly to heterozygosity rates, both minor allele frequencies and deviations from Hardy-Weinberg will be affected by differ- ences in ancestry and should be considered in a population-specific way."
        #therefore I think we should check first PCA
    #note that in that paragraph they are also talking about missing rate per SNP, and they do not mention it here, so I guess it is ok to do it before the PCA, it should not be affected by the pop differen


#sex (--check-sex) and hetero (--het) should be checked after PCA because accorindg to plink info, it can be problems if we have a sample with most of samples from one ancestry and then a few from another ancestry
    #https://www.cog-genomics.org/plink/1.9/basic_stats

    #R log is present in our data, so we could check prob intensity in X for a full detail sex determination, think about it.

    #with sex check done, tell David about mismatches

    #Warning: 30589 het. haploid genotypes present (see
    #./04_inspect_snp_dup/01_remove_dup/ILGSA24-17873_merged_data_no_snp_dup.hh );
    #many commands treat these as missing.
    #Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.




    #think why you have a missing sample (2397LDJA), but when you check the number of samples between illumina and pheno_data, you only have 1 of difference. In the email to Bishop, you said that the missing sample could be the AGO... that was duplicated in illumina but not in pheno_data, so we should have 2 less samples, not 1. look the empty row between data and NAs... 

#filters I have not seen in ritchie's paper
    #filter by chromosome
        #check strange chromosome numbers (non-autosomal)
        #also check that no genetic position is added

    #check if indels!
        #1:207754848-GATAA-G
        #plink has flag  --snps-only to keep snps
            #https://www.biostars.org/p/378475/

#FOR QUALITY CONTROL, YOU COULD USE MIGHIGGAN SERVER, WHICH ALREADY APPLIES MULTIPLE FILTERS like duplicates AND THEN ADD A FEW WITH PLINK, AUGUSTO DID THAT
    #this michigan sever can be used also for imputing
        #https://www.mdpi.com/2073-4425/14/2/248
        #https://imputationserver.readthedocs.io/en/latest/pipeline/





#duplicates
#you can merge snps with the same position merging the fileset with itself and using --merge-equal-pos
    #If two variants have the same position, PLINK 1.9's merge commands will always notify you. If you wish to try to merge them, use --merge-equal-pos. (This will fail if any of the same-position variant pairs do not have matching allele names.) Unplaced variants (chromosome code 0) are not considered by --merge-equal-pos.
    #Note that you are permitted to merge a fileset with itself; doing so with --merge-equal-pos can be worthwhile when working with data containing redundant loci for quality control purposes.












#you could end this script by creating the missingness-heterogizosity rate plot and the PCA
#then use the information of both approaches to remove samples because different ancestry (outlier PCA), contamination (high hetero) or inbreeding (low hetero)
#yes, hetero can be influenced by ancestry, but I am going to use the threshold of the tutorial and if this is not ok for other ancestries, we can check if the outlier hetero are from other ancestries
    #see tutorials

print_text("start PCA to detect individuals that are outliers", header=2)

print_text("prepare folder for PCA", header=3)
run_bash(" \
    cd ./data/genetic_data/quality_control/; \
    ls -l; \
    mkdir \
        --parents \
        " + batch_name + "/pca; \
    ls -lR")
        #ls -R is for recursive


#read the two summary illumina reports once you have both and understand the quality checks done looking at biostars
    #ILGSA24-17303 Report v1.0.pdf


#if you are going to use pheno data, clean the file using a different script!!!
    #change the name of the 2399LDJA for 2397LDJA in the excel file with phenotype data.
        #if we change the ID in illumina, we would have to change the FinalReport, SampleMap, sample_sheet.... and the new data we receive from illumina (IDAT files of the first batch) would need to be changed also. Therefore, it is much more complicated.
    #check if the fact you did not change dtype of week 8 test beep to float in script 1, could be a problem
        #save pheno_data after the cleaning?
    #CHECK THE QUESTIONS TO DAVID about pheno

print_text("decompress merged data", header=2)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/; \
    gunzip \
        --keep \
        --force merged_batches*.gz")


print_text("create folder to save pca", header=2)
run_bash(" \
    cd ./data/genetic_data/; \
    rm -rf quality_control/pca; \
    mkdir -p quality_control/pca")


print_text("run PCA to detect clusters of samples", header=2)
run_bash("\
    cd ./data/genetic_data/; \
    plink \
        --pca \
            'tabs' \
            'header' \
        --bfile ./plink_bed_files/merged_batches/merged_plink_files/merged_batches \
        --out ./quality_control/pca/first_pca")
        #dimensionality reduction
            #PLINK 1.9 provides two dimension reduction routines: --pca, for principal components analysis (PCA) based on the variance-standardized relationship matrix, and --mds-plot, for multidimensional scaling (MDS) based on raw Hamming distances. 
            #Top principal components are generally used as covariates in association analysis regressions to help correct for population stratification, while MDS coordinates help with visualizing genetic distances.
            #By default, --pca extracts the top 20 principal components of the variance-standardized relationship matrix; you can change the number by passing a numeric parameter. Eigenvectors are written to plink.eigenvec, and top eigenvalues are written to plink.eigenval. 
            #The 'header' modifier adds a header line to the .eigenvec file(s), and the 'tabs' modifier makes the .eigenvec file(s) tab- instead of space-delimited.
            #https://www.cog-genomics.org/plink/1.9/strat#pca



##USE OTHER TECHINES TO CHECK FOR OUTLIERS AND POP STRATIFICATION?
    #look tutorials
    #https://www.cog-genomics.org/plink/1.9/strat


#load the eigenvec file generated
import pandas as pd
first_pca=pd.read_csv(\
    "./data/genetic_data/plink_bed_files/merged_batches/merged_batches_pca/first_pca.eigenvec", \
    sep="\t", \
    header=0, \
    low_memory=False)
        #Produced by --pca. Accompanied by an .eigenval file, which contains one eigenvalue per line.
        #The .eigenvec file 
            #is, by default, a space-delimited text file with no header line and 2+V columns per sample, where V is the number of requested principal components. 
            #The --pca 'header' modifier causes a header line to be written, and the 'tabs' modifier makes this file tab-delimited.
            #The first two columns are the sample's FID/IID, and the rest are principal component weights in the same order as the .eigenval values (if the header line is present, these columns are titled 'PC1', 'PC2', ...).
        #https://www.cog-genomics.org/plink/1.9/formats#eigenvec

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

#make interactive plot of the two first PCAs, so you can zoom in
import plotly.express as px
fig = px.scatter(\
    data_frame=pca_pheno, \
    x="PC1", \
    y="PC2", \
    color="FID",
    hover_data=[\
        "IID", \
        "body_mass_diff", \
        "beep_test_diff", \
        "vo2max_diff"])
        #you can use columns of DF to add axis data, but also modify color, size, and show data per sample in a desplegable box
        #https://plotly.com/python/line-and-scatter/
        #https://plotly.com/python-api-reference/generated/plotly.express.scatter.html
#fig.show()
fig.write_html("./data/genetic_data/plink_bed_files/merged_batches/merged_batches_pca/pca_plot.html")
    #https://plotly.com/python/interactive-html-export/





##################################
###### check sex mismatches ######
##################################
print_text("check sex mismatches", header=1)




    #check sex uses allele frequencies, so you can have problems with different ancestries, yo have to prepare the data... so maybe is bettter to see the PCA for outliers and batch effects, filtering and then go to see sex once we have cleaner data
    #https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex

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












    

















#open a pool
if n_cores==None:
    pool_processors = None
        #This uses all cores available
else:
    pool_processors = int(n_cores)
import multiprocessing as mp
pool = mp.Pool(processes=pool_processors)




########################################################################
#DO NOT FORGET TO ASK DAVID QUESIONS ABOUT PHENO IN TODO.MD
########################################################################
