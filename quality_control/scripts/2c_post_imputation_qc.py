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
######## MERGE BATCHES AND ASSESS BATCH EFFECTS ########
########################################################

#This script will merge the plink binary files of both batches and then perform post imputation QC



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
print_text("checking function to print nicely", header=1)
print_text("checking function to print nicely", header=2)



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



#############################
###### merging batches ######
#############################
print_text("merging batches", header=1)


##IMPORTANT: We are merging batches after imputation, I think that merging can be done after imputation and then perform post-imputation QC in the merged dataset, as Ritchi did in this paper (https://www.frontiersin.org/articles/10.3389/fgene.2014.00370/full). See the protocol we are following to do post-imputation QC
    #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view?usp=sharing
    #https://github.com/RitchieLab/GWAS-QC

##READ the second ritchie paper to see how they merge and then what QC approaches used, the same than those of the protocol?
    #https://www.frontiersin.org/articles/10.3389/fgene.2014.00370/full

##IMPORTANT, wen merging different batches, check same strand
    #same strand in both batches?
        #we are using forward right?
    #https://www.biostars.org/p/310290/


#you have to use the correct path and files to merge, becuase this coded uses the previous non-cleaned data as input, change that

print_text("check we have 216 and 1242 samples for batch 1 and 2 (batch 2 lost 6 samples that were duplicated), respectively looking at the merged fam file of each batch and the list of samples to be merged", header=2)
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


print_text("create new folders to store files", header=2)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/; \
    rm -rf merged_batches/data_to_merge; \
    mkdir -p merged_batches/data_to_merge")
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/; \
    rm -rf merged_batches/merged_plink_files; \
    mkdir -p merged_batches/merged_plink_files")


print_text("copy the plink files of both batches", header=2)
run_bash("\
    cd ./data/genetic_data/plink_bed_files; \
    cp ./ILGSA24-17303/03_merged_data/ILGSA24-17303_merged_data* ./merged_batches/data_to_merge/; \
    cp ./ILGSA24-17873/03_merged_data/ILGSA24-17873_merged_data* ./merged_batches/data_to_merge/")


print_text("create a .txt with the name of the plink inputs for each batch", header=2)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/; \
    echo ILGSA24-17303_merged_data > list_files_to_merge.txt; \
    echo ILGSA24-17873_merged_data >> list_files_to_merge.txt")
        #for --merge-list
            #If a line contains only one name, it is assumed to be the prefix for a binary fileset
            #https://www.cog-genomics.org/plink/1.9/data#merge_list
        #">" creates a new file where the output of the command it is saved
        #">>" appends the output of the command to an existing file
            #https://unix.stackexchange.com/questions/159513/what-are-the-shells-control-and-redirection-operators
            #https://unix.stackexchange.com/questions/77277/how-to-append-multiple-lines-to-a-file


print_text("decompress the plink inputs", header=2)
run_bash("\
    gunzip --keep --force ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/ILGSA24-17303_merged_data*.gz; \
    gunzip --keep --force ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/ILGSA24-17873_merged_data*.gz")
        #-keep to keep the original compressed file and --force to overwrite if the decompressed file exists


print_text("merge the two batches using two different approaches in plink: --merge-list and --bmerge", header=2)
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


print_text("check that merging the two batches with --merge-list and --bmerge gives the same", header=2)
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


print_text("remove the second fileset created only for the check", header=2)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    rm merged_batches_2*; \
    n_files=$(ls | wc -l); \
    if [[ $n_files -eq 4 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")


print_text("check that the number of unique IDs in the merged fam file is equal to the total number of samples, i.e., no ID is duplicated", header=2)
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


#more info about merging
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


print_text("check we do NOT have a .missnp file, where variants with more than two alleles are indicated. We should not have this file")
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files; \
    count_miss_file=$(ls -l *.missnp | wc -l); \
    if [[ $count_miss_file -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #list files with extension ".missnp" and count them. This should be zero.


print_text("check that the bim files (SNP maps) of both batches and the merged file have the same SNPs, although the minor can be different", header=2)
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

#note about the differences in minor/major alleles between batches
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


print_text("see if the first SNPs that have different minor/major allele between batches are indeed different in terms of frequency and alleles according to plink --freq")
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


print_text("remove freq files", header=2)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    rm plink.frq; \
    rm plink.hh; \
    rm plink.log; \
    rm plink.nosex")


print_text("compress the merged files", header=2)
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/; \
    gzip --force ./merged_batches.bed; \
    gzip --force ./merged_batches.bim; \
    gzip --force ./merged_batches.fam")
        #--force: force overwrite of output file and compress links


print_text("remove the decompressed input files", header=2)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/data_to_merge/;\
    rm ./ILGSA24-17303_merged_data.bed; \
    rm ./ILGSA24-17303_merged_data.bim; \
    rm ./ILGSA24-17303_merged_data.fam;\
    rm ./ILGSA24-17873_merged_data.bed; \
    rm ./ILGSA24-17873_merged_data.bim; \
    rm ./ILGSA24-17873_merged_data.fam")

print_text("check again we have the exact sum of samples from both batches", header=2)
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
    if [[ $n_samples_merged -eq $n_samples_first_batch+$n_samples_second_batch ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")


print_text("check again after merging we still have 654027 snps, and this is the number of variants than in each of the batches", header=2)
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



###################################
###### explore batch effects ######
###################################
print_text("explore batch effects", header=1)


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


#por aqui




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




