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

#This script will merge the plink binary files of both batches and then check sex, batch effects and then perform QCC



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



#######################
###### check sex ######
#######################
print_text("check sex", header=1)


print_text("decompress merged data", header=2)
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/merged_plink_files/; \
    gunzip \
        --keep \
        --force merged_batches*.gz")


##por aquiii













###################################
###### explore batch effects ######
###################################

print_text("create folder to save pca")
run_bash(" \
    cd ./data/genetic_data/; \
    rm -rf quality_control/pca; \
    mkdir -p quality_control/pca")


###CHANGE THE PATH lines below


#
print("\n#######################################\n#######################################")
print("run PCA to detect clusters of samples")
print("#######################################\n#######################################")
run_bash("\
    cd ./data/genetic_data/plink_bed_files/merged_batches/; \
    plink \
        --pca \
            'tabs' \
            'header' \
        --bfile ./merged_plink_files/merged_batches \
        --out ./merged_batches_pca/first_pca")
        #dimensionality reduction
            #PLINK 1.9 provides two dimension reduction routines: --pca, for principal components analysis (PCA) based on the variance-standardized relationship matrix, and --mds-plot, for multidimensional scaling (MDS) based on raw Hamming distances. 
            #Top principal components are generally used as covariates in association analysis regressions to help correct for population stratification, while MDS coordinates help with visualizing genetic distances.
            #By default, --pca extracts the top 20 principal components of the variance-standardized relationship matrix; you can change the number by passing a numeric parameter. Eigenvectors are written to plink.eigenvec, and top eigenvalues are written to plink.eigenval. 
            #The 'header' modifier adds a header line to the .eigenvec file(s), and the 'tabs' modifier makes the .eigenvec file(s) tab- instead of space-delimited.
            #https://www.cog-genomics.org/plink/1.9/strat#pca

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


##POR AQUI

    '''
REQUEST DAVID, 7800AGSO and the list of sex mismatches

Hi Diego,

Thanks for being so thorough with this!

I caught up with my post-doc today and I was able to clarify some of your queries (below in red).

If it would help to have a quick chat on Zoom about this just let me know.

Repeated IDs in illumina reports
There are 6 illumina reports with almost the same ID, just differing in the suffix:
1100JHJM_1
1100JHJM_2
1200JPJM_1
1200JPJM_2
7800AGSO_1
7800AGSO_2
I have checked that each pair of samples (e.g., 1100JHJM_1 and 1100JHJM_2) is not identical (i.e., different genotypes), which is the case. So they are indeed different samples, but they share almost the same ID.
I have found out that both 1100JHJM and 1200JPJM appear two times in the AGRF code column of the excel with phenotypic data ("combact gene DNA GWAS 23062022.xlsx"). For example, row 157 shows a 21 years old male, while row 174 shows a 18.82 years old male. Both have the same ID (1200JPJM) but different values for the phenotypic variables.
This explains why we have 4 illumina reports with these IDs changing the suffix.
I do not think we can identify which phenotypic data belong to each of these final reports, so I would remove these 4 samples.
    - This is correct. In session 844, there are two individuals with the same code (1100JHJM, both male), and in session 845, there are two individuals with the same code (1200 JPJM, both male). I agree that we will need to remove these from our analysis.
        - Regarding 1100JHJM and 1200JPJM, I will removed these samples.

In the case of 7800AGSO, there is no duplication in the excel file, but there are still two illumina reports with this ID (_1 and _2). These two reports have different genotypes, so they seem to come from different samples, but I do not know where the second ID comes from. My guess is that some error occured during the genotyping, but I am not sure.
This explains another thing I wanted to mention. The total number of ID samples in the excel file is 1463. However, the number of illumina reports I have is 216 (first batch) + 1248 (second batch), totalling 1464 samples. There is a missing sample that could be this repeated 7800AGSO report.
I do not think we can be sure which one (7800AGSO_1 or 7800AGSO_2) represents sample 7800AGSO in the excel file. Therefore, I would remove this sample (7800AGSO) and the corresponding illumina reports (7800AGSO_1 and 7800AGSO_2).
    - There is only one individual (male) with code 7800 from session 878, so there appears to have been a labelling error while preparing the DNA samples. Are you able to tell us which is the sample from the excel file that doesn’t appear to have a DNA sample? That might help us to work this out. If not, I agree that we will need to remove this sample also.
        - Regarding 7800AGSO_1 and 7800AGSO_2, I do not think I can determine which of these genotypes correspond with 7800AGSO in the excel file, so we should also remove this sample, as you said.

Another missing sample
I have an Illumina report for one sample (2397LDJA; ILGSA24-17303) that is not present in the excel file.
In contrast, the sample ID 2399LDJA is present in the excel file but there is no illumina report for this ID. 
Both IDs are the same except for the 4th digit.
Maybe they are the same sample, but it would probably be a good idea to remove both IDs.
Yes! I had exactly the same note, I think it is a mislabelling of the last digit of the number (the labelling was very hard to read on some of the blood samples). So, I think 2397LDJA; ILGSA24-17303 is 2399LDJA in the excel file.
    - Ok, so I will change the name of 2397LDJA for 2399LDJA in the genotype data and maintain this sample.

Non-matching sex
I have detected some inconsistencies between the sex reported in the excel and the one reported by illumina and based on genetic data.
There are 10 samples that are considered as sex=unknown by illumina, while they have defined sex in the excel file. My guess is this a problem from Illumina?
I have also detected some samples that are defined as Male but have two copies for SNPs in the X chromosomes that should have only 1 copy (regions of the X chromosome that have no correspondence with the Y chromosome).
There are also 9 samples for which the illumina reported sex is opposite to that shown in the phenotypic data.
I still have to run more checks about this using specific bioinformatic tools, but still this is something to consider. My guess is that this is a problem from Illumina, but just to be sure, the "Gender" data in the excel file considers biological sex, right? 
    - Yes, this is biological sex (I’ll not to check if Defence allows recruits to self-identify).
    -  If I can get a list of the samples with non-matching sex, I can double check in the raw data and also check with the defence files.



    '''


    #illumina report is 1 based? look first script.

    #same strand in both batches?
        #we are using forward right?

    #check if the fact you did not change dtype of week 8 test beep to float in script 1, could be a problem

    #check a bit more plotting script
    
    #think why you have a missing sample (2397LDJA), but when you check the number of samples between illumina and pheno_data, you only have 1 of difference. In the email to Bishop, you said that the missing sample could be the AGO... that was duplicated in illumina but not in pheno_data, so we should have 2 less samples, not 1. look the empty row between data and NAs... 

    #there are some samples that are very far away, but both from batch 1 and 2. when you zoom in the bulk of the samples, samples of both batches are evenly distributed.

    #save pheno_data after the cleaning?

    ##LOOK ALSO CLUSTERING AND OUTLIER DETECTION
        #https://www.cog-genomics.org/plink/1.9/strat







###THIIIS
    #Genome-wide association studies
        #https://www.nature.com/articles/s43586-021-00056-9
        #lapalainen 2021
        #they talk about using michinga server to impute as a standard proceedure
    #Quality Control Procedures for Genome-Wide Association Studies
        #Ritche 2011 (la versión 2022 no está disponible para mi)
        #esta gente dice de hacer un case/control analysis con los batch to detect effects, just like in here (https://www.biostars.org/p/388300/). 
        #they say that in the event of batch effect, you can use the same techinques used for population stratification.
        #lee el paper antiguo pero mira el tutorial de github de la versión mas reciente
            #https://github.com/RitchieLab/GWAS-QC
    #Identifying and mitigating batch effects in whole genome sequencing data



#you can merge snps with the same position merging the fileset with itself and using --merge-equal-pos
    #If two variants have the same position, PLINK 1.9's merge commands will always notify you. If you wish to try to merge them, use --merge-equal-pos. (This will fail if any of the same-position variant pairs do not have matching allele names.) Unplaced variants (chromosome code 0) are not considered by --merge-equal-pos.
    #Note that you are permitted to merge a fileset with itself; doing so with --merge-equal-pos can be worthwhile when working with data containing redundant loci for quality control purposes.


#FOR QUALITY CONTROL, YOU COULD USE MIGHIGGAN SERVER, WHICH ALREADY APPLIES MULTIPLE FILTERS like duplicates AND THEN ADD A FEW WITH PLINK, AUGUSTO DID THAT
#this michigan sever can be used also for imputing
    #https://www.mdpi.com/2073-4425/14/2/248
    #https://imputationserver.readthedocs.io/en/latest/pipeline/


#check if indels!
    #1:207754848-GATAA-G
    #plink has flag  --snps-only to keep snps
        #https://www.biostars.org/p/378475/




#CHECK TUTORIALs so you are sure you are follwing the correct steps for filtering...
    #In particular, we are going to use the a paper about QC by Ritchie. There is a first version 2011 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/) and a second version in 2022 (https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.603).

#CHECK THE PDF REPORTs FROM ILLUMINA FOR EACH BATCH!
    #check folder called data in the second batch?


##IMPORTANT, wen merging different batches, check same strand
    #https://www.biostars.org/p/310290/


#check there is not overlap in individuals between cbatches, no indivuals is both batches

#filter by chromosome
    #check strange chromosome numbers (non-autosomal)
    #also check that no genetic position is added


##check snp maps after merging






##remove these duplicates
#duplicated positions should be merged or removed. In our case, we are talking about 1% of the snps, so it should not be a problem.

#if there are more than 2% of duplicates, stop
if (n_duplicates_plink/n_genotypes)*100 > 2:
    raise ValueError("ERROR! WE HAVE MORE THAN 2% OF SNPS WITH DUPLICATED POSITION")

#open a folder to save filtered dataset
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    rm -rf ./04_inspect_snp_dup/01_remove_dup/; \
    mkdir -p ./04_inspect_snp_dup/01_remove_dup/")

#filter these snps
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    plink \
        --bfile ./03_merged_data/" + batch_name + "_merged_data \
        -exclude ./04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates.dupvar \
        --out  ./04_inspect_snp_dup/01_remove_dup/" + batch_name + "_merged_data_no_snp_dup \
        --make-bed")
        #-exclude a list with the SNP names as input to remove snps

#check again for duplicates
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/; \
    plink \
        --bfile ./" + batch_name + "_merged_data_no_snp_dup \
        --list-duplicate-vars suppress-first ids-only\
        --out  ./" + batch_name + "_duplicates")

#count number of duplicates
print("\n#####################\n#####################")
print("Do we have zero duplicated positions after filtering? THIS CHECK CAN BE PRINTED BEFORE THIS LINE")
print("#####################\n#####################")
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/ \
    n_lines=wc -l " + batch_name + "_duplicates.dupvar; \
    FILE=" + batch_name + "_duplicates.dupvar; \
    if [ -f $FILE ] && [ $n_lines==0 ]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
    #count the number of lines in the duplicates list
    #save the name of that file
    #if the file exists and the number of lines is zero, perfect because there no snp duplicated by position
    #else False

#remove the files we are not interested in
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup/; \
    rm " + batch_name + "_duplicates.dupvar")

#remove hh files only if they are present
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/01_remove_dup; " \
    "n_hh_files=$(ls *.hh | wc -l); \
    if [ $n_hh_files -gt 0 ]; then \
        rm ./" + batch_name + "*.hh; \
    fi")
    #count the number of files with the "hh" extension
    #if that number is greater than 0, then remove all the hh files

#compress the bed/bim/fam files
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    gzip ./03_merged_data/" + batch_name + "_merged_data.bed; \
    gzip ./03_merged_data/" + batch_name + "_merged_data.bim; \
    gzip ./03_merged_data/" + batch_name + "_merged_data.fam; \
    gzip ./04_inspect_snp_dup/01_remove_dup/" + batch_name + "_merged_data_no_snp_dup.bed; \
    gzip ./04_inspect_snp_dup/01_remove_dup/" + batch_name + "_merged_data_no_snp_dup.bim; \
    gzip ./04_inspect_snp_dup/01_remove_dup/" + batch_name + "_merged_data_no_snp_dup.fam")





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



#####################
###### SEX CASES #####

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

#to remove the dir
temp_dir.cleanup()
    #https://stackoverflow.com/questions/3223604/how-to-create-a-temporary-directory-and-get-its-path-file-name


#DO NOT FORGET TO ASK DAVID QUESIONS ABOUT PHENO IN TODO.MD



########################################################################
#DO NOT FORGET TO ASK DAVID QUESIONS ABOUT PHENO IN TODO.MD
########################################################################

#- Deep learning-based polygenic risk analysis for Alzheimer’s disease prediction
#    - https://www.nature.com/articles/s43856-023-00269-x
    #they select variants from a previos GWas and then train models in their cohorts
        #they select between 8000 and 2000 snps and use them as inputs in the different models. 
        #we have a big problem here because we cannot use a previous gwas to select snps, so we have 600K snps and onlt 1400 samples
            #maybe gwas cardiorespiratory fitness
        #this sounds like typical gwas would be better option, but we do not have validation cohort
        #maybe do like zaragoza people? select genes/snps from literature (gwas potente de cardiorespitratory fitness) and use that to make the score that in turn would be tested by splitting in training and evaluation?
            #right now this sounds like the best option.
        #BEST OPTION
            #select variants from literature (GWAS VO2 max?) and use them to create a genetic algorithm to predict response to train in terms of VO2 max.
            #indeed, they do for the chinese analyses.
                #they took seveeral GWAS across different ethinc groups, obtain 200 gwas hits, 
                #analyze the assocition of these variants with Alzhemier in the two chinese cohort
                #calculate meta-pvalue for the association controlling for PCA, obtainint 37 variants significant
                #from this associatoons extract effect size to develop PRS with these 37 variants
                #for lasso and DNN use the 37 variants
            #maybe we can do a meta from differene erupean gwas? maybe the filter again in their cohort because the etnicity was different...
            #do CV on 1000 samples and then test in 400.
            #Yes, you are losing variants that are not known, but if you get a good performance in neven seen data, then you are good.
            #remember paper of Jones, where they select genes bases on literature. Yes, it would be ideal to have another cohort like them, but we have a very controlled study... I think we have a very good study to test polygenic scores, even not having a second cohort.
                #you should read about what is the situation of poly scores in fitness research. if there are not frequent, this could be a hit.
            #then we could extract feature importance and see what variants are more important, do enrichment analyses in functions..
            #even though think about the alternative of doing gwas, you propose new genes and then what? do no have more accuacy or anything, and no another cohort
                #THINK
        #EVEN MORE BETTER?
            #selet gwas hitas from previous gwas
            #do asssociaiton analysis in our cohort and select those hits with signal in our experiment
            #predic VO2 max wihth PRS, lasso, DNN using CV and test set
            #then repeat in publich gwas
            #so we can use our great design to propose a genetic lgorithm and we test it in this and other cohorts?!!!
    #they do typcal CV with the three cohorts, but also repeat with trainin in two and evaluation in the third
        #we will probably have to do training and evalution.
        #
    #the approch is very interesting combiniing genetics and good deep learning for prediction
    #the last author is H-index 90, 
    #they used plink to detect duplicated samples between sets and remove them
    #network architecture basd on genetics!!
    #thye check model performance in subset of etnicity and sex
        #I guess they combined different ethinicities thanks to add the PCA axes, can we do that? most would be european..
        #maybe we can do it between age groups... factors influencif fitness
    #DNN works better than PRS and lasso with a small number of predictors, just 37!!!
    #INTERESTING
        #they used an R function to do botstrap and calculate CI of AUC and p-value for the differences in AUC between DNNs, lasso and PRS!!!
    #i undesrtand they got the scores from the DNN and use them to stratify indiviausl in low, high and meidum risk of AD, they also correlate these scores with multiple phenotipic variables to understand the mechiansms through which these snps are acting.
    #in patients without AD
        #they correlated variants with different phenotypes that are AD markers, so these variants are related to AD even in patients without AD, so they could be useful for screening
        #they use proteomics  to check correlation between variants and different protein levels, then go enricument analysis and foudn that inmmune stuff is implicate.
    #they also used the arcihecture of the networks, instead of using the unique score of the last layer, they used the 5 scores from the penultimate layer. These scores are not correlated between them, but correlate each one with different plasma proteins, suggesting the network is targeting different metabolic pathways implicated in the pathogenesis of AD. They also used the scores of these 5 nodes of the penultimate layer to classify patients with AD as in risk or not.

    #then they correlated variants with the score of the DNN so you can see the impact of each snp and look for functional impact of that variant.

    #graph neural network model
        #I think they present variants as nodes and connect them according to their correlation, but not sure how the CNN works here
        #it seems that this analyses tries to maximize or minimize something from one point of the network to the other, for example, flight routes minimixing fuel...

    #you could work on predicting obesity status or fat percentage in HELENA using this approach given we have a loooot of measurements in blood and good obesity markers, we could select markers by a meta-gwas (mira loos) y luego ir filtrando en la misma cohorte. Tenemos size para split y luego podriamos hablar con chiqui para validation.
        #se podria mirar incluso si hacer el proteomic, there are a lot of peolple in helena with lab experience.

    #you could also apply this to the UK biobank! they have good imaging data bout different adiposity compartiments. Visceral seems to be better indicator of complications, see nature paper of Pedro
        #Obesity and the risk of cardiometabolic diseases


    #todo esto es muy aplicado!! Data science portfolio!

#maybe these plots woth R2 and gene names of the variants and p-value could be useful?  
    #Page 46 of "aau1043-halldorsson-sm-revision1.pdf"