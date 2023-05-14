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

#In our particular case, this should not be very important because the same SNPs are genotyped in both batches, so it is very inlikely that all samples of a batch have missing for a given SNP so that SNP is not present in one batch but it is present in the other one. 

#I going to do imputation first just in case, to avoid any problems. Also, remember that the protocol mentions the possibility to merge after imputation the batch section, so it does not seem very strange to merge batches after imputation. Indeed, I have seen a post in biostarts doing PCAs and other stuff in different batches and then in the merged dataset.
    #https://www.biostars.org/p/438079/




##por aquiii

#check the text I have in the next lines and then go directly to the protocol, follow it, except not starting with sex and hetero, see next comment.
    #https://github.com/RitchieLab/GWAS-QC

#sex (--check-sex) and hetero (--het) should be checked after PCA because accorindg to plink info, it can be problems if we have a sample with most of samples from one ancestry and then a few from another ancestry
    #https://www.cog-genomics.org/plink/1.9/basic_stats

    #R log is present in our data, so we could check prob intensity in X for a full detail sex determination, think about it.


    #IMPORTANT INFORMATION about integrity of the files
        #I get a warning when unzipping the zip of the second batch (CAGRF20093767.zip) with unzip, and I cannot access the data
        #with 7z I can access the data but I get a warning "headers error".

        #I am checking the integrity of the files using checksums:
            #https://www.tecmint.com/generate-verify-check-files-md5-checksum-linux/
            #A checksum is a digit which serves as a sum of correct digits in data, which can be used later to detect errors in the data during storage or transmission. MD5 (Message Digest 5) sums can be used as a checksum to verify files or strings in a Linux file system.
            #MD5 Sums are 128-bit character strings (numerals and letters) resulting from running the MD5 algorithm against a specific file. The MD5 algorithm is a popular hash function that generates 128-bit message digest referred to as a hash value, and when you generate one for a particular file, it is precisely unchanged on any machine no matter the number of times it is generated.
            #It is normally very difficult to find two distinct files that results in same strings. Therefore, you can use md5sum to check digital data integrity by determining that a file or ISO you downloaded is a BIT-FOR-BIT COPY OF THE REMOTE FILE or ISO.
            #you can do "md5sum file.txt" and get the 128-bit message. If you modify just a line, anything, from the file, you get a new 128-bit message.
            #you can obtain the message just piping:
                #echo "eso" | md5sum
                #echo "eso1" | md5sum #they have different message
            #This checks the content, but not the file name. So this would get the same checksum even having different names
                #echo "eso" > file_1; md5sum file_1
                #echo "eso" > file_2; md5sum file_2
            #You can redirect the hash value(s) of a file(s) into a text file and store, share them with others. For the two files above, you can issues the command below to redirect generated hash values into a text file for later use:
                #md5sum file_1 file_2 > checksum.md5
                #cat checksum.md5 # you get two lines, one per file, with checksum and file name
            #To check that the files have not been modified since you created the checksum, run "md5sum --check checksum.md5". You should be able to view the name of each file along with “OK”.
                #md5sum --check checksum.md5
                    #file_1: OK
                    #file_2: OK
            #Remember that after creating the checksum, you can not rename the files or else you get a “No such file or directory” error, when you try to verify the files with new names.
                #mv file_1 file_3
                #md5sum --check checksum.md5
                    #file_1: FAILED open or read
                    #file_2: OK
            #This is exactly what the PDFs about data integrity in both batches says
            #Our results
                #first batch: 
                    #in "CombatGenes" run "md5sum --check checksums.md5"
                    #checksums.md5: FAILED
                    #example_ILGSA24-17303_FinalReport1.txt: OK
                    #example_Sample_Map.txt: OK
                    #example_SNP_Map.txt: OK
                    #ILGSA24-17303.zip: OK
                    #md5sum: WARNING: 1 computed checksum did NOT match
                #second batch
                    #in "CombatGenes/17873" run "md5sum --check checksums.md5"
                    #205771890087.zip: OK
                    #205771890120.zip: OK
                    #205771890129.zip: OK
                    #205771890173.zip: OK
                    #205785500018.zip: OK
                    #205785500033.zip: OK
                    #205785500038.zip: OK
                    #205785500075.zip: OK
                    #205857150085.zip: OK
                    #205857150090.zip: OK
                    #205857150136.zip: OK
                    #205857150149.zip: OK
                    #205955840060.zip: OK
                    #205955840105.zip: OK
                    #205955840108.zip: OK
                    #205955840137.zip: OK
                    #205960020112.zip: OK
                    #206023350028.zip: OK
                    #206023350029.zip: OK
                    #206023350043.zip: OK
                    #206023350156.zip: OK
                    #206036460037.zip: OK
                    #206036460040.zip: OK
                    #206036460041.zip: OK
                    #206036460042.zip: OK
                    #206036460091.zip: OK
                    #206036460093.zip: OK
                    #206036460121.zip: OK
                    #206036460124.zip: OK
                    #206036460127.zip: OK
                    #206036460128.zip: OK
                    #206036460129.zip: OK
                    #206036460164.zip: OK
                    #206036460165.zip: OK
                    #206036460169.zip: OK
                    #206053690035.zip: OK
                    #206063100046.zip: OK
                    #206063100048.zip: OK
                    #206063100056.zip: OK
                    #206063100059.zip: OK
                    #206063100076.zip: OK
                    #206063100077.zip: OK
                    #206063100078.zip: OK
                    #206063100079.zip: OK
                    #206063100080.zip: OK
                    #206063100081.zip: OK
                    #206063100119.zip: OK
                    #206063100120.zip: OK
                    #206063100131.zip: OK
                    #206063100132.zip: OK
                    #206123430018.zip: OK
                    #md5sum: 206123430033.zip: No such file or directory
                    #206123430033.zip: FAILED open or read
                    #CAGRF20093767_CNMetrics.csv: OK
                    #CAGRF20093767_DNAReport.csv: OK
                    #CAGRF20093767_Reproducibility and Heritability Report.csv: OK
                    #CAGRF20093767_SampleSheet.csv: OK
                    #CAGRF20093767.zip: OK
                    #checksums.md5: FAILED
                    #GSA-24v3-0_A1_ClusterFile.egt: OK
                    #GSA-24v3-0_A2.bpm: OK
                    #PLINK_030222_0457.zip: OK
                    #md5sum: WARNING: 1 listed file could not be read
                    #md5sum: WARNING: 1 computed checksum did NOT match
                #In both batches we get some errors for some files, but the important thing is that the two zips we are using are EXACTLY THE SAME than when they were created in first place by the sequencing center. Therefore, all our problem transfering the data has not affected the data.
                    #in the first batch we get Fail for checksum file itself (checksums.md5), but the rest of files that the checksum file targets, which are the one we are interested, are OK.
                    #in the second bath we have problems again with the checksum file but also with 1 compressed file: 206123430033.zip is not present.
                        #I think these files include the image data used to do the sequencing. It seems that new illumina sequencing technology uses two color channels (red and green, which are the colors present in these zips) to determine which is present in a given position. A priori, we are not going to use this information, and if we need prob intensity, I think we can just use log R, so we should be fine.
                            #https://www.ogc.ox.ac.uk/wp-content/uploads/2017/09/techspotlight_two-channel_sbs.pdf

        #the check with sumchecks is great because we now know that the data I am using is exactly the same as created, but still we get the header warning in the zip of the second batch. As I said, that problem was not created by me nor David because the sumchecks are the same, we did not change anything, but maybe there was a problem when compressing the file by the sequencing center.


    #CHECK that log R is prob intensity again

    #do the check of integrity with 7z?

    #heritability excel, they talk about the duplicates!! these were checks about reproducibility? but why we have two times in the pheno excel file

    #look for summary PDF for both batches, I think I got the summary of the frist batch, but in summer, and it is not in the compressed file
        #preguntar david for this PDF in second batch?

    #check folder called data in the second batch?

    #check if the fact you did not change dtype of week 8 test beep to float in script 1, could be a problem
        #save pheno_data after the cleaning?

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










###################################
###### PCA ######
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












    

















#to remove the dir
temp_dir.cleanup()
    #https://stackoverflow.com/questions/3223604/how-to-create-a-temporary-directory-and-get-its-path-file-name




########################################################################
#DO NOT FORGET TO ASK DAVID QUESIONS ABOUT PHENO IN TODO.MD
########################################################################
