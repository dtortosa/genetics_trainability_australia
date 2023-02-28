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



#########################################################
######## CONVERT ILLUMINA REPORT TO PLINK FORMAT ########
#########################################################



##################
#### Starting ####
##################


######
# wd #
######

#import os
#print(os.getcwd())

#you can also use the "!" character, but only whithin ipython, it does not work if you run the script as a program
#message="hola mundo"
#!echo {message}
#!pwd
    #https://jakevdp.github.io/PythonDataScienceHandbook/01.05-ipython-and-shell-commands.html


#####################################
# General info about available data #
#####################################

#Batches
    #Within combat_genes.zip, we have two batches (ILGSA24-17303 and ILGSA24-17873), being the data separated in these two.
    #In ILGSA24-17303.zip, we have the final reports for each 216 samples, along with the sample and snp maps.
        #In the initial_stuff folder there is a zip called "ILGSA24-17303.zip" that I may have downloaded from the initial location where this data was stored in summer 2022. There are Plink files, but I am not sure this is the correct data and I cannot find the final_report files.
    #In 17873, we have the IDAT files with probs intensity from the microarrays used to genotype (first zips), the final reports+sample/snp maps (CAGRF20093767.zip) and a inputs for plink. But all of this only for 1248 individuals, not the whole cohort.
    #the phenotype csv has 1464 samples, not having phenotype data for 41 of them.
        #1248+216=1464

#I have received a Illumina report with 3 files ([link](https://www.biostars.org/p/51928/)):
    #The "FinalReport.txt" for Illumina raw genotype data generated from Genome Bead Studio for 2.5M (GSGT Version  2.0.4). This includes a header with number of SNPs, samples.... and then the data with sample index, sample names, alleles... the first row includes the column names. This is a BPM file.
        #From this file, we can obtain the a lgen file. It is just taking the first few columns of the FinalReport. Plink has an option to read the .lgen format and convert it to PED file (see below; [link](https://www.biostars.org/p/13302/))
    #A SNP map file with the physical positions of the snps.
    #A sample map file with info for each individual, like the sex, ID..

#It is usually recommended to prepare this data to create a ped file with Plink, which is a tool to process genetic data ([link](https://www.cog-genomics.org/plink/)), perform some filtering and QC analyses and then export as VCF ([link](https://www.biostars.org/p/210516/), [link](https://www.biostars.org/p/135156/), [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/)) and use packages like Hail in python.
    #https://hail.is/docs/0.2/tutorials/01-genome-wide-association-study.html

#There is an interesting alternative, gtc2vcf ([link](https://github.com/freeseek/gtc2vcf), [link](https://software.broadinstitute.org/software/gtc2vcf/)), which can directly transform Ilumina reports into VCF files from command line. We are going to use Plink, though, because it is much more widely used and there are tutorials and best practice papers using it.

#In particular, we are going to use the a paper about QC by Ritchie. There is a first version 2011 ([link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/)) and a second version in 2022 ([link](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.603)).


#######################################
# Passing arugments of python program #
#######################################

#define input arugments to be passed when running this script in bash
#we will use a bash script to run this python program two times instead of doing a function and parallelize that function. In this way, the same python script is run for each batch and if we want to run again one batch, we just need to modify the bash script, not the python program. This makes sense because we have only two batches, and the parallelization will occur inside each batch. Importantly, we can also have separated .out files, so we can look at the checks separately.
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--batch_name", type=str, default="ILGSA24-17873", help="Name of the batch used as input")
parser.add_argument("--n_cores", type=int, default=5, help="Number of cores requested")
parser.add_argument("--n_samples", type=int, default=2, help="Number of samples to be analyzed")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
batch_name = args.batch_name
n_cores = args.n_cores
n_samples = args.n_samples

#starting
print("#################################################################################################################################\n#################################################################################################################################")
print("############################################# STARTING BATCH NUMBER: " + batch_name + " ##############################################")
print("#################################################################################################################################\n#################################################################################################################################")



#########################################
# Unzip files of the batch in temp dict #
#########################################

#set the name of the zip file to be unzzipped for each batch
if batch_name=="ILGSA24-17303":
    zip_name="ILGSA24-17303"
elif batch_name=="ILGSA24-17873":
    zip_name="CAGRF20093767"

#create temporary folder to save
import tempfile
temp_dir = tempfile.TemporaryDirectory()
#print(temp_dir.name)

#to remove the dir
#temp_dir.cleanup()
    #https://stackoverflow.com/questions/3223604/how-to-create-a-temporary-directory-and-get-its-path-file-name

#read only the zip and get list of files in it
import numpy as np
import zipfile
zipdata = zipfile.ZipFile("data/genetic_data/illumina_batches/" + zip_name + ".zip")
zipinfos = zipdata.infolist()

#if we have selected only a subset of samples, then we have to make the subset of samples and then manually select the SNP and Sample maps
if (n_samples != None) and (batch_name=="ILGSA24-17303" and int(n_samples)<216) | (batch_name=="ILGSA24-17873" and int(n_samples)<1248):
    
    #select a subset of the samples
    zipinfos = zipinfos[0:n_samples+1]

    #add the paths for SNP and sample maps
    zipinfos.append(zipdata.infolist()[np.where([element.filename == zip_name + "/SNP_Map.txt" for element in zipdata.infolist()])[0][0]])
    zipinfos.append(zipdata.infolist()[np.where([element.filename == zip_name + "/Sample_Map.txt" for element in zipdata.infolist()])[0][0]])

#get the name of each zip file
names_files = [zipinfo.filename for zipinfo in zipinfos]

#iterate across files and get only final reports
print("\n#######################################\n#######################################")
print("Unzipping data: ")
print("#######################################\n#######################################")
#we are using a loop because it is a bit faster than zipdata.extractall. Parallelizing does not seems to be a good solution for unzipping and I got problems with pool. This takes a few minutes anyways.
    #parallel slower
        #https://stackoverflow.com/questions/43313666/python-parallel-processing-to-unzip-files
    #count time
        #https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution
import numpy as np
#zipinfo=zipinfos[1]
#zipinfo=zipinfos[np.where([element == zip_name + "/SNP_Map.txt" for element in names_files])[0][0]]
#zipinfo=zipinfos[np.where([element == zip_name + "/Sample_Map.txt" for element in names_files])[0][0]]
import os
for zipinfo in zipinfos:

    #fir the file name starts with the adequate zip (batch) name and FinalReport
    if (zipinfo.filename.startswith(zip_name + "/" + zip_name + "_FinalReport")) | (zipinfo.filename.startswith(zip_name + "/SNP_Map.txt")) | (zipinfo.filename.startswith(zip_name + "/Sample_Map.txt")):

        #rename by removing the name of the parent folder, so we only get the file not the parent folder
        zipinfo.filename = zipinfo.filename.split(zip_name+"/")[1]
        print(zipinfo.filename)

        #extract the file in the temp dict
        zipdata.extract(zipinfo, temp_dir.name)
            #https://stackoverflow.com/questions/44079913/renaming-the-extracted-file-from-zipfile/56362289#56362289
            #https://stackoverflow.com/questions/13765486/how-to-unzip-specific-folder-from-a-zip-with-python

        #if the selected file is a FinalReport
        if zipinfo.filename.startswith(zip_name + "_FinalReport"):

            #remove the first 10 lines of the file, which is the header
            os.system("cd " + temp_dir.name + "; tail -n +11 " + zipinfo.filename + " > tmp.txt && mv tmp.txt " + zipinfo.filename)
                #if you use tail with "+number_line" you can get all lines starting from the selected number of line
                    #tail is faster than sed
                #then save the result in a temporal file and overwrite the original file. The && will make sure that the file doesn't get overwritten when there is a problem.
                    #https://stackoverflow.com/questions/339483/how-can-i-remove-the-first-line-of-a-text-file-using-bash-sed-script
                    #https://www.baeldung.com/linux/remove-first-line-text-file

#get the number of files extract
import subprocess
#run the command to count
get_n_files = subprocess.Popen(
    args="ls " + temp_dir.name + " | wc -l", 
    shell=True, #If true, the command will be executed through the shell.
    stdout=subprocess.PIPE).stdout
    #https://unix.stackexchange.com/questions/418616/python-how-to-print-value-that-comes-from-os-system
#read and transform to int
n_files = int(get_n_files.read())
#check
print("\n#######################################\n#######################################")
print("the number of files extracted is correct?: ")
print("#######################################\n#######################################")
if batch_name == "ILGSA24-17303":
    correct_number_files = 216 + 2
elif batch_name == "ILGSA24-17873":
    correct_number_files = 1248 + 2
    #plus 2 because we also extracted sample and snp maps
print(n_files == correct_number_files)

#list files present in temp
import glob
list_files = glob.glob(temp_dir.name + "/*")

#natural sorting, 1,10, 21.. that works with numbers + strings like 1b, 5c...
from natsort import natsorted
list_files = natsorted(list_files)
    #https://github.com/SethMMorton/natsort#installation

#check
print("\n#######################################\n#######################################")
print("the number of listed files is correct?: ")
print("#######################################\n#######################################")
print(len(list_files) == correct_number_files)

#extract the names of the final reports only, which start with the zip name
list_files_samples = [file for file in list_files if file.startswith(temp_dir.name + "/" + zip_name)]



#############################################
# read illumina reports and maps with spark #
#############################################

print("\n#####################\n#####################")
print("read illumina reports and maps with spark")
print("#####################\n#####################")

#spark gives us the possibility to analyze big datasets with much less ram and relatively fast. You can work with multiple final reports, being all connected, making possible to make queries on them, but not being all loaded in memory

#the header of the final reports has been remove, so now the first row has the column names of the genotype data, but no illumina info about the report

#open a Spark session
#sparkcontext, which is used by TDI, is going to be deprecated
from pyspark.sql import SparkSession
spark = SparkSession.builder \
    .appName("working with illumina reports") \
    .master("local[" + str(n_cores) + "]") \
    .config("spark.driver.memory", "6g") \
    .getOrCreate()
    #master: local (not cluster) and selecting number cores instead of using all "*"
        #https://sparkbyexamples.com/spark/what-does-setmaster-local-mean-in-spark/#:~:text=What%20does%20setMaster(local%5B*%5D)%20in%20Spark%3F,a%20SparkSession%20or%20SparkConf%20object.&text=Here%2C%20setMaster()%20denotes%20where,URL%20for%20a%20distributed%20cluster.
    #.config() to change configuration
        #https://stackoverflow.com/questions/41886346/spark-2-1-0-session-config-settings-pyspark 
    #increase the memory for executor, if not, we get refused connection error, probably because java is killed due to out of memory problems
        #https://stackoverflow.com/questions/49995044/bizarre-exception-from-pyspark-of-local-mode-in-docker-container
    #in case "refused connection" error, you have to commenting the first two lines "/etc/hosts", the ones with 127.0.0.1 IP. Indeed sparks "prefers" not to use this IP
        #https://stackoverflow.com/questions/24881757/apache-spark-connection-refused-for-worker
    #CHECK THIS
        #select number of cores
            #https://stackoverflow.com/questions/24622108/apache-spark-the-number-of-cores-vs-the-number-of-executors

#we are going to work with data frames
    #https://spark.apache.org/docs/2.2.0/sql-programming-guide.html#spark-sql-dataframes-and-datasets-guide

#we need the same column order in all illumina reports in order to use all of them with spark. For that, extract the column names for each report using pandas. Only header
# Get first the paths for each final report
list_reports_files_full_path = glob.glob(temp_dir.name+ "/" + zip_name + "_FinalReport*.txt") 
    #I prefer using the glob module, as it does pattern matching and expansion.
        #https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
#read header of each report
import pandas as pd
list_df_reports = [
    pd.read_csv(
        path, 
        nrows=0, #no rows
        header=0, #only header
        sep='\t',
        low_memory=False) #to avoid problems with mixed types
    for path in list_reports_files_full_path]
print("\n#####################\n#####################")
print("all the columns in each report are the same than in the first report?")
print("#####################\n#####################")
print(all([(np.array(report.columns) == np.array(list_df_reports[0].columns)).all() for report in list_df_reports]))

#define the scheme, i.e., the type of the data
from pyspark.sql.types import StructType, IntegerType, StringType, DoubleType, BooleanType
schema_reports = StructType() \
    .add("Sample Index",IntegerType(), nullable=True) \
    .add("Sample ID",StringType(), nullable=True) \
    .add("Sample Name",StringType(), nullable=True) \
    .add("SNP Index",IntegerType(), nullable=True) \
    .add("SNP Name",StringType(), nullable=True) \
    .add("Chr",StringType(), nullable=True) \
    .add("Position",IntegerType(), nullable=True) \
    .add("GT Score",DoubleType(), nullable=True) \
    .add("GC Score",DoubleType(), nullable=True) \
    .add("Allele1 - AB",StringType(), nullable=True) \
    .add("Allele2 - AB",StringType(), nullable=True) \
    .add("Allele1 - Top",StringType(), nullable=True) \
    .add("Allele2 - Top",StringType(), nullable=True) \
    .add("Allele1 - Forward",StringType(), nullable=True) \
    .add("Allele2 - Forward",StringType(), nullable=True) \
    .add("Allele1 - Design",StringType(), nullable=True) \
    .add("Allele2 - Design",StringType(), nullable=True) \
    .add("Theta",DoubleType(), nullable=True) \
    .add("R",DoubleType(), nullable=True) \
    .add("X Raw",IntegerType(), nullable=True) \
    .add("Y Raw",IntegerType(), nullable=True) \
    .add("X",DoubleType(), nullable=True) \
    .add("Y",DoubleType(), nullable=True) \
    .add("B Allele Freq",DoubleType(), nullable=True) \
    .add("Log R Ratio",DoubleType(), nullable=True) \
    .add("SNP Aux",DoubleType(), nullable=True) \
    .add("SNP",StringType(), nullable=True) \
    .add("ILMN Strand",StringType(), nullable=True) \
    .add("Top Genomic Sequence",StringType(), nullable=True) \
    .add("Customer Strand",StringType(),nullable=True)
        #Chr has to be string because in case we have Y, X or other stuff
        #nullable means you can have null values
        #https://spark.apache.org/docs/3.1.3/api/python/reference/api/pyspark.sql.types.StructType.html

#create a spark DF with all the samples
df_samples = spark.read.option("delimiter", "\t").option("header", True).schema(schema_reports).csv(temp_dir.name+ "/" + zip_name + "_FinalReport*.txt")
    #Note this DF is NOT in memory, but you can make queries to go through the whole data while using few RAM
    #read using tab, using first row as header and then go to specific folder and select all files with the batch name and being final reports
    #use the previously defined schema
    #IMPORTANT, all the files have to have the same schema to work in this way. We already checked that all final reports are the same.
        #https://stackoverflow.com/questions/69350586/pyspark-read-multiple-csv-files-at-once

#check
print("\n#####################\n#####################")
print("check the column names of the reports DF is correct")
print("#####################\n#####################")
print(all(df_samples.schema.names == list_df_reports[0].columns))

#see schema and column names separately
print(df_samples.printSchema())
print(df_samples.schema.names)
    #https://stackoverflow.com/questions/39746752/how-to-get-name-of-dataframe-column-in-pyspark

#you can do operations on the columns
df_samples.select(df_samples["Position"]+1).show()

#select only columns of interest
df_samples_subset_raw = df_samples.select(["Sample Index", "Sample ID", "SNP Index", "SNP Name", "Chr", "Position", "Allele1 - Forward", "Allele2 - Forward"])
df_samples_subset_raw.show()

#add a column with the input file name
from pyspark.sql.functions import input_file_name
df_samples_subset = df_samples_subset_raw.withColumn("input_file", input_file_name())
df_samples_subset.show()

#load the maps
snp_map = pd.read_csv(temp_dir.name+ "/" + "SNP_Map.txt",
    delimiter="\t",
    header=0,
    low_memory=False)
sample_map = pd.read_csv(temp_dir.name+ "/" + "Sample_Map.txt",
    delimiter="\t",
    header=0,
    low_memory=False)

#you can also do SQL queries
#SQL vs pandas in spark, not great differences
    #https://www.confessionsofadataguy.com/dataframes-vs-sparksql-which-one-should-you-choose/
# Register the DataFrame as a SQL temporary view
df_samples_subset.createOrReplaceTempView("df_samples_subset_sql")
#do the query, asking for rows with chromosome 1
spark \
    .sql("SELECT * FROM df_samples_subset_sql WHERE Chr=1") \
    .show()

#reorder using chromosome inverse
df_samples_subset.orderBy("Chr", ascending=False).show()
    #to use in specific order
        #https://stackoverflow.com/questions/54071665/pyspark-read-multiple-csv-files-into-a-dataframe-in-order

#calculate the number of reports, i.e., samples as the number of unique input_file names
n_unique_samples = df_samples_subset.select("input_file").distinct().count()

#check
print("\n#####################\n#####################")
print("We have as many unique samples in illumina reports as samples in sample map?")
print("#####################\n#####################")
print(n_unique_samples == sample_map.shape[0])
print(len(list_files_samples) == sample_map.shape[0])

#count the number of genotypes, i.e., total number of rows
n_genotypes = df_samples_subset.count()

#check
print("\n#####################\n#####################")
print("We have as many genotypes as samples*snps in the maps")
print("#####################\n#####################")
print(n_genotypes == snp_map.shape[0]*sample_map.shape[0])

#checks
print("\n#####################\n#####################")
print("Multiple checks within each sample")
print("#####################\n#####################")

#check that the sample ID is the same than that showed in the input file name
import pyspark.sql.functions as F
#define the two columns to be compared
col_input_name_check = F.split(F.col("input_file"), "file://" + temp_dir.name + "/" + zip_name + "_").getItem(1)
    #prepare a column from splitting the input_file column and get the second item, which is FinalReportXX.txt, where XX is the index of the sample. This column was previously created.
col_index_check = F.concat(F.lit("FinalReport"), F.col("Sample Index"), F.lit(".txt"))
    #prepare a column concatenating FinalReport and .txt as literal values to the sample index, so we have the same format than in the input file name
df_samples_subset \
    .withColumn("check_1", col_index_check == col_input_name_check) \
    .select("check_1") \
    .distinct() \
    .show()
    #to the orignal data.frame, add a new column comparing the two previously defined columns, select the column with the check, get the distinct values and show.
    #Only true should be present.

#check we have the correct number of samples
n_distinct_samples = df_samples_subset \
    .select(df_samples_subset["Sample ID"]) \
    .distinct() \
    .count()
print("CHECK 2:") 
print(n_distinct_samples == sample_map.shape[0])
print(n_distinct_samples == len(list_files_samples))
    #list_files is a list with the file names of the final reports

#check we do not have duplicated SNP names
n_distinct_snp_names = df_samples_subset \
    .select(df_samples_subset["SNP Name"]) \
    .distinct() \
    .count()
print("CHECK 2b:") 
print(n_distinct_snp_names == snp_map.shape[0])

#check we do not have duplicated SNP positions
n_distinct_snp_positions = df_samples_subset \
    .select(df_samples_subset["Position"]) \
    .distinct() \
    .count()
print("CHECK 2c:")
diff_n_snps_distinct_position = snp_map.shape[0] - n_distinct_snp_positions
if diff_n_snps_distinct_position != 0:
    print("IMPORTANT: DO WE HAVE DUPLICATED POSITIONS")
    print(f"The difference between the number of SNPs in the map and the distinct positions is {diff_n_snps_distinct_position}")

#check that each sample has the same number of rows, i.e., SNPs
df_samples_subset \
    .groupBy("Sample Index") \
    .count() \
    .toDF("Sample Index", "count") \
    .withColumn("check_2d", F.col("count") == snp_map.shape[0]) \
    .select("check_2d") \
    .distinct() \
    .show()
    #group rows per sample, count the number of rows per sample, convert to DF, then create a new column checking whether count is equal to the number of snps in the map, select that column and see if only true

#do more checks
#define function to do the checks. this will work on each final report as a pandas DF
#pdf = pd.read_csv(temp_dir.name+ "/" + zip_name + "_FinalReport24.txt", delimiter="\t", header=0, low_memory=False)
def calc_checks(pdf):

    #create a pandas DF with the column names of the checks
    df = pd.DataFrame(columns=["index", "check_3", "check_4", "check_5", "check_6", "check_7", "check_8", "check_9"])

    #fill the different columns
    df["index"] = pdf["Sample Index"]
        #get the index of the sample
    df["check_3"] = len(pdf["SNP Index"].unique()) == snp_map.shape[0]
        #the number of unique SNP indexes is exactly the number of SNPs in the map?
    df["check_4"] = all(pdf["SNP Index"].values == snp_map["Index"].values)
    df["check_5"] = all(pdf["SNP Name"].values == snp_map["Name"].values)
    df["check_6"] = all(pdf["Chr"].values == snp_map["Chromosome"].values)
    df["check_7"] = all(pdf["Position"].values == snp_map["Position"].values)
        #each row (i.e., genotype) has exactly the same data regarding SNP, position... than int eh snp map?
    df["check_8"] = np.isin(pdf["Sample ID"].values[0], sample_map["ID"].values)
    df["check_9"] = np.isin(pdf["Sample Index"].values[0], sample_map["Index"].values)
        #the sample ID (all rows of a report belongs to the same sample) is included in the sample map?
    return df
#define schema of resulting DF
schema_checks = StructType() \
    .add("index",IntegerType(), nullable=True) \
    .add("check_3", BooleanType(), nullable=True) \
    .add("check_4", BooleanType(), nullable=True) \
    .add("check_5", BooleanType(), nullable=True) \
    .add("check_6", BooleanType(), nullable=True) \
    .add("check_7", BooleanType(), nullable=True) \
    .add("check_8", BooleanType(), nullable=True) \
    .add("check_9", BooleanType(), nullable=True)
#define a function to calculate conditional count only when condition is satisfied
cnt_cond = lambda cond: F.sum(F.when(cond, 1).otherwise(0))
#run the checks
checks_raw = df_samples_subset \
    .groupby("Sample Index") \
    .applyInPandas(calc_checks, schema_checks) \
    .agg(
        cnt_cond(F.col("check_3") == True).alias('check_3_cnt'), 
        cnt_cond(F.col("check_4") == True).alias('check_4_cnt'), 
        cnt_cond(F.col("check_5") == True).alias('check_5_cnt'), 
        cnt_cond(F.col("check_6") == True).alias('check_6_cnt'), 
        cnt_cond(F.col("check_7") == True).alias('check_7_cnt'), 
        cnt_cond(F.col("check_8") == True).alias('check_8_cnt'), 
        cnt_cond(F.col("check_9") == True).alias('check_9_cnt')) \
    .collect()
        #group by sample, and within each sample calculate the checks and generate a DF. This DF has a row per genotype, so all genotypes of the same sample are repeated, because of this we aggreagte
        #count only those cases that are true for each check, so we get just 1 row and 5 columns with the total number of Trues
            #https://stackoverflow.com/questions/49021972/pyspark-count-rows-on-condition

#The total number of trues in each checks should be equal to the number of genotypes, although we did operations within sample
for index, check in enumerate(checks_raw[0]):
    print("CHECK " + str(index+3) + ": " + str(check == n_genotypes))

#check that the distinct sample IDs are the same than the IDs in the sample map
distinct_sample_ids = df_samples_subset \
    .select(df_samples_subset["Sample ID"]) \
    .distinct() \
    .collect()
print("CHECK 10")
#do the check only if we have all the samples (n_samples equals to None or the total number of samples), because you cannot compare two arrays with different size and the distinct sample id will be lower than the cases in sample map if we use a subset of samples
if (n_samples == None):
    print(np.equal(np.ndarray.flatten(np.array(distinct_sample_ids)), sample_map["ID"].values))
elif (batch_name=="ILGSA24-17303" and int(n_samples)==216) | (batch_name=="ILGSA24-17873" and int(n_samples)==1248):
    print(np.equal(np.ndarray.flatten(np.array(distinct_sample_ids)), sample_map["ID"].values))

#see unique alleles
print("\n#####################\n#####################")
print("see unique alleles")
print("#####################\n#####################")
df_samples_subset \
    .select(df_samples_subset["Allele1 - Forward"]) \
    .distinct() \
    .show()
df_samples_subset \
    .select(df_samples_subset["Allele2 - Forward"]) \
    .distinct() \
    .show()


######################
# Plink Installation #
######################

#I have downloaded Plink ([link](https://www.cog-genomics.org/plink/)) and copied the executable ("plink") to `bin`, so we can use Plink from any place just typing `plink`. We are using Plink version 1.9 (see singularity recipe for the version used).
#os.system("plink --version")

#Note that there is Plink 1.9 ([link](https://www.cog-genomics.org/plink/1.9/)) and Plink 2.0 ([link](https://www.cog-genomics.org/plink/2.0/)), these are not connected but different programs. 
    #This [threat](https://www.biostars.org/p/299855/#:~:text=The%20main%20difference%20is%20that,for%20a%20while%20to%20come.) of biostars explains the differences:
        # "The main difference is that plink 1.9 is essentially finished, while plink 2.0 is an alpha-stage program which will have significant unfinished components for a while to come. As a consequence, current development priorities for plink 2.0 are centered around things which are impossible to do with plink 1.9, such as handling multiallelic/phased variants and dosage data and reliably tracking REF/ALT alleles; while things that plink 1.9 already handles perfectly well, such as working with .ped/.map file pairs, have been deliberately deprioritized for now.
        # So, **you should stick to 1.9 as long as it's good enough for the jobs you need to perform. But once you need to do something outside 1.9's scope, you're likely to find that 2.0 already has the additional feature you need** (or it'll be added quickly after you ask for it)"

#see the dev version of this script for details about the final report like allele strands, plink usage, etc....


####################################################
#### Conversion of final report to Plink inputs ####
####################################################

# We are going to convert the illumina report (BPM file) to lgen file that can be used as input in Plink. This seems to be trivial ([link](https://www.biostars.org/p/51928/)), but we are going to use a tutorial just in case ([link](https://www.youtube.com/watch?v=5_pNby7l2dE&t=1s&ab_channel=GenomicsBootCamp), [link](https://pastebin.com/pzgw7JVp)).

# The final report has 1 line for the SNP of one individual. We have separated final reports for each individual. This is close to the lgen format of plink.

# Therefore, our goal is to create lgen, also the fam and map files required to load the lgen in Plink. Some information might be missing in the final report, so you need to replace them.


#############
# lgen file #
#############

print("\n#####################\n#####################")
print("Preparing lgen file")
print("#####################\n#####################")

# **lgen file [plink info](https://www.cog-genomics.org/plink/1.9/formats#lgen)**
    # A text file with no header line, and one line per genotype call (or just not-homozygous-major calls if 'lgen-ref' was invoked) usually with the following five fields:
        #Family ID
        #Within-family ID
        #Variant identifier
        #Allele call 1 ('0' for missing)
        #Allele call 2

#subset all_reports selecting some columns. The non-selected columns were used for checks
lgen_file = df_samples_subset.select(["Sample ID", "SNP Name", "Allele1 - Forward", "Allele2 - Forward"])
lgen_file.show()

# Change to "0" those genotype calls with "--" to match the format of Plink
#allele 1
lgen_file = lgen_file \
    .withColumn("Allele1 - Forward_final", \
        F.when( \
            (lgen_file["Allele1 - Forward"] == "-") | (lgen_file["Allele1 - Forward"] == "--"), \
            "0") \
            .otherwise(lgen_file["Allele1 - Forward"])) \
    .withColumn("Allele2 - Forward_final", \
        F.when( \
            (lgen_file["Allele2 - Forward"] == "-") | (lgen_file["Allele2 - Forward"] == "--"), \
            "0") \
            .otherwise(lgen_file["Allele2 - Forward"])) \
    .drop("Allele1 - Forward") \
    .drop("Allele2 - Forward")
        #https://stackoverflow.com/questions/44773758/how-to-conditionally-replace-value-in-a-column-based-on-evaluation-of-expression
        #https://stackoverflow.com/questions/29600673/how-to-delete-columns-in-pyspark-dataframe
lgen_file.show()

#check
print("\n#####################\n#####################")
print("Check that all SNPs do NOT have '-' or '--' for allele 1 and 2")
print("#####################\n#####################")
(lgen_file \
    .filter( \
        (lgen_file["Allele1 - Forward_final"].isin(["-", "--"])) | (lgen_file["Allele2 - Forward_final"].isin(["-", "--"]))) \
    .count()) == 0
        #https://stackoverflow.com/questions/63330350/pyspark-dataframe-filter-column-contains-multiple-value

#Add additional columns that are required for lgen files
lgen_file = lgen_file \
    .withColumn("FID", F.lit("combat_" + batch_name))

# Reorder the columns
lgen_file = lgen_file \
    .select(["FID", "Sample ID", "SNP Name", "Allele1 - Forward_final", "Allele2 - Forward_final"])
        #https://stackoverflow.com/questions/42912156/python-pyspark-data-frame-rearrange-columns

#change the column of sample ID
lgen_file = lgen_file \
    .withColumnRenamed("Sample ID", "sample_id")
        #this column name will be used when writing the data per sample, so it is better to avoid the spaces

#duplicate sample_id column because when you partition by that column during writing, it gets removed
lgen_file = lgen_file \
    .withColumn("sample", F.col("sample_id"))
    #https://stackoverflow.com/questions/44575911/spark-how-to-prevent-dataframewriter-from-dropping-the-partitioning-columns-on

#look
lgen_file.show()

#save the data

#we can data in binary plink format split across many files, this can be then merged, so we can just save the results of spark in different csv file
    #book plink, tutorial chirstofer
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3

#prepare folders to save the lgen files
os.system("rm -rf data/genetic_data/plink_inputs/" + batch_name + "; mkdir -p data/genetic_data/plink_inputs/" + batch_name)

#you can save each final report separately by ID
lgen_file \
    .write \
        .partitionBy("sample") \
        .option("header", False) \
        .option("delimiter", "\t") \
        .option("compression", "gzip") \
        .csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/")
        #https://stackoverflow.com/questions/47780397/saving-dataframe-records-in-a-tab-delimited-file
        #https://stackoverflow.com/questions/40163996/how-to-save-a-dataframe-as-compressed-gzipped-csv
        #https://sparkbyexamples.com/pyspark/pyspark-partitionby-example/

#NOTE about cores and partitions
    #Note that the number of files generated depends on the size of each file and the number of cores. So if you have 5 cores, i.e., local[5], at least you will get 5 partions of the data.
    #you can check the number of partitions with "lgen_file.rdd.getNumPartitions()".
    #I have been doing analyses with the partitions that my local/clusterlet me, i.e., 5 in my laptop and tens in the HPC.
    #Then, when writing, I specify that I want the data of each sample separated in different folders with "partitionBy". I do not care if several files are saved in the same individual if all belong to the same Sample ID, later i will merge the data
        #https://sparkbyexamples.com/pyspark/pyspark-partitionby-example/
        #https://stackoverflow.com/questions/48143159/spark-write-to-disk-with-n-files-less-than-n-partitions

#or a single file you can use PARTITION OR COALESCEN
#lgen_file \
#    .coalesce(1) \
#    .write \
#        .option("header", True) \
#        .option("delimiter", "\t") \
#        .option("compression", "gzip") \
#        .csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/")
    #https://sparkbyexamples.com/spark/spark-write-dataframe-single-csv-file/
    #https://sparkbyexamples.com/spark/spark-repartition-vs-coalesce/#dataframe-%20coalesce

#remove the duplicated sample column
lgen_file = lgen_file \
    .drop("sample")
lgen_file.show()


############
# fam file #
############

print("\n#####################\n#####################")
print("Preparing fam file")
print("#####################\n#####################")

# **Fam file ([PLINK sample information file](https://www.cog-genomics.org/plink/1.9/formats#fam))**
    # A text file with no header line, and one line per sample with the following six fields:
        # - Family ID ('FID')
        # - Within-family ID ('IID'; cannot be '0')
        # - Within-family ID of father ('0' if father isn't in dataset)
        # - Within-family ID of mother ('0' if mother isn't in dataset)
        # - Sex code ('1' = male, '2' = female, '0' = unknown)
        # - Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

##check sex variable
#before creating the fam file, we need to see what is going on with the sex, because there are differences between phenotype dataset and the sample map of the first batch

print("\n#####################\n#####################")
print("Check sex in sample map and phenotype data")
print("#####################\n#####################")

#load pheno data, this include reported sex and VO2 max data. I have checked that the data is the same directly reading from excel than converting to csv
pheno_data = pd.read_excel(
    "data/pheno_data/combact gene DNA GWAS 23062022.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data)

#there are several phenotypes, see dev version for details about possible errors.
#here we are only going to modify an entry that is clearly wrong and do some checks with the sample map. Also, we are going to use the sex indicated in pheno_data because there are some samples with differences between illumina estimated sex and the one showed in the pheno_data

#beep test 8 has an error, "o" letter instead "0" number in one sample
print(pheno_data.loc[pheno_data["Week 8 beep test"] == "11.1O", "Week 8 beep test"])

#change 11.1O for 11.10
pheno_data.loc[pheno_data["Week 8 beep test"] == "11.1O", "Week 8 beep test"] = 11.10

#check
print("\n#####################\n#####################")
print("All samples with genetic data are included in pheno_data?")
print("#####################\n#####################")
print(sum(sample_map["ID"].isin(pheno_data["AGRF code"])) == sample_map.shape[0])

print("\n#####################\n#####################")
print("How many samples with genetic data are included in pheno_data?")
print("#####################\n#####################")
print(f'{sum(sample_map["ID"].isin(pheno_data["AGRF code"]))} out {sample_map.shape[0]}')

#get sample IDs and gender from sample map, then modify using pheno_data
fam_file = sample_map.loc[:, ["ID", "Gender"]]

#remove those samples without phenotyipic data
#IMPORTANT HERE, we lose samples with genetic data because they do not have pheno data and hence, no sex data from the study can be obtained, only illumina
#fam_file_raw = fam_file_raw.loc[fam_file_raw["ID"].isin(pheno_data["AGRF code"]), :]
    #IN CASE YOU NEED TO REMOVE

# Codify the sex variable following plink notation and considering the sex in pheno_data
fam_file.loc[
    (fam_file["ID"].isin(pheno_data.loc[pheno_data["Gender"].isna(), "AGRF code"])) | 
    (~fam_file["ID"].isin(pheno_data["AGRF code"])), 
    "Gender"] = "0"
    #set as 0 the sex of those samples whose
        #ID is associated to NA gender in pheno_data
        #OR
        #ID is NOT included in pheno_data
fam_file.loc[
    fam_file["ID"].isin(pheno_data.loc[pheno_data["Gender"] == "M", "AGRF code"]), 
    "Gender"] = "1"
    #set as 1 the sex of those samples whose ID is associated with Gender=M in pheno_data
fam_file.loc[
    fam_file["ID"].isin(pheno_data.loc[pheno_data["Gender"] == "F", "AGRF code"]), 
    "Gender"] = "2"

#check
print("\n#####################\n#####################")
print("The new sex variable is correctly coded?")
print("#####################\n#####################")
print((len(fam_file["Gender"].unique()) == 3) & (np.isin(fam_file["Gender"].unique(), test_elements=["1", "2", "0"])).all())
    #only three possible sexes, being 1, 2 and 0.

# Add the family variables and the phenotype (not added for now)
fam_file["FID"] = "combat_" + batch_name #ID for the whole study
fam_file["IID_father"] = "0" #'0' if father isn't in dataset
fam_file["IID_mother"] = "0" 
fam_file["phenotype"] = -9 #this is no data for phenotype

# Reorder
fam_file = fam_file[["FID", "ID", "IID_father", "IID_mother", "Gender", "phenotype"]]

print("\n#####################\n#####################")
print("See fam file:")
print("#####################\n#####################")
print(fam_file)

#save without header:
fam_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + ".fam.gz",
    sep="\t",
    header=None,
    compression='gzip',
    index=False)


############
# map file #
############

print("\n#####################\n#####################")
print("Preparing map file")
print("#####################\n#####################")

# **Map file [Plink info](https://www.cog-genomics.org/plink/1.9/formats#map)**

# A text file with no header line, and one line per variant with the following 3-4 fields:
    # - Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
    # - Variant identifier
    # - Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
    # - Base-pair coordinate

# We can get this information from the SNP map.
# Make a copy with copy(), so modifications to the data or indices of the copy will not be reflected in the original object
map_file = snp_map.copy()
print(map_file)

#set genetic position as dummy
map_file["centimorgans"] = 0 #dummy value of zero

#select columns
map_file = map_file.loc[:, ["Chromosome", "Name", "centimorgans", "Position"]]
print(map_file)

# Save
map_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + ".map.gz",
    sep="\t",
    header=None,
    compression='gzip',
    index=False)


#################################
# check with lgen and map files #
#################################

print("\n#####################\n#####################")
print("checks with lgen, map and fam files")
print("#####################\n#####################")
print("SNPs are the same in lgen and map files?")
unique_snp_names = lgen_file \
    .select("SNP Name") \
    .distinct() \
    .collect()
np.array_equal(np.array(unique_snp_names), map_file["Name"].values)
print("samples are the same in lgen and fam files?")
unique_sample_ids = lgen_file \
    .select("sample_id") \
    .distinct() \
    .collect()
np.array_equal(np.array(unique_sample_ids), fam_file["ID"].values)



##########################################
#### create plink binaries per sample ####
##########################################

print("\n#####################\n#####################")
print("create plink binaries per sample")
print("#####################\n#####################")

#list the of folders with the lgen data per sample
import glob
lgen_full_paths = glob.glob("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/sample=*", recursive=False)
#check
print("\n#####################\n#####################")
print("We have the correct number of lgen folders? i.e., one per sample?")
print("#####################\n#####################")
print(len(lgen_full_paths) == fam_file.shape[0])

#get the extensions of the files in these folders
list_files = glob.glob("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/**/*.csv.gz", recursive=True)
list_files = [path.split("/")[-1] for path in list_files]
list_extensions = ["-".join(path.split("-")[2:]) for path in list_files]
    #load paths of all compressed partitions of the lgen files for all samples
    #split each one by "/" and select the last one to get the file name
    #split each one by "-" and select all except the first two, then join with "-" so we get only the extension avoiding the partition name
    
#then check all files have the same extension. They should even if the belong to different samples
print("\n#####################\n#####################")
print("lgen files have the correct extension?")
print("#####################\n#####################")
if all([file == list_extensions[0] for file in list_extensions]):

    #perfect
    print("TRUE")

    #extract the extension with file type
    lgen_extension = list_extensions[0]
    lgen_extension = lgen_extension.split(".csv.gz")[0]
else:
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTENSION OF THE LGEN FILES")

#natural sorting, 1,10, 21.. that works with numbers + strings like 1b, 5c...
from natsort import natsorted
lgen_full_paths = natsorted(lgen_full_paths)
    #https://github.com/SethMMorton/natsort#installation

#create a folder to save the bed files for the selected batch
import os
os.system(
    "rm -rf data/genetic_data/plink_bed_files/" + batch_name + "; \
    mkdir -p data/genetic_data/plink_bed_files/" + batch_name)
    #mkdir
        #p flag: A flag which enables the command to create parent directories as necessary. If the directories exist, no error is specified

#define function to calculate the bed file for each lgen file
#lgen_full_path = lgen_full_paths[1]
def plink_inputs_prep(lgen_full_path):

    sample_id_file = lgen_full_path.split("/")[-1]

    #extract the name of the lgen file and path
    #lgen_file_name = lgen_full_path.split("/")[-1]
    #lgen_file_name_no_ext = lgen_file_name.split(".csv.gz")[0]
    #sample_number = lgen_file_name.split("-")[1]
    #lgen_path = "/".join(lgen_full_path.split("/")[0:-1])

    #create a folder to temporary save the plink inputs
    #os.system(
    #    "cd " + lgen_full_path + "; rm -rf inputs_sample_" + sample_number + "; \
    #    mkdir -p inputs_sample_" + sample_number)

    #decompress the different partitions of the lgen file but remove just in case before that
    os.system(
        "cd " + lgen_full_path + "; \
        rm *" + lgen_extension + ".csv; \
        gunzip -k *" + lgen_extension + ".csv.gz")



    files_selected_sample = glob.glob(lgen_full_path + "/*.csv", recursive=False)
    files_selected_sample = natsorted(files_selected_sample)
        #to have the same order of snps than in the snp map 
    list_dfs = [
        pd.read_csv(
            path, 
            header=None,
            sep='\t',
            low_memory=False) #to avoid problems with mixed types
        for path in files_selected_sample]
    len(list_dfs) == len(files_selected_sample)

    # Concatenate, because we have one row per SNP and sample, so we can just concateneate DFs.
    lgen_file_selected_sample = pd.concat(
        objs=list_dfs, #list DFs
        axis=0, #concatenate along rows
        ignore_index=True) #clear the existing index and reset it
    
    # We should have as many rows as the total sum of rows across the list of DFs
    print(f' {lgen_file_selected_sample.shape[0] == sum([element.shape[0] for element in list_dfs])}')


    selected_sample_id = np.unique(lgen_file_selected_sample[1])

    


    #check whether 
        #the number of unique Sample IDs is equal to 1 AND
        #the ID in the file is the same than the ID in the name of the folder
        #the number of unique SNP names is not equal than the number of SNPs in the snp map OR
        #the SNP names are not exactly the same than in the snp map
    sample_checks = \
        (len(selected_sample_id) == 1) and \
        (selected_sample_id == sample_id_file.split("=")[1]) and \
        (len(np.unique(lgen_file_selected_sample[2])) == snp_map.shape[0]) and \
        all(lgen_file_selected_sample[2] == snp_map["Name"])
    if sample_checks: #if all True
        print("WE PASS SAMPLE CHECKS FOR LGEN FILE OF SAMPLE " + selected_sample_id)
    else:
        raise ValueError("ERROR! FALSE! WE HAVE IMPORTANT ERRORS IN THE LGEN FILE OF SAMPLE")

    #save 
    lgen_file_selected_sample.to_csv(lgen_full_path + "/" + sample_id_file + ".lgen", 
        sep='\t',
        header=False, 
        index=False)

    #get the row of the fam file corresponding to the selected sample
    #if we use the whole fam file, we would get a ped file with as many rows as sample, but as we are only using the final report of one sample, all the rest of rows would be zero.
    #we need a a fma file with one row, so we can generate a ped file with one row for the sample selected.
    fam_file_to_subset = pd.read_csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + ".fam.gz", 
        delimiter="\t", 
        header=None, 
        low_memory=False)
    select_fam = fam_file_to_subset.loc[fam_file_to_subset[1] == selected_sample_id[0], :]
    
    #additional check about the number of samples in the lgen file
    if select_fam.shape[0] != 1:
        raise ValueError("ERROR! FALSE! WE DO NOT HAVE 1 SAMPLE IN THE FAM FILE")

    #save that as the fam file for this sample
    select_fam.to_csv(lgen_full_path + "/" + sample_id_file + ".fam", 
        sep='\t',
        header=False, 
        index=False)

    #remove the lgen file
    del(lgen_file_selected_sample)

    #copy the whole map file, in this case, we have the same snps, i.e., the same map in all samples. Indeed, illumina gives one map for the whole batch. Save the file in the folder of the sample
    os.system(
        "cd data/genetic_data/plink_inputs/" + batch_name + "; \
        gunzip -c " + batch_name + ".map.gz > " + batch_name + "_lgen_files/" + sample_id_file + "/" + sample_id_file + ".map")
        #c flag: writes the output stream to stdout and then you can use > to redirect to another folder. This will leave the compressed file untouched.
            #https://superuser.com/questions/45650/how-do-you-gunzip-a-file-and-keep-the-gz-file

    #create a folder to save the bed file of the sample
    os.system(
        "rm -rf data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file + "; \
        mkdir -p data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file)

    #create ped files file using the lgen files
    os.system(
        "cd ./data/genetic_data/; \
        plink \
            --lfile ./plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/" + sample_id_file + " \
            --out ./plink_bed_files/" + batch_name + "/01_bed_per_sample/sample_" + sample_number + "/" + batch_name + "_sample_" + sample_number + " \
            --recode")
        #go to the folder with plink inputs
        #--lfile for loading the lgen file, which should be accompanied by a .fam and .map files having the same name except the extension.
            #https://www.cog-genomics.org/plink/1.9/formats
        #--recode creates a new text fileset, after applying sample/variant filters and other operations. By default, the fileset includes a .ped and a .map file, readable with --file.
            #https://www.cog-genomics.org/plink/1.9/data#recode
        #--out for the output name

    #the result is ped and map files 
    #.ped
        #Original standard text format for sample pedigree information and genotype calls. Normally must be accompanied by a .map file; Haploview requires an accompanying .info file instead. Loaded with --file, and produced by --recode.
        #Contains no header line, and one line per sample with 2V+6 fields (columns) where V is the number of variants (as many columns as variants). The first six fields are the same as those in a .fam file. The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.
        #If all alleles are single-character, PLINK 1.9 will correctly parse the more compact "compound genotype" variant of this format, where each genotype call is represented as a single two-character string. This does not require the use of an additional loading flag. You can produce such a file with "--recode compound-genotypes".    
    #map.
        #Variant information file accompanying a .ped text pedigree + genotype table. Also generated by "--recode rlist".
        #A text file with no header line, and one line per variant with the following 3-4 fields:
            #Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
            #Variant identifier
            #Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
            #Base-pair coordinate
        #All lines must have the same number of columns (so either no lines contain the morgans/centimorgans column, or all of them do)

    #convert the ped file to bed
    os.system(
        "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/sample_" + sample_number + "; \
        plink \
            --file ./" + batch_name + "_sample_" + sample_number + " \
            --out ./" + batch_name + "_sample_" + sample_number + " \
            --make-bed")
        #--file let you load ped files
            #https://www.cog-genomics.org/plink/1.9/formats#ped
        #--make-bed creates bed file
            #https://www.cog-genomics.org/plink/1.9/formats#bed
            #https://zzz.bwh.harvard.edu/plink/binary.shtml
        #you could just do these two step in one loading the lgen file and using make-bed
            #https://icb.med.cornell.edu/wiki/index.php/Plink/howto
        #bed file is NOT the bed format of USCS!! It is a plink format more compacted and it seems is faster to work with. It is in hexadecimal, meaning that the genotypes are not stored as number but in just two letters. For example, in 0xdc, the two last letters are storing the genotype of several individuals. dc in binary is 11011100. This means that the first sample is is 00, i.e., homozygous for the first allele of this snp in the bim file, the second sample is mozygous for the second allele of this snp in the bim file, and so on... If you more samples, another hexa code is added until all samples are coded and then it moves to the next snp...
            #https://www.cog-genomics.org/plink/1.9/formats#bed
            #https://www.biostars.org/p/113166/
            #https://coolconversion.com/math/binary-octal-hexa-decimal/Convert_binary_number_11011100_in_decimal_

    #remove the files we are not interested in
    os.system(
        "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/sample_" + sample_number + "; \
        rm ./" + batch_name + "_sample_" + sample_number + ".ped; \
        rm ./" + batch_name + "_sample_" + sample_number + ".map")
        #the hh file says there is a genotype that is diplodi but should be haplo, like Y in male
        #it is created when doing other operations like calculate the frequency

    #remove hh files only if exists
    os.system(
        "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/sample_" + sample_number + "; " \
        "FILE=./" + batch_name + "_sample_" + sample_number + ".hh; \
        if test -f $FILE; then \
            rm ./" + batch_name + "_sample_" + sample_number + ".hh; \
        fi")
        #https://linuxize.com/post/bash-check-if-file-exists/

    #remove also the folder with the inputs for the sample
    os.system("rm -rf data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/inputs_sample_" + sample_number)

#open a pool
import multiprocessing as mp
pool = mp.Pool(processes=n_cores)

#apply the function across the pool
results_pararllelize = [pool.apply_async(plink_inputs_prep, args=(lgen_file,)) for lgen_file in lgen_full_paths]

#close the pool
pool.close()

#wait for the completion of all scheduled jobs
pool.join() #without this, the script is finished without that all samples have been finished.
    #https://stackoverflow.com/questions/44896547/python-apply-async-doesnt-execute-function

#NOTE: 
    #when we have one sample per bed file, the bim file, which is a map file with also the alele options, will have partial information. for example, if the unique sample is homozigous for a SNP, e.g., AA, the map in that position will have A and 0, because there are no more alleles.
        #no problem, when a sample with the other allele is included, the bim file adds the other allele. I have checked several cases before and after mergin all bed/bim files into one file.



#######################################
#### merge binaries of all samples ####
#######################################

print("\n#####################\n#####################")
print("merge binaries of all samples")
print("#####################\n#####################")

#create a folder to save the binary files of all samples and then merge
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "; \
    rm -rf ./02_data_to_merge; mkdir -p ./02_data_to_merge")

#copy the bed/bim/fam files of all samples to the same folder
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "; \
    find ./01_bed_per_sample/  | \
        egrep '" + batch_name + ".*(bed|bim|fam)' | \
        while read file_path; \
            do cp $file_path ./02_data_to_merge/; \
        done")
    #cd to the folder with bed files for the selected batch
    #find in a subfolder with bed files per sample
        #any file starting with batch name and ending with bed/bim/fam extensions
        #read the paths of these files
            #copy the file to the folder for merging
        #https://stackoverflow.com/questions/15617016/copy-all-files-with-a-certain-extension-from-all-subdirectories
        #https://unix.stackexchange.com/questions/189210/find-command-multiple-conditions

#get all bed files paths and names
bed_full_paths = glob.glob("data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/" + batch_name + "*.bed", recursive=False)
    #https://www.geeksforgeeks.org/how-to-use-glob-function-to-find-files-recursively-in-python/
bed_file_names = [path.split("/")[-1].split(".")[0] for path in bed_full_paths]

#save the names in a txt file
with open(r"./data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/list_to_merge.txt", "w") as fp:
    fp.write("\n".join(bed_file_names))
        #each name in a different line so we have to add "\n" to the name
        #https://pynative.com/python-write-list-to-file/

#create folder to save merged data
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "; \
    rm -rf 03_merged_data; \
    mkdir -p 03_merged_data")

#merge with plink
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/; \
    plink \
        --merge-list ./list_to_merge.txt \
        --out ../03_merged_data/" + batch_name + "_merged_data")
    #run plink on 02_data_to_merge where all the bed files are located
    #then save the output in a folder outside this folder ("../")    
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3



#############################
#### look for duplicates ####
#############################

print("\n#####################\n#####################")
print("look for duplicates")
print("#####################\n#####################")

#create general folder to do operations related to duplicates
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    rm -rf ./04_inspect_snp_dup; \
    mkdir -p ./04_inspect_snp_dup")

#there are duplicated positions in the original snp map, because of this we see a warning when merging. This is common in illumina data
    #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
    #https://www.biostars.org/p/281276/

#create a folder to list duplicates
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/; \
    rm -rf ./00_list_dup; \
    mkdir -p ./00_list_dup")

#list duplicates
os.system(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    plink \
        --bfile ./03_merged_data/" + batch_name + "_merged_data \
        --list-duplicate-vars suppress-first ids-only\
        --out  ./04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates")
        #list-duplicate-vars to list duplicates by POSITION
            #suppress-first prevents the first variant in each group from being reported (since, if you're removing duplicates, you probably want to keep one member of each group).
            #ids-only modifier removes the header and the position/allele columns, so the generated list can be used as input for --exclude
                #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
        #in plink 2 you can also look for duplicates in ID, but we already check SNP names duplicates with spark and it is seems is ok

#load the duplicate list
duplicate_cases = pd.read_csv(
    "./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates.dupvar", 
    sep="\t", 
    header=None,
    low_memory=False)

#check
print("\n#####################\n#####################")
print("the number of duplicates according with plink matches my previous estimation with spark?")
print("#####################\n#####################")
print(duplicate_cases.shape[0] == diff_n_snps_distinct_position)


##remove these duplicates
#duplicated positions should be merged or removed. In our case, we are talking about 1% of the snps, so it should not be a problem.

#if there are more than 2% of duplicates, stop
if (duplicate_cases.shape[0]/n_genotypes)*100 > 2:
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



######################
#### finishing up ####
######################

#remove the temp dir
temp_dir.cleanup()

#stop spark env
spark.stop()


#the second batch fails because it produces too much splits, it should be 2 only but we see 4
    #SOLVED, but do checks about the samples separations

#search for duplicates snp positions. there is a link to convert gsa names to rsi!
    #https://www.biostars.org/post/search/?query=duplicated+snp+names+illumina
    #https://www.biostars.org/p/277737/


#CHECK THE DIFFERENT FORMATS GENErated in ped/bim...
