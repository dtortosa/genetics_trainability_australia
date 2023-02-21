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
parser.add_argument("--batch_name", help="Name of the batch used as input")
parser.add_argument("--n_cores", type=int, default=5, help="Number of cores requested")
parser.add_argument("--n_samples", type=int, default=None, help="Number of samples to be analyzed")
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
#batch_name = args.batch_name
#n_cores = args.n_cores
#n_samples = args.n_samples

#for debugging
#batch_name = "ILGSA24-17873"
batch_name = "ILGSA24-17303"
n_cores = 7
n_samples = 10

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
zipinfos = zipdata.infolist()[0:n_samples+1]

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

#if we set n_samples, then reduce the number of samples using this arugment
if n_samples != None:
    list_files_samples = list_files_samples[0:n_samples]


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
list_reports_files_full_path = glob.glob(temp_dir.name+ "/" + batch_name + "_FinalReport*.txt") 
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
df_samples = spark.read.option("delimiter", "\t").option("header", True).schema(schema_reports).csv(temp_dir.name+ "/" + batch_name + "_FinalReport*.txt")
    #Note this DF is NOT in memory, but you can make queries to go through the whole data while using few RAM
    #read using tab, using first row as header and then go to specific folder and select all files with the batch name and being final reports
    #use the previously defined schema
    #IMPORTANT, all the files have to have the same schema to work in this way. We already checked that all final reports are the same.
        #https://stackoverflow.com/questions/69350586/pyspark-read-multiple-csv-files-at-once

print("\n#####################\n#####################")
print("inspect the data and do some operations to see what spark can do")
print("#####################\n#####################")


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
#create schemas
schema_snp_map = StructType() \
    .add("Index",IntegerType(), nullable=True) \
    .add("Name",StringType(), nullable=True) \
    .add("Chromosome",StringType(), nullable=True) \
    .add("Position",IntegerType(), nullable=True) \
    .add("GenTrain Score",DoubleType(), nullable=True) \
    .add("SNP",StringType(), nullable=True) \
    .add("ILMN Strand",StringType(), nullable=True) \
    .add("Customer Strand",StringType(),nullable=True) \
    .add("NormID",StringType(),nullable=True)
        #Chr has to be string because in case we have Y, X or other stuff
        #nullable means you can have null values
        #https://spark.apache.org/docs/3.1.3/api/python/reference/api/pyspark.sql.types.StructType.html
schema_sample_map = StructType() \
    .add("Index",IntegerType(), nullable=True) \
    .add("Name",StringType(), nullable=True) \
    .add("ID",StringType(), nullable=True) \
    .add("Gender",StringType(), nullable=True) \
    .add("Plate",StringType(), nullable=True) \
    .add("Well",StringType(), nullable=True) \
    .add("Group",StringType(), nullable=True) \
    .add("Parent1",StringType(), nullable=True) \
    .add("Parent2",StringType(), nullable=True) \
    .add("Replicate",StringType(), nullable=True) \
    .add("SentrixPosition",StringType(), nullable=True) \
        #nullable means you can have null values
        #https://spark.apache.org/docs/3.1.3/api/python/reference/api/pyspark.sql.types.StructType.html
#read
snp_map = spark.read.option("delimiter", "\t").option("header", True).schema(schema_snp_map).csv(temp_dir.name+ "/" + "SNP_Map.txt")
sample_map = spark.read.option("delimiter", "\t").option("header", True).schema(schema_sample_map).csv(temp_dir.name+ "/" + "Sample_Map.txt")
print(snp_map.printSchema())
print(sample_map.printSchema())

#check the column names of the three files are correct
snp_map_pandas = pd.read_csv(temp_dir.name+ "/" + "SNP_Map.txt",
    delimiter="\t",
    header=0,
    low_memory=False)
sample_map_pandas = pd.read_csv(temp_dir.name+ "/" + "Sample_Map.txt",
    delimiter="\t",
    header=0,
    low_memory=False)
print("\n#####################\n#####################")
print("check the column names of the three files are correct")
print("#####################\n#####################")
print(all(df_samples.schema.names == list_df_reports[0].columns))
print(all(snp_map.schema.names == snp_map_pandas.columns))
print(all(sample_map.schema.names == sample_map_pandas.columns))

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
print(n_unique_samples == sample_map_pandas.shape[0])
print(len(list_files_samples) == sample_map_pandas.shape[0])

#count the number of genotypes, i.e., total number of rows
n_genotypes = df_samples_subset.count()

#check
print("\n#####################\n#####################")
print("We have as many genotypes as samples*snps in the maps")
print("#####################\n#####################")
print(n_genotypes == snp_map.count()*sample_map.count())

#checks
print("\n#####################\n#####################")
print("Multiple checks within each sample")
print("#####################\n#####################")

#check that the sample ID is the same than that showed in the input file name
from pyspark.sql.functions import split, concat, col, lit
#define the two columns to be compared
col_input_name_check = split(col("input_file"), "file:" + temp_dir.name + "/" + batch_name + "_").getItem(1)
    #prepare a column from splitting the input_file column and get the second item, which is FinalReportXX.txt, where XX is the index of the sample. This column was previously created.
col_index_check = concat(lit("FinalReport"), col("Sample Index"), lit(".txt"))
    #prepare a column concatenating FinalReport and .txt as literal values to the sample index, so we have the same format than in the input file name
df_samples_subset \
    .withColumn("check_1", col_index_check == col_input_name_check) \
    .select("check_1") \
    .distinct() \
    .show()
    #to the orignal data.frame, add a new column comparing the two previously defined columns, select the column with the check, get the distinct values and show.
    #Only true should be present.

#check that each sample has the same number of rows
df_samples_subset \
    .groupBy("Sample Index") \
    .count() \
    .toDF("Sample Index", "count") \
    .withColumn("check_2", col("count") == snp_map.count()) \
    .select("check_2") \
    .distinct() \
    .show()
    #group rows per sample, count the number of rows per sample, convert to DF, then create a new column checking whether count is equal to the number of snps in the map, select that column and see if only true


#check that index, name, chromosome and position of all SNPs are the same than in the snp map. They should be in the exact same order in all files

#pdf = pd.read_csv(temp_dir.name+ "/" + batch_name + "_FinalReport24.txt", delimiter="\t", header=0, low_memory=False)


def calc_checks(pdf):

    df = pd.DataFrame(columns=["index", "check_3", "check_4", "check_5", "check_6"])

    df["index"] = pdf["Sample Index"]
    df["check_3"] = all(np.array(pdf["SNP Index"]) == np.array(snp_map_pandas["Index"]))
    df["check_4"] = all(np.array(pdf["SNP Name"]) == np.array(snp_map_pandas["Name"]))
    df["check_5"] = all(np.array(pdf["Chr"]) == np.array(snp_map_pandas["Chromosome"]))
    df["check_6"] = all(np.array(pdf["Position"]) == np.array(snp_map_pandas["Position"]))
    df["check_7"] = np.isin(np.array(pdf["Sample ID"][0]), np.array(sample_map_pandas["ID"]))
    df["check_8"] = np.isin(np.array(pdf["Sample Index"][0]), np.array(sample_map_pandas["Index"]))
    return df


schema_checks = StructType() \
    .add("index",IntegerType(), nullable=True) \
    .add("check_3", BooleanType(), nullable=True) \
    .add("check_4", BooleanType(), nullable=True) \
    .add("check_5", BooleanType(), nullable=True) \
    .add("check_6", BooleanType(), nullable=True) \
    .add("check_7", BooleanType(), nullable=True) \
    .add("check_8", BooleanType(), nullable=True)

import pyspark.sql.functions as F


cnt_cond = lambda cond: F.sum(F.when(cond, 1).otherwise(0))

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
        F.count("*").alias("check_9")) \
    .collect()
        #count is for the total number of genotypes, but maybe we need number of samples?


for index, check in enumerate(checks_raw[0]):
    print("CHECK " + str(index) + ": " + str(check == n_genotypes))



#aggreate grouups by counting only true cases?
    #https://stackoverflow.com/questions/49021972/pyspark-count-rows-on-condition





    #make a tuple with all checks and append it to empty list
    check_across_reports_list.append(tuple([check_1, check_2, check_3, check_4, check_5, check_6, check_7, check_8]))

#convert to DF and add column names
check_across_reports_df = pd.DataFrame(check_across_reports_list)
check_across_reports_df.columns = ["check_" + str(number) for number in np.arange(1, check_across_reports_df.shape[1]+1, 1)]
    #as many checks as columns we have in the DF + 1, because the end of the range is NOT included.

#see if all true across all rows of each column
print(check_across_reports_df.all(axis=0))

#check
print("\n#####################\n#####################")
print("The multiple checks are done across all samples?")
print("#####################\n#####################")
print(check_across_reports_df.shape[0] == len(list_df_samples))

#concat
print("\n#####################\n#####################")
print("concatenate all DFs into one single DF, i.e., all_reports")
print("#####################\n#####################")
all_reports = pd.concat(
    list_df_samples, 
    axis=0, #concatenate along rows
    ignore_index=True) #clear the existing index and reset it, so do we have different indexes for each combination of snp and sample
print(all_reports)

#check
print("\n#####################\n#####################")
print("all_reports has the correct number of rows?")
print("#####################\n#####################")
print(all_reports.shape[0] == sample_map.shape[0]*snp_map.shape[0])
    #number of SNPs should be the number of samples times the number of snps per sample
print(all_reports.shape[0] == sum([element.shape[0] for element in list_df_samples]))
    #we should have as many rows as the total sum of rows across the list of DFs

#checks samples
print("\n#####################\n#####################")
print("all_reports has the correct samples")
print("#####################\n#####################")
print(np.array_equal(all_reports["Sample ID"].unique(), sample_map["ID"].values))
print(np.array_equal(all_reports["Sample Index"].unique(), sample_map["Index"].values))

#checks snps
print("\n#####################\n#####################")
print("all_reports has the correct snps")
print("#####################\n#####################")
print(np.array_equal(all_reports["SNP Index"].unique(), snp_map["Index"].values))
print(np.array_equal(all_reports["SNP Name"].unique(), snp_map["Name"].values))
print(np.array_equal(all_reports["Chr"].unique(), snp_map["Chromosome"].unique()))
print(np.array_equal(all_reports["Position"].unique(), snp_map["Position"].unique()))

#see unique alleles
print("\n#####################\n#####################")
print("see unique alleles")
print("#####################\n#####################")
print(all_reports["Allele1 - Forward"].unique())
print(all_reports["Allele2 - Forward"].unique())

#remove the list of DFs
import gc
del(list_df_samples)
gc.collect()

#save all_reports
#all_reports.to_csv("data/genetic_data/illumina_batches/"+batch_name+"_all_samples.txt",
#    sep="\t",
#    index=False)


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
lgen_file = all_reports.loc[:, ["Sample ID", "SNP Name", "Allele1 - Forward", "Allele2 - Forward"]]

#remove all_reports
del(all_reports)
gc.collect()

# Change to "0" those genotype calls with "--" to match the format of Plink
#allele 1
lgen_file.loc[
    (lgen_file["Allele1 - Forward"] == "-") | 
    (lgen_file["Allele1 - Forward"] == "--"), 
    "Allele1 - Forward"] = "0"
#allele 2
lgen_file.loc[
    (lgen_file["Allele2 - Forward"] == "-") | 
    (lgen_file["Allele2 - Forward"] == "--"), 
    "Allele2 - Forward"] = "0"

#check
print("\n#####################\n#####################")
print("Check that all SNPs do NOT have '-' or '--' for allele 1 and 2")
print("#####################\n#####################")
all(~np.isin(lgen_file["Allele1 - Forward"].unique(), test_elements=["-", "--"]))
all(~np.isin(lgen_file["Allele2 - Forward"].unique(), test_elements=["-", "--"]))
    #"~" to negate 

#Add additional columns that are required for lgen files
lgen_file["FID"] = "combat"

# Reorder the columns
lgen_file = lgen_file[["FID", "Sample ID", "SNP Name", "Allele1 - Forward", "Allele2 - Forward"]]

#look
print(lgen_file)

#Save without header:
lgen_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + ".lgen.gz",
    sep="\t",
    header=None,
    compression='gzip',
    index=False)



#you can save each final report separately
#df_subset_2.write.option("header", True).option("delimiter", "\t").csv(temp_dir.name + "/" + batch_name + "full.txt")
    #https://stackoverflow.com/questions/47780397/saving-dataframe-records-in-a-tab-delimited-file

#or a single file
df_subset_2.coalesce(1).write.option("header", True).option("delimiter", "\t").csv(temp_dir.name + "/" + batch_name + "_full")
    #https://sparkbyexamples.com/spark/spark-write-dataframe-single-csv-file/
    #https://sparkbyexamples.com/spark/spark-repartition-vs-coalesce/#dataframe-%20coalesce

#https://spark.apache.org/docs/2.2.0/sql-programming-guide.html#creating-dataframes




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
fam_file["FID"] = "combat" #ID for the whole study
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
fam_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + ".fam.gz",
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
map_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + ".map.gz",
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
print(np.array_equal(lgen_file["SNP Name"].unique(), map_file["Name"].values))
print("samples are the same in lgen and fam files?")
print(np.array_equal(lgen_file["Sample ID"].unique(), fam_file["ID"].values))


#################################
# check with lgen and map files #
#################################

del([lgen_file, map_file, fam_file])
gc.collect()


#######################
# Convert to ped file #
#######################

#calculate the ped file using the lgen, map and fam files
import os
os.system("cd data/genetic_data/plink_inputs; gunzip -k " + batch_name + ".lgen.gz")
os.system("cd data/genetic_data/plink_inputs; gunzip -k " + batch_name + ".map.gz")
os.system("cd data/genetic_data/plink_inputs; gunzip -k " + batch_name + ".fam.gz")


os.system("cd data/genetic_data; plink --lfile plink_inputs/" + batch_name + " --recode --out plink_ped_files/" + batch_name + "_plink")
    #go to the folder with plink inputs
    #--lfile for loading the lgen file, which should be accompanied by a .fam and .map files having the same name except the extension.
        #https://www.cog-genomics.org/plink/1.9/formats
    #--recode creates a new text fileset, after applying sample/variant filters and other operations. By default, the fileset includes a .ped and a .map file, readable with --file.
        #https://www.cog-genomics.org/plink/1.9/data#recode
    #--out for the output name

os.system("cd data/genetic_data/plink_inputs; rm " + batch_name + ".lgen")
os.system("cd data/genetic_data/plink_inputs; rm " + batch_name + ".map")
os.system("cd data/genetic_data/plink_inputs; rm " + batch_name + ".fam")

os.system("cd data/genetic_data/plink_ped_files; gzip " + batch_name + "_plink.map")
os.system("cd data/genetic_data/plink_ped_files; gzip " + batch_name + "_plink.ped")


###compress with pigz?
#memory issiues, we woulld need between 5-10 times of RAM than data we have, so we have 130GB, requiring between 650 and 1300GB. We are using more in the HPC with optuna across 50 cores...
    #https://wesmckinney.com/blog/apache-arrow-pandas-internals/

#pyspark used in TDI?



#This creates several files with the name used in --out, i.e., "batch1_example_plink"

##.ped
    #Original standard text format for sample pedigree information and genotype calls. Normally must be accompanied by a .map file; Haploview requires an accompanying .info file instead. Loaded with --file, and produced by --recode.
    #Contains no header line, and one line per sample with 2V+6 fields (columns) where V is the number of variants (as many columns as variants). The first six fields are the same as those in a .fam file. The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.
    #If all alleles are single-character, PLINK 1.9 will correctly parse the more compact "compound genotype" variant of this format, where each genotype call is represented as a single two-character string. This does not require the use of an additional loading flag. You can produce such a file with "--recode compound-genotypes".
#.ped too many columns to be loaded

##map.
    #Variant information file accompanying a .ped text pedigree + genotype table. Also generated by "--recode rlist".
    #A text file with no header line, and one line per variant with the following 3-4 fields:
        #Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
        #Variant identifier
        #Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
        #Base-pair coordinate
    #All lines must have the same number of columns (so either no lines contain the morgans/centimorgans column, or all of them do)

#see the map file
map_plink = pd.read_csv("data/plink_inputs_example/batch1_example_plink.map",
    names=["chromosome", "snp_name", "genetic_position", "pair_base"],
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
print(map_plink)

#chromosomes
print("see the chromosome names and compare with the original map")
print(np.unique(map_plink["chromosome"]))
print(np.unique(snp_map["Chromosome"]))
    #we have 22 autosomals, then X, Y, XY, MT, i.e., a total of 26 chromosomes

#variant names
print("check that the snp names are the same than those in the original map")
print(all(map_plink["snp_name"].isin(snp_map["Name"])))
print(all(snp_map["Name"].isin(map_plink["snp_name"])))

#genetic position
print("genetic position is always zero?")
print(np.unique(map_plink["genetic_position"]) == 0)
    #all zero, which is considered as none

##nosex
    #List of samples with ambiguous sex codes
    #https://www.cog-genomics.org/plink/1.9/output



###maybe use pigz to compress the lgen file? ofr compressing seems to be much more faster because uses several cores, it is less useful for decompressing
    #https://stackoverflow.com/questions/12313242/utilizing-multi-core-for-targzip-bzip-compression-decompression




###ELIMINNA TEMP FOLDER

################################

#remove the temp dir
temp_dir.cleanup()

#stop spark env
spark.stop()







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

##.log
    #info about the session, random seed, input files..


#################
# Load ped file #
#################

#load the ped file in plink and calculate frequency of snps
os.system("cd data/plink_inputs_example; plink --file batch1_example_plink --freq --out  batch1_example_plink_analysis")
    #--file indicate name of the ped and map files.
    #--freq generates an allele frequency report.
    #--output is the prefix of the output files generated
    #https://www.cog-genomics.org/plink/1.9/general_usage

##frq file
    #Produced by --freq. Valid input for --read-freq.
    #A text file with a header line, and then one line per variant with the following six fields:
        #CHR Chromosome code
        #SNP Variant identifier
        #A1  Allele 1 (usually minor)
        #A2  Allele 2 (usually major)
        #MAF Allele 1 frequency
        #NCHROBS Number of allele observations

#I cannot load this file due to different delimiters in the file

##hh_file
    #Produced automatically when the input data contains heterozygous calls where they shouldn't be possible (haploid chromosomes, male X/Y), or there are nonmissing calls for nonmales on the Y chromosome.
    #A text file with one line per error (sorted primarily by variant ID, secondarily by sample ID) with the following three fields:
        #Family ID
        #Within-family ID
        #Variant ID
    #https://www.cog-genomics.org/plink/1.9/formats#hh

#see the file
hh_plink_analysis = pd.read_csv("data/plink_inputs_example/batch1_example_plink_analysis.hh",
    names=["FID", "ID", "snp_name"],
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
print(hh_plink_analysis)

#check
print("check that hh_file generated in analysis is identical to the hh_file generated when creating the .ped file")
print(hh_plink_analysis.equals(hh_plink))

##nosex
    #List of samples with ambiguous sex codes
    #https://www.cog-genomics.org/plink/1.9/output

#see the file
nosex_plink_analysis = pd.read_csv("data/plink_inputs_example/batch1_example_plink_analysis.nosex",
    names=["FID", "ID"],
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
print(nosex_plink_analysis)

#check
print("check that nosex_file generated in analysis is identical to the nosex_file generated when creating the .ped file")
print(nosex_plink_analysis.equals(nosex_plink))

##.log
    #info about the session, random seed, input files..




'''