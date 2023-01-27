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

#
print("\n#########print the working directory#########")
import os
print(os.getcwd())

#you can also use the "!" character, but only whithin ipython, it does not work if you run the script as a program
#message="hola mundo"
#!echo {message}
#!pwd
    #https://jakevdp.github.io/PythonDataScienceHandbook/01.05-ipython-and-shell-commands.html


#######################################
# Passing arugments of python program #
#######################################

#define input arugments to be passed when running this script in bash
#we will use a bash script to run this python program two times instead of doing a function and parallelize that function. In this way, the same python script is run for each batch and if we want to run again one batch, we just need to modify the bash script, not the python program. This makes sense because we have only two batches, and the parallelization will occur inside each batch. Importantly, we can also have separated .out files, so we can look at the checks separately.
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--batch_name", help="Name of the batch used as input")
args=parser.parse_args()

#get the arguments of the function that have been passed through command line
#batch_name = args.batch_name
batch_name = "ILGSA24-17303"
#batch_name = "ILGSA24-17873"
    #for debugging

#
print("########################################\n########################################")
print("Starting batch number: " + batch_name)
print("########################################\n########################################")


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
import zipfile
zipdata = zipfile.ZipFile("data/genetic_data/illumina_batches/" + zip_name + ".zip")
zipinfos = zipdata.infolist()

#get the name of each zip file
names_files = [zipinfo.filename for zipinfo in zipinfos]

#iterate across files and get only final reports
print("#####################\n#####################")
print("Unzipping data: ")
print("#####################\n#####################")
#we are using a loop because it is a bit faster than zipdata.extractall. Parallelizing does not seems to be a good solution for unzipping and I got problems with pool. This takes a few minutes anyways.
    #parallel slower
        #https://stackoverflow.com/questions/43313666/python-parallel-processing-to-unzip-files
    #count time
        #https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution

import numpy as np
#zipinfo=zipinfos[np.where([element == zip_name + "/SNP_Map.txt" for element in names_files])[0][0]]
#zipinfo=zipinfos[np.where([element == zip_name + "/Sample_Map.txt" for element in names_files])[0][0]]
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
print("#####################\n#####################")
print("the number of files extracted is the number of samples?: ")
print("#####################\n#####################")
if batch_name == "ILGSA24-17303":
    correct_number_files = 216 + 2
elif batch_name == "ILGSA24-17873":
    correct_number_files = 1248 + 2
    #plus 2 because we also extracted sample and snp maps
print(n_files == correct_number_files)


import glob
list_files = glob.glob(temp_dir.name + "/*")
list_files

#natural sorting, 1,10, 21.. that works with numbers + strings like 1b, 5c...
from natsort import natsorted
list_files = natsorted(list_files)
    #https://github.com/SethMMorton/natsort#installation



len(list_files) == correct_number_files

list_files_samples = [file for file in list_files if file.startswith(temp_dir.name + "/" + zip_name)]
list_files_samples

import pandas as pd
sample_map = pd.read_csv(temp_dir.name + "/Sample_Map.txt",
    header=0,
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
sample_map

len(np.unique(sample_map["ID"])) == len(list_files_samples)


snp_map = pd.read_csv(temp_dir.name + "/SNP_Map.txt",
    header=0,
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
snp_map

print("#####################\n#####################")
print("Reading the final reports into a list")
print("#####################\n#####################")

#sample_path = list_files_samples[0]
def read_final_reports(sample_path):

    print(sample_path.split(temp_dir.name+"/")[1])

    final_report = pd.read_csv(
        sample_path,
        skiprows=10, #skip rows with the header
        delimiter="\t", 
        low_memory=False)
        #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
    return final_report

#read_final_reports(list_files_samples[0])

import multiprocessing as mp

#open the pool with 5 cores
pool = mp.Pool(5)

#run the function to calculate distance of each gene to the closest
#selective pressure gene. This takes the gene IDs as inputs
#and the function will calculate the distance for each gene id
#https://stackoverflow.com/questions/64763867/parallel-processing-of-each-row-in-pandas-iteration
list_df_samples = pool.map(read_final_reports, list_files_samples[0:20])
#close the pool
pool.close()


len(list_df_samples) == sample_map.shape[0]

check_across_reports_list = []
#index=0; df_sample = list_df_samples[0]
for index, df_sample in enumerate(list_df_samples):

    check_1 = int(list_files_samples[index].split(temp_dir.name+"/"+zip_name+"_FinalReport")[1].split(".txt")[0]) == df_sample["Sample Index"].unique()[0]

    check_2 = list_df_samples[0].columns.equals(df_sample.columns)

    check_3 = snp_map["Name"].equals(df_sample["SNP Name"])
    check_4 = snp_map["Chromosome"].equals(df_sample["Chr"])
    check_5 = snp_map["Position"].equals(df_sample["Position"])
    check_6 = df_sample["Sample ID"].unique()[0] == sample_map.iloc[index, :].loc["ID"]
    check_7 = df_sample["Sample Index"].unique()[0] == sample_map.iloc[index, :].loc["Index"]
        #we have ordered the final_reports so first goes the report of sample 1, then 2...

    check_across_reports_list.append(tuple([check_1, check_2, check_3, check_4, check_5, check_6, check_7]))

check_across_reports_df = pd.DataFrame(check_across_reports_list)

check_across_reports_df.columns = ["check_" + str(number) for number in np.arange(1, check_across_reports_df.shape[1]+1, 1)]


check_across_reports_df.all(axis=0)

all_reports = pd.concat(list_df_samples, axis=0)


#check number of rows
all_reports.shape[0] == sample_map.shape[0]*snp_map.shape[0]


len(all_reports["Sample ID"].unique()) == sample_map.shape[0]
all(all_reports["Sample ID"].isin(sample_map["ID"]))
all(all_reports["Sample Index"].unique() == sample_map["Index"])


len(all_reports["SNP Name"].unique()) == snp_map.shape[0]
all(all_reports["SNP Name"].unique() == snp_map["Name"])

all(all_reports["Chr"].unique() == snp_map["Chromosome"].unique())
all(all_reports["Position"].unique() == snp_map["Position"].unique())




###########################################
#### General info about available data ####
###########################################

#Batches
    #Within combat_genes.zip, we have two batches (ILGSA24-17303 and ILGSA24-17873), being the data separated in these two.
    #In ILGSA24-17303.zip, we have the final reports for each 216 samples, along with the sample and snp maps.
        #In the initial_stuff folder there is a zip called "ILGSA24-17303.zip" that I may downloaded from the initial location where this data was stored in summer 2022. There are Plink files, but I am not sure this is the correct data and I cannot find the final_report files.
    #In 17873, we have the IDAT files with probs intensity from the microarrays used to genotype (first zips), the final reports+sample/snp maps (CAGRF20093767.zip) and a inputs for plink. But all of this only for 1248 individuals, not the whole cohort.
    #the phenotype csv has 1264 samples, not having phenotype data for 41 of them.
        #1248+216=1464

    #WARNING ABOUT UNZIPPING TO CHECK
        #warning [CAGRF20093767.zip]:  32332459452 extra bytes at beginning or within zipfile (attempting to process anyway)
        #error [CAGRF20093767.zip]:  start of central directory not found; zipfile corrupt. (please check that you have transferred or created the zipfile in the appropriate BINARY mode and that you have compiled UnZip properly)
        #CHECK WARNING
        #solution for the error: https://askubuntu.com/questions/54904/unzip-error-end-of-central-directory-signature-not-found



#I have received a Illumina report with 3 files ([link](https://www.biostars.org/p/51928/)):
    #The "FinalReport.txt" for Illumina raw genotype data generated from Genome Bead Studio for 2.5M (GSGT Version	2.0.4). This includes a hader with number of SNPs, samples.... and then the data with sample index, sample names, alleles... the first row includes the column names. This is a BPM file.
        #From this file, we can obtain the a lgen file. It is just taking the first few columns of the FinalReport. Plink has an option to read the .lgen format and convert it to PED file (see below; [link](https://www.biostars.org/p/13302/))
    #A SNP map file with the physical positions of the snps.
    #A sample map file with info for each individual, like the sex, ID..

#It is usually recommended to prepare this data to create a ped file with Plink, which is a tool to process genetic data ([link](https://www.cog-genomics.org/plink/)), perform some filtering and QC analyses and then export as VCF ([link](https://www.biostars.org/p/210516/), [link](https://www.biostars.org/p/135156/), [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/)) and use packages like Hail in python.
    #https://hail.is/docs/0.2/tutorials/01-genome-wide-association-study.html

#There is an interesting alternative, gtc2vcf ([link](https://github.com/freeseek/gtc2vcf), [link](https://software.broadinstitute.org/software/gtc2vcf/)), which can directly transform Ilumina reports into VCF files from command line. We are going to use Plink, though, because it is much more widely used and there are tutorials and best practice papers using it.

#In particular, we are going to use the a paper about QC by Ritchie.There is a first version 2011 ([link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/)) and a second version in 2022 ([link](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.603)).








############################
#### Plink Installation ####
############################

#I have downloaded Plink ([link](https://www.cog-genomics.org/plink/)) and copied the executable ("plink") to `bin`, so we can use Plink from any place just typing `plink`. We are using Plink version 1.9 (see singularity recipe for the version used).
os.system("plink --version")

#Note that there is Plink 1.9 ([link](https://www.cog-genomics.org/plink/1.9/)) and Plink 2.0 ([link](https://www.cog-genomics.org/plink/2.0/)), these are not connected but different programs. 
    #This [threat](https://www.biostars.org/p/299855/#:~:text=The%20main%20difference%20is%20that,for%20a%20while%20to%20come.) of biostars explains the differences:
        # "The main difference is that plink 1.9 is essentially finished, while plink 2.0 is an alpha-stage program which will have significant unfinished components for a while to come. As a consequence, current development priorities for plink 2.0 are centered around things which are impossible to do with plink 1.9, such as handling multiallelic/phased variants and dosage data and reliably tracking REF/ALT alleles; while things that plink 1.9 already handles perfectly well, such as working with .ped/.map file pairs, have been deliberately deprioritized for now.
        # So, **you should stick to 1.9 as long as it's good enough for the jobs you need to perform. But once you need to do something outside 1.9's scope, you're likely to find that 2.0 already has the additional feature you need** (or it'll be added quickly after you ask for it)"


##Plink dummy example

# Follow the general usage page ([link](https://www.cog-genomics.org/plink/1.9/general_usage)) of plink and run the toy example found in the downloaded folder.

#os.system("cd /media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/software/plink/plink_linux_x86_64_20221210; plink --file toy --freq --out toy_analysis')


# Explanation of the arguments
    # - "--file toy" tells PLINK to use the genomic data in the text files toy.ped and toy.map. You'll see several other ways to specify input data on the next page.
        #toy.map has two rows, one for each SNP ([link](https://www.cog-genomics.org/plink/1.9/formats#map)).
            #- Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
            #- Variant identifier
            #- Position in morgans or centimorgans (optional; also safe to use dummy value of          '0')
            #- Base-pair coordinate 
            #- All lines must have the same number of columns (so either no lines contain the morgans/centimorgans column, or all of them do).
        #- toy.ped. Contains no header line, and one line per sample with 2V+6 fields where V is the number of variants. 
            #- The first six fields are the same as those in a .fam file:
                #Family ID ('FID')
                #Within-family ID ('IID'; cannot be '0')
                #Within-family ID of father ('0' if father isn't in dataset)
                #Within-family ID of mother ('0' if mother isn't in dataset)
                #Sex code ('1' = male, '2' = female, '0' = unknown)
                #Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            #The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.
        #If all alleles are single-character, PLINK 1.9 will correctly parse the more compact "compound genotype" variant of this format, where each genotype call is represented as a single two-character string. This does not require the use of an additional loading flag. You can produce such a file with "--recode compound-genotypes". 
        #It is also possible to load .ped files missing some initial fields.
            # - --freq tells PLINK to generate an allele frequency report. The full range of supported calculations is summarized under "Main functions" in the sidebar, and the formats of all reports are described in the file formats appendix.
            # - --out is for the output file prefix.
 
# So this particular combination makes PLINK calculate allele frequencies in toy.ped + toy.map, and write a report to toy_analysis.frq.
    # Here ([link](https://www.cog-genomics.org/plink/1.9/general_usage)) you can find information about the flag usage in plink. As an example, we are going to work with --out flag.
    # By default, the output files generated by PLINK all have names of the form 'plink.<one of these extensions>'. This is fine for a single run, but as soon as you make more use of PLINK, you'll start causing results from previous runs to be overwritten.
    # Therefore, you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide. E.g. in the example above, "--out toy_analysis" caused PLINK to create a file named toy_analysis.frq instead of plink.frq.
    # Since the prefix is a required parameter, invoking --out without it will cause PLINK to quit during command line parsing:



#######################################
#### Exploration of a final report ####
#######################################

#Explore the final report file of one sample. 

#We are going to avoid the header, so we skip the first 10 rows leaving a total of 40 rows, i.e., 39 SNPs plus the columns names.

# In[2]:


import pandas as pd

final_report1 = pd.read_csv("data/example_data/ILGSA24-17303_FinalReport1.txt",
    skiprows=10, #skip rows with the header
    delimiter="\t", 
    low_memory=False) 
    #low_memory: Internally process the file in chunks, resulting in lower memory use while parsing, but possibly mixed type inference. To ensure no mixed types either set False, or specify the type with the dtype parameter. 
final_report1

#See columns names (see this [tutorial](https://www.youtube.com/watch?v=_-Gu9hO8aFM&ab_channel=GenomicsBootCamp)):
final_report1.columns


# **Sample identification**
# The first three columns are about the sample. At least in our case, we have 1 final report for each sample, i.e., each individual (see below), so these columns should have only one value in the entire file. The same sample index and ID. We do not have sample name.

import numpy as np

print(f'Unique cases are 1? {len(np.unique(final_report1["Sample Index"]))==1}')
print(f'Unique cases are 1? {len(np.unique(final_report1["Sample ID"]))==1}')
print(f'Unique cases are 1? {len(np.unique(final_report1["Sample Name"]))==1}')


# **SNPs identification**
# 
# Then, we have the SNP indexes and names. We should have a different value for each row, i.e., the number of unique cases should be the same than the number of rows, i.e., no duplicates. Therefore, **we have as many rows as SNPs**.

print(f'Unique cases are the number of rows? {len(np.unique(final_report1["SNP Index"]))==final_report1.shape[0]}')
print(f'Unique cases are the number of rows? {len(np.unique(final_report1["SNP Name"]))==final_report1.shape[0]}')


# **Chromosome data**
# 
# Then, we have the chromosome. We have some extrange cases like chromosome 0, XY and MT
np.unique(final_report1["Chr"])


# We have 210 SNPs out 600K that have chromosome zero, and they also have position zero. These seem to be outdated SNPs as they do not have any physical position.
f'SNPs with chromosome zero are less than 1000? {final_report1.loc[final_report1["Chr"] == "0"].shape[0] < 1000}'
f'Chromosome zero SNPs have also position zero? {np.unique(final_report1.loc[final_report1["Chr"] == "0", "Position"]) == 0}'

# We have also SNPs in chromosome XY.
print(final_report1.loc[final_report1["Chr"] == "XY"])
f'The number of XY SNPs is less than 1000: {final_report1.loc[final_report1["Chr"] == "XY"].shape[0]<1000}'
    #These seems to be pseudoautosomal regions (PAR), i.e., homologous sequences of nucleotides on the X and Y chromosomes. Therefore, this should be ok, the error is to put them in X or XY while the convention is to set them as XY ([link; section Sex-chromosome marker annotation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100670/)).

# We have also mitochondrial SNPs.
print(final_report1.loc[final_report1["Chr"] == "MT"])
f'The number of mito SNPs is less than 1500: {final_report1.loc[final_report1["Chr"] == "MT"].shape[0]<1500}'


# **Position**
# 
# Then, we have the position, which is an integer.
print(f'The type of Position is integer? {final_report1["Position"].dtype == "int64"}')
print(f'There are no NAs? {final_report1.loc[final_report1["Position"].isna()].shape[0] == 0}')


# **Quality measures**
# 
# Then, we have the GT and GC scores, which are quality measures.

# The most important QC parameter is the GenTrain score. The GenTrain score is computed from the GenTrain 2.0 clustering algorithm. It is a measurement of SNP calling quality, ranging from 0 to 1, with higher value meaning better quality ([link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/)).

print(f'The type of GT score is integer? {final_report1["GT Score"].dtype == "float64"}')
print(f'There are no NAs? {final_report1.loc[final_report1["GT Score"].isna()].shape[0] == 0}')
print(f'Min is equal of higher than 0? {np.min(final_report1["GT Score"]) >= 0}')
print(f'Max is equal or lower than 1? {np.max(final_report1["GT Score"]) <= 1}')


# The second most important QC parameter is the cluster separation score, which measures how well the AA, AB and BB clusters are separated. The cluster separation score also ranges from 0 to 1, with higher meaning better (more separation; [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/)). 
# 
# In the link cited, in figure 4, you can see several cases for specific snps where the separation between clusters (AA, AB and BB) is not clear (left panels). Meaning that there individuals are not clearly differentiated for this SNP. With realignement they get a better separation and a better score.

# NOT PRESENT IN MY DATASET

# The third most important QC parameter is call frequency, which measures the percentage of samples with successful calls for that SNP. The call frequency also ranges from 0 to 1, with higher meaning more samples have successful calls for this SNP ([link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/)).

print(f'The type of GT score is integer? {final_report1["GC Score"].dtype == "float64"}')
print(f'There are no NAs? {final_report1.loc[final_report1["GC Score"].isna()].shape[0] == 0}')
print(f'Min is equal of higher than 0? {np.min(final_report1["GC Score"]) >= 0}')
print(f'Max is equal or lower than 1? {np.max(final_report1["GC Score"]) <= 1}')


# **Log R ratio and BAF**
# 
# These metrics could be used to detect copy number variation (e.g., deletions) without having whole genome sequencing ([Nandolo et al 2018](https://www.youtube.com/watch?v=_-Gu9hO8aFM&ab_channel=GenomicsBootCamp)).


# **Alleles**
# 
# One of the biggest weakness of the Illumina genotyping array design is Illumina’s definition of strand. As DNA is double-stranded, significant SNPs from GWAS need to be presented with their strand information to properly report the risk alleles. This unfortunately has not been a standard practice. The most intuitive definition of strand is to use the human genome reference as the forward strand. Defying logic, Illumina introduced a more convoluted definition of strand: top and bottom [25], which has caused great confusion with reference to the forward and reverse strand [26, 27] ([link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/)).
# 
# DNA strand designations ([link](https://www.biostars.org/p/4885/)):
# - Customer Strand:Same as the Source strand. For custom content, it is the strand submitted by the customer for probe design.
# - ILMN Strand, a.k.a. Design Strand: The strand used by Illumina to design probes based on thermodynamic stability and locus specificity according to NCBI BLAST. For this reason, it can differ from the Customer/ Source strand.
# - Forward/Reverse (Fwd/Rev) Strand: Used by dbSNP, **Fwd/Rev designations can change with NCBI Genome Build updates, so Genome Build must be specified when reporting Fwd/Rev strands**. 1. For SNPs in standard array products, Fwd strand = Source strand, and originates from dbSNP. 2. For custom array product SNPs without rsid’s, the customer can identify the Source strand as Fwd or Rev, based on their own criteria. Illumina custom product files use the customer’s Fwd/Rev designations. Note: The Fwd strand, as identified in Illumina standard product files, should not be confused with Plus (+) strand, which HapMap interchangeably calls the “forward strand.”
# - Plus/Minus (+/-) Strand: The standard designation for all eukaryotic organisms used by HapMapand 1000 Genomes Project. The 5′ end of the (+) strand is at the tip of the short arm (p arm) of the chromosome and the 5′ end of the (-) strand is at the tip of the long arm (q arm). (+/-) designations can change with NCBI Genome Build updates, so Genome Build must be specified when reporting (+/-) strands.
# - Source Strand: Same as the Customer strand. The strand submitted to the Illumina designer for probe design. 1. For standard SNPs, it is the Fwd strand as reported in the source database (i.e., dbSNP). 2. Custom content can be reported as rsid’s or as the DNA sequences or chromosomal regions, depending on the format submitted by the customer.
# - Top/Bottom (Top/Bot) Strand:Top/Bot nomenclature was developed by Illumina using sequence-based context to assign strand designations that does not change regardless of database or genome assembly used. (e.g., depending on the NCBI Genome Build referenced, strand and allele designations can change). Top/Bot is not directly related to Fwd/Rev or (+/-).Top/Bot strand is determined by examining the SNP and the surrounding DNA sequence and it only applies to SNPs with two possible alleles. See the Top/Bot A/B Allele bulletin for more details ([link](https://www.illumina.com/documents/products/technotes/technote_topbot.pdf), [link](https://support.illumina.com/bulletins/2017/06/how-to-interpret-dna-strand-and-allele-information-for-infinium-.html)).

# See the different allele names for some SNPs with rs number:
final_report1.loc[final_report1["SNP Name"].apply(lambda x: "rs" in x), ['SNP Name', 'Allele1 - AB', 'Allele2 - AB', 'Allele1 - Top', 'Allele2 - Top', 'Allele1 - Forward','Allele2 - Forward','Allele1 - Design', 'Allele2 - Design','ILMN Strand', 'SNP']]
    #get several allele columns for tose rows having rs in the SNP name

# For example, rs9999844 has for foward, the aleles indicated in dbSNP. But this does not work always. For example, according to dbSNP, rs999994 is C/T, but we have G.
# 
# We will use Allele-Foward in general, but we will have to check the alleles at some point, maybe after selecting SNPs for the polygenic score?


# **Other columns**
# 
# There still other columns, probably related to the plots in the report of illumina, but we are not going to explore them for now unless the QC tutorials make use of them.
final_report1.columns


## Load and explore the sample and snp maps
# We need these maps in order to generate the inputs for Plink (see below).

# **SNP map**
snp_map = pd.read_csv("data/example_data/SNP_Map.txt",
    delimiter="\t", 
    header=0,
    low_memory=False) 
    #low_memory: Internally process the file in chunks, 
    #resulting in lower memory use while parsing, 
    #but possibly mixed type inference. To ensure no mixed 
    #types either set False, or specify the type with the dtype parameter. 
snp_map


# Check we have the same number of SNPs in the final report and the SNP map
f'Number of SNP matches between SNP map and final report? {snp_map.shape[0] == final_report1.shape[0]}'


# Therefore, we have a row per SNP, having each SNP an Index and a Name
print(f'Unique cases are the number of rows? {len(np.unique(snp_map["Index"]))==snp_map.shape[0]}')
print(f'Unique cases are the number of rows? {len(np.unique(snp_map["Name"]))==snp_map.shape[0]}')


# Each SNP has its chromose, having autosomals but also MT, zero cases and sex chromosomes..
np.unique(snp_map["Chromosome"])


# We have SNPS that have chromosome zero, and they also have position zero. They seem to be outdated SNPs.
print(f'SNPs with chromosome zero are less than 1000? {snp_map.loc[snp_map["Chromosome"] == "0"].shape[0] < 1000}')
print(f'Chromosome zero SNPs have also position zero? {np.unique(snp_map.loc[snp_map["Chromosome"] == "0", "Position"]) == 0}')


# We have also SNPs in chromosome XY.
print(snp_map.loc[snp_map["Chromosome"] == "XY"])
print(f'The number of XY SNPs is less than 1000: {snp_map.loc[snp_map["Chromosome"] == "XY"].shape[0]<1000}')
    # These seems to be pseudoautosomal regions (PAR), i.e., homologous sequences of nucleotides on the X and Y chromosomes. Therefore, this should be ok, the error is to put them in X or XY while the convention is to set them as XY ([link; section Sex-chromosome marker annotation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100670/)).

# We have also mitochondrial SNPs.
print(snp_map.loc[snp_map["Chromosome"] == "MT"])
print(f'The number of mito SNPs is less than 1500: {snp_map.loc[snp_map["Chromosome"] == "MT"].shape[0]<1500}')


# Then, we have GenTrain score for each SNP.
print(f'There is no NA? {all(snp_map["GenTrain Score"].notna())}')
print(f'Float? {snp_map["GenTrain Score"].dtype == "float64"}')
print(snp_map["GenTrain Score"].describe())


# There is also allele and strand information. Not sure what strand is shown, but we are going to use foward in final report anyways.

# The SNP I checked with mismatch between dbSNP and my database has also mismatch here.
snp_map.loc[snp_map["Name"] == "rs999994"]


# **Sample map**

sample_map = pd.read_csv("data/example_data/Sample_Map.txt",
    delimiter="\t", 
    header=0,
    low_memory=False) 
    #low_memory: Internally process the file in chunks, 
    #resulting in lower memory use while parsing, 
    #but possibly mixed type inference. To ensure no mixed 
    #types either set False, or specify the type with the dtype parameter. 
sample_map


# The number of samples should be 216 because this is the number of final reports present in the first batch
sample_map.shape[0] == 216


# We have index and ID for each sample
print(f'Unique cases are the number of rows? {len(np.unique(sample_map["Index"]))==sample_map.shape[0]}')
print(f'Unique cases are the number of rows? {len(np.unique(sample_map["ID"]))==sample_map.shape[0]}')


# We have gender information, whith male, female and Unknown. This is estimated from illumina, not from the study!!
print(np.unique(sample_map["Gender"]))


# There are three individuals with unknown sex according to Illumina data. These are indeed the three problematic individuals according to the logR value, which is above the Illumina expecation of 0.3.
print(sample_map.loc[sample_map["Gender"] == "Unknown"])
print(f'Unknown sex is below 10? {sample_map.loc[(sample_map["Gender"] == "Unknown")].shape[0] < 10}')

# Then, we have columns for ID of the parents, but this is not the case in our study.



################################
#### Explore phenotype data ####
################################

#This include reported sex and VO2 max data.
pheno_data = pd.read_csv("data/pheno_data/combact gene DNA GWAS 23062022_all_dna_samples.csv",
    delimiter=",", 
    header=0,
    low_memory=False) 
    #low_memory: Internally process the file in chunks, 
    #resulting in lower memory use while parsing, 
    #but possibly mixed type inference. To ensure no mixed 
    #types either set False, or specify the type with the dtype parameter. 
pheno_data

# We have gender (F/M) with some NAs
print(pheno_data["Gender"].describe())
print(f'NAs for sex: {pheno_data[pheno_data["Gender"].isna()].shape[0]}')

#the NAs for sex have no data at all
print(pheno_data.loc[pheno_data["Gender"].isna(),:])

# Also age, which is a float, as shown in the paper
pheno_data["Age"].describe()

# We also have the ID of the sample
print(pheno_data["AGRF code"])

# This code is the ID of the samples in the sample map of illumina.
# We have almost all samples of the first batch included in the pheno data, we only lack one individual.
print(f'Number of samples in sample map: {sample_map.shape[0]}')
print(f'Number of samples of pheno_data included in sample map')
print(pheno_data[pheno_data["AGRF code"].isin(sample_map["ID"])].shape[0])


# Then, different body mass and cardiorespiratory fitness variables.
# 
# NOT CHECKED FOR NOW.

#note that we have last rows with individuals of the study that have no phenotype data
print("Number of samples without phenotype data")
print(sum(pheno_data["Gender"].isna()))



####################################################
#### Conversion of final report to Plink inputs ####
####################################################

# We are going to convert the illumina report (BPM file) to lgen file that can be used as input in Plink. This seems to be trivial ([link](https://www.biostars.org/p/51928/)), but we are going to use a tutorial just in case ([link](https://www.youtube.com/watch?v=5_pNby7l2dE&t=1s&ab_channel=GenomicsBootCamp), [link](https://pastebin.com/pzgw7JVp)).

# The final report has 1 line for the SNP of one individual. We have separated final reports for each individual. This is close to the lgen format of plink.

# Therefore, our goal is to create lgen, also the fam and map files required to load the lgen in Plink. Some information might be missing in the final report, so you need to replace them.


#############
# lgen file #
#############

# **lgen file [plink info](https://www.cog-genomics.org/plink/1.9/formats#lgen)**
    # A text file with no header line, and one line per genotype call (or just not-homozygous-major calls if 'lgen-ref' was invoked) usually with the following five fields:
        #Family ID
        #Within-family ID
        #Variant identifier
        #Allele call 1 ('0' for missing)
        #Allele call 2

# As we have each sample in a separated final report, we need to bind them in order to get all the genotype calls

# Get first the paths for each final report
import glob
list_reports_files_full_path = glob.glob("data/example_data/ILGSA24-17303_FinalReport*") 
    #I prefer using the glob module, as it does pattern matching and expansion.
        #https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
list_reports_files_full_path

# Load the data.frame selecting only the colums needed for the lgen file
#using loop for now, BUT THIS SHOULD BE PARALLELIZED
list_df_reports = [
    pd.read_csv(
        path, 
        skiprows=10, #skip first 10 rows to get only the data
        usecols=[
            "Sample ID", #load only a subset of the columns
            "SNP Name", 
            "Allele1 - Forward", 
            "Allele2 - Forward"],                               
        header=0,
        sep='\t',
        low_memory=False) #to avoid problems with mixed types
    for path in list_reports_files_full_path]

# You can see how the selection of less columns decreases the size using the first final_report as example. We reduce the size almost 6 times.
import sys
print(sys.getsizeof(final_report1)/(10**6))
print(sys.getsizeof(final_report1[["Sample ID", "SNP Name", "Allele1 - Forward", "Allele2 - Forward"]])/(10**6))
    #get the size of final_report1 as a whole and then to only the selected columns
    #divide by 10^6 to get MBs

# The rest of necessary columns like sample ID or sample SNP and position/chromsome are just repeated across final_reports but in the sample and SNP map are not duplicated, so we can just use these smaller files to get these required columns. See next sections.

# We get a list with one DF per final report
print(list_df_reports[0])
print(f'Do we have as many DFs as paths we got? {len(list_df_reports) == len(list_reports_files_full_path)}')

# Concatenate, because we have one row per SNP and sample, so we can just concateneate DFs.
full_final_report = pd.concat(
    objs=list_df_reports, #list DFs
    axis=0, #concatenate along rows
    ignore_index=True) #clear the existing index and reset it

# We should have as many rows as the total sum of rows across the list of DFs
print(f'Do we have all the genotype calls? {full_final_report.shape[0] == sum([element.shape[0] for element in list_df_reports])}')

#take a look
print(full_final_report)

# We have I and D as genotype calls, we will have to check that
print(np.unique(full_final_report["Allele1 - Forward"]))

#Make a copy, using copy(), so modifications to the data or indices of the copy will not be reflected in the original object
lgen_file_raw = full_final_report.copy()

# Change to "0" those genotype calls with "--" to match the format of Plink
#allele 1
lgen_file_raw.loc[
    (lgen_file_raw["Allele1 - Forward"] == "-") | 
    (lgen_file_raw["Allele1 - Forward"] == "--"), 
    "Allele1 - Forward"] = "0"
#allele 2
lgen_file_raw.loc[
    (lgen_file_raw["Allele2 - Forward"] == "-") | 
    (lgen_file_raw["Allele2 - Forward"] == "--"), 
    "Allele2 - Forward"] = "0"

# Check that all SNPs do NOT have "-" or "--" for allele 1 and 2
all(~lgen_file_raw["Allele1 - Forward"].isin(["-", "--"]))
all(~lgen_file_raw["Allele1 - Forward"].isin(["-", "--"]))
    #"~" to negate 

# Add additional columns that are required for lgen files
lgen_file_raw["FID"] = "combat"

# Reorder the columns
lgen_file = lgen_file_raw[["FID", "Sample ID", "SNP Name", "Allele1 - Forward", "Allele2 - Forward"]]

#look
print(lgen_file)

#Save without header:
lgen_file.to_csv("data/plink_inputs_example/batch1_example.lgen",
    sep="\t",
    header=None,
    index=False)


############
# fam file #
############

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

#first check if we have samples with genetic but no phenotipic data
print("How many samples have genetic but no phenotipic data:")
print(f'{sum(~sample_map["ID"].isin(pheno_data["AGRF code"]))} out {sample_map.shape[0]} samples in the batch')

#merge the sample map and pheno_data retaining only those samples with ID in the batch
test_sample_map = sample_map.merge(
    pheno_data, 
    left_on="ID", #use the column with IDs in sample_map
    right_on="AGRF code", #use the column with IDs in pheno_data
    suffixes=["_sample_map", "_pheno_data"], #set the suffix for repeated columns (e.g., Gender)
    how="left") #only IDs from the sample_map

#recode sex column in sample_map to match coding of pheno_data
test_sample_map.loc[test_sample_map["Gender_sample_map"] == "Male", "Gender_sample_map"] = "M"
test_sample_map.loc[test_sample_map["Gender_sample_map"] == "Female", "Gender_sample_map"] = "F"
    #Male instead of M and Female instead of F

#see samples with different sex between datasets
unmatched_sex = test_sample_map.loc[
    (test_sample_map["Gender_pheno_data"] != test_sample_map["Gender_sample_map"]) & 
    ~test_sample_map["Gender_pheno_data"].isna(), 
    ["Gender_pheno_data", "Gender_sample_map", "ID", "AGRF code"]]

#we have some cases with problematic sex
print("how many sample have different sex between sample map (batch data) and phenotipic data?")
print(unmatched_sex.shape[0])
print("See these cases")
print(unmatched_sex)

#FOR NOW, we are going to use the sex reported in the phenotipic data, but we have to ask to David. We will add this sex to the sample_map

# We can get the self-reported sex and sample ID from the pheno data.
fam_file_raw = sample_map.loc[:, ["ID", "Gender"]]
print(fam_file_raw)

#remove those samples without phenotyipic data
#IMPORTANT HERE, we lose samples with genetic data because they do not have pheno data and hence, no sex data from the study can be obtained, only illumina
#fam_file_raw = fam_file_raw.loc[fam_file_raw["ID"].isin(pheno_data["AGRF code"]), :]

# Codify the sex variable following plink notation and considering the sex in pheno_data
fam_file_raw.loc[
    (fam_file_raw["ID"].isin(pheno_data.loc[pheno_data["Gender"].isna(), "AGRF code"])) | 
    (~fam_file_raw["ID"].isin(pheno_data["AGRF code"])), 
    "Gender"] = "0"
    #set as 0 the sex of those samples whose
        #ID is associated to NA gender in pheno_data
        #OR
        #ID is NOT included in pheno_data
fam_file_raw.loc[
    fam_file_raw["ID"].isin(pheno_data.loc[pheno_data["Gender"] == "M", "AGRF code"]), 
    "Gender"] = "1"
    #set as 1 the sex of those samples whose ID is associated with Gender=M in pheno_data
fam_file_raw.loc[fam_file_raw["ID"].isin(pheno_data.loc[pheno_data["Gender"] == "F", "AGRF code"]), "Gender"] = "2"
print(np.unique(fam_file_raw["Gender"]))

# Add the family variables and the phenotype (not added for now)
fam_file_raw["FID"] = "combat" #ID for the whole study
fam_file_raw["IID_father"] = "0" #'0' if father isn't in dataset
fam_file_raw["IID_mother"] = "0" 
fam_file_raw["phenotype"] = -9 #this is no data for phenotype

# Reorder
fam_file = fam_file_raw[["FID", "ID", "IID_father", "IID_mother", "Gender", "phenotype"]]
print(fam_file)

# Save without header:
fam_file.to_csv("data/plink_inputs_example/batch1_example.fam",
    sep="\t",
    header=None,
    index=False)


############
# map file #
############

# **Map file [Plink info](https://www.cog-genomics.org/plink/1.9/formats#map)**

# A text file with no header line, and one line per variant with the following 3-4 fields:
    # - Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
    # - Variant identifier
    # - Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
    # - Base-pair coordinate

# We can get this information from the SNP map.
    # We could also get this information from the final report but, as shown in previous lines, we have only loaded the columns required for the lgen file because we have to load hundreds of final reports and concatenate them, so it will use less memory with less columns. The Indexes and positions are not duplicated in SNP map, only 1 row per SNP, being a much smaller file.
print(final_report1)

# Make a copy with copy(), so modifications to the data or indices of the copy will not be reflected in the original object
map_file_raw = snp_map.copy()
print(map_file_raw)

#set genetic position as dummy
map_file_raw["centimorgans"] = 0 #dummy value of zero

#select columns
map_file = map_file_raw.loc[:, ["Chromosome", "Name", "centimorgans", "Position"]]
print(map_file)

# Save
map_file.to_csv("data/plink_inputs_example/batch1_example.map",
    sep="\t",
    header=None,
    index=False)


############
# map file #
############

#checks map and lgen_file
print("All snps of the map_file are included in the lgen_file?")
print(all(map_file["Name"].isin(lgen_file["SNP Name"])))
print("All snps of the lgen_file are included in the map_file?")
print(all(lgen_file["SNP Name"].isin(map_file["Name"])))

#checks fam_file and lgen_file
print("All samples of the fam_file are included in the lgen_file?")
print(all(fam_file["ID"].isin(lgen_file["Sample ID"])))
    #this can be FALSE because we have used only two final_reports, i.e., 2 individuals
print("All samples of the lgen_file are included in the fam_file?")
print(all(lgen_file["Sample ID"].isin(fam_file["ID"])))


#######################
# Convert to ped file #
#######################

#calculate the ped file using the lgen, map and fam files
os.system("cd data/plink_inputs_example; plink --lfile batch1_example --recode --out batch1_example_plink")
    #go to the folder with plink inputs
    #--lfile for loading the lgen file, which should be accompanied by a .fam and .map files having the same name except the extension.
        #https://www.cog-genomics.org/plink/1.9/formats
    #--recode creates a new text fileset, after applying sample/variant filters and other operations. By default, the fileset includes a .ped and a .map file, readable with --file.
        #https://www.cog-genomics.org/plink/1.9/data#recode
    #--out for the output name

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




###ELIMINNA TEMP FOLDER

################################

#remove the temp dir
temp_dir.cleanup()