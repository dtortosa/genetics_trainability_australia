#!/usr/bin/env python
# coding: utf-8


#########################################################
######## CONVERT ILLUMINA REPORT TO PLINK FORMAT ########
#########################################################



##################
#### Starting ####
##################

#set the working directory
import os
os.chdir("/home/dftortosa/singularity/australian_army_bishop/quality_control/")
os.getcwd()



###########################################
#### General info about available data ####
###########################################

#I have received a Illumina report with 3 files ([link](https://www.biostars.org/p/51928/)):
    #The "FinalReport.txt" for Illumina raw genotype data generated from Genome Bead Studio for 2.5M (GSGT Version	2.0.4). This includes a hader with number of SNPs, samples.... and then the data with sample index, sample names, alleles... the first row includes the column names. This is a BPM file.
        #From this file, we can obtain the a lgen file. It is just taking the first few columns of the FinalReport. Plink has an option to read the .lgen format and convert it to PED file (see below; [link](https://www.biostars.org/p/13302/))
    #A SNP map file with the physical positions of the snps.
    #A sample map file with info for each individual, like the sex, ID..

#It is usually recommended to prepare this data to create a ped file with Plink, which is a tool to process genetic data ([link](https://www.cog-genomics.org/plink/)), perform some filtering and QC analyses and then export as VCF ([link](https://www.biostars.org/p/210516/), [link](https://www.biostars.org/p/135156/), [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/)) and use packages like Hail in python.

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

# We can get the self-reported sex and sample ID from the pheno data.
fam_file_raw = pheno_data.loc[:, ["Gender", "AGRF code"]]
print(fam_file_raw)

# Codify the sex variable following plink notation
fam_file_raw.loc[fam_file_raw["Gender"].isna(), "Gender"] = "0"
fam_file_raw.loc[fam_file_raw["Gender"] == "M", "Gender"] = "1"
fam_file_raw.loc[fam_file_raw["Gender"] == "F", "Gender"] = "2"
print(np.unique(fam_file_raw["Gender"]))

# Add the family variables and the phenotype (not added for now)
fam_file_raw["FID"] = "combat" #ID for the whole study
fam_file_raw["IID_father"] = "0" #'0' if father isn't in dataset
fam_file_raw["IID_mother"] = "0" 
fam_file_raw["phenotype"] = -9 #this is no data for phenotype

# Reorder
fam_file = fam_file_raw[["FID", "AGRF code", "IID_father", "IID_mother", "Gender", "phenotype"]]
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



#POR AQUIIII


# CHECK THE THREE FILES

# Remove SNPs/samples not included in the lgen file!!! 

# In[ ]:





# Change to ped file with PLINK 

# In[181]:


import os 

os.system("ls data/plink_inputs_example")


# In[191]:


import os 

os.system("cd data/plink_inputs_example; plink --nonfounders --lfile batch1_example --recode --out batch1_example_plink")

#I have not checked this code


# This error is caused because in the phenotype file, row 1423 is empty and then we have rows with sample ID but no phenotype, so we should remove all rows from 1423, 

# In[ ]:


fam_file_raw.iloc[1422,:]


# In[ ]:





# Full tutorial for plink
# 
# https://genomicsbootcamp.github.io/book/

# In[ ]:





# In[ ]:





# Things to check when QC:
# 
# - Check the genome build. I have manually checked some SNPs and they have the position of hg38.
#     - YOU ARE NOT ADDED CENTIMORGANS, SO BE CAREFUL WITH WHAT PLINK DO ABOUT LD.
# - Remove zero chromosomes
#     - also MT, X, Y and XY?
# - Check strand, because the foward strand of illumina not always matches that of dbSNP
#     - StrandScript is a software to check and solve that
# - Check the illumina pdf report for the first batch (zip file in initial stuff), because they say there are three problematic samples.
#     - these three sample also have unknown sex according to illumina data!!
# - check sex between illumina and real sex data (pheno_data)
#     - remove uknown sex individuals?
# - check genotype calls that are I or D
#     - how plink deals with this?
# - ask David about the samples without phenotype in the excel file?
    #the last 42
#Ask david about the sample included in first bath but without phenotype data

# In[ ]:





# In[ ]:





# In[ ]:


get_ipython().run_cell_magic('bash', '', '\ncd /home/dftortosa/Desktop\n\nplink --file ILGSA24-17873 --freq --out ILGSA24-17873_analysis\n')


# In[ ]:





# In[ ]:


get_ipython().run_cell_magic('bash', '', '\ncd data/example_data\n\n#https://www.cog-genomics.org/plink/1.9/formats#lgen\n\n#plink --lfile example_ILGSA24-17303_FinalReport1.txt --recode\n\n#https://www.biostars.org/p/51928/\n')


# In[ ]:





# In[ ]:





# There are two batches (ILGSA24-17303 and ILGSA24-17873), being the data separated for these. 
# 
# - In ILGSA24-17303.zip, we have the final reports for each 216 samples, along with the sample and snp maps.
#     - In the initial_stuff folder there is a zip called "ILGSA24-17303.zip" that I may downloaded from the initial location where this data was stored in summer 2022. There are Plink files, but I am not sure this is the correct data and I cannot find the final_report files.
# 
# - In 17873, we have the IDAT files with probs intensity from the microarrays used to genotype (first zips), the final reports (CAGRF20093767.zip) and a inputs for plink. But all of this only for 1248 individuals, not the whole cohort.
#     - CAGRF20093767.zip includes the final reports of 1248 individuals, along with the sample and snp maps.

# warning [CAGRF20093767.zip]:  32332459452 extra bytes at beginning or within zipfile
#   (attempting to process anyway)
# error [CAGRF20093767.zip]:  start of central directory not found;
#   zipfile corrupt.
#   (please check that you have transferred or created the zipfile in the
#   appropriate BINARY mode and that you have compiled UnZip properly)
# 
# 
# CHECK WARNING
# 
# solution for the error: https://askubuntu.com/questions/54904/unzip-error-end-of-central-directory-signature-not-found

# In[ ]:





# In[ ]:





# - Quality Control Procedures for Genome-Wide Association Studies
# - A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis
# - Genetic prediction of complex traits with polygenic scores: a statistical review
# - Addressing the challenges of polygenic scores in human genetic research
