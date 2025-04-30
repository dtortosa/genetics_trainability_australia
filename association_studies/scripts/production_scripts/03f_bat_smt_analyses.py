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



###################################################################
######## ANALYZE ASSOCIATIONS WITHIN RELEVANT SET OF GENES ########
###################################################################

#We are going to analyze the associations within relevant set of genes like the BAT connectome.






########################
# region SUMAMRY STEPS #
########################

#FOR BAT ANALYSES
    #this would an additional step in this project that would be outside of the paper
    #take the 1000kb gene windows for all coding genes, liftover to hg38. If the USCS tool accepts genomic ranges, just use them as input
    #for each phenotype (VO2, beep....), calculate the average effect size for the association of SNPs inside each gene. better than median because want influence of outliers within gene like in iHS, if a SNPs is veery important in a gene that should influence the info about the whole gene
    #then, calculate 100000 random sets of genes, within each set, calculate the median effect of all genes inside the set and compare with the BAT set to obtain a distribution and empirical p-value (is effect higher in BAT). Here I want median because I do not want a gene outliser change things, I want the overall impact of BAT genes in general, not just a few genes.

# endregion






#######################################
# region INITIAL ANOTATIONS AND STEPS #
#######################################


###########
# imports #
###########

import pandas as pd
import numpy as np
import argparse
from multiprocessing import Pool
from functools import partial



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
    elif ("warning" in complete_process.stderr):
        
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
parser=argparse.ArgumentParser()
parser.add_argument("--pressure_name", type=str, default="bat", help="Selective pressure to analyze. String always, None does not work!")
parser.add_argument("--connectome_percentile", type=int, default=10, help="Percentile of top genes closest to the core gene of the corresponding connectome (BAT or SMT). Integer always, None does not work!")
parser.add_argument("--n_iterations", type=int, default=100000, help="Number of iterations to calculate the empirical p-value about differences between set of genes and the random expectation. Integer always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
pressure_name = args.pressure_name
connectome_percentile = args.connectome_percentile
n_iterations = args.n_iterations

# endregion





######################################################################
# region CALCULATE THE AVERAGE SNP EFFECT PER GENE ACROSS PHENOTYPES #
######################################################################
print_text("Calculate the average SNP effect per gene across phenotypes", header=1)
print_text("data preparation", header=2)
print_text("load the gene coordinates file in hg19", header=3)
gene_coords = pd.read_table(\
    "../../../../../postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt", \
    sep="\t", \
    header=0, \
    low_memory=False)
print(gene_coords)

print_text("process the file", header=3)
print_text("Select only one row per gene because in this file, each row is an exon, but all the rows of the same gene have the same middle gene position and gene windows as all the rows belongs to the same gene", header=4)
gene_coords_no_duplicated = gene_coords[~gene_coords["gene_id"].duplicated(keep="first")]
    #set as True all duplicates except the first occurrence.
    #then negate with "~", so all duplicates are False except the first occurrence, leading to select only that first occurrence
print(gene_coords_no_duplicated)

print("Check that the subset has as many rows as unique gene IDs in the original file. If True, remove the original gene coirdinate file")
check_gene_coords_subset = len(gene_coords["gene_id"].unique()) == gene_coords_no_duplicated.shape[0]
if(check_gene_coords_subset):
    del(gene_coords)
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE SUBSET OF THE GENE COORDINATE FILE")

print_text("select only the columns we need", header=4)
gene_coords_no_duplicated_subset = gene_coords_no_duplicated[[ 
    "chromosome_name", \
    "gene_id", \
    "hgnc_symbol", \
    "gene_start", \
    "gene_end", \
    "middle_point", \
    "lower_end_window_1000kb", \
    "upper_end_window_1000kb" \
]]
print(gene_coords_no_duplicated_subset)
    #We are only using 1000kb windows as these where we assessed positive selection with flexsweep and they are the ones with the most power

print_text("check that the gene_id has no NA", header=4)
if gene_coords_no_duplicated_subset["gene_id"].isna().sum()!=0:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE SUBSET OF THE GENE COORDINATE FILE")

print_text("remove genes with NA for the 1000kb windows", header=4)
gene_coords_no_duplicated_subset = gene_coords_no_duplicated_subset.dropna(
    subset=["lower_end_window_1000kb", "upper_end_window_1000kb"]
)

print_text("add the string 'chr' to the chromosome numbers", header=4)
gene_coords_no_duplicated_subset["chromosome_name"] = gene_coords_no_duplicated_subset["chromosome_name"].astype(str)
gene_coords_no_duplicated_subset.loc[:, "chromosome_name"] = "chr" + gene_coords_no_duplicated_subset["chromosome_name"]
    #before adding the string you need to convert the column to string (it was initially a integer)

print_text("convert the coordinates of the windows to 0-based to match the format of BED files in USCS", header=4)
gene_coords_no_duplicated_subset["lower_end_window_1000kb"] = gene_coords_no_duplicated_subset["lower_end_window_1000kb"] - 1
    #In liftover, supported formats include BED (e.g. "chr4 100000 100001", 0-based) and position box ("chr4:100,001-100,001", 1-based).
        #As you can see, in BED files the start is 0 while the end is the same than in 1-based coordinates. So we have to remove 1 from the start to convert it to 0-based.
        #https://genome.ucsc.edu/cgi-bin/hgLiftOver
    #I am assuming that our gene windows (calculated by me) are 1-based as in conserved elements density script we calculated the average density inside the windows after summing 1 to the start of the uscs data
        #see "/home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/method_deep/scripts/genomic_features_calculations/gene_windows/cons_elements_density/cons_elements_density_v3.R"

print_text("check the gene windows are float with .0 only, no more decimals", header=4)
if (\
    (sum(gene_coords_no_duplicated_subset["lower_end_window_1000kb"] % 1 != 0)!=0) |
    (sum(gene_coords_no_duplicated_subset["upper_end_window_1000kb"] % 1 != 0)!=0)
):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE COORDINATES OF THE GENE WINDOWS. THEY ARE NOT INTEGERS")
    #calculate the remainder of each value when divided by 1, that should be always 0, if not, stop

print_text("convert the coordinates of the windows to int", header=3)
gene_coords_no_duplicated_subset["lower_end_window_1000kb"] = gene_coords_no_duplicated_subset["lower_end_window_1000kb"].astype(int)
gene_coords_no_duplicated_subset["upper_end_window_1000kb"] = gene_coords_no_duplicated_subset["upper_end_window_1000kb"].astype(int)
    #they originally had a ".0" at the end because they were float numbers, but we need integers for liftover

print_text("print the updated data frame", header=3)
print(gene_coords_no_duplicated_subset[["chromosome_name", "gene_id", "hgnc_symbol", "lower_end_window_1000kb", "upper_end_window_1000kb"]])

print_text("save the DF to be used as input in liftover", header=3)
run_bash(" \
    mkdir \
        -p \
        ./data/bat_smt_analyses/ \
")
gene_coords_no_duplicated_subset[["chromosome_name", "lower_end_window_1000kb", "upper_end_window_1000kb", "gene_id"]].to_csv( \
    "./data/bat_smt_analyses/gene_windows_hg19.bed", \
    sep="\t", \
    header=False, \
    index=False \
)
    #The BED format is a tab-separated file with the following columns:
        #chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671). Many assemblies also support several different chromosome aliases (e.g. '1' or 'NC_000001.11' in place of 'chr1').
        #chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        #chromEnd - The ending position of the feature in the chromosome or scaffold.
        #name - Defines the name of the BED line. In our case, the gene ID.
        #https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#Liftover
        #https://genome.ucsc.edu/FAQ/FAQformat.html#format1

print_text("do the liftover from hg19 to hg38", header=3)
print_text("run liftover with the following parameters", header=4)
#From hg19 to hg38
#Minimum ratio of bases that must remap: 0.95 (default)
    #The minimum ratio of basepairs of the input region covered by an alignment. Regions scoring lower than this will not be lifted at all.
#Regions defined by chrom:start-end (BED 4 to BED 6; our case):
    #Keep original positions in output: NO
        #Lifted items for BED4 and up will include their original positions as part of their output names to assist in determining what got mapped where (in case multiple items have the same name in the input). Coordinates are 1-based fully closed, so the BED entry "chr1 100 150 item1" will be labeled "chr1:101-150:item1".
        #We are adding the gene ID as the name field after the coordinates, so we do not need this
    #Allow multiple output regions: NO
        #By default, input regions that map to multiple regions will not be lifted at all. When this option is checked, all targets are output.
        #We do not this as it will cause confusion, we just one regions with one to one match
#input file: gene_windows_hg19.bed
    #see above for the columns and the BED format
#https://genome.ucsc.edu/cgi-bin/hgLiftOver

print_text("load results liftover", header=4)
gene_windows_hg38 = pd.read_csv(\
    "./data/bat_smt_analyses/gene_windows_hg38.bed", \
    sep="\t", \
    header=None, \
    low_memory=False \
)

print_text("see lost genes to lack of alignemnt between ensembles", header=4)
n_lost_genes = gene_coords_no_duplicated_subset.shape[0] - gene_windows_hg38.shape[0]
if(n_lost_genes/gene_coords_no_duplicated_subset.shape[0]*100 > 10):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE LIFTOVER. MORE THAN 10% OF THE GENES WERE LOST")

print_text("set column names", header=4)
gene_windows_hg38.columns = ["chromosome_name", "lower_end_window_1000kb", "upper_end_window_1000kb", "gene_id"]
if(sum(gene_windows_hg38["lower_end_window_1000kb"] >= gene_windows_hg38["upper_end_window_1000kb"])!=0):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE LIFTOVER. THE LOWER END OF THE WINDOW IS GREATER THAN THE UPPER END FOR SOME GENES")

print_text("convert the start of the windows to 1-based like the rest of our data", header=4)
gene_windows_hg38["lower_end_window_1000kb"] = gene_windows_hg38["lower_end_window_1000kb"] + 1
    #Remember that in the liftover we are working with 0-based coordinates, but in the rest of our data we are working with 1-based coordinates.

print_text("define function to calculate the average effect size obtained from the GWAS per phenotype and gene", header=3)
#phenotype_name = "weight_change"
def process_phenotype(phenotype_name):

    #loadd the assoc file
    assoc_results = pd.read_csv( \
        "./results/final_results/analysis_full_data/" + phenotype_name + "/" + phenotype_name + "_small_set_predictors_set_linear_raw.assoc.gz", \
        sep="\t", \
        header=0, \
        low_memory=False, \
    )
    if(assoc_results.shape[1] != 13):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE ASSOC FILE. IT HAS NOT 13 COLUMNS")

    #add the string "chr" to the chromosome numbers
    assoc_results["Chromosome"] = "chr" + assoc_results["Chromosome"].astype(str)

    #empty list to store the results
    results = []

    #Group genes by chromosome for efficient processing
    #chromosome, gene_data = [i for i in gene_windows_hg38.groupby("chromosome_name")][0]
    for chromosome, gene_data in gene_windows_hg38.groupby("chromosome_name"):
        
        #filter assoc_results for the current chromosome
        snps_in_chromosome = assoc_results[assoc_results["Chromosome"] == chromosome]
        if(snps_in_chromosome["Chromosome"].unique() != chromosome):
            raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CHROMOSOME NAME IN THE ASSOC FILE. IT IS NOT THE SAME AS IN THE GENE COORDINATE FILE")

        #iterate over each gene in the chromosome within the coordinate file
        #_, gene_row = [i for i in gene_data.iterrows()][0]
        for _, gene_row in gene_data.iterrows():
            
            #extract the genomic window for the gene
            lower_bound = gene_row["lower_end_window_1000kb"]
            upper_bound = gene_row["upper_end_window_1000kb"]

            #filter SNPs within the genomic window
            snps_in_window = snps_in_chromosome[
                (snps_in_chromosome["Basepair"] >= lower_bound) &
                (snps_in_chromosome["Basepair"] <= upper_bound)
            ]

            #calculate the absolute effect size
            effect_sizes_abs = snps_in_window["Effect"].abs()
                #we do not care about the sense (ref or alt) of the effect, just care about the absolute value, the difference respect to zero

            #calculate mean and SD of the absolute effect sizes column
            mean_effect = effect_sizes_abs.mean() if not snps_in_window.empty else np.nan
            sd_effect = effect_sizes_abs.std() if not snps_in_window.empty else np.nan
                #prefer mean over median because we want to take into account the outliers, as in iHS, if a SNP is very important in a gene that should influence the info about the whole gene

            #create a result dictionary for the gene
            results_row = {
                "phenotype": phenotype_name,
                "chromosome_name": chromosome,
                "gene_id": gene_row["gene_id"],
                "mean_effect": mean_effect,
                "sd_effect": sd_effect
            }
            results.append(results_row)

    #check we have analyzed all the genes in the gene windows file
    if len(results) != gene_windows_hg38.shape[0]:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF GENES IN THE RESULTS")

    return results

print_text("run the function", header=3)
print_text("get the list of phenotypes", header=4)
phenotype_names = ["distance_change", "vo2_change", "beep_change", "weight_change"]

print_text("use multiprocessing Pool to process each phenotype in parallel", header=4)
with Pool(len(phenotype_names)) as pool:
    all_results = pool.map(process_phenotype, phenotype_names)
if(len(all_results) != len(phenotype_names)):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF PHENOTYPES IN THE RESULTS")

print_text("flatten the list of lists into a single list", header=4)
flattened_results = [item for sublist in all_results for item in sublist]

print_text("convert the results into a DataFrame", header=4)
results_df = pd.DataFrame(flattened_results)
if(results_df.shape[0] != gene_windows_hg38.shape[0]*len(phenotype_names)):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF GENES IN THE RESULTS")
else:
    print(results_df)

print_text("merge the results with the original coordinate files having gene symbols", header=4)
results_df_symbol = results_df.merge(gene_coords_no_duplicated_subset[["gene_id", "hgnc_symbol"]], how="left", on="gene_id")
    #we only retain those genes present in the results_df
        #left: left: use only keys from left frame, similar to a SQL left outer join; preserve key order
    #merge based on "gene_id" column
        #on: Column or index level names to join on. These must be found in both DataFrames.

# endregion






######################################################################
# region calculate enrichment of effect size in interest set of genes #
######################################################################
print_text("calculate enrichment of effect size in interest set of genes", header=2)
#calculate the median effect size across genes within a set of interest genes and then test whether it is higher than the median effect size across genes in random sets of genes
print_text("first extract the connectome of the interest gene", header=3)
if(pressure_name == "bat"):
    core_gene = "UCP1"
elif(pressure_name == "smt"):
    core_gene = "SLN"
elif(pressure_name == "metabolic"):
    core_gene = "NA"

print_text("create folder for the selective pressure", header=3)
run_bash(" \
    mkdir \
        -p \
        ./data/bat_smt_analyses/" + pressure_name + "_distance/; \
")

print_text("do operations specific of the connectome sets", header=3)
if pressure_name in ["bat", "smt"]:
    
    print_text("load the connectome having " + core_gene + " as core gene", header=4)
    run_bash(" \
        unzip \
            -p \
            ./data/bat_smt_analyses/all_human_gene-specific_connectomes_122015.zip \
            " + core_gene + ".txt > \
        ./data/bat_smt_analyses/" + pressure_name + "_distance/" + core_gene + ".txt \
    ")
        #extract the connectome of the gene from a Zip fila with ALL individual connectomes
            #downloaded from 
                #https://lab.rockefeller.edu/casanova/HGC
                #Human gene-specific connectomes - download specific genes
                #last file:
                    #all_human_gene-specific_connectomes_122015
        #unzip -p: 
            #extract files to pipe, no messages
            #we can get a specific file within the zip
                #https://unix.stackexchange.com/a/14125
        #We get a warning with both UCP1 and SNL
            #warning [./data/all_human_gene-specific_connectomes_122015.zip]:  4294967296 extra bytes at beginning or within zipfile (attempting to process anyway)
            #the returncode is equal to 1.
                #According to the man, this is not serious: 
                    #one or more warning errors were encountered, but processing completed successfully anyway. This includes zipfiles where one or more files was skipped due to unsupported compression method or encryption with an unknown password.
                        #https://linux.die.net/man/1/unzip
                #this is not error but warning, error would be returncode 2
                    #a generic error in the zipfile format was detected. Processing may have completed successfully anyway; some broken zipfiles created by other archivers have simple work-arounds.
            #We are unzipping just 1 file, and the unzipping was successful, also this is only a warning, so I think we are good here.
                #I am going to check below in the case of UCP1 is the file is identical to the original used for the BAT analyses. If so, the fact that we get a warning is not affecting.

    print_text("load the connectome having " + core_gene + " as core gene", header=4)
    pressure_conn = pd.read_csv( \
        "./data/bat_smt_analyses/" + pressure_name + "_distance/" + core_gene + ".txt", \
        sep="\t", \
        header=0, \
        low_memory=False)
    print(pressure_conn)
    if(pressure_conn.shape[1] != 10):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CONNECTOME FILE. IT HAS NOT 10 COLUMNS")

    print_text("check we have selected the correct connectome, i.e., the one with the correct core gene", header=4)
    if(pressure_conn["Target"].iloc[0] != core_gene):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CONNECTOME FILE. IT IS NOT THE ONE WITH THE CORRECT CORE GENE")

    print_text("if UCP1 is the core gene, check that the connectome extracted from the zip is the same as the original used for BAT analyses", header=4)
    if(pressure_name == "bat"):
        old_ucp1_conn = pd.read_csv( \
            "../../../../human_genome_connectome/bat_connectome/data/human_connectome/UCP1.txt", \
            sep="\t", \
            header=0, \
            low_memory=False)
        ucp1_2020 = pd.read_csv( \
            "../../../../human_genome_connectome/bat_connectome/data/human_connectome/UCP1_downloaded_2020.txt", \
            sep="\t", \
            header=0, \
            low_memory=False)
            #ucp1 connectome again but downloaded in 2020 (19/06/2020
        if( \
            (not pressure_conn.equals(old_ucp1_conn)) |
            (not pressure_conn.equals(ucp1_2020)) \
        ):
            raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CONNECTOME FILE. IT IS NOT THE ONE WITH THE CORRECT CORE GENE")
    else:
        print("We are not working with BAT, but with " + pressure_name)

    print_text("if UCP1 is the core gene, check that file with BAT relationships has the same genes than the connectome 1% obtained now", header=4)
    if(pressure_name == "bat"):
        bat_relationship = pd.read_csv( \
            "../../../../human_genome_connectome/bat_connectome/results/connectome_results/tables/appendix_S1_ordered.csv", \
            sep=",", \
            header=0, \
            low_memory=False)
        from natsort import natsort_keygen
        if(not pressure_conn \
            .loc[ \
                pressure_conn["Target_in_source_P-value(percentile)"] < 0.01, \
                "Target"] \
            .sort_values( \
                axis=0, 
                key=natsort_keygen(), \
                ignore_index=True) \
            .equals( \
                bat_relationship["Genes"]\
                .sort_values( \
                    axis=0, 
                    key=natsort_keygen(), \
                    ignore_index=True)) \
        ):
            raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CONNECTOME FILE. IT IS NOT THE ONE WITH THE CORRECT CORE GENE")
                #from the UCP1 connectome
                    #select those rows for which the p-value percentile of the gene is below 1% and get the gene name
                #natural sort the gene names
                    #"by" is not needed here because we only have 1 column
                    #axis=0: rows
                    #key: Apply the key function to the values before sorting.
                        #natsort_keygen()
                            #Generate a key to sort strings and numbers naturally.
                    #ignore_index:
                        #If True, the resulting axis will be labeled 0, 1, …, n - 1
                        #so you avoid the previous index
                #check that the resulting series is equals to the Genes included in BAT relationships after sorting in the same way
                    #https://stackoverflow.com/a/63890954/12772630
    else:
        print("We are not working with BAT, but with " + pressure_name)

    print_text("Percentile " + str(connectome_percentile) + "% :  select the interest genes", header=4)
    selected_connectome_genes = pressure_conn \
        .loc[ \
            pressure_conn["Target_in_source_P-value(percentile)"] < (connectome_percentile/100), \
            "Target"]
    print(selected_connectome_genes)

elif(pressure_name == "metabolic"):

    print_text("copy the file with the list of metabolic genes")
    run_bash(" \
        cp \
            ../../../../../postdoc_enard_lab/projects/ancient_selection_dating/data/metabolic_genes/metabolic_gene_list.txt.gz \
            ./data/bat_smt_analyses/" + pressure_name + "_distance/; \
    ")

    print_text("read the metabolic genes file", header=4)
    metabolic_genes = pd.read_csv( \
        "./data/bat_smt_analyses/" + pressure_name + "_distance/metabolic_gene_list.txt.gz", \
        sep="\t", \
        header=0, \
        low_memory=False \
    )

    print_text("check the column names", header=4)
    metabolic_genes.columns=["hgnc_symbol", "gene_id"]

    print_text("extract the gene symbols", header=4)
    selected_connectome_genes = metabolic_genes["hgnc_symbol"]
    print(selected_connectome_genes)

print_text("check we have not lost many genes for not having hg38 coordinates", header=4)
n_connectome_genes_in = sum(np.isin(results_df_symbol["hgnc_symbol"].unique(), selected_connectome_genes))
if (len(selected_connectome_genes)-n_connectome_genes_in)/len(selected_connectome_genes)*100 > 15:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CONNECTOME GENES. MORE THAN 15% OF THEM ARE NOT IN THE RESULTS")

print_text("calculate metrics in the interest set of genes", header=3)
print_text("empty list for results", header=4)
connectome_reference = []

print_text("iterate over phenotypes", header=4)
#pheno=phenotype_names[0]
for pheno in phenotype_names:
    
    #select the mean effects for the selected phenotype and within the set of interest genes
    results_df_symbol_subset = results_df_symbol.loc[ \
        (results_df_symbol["phenotype"]==pheno) & \
        (results_df_symbol["hgnc_symbol"].isin(selected_connectome_genes)), \
    :]
    if(results_df_symbol_subset.shape[0]!=n_connectome_genes_in):
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM CALCULATING AVERAGE EFFECT FOR THE INTEREST SET OF GENES")
        #check we have the final expected number of connectome genes after counting only those with hg38 coordinates

    #make a disct with the results
    results_connectome = { \
        "pheno": pheno, \
        "median_avg_effect_connectome": results_df_symbol_subset["mean_effect"].median(), \
        "median_std_effect_connectome": results_df_symbol_subset["sd_effect"].median() \
    }

    #append to the list
    connectome_reference.append(results_connectome)

print_text("convert results to DF", header=4)
connectome_reference_df=pd.DataFrame(connectome_reference)
connectome_reference_df

print_text("calculate empirical p-value about the interest set of genes having a higher effect compared to random set of genes", header=3)
print_text("define function", header=4)
#pheno=phenotype_names[0]
def process_phenotype_randomization(pheno, n_iterations):
    
    #subset the results for the current phenotype and NOT within the interest set of genes
    phenotype_results = results_df_symbol[ \
        (results_df_symbol["phenotype"] == pheno) & \
        (~results_df_symbol["hgnc_symbol"].isin(selected_connectome_genes)) \
    ]

    #get the median average effect for the connectome genes for this phenotype
    connectome_median = connectome_reference_df.loc[
        connectome_reference_df["pheno"] == pheno, "median_avg_effect_connectome"
    ].values[0]

    #initialize a counter for how many random medians are higher than the connectome median
    count_higher = 0

    #empty list to save the median value in each random set
    random_median_list = []

    #perform iterations
    #iteration=[i for i in range(n_iterations)][0]
    for iteration in range(n_iterations):
        
        #randomly sample rows the phenotype results
        random_sample = phenotype_results.sample( \
            n=n_connectome_genes_in, \
            replace=False, \
            random_state=iteration \
        )
            #we take as many rows, i.e., genes, as genes we have in the interest set of genes after removing those without hg38 coordinates
            #replace=False to avoid the same gene to fall within the same set two times
            #set the seed using the iteration number

        #calculate the median of "mean_effect" for the random sample
        random_median = random_sample["mean_effect"].median()

        #check if the random median is higher than the connectome median
        if random_median >= connectome_median:
            count_higher += 1

        #save the median value
        random_median_list.append(random_median)

    #calculate the empirical p-value
    empirical_p_value = count_higher / n_iterations
        #if 10 out 100 random set are above the interest set, that would be a p-value of 0.1, we would need in less than 5 out 100, i.e., 0.05.

    #calculate the 95CI for the median effect across random sets
    random_lower_95CI = pd.Series(random_median_list).quantile(0.025)
    random_upper_95CI = pd.Series(random_median_list).quantile(0.975)

    # Return the result for this phenotype
    return { \
        "pheno": pheno, \
        "random_lower_95CI": random_lower_95CI, \
        "random_upper_95CI": random_upper_95CI, \
        "connectome_median": connectome_median, \
        "empirical_p_value": empirical_p_value \
    }

print_text("prepare parallelization creating an instance of the function to be run specifically with the required number of iterations", header=4)
process_with_iterations = partial( \
    process_phenotype_randomization, \
    n_iterations=n_iterations \
)

print_text("Use multiprocessing Pool to process each phenotype in parallel", header=4)
with Pool(len(phenotype_names)) as pool:
    empirical_p_values = pool.map(process_with_iterations, phenotype_names)

print_text("convert the results into a DataFrame", header=4)
empirical_p_values_df = pd.DataFrame(empirical_p_values)
pd.set_option("display.max_columns", None)
    #ensure all columns are displayed when printing the DataFrame
print(empirical_p_values_df)
pd.reset_option("display.max_columns")
    #reset the display settings to default after printing

# endregion







print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/association_studies/
#chmod +x ./scripts/production_scripts/03f_bat_smt_analyses.py
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03f_bat_smt_analyses.py --pressure_name="bat" --connectome_percentile=10 --n_iterations=1000000 > ./03f_bat_smt_analyses_bat.out 2>&1
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03f_bat_smt_analyses.py --pressure_name="smt" --connectome_percentile=10 --n_iterations=1000000 > ./03f_bat_smt_analyses_smt.out 2>&1
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03f_bat_smt_analyses.py --pressure_name="metabolic" --connectome_percentile=10 --n_iterations=1000000 > ./03f_bat_smt_analyses_metabolic.out 2>&1
#grep -Ei 'error|false|fail' ./03f_bat_smt_analyses_bat.out
#grep -Ei 'error|false|fail' ./03f_bat_smt_analyses_smt.out
#grep -Ei 'error|false|fail' ./03f_bat_smt_analyses_metabolic.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.
