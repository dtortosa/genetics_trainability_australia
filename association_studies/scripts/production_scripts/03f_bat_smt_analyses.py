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
    #take the 1000kb gene windows for all coding genes, liftover to hg38. If the USCS tool accepts genomic ranges, just use them as input, if not, split in two datasets the start and the end of the gene windows
    #for each phenotype (VO2, beep....), calculate the average (better than median because want influence of outliers within gene like in iHS, if a SNPs is veery important in a gene that should influence the info about the whole gene) effect size for the association of SNPs inside each gene
    #then, calculate 1000 random sets of genes, within each set, calculate the median effect of all genes inside the set and compare with the BAT set to obtain a distribution and empirical p-value (is association lower in BAT? LOOF BAT PAPER SCRIPTS FOR THIS). Here I want median because i do not want a gene outliser change things, I want the overall impact of BAT genes in general, not just a few genes.

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





gene_coords = pd.read_table(\
    "../../../../../postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt", \
    sep="\t", \
    header=0, \
    low_memory=False)
print(gene_coords)


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
    "lower_end_window_50kb", \
    "upper_end_window_50kb",  \
    "lower_end_window_100kb", \
    "upper_end_window_100kb",  \
    "lower_end_window_200kb", \
    "upper_end_window_200kb", \
    "lower_end_window_500kb", \
    "upper_end_window_500kb", \
    "lower_end_window_1000kb", \
    "upper_end_window_1000kb"]]
        #1000kb windows better, see our papers


    #assuming windows are 1-based as in conserved elements density script we calculate the average density inside the windows after summing 1 to the start of the uscs data
    #this means we can use the format "chr4:100,001-100,001" for liftover as this makes it assume 1-based coordinates
        #"If you submit data to the browser in position format (chr#:##-##), the browser assumes this information is 1-based."
    #see "/home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/method_deep/scripts/genomic_features_calculations/gene_windows/cons_elements_density/cons_elements_density_v3.R"


if gene_coords_no_duplicated_subset["gene_id"].isna().sum() !=0:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE SUBSET OF THE GENE COORDINATE FILE")

gene_coords_no_duplicated_subset = gene_coords_no_duplicated_subset.dropna(
    subset=["lower_end_window_1000kb", "upper_end_window_1000kb"]
)

gene_coords_no_duplicated_subset["chromosome_name"] = "chr" + gene_coords_no_duplicated_subset["chromosome_name"].astype(str)


gene_coords_no_duplicated_subset["lower_end_window_1000kb"] = gene_coords_no_duplicated_subset["lower_end_window_1000kb"] - 1

gene_coords_no_duplicated_subset["lower_end_window_1000kb"] = gene_coords_no_duplicated_subset["lower_end_window_1000kb"].astype(int)
gene_coords_no_duplicated_subset["upper_end_window_1000kb"] = gene_coords_no_duplicated_subset["upper_end_window_1000kb"].astype(int)




# Print the updated DataFrame
print(gene_coords_no_duplicated_subset[["chromosome_name", "lower_end_window_1000kb", "upper_end_window_1000kb", "gene_id", "hgnc_symbol"]])

run_bash(" \
    mkdir \
        -p \
        ./data/bat_smt_analyses/ \
")

gene_coords_no_duplicated_subset[["chromosome_name", "lower_end_window_1000kb", "upper_end_window_1000kb", "gene_id"]].to_csv( \
    "./data/bat_smt_analyses/gene_windows_hg19.tsv", \
    sep="\t", \
    header=False, \
    index=False \
)


gene_windows_hg38 = pd.read_csv(\
    "./data/bat_smt_analyses/gene_windows_hg38.bed", \
    sep="\t", \
    header=None, \
    low_memory=False \
)


gene_windows_hg38.columns = ["chromosome_name", "lower_end_window_1000kb", "upper_end_window_1000kb", "gene_id"]

gene_windows_hg38["lower_end_window_1000kb"] = gene_windows_hg38["lower_end_window_1000kb"] + 1







# Function to process a single phenotype
#phenotype_name = "vo2_change"
def process_phenotype(phenotype_name):

    run_bash(" \
        gunzip \
            --keep \
            --force \
            ./results/final_results/analysis_full_data/" + phenotype_name + "/" + phenotype_name + "_small_set_predictors_set_linear_raw.assoc.gz \
    ")

    assoc_results = pd.read_csv( \
        "./results/final_results/analysis_full_data/" + phenotype_name + "/" + phenotype_name + "_small_set_predictors_set_linear_raw.assoc", \
        sep="\t", \
        header=0, \
        low_memory=False, \
    )

    assoc_results["Chromosome"] = "chr" + assoc_results["Chromosome"].astype(str)

    results = []

    #Group genes by chromosome for efficient processing
    #chromosome, gene_data = [i for i in gene_windows_hg38.groupby("chromosome_name") ][0]
    for chromosome, gene_data in gene_windows_hg38.groupby("chromosome_name"):
        
        

        # Filter assoc_results for the current chromosome
        snps_in_chromosome = assoc_results[assoc_results["Chromosome"] == chromosome]

        # Iterate over each gene in the chromosome
        #_, gene_row = [i for i in gene_data.iterrows()][0]
        for _, gene_row in gene_data.iterrows():
            
            # Extract the genomic window for the gene
            lower_bound = gene_row["lower_end_window_1000kb"]
            upper_bound = gene_row["upper_end_window_1000kb"]

            # Filter SNPs within the genomic window
            snps_in_window = snps_in_chromosome[
                (snps_in_chromosome["Basepair"] >= lower_bound) &
                (snps_in_chromosome["Basepair"] <= upper_bound)
            ]

            effect_sizes_abs = snps_in_window["Effect"].abs()
                #do not care of the sense (ref or alt) of the effect, just care about the absolute value

            # Calculate mean and SD of the "Effects" column
            mean_effect = effect_sizes_abs.mean() if not snps_in_window.empty else np.nan
            sd_effect = effect_sizes_abs.std() if not snps_in_window.empty else np.nan
                #prefer mean over median because we want to take into account the outliers, as in iHS, if a SNP is very important in a gene that should influence the info about the whole gene

            # Create a result dictionary for the gene
            result = {
                "phenotype": phenotype_name,
                "chromosome_name": chromosome,
                "gene_id": gene_row["gene_id"],
                "mean_effect": mean_effect,
                "sd_effect": sd_effect
            }
            results.append(result)

    if len(results) != gene_windows_hg38.shape[0]:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF GENES IN THE RESULTS")

    run_bash(" \
        rm ./results/final_results/analysis_full_data/" + phenotype_name + "/" + phenotype_name + "_small_set_predictors_set_linear_raw.assoc \
    ")

    return results

# Get the list of phenotypes
phenotype_names = ["distance_change", "vo2_change", "beep_change", "weight_change"]

# Use multiprocessing Pool to process each phenotype in parallel
with Pool(len(phenotype_names)) as pool:
    all_results = pool.map(process_phenotype, phenotype_names)

# Flatten the list of lists into a single list
flattened_results = [item for sublist in all_results for item in sublist]

# Convert the results into a DataFrame
results_df = pd.DataFrame(flattened_results)

# Print or save the results
print(results_df)




results_df_symbol = results_df.merge(gene_coords_no_duplicated_subset[["gene_id", "hgnc_symbol"]], how="left", on="gene_id")





print_text("first extract the connectome of the interest gene", header=4)
if(pressure_name == "bat"):
    core_gene = "UCP1"
elif(pressure_name == "smt"):
    core_gene = "SLN"
elif(pressure_name == "metabolic"):
    core_gene = "NA"




if pressure_name in ["bat", "smt"]:
    run_bash(" \
        mkdir \
            -p \
            ./data/bat_smt_analyses/" + pressure_name + "_distance/; \
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


    print_text("check we have selected the correct connectome, i.e., the one with the correct core gene", header=4)
    print(pressure_conn["Target"].iloc[0] == core_gene)


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
        print(pressure_conn.equals(old_ucp1_conn))
        print(pressure_conn.equals(ucp1_2020))
            #both files are the same, suggesting that the warning obtained when zipping is not a problem
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
        print(pressure_conn \
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
                    ignore_index=True)))
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



    print_text("Percentile 10% :  select the interest genes genes", header=3)
    selected_connectome_genes = pressure_conn \
        .loc[ \
            pressure_conn["Target_in_source_P-value(percentile)"] < (connectome_percentile/100), \
            "Target"]
    print(selected_connectome_genes)

elif(pressure_name == "metabolic"):
    metabolic_genes = pd.read_csv("../../../../../postdoc_enard_lab/projects/ancient_selection_dating/data/metabolic_genes/metabolic_gene_list.txt.gz", sep="\t", header=0, low_memory=False)

    metabolic_genes.columns=["hgnc_symbol", "gene_id"]

    selected_connectome_genes = metabolic_genes["hgnc_symbol"]


selected_connectome_genes


n_connectome_genes_in = sum(np.isin(results_df_symbol["hgnc_symbol"].unique(), selected_connectome_genes))



if (len(selected_connectome_genes)-n_connectome_genes_in)/len(selected_connectome_genes)*100 >15:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE CONNECTOME GENES. MORE THAN 10% OF THEM ARE NOT IN THE RESULTS")
    #more lost due to liftover?



connectome_reference = []

#pheno=phenotype_names[0]
for pheno in phenotype_names:
    
    results_df_symbol_subset = results_df_symbol.loc[ (results_df_symbol["hgnc_symbol"].isin(selected_connectome_genes)) & (results_df_symbol["phenotype"]==pheno), :]

    connectome_reference.append({"pheno": pheno, "median_avg_effect_connectome": results_df_symbol_subset["mean_effect"].median(), "median_std_effect_connectome": results_df_symbol_subset["sd_effect"].median()})


connectome_reference_df=pd.DataFrame(connectome_reference)
connectome_reference_df



from functools import partial



#pheno=phenotype_names[0]

# Function to process a single phenotype
#pheno="weight_change"
def process_phenotype_randomization(pheno, n_iterations):
    # Subset the results for the current phenotype
    phenotype_results = results_df_symbol[(results_df_symbol["phenotype"] == pheno) & (~ results_df_symbol["hgnc_symbol"].isin(selected_connectome_genes))]

    # Get the median average effect for the connectome genes for this phenotype
    connectome_median = connectome_reference_df.loc[
        connectome_reference_df["pheno"] == pheno, "median_avg_effect_connectome"
    ].values[0]

    # Initialize a counter for how many random medians are higher than the connectome median
    count_higher = 0

    random_median_list = []


    # Perform 100 randomizations
    for iteration in range(n_iterations):
        # Randomly sample 167 rows from the phenotype results
        random_sample = phenotype_results.sample(n=n_connectome_genes_in, replace=False, random_state=iteration)

        # Calculate the median of "mean_effect" for the random sample
        random_median = random_sample["mean_effect"].median()

        # Check if the random median is higher than the connectome median
        if random_median >= connectome_median:
            count_higher += 1

        #save
        random_median_list.append(random_median)

    # Calculate the empirical p-value
    empirical_p_value = count_higher / n_iterations

    random_lower_95CI = pd.Series(random_median_list).quantile(0.025)
    random_upper_95CI = pd.Series(random_median_list).quantile(0.975)

    # Return the result for this phenotype
    return {"pheno": pheno, "random_lower_95CI": random_lower_95CI, "random_upper_95CI": random_upper_95CI, "connectome_median": connectome_median, "empirical_p_value": empirical_p_value}


process_with_iterations = partial(process_phenotype_randomization, n_iterations=n_iterations)


# Use multiprocessing Pool to process each phenotype in parallel
with Pool(len(phenotype_names)) as pool:
    empirical_p_values = pool.map(process_with_iterations, phenotype_names)

# Convert the results into a DataFrame
empirical_p_values_df = pd.DataFrame(empirical_p_values)

# Print or save the results
print(empirical_p_values_df)


