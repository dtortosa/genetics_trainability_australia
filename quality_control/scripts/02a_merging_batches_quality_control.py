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

#QC BEFORE OR AFTER MERGINIG?
    #Identifying and mitigating batch effects in whole genome sequencing data

#This script will merge the plink binary files of both batches and then check if there is any batch effect and then perform QCC



#FOR QUALITY CONTROL, YOU COULD USE MIGHIGGAN SERVER, WHICH ALREADY APPLIES MULTIPLE FILTERS like duplicates AND THEN ADD A FEW WITH PLINK, AUGUSTO DID THAT
    #https://www.mdpi.com/2073-4425/14/2/248
    #https://imputationserver.readthedocs.io/en/latest/pipeline/





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
