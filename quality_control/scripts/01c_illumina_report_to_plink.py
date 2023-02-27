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

#This script will merge the plink binary files of both batches and then check if there is any batch effect



#CHECK TUTORIALs so you are sure you are follwing the correct steps for filtering...



##IMPORTANT, wen merging different batches, check same strand
    #https://www.biostars.org/p/310290/


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
