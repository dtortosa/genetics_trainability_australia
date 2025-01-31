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



#######################################
######## PHENOTYPE PREPARATION ########
#######################################

#This script will prepare the phenotypes for the association analyses. 

#This and previous scripts are based on the following tutorials:
    #1. Quality Control Procedures for Genome-Wide Association Studies
        #https://github.com/RitchieLab/GWAS-QC
        #https://drive.google.com/file/d/1kxV3j_qCF_XMX47575frXhVwRMhzjqju/view
    #2. Genome-wide association studies
        #https://www.nature.com/articles/s43586-021-00056-9
    #3. Data Management and Summary Statistics with PLINK
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
    #4. Genomics Boot Camp
        #https://genomicsbootcamp.github.io/book/
    #5. Tutorial: a guide to performing polygenic risk score analyses
        #https://www.nature.com/articles/s41596-020-0353-1
        #https://choishingwan.github.io/PRS-Tutorial/
    #6. Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
    #7. Omics Data Preprocessing for Machine Learning: A Case Study in Childhood Obesity
        #https://www.mdpi.com/2073-4425/14/2/248
    #8. LDAK
        #https://dougspeed.com/
        #https://www.nature.com/articles/s41467-021-24485-y



#######################################
# region INITIAL ANOTATIONS AND STEPS #
#######################################

###########
# imports #
###########

import pandas as pd
import numpy as np



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
    else:
        #print the standard error and stop
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

#test it
print_text("check behaviour run_bash", header=1)
print_text("see working directory", header=2)
run_bash("pwd")
print_text("list files/folders there", header=2)
run_bash("ls")



######################
# folder preparation #
######################
run_bash(" \
    mkdir -p ./data/pheno_data; \
    mkdir -p ./data/plink_inputs; \
    mkdir -p ./results/final_results/distribution_hist; \
    ls -l")

# endregion






#####################################################################
# region process original excel file with phenotypes and covariates #
#####################################################################
print_text("combine the excel file, FAM file and the PCAs", header=1)
print_text("load the excel file", header=2)
#This include reported sex and VO2 max data. I have checked that the data is the same directly reading from excel than converting to csv
run_bash(" \
    cp \
        ../quality_control/data/pheno_data/'Combat gene DNA GWAS 23062022_v2.xlsx' \
        ./data/pheno_data/pheno_data.xlsx; \
    echo 'pheno_data.xlsx comes from ../quality_control/data/pheno_data/'Combat gene DNA GWAS 23062022_v2.xlsx', which is the second raw excel I got from D. Bishop after they add beep distance' > ./data/pheno_data/README.txt; \
    ls -l ./data/pheno_data/")
pheno_data = pd.read_excel(
    "./data/pheno_data/pheno_data.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data)


print_text("do some checks about the variables", header=2)
print_text("Check that the dtype of the new beep distance columns is float64", header=3)
if (pheno_data["Week 1 Distance (m)"].dtype == "float64") & (pheno_data["Week 8 Distance (m)"].dtype == "float64"):
    print("YES! GOOD TO GO!!")
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE DTYPE OF THE NEW BEEP DISTANCE COLUMNS")

print_text("Check that the new version of the pheno_data is exactly the same than the previous one if we remove the new beep distance variables", header=3)
print_text("load the original phenotype data", header=4)
pheno_data_no_beep_distance = pd.read_excel(
    "../quality_control/data/pheno_data/combact gene DNA GWAS 23062022.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data_no_beep_distance)

print_text("make the check", header=4)
check_pheno_versions = \
    pheno_data \
        .loc[:, ~pheno_data.columns.isin(["Week 1 Distance (m)", "Week 8 Distance (m)"])] \
        .equals(pheno_data_no_beep_distance)
if check_pheno_versions:
    print("YES! GOOD TO GO!!")
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NEW VERSION OF THE PHENO DATA")

print_text("Week 8 beep test for one of the samples is 11.1O, i.e., letter O instead number 0", header=3)
print_text("see the problem", header=4)
index_problematic_sample = np.where(pheno_data["Week 8 beep test"] == "11.1O")[0][0]
index_problematic_column = np.where(pheno_data.columns == "Week 8 beep test")[0][0]
print(pheno_data.iloc[index_problematic_sample, index_problematic_column])
print(pheno_data.iloc[index_problematic_sample,:])

print_text("change 11.1O for 11.10", header=4)
pheno_data.iloc[index_problematic_sample, index_problematic_column] = 11.1
print(pheno_data.iloc[index_problematic_sample,:])

print_text("convert the dtype of this column from string to float", header=4)
print("old dtype: " + str(pheno_data["Week 8 beep test"].dtype))
pheno_data["Week 8 beep test"] = pheno_data["Week 8 beep test"].astype(dtype="float64")
print("new dtype: " + str(pheno_data["Week 8 beep test"].dtype))
if(pheno_data["Week 8 beep test"].dtype == "float64"):
    print("GOOD TO GO: we have correctly change the dtype of the column with beep test data week 8")
else:
    raise ValueError("ERROR: FALSE! The dtype of the beep test week 8 column is not correct")
print("Sample 8244FGNJ has NA for Week 8 Distance (m). I detected a typo in the Week 8 beep test for that same sample, having 11.1O instead of 11.10, i.e., there is letter O instead of the zero number")


print_text("load fam file", header=2)
print_text("load the last fam file", header=3)
fam_file = pd.read_csv( \
    "../quality_control/data/genetic_data/quality_control/21_post_imputation_qc/03_third_qc_step/merged_3_geno.fam", \
    sep=" ", \
    header=None, \
    low_memory=False, \
    names=["family_id", "AGRF code", "ID_father", "ID_mother", "sex_code", "phenotype_value"])
    #add before hand the column names. The order is indicated in plink docs
        #Family ID ('FID')
        #Within-family ID ('IID'; cannot be '0')
        #Within-family ID of father ('0' if father isn't in dataset)
        #Within-family ID of mother ('0' if mother isn't in dataset)
        #Sex code ('1' = male, '2' = female, '0' = unknown)
        #Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            #https://www.cog-genomics.org/plink/1.9/formats#fam
print(fam_file)

print_text("check the fam file do not have any duplicated IDs", header=3)
#the plink fam file is going to be the reference because this is the result of all QC we applied, so any sample no present here, should NOT be considered in futher analyses
print_text("see duplicated IDs in pheno data", header=4)
dup_ids_pheno_data=pheno_data.loc[pheno_data["AGRF code"].duplicated(), "AGRF code"]
print(dup_ids_pheno_data)

print_text("check these samples are not included in the plink FAM file", header=4)
if( \
    (fam_file.loc[fam_file["AGRF code"].duplicated(),:].shape[0] != 0) | \
    (fam_file.loc[fam_file["AGRF code"].isin(dup_ids_pheno_data), :].shape[0] != 0) \
):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH DUPLICATED IDS")
else:
    print("GOOD TO GO: we do not have duplicated IDs in the plink FAM file, SO WE CAN MERGE PHENO DATA AND FAM FILE USING ONLY THE AGRF CODE")


print_text("merge and the excel file with the FAM file", header=2)
print_text("do the merge", header=3)
merged_data_raw_1 = pheno_data.merge( \
    fam_file, \
    on="AGRF code", \
    how="right" \
)
    #how="right":
        #use only keys from right frame, similar to a SQL right outer join
        #we use right, i.e., fam file (genetic data), for the keys because we only want to retain those samples that remained after the filters done in plink. We want to take advantage of the filtering work already done.
#check
if(merged_data_raw_1.shape[0]!=1203):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE MERGE")

#make a deep copy to leave the raw merged data intact
merged_data_raw_2 = merged_data_raw_1.copy(deep=True)
    #deep=True
        #Modifications to the data or indices of the copy will not be reflected in the original object (see notes below).
print(merged_data_raw_2)
print(merged_data_raw_2.describe().T)

print_text("total sum of NaNs or zero cases per column", header=4)
nan_zeros_1 = merged_data_raw_2.apply(lambda x: sum((x==0) | (x.isna())), axis=0)
print(nan_zeros_1)

print_text("check and remove cases with zero in phenotypes", header=3)
print_text("sum the number of Zeros per column", header=4)
print(merged_data_raw_2.apply(lambda x: sum(x==0), axis=0))
    #note that the number of zeros does NOT have to match that of the pheno_data because we have already filtered many samples in plink
    #for example, ["8502JGMJ", "9601AVJM"] have zero for "Week 1 Pred VO2max", but they are not included in the fam file
        #merged_data.loc[merged_data["AGRF code"].isin(["8502JGMJ", "9601AVJM"]),:]
        #fam_file.loc[fam_file["AGRF code"].isin(["8502JGMJ", "9601AVJM"]),:]

print_text("Weight: Given that weight is going to be used for the rest of variables as covariate, we will set as missing the weight of these samples. Importantly, a weight of zero does not make sense.", header=4)
merged_data_raw_2.loc[merged_data_raw_2["Week 1 Body Mass"]==0, "Week 1 Body Mass"] = np.nan
merged_data_raw_2.loc[merged_data_raw_2["Week 8 Body Mass"]==0, "Week 8 Body Mass"] = np.nan

print_text("VO2 max week 1: In the case of VO2 max, this is probably caused by the equation. It is indeed a very strange value because the beep test data is ok and the VO2 max of the 8th week is also ok. We are going to set this as nan", header=4)
print(merged_data_raw_2.loc[merged_data_raw_2["Week 1 Pred VO2max"]==0, ["Week 1 Beep test", "Week 8 beep test", "Week 1 Pred VO2max", "Week 8 Pred VO2max",]])
merged_data_raw_2.loc[merged_data_raw_2["Week 1 Pred VO2max"]==0, "Week 1 Pred VO2max"] = np.nan

print_text("Sex_code: We initially had missing for some samples, but these have been removed during the QC, also 2397LDJA was updated to 2399LDJA, so we should not have any zero for sex", header=4)
if (merged_data_raw_2.loc[merged_data_raw_2["sex_code"]==0, :].shape[0]!=0):
    raise ValueError("ERROR: FALSE! WE HAVE EMPTY SEX VALUES")

print_text("Parents ID: We do not have any parent in the cohort, so they all should be zero", header=4)
if(merged_data_raw_2.loc[(merged_data_raw_2["ID_father"]!=0) | (merged_data_raw_2["ID_mother"]!=0), :].shape[0] != 0):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE PARENTS IDs, THEY SHOULD BE ZERO")

print_text("sum again the number of Zeros per column", header=4)
sum_zeros_columns = merged_data_raw_2.apply(lambda x: sum(x==0), axis=0)
print(sum_zeros_columns)
check_zeros = np.array_equal( \
    sum_zeros_columns.loc[sum_zeros_columns!=0].index.to_numpy(), \
    np.array(['ID_father', 'ID_mother']) \
)
if (not check_zeros):
    raise ValueError("We have zeros in columns that should not have zeros, i.e., not only in ID_father and ID_mother")


print_text("deal with 2397LDJA/2399LDJA", header=3)
print("This is the misslabeled sample, it has one ID in genetic and other different ID in pheno. As the postdoc of David said: I think it is a mislabelling of the last digit of the number (the labelling was very hard to read on some of the blood samples). So, I think 2397LDJA; ILGSA24-17303 is 2399LDJA in the excel file. We already changed the ID during the QC, so only 2399LDJA should be present")
if ( \
    merged_data_raw_2.loc[merged_data_raw_2["AGRF code"] == "2397LDJA"].shape[0]!=0 | \
    merged_data_raw_2.loc[merged_data_raw_2["AGRF code"] == "2399LDJA"].shape[0]!=1 \
):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH 2397LDJA/2399LDJA")


print_text("bind PCA axes", header=2)
print_text("load axes", header=3)
smartpca_eigenvectors = pd.read_csv(\
    "../quality_control/data/genetic_data/quality_control/14_pop_strat/01_pca/eigen_out/loop_maf_missing_2_ldprunned_autosome_pca_tab.eigs", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(smartpca_eigenvectors)
#check we have the correct number of columns
if(smartpca_eigenvectors.shape[1]-1!=20):
    raise ValueError("ERROR! FALSE! We are not considering all the interesting axes in smartpca")

print_text("split IDs", header=3)
smartpca_eigenvectors[["family_id", "AGRF code"]] = smartpca_eigenvectors.iloc[:,0].str.split(":", expand=True)
    #splits the column at each ":" and expands the result into separate columns.
    #df[['new_col1', 'new_col2']] assigns the split columns to new columns named new_col1 and new_col2.
print(smartpca_eigenvectors)
#check
if(not smartpca_eigenvectors[0].equals(smartpca_eigenvectors["family_id"] + ":" + smartpca_eigenvectors["AGRF code"])):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM SPLITTING THE IDS IN THE PCA FILE")

print_text("change 2397LDJA by 2399LDJA. 2399LDJA is the correct ID showed in the pheno data, and we changed in the genetic data but after the PCA, so the error remains in the PCA data", header=3)
if(
    (smartpca_eigenvectors.loc[smartpca_eigenvectors["AGRF code"] == "2397LDJA", :].shape[0]!=1) |
    (smartpca_eigenvectors.loc[smartpca_eigenvectors["AGRF code"] == "2399LDJA", :].shape[0]!=0)
):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH 2397LDJA/2399LDJA")
else:
    smartpca_eigenvectors.loc[smartpca_eigenvectors["AGRF code"] == "2397LDJA", "AGRF code"] = "2399LDJA"

print_text("select the PCAs we are interested in", header=3)
#We are going to use only the 2 first PCAs:
    #There is a clear reduction of explained variability from 1 to and from 2 to 3, but from 3 there are small differences in explained variability. 
    #Admixture considers that we only have ONE continental ancestry in our dataset. Moreover, the PCA plots show that if there is an axis that is spreading a bit the samples is the first one, the rest show very cloudy patterns, almost perfect.
    #Yes, according to smartpca, we have 6 significant PCAs, but if we used that criteria we would include PCA5 which show a little bit of pattern, the loadings tend to increase a lot in one positon of chromosome 5 for this PCA. We may see a bit of that pattern in the same region for PCA 4 and 3, so I would not include any of these PCAs and just focus on the two most important.
    #We are going to stick to 1 and 2.
smartpca_eigenvectors_subset = smartpca_eigenvectors[["family_id", "AGRF code", 1, 2]]
print(smartpca_eigenvectors_subset)
    #I have plotted PCA1 against PC2 and I get the same plot we saw in "eigen_vectors_pairwise.png" for PC1 agaist PC2. So, as expected, the first two columns correspond with the two first axes.
#change the column names of the PCAs
smartpca_eigenvectors_subset = smartpca_eigenvectors_subset.rename({1: "PCA1", 2: "PCA2"}, inplace=False, axis=1)
print(smartpca_eigenvectors_subset)


print_text("merge the PCA with the rest of phenotypes", header=2)
print_text("do the merge", header=3)
merged_data=smartpca_eigenvectors_subset.merge( \
    merged_data_raw_2, \
    on=["family_id", "AGRF code"], \
    how="right" \
)
    #how="right":
        #use only keys from right frame, similar to a SQL right outer join
        #we use right, i.e., fam file (genetic data), for the keys because we only want to retain those samples that remained after the filters done in plink. We want to take advantage of the filtering work already done.
        #Remember that we removed additional samples after the PCA, mainly because of sex inconsistences...
#check
if(merged_data.shape[0]!=1203):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE MERGE")



print_text("process the merged file", header=1)
print_text("create new variables for the change before and after", header=2)
print_text("make the operations", header=3)
merged_data["weight_change"] = merged_data["Week 8 Body Mass"]-merged_data["Week 1 Body Mass"]
merged_data["beep_change"] = merged_data["Week 8 beep test"]-merged_data["Week 1 Beep test"]
merged_data["distance_change"] = merged_data["Week 8 Distance (m)"]-merged_data["Week 1 Distance (m)"]
merged_data["vo2_change"] = merged_data["Week 8 Pred VO2max"]-merged_data["Week 1 Pred VO2max"]
    #NA remains as NA

print_text("check we have the correct NANs", header=3)
print(merged_data.isna().apply(lambda x: sum(x), axis=0))
print("Weight change")
print( \
    sum(merged_data["weight_change"].isna()) == \
    sum( \
        (merged_data["Week 1 Body Mass"].isna()) | \
        (merged_data["Week 8 Body Mass"].isna())))
print("beep change")
print( \
    sum(merged_data["beep_change"].isna()) == \
    sum( \
        (merged_data["Week 1 Beep test"].isna()) | \
        (merged_data["Week 8 beep test"].isna())))
print("distance change")
print( \
    sum(merged_data["distance_change"].isna()) == \
    sum( \
        (merged_data["Week 1 Distance (m)"].isna()) | \
        (merged_data["Week 8 Distance (m)"].isna())))
print("VO2 max change")
print( \
    sum(merged_data["vo2_change"].isna()) == \
    sum( \
        (merged_data["Week 1 Pred VO2max"].isna()) | \
        (merged_data["Week 8 Pred VO2max"].isna())))
    #the total number of NAN for the change variables should be the same than the sum of cases that are NAN for week 1 or week 2 of the corresponding variable.

print_text("Check we only have NA", header=3)
print("Remember that LDAK consider missing pheno as NA, not -9, so we do not need to change the missing value. We do not have problems with the negative values because now -9 is not going to be NA, so a person losing 9Kg is going to be considered a value, not missing")
#Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). Binary phenotypes should only take values 0 (control), 1 (case) or NA (missing). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype).
    #http://dougspeed.com/phenotypes-and-covariates/

print_text("check if we have -9 cases aside from the 'change variables'", header=4)
minus_nine_cases=merged_data.apply(lambda x: x==-9, axis=1).sum()
print(minus_nine_cases)
    #sum number of -9 cases per column
#check if the columns with -9 are only "phenotype_value" coming from the FAM file and "weight_change", which is suppose to have zero values. Also check whether phenotype_value has only -9
check_minus_nine = ( \
    np.array_equal( \
        minus_nine_cases.loc[minus_nine_cases>0].index.to_numpy(), \
        np.array(["phenotype_value", "weight_change"]) \
    ) & \
    (bool(minus_nine_cases["phenotype_value"]==merged_data.shape[0])) \
)
if(check_minus_nine):
    #if that is the case, then remove the phenotype_value column because this is basically a column filled with NAs
    merged_data = merged_data.drop("phenotype_value", axis=1, inplace=False)
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH -9 VALUES")
print(merged_data)  
#check that after cleaning the only column with -9 is the weight_change
minus_nine_cases_after = merged_data.apply(lambda x: x==-9, axis=1).sum()
if( \
    np.array_equal( \
        minus_nine_cases_after[minus_nine_cases_after>0].index.to_numpy(), \
        np.array(["weight_change"]) \
    )
):
    print("GOOD TO GO: we only have -9 cases in the 'change variables'")
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH -9 VALUES")

print_text("We do NOT have to check for zero! cases", header=4)
    #We have calculated the difference of weight and cardiorespiratory fitness so a person that did not change his/her weight would have zero, totally ok. Also we can have 0 values for the PCA, if you check the PCA plots, both PCA1 and PC2 have values around zero, so ok.
    #We checked for zeros in the pheno file because before calculating any "change variable", we had zeros for weight and VO2, which is not possible, so we converted them to NaN. We should also have NAs for mother/father IDs
zero_cases_after = merged_data.apply(lambda x: x==0, axis=1).sum()
#check that the columns with zero are all the ones previously explained. Also father and mother IDs are all zero.
if( \
    bool(zero_cases_after[zero_cases_after>0].index.isin(["PCA1", "PCA2", "ID_father", "ID_mother", "weight_change", "beep_change", "distance_change", "vo2_change"]).all()) & \
    bool(zero_cases_after["ID_father"]==merged_data.shape[0]) & \
    bool(zero_cases_after["ID_mother"]==merged_data.shape[0]) \
):
    print("GOOD TO GO: we only have 0 cases in the 'change' variables, PCA and father/mother IDs")
    merged_data = merged_data.drop(["ID_father", "ID_mother"], axis=1, inplace=False)
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH 0 VALUES")


print_text("perform some checks", header=2)
print_text("check we have the correct index in order to filter by index and iloc, i.e., no gaps in the indexes", header=3)
print(np.array_equal(merged_data.index, range(0,merged_data.shape[0])))

print_text("check the problematic samples are NOT included", header=3)
print(sum(merged_data["AGRF code"].isin(["1100JHJM", "1200JPJM", "7800AGSO", "2397LDJA"]))==0)

print_text("check we do not have duplicated samples and no sample with NaN", header=3)
print(sum(merged_data["AGRF code"].duplicated(keep=False)) == 0)
    #keep=False
        #Mark all duplicates as True
print(sum(merged_data["AGRF code"].isna()) == 0)

print_text("check we do have two sex columns", header=3)
if(bool(merged_data.columns.isin(["Gender", "sex_code"]).sum()==2)):
    print("GOOD TO GO: we have the two expected Sex columns. Now update the names of gender to self_reported_sex to avoid confusion with plink sex")
    merged_data = merged_data.rename({"Gender": "self_reported_sex"}, axis=1, inplace=False)
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE SEX COLUMNS")


###por aquii

print_text("final exploration of the data", header=2)
print_text("count missing per column", header=3)
count_missing = merged_data.isna().sum()
print(count_missing)

print_text("check we have the correct number of missing based on the raw merged data", header=3)
print(count_missing["self_reported_sex"] == sum(merged_data_raw_2["Gender"].isna()))
print(count_missing["Age"] == sum(merged_data_raw_2["Age"].isna()))
#selected_phenotype="Week 8 Distance (m)"
for selected_phenotype in ["Week 1 Body Mass", "Week 8 Body Mass", "Week 1 Beep test", "Week 8 beep test", "Week 1 Distance (m)", "Week 8 Distance (m)", "Week 1 Pred VO2max", "Week 8 Pred VO2max"]:
    print(selected_phenotype)
    print(count_missing[selected_phenotype] == merged_data_raw_2[selected_phenotype].isna().sum())

print("In the case of family/parents Ids and sex code, we should not have any missing, but 0 for unknown, so we do not need to subtract 1 for the mislabelled sample (2397LDJA/2399LDJA)")   
print(count_missing["family_id"] == sum(merged_data_raw_2["family_id"].isna()))
print(count_missing["sex_code"] == sum(merged_data_raw_2["sex_code"].isna()))
print("check also that the family/parents IDs and the sex code are exactly the same")
print(merged_data_raw_2[["family_id", "sex_code"]].equals(merged_data[["family_id", "sex_code"]]))

print_text("check that the number to missing is equal to the cases of NaN OR 0 in the initial merged dataset", header=4)
print( \
    nan_zeros_1 \
        [["Age", "Week 1 Body Mass", "Week 8 Body Mass", "Week 1 Beep test", "Week 8 beep test", "Week 1 Distance (m)", "Week 8 Distance (m)", "Week 1 Pred VO2max", "Week 8 Pred VO2max"]] == \
    (count_missing \
        [["Age", "Week 1 Body Mass", "Week 8 Body Mass", "Week 1 Beep test", "Week 8 beep test", "Week 1 Distance (m)", "Week 8 Distance (m)", "Week 1 Pred VO2max", "Week 8 Pred VO2max"]]))
    #nan_zeros_1 was obtained at the begining after merging using two conditions, isna() and ==0
    #Gender is not included because we have change it its name, so ti does not work

print_text("see the final data", header=4)
print(merged_data)

print_text("see summary statistics of non_pheno columns", header=4)
print(merged_data.describe().T)

print_text("see summary statistics of pheno columns", header=4)
#selected_phenotype="Week 8 Distance (m)"
for selected_phenotype in ["PCA1", "PCA2", "Age", "Week 1 Body Mass", "Week 8 Body Mass", "Week 1 Beep test", "Week 8 beep test", "Week 1 Distance (m)", "Week 8 Distance (m)", "Week 1 Pred VO2max", "Week 8 Pred VO2max", "weight_change", "beep_change", "distance_change", "vo2_change"]:
    print("\n" + selected_phenotype)
    non_missing_subset = merged_data.loc[~merged_data[selected_phenotype].isna(), selected_phenotype]
        #select values of the selected phenotype that are NOT missing
    print("mean: " + str(np.mean(non_missing_subset)))
    print("percentile 2.5: " + str(np.percentile(non_missing_subset, 2.5)))
    print("percentile 25: " + str(np.percentile(non_missing_subset, 25)))
    print("percentile 50: " + str(np.percentile(non_missing_subset, 50)))
    print("percentile 75: " + str(np.percentile(non_missing_subset, 75)))
    print("percentile 97.5: " + str(np.percentile(non_missing_subset, 97.5)))
print("The values of the phenotypes make sense: ")
print("A tendency for lower values in the highest percentiles of weight at week 8, i.e., more decrease for people starting with more weight. There is a tendency of increased weight for lower percentile. The median change of weight is 0, meaning that half of sample did not lose weight. The distribution is centered around zero.")
print("A tendency of increase of beep test, distance and VO2 max in all percentiles at week 8. The higher percentiles of changes tend to be more separated from zero than the lower percentile, i.e., the distribution is displaced to positive values, meaning increase of cardiorespiratory fitness.")

print_text("save the data", header=3)
merged_data.to_csv("./data/pheno_data/pheno_data_cleaned.tsv",
    sep="\t",
    header=True,
    index=False, 
    na_rep="NA")
    #na_rep="NA"
        #This parameter specifies the string representation of NaN values in the output file. "NA" means that any NaN values in the DataFrame will be written as "NA" in the CSV file.

print_text("see some columns of the file", header=3)
run_bash(" \
    cd ./data/pheno_data/; \
    awk \
        'BEGIN{FS=\"\t\"}{if(NR<=10){print $1, $2, $3, $4}}' \
        pheno_data_cleaned.tsv; \
    ls -l")


print_text("create objects with phenotypes and covariates", header=2)
print_text("create multi pheno files", header=3)
print("According to plink man: \
    '--pheno causes phenotype values to be read from the 3rd column of the specified space- or tab-delimited file, instead of the .fam or .ped file. The FIRST AND SECOND COLUMNS OF THAT FILE MUST CONTAIN FAMILY AND WITHIN-FAMILY IDS, respectively. \
    In combination with --pheno, --mpheno lets you use the (n+2)th column instead of the 3rd column, while --PHENO-NAME LETS YOU SELECT A COLUMN BY TITLE. (In order to use --pheno-name, there must be a header row with first two entries 'FID' and 'IID'.) The new --pheno-merge flag tells PLINK to use the phenotype value in the .fam/.ped file when no value is present in the --pheno file; without it, the phenotype is always treated as missing in this case.'")
    #https://www.cog-genomics.org/plink/1.9/input
print("Therefore, we need family and within family ID and then the phenotypes in a tab-delimited file with a header, including FID and IID as the two first columns.")
multi_pheno_file = merged_data[["family_id", "AGRF code", "weight_change", "beep_change", "distance_change", "vo2_change"]]
multi_pheno_file.columns=["FID", "IID", "weight_change", "beep_change", "distance_change", "vo2_change"]
    #select the interest columns and set the required column names
multi_pheno_file.to_csv( \
    "./data/plink_inputs/pheno_file.tsv", \
    sep="\t", \
    index=False, \
    header=True)
run_bash(" \
    cd ./data/plink_inputs/; \
    awk \
        'BEGIN{FS=\"\t\"}{if(NR<=10){print $0}}' \
        pheno_file.tsv")


print_text("create pheno/covariate files", header=3)
print("Plink man: ' \
    --covar designates the file to load covariates from. The file format is the same as for --pheno (optional header line, FID and IID in first two columns, covariates in remaining columns). BY DEFAULT, THE MAIN PHENOTYPE IS SET TO MISSING IF ANY COVARIATE IS MISSING; you can disable this with the 'keep-pheno-on-missing-cov' modifier. \
    --covar-name lets you specify a subset of covariates to load, by column name; separate multiple column names with spaces or commas, and use dashes to designate ranges. (Spaces are not permitted immediately before or after a range-denoting dash.) --covar-number lets you use column numbers instead.'")
    #https://www.cog-genomics.org/plink/1.9/input
covar_file = merged_data[["family_id", "AGRF code", "Age", "Week 1 Body Mass", "Week 1 Beep test", "Week 1 Distance (m)", "Week 1 Pred VO2max"]]
covar_file.columns=["FID", "IID", "age", "week_1_weight", "week_1_beep", "week_1_distance", "week_1_vo2"]
    #select the baseline for each predictor, so VO2 max base line for 8 week baseline, etc... plus age and set the correct column names
covar_file.to_csv( \
    "./data/plink_inputs/covar_file.tsv", \
    sep="\t", \
    index=False, \
    header=True)
run_bash(" \
    cd ./data/plink_inputs/; \
    awk \
        'BEGIN{FS=\"\t\"}{if(NR<=10){print $0}}' \
        covar_file.tsv")


print_text("perform covariate selection based on association levels", header=3)
print("In the HINT study, they did this: 'the Baseline VO2peak, the individual study and the PC6 from the principal component analysis were found to be significantly associated with VO2peak response and were included as covariates in analysis. Age and sex were not associated with the trait. Our findings did not change when age and sex were also included in the association analysis. Thus, we included covariates based on a posteriori instead of a priori knowledge'")
        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
print_text("create a dict with correct, simplified names for the covariates", header=4)
dict_names_covs = {
    "Age": "age", \
    "sex_code": "sex_code", \
    "Week 1 Body Mass": "week_1_weight", \
    "Week 1 Beep test": "week_1_beep", \
    "Week 1 Distance (m)": "week_1_distance", \
    "Week 1 Pred VO2max": "week_1_vo2"}
print(dict_names_covs)

print_text("plot distribution of each change phenotype to visually check normality", header=4)
import matplotlib.pyplot as plt
from scipy.stats import shapiro
from sklearn.preprocessing import quantile_transform
#change_pheno="distance_change"
for change_pheno in ["weight_change", "beep_change", "distance_change", "vo2_change"]:
    print("\n \n Plotting the distribution of " + change_pheno)

    #remove rows with missing for the phenotype of interest
    modeling_data = merged_data \
        .loc[:,[change_pheno]] \
        .apply(lambda x: np.nan if (x[change_pheno]==new_nan_value) else x, axis=1) \
        .dropna()
        #select the interest column
        #for each row
            #set NA if the column is the new_nan_value, i.e., missing in plink.
            #else do not change anything of the row
        #drop rows with NaN

    print("The number of missing removed is correct?")
    check_na = merged_data.shape[0] - modeling_data.shape[0] == count_missing[change_pheno]
        #the number of rows lost in modeling_data should be the same than the number of missing for the selected pheno
    if check_na:
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NAN IN " + change_pheno)

    print("performing quantile transformation")
    pheno_transformed = quantile_transform( \
        X=modeling_data[change_pheno].values.reshape(-1, 1), \
        n_quantiles=500, 
        output_distribution="normal")
        #This method transforms the features to follow a uniform or a normal distribution. Therefore, for a given feature, this transformation tends to spread out the most frequent values. It ALSO REDUCES THE IMPACT OF (MARGINAL) OUTLIERS: this is therefore a robust preprocessing scheme. Note that this transform is non-linear. IT MAY DISTORT LINEAR CORRELATIONS BETWEEN VARIABLES MEASURED AT THE SAME SCALE BUT RENDERS VARIABLES MEASURED AT DIFFERENT SCALES MORE DIRECTLY COMPARABLE.
            #This should not be a problem because our phenotypes and covariates are not in the same scale, some are close to zero, other closer to 10, 100, 1000... so we have different scales.
        #this would be another cause of data leak
            #n_quantiles
                #Number of quantiles to be computed. It corresponds to the number of landmarks used to discretize the cumulative distribution function. If n_quantiles is larger than the number of samples, n_quantiles is set to the number of samples as a larger number of quantiles does not give a better approximation of the cumulative distribution function estimator.
            #output_distribution
                #Marginal distribution for the transformed data. The choices are ‘uniform’ (default) or ‘normal’.
            #https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.quantile_transform.html

    print("perform shapiro test: statistic close 1 and p-value>0.05 indicates normality")
    print(shapiro(pheno_transformed))
        #The Shapiro-Wilk test tests the null hypothesis that the data was drawn from a normal distribution.
            #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html

    print("plotting the raw vs. transformed phenotype and the histogram of the transformed phenotype")
    fig, (ax1, ax2) = plt.subplots(2, 1)
    #make more vertical space between subplots
    fig.subplots_adjust(hspace=0.4)
            #https://stackoverflow.com/a/5159405/12772630
    #correlate the phenotype vs transformed phenotype
    ax1.scatter(x=modeling_data[change_pheno], y=pheno_transformed, s=0.5)
    ax1.set_ylabel(change_pheno + " - quantile trans")
    ax1.set_xlabel(change_pheno + " - raw")
    #plot histogram of the transformed phenotype
    ax2.hist(pheno_transformed, 50, density=True, alpha=0.4, label="Distribution of " + change_pheno)
    ax2.set_xlabel("Histrogram of " + change_pheno + " after the trasformation")
    plt.savefig( \
        fname="./results/prelim_results/distribution_hist/" + change_pheno + "_trans_hist.png")
    plt.close()
print("IMPORTANT: \
    Our phenotypes are in general more or less normal, but in order to make them more normal, we have performed a quantile transformation. The W statistic of shapiro is almost 1, and the p-value is below 0.05 but not veery small. The HIIT study used the p-value of shapiro test as indication of normality. \
    Given these results, we are going to pre-process phenotypes and covariates directly in plink. It is interesting to apply his preprocesing to the covariates also so we have all variables in the same scale? \
    THINK MORE ABOUT THIS IN THE FUTURE")

print_text("run the associations between phenotypes and covariates", header=4)
#empty dict to save covariates per phenotype
dict_pvalue_covar = { 
    "weight_change": [], \
    "beep_change": [], \
    "distance_change": [], \
    "vo2_change": []}
#define a function to calculate the association between a pair of phenotype and covariate
from scipy.stats import linregress
#change_pheno="distance_change"; covar="Age"
def fun_assoc(change_pheno, covar):

    #remove rows with missing for any of the interest columns
    modeling_data = merged_data \
        .loc[:,[change_pheno, covar]] \
        .apply(lambda x: np.nan if (x[change_pheno]==new_nan_value) | (x[covar]==new_nan_value) else x, axis=1) \
        .dropna()
        #select the interest columns
        #for each row
            #set NA if any of the two columns is new_nan_value, i.e., missing in plink. The np.nan goes to all columns in that row
            #else do not change anything of the row
        #drop rows with NaN

    #model
    model_results = linregress( \
        x=modeling_data[covar], \
        y=modeling_data[change_pheno], \
        alternative='two-sided')
        #the phenotype is the response
        #alternative
            #'two-sided': the slope of the regression line is nonzero (greater or lower than zero, i.e., negative or positive association)
    print(model_results)

    #add the covariate to the list of covariates of the selected pheno if the p_value is below 0.1. We use the correct, simplified name of the covariate
    if model_results.pvalue < 0.1:
        dict_pvalue_covar[change_pheno].append(dict_names_covs[covar])
print("IMPORTANT: \
    We are analyzing here sex_code (0,1,2) as a continuous variable, but we should use a logistic regression. \
    This is not a big deal because we are using these associations just to include or not a covariate, and sex is likely very important in all the phenotypes studied. \
    In the future use statsmodels module to perform logistic regression.")
    #https://www.statsmodels.org/stable/discretemod.html

#for each phenotype
#change_pheno="distance_change"; covar="Age"
for change_pheno in ["weight_change", "beep_change", "distance_change", "vo2_change"]:
    print("\n \n #" + change_pheno + "#")
    #test the association with each covariate
    for covar in ["Age", "sex_code", "Week 1 Body Mass"]:
        print("\n" + covar)
        fun_assoc(change_pheno=change_pheno, covar=covar)
    #if the phenotype is beep_change, we also want to check the association with baseline beep
    if change_pheno == "beep_change":
        covar="Week 1 Beep test"
        print("\n " + covar)
        fun_assoc(change_pheno=change_pheno, covar=covar)
    #if the phenotype is distance_change, we also want to check the association with baseline distance
    elif change_pheno == "distance_change":
        covar = "Week 1 Distance (m)"
        print("\n " + covar)
        fun_assoc(change_pheno=change_pheno, covar=covar)
    #if the phenotype is vo2_change, we also want to check the association with baseline vo2 max
    elif change_pheno == "vo2_change":
        covar = "Week 1 Pred VO2max"
        print("\n " + covar)
        fun_assoc(change_pheno=change_pheno, covar=covar)
print("All covariates except Age for beep_test and distance which is marginally significant (p=0.06). See dict with covariates")
print(dict_pvalue_covar)

print_text("add to the dict of covariates a new entry to indicate to plink if sex is going to be used or not as covariate", header=4)
#make deep copy of the previous dict (no conection between source and new copy)
import copy
dict_pvalue_covar_final = copy.deepcopy(dict_pvalue_covar)
#for each phenotype
#key="distance_change"
for key in dict_pvalue_covar_final:
    #join all covariates of the selected phenotype in a single string and save as a list of one item
    dict_pvalue_covar_final[key] = [",".join(dict_pvalue_covar_final[key])]
    
    #remove sex and add it as a separate string if it is present
    if "sex_code" in dict_pvalue_covar_final[key][0]:
        dict_pvalue_covar_final[key][0] = dict_pvalue_covar_final[key][0].replace("sex_code,", "")
            #SEX is included as an argument in plink!!
        dict_pvalue_covar_final[key].append("sex")
    else:
        dict_pvalue_covar_final[key].append("")
print("See the new dict")
print(dict_pvalue_covar_final)



#USA LA VARIABLE  DE SEXO DE PLINK
#<datastem>.fam has one row per individual and six columns, which provide the Individual ID, the Family ID, as well as Maternal and Paternal IDs, Sex and Phenotype. Note that LDAK only uses the first two IDs; the remaining four columns are ignored (so to use sex as a covariate or to provide phenotypic values, these need to be supplied separately)

#UN ARHIVO.PHENO PARA CADA PHENO!
    #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). Binary phenotypes should only take values 0 (control), 1 (case) or NA (missing). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype).
    #http://dougspeed.com/phenotypes-and-covariates/




#ask David the last questions about the phenotypes
#explain the selected covariates
#decide if we are going to analyze the three training phenotypes
#more questions
    #- ask david that from sample 1161 to 1376, age is integer, not float, in contrast with almost all the rest samples. This is ok?
    #- in some phenotypes, some some samples have value of 0 and others have no value. I guess zero should be NA, right?
    #    - body mass week 1 and 8
    #    - VO2 max week 1
    #- sample 1194, the value for week 8 beep includes a letter: 11.1O. I guess I can safely change that "o" letter by zero.
    #- I guess that the sheet "DNA with only wk1" includes genotyped samples with only data for the first week, not week 8. So I should only use the sheet "All DNA samples" and discard the 42 samples at the bottomn with NA for all columns except the AGRF code.
    #- some rows are coloured, there is something special about these samples it could be relevant for the analysis?