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



######################################
######## ASSOCIATION ANALYSIS ########
######################################

#This script will perform association analyses. 

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
    mkdir -p ./results/prelim_results; \
    ls -l")



#################################
# approach to calculate the PRS #
#################################

#Data and validation scheme
    #Make a Polygenic predictor score (PPS) of VO2 max trainibility
        #In order to calculate a PPS, we need effect sizes and p-values for the association between SNPs and the trait obtained from a different study than the one used to train the scores, so we avoid overestimating its predictive power.
            #"Accuracy measurements can be inflated if the discovery GWAS and the target cohort share individuals"
                #Genome-wide association studies
                    #https://www.nature.com/articles/s43586-021-00056-9
            #if the score has a small effect size but it is associated and robust to confounders, you can at least say that there is a connection between genetics and the trait
                #"However, if the results are shown to be robust to confounding (see Population genetic structure and the generalizability of PRSs), then the effect size is not important if the aim is only to establish whether an association exists, which may provide etiological insight."
                    #Interpretation of PRS-trait associations
                        #https://www.nature.com/articles/s41596-020-0353-1
        #External validation
            #I have checked previous GWAS about VO2 max trainability and they all have much lower sample size than this study and hence, less power to detect true associations than us. Therefore, I think makes sense to use this study to associate genes and trainability rather than use p-values from previous GWAS to calculate the scores.
            #The problem is that we then need at least a second cohort to train the score and see how they predict. 
                #In the PRS tutorial, they say that you need at least two cohort, the base to obtain GWAS summary statistics and then the target to do cross-validation with the PRS. 
                #Of course, a test set would be then necessary for truly real predictive power, but in the absence of that, we could just use CV.
                    #Overfitting in PRS-trait association testing
                        #https://www.nature.com/articles/s41596-020-0353-1
                #maybe Improve-HIIT could be used as a test set?
            #David, I have seen you are co-author of one of the previous GWAS on V02 max trainability, the HIIT-Predict study. Do you think there would have been any possibility to talk with them about using their data for training polygenic scores purposes? They say in the paper data is avaiable upon reasonable request.
                #we should check that samples from both cohorts are not related, as both come from Australia
                #maybe Improve-HIIT could be used as a test set?
                #Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
                    #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7#Sec21
            #maybe we could even use our PPS developed in the cohort to test whether it interacts with the intervention type in HIIT?
                #In other words, try also models with the interaction PPSxTrainingType
                #maybe high PPS values associate with higher increase in VO2 max only in high volumen HIIT, while not in constant training. This would be a small window to training personalization based on genetics.
                #I am just speculating in this point because I am not sure if there is enough sample size in Predict-HIIT to do test interaction and also split data for training-evaluation.
            #In summary, we would take advantage of your experiment and machine learning to develop a new PPS of trainability that would be then tested.
        #An alternative would be 
            #just ask for p-values/betas from previous GWAS on VO2 max change like the HERITAGE of Predict-HIIT studies. 
                #This is not genotypes, just the summary statistics.
                #the problem is that we would be using studies with smaller sample size to discover variants. Indeed, the Predict-HIIT study did not have any SNP with a significant p-value after correction and the paper say that they are probably underpowered. Their power calculations indicates that "a cohort of 2960 samples would have 80% power to detect a quantitative trait with a true heritability of 30%.". Therefore, it would be better to use our study to discover, given its larger sample size.
                #if we have data from different gwas we can obtain meta-p-values or something like that
                    #Deep learning-based polygenic risk analysis for Alzheimer’s disease prediction
                #previos GWAS about trainability
                    #some GWAS of VO2 max trainability
                        #Genes to predict VO2max trainability: a systematic review
                            #https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4192-6
                        #Genome wide association study of response to interval and continuous exercise training: the Predict-HIIT study
                            #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7
                        #Genotype-Phenotype Models Predicting VO2max Response to High-Intensity Interval Training in Physically Inactive Chinese
            #Obtain p-values/betas from GWAS of baseline VO2 max in the UK Biobank. 
                #This biobank includes estimates of VO2 max based on a submaximal cycle ramp test and it has a large sample size, which is good for power in the base cohort.
                #there is already a study that have done a GWAS on this data, so we could just ask them for the summary statistics. In the worst case scenario, we could do the GWAS on our own.
                #then we can test the polygenic scores in our cohort.
                    #The genetic case for cardiorespiratory fitness as a clinical vital sign and the routine prescription of physical activity in healthcare
                        #https://www.medrxiv.org/content/10.1101/2020.12.08.20243337v2
                #problems
                    #we would use a correlate of VO2 change, not the trait itself. this can be done, but it is not optimal
                        #See "Predicting different traits and exploiting multiple PRSs" 
                            #https://www.nature.com/articles/s41596-020-0353-1
                    #we would not take advantage of the greater sample size we have in this study compared to previous interventions in order to derive the polygenic score.
            #Make a Polygenic predictor score of VO2 max
                #It is less interesting, but it would be easier because we could use the previous GWAS studies to obtain effect sizes, then train in our dataset and even validate in a third dataset, as the UK biobank has VO2 max data along with genetics.
        #THE OPTION SELECTED BY DAVID AND JONATAN
            #Just do GWAS of trainability with our cohort
                #I would not got for this option as we do not have a validation cohort, only the discovery cohort.
                #In some cases like studies of rare diseases, they use interval validation of the GWAS, i.e., using the same cohort but applying resampling methods. 
                    #This is an option but it would be a clear limitation of the study.
                    #the gold standard in this point is to use an independent and diverse cohort for validation
                #We also have a good point in the fact we have a very specific cohort. It is very difficult to find a population with similar characteristics. This is just like having a cohort of a rare disease.
                #Data leaking
                    #we remove SNPs based on the allele frequencies of all samples, including those of the validation dataset.
                    #We also remove samples based on relatedness and PCA considering test and evaluation samples, but you have to related individuals and as you have more data you can do it better...
                    #the most problematic point would be to use PCAs as predictors...
                #info
                    #Evaluation of a two-step iterative resampling procedure for internal validation of genome-wide association studies
                        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4859941/
                    #In "Genome-wide association studies", nature review methods, they say that internal validation is a possibility, but no the gold standard
                        #https://www.nature.com/articles/s43586-021-00056-9
                    #David sent a paper also using this approach
                        #Personalized Nutrition by Prediction of Glycemic Responses
    #mail to david
        #Also, I have been thinking and reading about the approach for doing the validation and I have some thoughts on this. I have found several papers talking about internal validation of GWAS hits for cohorts where replication is not possible like in rare disease (REFS). As we discussed, this is also our case, because our population has very specific characteristics.
        #I have found, however, more insistence for external validation when reading about poligenic scores due to the risk of overfitting. For example, in this Nature Protocol (REFS), the authors insists in the necessity of having complete independent cohorts when doing the scores. 
        #As we already discussed, we could just split our dataset in training and evaluation, repeating the associations and scores hundreds of times, but here the risk of "data leak" is extremely high in my opinion. In other words, when we use the model on the evaluation dataset, this data is not completely new, the model has already received information of that dataset when doing the training in the train dataset.
        #For example, 
            #we are removing SNPs based on the minor allele frequency across all samples, thus we are pruning our genetic information using data from training and evaluation samples. 
            #In the same vein, we are going to do imputation, which is the "estimation" of genotypes of SNPs we have not genotyped using the information of SNPs around we did genotype and the information about structural variation in the genomes of the ancestry (~race) of our studied population. So we are extending our data using again all samples.
            #we are removing samples that are related (e.g., cousins) and, of course, we need to use the whole set of samples.
            #Last but not least, we need to include a PCA with ancestry information in the models. That PCA is calculated using all samples, so even the predictors already contain information of both training and evaluation.
        #Therefore, we do not have truly unseen (evaluation) data, so the models are going to be training knowing before the data that will be used for evaluate them later. This can artificially inflate a lot the estimation of the predictive power of the models, so we cannot be confident that the model would be the equally good in other population of the same characteristics.
        #We could try to do the quality control separated between training and evaluation, but then we would have to repeat the whole analyses (not only associations) hundreds of times and, more important, the QC does not work with a low number of samples, like the one we would have in the smaller evaluation dataset.
        #Finally, the problem of overfitting does not end here. We will likely need to use models where some hyperparameters need to be tuned or optimized (e.g., degree of correction of p-values in a lasso regression). Therefore, the evaluation set should be split again in training-evaluation, so we train the model with a different set of parameters and then check the performance of these different set of parameters in the new evaluation set. This is a common approach in machine learning, but it is also specifically recommended for polygenic scores, because the lack of parameter optimization can greatly reduce the power of the scores.
        #In summary, using only our data would lead to split the data several times leading to low sample sizes. The risk of overfitting would be really great, limiting the ability of generalize of our results, i.e., we cannot know how powerful the score is for other individuals (even of the same race). Note that the main motivation/necessity to do training/evaluation splits is to obtain less biased estimations of predictive power, but a data leakage greatly difficult this.
        #We can, of course, do what we discussed and go for the internal validation, saying that we could not find a sample with similar characteristics and making e a clear disclaimer about the potential inflation of the predictive power estimation. But, if possible, it would be much better to avoid this. 
        #there is no rush to talk about this as before anything. As we indeed decided, I have to first do the QC and run the GWAS on our dataset. Then, we can see what to do. 

#Specific section for DNNs becuase it can be interesting, but in general use a battery of ML models
    #alzheimer paper that uses a very cool approach to predict using polygenic scores
        #https://www.nature.com/articles/s43856-023-00269-x
    #other option
        #SNPRS:Stacked Neural network for predicting Polygenic Risk Score

#things to do about associations
    #add body mass baseline (week 1) for the associations of vO2 max
    #stratified analysis by sex
        #According to Jonatan, this is a big deal right now, being a movement to report the differences between sexes if possible.
        #Remember we have a great imbalance, you have 1057 men! so around 300 females!
        #We can do this for the association analyses to see what happens
        #I do not think it is worth it for poligenic scores, because there is big big reduction of sample size with the consequent reduction of power to detect associations and develop the score
            #think?
    #impact of body weight
        #test removing change weight as a covariate for vO2 max
        #check body mass impact on predicitive power VO2 max
        #not sure about this, reduction of power across hundreds of SNPs? or in the final polygenic score? if we do it for this covariate, why not for the other three?
            #I indeed plant to use lasso as one of the tested models, and this model penalize predictors, so maybe it is not necessary?
                #lasso penalyze all predictors in the same degree?
    #think about using the batch as covariate
        #https://academic.oup.com/cardiovascres/article/118/Supplement_1/cvac066.013/6605381
    #we can leave X, Y, PAR and MT for now
        #but CHECK WHAT TO DO IN POLYGENIC SCORES
        #they can impede the score?


    #check manhatan plots
        #there is a strange gap in VO2 max for one of the first chromosomes




######################
# prepare covariates #
######################
print_text("prepare covariates", header=1)



print_text("load pheno data", header=2)
print_text("This include reported sex and VO2 max data. I have checked that the data is the same directly reading from excel than converting to csv", header=3)
run_bash(" \
    cp \
        ../quality_control/data/pheno_data/'combact gene DNA GWAS 23062022.xlsx' \
        ./data/pheno_data/pheno_data.xlsx; \
    echo 'pheno_data.xlsx comes from ../quality_control/data/pheno_data/combact gene DNA GWAS 23062022.xlsx, which is the raw excel and got from D. Bishop' > ./data/pheno_data/README.txt; \
    ls -l ./data/pheno_data/")
import pandas as pd
import numpy as np
pheno_data = pd.read_excel(
    "./data/pheno_data/pheno_data.xlsx",
    header=0,
    sheet_name="All DNA samples")


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

print_text("look the pheno data", header=4)
print(pheno_data)



print_text("load fam file", header=2)
print_text("load the fam file of the current steps I am working on for QC (date 08/24/2023", header=3)
fam_file = pd.read_csv( \
    "../quality_control/data/genetic_data/quality_control/08_loop_maf_missing/loop_maf_missing_2.fam", \
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

print_text("see the fam file", header=3)
print(fam_file)



print_text("merge and create file with the final set of samples and phenotypes", header=2)
print_text("merge pheno data and fam file", header=3)
merged_data_raw = pheno_data.merge( \
    fam_file, \
    on="AGRF code", \
    how="right")
    #how="right":
        #use only keys from right frame, similar to a SQL right outer join
        #we use right, i.e., fam file (genetic data), for the keys because we only want to retain those samples that remained after the filters done in plink. We want to take advantage of the filtering work already done.
        #the only sample we need to recover from pheno_data is 2399LDJA, because this sample is 2397LDJA in plink and it is just a mislabel problem so we can recover its phenotypic data.
#make a deep copy to leave the raw merged data intact
merged_data = merged_data_raw.copy(deep=True)
    #deep=True
        #Modifications to the data or indices of the copy will not be reflected in the original object (see notes below).
print(merged_data)
print(merged_data.describe().T)

print_text("total sum of NaNs or zero cases per column", header=4)
nan_zeros_1 = merged_data.apply(lambda x: sum((x==0) | (x.isna())), axis=0)
print(nan_zeros_1)


print_text("check and remove cases with zero in phenotypes", header=3)
print_text("sum the number of Zeros per column", header=4)
print(merged_data.apply(lambda x: sum(x==0), axis=0))
    #note that the number of zeros does NOT have to match that of the pheno_data because we have already filtered many samples using maf and call rate in plink
    #for example, ["8502JGMJ", "9601AVJM"] have zero for "Week 1 Pred VO2max", but they are not included in the fam file
        #merged_data.loc[merged_data["AGRF code"].isin(["8502JGMJ", "9601AVJM"]),:]
        #fam_file.loc[fam_file["AGRF code"].isin(["8502JGMJ", "9601AVJM"]),:]
print_text("Weight: we have zeros for weight and VO2 max. Given that weight is going to be used for the rest of variables as covariate, we will set as missing the weight of these samples. Importantly, a weight of zero does not make sense.", header=4)
merged_data.loc[merged_data["Week 1 Body Mass"]==0, "Week 1 Body Mass"] = np.nan
merged_data.loc[merged_data["Week 8 Body Mass"]==0, "Week 8 Body Mass"] = np.nan

print_text("VO2 max: In the case of VO2 max, this is probably caused by the equation. It is indeed a very strange value because the beep test data is ok and the VO2 max of the 8th week is also ok. We are going to set this as nan", header=4)
print(merged_data.loc[merged_data["Week 1 Pred VO2max"]==0, ["Week 1 Beep test", "Week 8 beep test", "Week 1 Pred VO2max", "Week 8 Pred VO2max",]])
merged_data.loc[merged_data["Week 1 Pred VO2max"]==0, "Week 1 Pred VO2max"] = np.nan

print_text("Sex_code: We have missing for 42 samples which are the ones with genetic data that are completely empty in the pheno_data PLUS 2397LDJA, the misslabeled sample with different Id geno and pheno data just by 1 number (see below). These are samples without phenotypic data, thus we cannot know the self-reported sex. Zero is ok here, so we do not do anything", header=4)
print("Do we have exactly 42 samples with sex_code=0?")
print(merged_data.loc[merged_data["sex_code"]==0, :].shape[0]==41+1)
print("has 2397LDJA zero as sex code?")
print(merged_data.loc[merged_data["AGRF code"] == "2397LDJA", "sex_code"]==0)

print_text("Parents ID: We do not have any parent in the cohort, so they all should be zero", header=4)
print(merged_data.loc[(merged_data["ID_father"]==0) & (merged_data["ID_mother"]==0), :].shape[0] == merged_data.shape[0])

print_text("sum again the number of Zeros per column", header=4)
print(merged_data.apply(lambda x: sum(x==0), axis=0))
print("We only have zeros in sex_code and parents IDs, which is ok, see above")


print_text("deal with 2397LDJA/2399LDJA", header=3)
print("I guess this is the misslabeled sample, it has one ID in genetic and other different ID in pheno. As the postdoc of David said: I think it is a mislabelling of the last digit of the number (the labelling was very hard to read on some of the blood samples). So, I think 2397LDJA; ILGSA24-17303 is 2399LDJA in the excel file")
mislabelled_sample = merged_data.loc[merged_data["AGRF code"].isin(["2397LDJA", "2399LDJA"]), :]
print("If we have 2397LDJA and is all missing")
if mislabelled_sample.shape[0] == 1:
    mislabel_all_na =  \
        sum( \
            mislabelled_sample.loc[:, ~mislabelled_sample.columns.isin(["AGRF code", "family_id", "ID_father", "ID_mother", "sex_code", "phenotype_value"])].isna() \
            .values[0]) == \
        sum( \
            ~mislabelled_sample \
                .columns \
                .isin(["AGRF code", "family_id", "ID_father", "ID_mother", "sex_code", "phenotype_value"]))
        #make the sum of phenotypic entries (not fam fields) of the mislabeled row that are missing
        #make the sum of the columns that are not phenotypic entries (not fam fields)
        #both sums should be same, meaning that all phenotypic entries (not fam fields) are missing
    if (mislabelled_sample["AGRF code"]=="2397LDJA").values[0] & (mislabel_all_na):
        print("change phenotypic values of 2397LDJA from missing to the values of 2399LDJA")
        merged_data \
            .iloc[ \
                np.where(merged_data["AGRF code"] == "2397LDJA")[0], \
                np.where(merged_data.columns.isin(pheno_data.columns))[0]] = \
        pheno_data \
            .loc[pheno_data["AGRF code"] == "2399LDJA", :] \
            .squeeze(axis=0)
            #in merged data
                #select the row of 2397LDJA
                #and the phenoypic columns. 
                    #we need to use iloc to select the row to be changed, we get an error with loc
                    #https://stackoverflow.com/questions/45241992/updating-a-row-in-a-dataframe-with-values-from-a-numpy-array-or-list
            #as new values use the row of 2399LDJA in pheno_data and the squeeze to pandas series
            #IMPORTANT: The pheno columns are in the same order in merged data and in pheno_data because the source DF is pheno_data, as we did pheno_data.merge(fam) We need to maintain this and not do 
        #extract the modified row
        modified_row = merged_data \
            .loc[ \
                merged_data["AGRF code"] == "2399LDJA", \
                merged_data.columns.isin(pheno_data.columns)] \
            .squeeze()
        #sex_code is not in pheno_data as it comes from plink, we need to update that field using the Gender info from the modified row
        if modified_row["Gender"]=="M":  
            merged_data.loc[merged_data["AGRF code"] == "2399LDJA", "sex_code"] = 1
        elif modified_row["Gender"]=="F":
            merged_data.loc[merged_data["AGRF code"] == "2399LDJA", "sex_code"] = 2
        else:
            merged_data.loc[merged_data["AGRF code"] == "2399LDJA", "sex_code"] = 0      
            #The sex code in plink is 1 for male, 2 for female and 0 for unknown
                #https://www.cog-genomics.org/plink/1.9/formats#fam

        print("Check the new rows is identical to 2399LDJA in pheno_data")
        print(modified_row.equals( \
            pheno_data \
                .loc[pheno_data["AGRF code"] == "2399LDJA", :] \
                .squeeze(axis=0)))
        print("Check NO NaN is in the new row")
        print(True not in modified_row.isna().values)
        print("Print the new row")
        print(modified_row)
        print("check the sex_code of the modified row is correct")
        sex_code_check = merged_data.loc[merged_data["AGRF code"] == "2399LDJA", "sex_code"]
        if (sex_code_check == 1).values:
            print(pheno_data.loc[pheno_data["AGRF code"] == "2399LDJA", "Gender"] == "M")
        elif (sex_code_check == 2).values:
            print(pheno_data.loc[pheno_data["AGRF code"] == "2399LDJA", "Gender"] == "F")
    else:
        raise ValueError("ERROR: FALSE! WE DO HAVE A PROBLEM WITH 2397LDJA/2399LDJA")
else:
    raise ValueError("ERROR: FALSE! WE DO HAVE A PROBLEM WITH 2397LDJA/2399LDJA")


print_text("create new variables for the change before and after", header=3)
print_text("make the operations", header=4)
merged_data["weight_change"] = merged_data["Week 8 Body Mass"]-merged_data["Week 1 Body Mass"]
merged_data["beep_change"] = merged_data["Week 8 beep test"]-merged_data["Week 1 Beep test"]
merged_data["vo2_change"] = merged_data["Week 8 Pred VO2max"]-merged_data["Week 1 Pred VO2max"]
    #NA remains as NA

print_text("check we have the correct NANs", header=4)
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
print("VO2 max change")
print( \
    sum(merged_data["vo2_change"].isna()) == \
    sum( \
        (merged_data["Week 1 Pred VO2max"].isna()) | \
        (merged_data["Week 8 Pred VO2max"].isna())))
    #the total number of NAN for the change variables should be the same than the sum of cases that are NAN for week 1 or week 2 of the corresponding variable.


print_text("Change the default NaN value to a negative value that is smaller than the maximum value in absolute value of the new change variables", header=3)
print("Note that plink considers -9 as missing, but we can have -9 as a real value because we are calculating differences between week 1 and 8. We have to deal with that and select a suitable missing value, probably a large negative number. Larger than the maximum value for the difference phenotypes")
        #https://www.cog-genomics.org/plink/1.9/input#pheno_encoding

print_text("calculate absolute value and then get the max value skipping NaNs for each variable", header=4)
max_weight_change = merged_data["weight_change"].abs().max(skipna=True)
max_beep_change = merged_data["beep_change"].abs().max(skipna=True)
max_vo2_change = merged_data["vo2_change"].abs().max(skipna=True)
print("Maximum different in weight of " + str(max_weight_change))
print("Maximum different in beep test of " + str(max_beep_change))
print("Maximum different in VO2 of " + str(max_vo2_change))

print_text("select the max value from all three variables, sum 20 and then convert to negative. This is the new missing value", header=4)
new_nan_value = int(-(np.round(np.max([max_weight_change, max_beep_change, max_vo2_change]))+20))
print(new_nan_value)
    #we save it as integer, although pandas will change it to float (adding .0), we will change back that using "new_nan_value" as integer.

print_text("fill NaN cases with the new missing value using 'fillna' of pandas", header=4)
merged_data = merged_data.fillna(new_nan_value)
print(merged_data)

print_text("fill -9 cases of phenotype_value with the new missing value", header=4)
merged_data.loc[merged_data["phenotype_value"]==-9, "phenotype_value"] = new_nan_value
print(merged_data)

print_text("check we do NOT have the row coming from the pheno excel with missing for all columns", header=4)
all_missing_row = np.where(merged_data.apply(lambda x: x == new_nan_value, axis=0).apply(np.sum, axis=1) == len(merged_data.columns))[0]
    #look for any value equal to the new missing across rows
    #then sum all True cases per row
    #check if any row has a many Trues as columns we have, i.e., a row with missing for all columns
    #get the index of the row and extract it
if len(all_missing_row)==0:
    print("We do not have the row with ALL missing from the pheno excel")
else:
    raise ValueError("ERROR: FALSE! We DO have the row with ALL missing from the pheno excel, but we should not have as we used only the keys of the fam file")

print_text("check we have the correct index in order to filter by index and iloc, i.e., no gaps in the indexes", header=4)
print(np.array_equal(merged_data.index, range(0,merged_data.shape[0])))

print_text("check the problematic samples are NOT included", header=4)
print(sum(merged_data["AGRF code"].isin(["1100JHJM", "1200JPJM", "7800AGSO", "2397LDJA"]))==0)

print_text("check we do not have duplicated samples and no sample with NaN", header=4)
print(sum(merged_data["AGRF code"].duplicated(keep=False)) == 0)
    #keep=False
        #Mark all duplicates as True
print(sum(merged_data["AGRF code"].isna()) == 0)


print_text("Removing the float point of the new missing", header=3)
print("According to plink man: 'Missing phenotypes are normally expected to be encoded as -9. You can change this to another integer with --missing-phenotype. IMPORTANTLY, floating point values are now disallowed due to rounding issues, and nonnumeric values such as 'NA' are rejected since they're treated as missing phenotypes no matter what.)'.")
    #https://www.cog-genomics.org/plink/1.9/input#pheno_encoding
print("Therefore, we need to change our missing value from float to integer. We originally created the new missing value as integer but pandas converts it to float when included in a float column by just adding '.0'.")
print("For each phenotype, we are to convert the column to 'object', which is a type that can be strings and floats. Then convert the cases of missing to integer effectively removing the float point.")
    #https://stackoverflow.com/a/34524168/12772630
#selected_phenotype="Week 1 Body Mass"
for selected_phenotype in merged_data.columns:

    #if the dtype of the selected phenotype starts with "float"
    if str(merged_data[selected_phenotype].dtype).startswith("float"):
        print(selected_phenotype)

        #convert that column to object
        merged_data[selected_phenotype] = merged_data[selected_phenotype].astype(object)

        #change to integer the missing cases of the selected phenotype
        merged_data.loc[merged_data[selected_phenotype] == new_nan_value, selected_phenotype] = int(new_nan_value)

        print("the new missing is now an integer?")
        if isinstance(merged_data.loc[merged_data[selected_phenotype] == new_nan_value, selected_phenotype].unique()[0], int):
                #https://stackoverflow.com/a/4541167/12772630
            print("TRUE \n")
        else:
            raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM CONVERTING TO INTEGER THE MISSING CASES \n")

print_text("reorder columns", header=4)
col_names = merged_data.columns
    #get the column names
col_index_ordered = np.concatenate(( \
    np.where(col_names == "AGRF code")[0], \
    np.where(col_names != "AGRF code")[0]))
    #make a numpy array with the index of AGRF code and another array with the index of the rest of columns
merged_data = merged_data.iloc[:, col_index_ordered]
    #reorder the columns
    #https://stackoverflow.com/questions/13148429/how-to-change-the-order-of-dataframe-columns

print_text("final exploration of the data", header=3)
print_text("count missing per column", header=4)
count_missing = merged_data.apply(lambda x: sum(x==new_nan_value), axis=0)
print(count_missing)

print_text("check we have the correct number of missing based on the raw merged data", header=4)
print("In general, we have to subtract 1 from the number of missing in the raw dataset because the misslabeled sample (2397LDJA/2399LDJA) was missing in that file but not in the new one, where we did the required corrections")
count_missing["Gender"] == sum(merged_data_raw["Gender"].isna())-1
count_missing["Age"] == sum(merged_data_raw["Age"].isna())-1
#selected_phenotype="Week 1 Body Mass"
print("In the case of phenotypes, we converted to missing no only NaN, but also zero values and -9, so we have to consider samples with these values in the raw dataset")
for selected_phenotype in ["Week 1 Body Mass", "Week 8 Body Mass",
       "Week 1 Beep test", "Week 8 beep test", "Week 1 Pred VO2max",
       "Week 8 Pred VO2max", "phenotype_value"]:
    print(selected_phenotype)
    if selected_phenotype != "phenotype_value":
        print( \
            count_missing[selected_phenotype] == \
            sum( \
                (merged_data_raw[selected_phenotype]==0) | \
                (merged_data_raw[selected_phenotype].isna()) | \
                (merged_data_raw[selected_phenotype]==-9))-1)
    elif selected_phenotype == "phenotype_value":
        print( \
            count_missing[selected_phenotype] == \
            sum( \
                (merged_data_raw[selected_phenotype]==0) | \
                (merged_data_raw[selected_phenotype].isna()) | \
                (merged_data_raw[selected_phenotype]==-9)))
            #the mislabeled sample always has missing, before and after the corrections, so we do not need to subtract 1.
print("In the case of family/parents Ids and sex code, we should not have any missing, but 0 for unknown, so we do not need to subtract 1 for the mislabelled sample (2397LDJA/2399LDJA)")   
count_missing["family_id"] == sum(merged_data_raw["family_id"].isna())
count_missing["ID_father"] == sum(merged_data_raw["ID_father"].isna())
count_missing["ID_mother"] == sum(merged_data_raw["ID_mother"].isna())
count_missing["sex_code"] == sum(merged_data_raw["sex_code"].isna())
print("check also that the family/parents IDs and the sex code are exactly the same if we do not consider the mislabelled sample (2397LDJA/2399LDJA), before and after the corrections we have non-missing for that sample in these variables. Unknown is zero in these variables")
print(merged_data_raw \
        .loc[ \
            merged_data_raw["AGRF code"]!="2397LDJA", \
            ["family_id", "ID_father", "ID_mother", "sex_code"]] \
        .equals( \
            merged_data \
                .loc[ \
                    merged_data["AGRF code"]!="2399LDJA", \
                    ["family_id", "ID_father", "ID_mother", "sex_code"]]))

print_text("check that the number to missing is equal to the cases of NaN OR 0 in the initial merged dataset PLUS 1. We have to add 1 because at the begining, the misslabeled sample (2397LDJA/2399LDJA) was not correctly included in our data", header=4)
print( \
    nan_zeros_1 \
        [["Gender", "Age", "Week 1 Body Mass", "Week 8 Body Mass", "Week 1 Beep test", "Week 8 beep test", "Week 1 Pred VO2max", "Week 8 Pred VO2max"]] == \
    (count_missing \
        [["Gender", "Age", "Week 1 Body Mass", "Week 8 Body Mass", "Week 1 Beep test", "Week 8 beep test", "Week 1 Pred VO2max", "Week 8 Pred VO2max"]]+1))
    #nan_zeros_1 was obtained at the begining after merging using two conditions, isna() and ==0

print_text("see the final data", header=4)
print(merged_data)

print_text("see summary statistics of non_pheno columns", header=4)
print(merged_data.describe().T)

print_text("see summary statistics of pheno columns", header=4)
#selected_phenotype="Week 1 Body Mass"
for selected_phenotype in ["Age", "Week 1 Body Mass", "Week 8 Body Mass",
       "Week 1 Beep test", "Week 8 beep test", "Week 1 Pred VO2max",
       "Week 8 Pred VO2max", "weight_change",
       "beep_change", "vo2_change"]:
    print("\n" + selected_phenotype)
    non_missing_subset = merged_data.loc[merged_data[selected_phenotype] != new_nan_value, selected_phenotype]
        #select values of the selected phenotype that are NOT missing
    print("mean: " + str(np.mean(non_missing_subset)))
    print("percentile 2.5: " + str(np.percentile(non_missing_subset, 2.5)))
    print("percentile 25: " + str(np.percentile(non_missing_subset, 25)))
    print("percentile 50: " + str(np.percentile(non_missing_subset, 50)))
    print("percentile 75: " + str(np.percentile(non_missing_subset, 75)))
    print("percentile 97.5: " + str(np.percentile(non_missing_subset, 97.5)))
print("The values of the phenotypes make sense: ")
print("A tendency for lower values in the highest percentiles of weight at week 8, i.e., more decrease for people starting with more weight. There is a tendency of increased weight for lower percentile. The median change of weight is 0, meaning that half of sample did not lose weight. The distribution is centered around zero.")
print("A tendency of increase of beep test and VO2 max in all percentiles at week 8. The higher percentiles of changes tend to be more separated from zero than the lower percentile, i.e., the distribution is displaced to positive values, meaning increase of cardiorespiratory fitness.")


print_text("save the data", header=3)
merged_data.to_csv("./data/pheno_data/pheno_data_cleaned_missing_as_minus_as_"+str(np.abs(new_nan_value))+".tsv",
    sep="\t",
    header=True,
    index=False)
print("see some columns of the file")
run_bash(" \
    cd ./data/pheno_data/; \
    awk \
        'BEGIN{FS=\"\t\"}{if(NR<=10){print $1, $2, $3, $4}}' \
        pheno_data_cleaned_missing_as_minus_51.tsv; \
    ls -l")



print_text("create objects with phenotypes and covariates", header=2)
print_text("create multi pheno files", header=3)
print("According to plink man: \
    '--pheno causes phenotype values to be read from the 3rd column of the specified space- or tab-delimited file, instead of the .fam or .ped file. The FIRST AND SECOND COLUMNS OF THAT FILE MUST CONTAIN FAMILY AND WITHIN-FAMILY IDS, respectively. \
    In combination with --pheno, --mpheno lets you use the (n+2)th column instead of the 3rd column, while --PHENO-NAME LETS YOU SELECT A COLUMN BY TITLE. (In order to use --pheno-name, there must be a header row with first two entries 'FID' and 'IID'.) The new --pheno-merge flag tells PLINK to use the phenotype value in the .fam/.ped file when no value is present in the --pheno file; without it, the phenotype is always treated as missing in this case.'")
    #https://www.cog-genomics.org/plink/1.9/input
print("Therefore, we need family and within family ID and then the phenotypes in a tab-delimited file with a header, including FID and IID as the two first columns.")
multi_pheno_file = merged_data[["family_id", "AGRF code", "weight_change", "beep_change", "vo2_change"]]
multi_pheno_file.columns=["FID", "IID", "weight_change", "beep_change", "vo2_change"]
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
covar_file = merged_data[["family_id", "AGRF code", "Age", "Week 1 Body Mass", "Week 1 Beep test", "Week 1 Pred VO2max"]]
covar_file.columns=["FID", "IID", "age", "week_1_weight", "week_1_beep", "week_1_vo2"]
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

print_text("create a dict the covariate for each phenotype based on their association", header=3)
print_text("create a dict with correct names for the covariates", header=4)
dict_names_covs = {
    "Age": "age", \
    "sex_code": "sex_code", \
    "Week 1 Body Mass": "week_1_weight", \
    "Week 1 Beep test": "week_1_beep", \
    "Week 1 Pred VO2max": "week_1_vo2"}
print(dict_names_covs)

print_text("perform covariate selection based on association levels", header=4)
print("In the HINT study, they did this: 'the Baseline VO2peak, the individual study and the PC6 from the principal component analysis were found to be significantly associated with VO2peak response and were included as covariates in analysis. Age and sex were not associated with the trait. Our findings did not change when age and sex were also included in the association analysis. Thus, we included covariates based on a posteriori instead of a priori knowledge'")
        #https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-021-00733-7

print_text("run the associations", header=4)
#empty dict to save covariates per phenotype
dict_pvalue_covar = { 
    "weight_change": [], \
    "beep_change": [], \
    "vo2_change": []}
#for each phenotype and covariate
from scipy.stats import linregress

#change_pheno="weight_change"; covar="Age"
def fun_assoc(change_pheno, covar):

    #remove rows with missing for any of the interest columns
    modeling_data = merged_data \
        .loc[:,[change_pheno, covar]] \
        .apply(lambda x: np.nan if (x[change_pheno]==-51) | (x[covar]==-51) else x, axis=1) \
        .dropna()
        #select the interest columns
        #for each row
            #set NA if any of the two columns is -51, i.e., missing in plink. The np.nan goes to all columns in that row
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

#apply the function
#quciki check normality of the three phenotypes?

#change_pheno="weight_change"; covar="Age"
for change_pheno in ["weight_change", "beep_change", "vo2_change"]:
    print("\n \n #" + change_pheno + "#")
    for covar in ["Age", "sex_code", "Week 1 Body Mass"]:
        print("\n" + covar)
        #function here!!
    if change_pheno == "beep_change":
        covar="Week 1 Beep test"
        print("\n " + covar)
        modeling_data = merged_data[[change_pheno, covar]].apply(lambda x: np.nan if (x[0]==-51) | (x[1]==-51) else x, axis=1).dropna()
        model_results = linregress(x=modeling_data[covar], y=modeling_data[change_pheno], alternative='two-sided')
        print(model_results)
        if model_results.pvalue < 0.1:
            dict_pvalue_covar[change_pheno].append(dict_names_covs[covar])
    elif change_pheno == "vo2_change":
        covar = "Week 1 Pred VO2max"
        print("\n " + covar)
        modeling_data = merged_data[[change_pheno, covar]].apply(lambda x: np.nan if (x[0]==-51) | (x[1]==-51) else x, axis=1).dropna()
        model_results=linregress(x=modeling_data[covar], y=modeling_data[change_pheno], alternative='two-sided')
        print(model_results)
        if model_results.pvalue < 0.1:
            dict_pvalue_covar[change_pheno].append(dict_names_covs[covar])
print("All covariates except Age for beep_test which is marginally significant (p=0.06)")



#add batch as covariate????
    #note sure if this can cause problems because we are already indicating the family ID (batch) in both covariates and phenos... In the final analyses, if you do not see batch effects, I think you can skip completely this.
    #maybe in the future...



##run the associations
#make a dict with the covariates for each phenotype
import copy
dict_pvalue_covar_final = copy.deepcopy(dict_pvalue_covar)

#key="weight_change"
for key in dict_pvalue_covar_final:
    print(key)
    dict_pvalue_covar_final[key] = [",".join(dict_pvalue_covar_final[key])]
    if "sex_code" in dict_pvalue_covar_final[key][0]:
        dict_pvalue_covar_final[key][0] = dict_pvalue_covar_final[key][0].replace("sex_code,", "")
            #SEX is included as an argument in plink!!
        dict_pvalue_covar_final[key].append("sex")
    else:
        dict_pvalue_covar_final[key].append("")


dict_pvalue_covar_final

#make dict for title plots
dict_titles = {
    "weight_change": "Change of body mass", \
    "beep_change": "Change in beep test", \
    "vo2_change": "Change in predicted VO2 max"}

#create folder to have results of each pheno separated

#run associations and plot each pheno
#pheno="beep_change"
for pheno in ["weight_change", "beep_change", "vo2_change"]:

    #print the phenotype name
    print(pheno)

    #get the covariate names
    covs = dict_pvalue_covar_final[pheno][0]
    sex_cov = dict_pvalue_covar_final[pheno][1]

    #run plink assoc
    run_bash(" \
        mkdir ./results/prelim_results/" + pheno + "/; \
        plink \
            --bfile .././quality_control/data/genetic_data/quality_control/08_loop_maf_missing/loop_maf_missing_2 \
            --linear " + sex_cov + " \
            --missing-phenotype " + str(new_nan_value) + " \
            --pheno ./pheno_file.tsv \
            --pheno-name " + pheno + "\
            --covar ./covar_file.tsv \
            --covar-name " + covs + " \
            --out ./results/prelim_results/" + pheno + "/" + pheno + "; \
        ls -l")
            #Given a quantitative phenotype and possibly some covariates (in a --covar file), --linear writes a linear regression report to plink.assoc.linear.
                #the result is a file with several rows per SNP, having the slope and P for each covariate and the ADD effect of the SNP
                #https://www.cog-genomics.org/plink/1.9/assoc
            #By default, when at least one male and one female is present, sex (male = 1, female = 0) is automatically added as a covariate on X chromosome SNPs, and nowhere else. The 'sex' modifier causes it to be added everywhere, while 'no-x-sex' excludes it.
                #https://www.cog-genomics.org/plink/1.9/assoc

    #use awk to select only rows for ADD effect of the SNP and the first row, then save as tsv.
    run_bash(" \
        cd ./results/prelim_results/" + pheno + "/; \
        awk \
            'BEGIN{FS=\" \"; OFS=\"\t\"}{if($5==\"ADD\" || NR==1){print $1, $2, $3, $4, $5, $6, $7, $8, $9}}'\
            " + pheno + ".assoc.linear > \
        " + pheno + ".assoc.linear.tsv")

    #make plot
    assoc_results = pd.read_csv("./results/prelim_results/" + pheno + "/" + pheno + ".assoc.linear.tsv", sep="\t")
        #use pheno to avoid problems with the delimiter
            #https://stackoverflow.com/questions/69863821/read-output-model-of-plink
    
    print("check we do not have duplicates by position")
    print(sum(assoc_results.duplicated(subset=["CHR", "BP"], keep=False))==0)

    #make manhattan plot with plotly
    import dash_bio
    import plotly.graph_objects as go
    fig = dash_bio.ManhattanPlot(
            dataframe=assoc_results, \
            chrm='CHR', \
            bp='BP', \
            p='P', \
            snp='SNP', \
            gene=None, \
            logp=True, \
            suggestiveline_value=-np.log10(0.05), \
            genomewideline_value=-np.log10(0.05/assoc_results.shape[0]), \
            title="Manhattan Plot - " + dict_titles[pheno])
                #the significance line is bonferroni, which is very stringent
                #https://plotly.com/python/manhattan-plot/
    fig.write_html("./results/prelim_results/" + pheno + "/" + pheno + ".html")
        #you can also save as an interactive html with "fig.write_html("path/to/file.html")"
            #https://plotly.com/python/interactive-html-export/


#271 people had missing value(s), and 1 person was not seen in the covariate file.
#???

##CHECK THE GAPS 
##CHECK THE OVERLAP BETWEEN CHROMSOME 25 AND 26 IN BEEP CHANGE
##very low p-values in chromosome 25 for beep change, maybe remove non-autosomal chromosomes?
##zeros in VO2 max? ask in the mail? Jonatan was indeed talking about not using VO2 max in the previous email.
