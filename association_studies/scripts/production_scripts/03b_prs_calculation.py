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
import argparse
from sklearn.model_selection import train_test_split
import pickle
from sklearn.preprocessing import quantile_transform


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



#######################################
# Passing arguments of python program #
#######################################

#define input arguments to be passed when running this script in bash
parser=argparse.ArgumentParser()
parser.add_argument("--iter_number", type=int, default=1, help="Number of the iteration run. Integer always, None does not work!")
parser.add_argument("--response_variable", type=str, default="distance_change", help="Phenotype to model. String always, None does not work!")
parser.add_argument("--covariate_dataset", type=str, default="small_set_predictors", help="Type of dataset used, short or long list of covariates. String always, None does not work!")
parser.add_argument("--manhattan_plot", type=bool, default=False, help="Ask for a Manhattan plot considering the whole dataset. Bool always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
iter_number = args.iter_number
response_variable = args.response_variable
covariate_dataset = args.covariate_dataset
manhattan_plot = args.manhattan_plot




######################
# folder preparation #
######################
run_bash(" \
    mkdir -p \
        ./data/train_test_sets/" + covariate_dataset + " \
        ./results/final_results/" + covariate_dataset + "; \
    ls -l \
")

# endregion






########################################
# region approach to calculate the PRS #
########################################

#The QC procedures are intentionally conservative, and particular care should be taken in performing them, because small errors can become inflated when aggregated across SNPs in PRS calculation. Researchers should follow established guidelines.
#Prior to performing heritability analysis, it is necessary to perform very careful quality control of samples and predictors. This is because heritability analyses tend to be very sensitive to genotyping errors. Estimates of heritability represent the sum of the contributions of many predictors, so even if the error at individual predictors is small, this can accumulate across predictors to produce a large overall error. Further, estimates of heritability are sensitive to inflation due to population structure and familial relatedness; the former can result in correlations between predictors and the phenotype that are due to ancestry, not causality, while the latter leads to long-range linkage disequilibrium (meaning that we can no longer be confident that the heritability estimates reflect only the heritability contributed by the predictors in the dataset).
#Quality control is also important when constructing prediction models. Genotyping errors, population structure and familial relatedness can all bias estimates of predictor effect sizes, meaning that the model will under perform when used to predict phenotypes for an independent set of samples.
    #O´Really Nature Protocol
        #Tutorial: a guide to performing polygenic risk score analyses. Nature protocol.
        #(e.g., refs. 28–30) — we recommend ref. 29— to perform standard GWAS QC on the base and target data. Since the option of performing QC on the base GWAS data will typically be unavailable, researchers should ensure that high-quality QC was performed on the GWAS data that they utilize. We recommend the following QC criteria for standard analyses: genotyping rate >0.99, sample missingness <0.02, Hardy-Weinberg Equilibrium P >1 × 10−6, heterozygosity within 3 standard deviations of the mean, minor allele frequency (MAF) >1% (MAF >5% if target sample N <1000) and imputation ‘info score’ >0.8.
        #Strand issues, duplicates, sex incosistences and removal of sex chromosomes, population structure and sample relatedness
            #if the aim of an analysis is to model autosomal genetics only, then we recom- mend that all X and Y chromosome SNPs are removed from the base and target data to eliminate the possibility of non- autosomal sex effects influencing results
            #A high degree of relatedness between individuals between the base and target data can also generate inflation of the associa- tion between the PRS and target phenotype. While population structure produces a correlation between genetics and envir- onmental risk factors that requires a broad solution, the pro- blem is exacerbated with inclusion of very close relatives, since they may share the same household environment as well (dis- cussed below). Thus, if genetic data from the relevant base data samples can be accessed, then any closely related individuals (e.g., first/second degree relatives) across base and target sam- ples should be removed to eliminate this risk. If this is not an option, then every effort should be made to select base and target data that are unlikely to contain highly related indivi- duals. However, statistical power can be compromised in ana- lyzing base and target samples from different populations, as discussed below, and so ideally base and target samples should be as similar as possible without risking inclusion of over- lapping or highly related samples.
    #Doug Speed
        #Most heritability analyses are designed for common predictors (typically MAF>0.01 or MAF>0.005), as it is not currently clear how best to include rare variants. We also generally restrict to common predictors when constructing prediction models, as it is very difficult to accurately estimate effect sizes for rare predictors. Further, both heritability and prediction analyses will usually focus on autosomal predictors, as it is not straightforward to include sex chromosomes.
        #Our preferred measure of SNP quality is the information score from imputation, which reflects how closely a SNP's genotypes match those expected given the neighbouring SNPs. We typically exclude SNPs with information score below 0.95 or 0.99. While it is reasonable to also filter SNPs based on other metrics, such as missingness or Hardy-Weinberg Equilibrium, in general we find this is not necessary (at least for common SNPs). Note that some imputation software, such as IMPUTE2, compute multiple information scores, in which case we would filter based on all scores.
    #We have applied all these filters
        #snp call rate >0.99 (--geno 0.01)
        #sample call rate >0.99 (--mind 0.01), this is more stringent than 0.02 because samples with a missinigness rate of 0.015 are removed with our threshold.
        #Hardy-Weinberg Equilibrium P >1 × 10−6 (--hwe 1e-6 midp)
        #heterozygosity
            #In our case, just 13 samples are above the mean hetero plus 3 standard deviations, while 42 are below mean hetero less 3 standard deviations. It seems we are ok because in both cases, the outliers are close to the limits (see below). At least closer compared to the plot of Ritchie´s review. In that figure you can see the limits are 0.16-0.18 but there are samples at 0.22 and 0.12.
        #MAF > 0.05 (--maf 0.05)
        #Imputation score >= 0.95.
            #This meets Doug and O´Reilly´s recommendations.
        #Mismatching SNPs due to strand
        #Duplicate SNPs
        #Sex incosistences
        #Population structure and sample relatedness

#PRS calculation:
    #There are several options in terms of how PRSs are calculated. GWASs are performed on finite samples drawn from particular subsets of the human population, and so the SNP effect size estimates are some combination of true effect and stochastic variation— producing ‘winner’s curse’ (see overfitting) among the top- ranking associations—and the estimated effects may not gen- eralize well to different populations (discussed below). The aggregation of SNP effects across the genome is also compli- cated by the correlation between SNPs—LD. Thus, key factors in the development of methods for calculating PRSs are: (i) the potential adjustment of GWAS estimated effect sizes via, for example, shrinkage, (ii) the tailoring of PRSs to target popu- lations and (iii) the task of accounting for LD
    #Shrinkage (reduction) of GWAS effect size estimates
        #PRS methods that perform shrinkage of all SNPs 19,20,45,46 generally exploit commonly used statistical shrinkage/ regularization techniques, such as LASSO or ridge regres- sion 19 , or Bayesian approaches that perform shrinkage via prior distribution specification20,45,46. Under different approaches or parameter settings, varying forms of shrinkage can be achieved: e.g., LASSO regression reduces small effects to zero, while ridge regression shrinks the largest effects more than LASSO but does not reduce any effects to zero. The most appropriate shrinkage to apply is dependent on the underlying mixture of null and true effect size distributions, which are probably a complex mixture of distributions that vary by trait. Since the optimal shrinkage parameters are unknown a priori, PRS prediction is typically optimized across a range of possible parameter values (see below for overfitting issues relating to this), which in the case of LDpred, for example, includes a parameter for the fraction of causal variants 45.
            #This is what LDAK does.
        #In the classic PRS calculation method 5,11,13, only those SNPs with a GWAS association P value below a certain threshold (e.g., P < 1 × 10−5 ) are included in the calculation of the PRS, while all other SNPs are excluded. This approach effectively shrinks all excluded SNPs to an effect size estimate of zero and performs no shrinkage on the effect size estimates of those SNPs included. Since the optimal P value threshold is unknown a priori, PRSs are typically calculated over a range of thresholds, association with the target trait is tested for each, and the prediction is optimized accordingly (see Overfitting in PRS-trait association testing). This process is analogous to tuning parameter optimization in the formal shrinkage methods. An alternative way to view this approach is as a parsimonious variable selection method, effectively performing forward selection ordered by GWAS P value, involving block-updates of predictors (SNPs), with size dependent on the increment between P value thresholds. Thus the ‘optimal threshold’ selected is defined as such only within the context of this forward selection process; a PRS computed from another subset of the SNPs could be more predictive of the target trait, but the number of possible subsets of SNPs is too large to feasibly test given that GWAS are based on millions of SNPs.
            #This is what HIIT-predict did.
    #Controlling for LD
        #If genetic association testing is performed using joint models of multiple SNPs47 , then independent genetic effects can be esti- mated despite the presence of LD. However, association tests in GWASs are typically performed one SNP at a time, which, combined with the strong correlation structure across the genome, makes estimating the independent genetic effects (or best proxies of these if not genotyped/imputed) extremely challenging. If independent effects were estimated in the GWAS or by subsequent fine-mapping, then PRS calculation can be a simple summation of those effects. If, instead, the investigator is using a GWAS based on one-SNP-at-a-time testing, then there are two main options for approximating the PRS that would be obtained from independent effect estimates: (i) SNPs are clumped (i.e., thinned, prioritizing SNPs at the locus with the smallest GWAS P value) so that the retained SNPs are largely independent of each other, and, thus, their effects can be summed, assuming additivity; and (ii) all SNPs are included, accounting for the LD between them. In the classic PRS cal- culation method 5,11,13 , option (i) is combined with P value thresholding and called the C+T (clumping + thresholding) method, while option (ii) is generally favored in methods that implement traditional shrinkage techniques 19,20,45,46 . The rela- tively similar performance of the classic approach to more sophisticated methods14,19,20 may be due to the clumping process capturing conditionally independent effects well; note that clumping does not merely thin SNPs by LD at random (like pruning) but preferentially selects SNPs most associated with the trait under study, and retains multiple SNPs in the same genomic region if there are multiple independent effects there: clumping does not simply retain only the most-associated SNP in a region. A criticism of clumping, however, is that researchers typically select an arbitrarily chosen correlation threshold 41 for the removal of SNPs in LD, and so while no strategy is without arbitrary features, this may be an area for future development of this approach. The key benefits of the classic PRS method are that it is relatively fast to apply and is more interpretable than present alternatives. Both clumping and LD modeling require estimation of the LD between SNPs.
        #HIIT-PREDICT used thresholding + clumping selecting SNPS by p-value and then selecting the most important (lead) SNPs in each group of correlated SNPs.

#The problem of calculating heritability in small datasets (i.e., <5000 samples)
    #from LDAK webpage
        #https://dougspeed.com/small-datasets/
        #To perform heritability analysis (e.g., estimate SNP heritability or partition heritability), requires that the samples are "unrelated" (in practice, this means at most distantly related, with no pair closer than, say, second cousins). This ensures that the heritability estimates reflect only the heritability contributed by predictors in the dataset (or predictors in local linkage disequilibrium with these). By contrast, if your dataset includes (substantial) relatedness, there will likely be long-range linkage disequilibrium (e.g., between predictors on different chromosomes), and you will end up with inflated estimates.
            #we removed first and second degree, i.e., parent-child, siblings, grandparents, grandchildren, aunts, uncles, nieces, nephews, and half-siblings...
            #we did not remove cousins, and we could have them. We used the threshold frequently used in GWAS studies. As you will see below, we are going to be limited for the calculation of heritability due to sample size.
        #Furthermore, heritability analyses generally require a large sample size. For example, to reliably estimate SNP heritability (standard deviation less than 5%) using a single kinship matrix typically needs at least 7,000 unrelated samples; if you wish to use multiple kinship matrices (i.e., perform Genomic Partitioning), the required number of samples is even higher.
    #According to Doug, the clumping and thresholding method technically assumes NO heritability model (it is classical - only bayesian methods assume a heritability model). In contrast, we will likely get noisy estimates of h2, making it hard to leverage the power of good heritability models when using 1.2K samples. So, yes, 1200 samples is typically considered small (if the trait was very simple, eg Mendelian, 1200 would be sufficient; but assuming the trait is complex, it will be challenging to get prediction with 1200 samples). Note, however, that in the same way, the sample size will affect all PRS tools. 
    #More important, in general, individual-level data tools (e.g., .LDAK) are always better than summary statistic tools (e.g., clumping and thresholding). So even for small sample sizes, it remains better to use individual-level data tools. While the main advantage is accuracy, ind-level methods also have the advantage of reporting measures of accuracy (e.g., LDKA uses internal CV, and will report accuracy from this).
        #Note that calculating p-values/effect sizes and then use them as input for the PRS is basically use to summary statistics, instead of taking the statistics from a previous paper, you calculate them yourself. 
        #The alternative is to use individual-level data that, according to Doug, have more accuracy than the summary statistic approach.
    #So if we have an acceptable training-evaluations scheme, we can just use the individual-level approach and then check the accuracy that should be higher anyways than traditional approaches.
        #Doug: Lastly, if doing CV, make sure the training and test samples are quasi independent (normally people with 1200 samples are analyzing animal or plant datasets, with high relatedness; however, it sounds like you have unrelated humans, so all ok)
        #Diego: Regarding the independence of the training and test set, yes, this was a study that included unrelated humans, at least on paper. Anyways, I checked for cryptic relatedness using KING and removed first and second-degree related pairs. Still, the datasets will not be completely independent as they come from the same study, but I think this is the best I can do given the particularities of my cohort.
        #Doug: Regarding the dataset, yes, that sounds absolutely fine (even if you had not removed the first / second degree relatives, I would still think it would be fine for cross validation)
        #Indeed, in their 2021 Nat Comm they use training and test from the same cohort! i.e., the UK Biobank
    #It is possible the small sample size will "break" LDAK, but this will be obvious from the screen output (e.g., if it is impossible to get a meaningful estimate of total heritability) - and if this does happen, it would be interesting to see.
    #Therefore, given the data I have (low sample size and complex trait), individual-level data is still a better option despite the limitations in the calculation of heritability. I going definitively to give a try to --elastic from LDAK and update here if the tools fails because of the sample size.

#LDAK is toolbox that includes many tools and, in particular, the calculation of polygenic scores for trait prediction. This has been developed during years by Doug Speed, a reseracher at Denmark that has published his new methods in journals like Genome Research or Nature Communications, having multiple citations. His methods are indeed cited by the PRS tutorial of O´Really.
        #https://dougspeed.com/
    #It takes care of the Shrinkage of GWAS effect size estimates by applying the srhinkage across all SNPs. As we will sue elastic net, I assume that this will be a combination of lasso (small weights go to zero) and ridge (no weight goes to zero while big weights are reduced).
    #But the great adventage over other methods is that it does not assume that all SNPs have the same influence on heritability, instead the impact of each individual SNPs is estimated consider "minor allele frequency (MAF), local levels of linkage disequilibrium and functional annotations".
        #We are not taking full advantage of this because of our sample size, but it is ok (see above).
    #Across +200 phenotypes, they found that with this new approach the proportion of phenotypic variance explained increases by on average 14%, which is equivalent to increasing the sample size by a quarter.
        #https://www.nature.com/articles/s41467-021-24485-y
    #The input are just plink fileset with the individual genotyping data (although you can also use summary statistics) and phenotype data (response and covariables) using plink format. 
    #This outputs: 1) the heritability estimate of trait. This can be used as reference to see how much of that heriability is explained by the PRS; 2) The predictive model is saved in a .effect file that can be then used as input for "Calculate Scores", another tool of LDAK that predict the trait in a new set samples.
    #The approach they use is to split the samples in two groups in training and test. The PRS is calculated in the training set using a CV approach so the hyperparameters of the model (prior distribution parameters) can be tuned. Then, the model is used to predict the trait in the test set and the correlation between predicted and obserbed is used to evaluate the predicitive power.

#The heritability model 
    #The first heritability analyses used the GCTA Model, which assumes that all SNP contributes equally. Since then, many different models have been developed, that try to more accurately describe how heritability varies across the genome. For example, we proposed the LDAK Model, where the expected heritability contributed by a SNP depends on its minor allele frequency (MAF) and local levels of linkage disequilibrium (LD).
        #http://dougspeed.com/technical-details/
    #The LDAK weightings are designed to account for the fact that levels of linkage disequilibrium vary across the genome. The LDAK weightings are designed to equalize the tagging of SNPs across the genome; SNPs in regions of high linkage disquilbirum (LD) will tend to get low weightings, and vice versa. This will influence the impact of each SNP on the overall heredability (see below for the heirtability formula below). If you have 4 SNPs very close and in linakge they all are considered to tag the same underlying variation having then the same weight than a single SNP that is tagging a single source of underlying variation. HOWEVER, they only advise using them when constructing the BLD-LDAK or BLD-LDAK+Alpha Models, our recommended Heritability Models when analysing summary statistics.
        #https://dougspeed.com/method-overview/
    #We have spent much time investigating the best Heritability Model. We now generally recommend using the Human Default Model. Although there are more realistic heritability models (e.g., those that take into account functional annotations), we recommend the Human Default Model because it has robust performance and is easy to use.
    #We now recommend using the Human Default Model, in which the expected heritability of a SNP depends only on its MAF. We first proposed this model in our paper LDAK-GBAT: fast and powerful gene-based association testing using summary statistics (American Journal of Human Genetics). Although there are more realistic heritability models (e.g., those that take into account functional annotations), we recommend the Human Default Model because it has robust performance and is easy to use.

#additional steps besides QC?
    #A quick follow-up question: If I understand correctly from LDAK documentation, after applying quality control procedures (imputation quality [R2], MAF, sample and SNP missingness, population structure, sample relatedness, sex inconsistency, etc....) I would just need to use the Plink filesets as input for --elastic along with the phenotype and covariates. Then I could use the .effect file as input for --calc-scores in the test set and also assess model performance with --jackknife. So my question, which is probably obvious given what you just explained but I just want to be sure, I do not need to do additional steps besides the QC (i.e., no clumping, thresholding etc...) as I am not using the classical approach, right?
    #Yes, your steps sound right - perform QC, then use the PLINK files with --elastic (adding --LOCO NO, because you only care about prediction) and then --calc-scores

#you could also use LDAK to create the regular approach
    #You may wish to do a classical PRS just for interest / comparison, in which case, you can run --linear (on the training samples), then the file with suffix .score gives classical PRS .corresponding to 7 p-value thresholds (that you can then provide to --calc-scores)
    #NOT SURE IF WE HAVE TO DO CLUMPING, ASK IN THE FUTURE TO DOUG IN CASE NEEDED. 

# endregion






#############################################################
# region define function to calculate PRS across iterations #
#############################################################
#iter_number=1; response_variable="distance_change"; covariate_dataset="small_set_predictors"
def prs_calc(iter_number, response_variable, covariate_dataset):

    print_text(f"For phenotype {response_variable}, and the {covariate_dataset}, starting iteration {iter_number}", header=1)
    print_text("split training and test", header=2)
    print_text("create folders for results", header=3)
    run_bash(" \
        mkdir -p \
            ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "; \
    ")

    print_text("load the phenotype data", header=3)
    pheno_subset = pd.read_csv( \
        "./data/pheno_data/" + covariate_dataset + "/" + response_variable + "_subset/" + response_variable + "_subset.tsv", \
        header=0, \
        sep="\t" \
    )
    print(pheno_subset)

    print_text("specify the covariates", header=3)
    selected_covariates = pheno_subset.columns[~pheno_subset.columns.isin(["family_id", "AGRF code", response_variable])]
    print(selected_covariates)
    #check we have correct covariates
    total_list_covariates = ["Age", "sex_code", "Week 1 Body Mass", "Week 1 Beep Test", "Week 1 Distance (m)", "Week 1 Pred VO2max"] + [f"PCA{i}" for i in range(1,21)]
    if(sum([1 for cov in selected_covariates if cov not in total_list_covariates])!=0):
        raise ValueError("ERROR: FALSE! WE HAVE COVARIATES THAT ARE NOT IN THE TOTAL LIST")

    print_text("transform the phenotypes", header=3)
    #We applied in the previous step the post-imputation QC for the last time after we have removed samples for the last time, this time due to missing phenotype data. These samples without weight are not going to be used anymore, so we should check how the genetic data changes.
    #This is different from the training-evaluation split, where we will train the models with a subset of the samples, but the rest of samples will be used in evaluation, we are not technically remove them, so we should not repeat the QC for each training-eval dataset. If we would do that, we would end up with different genetic predictors (i.e., SNPs) between training-evaluation partitions, and we do not want that.
    #To limit data leakage, we are going to apply the transformation of the phenotypes separataley in trainining and test. Remember what we did in the niche paper, we preprocessed the occurrences, selected the predictors and then we split in training and evaluation just before modeling. This is equivalent with what we have done here.

    print_text("create a copy to save transformed variables", header=4)
    pheno_subset_transform = pheno_subset.copy(deep=True)
        #deep=True: This creates a deep copy of the DataFrame. A deep copy means that a new DataFrame object is created, and all the data is copied. Changes to the new DataFrame will not affect the original DataFrame, and vice versa.

    #select the variable to be transformed
    variables_to_transform = [cov for cov in selected_covariates if cov != "sex_code"] + [response_variable]
    #all covariates (except sex_code) and the response variable
    #Binary phenotypes should only take values 0 (control), 1 (case) or NA (missing). So do NOT transform!
        #http://dougspeed.com/phenotypes-and-covariates/

    #transform all the continuous variables and overwrite
    #pheno_to_transform=variables_to_transform[0]
    for pheno_to_transform in variables_to_transform:
        pheno_subset_transform[pheno_to_transform] = quantile_transform( \
            X=pheno_subset_transform[pheno_to_transform].values.reshape(-1, 1), \
            n_quantiles=int(pheno_subset_transform.shape[0]*0.5), 
            output_distribution="normal" \
        )
            #we select half of the samples to calculate the quantiles. See 03a_phenotype_prep.py for details.
            #Note that all variables numeric and continuous except sex
                #sex_code is numeric but categorical so it should not be transformed!
                #beep test is numeric and continuous, it has float values
            #we use quantile_transform, which is the quivalent function without the estimator API. If you need to invert the transformation, use the transformer: transformer.inverse_transform()
                #https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.QuantileTransformer.html
            #see 03a_phenotype_prep.py for details about the transformation
        #check
        if (\
            (pheno_subset_transform[pheno_to_transform].max() > 6) | \
            (pheno_subset_transform[pheno_to_transform].min() < -6) \
        ):
            raise ValueError(f"ERROR: FALSE! {pheno_to_transform} is not transformed correctly")
    
    #check sex_code has not been transformed
    if (pheno_subset_transform["sex_code"].unique().size != 2):
        raise ValueError("ERROR: FALSE! WE HAVE TRANSFORMED SEX_CODE")

    print_text("save the full transformed dataset", header=3)
    if(iter_number==1):
        run_bash(" \
            mkdir \
                -p \
                ./data/full_set/" + covariate_dataset + "/" + response_variable + " \
                ./results/final_results/" + covariate_dataset + "/" + response_variable + "/full_dataset/; \
        ")
        pheno_subset_transform.to_csv( \
            "./data/full_set/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )

    print_text("create dict to change names of covariates that are problematic for LDAK", header=3)
    dict_change_names={
        "family_id": "FID",
        "AGRF code": "IID",
        "Age": "age",
        "Week 1 Body Mass": "week_1_weight",
        "Week 1 Beep Test": "week_1_beep",
        "Week 1 Distance (m)": "week_1_distance",
        "Week 1 Pred VO2max": "week_1_vo2"
    }

    print_text("split", header=3)
    train_df, test_df = train_test_split( \
        pheno_subset_transform, \
        test_size=0.25, \
        random_state=iter_number \
    )
        #25% for the test set, using as random state the iteration number, so each iteration will have a different training-test set

    print_text("define function to prepare LDAK inputs", header=2)
    #type_df="training"
    #type_df="test"
    def ldak_input_prep(type_df):

        print_text("select the input datasets, training, set or full", header=3)
        if (type_df=="training"):
            input_df = train_df
        elif (type_df=="test"):
            input_df = test_df

        print_text("open a folder for the set", header=3)
        run_bash(" \
            mkdir -p ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/ \
        ")

        print_text("save the input dataframe, test or training", header=3)
        input_df.to_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_transform.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )

        #save the response variable after changing the names
        input_df[["family_id", "AGRF code", response_variable]].rename(columns=dict_change_names).to_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_transform_subset_response.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )

        #save the covariates that are factors
        if ("sex_code" in selected_covariates):

            #save sex_code
            input_df[["family_id", "AGRF code", "sex_code"]].rename(columns=dict_change_names).to_csv( \
                "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_transform_subset_covars_factors.tsv", \
                sep="\t", \
                header=True, \
                index=False, \
                na_rep="NA" \
            )

        #save the covariates that are continuous
        selected_covariates_cont = [cov for cov in selected_covariates if cov != "sex_code"]
        input_df[["family_id", "AGRF code"] + selected_covariates_cont].rename(columns=dict_change_names).to_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_transform_subset_covars_cont.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )

        #create a file with the samples of the selected set
        input_df[["family_id", "AGRF code"]].to_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_samples_in.tsv", \
            sep="\t", \
            header=False, \
            index=False, \
            na_rep="NA" \
        )

        #use that file to select the corresponding sample from the plink fileset
        run_bash(" \
            plink \
                --bfile ./data/plink_filesets/" + covariate_dataset + "/" + response_variable + "_filesets/" + response_variable + "_subset_missing_clean_maf_hwe_sample_snp_missing \
                --keep ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_samples_in.tsv \
                --make-bed \
                --out ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_plink_fileset \
        ")
            #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.

        #Note abut QC:
            #We are NOT doing anymore QC.
            #We did the last one when samples were removed due to NAs in weight data and other covariates. We ensured that for each phenotype, after removing all NAs for its selected covariates, the remaining samples have the correct missingness along with maf, HWE and SNP missingness
            #If we do the QC filtering inside each training/test, i.e., after splitting the dataset, we could end up with different SNPs across training sets!!
            #This is just like we did for the niche paper, preprocessing of the data together, and then split into training and evaluation, not doing anything else to the data

        #check we have selected the correct samples
        subset_fam = pd.read_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_plink_fileset.fam", \
            sep=" ", \
            header=None \
        )
        #if the family_id or the sample IDs are not identical between our transformed dataset and the plink fileset, we have a problem
        input_df_transform_sorted = input_df.sort_values(by=["family_id", "AGRF code"]).reset_index(drop=True)
            #sort the input df to have the same order than the FAM file as it seems that plink reorder columns when preparing the output
        if ( \
            (not subset_fam[0].rename("family_id").equals(input_df_transform_sorted["family_id"])) | \
            (not subset_fam[1].rename("AGRF code").equals(input_df_transform_sorted["AGRF code"])) \
        ):
            raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM SELECT THE SAMPLES WITHOUT NA IN THE PLINK FILESETS")


    print_text("prepare LDAK inputs", header=2)
    ldak_input_prep(type_df="training")
    ldak_input_prep(type_df="test")

    print_text("run LDAK on the training set", header=2)
    run_bash(" \
        mkdir \
            -p ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/; \
        ldak6.1.linux \
          --elastic ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_elastic \
          --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_plink_fileset \
          --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_response.tsv \
          --covar ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_cont.tsv \
          --factors ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_factors.tsv \
          --LOCO NO \
    ")
    #ELASTIC:
        #https://dougspeed.com/elastic-net/


    ###use -mpheno!!!! to select the pheno you want
        #http://dougspeed.com/phenotypes-and-covariates/

    ###NO USES MPHENO!!! UN FILE POR PHENO ASI NO TE HACE LA MEDIA
    ##Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). Binary phenotypes should only take values 0 (control), 1 (case) or NA (missing). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype).

    run_bash(" \
        mkdir \
            -p ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/; \
        ldak6.1.linux \
            --calc-scores ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation \
            --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_plink_fileset \
            --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_transform_subset_response.tsv \
            --scorefile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_elastic.effects \
            --power 0 \
    ")

    #CALCULATE SCORES
        #https://dougspeed.com/profile-scores/

    run_bash(" \
        ldak6.1.linux \
            --jackknife ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_jacknife_eval \
            --profile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation.profile \
            --num-blocks 200 \
    ")

    #jacknife
        #https://dougspeed.com/jackknife/


    #You may wish to do a classical PRS just for interest / comparison, in which case, you can run --linear (on the training samples), then the file with suffix .score gives classical PRS .corresponding to 7 p-value thresholds (that you can then provide to --calc-scores)
    #get also from here the manhtan plot? we need to check inflation? maybe just bonferroni and show no snp is signifncat, so PRSs are useful

        #--linear: https://dougspeed.com/single-predictor-analysis/

    run_bash(" \
        ldak6.1.linux \
            --linear ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_linear \
            --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_plink_fileset \
            --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_response.tsv \
            --covar ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_cont.tsv \
            --factors ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_factors.tsv \
            --permute YES \
    ")

    run_bash(" \
        mkdir \
            -p ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/; \
        ldak6.1.linux \
            --calc-scores ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation_classical \
            --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_plink_fileset \
            --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_transform_subset_response.tsv \
            --scorefile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_linear.score \
            --power 0 \
    ")

    run_bash(" \
        ldak6.1.linux \
            --jackknife ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_classical_jacknife_eval \
            --profile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation_classical.profile \
            --num-blocks 200 \
    ")



    #Error, ./data/train_test_sets/distance_change/train_test_iter_1/training_set/distance_change_training_set_transform_subset_response.tsv contains multiple phenotypes, so you must specify one using "--mpheno", or use "--mpheno ALL" to analyse all phenotypes

        #Yes, your steps sound right - perform QC, then use the PLINK files with --elastic (adding --LOCO NO, because you only care about prediction) and then --calc-scores
        #Select the Heritability model
            #The GCTA Model assumes E[h2j] = tau1 (i.e., that expected heritability is constant across SNPs). Note that this is the model assumed by any method that first standardizes SNPs, then assigns to each the same prior distribution or penalty function. In LDAK, this model is specified by adding the options --ignore-weights YES and --power -1 when Calculating Kinships or Calculating Taggings.
            #The Human Default Model assumes E[h2j] = tau1 [fj(1-fj)]0.75, where fj is the minor allele frequency of SNP j. Therefore, the Human Default Model assumes that more common SNPs have higher E[h2j] than less common SNPs. In LDAK, this model is specified by adding the options --ignore-weights YES and --power -0.25 when Calculating Kinships or Calculating Taggings.
                #THIS IS THE RECOMMENDED MODEL IN GENERAL BY DOUG
            #The LDAK Model assumes E[h2j] = tau1 wj [fj(1-fj)]0.75, where wj are the LDAK weightings; these will tend to be lower for SNPs in regions of high linkage disequilibrium (LD), and vice versa. Therefore, the LDAK Model assumes that E[h2j] is higher for SNPs in regions of lower LD and for those with higher MAF. In LDAK, this model is specified by adding the options --weights <weightsfile> and --power -0.25 when Calculating Kinships or Calculating Taggings, where <weightsfile> provides the LDAK Weightings.
                #Here is where the weigthings enters in action, considering not only the MAF but also the LD level in each genomic region.





    #you can use --elasti en el 75%, ahí te hace automaticamente CV para seleccionar hiperparamteros usando plink filset como inputs. El output se puede usar como input para otra funcion que calcula los scores para nuevos individuaos, el 25% restante. Do function that do this for a given iteration so you can use 10 cores to do 100 iterations in 1 day. you can also calculate the regulra PRS in each iteration. at the end of paralelization, in a different script, you will combine R2, accurcy of all iteratons and calculate media and CI. Caclulate the median difference respect to the original PRS approach. you could also calculate the degree of overlap across iterations for the top 10% responders. as a previous step in a different step, use LDAK to do manhatan plot (LDAK-KVIK).



    #if you compare with trad-PRS, you can use log?
        #https://dougspeed.com/compare-models/

    #David wants to show also manhatann plots (maybe just use plink or LDAK?), also repeat 10 times the 75-25% and see if the respodners are the same?


    ##use simple PRS (linear) as reference as Dr. Speed suggested in your github conversation?

    #you could take from here the p-values to calculate the manhattan plot and the p-values could be used later for the BAT analyses



#run the function across the three phenotypes
#ALSO CONSDIEIRNG REDUCED AND FULL COVARIATE SET
prs_calc(iter_number, response_variable, covariate_dataset, manhattan_plot)





#REMOVE PLINK FILES AFTER, IF NOT CRAZY AMOUNT OF SPACE


#MAKE A BASELINE regular approach?





###PUT THIS IN A BASH SCRIPT INSIDE LDAK FOLDER
"""
run_bash(" \
    mkdir -p ../ldak_versions/; \
    cd ../ldak_versions/; \
    wget --output-document ldak6.1.linux https://github.com/dougspeed/LDAK/raw/refs/heads/main/ldak6.1.linux; \
    chmod +x ./ldak6.1.linux; \
")
    #http://dougspeed.com/downloads2/
    #https://github.com/dougspeed/LDAK
"""


##WHEN INTERPRETING THE RESULTS OF THE PRS, LOOK SLIDES FROM DOUG
    #https://dougspeed.com/short/
###WHEN DONE, YOU CAN CHECK THE SECTION OF INTERPRETATION FROM ORRELLY
    #Interpretation and presentation of results



# endregion



##parallelize
##PARSE ITERATION AS A COMMAND SO YOU CAN RUN EACH ITERATION SEPARATETLY FOR EACH PHENO
#WE COULD SEND 100 JOBS OF JUST 1 CORE AND 1 HOUR AND THEY SHOULD BE RUNNING FAST
#ALSO INCLUDE THE PHENO AS A ARGUMENT, SO WE CAN DO THIS FOR EACH OF THE THREE PHENOTYPES





#FOR BAT ANALYSES
    #this would an additional step in this project that would be outside of the paper
    #take the 1000kb gene windows for all coding genes, liftover to hg38. If the USCS tool accepts genomic ranges, just use them as input, if not, split in two datasets the start and the end of the gene windows
    #for each phenotype (VO2, beep....), calculate the average (better than median because want influence of outliers within gene like in iHS, if a SNPs is veery important in a gene that should influence the info about the whole gene) p-value for the association of SNPs inside each gene
    #then, calculate 1000 random sets of genes, within each set, calculate the median association of all genes inside the set and compare with the BAT set to obtain a distribution and empirical p-value (is association lower in BAT? LOOF BAT PAPER SCRIPTS FOR THIS). Here I want median because i do not want a gene outliser change things, I want the overall impact of BAT genes in general, not just a few genes.





#Data and validation scheme
    #Integration of risk factor polygenic risk score with disease polygenic risk score for disease prediction
    #Estimating Disorder Probability Based on Polygenic Prediction Using the BPC Approach
    #Exploring the application of deep learning methods for polygenic risk score estimation
    #Machine Learning Models of Polygenic Risk for Enhanced Prediction of Alzheimer Disease Endophenotypes
    #Genome-wide association studies and polygenic risk score phenome-wide association studies across complex phenotypes in the human phenotype project
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




''''
##you can use the following code to create manhatan plots.

print_text("run the associations between SNPs and phenotypes", header=3)
print_text("make dict for title plots", header=4)
dict_titles = {
    "weight_change": "Change of body mass", \
    "beep_change": "Change in beep test", \
    "distance_change": "Change in distance", \
    "vo2_change": "Change in predicted VO2 max"}
print(dict_titles)

print_text("run the associations and make the Manhattan plots", header=4)
#pheno="distance_change"
for pheno in ["weight_change", "beep_change", "distance_change", "vo2_change"]:
    print("\n \n" + pheno)

    #get the covariate names
    covs = dict_pvalue_covar_final[pheno][0]
    sex_cov = dict_pvalue_covar_final[pheno][1]

    print("run plink2 assoc")
    #we use plink2 because it can directly preprocess phenotypes and covaraites using the quantile transformation (see above for results of the transformation in the phenotypes)
    run_bash(" \
        mkdir -p ./results/prelim_results/" + pheno + "/; \
        plink2 \
            --bfile .././quality_control/data/genetic_data/quality_control/08_loop_maf_missing/loop_maf_missing_2 \
            --linear " + sex_cov + " \
            --quantile-normalize \
            --input-missing-phenotype " + str(new_nan_value) + " \
            --pheno ./data/plink_inputs/pheno_file.tsv \
            --pheno-name " + pheno + "\
            --covar ./data/plink_inputs/covar_file.tsv \
            --covar-name " + covs + " \
            --out ./results/prelim_results/" + pheno + "/first_assoc; \
        ls -l ./results/prelim_results/" + pheno + "/")
            #Given a quantitative phenotype and possibly some covariates (in a --covar file), --linear writes a linear regression report to pheno_name.glm.linear.
                #sex
                    #First, sex (as defined in the .fam/.psam input file) is normally included as an additional covariate. If you don't want this, add the 'no-x-sex' modifier. Or you can add the 'sex' modifier to include .fam/.psam sex as a covariate everywhere. WHATEVER YOU DO, DON'T INCLUDE SEX FROM THE .FAM/.PSAM FILE AND THE --COVAR FILE AT THE SAME TIME; OTHERWISE THE DUPLICATED COLUMN WILL CAUSE THE REGRESSION TO FAIL
                #quantile-normalize
                    #--quantile-normalize forces named quantitative phenotypes and covariates to a N(0, 1) distribution, PRESERVING ONLY THE ORIGINAL RANK ORDERS; if no parameters are provided, all quantitative phenotypes and covariates are affected. --pheno-quantile-normalize does the same for just quantitative phenotypes, while --covar-quantile-normalize does this for just quantitative covariates.
                #input--missing-phenotype
                    #Missing case/control or quantitative phenotypes are expected to be encoded as 'NA'/'nan' (any capitalization) or -9. By default, other strings which don't start with a number are now interpreted as categorical phenotype/covariate values; to force them to be interpreted as missing numeric values instead, use --no-categorical.
                    #You can change the numeric missing phenotype code to another integer with --input-missing-phenotype, or just disable -9 with --no-input-missing-phenotype.
                #See above for details about --pheno/--covar and the corresponding names
                #https://www.cog-genomics.org/plink/2.0/assoc
                #https://www.cog-genomics.org/plink/2.0/assoc#sex
                #https://www.cog-genomics.org/plink/2.0/data#quantile_normalize
                #https://www.cog-genomics.org/plink/2.0/input

    print("IMPORTANT: \
        There are multiple disclamers in the --glm info page, like problems with low minor allele count or how to prune by relatedness, include PCAs for population stratification, etc... \
        TAKE A LOOK AT THIS!!")
        #https://www.cog-genomics.org/plink/2.0/assoc

    print("IMPORTANT: \
        It seems that SEX is not subjected to quantile transformation, only the covariates included in the covar file. \
        Be careful in the future if you add discrete covariates in the covar file, like batch. \
        Also you should check how the distribution of the covariates change with the quantile normal transformation, you only did it for phenotypes.")

    print("use awk to see the assoc file. Then use it select only rows for ADD effect of the SNP and the first row, then save as tsv")
    run_bash(" \
        cd ./results/prelim_results/" + pheno + "/; \
        echo 'See 7 first columns of the original assoc file'; \
        awk \
            'BEGIN{FS=\" \"}{if(NR<=10){print $1, $2, $3, $4, $5, $6, $7}}END{print \"The number of columns is \"; print NF}' \
            first_assoc." + pheno + ".glm.linear; \
        echo '\nProcess with awk to specify the delimiter and then take a look to the file'; \
        awk \
            'BEGIN{FS=\" \"; OFS=\"\t\"}{if($7==\"ADD\" || NR==1){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}}'\
            first_assoc." + pheno + ".glm.linear | \
        gzip \
            --stdout > first_assoc." + pheno + ".glm.linear.tsv.gz; \
        gunzip \
            --stdout \
            first_assoc." + pheno + ".glm.linear.tsv.gz | \
        awk \
            'BEGIN{FS=\"\t\"}{if(NR<=10){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}}'")

    #extract the unique cases in the TEST column in order to check if we have selected only ADD
    unique_test = run_bash(" \
        cd ./results/prelim_results/" + pheno + "/; \
        gunzip \
            --stdout \
            first_assoc." + pheno + ".glm.linear.tsv.gz | \
        awk \
            'BEGIN{FS=\"\t\"}{if(NR>1 && !a[$7]++){print $7}}'", return_value=True).strip()
        #this is the way to get uniq cases from a colum in awk
            #if the value in $7 has not been previously included in "a" before, get True (without the "!" you would get false). Then, print the value of $7 in those cases, i.e., the unique value of the TEST column (7th column)
            #https://stackoverflow.com/questions/1915636/is-there-a-way-to-uniq-by-column
    print("Do we have correctly selected the ADD rows only (i.e., additive effect of SNPs)?")
    if (type(unique_test)==str) & (unique_test=="ADD"):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM FILTERING THE ROWS WITH THE ADDITIVE EFFECT OF THE SNPS")

'''

print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/association_studies/
#chmod +x ./scripts/production_scripts/03b_prs_calculation.py
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03b_prs_calculation.py --iter_number=1 --response_variable=distance_change --covariate_dataset="large_set_predictors" > ./03b_prs_calculation_iter_1_distance_change.out 2>&1
    #We should run the 8 dataset_size and phenotype combinations separtaely in the cluster, each combination run in serial 100 iterations, alwaysing requiering just 20GB, each iteration is around 10 minutes, so in 16-20 hours should be done.
#grep -Ei 'error|false|fail' ./03b_prs_calculation.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.
