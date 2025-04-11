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

#This script will calculate the Polygenic Risk Scores. 

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
parser.add_argument("--total_iter", type=int, default=100, help="Total number of iterations. Integer always, None does not work!")
parser.add_argument("--response_variable", type=str, default="distance_change", help="Phenotype to model. String always, None does not work!")
parser.add_argument("--covariate_dataset", type=str, default="small_set_predictors", help="Type of dataset used, short or long list of covariates. String always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
total_iter = args.total_iter
response_variable = args.response_variable
covariate_dataset = args.covariate_dataset

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
            #we did not remove cousins, and we could have them. We used the threshold frequently used in GWAS studies. As you will see below, we are going to be limited for the calculation of heritability due to sample size anyways.
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
        #The GCTA Model assumes E[h2j] = tau1 (i.e., that expected heritability is constant across SNPs). Note that this is the model assumed by any method that first standardizes SNPs, then assigns to each the same prior distribution or penalty function.
        #http://dougspeed.com/technical-details/
    #The LDAK weightings are designed to account for the fact that levels of linkage disequilibrium vary across the genome. The LDAK weightings are designed to equalize the tagging of SNPs across the genome; SNPs in regions of high linkage disquilbirum (LD) will tend to get low weightings, and vice versa. This will influence the impact of each SNP on the overall heredability (see below for the heirtability formula below). If you have 4 SNPs very close and in linakge they all are considered to tag the same underlying variation having then the same weight than a single SNP that is tagging a single source of underlying variation. HOWEVER, they only advise using them when constructing the BLD-LDAK or BLD-LDAK+Alpha Models, our recommended Heritability Models when analysing summary statistics.
        #The LDAK Model assumes E[h2j] = tau1 wj [fj(1-fj)]0.75, where wj are the LDAK weightings; these will tend to be lower for SNPs in regions of high linkage disequilibrium (LD), and vice versa. Therefore, the LDAK Model assumes that E[h2j] is higher for SNPs in regions of lower LD and for those with higher MAF. Here is where the weigthings enters in action, considering not only the MAF but also the LD level in each genomic region.
        #https://dougspeed.com/method-overview/
    #We have spent much time investigating the best Heritability Model. We now generally recommend using the Human Default Model. Although there are more realistic heritability models (e.g., those that take into account functional annotations), we recommend the Human Default Model because it has robust performance and is easy to use.
        #The Human Default Model assumes E[h2j] = tau1 [fj(1-fj)]0.75, where fj is the minor allele frequency of SNP j. Therefore, the Human Default Model assumes that more common SNPs have higher E[h2j] than less common SNPs.
    #We now recommend using the Human Default Model, in which the expected heritability of a SNP depends only on its MAF. We first proposed this model in our paper LDAK-GBAT: fast and powerful gene-based association testing using summary statistics (American Journal of Human Genetics). Although there are more realistic heritability models (e.g., those that take into account functional annotations), we recommend the Human Default Model because it has robust performance and is easy to use.

#additional steps besides QC?
    #A quick follow-up question: If I understand correctly from LDAK documentation, after applying quality control procedures (imputation quality [R2], MAF, sample and SNP missingness, population structure, sample relatedness, sex inconsistency, etc....) I would just need to use the Plink filesets as input for --elastic along with the phenotype and covariates. Then I could use the .effect file as input for --calc-scores in the test set and also assess model performance with --jackknife. So my question, which is probably obvious given what you just explained but I just want to be sure, I do not need to do additional steps besides the QC (i.e., no clumping, thresholding etc...) as I am not using the classical approach, right?
    #Yes, your steps sound right - perform QC, then use the PLINK files with --elastic (adding --LOCO NO, because you only care about prediction) and then --calc-scores

#you could also use LDAK to create the regular approach
    #You may wish to do a classical PRS just for interest / comparison, in which case, you can run --linear (on the training samples), then the file with suffix .score gives classical PRS .corresponding to 7 p-value thresholds (that you can then provide to --calc-scores)
    #NOT SURE IF WE HAVE TO DO CLUMPING, ASK IN THE FUTURE TO DOUG IN CASE NEEDED. 

#CV scheme and independent dataset
            #In some cases like studies of rare diseases, they use interval validation of the GWAS, i.e., using the same cohort but applying resampling methods. 
                #This is an option but it would be a clear limitation of the study.
                #the gold standard in this point is to use an independent and diverse cohort for validation
            #We also have a good point in the fact we have a very specific cohort. It is very difficult to find a population with similar characteristics. This is just like having a cohort of a rare disease.
            #Data leaking
                #we remove SNPs based on the allele frequencies of all samples, including those of the validation dataset.
                #We also remove samples based on relatedness and PCA considering test and evaluation samples, but you have to related individuals and as you have more data you can do it better...
                #the most problematic point would be to use PCAs as predictors...
            #We have explained all of this to Dr.Speed and he does not seem to be specially concern regarding the independence of the samples for evaluation.
            #info
                #Evaluation of a two-step iterative resampling procedure for internal validation of genome-wide association studies
                    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4859941/
                #In "Genome-wide association studies", nature review methods, they say that internal validation is a possibility, but no the gold standard
                    #https://www.nature.com/articles/s43586-021-00056-9
                #David sent a paper also using this approach
                    #Personalized Nutrition by Prediction of Glycemic Responses

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

    #check we do not have NAs
    if(pheno_subset_transform.isna().sum().sum()!=0):
        raise ValueError("ERROR: FALSE! WE HAVE NAs IN THE TRANSFORMED DATASET")

    print_text("save the full transformed dataset", header=3)
    if(iter_number==1):
        run_bash(" \
            mkdir \
                -p \
                ./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "; \
        ")
        pheno_subset_transform.to_csv( \
            "./data/full_set_transform/" + covariate_dataset + "/" + response_variable + "/" + response_variable + "_full_set_transform.tsv", \
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

        print_text("save the input dataframe for test or training", header=3)
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
            #In LDAK, phenotype files should be in PLINK format. The first two columns should provide sample IDs, with subsequent columns providing values for each phenotype. The files can contain a header row, but then the first two elements must be named either FID & IID or ID1 & ID2.
            #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype). We do not have NAs anyways.
            #Also, we are going to have just 1 response variable, so no problem with filling NA cases with the average. We have just 1 response and it is free of NAs.

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
                #In LDAK, phenotype files should be in PLINK format. The first two columns should provide sample IDs, with subsequent columns providing values for each phenotype. The files can contain a header row, but then the first two elements must be named either FID & IID or ID1 & ID2.
                #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype). We do not have NAs anyways.

        #save the covariates that are continuous
        selected_covariates_cont = [cov for cov in selected_covariates if cov != "sex_code"]
        input_df[["family_id", "AGRF code"] + selected_covariates_cont].rename(columns=dict_change_names).to_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_transform_subset_covars_cont.tsv", \
            sep="\t", \
            header=True, \
            index=False, \
            na_rep="NA" \
        )
            #In LDAK, phenotype files should be in PLINK format. The first two columns should provide sample IDs, with subsequent columns providing values for each phenotype. The files can contain a header row, but then the first two elements must be named either FID & IID or ID1 & ID2.
            #Missing phenotypic values should be denoted by NA (note that while PLINK also treats -9 as missing, this is not the case in LDAK). In general, LDAK excludes samples with missing phenotypes (the exception is when analyzing multiple phenotypes, in which case LDAK generally replaces missing values with the mean value of the corresponding phenotype). We do not have NAs anyways.

        #create a file with the samples of the selected set
        input_df[["family_id", "AGRF code"]].to_csv( \
            "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/" + type_df + "_set/" + response_variable + "_" + type_df + "_set_samples_in.tsv", \
            sep="\t", \
            header=False, \
            index=False, \
            na_rep="NA" \
        )
            #This is the input for plink´s --keep argument. This accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.

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

    print_text("create and evaluate a PRS using elastic net", header=2)
    print_text("run LDAK on the training using elastic", header=3)
    run_bash(" \
        mkdir \
            -p \
            ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/; \
        ldak6.1.linux \
          --elastic ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_elastic \
          --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_plink_fileset \
          --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_response.tsv \
          --covar ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_cont.tsv \
          --factors ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_factors.tsv \
          --LOCO NO \
    ")
        #Dr. Speed: Yes, your steps sound right - perform QC, then use the PLINK files with --elastic (adding --LOCO NO, because you only care about prediction) and then --calc-scores
            #We DO NOT CARE ABOUT HERIABILITY BECAUSE WE CANNOT CALCULATE IT PROPERLY!! JUST PREDICTION.
        #Here we explain how to construct an elastic net prediction model. These instructions assume you are analysing individual-level data.
            #So I understand we are combining here ridge and lasso, which is the elastic net.
        #construct an elastic net prediction model by analysing individual-level data (if instead you are analysing summary statistics, you should use MegaPRS).
            #--bfile/--speed <datastem> or --bgen <datafile> - to specify the genetic data files (see File Formats). If your genetic data are in a different format, you should first Make Data.
                #LDAK accepts genetic data in many formats. When analyzing SNP data, it is easiest to use bed format (binary PLINK) which is designed to store hard-coded SNP genotypes (all values must be 0, 1, 2 or NA).
                #LDAK considers the last two columns of the bim file are the A1 and A2 alleles, respectively. During, QC, we checked for strand issues and flip alleles if needed. Then, we select REF and ALT from the VCF files generated from the imputation and converted them to plink format using plink.
                #Note that plink swaps REF by ALT if the REF is less frequent, but this should not be a problem because we are not going to use the alleles, just PRS. Anyways, if REF is "A" and ALT is "T", as long as we are using the foward strand, then we can compare with other studies, because "T" is indeed "T", so we could see if that allele associate with a given phenotype in the same way than other studies. And we have ensured this during the quality control also removing problematic palindromic cases...
                    #If needed, you can use "--ref-allele" and input the REF allele list from the source VCF files to correct this.
            #--pheno <phenofile> - to specify phenotypes (in PLINK format; see above). Samples without a phenotype will be excluded. If <phenofile> contains more than one phenotype, specify which should be used with --mpheno <integer>, or use --mpheno ALL to analyze all phenotypes.
                #We are using just one phenotype, so we do not need to specify --mpheno
            #You can use --covar <covarfile> or --factors <factorfile> to provide quantitative or categorical covariates (in PLINK format; see above); the phenotype will be regressed on these prior to estimating effect sizes.
                #Note that if a categorical covariate has U unique values, LDAK will (internally) replace it with U-1 indicator variables (LDAK will give an error if the total number of indicator variables is greater than half the sample size).
                #For example, sex code is 1 and 2, so in LDAK would be 0 and 1.
            #--LOCO NO - to tell LDAK to focus on creating the genome-wide prediction model (instead of creating leave-one-chromosome-out models for use with LDAK-KVIK).
                #This was the recommendation of Dr. Speed, and we are not using the LDAK-KVIK model anyways.
            #By default, LDAK will estimate the heritability and the power parameter alpha; to instead specify their values use --power <float> and --her <float> (note that if you use --her, you must also use --power).
            #By default, LDAK will use 90%/10% cross-validation to determine suitable prior distribution parameters. You can change the fraction of test samples uing --cv-proportion <float>,  specify the test samples using --cv-samples <cvsampsfile>, or turn off cross-validation, using --skip-cv YES (LDAK will then output multiple models, each trained using 100% of samples).
                #So it uses CV internally, which is great. We do CV to tune the hyperparameters of the model, and then we use the best model to predict the test set.
            #By default, LDAK will assign all predictors (i.e., SNPs) weighting one (equivalent to using --ignore-weights YES). If you prefer to provide your own weightings, use --weights <weightsfile> or --ind-hers <indhersfile> (note that if using --ind-hers, you can not use --her or --power).
                #we are just going to use the default weights. To my understanding, LDAK will give different weights to the SNPs depending on their MAF, so we do not need to worry about that.
            #You can use --keep <keepfile> and/or --remove <removefile> to restrict to a subset of samples, and --extract <extractfile> and/or --exclude <excludefile> to restrict to a subset of predictors (for more details, see Data Filtering).
            #output:
                #The estimated prediction model is saved in <outfile>.effects. Usually, this file has five columns, providing the predictor name (SNP), its A1 and A2 alleles, the average number of A1 alleles, then its estimated effect (relative to the A1 allele). If you used --skip-cv YES, there will be effect sizes for each of the different prior parameters. This file is ready to be used for Calculating Scores (i.e., to predict the phenotypes of new samples).
            #https://dougspeed.com/elastic-net/

    print_text("calculate the PRS in the test set", header=3)
    run_bash(" \
        mkdir \
            -p \
            ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/; \
        ldak6.1.linux \
            --calc-scores ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation \
            --scorefile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_elastic.effects \
            --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_plink_fileset \
            --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_transform_subset_response.tsv \
            --power -1 \
    ")
        #calculate linear combinations of predictor values (the linear projection of genetic data onto predictor effect sizes). These are most commonly used for creating polygenic risk scores and testing the performance of prediction models. Please note that if you include a phenotype when calculating scores, then the resulting profile file can be used with Jackknife (in order to compute additional measures of accuracy and corresponding estimates of precision).
        #--scorefile <scorefile> - to provide the predictor effect sizes. The score file should have at least five columns. The first 4 columns provide the name of each predictor, its two alleles (test allele followed by reference allele), then its centre (the mean number of test alleles); the remaining columns provide sets of predictor effect sizes. The file should have a header row, whose first element must be "Predictor" or "SNP". Note that if centre is "NA" for a predictor, then LDAK will centre values based on the mean number of test alleles in the genetic data.
        #–bfile/–gen/–sp/–speed <prefix> or --bgen <datafile> - to specify genetic data files (see File Formats)
        #-pheno <phenofile> - to specify phenotypes (in PLINK format; see above).
            #Typically, the sets of effect sizes in the score file will correspond to different prediction models, and we are interested in measuring their accuracy. We have phenotypes for samples in the genetic data (i.e., you have individual-level validation data), in which case you should provide these using --pheno <phenofile>. LDAK will then calculate the correlation between scores and phenotypic values for the samples in the genetic data.
        #Sometimes you will have two sets of prediction models, for example, one computed using training samples and one computed using all samples. You can then use --scorefile <scorefile> to provide the first set of prediction models, --pheno <phenofile> or --summary <sumsfile> to provide phenotypic values or summary statistics, and --final-effects <finaleffectsfile> to provide the second set of prediction models. LDAK will save to <outfile>.effects.best the prediction model from the second set that corresponds to the most accurate model from the first set (for example, if Model 5 from the first set has highest correlation, LDAK will save Model 5 from the second set).
        #--power <float> - to specify how predictors are scaled (see below). Usually, the score file contains raw effects, so you should use --power 0.
            #When power equals one, predictors are divided by their expected standard deviation (i.e., are standardised). Therefore, you should use --power=-1 when the score file contains standardised effect sizes. By contrast, when power equals zero, predictors are no longer scaled. Therefore, you should use --power=0 when the score file contains raw effect sizes.
            #I think our effect sizes coming from the elastic net are standardised, so we should use --power -1, but I have to ask Dr. Speed.
        #output:
            #The profile scores will be saved in <outfile>.profile, the (estimated) correlation between scores and phenotypic values will be saved in <outfile>.cors.
        #https://dougspeed.com/profile-scores/

    print_text("calculate evaluation metrics", header=3)
    run_bash(" \
        ldak6.1.linux \
            --jackknife ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_jacknife_eval \
            --profile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation.profile \
            --num-blocks 200 \
    ")
        #The jackknife function measures the similarity between pairs of vectors containing predicted and observed values. Specifically, it computes the correlation, correlation squared (on the observed scale), mean squared error and mean absolute error, as well as corresponding estimates of standard deviation. The jackknife function was designed for computing the accuracy of polygenic risk scores, such as those created by the Prediction tools.
            #--data-pairs <datapairs> or --profile <profile> - to provide pairs of predicted and observed values.  If using --data-pairs, then <datapairs> should have either two or three columns (with no headers), containing predicted values, observed values, and (if provided) regression weights. If using --profile, then <profile> should be the output from Calculation Scores.
            #--pheno (in PLINK format)
                #Because we included the option --pheno in the previous step (--calc-scores), the resulting profile contains both scores and the phenotypes. Therefore, we can measure how well the scores predict the phenotypes by running
            #--num-blocks <integer> - to specify the number of jackknife blocks (usually 200 is a good choice).
            #If the observed values are binary, you can add --AUC YES, and LDAK will additionally compute the area under curve. You can also add --prevalence <float> to specify the proportion of cases in the population, and then LDAK will also report correlation squared on the liability scale.
            #output:
                #The accuracy estimates will be saved in <outfile>.jack.
        #https://dougspeed.com/jackknife/


    print_text("calculate and evaluate a PRS using the classical approach", header=2)
    #Dr. Speed: You may wish to do a classical PRS just for interest / comparison, in which case, you can run --linear (on the training samples), then the file with suffix .score gives classical PRS .corresponding to 7 p-value thresholds (that you can then provide to --calc-scores)
    print_text("create Classical PRS in training dataset", header=3)
    run_bash(" \
        ldak6.1.linux \
            --linear ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_linear \
            --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_plink_fileset \
            --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_response.tsv \
            --covar ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_cont.tsv \
            --factors ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_transform_subset_covars_factors.tsv \
            --permute NO \
    ")
        #perform one-predictor-at-a-time analysis, using linear regression (classical model, i.e., not mixed model).
            #--bfile/--gen/--sp/--speed <datastem> or --bgen <datafile> - to specify the genetic data files (see File Formats).
            #--pheno <phenofile> - to specify phenotypes (in PLINK format). Samples without a phenotype will be excluded. If <phenofile> contains more than one phenotype, specify which should be used with --mpheno <integer>, or use --mpheno ALL to analyze all phenotypes.
            #You can use --covar <covarfile> or --factors <factorfile> to provide quantitative or categorical covariates (in PLINK format) as fixed effect in the regression.
            #If you add --spa-test YES, LDAK will recompute p-values for the most associated predictors using a SaddlePoint Approximation.
            #It remains possible to  add  --grm <kinfile>, in which case LDAK will perform mixed-model linear regression. However, we instead recommend using LDAK-KVIK, which is usually faster and more powerful.
            #To perform a within-family analysis, add --families YES (for more details of this analysis, see Howe et al.). Note that LDAK will infer families based on the 1st column of the fam file (the FID).
            #To perform a trio analysis, add --trios YES. Note that LDAK will infer trios based on the 3rd and 4th column of the fam file (the PID and MID).
            #To perform weighted linear regression, use --sample-weights <sampleweightfile>. The file <sampleweightfile> should have three columns, where each row provides two sample IDs followed by a positive float. Note that by default, LDAK will use the sandwich estimator of the effect size variance (see this page for an explanation); to instead revert to the standard estimator of variance, add --sandwich NO.
            #If you add --permute YES - the phenotypic values will be shuffled. This is useful if wishing to perform permutation analysis to see the distribution of p-values or test statistics when there is no true signal.
                #SO WHEN YOU SET THIS AS YES, YOU ARE GETTING THE P-VALUE WHEN THERE IS NO SIGNAL
            #output:
                #When performing a standard analysis, LDAK produces five output files: <outfile>.assoc contains the main results; <outfile>.summaries contains summary statistics (in the format required for use with SumHer, and MegaPRS); <outfile>.pvalues contains p-values (useful if you wish to Thin Predictors); <outfile>.coeff contains estimates of the fixed effects; <outfile>.score contains simple prediction models corresponding to six different p-value thresholds.
            #https://dougspeed.com/single-predictor-analysis/

    print_text("calculate the classical PRS in the test set", header=3)
    run_bash(" \
        ldak6.1.linux \
            --calc-scores ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation_classical \
            --scorefile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_linear.score \
            --bfile ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_plink_fileset \
            --pheno ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_transform_subset_response.tsv \
            --power -1 \
    ")
        #see above for the arguments of --calc-scores, just remember
            #--scorefile <scorefile> - to provide the predictor effect sizes. The score file should have at least five columns. The first 4 columns provide the name of each predictor, its two alleles (test allele followed by reference allele), then its centre (the mean number of test alleles); the remaining columns provide sets of predictor effect sizes. In the case of the output of --linear, the fifth columns is ALL, and it seems to have effect sizes.
            #power -1 assumes that the effect sizes of the SNPs in the score file are standardized. 

    print_text("calculate evaluation metrics", header=3)
    run_bash(" \
        ldak6.1.linux \
            --jackknife ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_classical_jacknife_eval \
            --profile ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation_classical.profile \
            --num-blocks 200 \
    ")
        #see above for the arguments of jackknife
        #in the case of --linear, as we have several p-value thresholds, we have several models and hence the metrics for each threshold


    print_text("check we have used the correct samples in both analyses", header=2)
    print_text("load the FAM files used for training and test", header=3)
    fam_file_training = pd.read_csv( \
        "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_set_plink_fileset.fam", \
        sep=" ", \
        header=None \
    )
    fam_file_test = pd.read_csv( \
        "./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_test_set_plink_fileset.fam", \
        sep=" ", \
        header=None \
    )
    
    print_text("samples in files generated by elastic net", header=3)
    prs_elastic = pd.read_csv( \
        "./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_elastic.prs", \
        sep=" ", \
        header=0 \
    )
    profile_prs_elastic = pd.read_csv( \
        "./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation.profile", \
        sep="\t", \
        header=0 \
    )
    if ( \
        (not fam_file_training[1].rename("IID").equals(prs_elastic["IID"])) | \
        (not fam_file_test[1].rename("ID2").equals(profile_prs_elastic["ID2"])) \
    ):
        raise ValueError("ERROR: FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLES FOR ELASTIC NET")

    print_text("samples in files generated by --linear", header=3)
    combined_linear = pd.read_csv( \
        "./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/training_set/" + response_variable + "_training_linear.combined", \
        sep=" ", \
        header=0 \
    )
    profile_prs_classic = pd.read_csv( \
        "./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/test_set/" + response_variable + "_prs_calculation_classical.profile", \
        sep="\t", \
        header=0 \
    )
    if ( \
        (not fam_file_training[1].rename("IID").equals(combined_linear["IID"])) | \
        (not fam_file_test[1].rename("ID2").equals(profile_prs_classic["ID2"])) \
    ):
        raise ValueError("ERROR: FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLES FOR ELASTIC NET")


    print_text("compress the results", header=2)
    print_text("bed and bim plink files", header=3)
    run_bash(" \
        cd ./data/train_test_sets/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/; \
        gzip \
            --force \
            ./training_set/" + response_variable + "_training_set_plink_fileset.bed \
            ./training_set/" + response_variable + "_training_set_plink_fileset.bim \
            ./test_set/" + response_variable + "_test_set_plink_fileset.bed \
            ./test_set/" + response_variable + "_test_set_plink_fileset.bim; \
    ")
    
    print_text("bed and bim plink files", header=3)
    run_bash(" \
        cd ./results/final_results/" + covariate_dataset + "/" + response_variable + "/train_test_iter_" + str(iter_number) + "/; \
        gzip \
            --force \
            ./training_set/" + response_variable + "_training_elastic.effects \
            ./training_set/" + response_variable + "_training_elastic.probs \
            ./training_set/" + response_variable + "_training_linear.assoc \
            ./training_set/" + response_variable + "_training_linear.pvalues \
            ./training_set/" + response_variable + "_training_linear.score \
            ./training_set/" + response_variable + "_training_linear.summaries; \
    ")

# endregion




#########################################
# region run function across iterations #
#########################################

#iterate from 1 to total_iter (add 1 because the last element in range is not considered)
#iter_number=1; response_variable="distance_change"; covariate_dataset="small_set_predictors"
for iter_number in range(1, total_iter+1):
    prs_calc(iter_number, response_variable, covariate_dataset)

# endregion





####################
# region questions #
####################

print_text("Questions", header=1)
print(""" \

1) I have applied a quantile transformation to normalize phenotypes and covariates, having all of them in the same range. Is this recommended? I have run just one example with unstandarized data and the correlation between the PRS and the phenotype seems to be similar between both approaches.
      
2) Related to that point, I have noticed that the phenotype values in the column "Adjusted_Phenotype" from the .prs file obtained from --elastic do not match the values of the input phenotype for each sample. I guess an additional transformation/processing is applied to the phenotypes by LDAK?

3) When using --elastic, I guess I should use the default and assign all predictors a weighting of one. Then, LDAK will consider MAF for the estimation of the impact of each SNP, right?

4) Should covariates also be used with --calc-scores in the test set? Given that the phenotype was regressed against the covariates in the training set prior to estimating effect sizes, I guess it is not necessary to include them in the test set, right? Still I made an attempt and added the continuous covariates (plus coefficients with --coeffsfile) to --calc-score and the correlation improved from 0.04 to 0.06... Although as you will see in a question below, this difference falls within the normal variation I have seen analyzing the same phenotype several times with the same parameters....

5) As you suggested, I am running the classical PRS approach (ussing --linear) as baseline. I just used as input the plink files (plus phenotype data and covariates) of the training set, all default. Should I add a clumping step with --thin-tops to get a reduced set of SNPs that then can be analyzed again with --linear? I took a quick look with one of the phenotypes and given that the p-values for the SNPs are very high, the resulting list after clumping had only 20-30 SNPs....

6) In --linear, --permute YES would give the p-value in the case o no association because the phenotype is suffled, right?

7) Regarding the power argument for --calc-scores, not sure if --power should be 0 or -1. The effect sizes in the .effects files coming out from --elastic seem to be standarized (e.g., quantile 0.25=-3.483200e-07 and 0.75=3.477500e-07 in one case or 0.25=-2.958000e-08 and 0.75=3.049600e-08 in another), so I am using --power -1. But could you confirm me this is the case when using --elastic --LOCO NO and rest default? and the same applies for --linear? I checked one case with --linear in the "ALL" column and the 25th and 75th percentile were -4.072700e-02 and 4.199300e-02, respectively, so I guess it is the same.

8) Also about --calc-scores, I guess it is recommended to use the default options for --calc-scores and NOT use "--hwe-stand NO", right?

9) I have noticed that the results slightly changed after changing the ldak executable. I was using one I downloaded at the end of 2024 and I am now using one I just downloaded today, both ldak6.1.linux from github. The correlation between PRS and phenotype changed from 0.0008 to 0.002. Then I repeated the analysis with the same executable and the correlation changed to 0.001704. I have ensured the same set of samples is used in all cases. This degree of variation is expected?

""")

# endregion






print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/association_studies/
#chmod +x ./scripts/production_scripts/03b_prs_calculation.py
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03b_prs_calculation.py --total_iter=100 --response_variable="distance_change" --covariate_dataset="small_set_predictors" > ./03b_prs_calculation_small_distance_change_100_iter.out 2>&1
    #We should run the 8 dataset_size and phenotype combinations separtaely in the cluster, each combination run in serial 100 iterations, alwaysing requiering just 20GB, each iteration is around 10 minutes, so in 16-20 hours should be done.
#grep -Ei 'error|false|fail' ./03b_prs_calculation.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.
