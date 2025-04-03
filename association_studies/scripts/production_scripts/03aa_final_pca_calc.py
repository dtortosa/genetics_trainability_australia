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
######## CALCULATE FINAL PCAS ########
######################################


#######################
# region INITIAL INFO #
#######################


#we are going to re-calculate PCAs
    #TLDR: 
        #we are going to do it with the final dataset after QC just before imputation as recommended by Dr. Speed. we have also going to include MORE thinning to avoid very high loadings in one specific region of chromosome 8 for PCA 5, 4 and 3.
    #Summary sent to co-authors
        #I asked for advice to Dr. Speed regarding the selection of covariates in the context of his PRS tool. He basically said the same as Jonatan: we could apply a sensible criteria (e.g., select covariates significantly associated) and then check what happens doing the opposite, i.e., increasing the number of covariates to check our results are robust. In this context, we ended up discussing the weights (influence) of specific SNPs on the different PCAs as I mentioned that, for some PCAs, there was an abnormal accumulation of SNPs with very high weights in a narrow region of chromosome 8. What I originally did in the QC was just to discard these PCAs as covariates for future modeling, but Dr. Speed recommended recalculating the PCAs removing SNPs in high linkage-disequilibrium regions as these cause this kind of problem. This opened a can of worms as it could imply I should have to repeat all the QC steps, including imputation!
        #Dr. Speed recommended me to not do that as long as the PCAs have not been used for anything else besides outlier removal. Even if the number of outliers slightly differ when using the new PCAs, that should be ok as imputation results should be robust to that. I have checked what happens if, in the QC, we would have removed samples using the "new PCAs" and the difference is that 8 additional samples would remain. So our current approach does not add problematic samples, just remove maybe a few more. Note that "saving" these samples would require repeating everything with the whole dataset.
        #In summary, I do not think it makes sense to repeat all the QC steps and imputation due to this matter. As suggested by Dr. Speed, we can just add an additional step before the PRS calculation where we recalculate the PCAs with the clean set of SNPs and samples after QC, obtaining more robust and clean PCAs for modeling. The initial approach I followed had as consequence the removal of a few additional samples, which should not be a problem given our sample size (1008 final samples). 
    #Conversation with Dr. Speed
        #Is there any specific recommendation about the use of covariates for LDAK besides the need of splitting factors vs continuous covariates? I have the usual covariates (age, sex, body mass, phenotype baseline, PCAs...) and I am using only those covariates that are individually associated with the phenotype of interest to limit model complexity.
            #I examined this question about 15 years ago, to try and find a more justified solution (instead of simply "include 5 PCAs, age and sex"). I remember I experimented with only including PCs that were significantly associated with the phenotype, or only PCs with significant eigenvalues (a software called EigenSoft proposed a test for significance) . Howveer, my conclusion was that it is generally not too important (albeit I had larger sample sizes).
            #Then about 10 years ago, I switched to using heritability analysis to decide number of PCs - see the left and right analyses on this page https://dougspeed.com/quality-control/ - basically, include enough PCs that the sum of estimates from left and right halves of genome equal estimate from whole genome.
            #Nowdays, I normally include a lot (eg 10 or 20 PCs, maybe 10 internal PCs and 10 from 1000 GP), and if I wanted to check, I would repeat the analysis with fewer. Now again, this is larger sample sizes - the impact is larger for smaller sample sizes. But perhaps it remains a good strategy - I suggest you first perform the analysis using the best choice you can make (whether that is using 5 or 10, because another paper did similar, or choosing the significant ones only - provided you have a sensible criterium). Then repeat using a different choice as a sensitivity analysis (so if you included a lot first time, repeat with fewer). Provided the results are not very different (e.g., you have a SNP that is highly significant in one, but not in the other), then your choices should be good
        #Thank you for your detailed answer! Yes, that makes a lot of sense. I was indeed considering the p-values from eigensoft, specifically those axes with a p-value<0.01 (first 5 in my case). Then I checked the distribution of SNPs loadings across the genome and found a suspicious increase of positive and negative loadings in a narrow region of chromosome 8 that was very marked for PC5, but was also visible for PC4 and PC3. So I decided to use only PC1 and PC2 for subsequent analyses. Then, as I noted before, I included these two PCAs if they were associated with the phenotype. Maybe this is a very stringent approach, but we are using clear criteria and it could make sense for a small dataset. I will follow your advice and check what happens including more axes (e.g., all significant axes according to eigensoft, 10, etc...)
            #Yes, very good to check the loadings (so be suspicous if you find a very dense patch - similar to what you observed on chr 8). Similarly, check the PCs and be suspicoius if you find groupings (e.g, divide into 2 or 3) that do not seem to correlate with ancestry. Both of these normally indicate problems due to LD - that the loadings are maiinly concentrated on a region of highLD, so instead of picking up genome-wide variation, they just record the genotypes for that small region
            #If this occurs, it is best to repeat with more thinning - and if human, I recommend excluding long range ld regions https://dougspeed.com/high-ld-regions/
        #Amazing all these insights! Yes, I initially checked the existence of clusters and thankfully no clear cluster was present, but I completely missed the option to remove the SNPs with the high loadings instead of the PCAs showing the problem. After increasing thinning (plink2 --indep-pairwise from 500kb 0.2 to 600kb 0.1) the peak completely disappears and the rest of regions/PCAs look perfect in this regard. At this point, after performing all QC steps and imputation, do you think it would necessary to re-run all my pipeline from this step? or it would be reasonable to just recalculate the PCAs with the increased thinning (to be used as covariates) and remove the high-LD from the plink fileset just before running LDAK?
            #Thanks, that is good to hear
            #In terms of your pipeline, my quick answer would be that you only need to repeat the steps that used the PCs - from looking above, you did the following: 1.2K samples and 3M of genotyped+imputed SNPs after QC - then make PRS using --elastic
            #So here I dont think the imputation depends on the PCs, so it sounds fine to leave, especially as I imagine that is the slowest step (well, maybe you used PCs to decide which individuals to exclude - but imputation results should be robust to changing the individuals slightly).
            #Then for --elastic, you can repeat with the new PCs. Note that it is fine to include high ld regions for --elastic (ie, remove them when computing the PCS, to avoid the uneven loadings, but you can retain them when making the PRS)
        #Exactly, I am only using the PCAs as covariates and to exclude outliers (iterative removal with smartpca). I have checked what happens doing the removal with higher thinning and the difference is that 8 samples are retained, while with less thinning they are removed. I understand this should not have a great impact on imputation considering the total size of my cohort. So yes, I think I will calculate the PCAs avoiding high ld regions and then use them as covariates for --elastic (retaining in that case high-ld regions). Would you recommend to calculate the PCAs with the fileset generated after imputation and QC? or with the fileset I originally used for the PCAs during pre-imputation QC?
            #I think it should not matter - normally, I compute PCs using only high quality (eg directly genotyped ) snps, rather than also include imputed
            #But I think general advice is to use only the final list of samples (and not the original list of samples). A simple case is if you had a few ethnic outliers, that were excluded during qc, then these would dominate PCs computed using all samples, so it would probably be best to use only the samples after qc
        #Ok, so I can calculate the PCAs using the data generated at the end of QC (when all sample and SNP filters have been applied) but before the imputation. Thank you so much for your insights!

# endregion






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
    mkdir -p ./data/final_pcas; \
    ls -l \
")

# endregion






#########################################
# region SELECT ONLY GENOTYPES VARIANTS #
#########################################
print_text("SELECT ONLY GENOTYPES VARIANTS", header=1)
print_text("create folders", header=2)
run_bash(" \
    mkdir \
        -p \
        ./data/final_pcas/00_typed_variants_after_qc/list_variants/ \
")


print_text("define function to extract the genotypes SNPs across chromosomes", header=2)
#chromosome=21
def genotyped_snps_prep(chromosome):

    #select only TYPED variants and update their IDs using position, then create plink files
    run_bash(" \
        cd ./data/final_pcas/00_typed_variants_after_qc/list_variants/; \
        chr_numb=" + str(chromosome) + "; \
        bcftools filter \
            --include 'INFO/R2>=0.95 && INFO/TYPED=1' \
            ../../../../../quality_control/data/genetic_data/quality_control/20_imputation_results/04_uncompressed_vcf_files/chr${chr_numb}.dose.vcf.gz | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2 \
            --min-alleles 2 | \
        bcftools view \
            --types snps | \
        bcftools annotate \
            --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' | \
        bcftools query \
            --format '%CHROM\_%POS\_%REF\_%ALT %INFO/TYPED %INFO/R2\n' | \
        awk \
            'BEGIN{FS=\" \"; OFS=\"\t\"}{print $1, $2, $3}' > ./chr${chr_numb}_snp_list.tsv \
    ")
    #for each chromosome
        #create a specific folder and move inside
        #include SNPs with TYPE flag
            #INFO/TYPED=1 filters for variants with the TYPED flag
            #INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped">
                #check 02bd_download_extract_imputation_results.sh in the script folder of QC
            #Therefore, we are discarding SNPs that were NOT genotyped
        #Remove SNPs with imputation quality < 0.95
            #This is more stringent than Ritchie´s but it is the recommendation of Doug Speed for using its PRS tool (0.95 or 0.99). Also remember Augusto used a filter of 0.9, not 0.7.
                #AUGUSTO POST-IMPTUATION: Once our genomic data were imputed, a second quality control analysis was performed with PLINK 1.9 software [13]. The second quality control exclusion criteria were: low imputation quality (𝑅2<0.9); variants that did not meet the Hardy–Weinberg equilibrium (HWE-P>10−6); and low minor allele frequency (MAF<0.01) [14].
            #Indeed, the removal of more SNPs using imputation quality has decreased the number of SNPs removed in the post-imputation quality control, suggesting we have targeting problematic SNPs when using the imputation quality.
            #In "02bd_download_extract_imputation_results.sh" you can see how R2 is the imputation quality.
            #Doug Speed: For example, when analysing the UK Biobank data, our information score threshold was 0.95. However, if in your data, almost all (very few) predictors exceed this threshold, then you should consider increasing (reducing) it. For some good advice on how to decide thresholds, see Mike Weale's Chapter on Quality Control.
                #In our case, we have reduced from 28M to 9M SNPs by increasing the imputation score from 0.7 to 0.95. This leaves 3M after QC.
                #If we increase the threshold to 0.99, the number of SNPs after QC reduces to 200K. This is basically the number of SNPs we have before imputation, so this is too much. 
                #We select a stringent threshold, but not the most stringent.
        #join biallelic SNPs into multiallelic records: We have SNPs with the same position but one allele different, suggesting we have multiallelic variants. We want to merge them and then remove.
            #--multiallelic +snps: 
                #this combines snps with the same position and at least equal REF or ALT. Therefore, this makes that a SNP with three alleles is within just one row, so we can filter it in the next step.
                #join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
        #remove variants with more than 2 alleles, as these are multiallelic (multiallelic SNPs were merged in the previous step)
            #applying this step and removing multiallelic variants makes all the merging errors to dissappear
                #PLINK cannot properly resolve genuine triallelic variants. We recommend exporting that subset of the data to VCF, using another tool/script to perform the merge in the way you want, and then importing the result. Note that, by default, when more than one alternate allele is present, --vcf keeps the reference allele and the most common alternate. (--[b]merge's inability to support that behavior is by design: the most common alternate allele after the first merge step may not remain so after later steps, so the outcome of multiple merges would depend on the order of execution.)
                #https://www.cog-genomics.org/plink/1.9/data#merge
        #select only SNPs 
            #some variants are indels, for example, chr21_10090492_AG_A (rs1211611360) or chr21_10088665_GC_G (rs1241342440) in chromosome 21.
            #Also, I have detected that these indels have a position that differ in 1 based respect to NCBI matching the "Canonical SPDI:". I do not fully understand this, but the point is that the non-indel SNPs have a coordinate matching exactly NCBI hg38. Given we are removing indels, I think we are good.
        #update the IDs of ALL variants using position and alleles: 
            #We have many SNPs without ID (".")
            #Also, in general, it is a good idea to remove rs IDs: "the rs ID system of naming variants is problematic because it's not curated sufficiently. A single rs ID can relate to multiple variants at the same position. Also, a single rs ID can relate to variants at more than 1 position in the genome. The fundamental issue being that rs IDs are not unique at all... It's too risky working with rs IDs on clinical data, where mix-ups / mess-ups just aren't allowed..."
                #https://www.biostars.org/p/281276/#281670
            #--set-id: assign ID on the fly. By default all existing IDs are replaced. If the format string is preceded by "+", only missing IDs will be set. For example, one can use: bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' file.vcf
            #Note that %FIRST_ALT gets the first alternate allele at the variant position, but we have already remove multiallelic SNPs, so this is not a problem, all SNPs left have only 1 ALT allele.
        #query the output to get the ID (based on chromosome, position and alleles), TYPED and R2 columns
            #note that we assesed and solved strand problems and ensure we had the same SNPs than in the imputation panel and even switch REF-ALT based on the imputation preparion tool, so the REF/ALT alleles should be ok.
        #then ensure that the columns are tab delimitted and save
        #we have applied the same steps for the last VCF files generated in this project, so we are sure we are using the last, cleaned data

    #check we correctly applied the TYPED flag and the imputation quality filter
    run_bash("\
        cd ./data/final_pcas/00_typed_variants_after_qc/list_variants/; \
        chr_numb=" + str(chromosome) + "; \
        n_snps_under_95_non_typed_imput=$( \
            awk \
                'BEGIN{FS=\"\t\"}{ \
                    if($2!=1 || $3<0.95){ \
                        count++ \
                    } \
                }END{print count}' \
                ./chr${chr_numb}_snp_list.tsv \
        ); \
        if [[ $n_snps_under_95_non_typed_imput -ne 0 ]]; then \
            echo \"ERROR! FALSE! WE HAVE SNPS WITH IMPUTATION QUALITY UNDER 0.95 \"; \
        fi; \
    ")
        #If the second field (TYPED) is not equal to 1 or the third field (R2) is less than 0.95, we count it
        #If the count is not equal to 0, we have SNPs with imputation quality under 0.95 or SNPs that were not genotyped


print_text("run function in parallel across chromosomes", header=2)
print_text("set chromosomes", header=3)
chromosomes = [str(i) for i in range(1, 23, 1)]
if (len(chromosomes) != 22):
    raise ValueError("ERROR! FALSE! NOT ALL CHROMSOMES ARE RUNNING")
else:
    print(chromosomes)

print_text("parallelize", header=3)
import multiprocessing as mp
pool = mp.Pool(len(chromosomes))
pool.map(genotyped_snps_prep, chromosomes)
pool.close()


print_text("merge all chromosomes", header=2)
import pandas as pd
import os

print_text("get a sorted list of all .tsv files in the folder with natural sorting (i.e., 1,2,3...", header=3)
from natsort import natsorted
import os
tsv_files = natsorted([f for f in os.listdir("./data/final_pcas/00_typed_variants_after_qc/list_variants/") if f.endswith("_snp_list.tsv")])
if(len(tsv_files)!=22):
    raise ValueError("ERROR! FALSE! NOT ALL CHROMSOMES ARE RUNNING")
else:
    print(tsv_files)

print_text("sequentially merge the .tsv files", header=3)
merged_df = pd.DataFrame()
#tsv_file=tsv_files[0]
for tsv_file in tsv_files:
    file_path = os.path.join("./data/final_pcas/00_typed_variants_after_qc/list_variants/", tsv_file)
    temp_df = pd.read_csv(file_path, sep="\t", header=None)
    merged_df = pd.concat([merged_df, temp_df], ignore_index=True)
merged_df.columns = ["ID", "TYPED", "R2"]

print_text("save the full DF and also just the IDs column", header=3)
merged_df.to_csv("./data/final_pcas/00_typed_variants_after_qc/merged_snps_list.tsv", sep="\t", index=False, header=True)
merged_df["ID"].to_csv("./data/final_pcas/00_typed_variants_after_qc/merged_snps_list_id_only.txt", index=False, header=False)

print_text("check we have all the SNPs", header=3)
#sum the number of rows across all the chromosomes
total_n_snps = run_bash(" \
    awk \
        'END{print NR}' \
        ./data/final_pcas/00_typed_variants_after_qc/list_variants/*_snp_list.tsv \
", return_value=True).strip()
#this should be equal to the number of rows in the merged file
if(total_n_snps != str(merged_df.shape[0])):
    raise ValueError("ERROR! FALSE! NOT ALL SNPS ARE MERGED")


print_text("prepare plink fileset", header=2)
print_text("check we have the correct versions of plink and plink2. We are using plink1.9 for mergning", header=3)
plink_version=run_bash("plink --version", return_value=True).strip()
plink2_version=run_bash("plink2 --version", return_value=True).strip()
if((plink_version!="PLINK v1.90b7 64-bit (16 Jan 2023)") | (plink2_version!="PLINK v2.00a4.2LM 64-bit Intel (31 May 2023)")):
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH PLINK VERSIONS")
else:
    print("PLINK VERSIONS ARE OK")

print_text("subset the plink filesets", header=3)
run_bash(" \
    plink \
        --bfile ../quality_control/data/genetic_data/quality_control/21_post_imputation_qc/03_third_qc_step/merged_3_geno \
        --extract ./data/final_pcas/00_typed_variants_after_qc/merged_snps_list_id_only.txt \
        --make-bed \
        --out ./data/final_pcas/00_typed_variants_after_qc/merged_3_geno_only_typed \
")
    #We use the last fileset just before the tests of batch effects (see line 322 in 02d_post_imputation_qc.py)
    #--extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.

print_text("check we have the correct number of SNPs", header=3)
run_bash(" \
    n_snps_pca=$( \
        awk \
            'BEGIN{FS=\"\t\"}{ \
        }END{print NR}' \
        ./data/final_pcas/00_typed_variants_after_qc/merged_3_geno_only_typed.bim \
    ); \
    if [[ $n_snps_pca < 242000 || $n_snps_pca > 243000 ]]; then \
        echo \"ERROR! FALSE! WE HAVE A PROBLEM WITH THE FINAL NUMBER OF SNPS FOR PCAS\"; \
    fi \
")

print_text("check we have the correct number of samples", header=3)
run_bash(" \
    n_samples_plink_qc=$( \
        awk \
            'BEGIN{FS=\"\t\"}{ \
        }END{print NR}' \
        ../quality_control/data/genetic_data/quality_control/21_post_imputation_qc/03_third_qc_step/merged_3_geno.fam \
    ); \
    n_samples_pca=$( \
        awk \
            'BEGIN{FS=\"\t\"}{ \
        }END{print NR}' \
        ./data/final_pcas/00_typed_variants_after_qc/merged_3_geno_only_typed.fam \
    ); \
    if [[ $n_samples_pca -ne $n_samples_plink_qc ]]; then \
        echo \"ERROR! FALSE! WE HAVE A PROBLEM WITH THE FINAL NUMBER OF SAMPLES FOR PCAS\"; \
    fi \
")

print_text("check we do NOT have a .missnp file", header=3)
run_bash(" \
    cd ./data/final_pcas/00_typed_variants_after_qc/; \
    if [[ -e merged_3_geno_only_typed.missnp ]]; then \
        echo 'ERROR! FALSE! FILTERING OF SNPS DID NOT WORK'; \
    else \
        echo 'FILTERING WORKED'; \
    fi \
")

print_text("see the total number of SNPs", header=3)
run_bash(" \
    cd ./data/final_pcas/00_typed_variants_after_qc; \
    awk \
        'END{print NR}' \
        ./merged_3_geno_only_typed.bim \
")
    
# endregion






#####################################
# region DO LD CLEANING OF VARIANTS #
#####################################
print_text("DO LD CLEANING OF VARIANTS", header=1)
print_text("create folders", header=2)
run_bash(" \
    mkdir \
        -p \
        ./data/final_pcas/01_ld_cleaning/ \
")

print_text("Filter out snps with high LD", header=2)
print_text("get a list of SNPs with low LD", header=3)
run_bash(" \
    cd ./data/final_pcas/; \
    plink2 \
        --bfile ./00_typed_variants_after_qc/merged_3_geno_only_typed \
        --indep-pairwise 600kb 1 0.1 \
        --out ./01_ld_cleaning/ld_prunning_filter; \
    ls -l ./01_ld_cleaning/")
        #This command produces a pruned subset of variants that are in approximate linkage equilibrium with each other, writing the IDs to plink2.prune.in (and the IDs of all excluded variants to plink2.prune.out). These files are valid input for --extract/--exclude in a future PLINK run; and, for backward compatibility, they do not affect the set of variants in the current run.
        #Since the only output of these commands is a pair of variant-ID lists, they now error out when variant IDs are not unique.
        #--indep-pairwise is the simplest approach, which only considers correlations between unphased-hardcall allele counts. We cannot use the alternative, which is --indep-pairphase and requires phased data. --indep-pairwise takes three parameters: 
            #a required window size in variant count or kilobase (if the 'kb' modifier is present) units, 
            #an optional variant count to shift the window at the end of each step (default 1, and now required to be 1 when a kilobase window is used).
                #I guess the window is shifted (moved forward) until the variant count changes in 1 unit?
            #a required r2 threshold.
        #At each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, and variants are greedily pruned from the window until no such pairs remain.
        #For example:
            #"--indep-pairwise 100kb 1 0.8"
            #removes SNPs so that no pair within 100 kilobases have squared-allele-count-correlation (r2) greater than 0.8, and saves the IDs of the remaining SNPs.
        #By default, when given a choice, the variant-pruning commands preferentially keep variants with higher nonmajor allele frequencies. However, if you provide a list of variant IDs to --indep-preferred, all variants in that list are prioritized over all variants outside it. (Allele frequencies will still be used for tiebreaking.)
            #it seems reasonable to select the SNP with the higher MAF within pairs of correlated SNPs.
        #On human data, some reasonable parameter settings are, in order of increasing strictness:
            #"--indep-pairwise 100kb 1 0.8"
            #"--indep-pairwise 200kb 1 0.5"
            #"--indep-pairwise 500kb 1 0.2"
        #As we get more strict, the number of selected SNPs is reduced. We will go for the most strict option, i.e., allowing a lower correlation between SNPs in larger windows. This means that we will have much less variants correlated across larger chunks of the genome. In this way, we will have a very clean set in terms of LD. This is important because we are going to do operations that are not LD-aware.
            #https://www.cog-genomics.org/plink/2.0/ld
            #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #IMPORTANT: we have gone even a bit further the most strict option. I did this because with "500kb 1 0.2" I was still getting a peak of loadings in chr8. After increasing the window to 600kb and the r2 to 0.1, this peak completely dissapeared and the rest of regions/PCAs look perfect in this regard. 

print_text("select only SNPs with low LD", header=3)
run_bash(" \
    cd ./data/final_pcas/; \
    plink \
        --bfile ./00_typed_variants_after_qc/merged_3_geno_only_typed \
        --extract ./01_ld_cleaning/ld_prunning_filter.prune.in \
        --make-bed \
        --out ./01_ld_cleaning/merged_3_geno_only_typed_ld_prunning;\
    ls -l ./01_ld_cleaning/ \
")
    #from the current fileset, select only those SNPs included in .prune.in
    #we use extract for that
        #--extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.
            #https://www.cog-genomics.org/plink/1.9/filter
    #make a new fileset

print_text("checks after filtering", header=3)
#In the Ritchie's tutorial, they ended up with 67,000 autosomal variants in linkage equilibrium in order to calculate IBD and pi_hat. I guess they used the same dataset for the PCA. They also say that "It is recommended that the user prune to approximately ∼100,000 SNPs for use in PCA analyses".
    #14_pop_strat/01_pca
#In the admixture tutorial and Admixture docs they use a more stringent LD filtering, but they say that if you want to control for population structure between related populations (i.e., within the same continent) then 10K is not enough and you should use 100K. 
        #--indep-pairwise 50 10 0.1
        #http://dalexander.github.io/admixture/admixture-manual.pdf
#It seems to be a good idea to have between 50K and 100K SNPs for PCA. For admixture you need more, but here the most important thing are the PCAs to control for population structure and in that case we are not very far off, around 50K while Ritchie used 67K.
run_bash(" \
    cd ./data/final_pcas/01_ld_cleaning/; \
    ld_snps_in=$( \
        awk \
            'BEGIN{FS=\" \"}END{print NR}' \
            ./merged_3_geno_only_typed_ld_prunning.bim); \
    printf 'The number of included autosomal SNPs in linkage equilibrium is %s\n' \"$ld_snps_in\"; \
    echo 'Is this number greater than 52K?'; \
    if [[ $ld_snps_in -gt 52000 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi \
")
    #you can print a variable with text using printf and %s for "string". Then you call the variable you want to add within "". You could use two variables: "$var1 $var2"
        #https://phoenixnap.com/kb/bash-printf

print_text("We are not removing SNPs by high-LD regions", header=3)
print("""
The high-LD regions indicated by Dr.Speed are in hg19 coordinates, while we are working with hg38. I tried liftovering the regions, but I was not able to do it and got extrange results like getting hg38 coordinates in the X chromosome when we do not have sex chromosomes.

Given that just incresing the prunning to 600kb and 0.1 made the SNP weigths peak completely dissapear. Importantly, if you compare the PCA plots between this subset and the PCA during QC, you can see they are really similar, it seems we are still detecting the same variability. So I think we are good.
      
This is what we have anyways. If we want to have SNPs with very low-LD, no peaks in the PCA loadings and, in general, the clean dataset, this is what is left after QC.
""")

# endregion






#######################
# region RUN SMARTPCA #
#######################
print_text("RUN SMART PCA", header=1)
print_text("create folders", header=3)
run_bash(" \
    mkdir \
        -p \
        ./data/final_pcas/02_smartpca/ \
")

print_text("PCA with smartpca", header=3)
#smartpca: general explanations
    #smartpca runs Principal Components Analysis on input genotype data and  outputs principal components (eigenvectors) and eigenvalues. We note that eigenvalue_k/(Sum of eigenvalues) is the proportion of variance explained by eigenvector_k. The method assumes that samples are unrelated. However, a small number of cryptically related individuals is usually  not a problem in practice as they will typically be discarded as outliers (we already removed related individuals).
    #The syntax of smartpca is "../bin/smartpca -p parfile"
        #The below for details about each argument

print_text("prepare fam file for smartpca", header=3)
#As a minor issue, smartpca ignores individuals in the .fam file if they are marked as missing in the phenotypes column. This awk command provides a new .fam file that will automatically include all individuals.
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
#copy also bim and bed files
run_bash(" \
    cd ./data/final_pcas/; \
    awk \
        '{print $1,$2,$3,$4,$5,1}' \
        ./01_ld_cleaning/merged_3_geno_only_typed_ld_prunning.fam > \
    ./02_smartpca/merged_3_geno_only_typed_ld_prunning.new_fam.fam; \
    cp ./01_ld_cleaning/merged_3_geno_only_typed_ld_prunning.bim ./02_smartpca; \
    cp ./01_ld_cleaning/merged_3_geno_only_typed_ld_prunning.bed ./02_smartpca \
")

print_text("decide the number of axes to output", header=3)
n_axes=20

print_text("smartpca parameter file", header=3)
run_bash(" \
    cd ./data/final_pcas/; \
    PREFIX=merged_3_geno_only_typed_ld_prunning; \
    echo genotypename: ./02_smartpca/$PREFIX.bed > ./02_smartpca/$PREFIX.par; \
    echo snpname: ./02_smartpca/$PREFIX.bim >> ./02_smartpca/$PREFIX.par; \
    echo indivname: ./02_smartpca/$PREFIX.new_fam.fam >> ./02_smartpca/$PREFIX.par; \
    echo snpweightoutname: ./02_smartpca/$PREFIX.snpeigs >> ./02_smartpca/$PREFIX.par; \
    echo evecoutname: ./02_smartpca/$PREFIX.eigs >> ./02_smartpca/$PREFIX.par; \
    echo evaloutname: ./02_smartpca/$PREFIX.eval >> ./02_smartpca/$PREFIX.par; \
    echo phylipoutname: ./02_smartpca/$PREFIX.fst >> ./02_smartpca/$PREFIX.par; \
    echo numoutevec: " + str(n_axes) + " >> ./02_smartpca/$PREFIX.par; \
    echo outliersigmathresh: 6 >> ./02_smartpca/$PREFIX.par; \
    echo numoutlieriter: 11 >> ./02_smartpca/$PREFIX.par; \
    echo numoutlierevec: 10 >> ./02_smartpca/$PREFIX.par; \
    echo outlieroutname: ./02_smartpca/$PREFIX.outliers >> ./02_smartpca/$PREFIX.par; \
    echo altnormstyle: YES >> ./02_smartpca/$PREFIX.par; \
    echo missingmode: NO >> ./02_smartpca/$PREFIX.par; \
    echo ldregress: 0 >> ./02_smartpca/$PREFIX.par; \
    echo noxdata: YES >> ./02_smartpca/$PREFIX.par; \
    echo nomalexhet: YES >> ./02_smartpca/$PREFIX.par; \
    echo newshrink: NO >> ./02_smartpca/$PREFIX.par \
")
    #WE ARE USING THE SAME PARAMETERS THAN IN THE QC
    #script for smartpca parameter file from ADMIXTURE TUTORIAL
        #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #genotypename:
        #contains genotype data for each individual at each SNP
    #snpname:
        #contains information about each SNP 
    #indivname:
        #contains information about each individual
        #accorind to the readme, the genotype and snp file can be just .bed and.bim files, respectively, from plink format, which is our format. 
        #according to the admixture tutorial, you can use the plink format, i.e., bed, bim and fam file as input. In the case of the fam file, the phenotypes have to be 1 for all samples, if not, they are not considered (see above)
            #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4
    #snpweightoutname:
        #output file containing SNP weightings of each principal component. Note that this output file does not contain entries for monomorphic SNPs from the input .snp file. 
    #evecoutname:
        #output file of eigenvectors. See numoutevec parameter below.
    #evaloutname:
        #output file of all eigenvalues
    #phylipoutname:
        #output file containing an fst matrix which can be used as input to programs in the PHYLIP package, such as the "fitch" program for constructing phylogenetic trees.
        #we are generating this just in case.
    #numoutevec:
        #number of eigenvectors to output.  Default is 10.
        #this should affect the results, for example, the axes used for outlier removal are indicated in another argument.
    #outlier removal. According to plink tutorial (https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3)
        #It is also a good idea to throw out gross outliers at this point; any sample which is more than, say, 8 standard deviations out on any top principal component is likely to have been genotyped/sequenced improperly; you can remove such samples by creating a text file with the bad sample IDs, and then using –remove +  –make-bed:
            #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3#Sec22
        #The phrase "8 standard deviations out" refers to a data point that falls 8 standard deviations away from the mean of a dataset. In statistics, the standard deviation is a measure of the amount of variation or dispersion in a set of values. A low standard deviation means that the values tend to be close to the mean, while a high standard deviation means that the values are spread out over a wider range. When a data point is said to be "8 standard deviations out", it means it's quite far from the mean. In a normal distribution, almost all data falls within 3 standard deviations of the mean. So, a data point that is 8 standard deviations away from the mean is extremely rare and could be considered an outlier. In the context of your sentence, it suggests that any sample which falls more than 8 standard deviations away on any top principal component might have been genotyped or sequenced improperly, and thus, it might be a good idea to remove these outliers from the analysis. This is because such extreme values can significantly skew the results and interpretations of the data analysis.
        #If you do this, follow it up by repeating the PCA, since the bad samples might have distorted the principal components. Occasionally, the new principal components will reveal another bad sample, and you have to repeat these two steps, etc. We are using the iterative process of smartpca as recommended in the tutorial
        #Steps
            #The program identifies individuals whose principal component scores deviate significantly from the majority of the population. This is typically done by setting a threshold for the number of standard deviations from the mean.
            #Iterative Process: Outlier removal is often an iterative process. The program may run PCA multiple times, each time removing individuals who are identified as outliers, until no more outliers are detected.
            #Impact on Analysis: Removing outliers helps to ensure that the PCA results are not skewed by individuals who have unusual genetic backgrounds. This can provide a clearer picture of the population structure and improve the accuracy of downstream analyses.
        #outliersigmathresh: 
            #We are using the default (6).
            #number of standard deviations which an individual must exceed, along one of the top (numoutlierevec) principal components, in order for that individual to be removed as an outlier.  Default is 6.0.
        #numoutlieriter:
            #maximum number of outlier removal iterations. Default is 5.  To turn off outlier removal, set this parameter to 0.
            #I have checked that outliers are removed until iteration 7, so we just using until 8.
        #numoutlierevec:
            #number of principal components along which to remove outliers during each outlier removal iteration.  Default is 10.
            #We are using the default. The PCA axes with a significant p-value (<0.05) are the first 11, but I have checked that using 11 instead of 10 for outlier removal gives the same outliers. In addition, I have checked the eigenvectors that detect the outliers and never 10 or 11 are included, so we do not need to consider more axes. We are sticking to the default.
        #outlieroutname:
            #output logfile of outlier individuals removed. If not specified, smartpca will print this information to stdout, which is the default.
    #altnormstyle:
        #Affects very subtle details in normalization formula. Default is YES (normalization formulas of Patterson et al. 2006). To match EIGENSTRAT (normalization formulas of Price et al. 2006), set to NO.
        #I have checked with and without altnormstyle and there are almost no differences.
    #missingmode:
        #If set to YES, then instead of doing PCA on # reference alleles, do PCA on whether each data point is missing or nonmissing.  Default is NO.
    #ldregress:
        #If set to a positive integer, then LD regression is turned on, and input to PCA will be the residual of a regression involving that many previous SNPs, according to physical location.  See Patterson et al. 2006. Default is 0 (no LD regression).  If desiring LD correction, we recommend 200.
        #this is done to correct for the correlation (i.e., the linkage disequilibrium) of the SNPs. We do NOT to do this because we have already prunned our data considering LD.
    #noxdata: 
        #if set to YES, all SNPs on X chr are excluded from the data set. The smartpca default for this parameter is YES, since different variances for males vs. females on X chr may confound PCA analysis.
        #Using YES, but not required anyway because we have only autosomal SNPs.
    #nomalexhet:
        #if set to YES, any het genotypes on X chr for males are changed to missing data. The smartpca default for this parameter is YES.
        #males should be homozigous for X (except pseudo-autosomic regions)
        #Using YES, but not required anyway because we have only autosomal SNPs.
    #lsqproject
        #PCA projection is carried out by solving least squares equations rather than an orthogonal projection step. This is approriate if PCs are calculated using samples with little missing data but it is desired to project samples with much missing data onto the top PCs. In other words, most of the samples have a lot of data, but some samples have a lot of missing, in that situation, instead of filling gaps with the average, this approach does something different that solves the problem. BUT, if the sample has few missing, then this works as the default orthogonal. I have checked that setting this to NO or YES does not change the results.
        #./EIG-8.0.0/POPGEN/lsqproject.pdf
    #shrinkmode/newshrink
        #A problem with smartpca is that samples used to calculate the PC axes "stretch" the axes. So that 2 populations in fact genetically identical (2 independent samples from the same underlying population) will appear different if one is used to compute axes, and one not.  shrinkmode: YES is an attempt to solve this problem.  Details to appear later, but this has been used successfully in the Reich lab.*** warning *** shrinkmode is slow and will greatly increase the runtime. (NEW) New version:  newshrink:  YES technical variation of shrinkmode,  should be (slightly). 
        #According to chatGTP, even within a single population, if there is significant internal structure, `shrinkmode` can help ensure that the PCA axes are not unduly influenced by subgroups within your population.
        #HOWEVER, I have run with and without shrink mode and the results are EXACTLY the same, with all decimals, while the run time goes up to 40 min from 5-10. We are not using this mode.
    #Multithreading (Code added by Chris Chang).
        #smartpca now supports multithreading but NOT with fastmode: YES. By default a (hopefully) system dependent number of threads is chosen. This can be overwritten by (for example) numthreads:   10
        #it seems by default it is using all the cores.

print_text("run smartpca", header=3)
#we use the version 8.0.0, which is the latest at the moment of writting (sep 2024). This was released in october 2022.
    #see the contianer receipte for further details about how install in the container
    #https://github.com/DReichLab/EIG/releases
run_bash(" \
    cd ./data/final_pcas/; \
    /bin/smartpca -p ./02_smartpca/merged_3_geno_only_typed_ld_prunning.par > ./02_smartpca/smart_pca_run.out \
")
    #smartpca also calculate the "Tracy-Widom statistics"
        #I have visually compared this table and the generated by twstats. They are the same, just a few decimals are different. Also that table has a column with eigenvalues and a p-value. Therefore, smartpca is also using twstats
            #see below for the code to ran twstats
        #The twstats program computes Tracy-Widom statistics to evaluate the statistical significance of each principal component identified by pca (Patterson et al. 2006).
        #it helps determine whether the observed eigenvalues (which correspond to the variance explained by each PC) are significantly larger than what would be expected by chance.
        #steps (from copilot):
            #Eigenvalue Calculation: After performing PCA, the eigenvalues corresponding to each principal component are calculated. These eigenvalues represent the amount of variance explained by each PC.
            #Comparison with Tracy-Widom Distribution: The largest eigenvalues are compared to the Tracy-Widom distribution. This comparison helps determine if the observed eigenvalues are significantly larger than those expected under the null hypothesis (i.e., no structure in the data, the axes are not absorbing any structure from the data).
            #Significance Testing: The program computes p-values for each eigenvalue based on the Tracy-Widom distribution. A low p-value indicates that the corresponding principal component explains a significant amount of variance, suggesting it captures meaningful structure in the data.
        #The twstats program assumes a random set of markers, and should not be used on data sets of ancestry-informative markers, as admixture-LD may violate its underlying assumptions. 
            #Ancestry-informative markers are those SNPs that significantly different between pops and are speficically used for infer ancestry of individuals.
            #This is not our case, as we have not specifically selected SNPs that differ between pops. 
            #Anyways, we are going to compare the significant PCA axes according to twstats and admixture program, so we are good here.

#in case you want to run twstats separately
'''
if False:
    run_bash(" \
        cd ./data/genetic_data/quality_control/14_pop_strat/01_pca/; \
        /usr/bin/twstats \
            -t ../../../../../../eigensoft_versions/EIG-8.0.0/POPGEN/twtable \
            -i ./eigen_out/loop_maf_missing_2_ldprunned_autosome_pca.eval \
            -o ./eigen_out/loop_maf_missing_2_ldprunned_autosome_pca_tracy.out")
'''

print_text("check the number of samples after outlier removal", header=3)
run_bash(" \
    cd ./data/final_pcas/02_smartpca; \
    n_outliers=$( \
        awk \
            'END{print NR}' \
            ./merged_3_geno_only_typed_ld_prunning.outliers \
    ); \
    if [[ $n_outliers -ne 0 ]]; then \
        echo 'ERROR! FALSE! WE HAVE A PROBLEM, WE HAVE NEW OUTLIERS BUT SHOULD NOT'; \
    fi \
")

print_text("see significant axes", header=3)
#extract the tracy table
run_bash(" \
    cd ./data/final_pcas/02_smartpca/; \
    awk \
        'BEGIN{ \
            OFS=\"\t\"; \
            start_table=0; \
            end_table=0; \
        }{ \
            if($0 ~ /^## Tracy-Widom statistics: rows:/){ \
                start_table=NR \
            }; \
            if($0 ~ /kurtosis/){ \
                end_table=1 \
            }; \
            if(start_table !=0 && NR > start_table && end_table==0){ \
                print $1,$2,$3,$4,$5,$6 \
            }; \
        }' \
        ./smart_pca_run.out > tracy_table.tsv \
")
    #force the output to be tab, and create two variables as zero
    #if the has the header of tracy table
        #save the number of the row to know when the table start
    #if the row includes kurtosis
        #we are at the end of the table, so activate ending
    #if start_table is not zero, and it is greater than the first row of the tracy table and we are not yet at the end of the table 
        #print all fields separately

print_text("load the table", header=4)
tracy_table = pd.read_csv(\
    "./data/final_pcas/02_smartpca/tracy_table.tsv", \
    sep="\t", \
    header=0, \
    low_memory=False)
print(tracy_table)

print_text("print significant axes", header=4)
tracy_signi_axes = tracy_table.loc[tracy_table["p-value"] < 0.05,:]
print(tracy_signi_axes)
print(f"""
Results:
We have {tracy_signi_axes.shape[0]} significant axes, meaining that they explain more genetic variance than expected by chance.
""")

print_text("plot smartpca eigenvectors", header=3)
print_text("process the eigenvector file to ensure is tab delimited", header=4)
run_bash(" \
    cd ./data/final_pcas/02_smartpca; \
    awk \
        'BEGIN{ \
            OFS=\"\t\" \
        }{ \
            if(NR>1){ \
                print " + ",".join(f"${i}" for i in range(1,n_axes+2)) + " \
            } \
        }'\
        ./merged_3_geno_only_typed_ld_prunning.eigs > ./merged_3_geno_only_typed_ld_prunning_tab.eigs \
")
    #force the OFS to be tabs
    #select all except the first 1, i.e., NR>1. The first one has eigenvalues of rach PCA axis. For each row
    #print the following columns
        #all columns from the first to the last one
        #the last is based on the number of PCA axes we have plus 1 because we have an ID column and 1 more because the last element of the range is not included in python
    #this effectively prints all columns as tab separated

#load the table
smartpca_eigenvectors = pd.read_csv(\
    "./data/final_pcas/02_smartpca/merged_3_geno_only_typed_ld_prunning_tab.eigs", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(smartpca_eigenvectors)
#check we have the correct number of columns
if(smartpca_eigenvectors.shape[1]-1!=n_axes):
    raise ValueError("ERROR! FALSE! We are not considering all the interesting axes in smartpca")

#plot the ifrst 5 axes againts each other
import matplotlib.pyplot as plt
import seaborn as sns
sns.pairplot(smartpca_eigenvectors.iloc[:,1:6])
plt.savefig("./data/final_pcas/02_smartpca/eigen_vectors_pairwise.png", dpi=300)
plt.close()

print_text("plot the explained variance using eigenvalues", header=3)
print_text("load the eigenvalues", header=4)
smartpca_eigenvalues = pd.read_csv(\
    "./data/final_pcas/02_smartpca/merged_3_geno_only_typed_ld_prunning.eval", \
    sep="\t", \
    header=None, \
    low_memory=False)
print(smartpca_eigenvalues)
    #I have checked these eigenvalues are the same than in the first row of the eigenvector table, which is named as "#eigvals"

print_text("check eigenvalues are in decreasing order", header=4)
is_decreasing = smartpca_eigenvalues.sort_values(by=0, ascending=False).equals(smartpca_eigenvalues)
if(is_decreasing == False):
    raise ValueError("ERROR! FALSE! The eigen values are not in decreasing order")

print_text("The proportion of total variance explained by each PC is a useful metric for understanding structure in a sample and for evaluating how many PCs one might want to include in downstream analyses (see Fig. 9). This can be computed as λi∕∑kλk, with λi being eigenvalues in decreasing order", header=4)
smartpca_eigenvalues[1] = smartpca_eigenvalues[0]/smartpca_eigenvalues[0].sum()
#smartpca_eigenvalues[2] = smartpca_eigenvalues.iloc[:, 1].cumsum(axis=0)
#check
if(smartpca_eigenvalues[1].sum() != 1.0):
    raise ValueError("ERROR! FALSE! We have not correctly calculated the proportion of explained variance")

print_text("change column names", header=4)
smartpca_eigenvalues = smartpca_eigenvalues.rename(columns={0:"eigenvalue", 1:"prop_variance"})
print(smartpca_eigenvalues)
    
print_text("select the eigenvalues of the relevant axes", header=4)
smartpca_eigenvalues_20 = smartpca_eigenvalues.iloc[0:n_axes, :]
print(smartpca_eigenvalues_20)

print_text("add a new column with the index (row numbers) and plot index vs the first and only column (eigenvalues)", header=4)
import matplotlib.pyplot as plt
smartpca_eigenvalues_20.reset_index().plot( \
    x='index', \
    y="prop_variance", \
    style='-o', \
    legend=None)
    #plot a line with a dot in eahc observation
    #reset_index() is called on the pca_eigenvalues DataFrame. This resets the index of the DataFrame, and the old index is added as a column named 'index'. This new DataFrame is then plotted using the plot function, with 'index' as the x-values and the first column (0) as the y-values.

#add labels
plt.title('Scatter Plot of PCA Eigenvalues')
plt.xlabel('PCs')
plt.ylabel('Proportion of variance explained')

#add xticks from the first to the last PC
plt.xticks( \
    ticks=smartpca_eigenvalues_20.index, \
    labels=range(1, smartpca_eigenvalues_20.shape[0]+1), \
    fontsize=8 \
)
    #the positions are just the index of the eigenvalues
    #the labels are 1 to 20, but add 1 at the end because the last element of the range is not included

#save and close
plt.savefig( \
    fname="./data/final_pcas/02_smartpca/explained_var.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()
    
#Results
    #PC1 explains 1.4% of variability, while PC2 explains 0.17%. From there, next PCAs until the number 20 explains around 0.11%. Therefore, the main reduction of explained variability occurs from PC1 to PC2 and a from PC2 to PC3. From there, the reductions are much smaller. Therefore, the main axes seems to be PC1 and PC2.
    #Note that the total explained variance is much less compared to the results of plink where variance explained by the first axis was of 40%! but it is in line what the admixture tutorial got with smartpca, so maybe this is a question of this approach and/or the removal of outliers
    #we are not calculating cumultaive percentage of explained variance because we have a small proportion explained by each axis, and the rule of Ritche of 80% would require to include hundreds of axes...

print_text("plot SNP weights", header=3)
print_text("process the file to ensure it is tab delimited", header=4)
run_bash(" \
    cd ./data/final_pcas/02_smartpca/; \
    awk \
        'BEGIN{ \
            OFS=\"\t\" \
        }{ \
            if(NR>0){ \
                print " + ",".join(f"${i}" for i in range(1,n_axes+4)) + " \
            } \
        }'\
        ./merged_3_geno_only_typed_ld_prunning.snpeigs > ./merged_3_geno_only_typed_ld_prunning_tab.snpeigs \
")
    #we to sum 1 because the last number is not included in the range
    #also 3 more because we have a column with the IDs and another with the chromosome and the phyisical position (I checked with the bim file, identical all cases)

print_text("load the file", header=4)
smartpca_snp_weights = pd.read_csv(\
    "./data/final_pcas/02_smartpca/merged_3_geno_only_typed_ld_prunning_tab.snpeigs", \
    sep="\t", \
    header=None, \
    low_memory=False)
smartpca_snp_weights
#check
if(smartpca_snp_weights.shape[1] != n_axes+3):
    raise ValueError("ERROR! FALSE! We have a problem with the number of axes in the SNP weights")
    #besides axes, we have 3 columns (SNP id, chromosome and position)

print_text("sort the DF by chromosome and position", header=4)
if(smartpca_snp_weights.sort_values(by=[1, 2]).equals(smartpca_snp_weights) != True):
    raise ValueError("ERROR! FALSE! Table with snp weights is not sorted by position and chromsome")

print_text("change column names", header=4)
smartpca_snp_weights.columns = ["rs_number", "chr", "pos"] + [f"pca_{i}" for i in range(1,n_axes+1)]
    #the axes are named as pca_1, pca_2... until pca_n_axes
print(smartpca_snp_weights)

print_text("melt the DataFrame to have one row for each combination of chromosome and PCA axis", header=4)
smartpca_snp_weights_melted = pd.melt( \
    smartpca_snp_weights, \
    id_vars=["rs_number", "chr", "pos"], \
    value_vars=[f"pca_{i}" for i in range(1,n_axes+1)], \
    var_name="pca_axis", \
    value_name="snp_weight" \
)
    #This function is useful to massage a DataFrame into a format where one or more columns are identifier variables (`id_vars`), while all other columns, considered measured variables (`value_vars`), are "unpivoted" to the row axis, leaving just two non-identifier columns, 'variable' and 'value'. In other words, pd.melt transforms the DataFrame from wide format to long format, with one row for each combination of chromosome, rs_numbre, position, and PCA axis. 
        #The id_vars parameter specifies the columns to use as identifier variables. 
            #The specific combination of chromosome, SNP and position
        #the value_vars parameter specifies the columns to unpivot.
            #we are taking all the PCA_axes columns and put them all together in the same column as different rows
        #The var_name and value_name parameters specify the names of the new columns.
print(smartpca_snp_weights_melted)

print_text("check we have all the rows", header=4)
if(smartpca_snp_weights[["rs_number", "chr", "pos"]].shape[0]*n_axes!=smartpca_snp_weights_melted.shape[0]):
    raise ValueError("ERROR! FALSE! Problem when melting the DF with SNP weights, we do not have the expected number of rows")
    #52989 SNPs times 20 PCA axes makes 1059780 rows, which is the number of rows of smartpca_snp_weights_melted

print_text("get unique PCA axes and chromosomes", header=4)
pca_axes = smartpca_snp_weights_melted['pca_axis'].unique()
chromosomes = smartpca_snp_weights_melted['chr'].unique()

print_text("create subplots", header=4)
#define dimensions
import matplotlib.pyplot as plt
fig, axes = plt.subplots( \
    nrows=len(pca_axes), \
    ncols=len(chromosomes), \
    figsize=(4*len(chromosomes), 4*len(pca_axes)))
    #define the number of rows/columns of the subplot grid.

print_text("make the plots", header=4)
#for each PCA axis (i starts at 0)
#pca_axis="pca_1"; chromosome=1; i=0; j=0
for i, pca_axis in enumerate(pca_axes):
    
    #for each chromosome (j starts at 0)
    for j, chromosome in enumerate(chromosomes):
        
        #select data for this PCA axis and chromosome
        snp_weights_subset = smartpca_snp_weights_melted[(smartpca_snp_weights_melted["pca_axis"] == pca_axis) & (smartpca_snp_weights_melted["chr"] == chromosome)]
        
        #create scatter plot
        axes[i, j].scatter(x=snp_weights_subset["pos"], y=snp_weights_subset["snp_weight"])
        
        #set x-axis range extending a bit from both sides
        axes[i, j].set_xlim(snp_weights_subset["pos"].min()*0.9, snp_weights_subset["pos"].max()*1.1)
        
        #set title
        axes[i, j].set_title(f"PCA Axis: {pca_axis}, Chromosome: {chromosome}")

#adjust layout
plt.tight_layout()

#save and close
plt.savefig( \
    fname="./data/final_pcas/02_smartpca/snps_weights_across_axes_chromosomes.png", \
    dpi=300) #dpi=300 for better resolution
plt.close()

#results
    #it is useful to inspect the PC loadings to ensure that they broadly represent variation across the genome, rather than one or a small number of genomic regions. SNPs that are selected in the same direction as genome-wide structure can show high loadings, but what is particularly pathological is if the only SNPs that show high loadings are all concentrated in a single region of the genome, as might occur if the PCA is explaining local genomic structure (such as an inversion) rather than population structure.
    #Everything clean

# endregion






############################
# region FINAL CONCLUSIONS #
############################
print_text("FINAL CONCLUSIONS", header=1)
"""
I have compared the results of the PCA with the ones from the QC and they are very similar. There are some differences regarding the number of significant axes, but if you check the the PCA plots, you can see the same patterns showing again there are no clear subgroups.

Because of this, we are not going to run Admixture again, we will just take these new PCAs for the PRS.
"""

# endregion





print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/other_projects/australian_army_bishop/heavy_analyses/australian_army_bishop/association_studies/
#chmod +x ./scripts/production_scripts/03aa_final_pca_calc.py
#singularity exec ./03a_association_analyses.sif ./scripts/production_scripts/03aa_final_pca_calc.py > ./03aa_final_pca_calc.out 2>&1