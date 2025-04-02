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
        ./data/final_pcas/typed_variants_after_qc/list_variants/ \
")


print_text("define function to extract the genotypes SNPs across chromosomes", header=2)
#chromosome=21
def genotyped_snps_prep(chromosome):

    #select only TYPED variants and update their IDs using position, then create plink files
    run_bash(" \
        cd ./data/final_pcas/typed_variants_after_qc/list_variants/; \
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
        #query the output to get the ID (based on chromosome, position and alleles), TYPED and R2 columns
            #note that we assesed and solved strand problems and ensure we had the same SNPs than in the imputation panel and even switch REF-ALT based on the imputation preparion tool, so the REF/ALT alleles should be ok.
        #then ensure that the columns are tab delimitted and save
        #we have applied the same steps for the last VCF files generated in this project, so we are sure we are using the last, cleaned data

    #check we correctly applied the TYPED flag and the imputation quality filter
    run_bash("\
        cd ./data/final_pcas/typed_variants_after_qc/list_variants/; \
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
tsv_files = natsorted([f for f in os.listdir("./data/final_pcas/typed_variants_after_qc/list_variants/") if f.endswith("_snp_list.tsv")])
if(len(tsv_files)!=22):
    raise ValueError("ERROR! FALSE! NOT ALL CHROMSOMES ARE RUNNING")
else:
    print(tsv_files)

print_text("sequentially merge the .tsv files", header=3)
merged_df = pd.DataFrame()
#tsv_file=tsv_files[0]
for tsv_file in tsv_files:
    file_path = os.path.join("./data/final_pcas/typed_variants_after_qc/list_variants/", tsv_file)
    temp_df = pd.read_csv(file_path, sep="\t", header=None)
    merged_df = pd.concat([merged_df, temp_df], ignore_index=True)
merged_df.columns = ["ID", "TYPED", "R2"]

print_text("save the full DF and also just the IDs column", header=3)
merged_df.to_csv("./data/final_pcas/typed_variants_after_qc/merged_snps_list.tsv", sep="\t", index=False, header=True)
merged_df["ID"].to_csv("./data/final_pcas/typed_variants_after_qc/merged_snps_list_id_only.txt", index=False, header=False)

print_text("check we have all the SNPs", header=3)
#sum the number of rows across all the chromosomes
total_n_snps = run_bash(" \
    awk \
        'END{print NR}' \
        ./data/final_pcas/typed_variants_after_qc/list_variants/*_snp_list.tsv \
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
        --extract ./data/final_pcas/typed_variants_after_qc/merged_snps_list_id_only.txt \
        --make-bed \
        --out ./data/final_pcas/typed_variants_after_qc/merged_3_geno_only_typed \
")
    #We use the last fileset just before the tests of batch effects (see line 322 in 02d_post_imputation_qc.py)
    #--extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.

print_text("check we have the correct number of SNPs", header=3)
run_bash(" \
    n_snps_pca=$( \
        awk \
            'BEGIN{FS=\"\t\"}{ \
        }END{print NR}' \
        ./data/final_pcas/typed_variants_after_qc/merged_3_geno_only_typed.bim \
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
        ./data/final_pcas/typed_variants_after_qc/merged_3_geno_only_typed.fam \
    ); \
    if [[ $n_samples_pca -ne $n_samples_plink_qc ]]; then \
        echo \"ERROR! FALSE! WE HAVE A PROBLEM WITH THE FINAL NUMBER OF SAMPLES FOR PCAS\"; \
    fi \
")

print_text("check we do NOT have a .missnp file", header=3)
run_bash(" \
    cd ./data/final_pcas/typed_variants_after_qc/; \
    if [[ -e merged_3_geno_only_typed.missnp ]]; then \
        echo 'ERROR! FALSE! FILTERING OF SNPS DID NOT WORK'; \
    else \
        echo 'FILTERING WORKED'; \
    fi \
")

print_text("see the total number of SNPs", header=3)
run_bash(" \
    cd ./data/final_pcas/typed_variants_after_qc; \
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
        ./data/final_pcas/ld_cleaning/ \
")

#do LD-subset increasing thinning
#also remove the areas indicated by Dr.Speed in his comment


# endregion