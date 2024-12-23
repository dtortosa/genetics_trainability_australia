#!/bin/bash 
	#to run this script: chmod +x script.sh; ./script.sh
	#!/bin/sh does not work with my terminal en msi of David.
	#if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
	#you can save the output and the errors
		#./script.sh > script.out #only output
		#./script.sh 2> error.out #only error
		#./script.sh > script.out 2> error.out #both in different files
		#./script.sh > script.out 2>&1 #both in the same file
		#https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/

##compare the final VCF files pre-imputation between versions

#we are assumming only autosomals
count_v1=$(ls -1 ./04_fifth_step_v1 | wc -l)
count_v2=$(ls -1 ./04_fifth_step_v2 | wc -l)
if [[ $count_v1 != 22 || $count_v2 != 22 ]]; then
    exit 1;
fi;

#for 1 to 22#!/bin/bash 
	#to run this script: chmod +x script.sh; ./script.sh
	#!/bin/sh does not work with my terminal en msi of David.
	#if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
	#you can save the output and the errors
		#./script.sh > script.out #only output
		#./script.sh 2> error.out #only error
		#./script.sh > script.out 2> error.out #both in different files
		#./script.sh > script.out 2>&1 #both in the same file
		#https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/

##compare the final VCF files pre-imputation between versions

#we are assumming only autosomals
count_v1=$(ls -1 ./04_fifth_step_v1 | wc -l)
count_v2=$(ls -1 ./04_fifth_step_v2 | wc -l)
if [[ $count_v1 != 22 || $count_v2 != 22 ]]; then
    exit 1;
fi;
    #if we do NOT have 22 files, i.e., 22 chromosomes, stop execution

#for 1 to 22
for i in {1..22}; do 
    cmp \
        --silent \
        <( \
            bcftools view \
                --no-header \
                ./04_fifth_step_v1/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated_flipped_sorted_chr${i}.vcf.gz \
        ) \
        <( \
            bcftools view \
                --no-header \
                ./04_fifth_step_v2/loop_maf_missing_2_pca_not_outliers_sex_full_clean_hwe_updated_chr_autosomals_miss_maf_snp_clean_sample_miss_clean-updated_flipped_sorted_chr${i}.vcf.gz \
        );
    if [[ $? -ne 0 ]]; then 
        echo "error" ${i}
    else
        echo "ok chr" ${i}
    fi;
done;
    #we have to compare only the genotypes because the header includes time stamps, so it changes with each run. Because of that we use bcftools view -H, to see only genotypes.
    #if the status after the run is not 0, means we have differences
