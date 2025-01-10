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

#Bash script to download files and decompress

#IMPORTANT: THIS IS A ONE TIME RUN SCRIPT! This worked while our imputation results were available in TOPMed, but this is no longer the case.

#set the wd
cd ./quality_control/data/genetic_data/quality_control/20_imputation_results/

#download the info files
curl
    -o ./01_qc_reports/
    -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1494026/bb2124a92934f68f104cce04b8c4316c40f30fa29a2afa0aa759187de2d2ff6c | bash
curl
    -o ./01_qc_reports/
    -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1494022/34d4333ac7af9d4c2d7c31a3546e71fe683ca5ad5994e7366e5f0562b122b856 | bash
    #-o: select the target directory where to save the output
    #-s: Silent mode. This option makes curl operate in silent mode, which means it won't show progress or error messages.
    #-L: If the server reports that the requested page has moved to a different location, this option will make curl redo the request on the new place.
    #|: Pipe the output of the curl command to the bash command
        #not sure what this means but I took this script from topmed download page, also used by Ritchie

#download the VCF and log files
curl
    -o ./03_compressed_vcf_files/
    -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1494028/d987347d570888509b07237a9bf205e3be2e4c732f6e4bc2a757f976425fe332 | bash
curl
    -o ./03_compressed_vcf_files/
    -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1494029/687e36150c419e7f12a307bf412289dc549406e3280551332141ab590d13f521 | bash

#decompress VCF files
for file in *.zip; do
    7z 
        e ./03_compressed_vcf_files/${file} 
        -p"8k&1zJRWVuPpdf" 
        -o"./04_uncompressed_vcf_files/"
        #e: This option stands for "extract". It extracts files from the archive without preserving the directory structure.
        #-p: This option specifies the password for the encrypted ZIP file. The password is "8k&1zJRWVuPpdf".
        #-o: This option specifies the output directory where the extracted files will be saved. In this case, the files will be extracted to the ./04_uncompressed_vcf_files/ directory.
done > ./04_uncompressed_vcf_files/output_decompression.txt

#check decompression
count_ok_decompression=$( \
    awk \
        'BEGIN{FS="\t"}{ \
            if($0 ~ /Everything is Ok/){ \
                count++ \
            } \
        }END{print count}' \
        ./04_uncompressed_vcf_files/output_decompression.txt \
)
if [[ $count_ok_decompression -eq 22 ]]; then
    echo "Decompression was successful"
else
    echo "ERROR! FALSE! Decompression failed"
fi
