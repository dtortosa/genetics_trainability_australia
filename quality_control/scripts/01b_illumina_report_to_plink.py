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



#########################################################
######## CONVERT ILLUMINA REPORT TO PLINK FORMAT ########
#########################################################



##################
#### Starting ####
##################



########################################
# define function to run bash commands #
########################################

#create a wrapper for subprocess.run in order to define a set of arguments and avoid typing them each time. We will ensure that we are using bash always and not sh.
from subprocess import run, PIPE
#command="ls"; return_value=False
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

    #if stderr is empty
    if complete_process.stderr=="":

        #print the standard output without "\n" and other characters
        print(complete_process.stdout)

        #return also the value if required
        if return_value==True:
            return complete_process.stdout
    elif ("Warning" in complete_process.stderr) | ("warning" in complete_process.stderr):

        #print the standard output without "\n" and other characters
        print(complete_process.stdout)

        #print the standard error without stopping
        print("WARNING: " + complete_process.stderr)

        #return also the value if required
        if return_value==True:
            return complete_process.stdout
    else:
        #print the standard error and stop
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

#test it
print("\n#######################################\n#######################################")
print("see working directory")
print("#######################################\n#######################################")
run_bash("pwd")
print("\n#######################################\n#######################################")
print("list files/folders there")
print("#######################################\n#######################################")
run_bash("ls")



######
# wd #
######

#import os
#print(os.getcwd())

#you can also use the "!" character, but only whithin ipython, it does not work if you run the script as a program
#message="hola mundo"
#!echo {message}
#!pwd
    #https://jakevdp.github.io/PythonDataScienceHandbook/01.05-ipython-and-shell-commands.html



#####################################
# General info about available data #
#####################################

#In general, I understand that the AGRF facility genotyped the DNA samples using illumina and then used GenomeStudio to analyze DNA quality and produce the FinalReports, being these reports the files we are using in our analyses.
    #https://www.illumina.com/Documents/products/technotes/technote_infinium_genotyping_data_analysis.pdf

#Batches
    #Within combat_genes.zip, we have two batches (ILGSA24-17303 and ILGSA24-17873), being the data separated in these two.
    #In ILGSA24-17303.zip, we have the final reports for each 216 samples, along with the sample and snp maps, while the IDAT files, DNA reports and other reports are in the folder ILGSA24-17303_idat_reports (I may have downloaded from the initial location where this data was stored in summer 2022).
    #In 17873 we have For 1248 individuals the final reports, other reports and IDAT files.

#Files we have
    #we have the IDAT files with probs intensity from the microarrays used to genotype (first zips):
        #I think these files include the image data used to do the sequencing. It seems that new illumina sequencing technology uses two color channels (red and green, which are the colors present in these zips) to determine which is present in a given position. In other words, I think these are IDAT files. A priori, we are not going to use this information, and if we need prob intensity, I think we can just use log R, so we should be fine. See next lines
            #https://www.ogc.ox.ac.uk/wp-content/uploads/2017/09/techspotlight_two-channel_sbs.pdf
    #CNMetrics (I guess copy number metrics)
        #This includes average log R and B allele frequency. This is useful for sex mismatches and copy number analyses, but we have this information in the final reports! This excel gives average and dev per sample and type of SNP (e.g., A/C), but you could calculate that information using the FinalReport as you have one value per sample and SNP. So you can do this for this and for the first batch.
            #log Ratio: 
                #We have a value per SNP, so we can get the mean/median in a subset of variants in the X chromosome to do sex checks
                    #"Examining the intensity of probe binding on the sex chromosomes will better resolve these cases. Illumina calls this intensity LogR ratio, and on Applied Biosystems (formerly Affymetrix systems), it is simply known as probe intensity. These metrics, once suitably normalized, are roughly linear in copy number. Because there are tens of thousands of loci on the X chromosome on modern platforms, it is appropriate to examine a subsample of markers and then take a measure of the central tendency of each sample, such as the median or mean intensity."
                        #Quality Control Procedures for Genome-Wide Association Studies
            #B allele frequency
                #This is a metric about the allele intensity ratio of two alleles (A and B), with 0 or 1 meaning the absence of one of the two alleles, and 0.5 means equal presence of both alleles. So this can be used to detect imbalances due to higher copy numbers due to duplications... so the frequency would not be 0, 0.5 or 1, but it would be in between.
                    #https://www.biostars.org/p/254848/#:~:text=%22The%20B%2DAllele%20Frequency%20is,alleles%20(e.g.%20AB).%22
    #DNA report
        #This file is generated with GenomeStudio prior the generation of the FinalReports. I think we can get all the information we need from the FinalReports, as we can calculate genotyping call, percentage of hetero. We even have the GenCall score (GC score) in the FinalReports so we can calculate calls vs no_calls, the percentiles of GC score...
            #https://www.illumina.com/Documents/products/technotes/technote_infinium_genotyping_data_analysis.pdf
        #No_Calls
            #Genotypes discarded
            #The total number of genotypes in each sample with a GenCall score below the no-call threshold as defined in the project options (default is 0.1500; see below). Genotypes that are not called are shown on the GenomeStudio SNP Graph as black points falling outside of the darkly shaded regions.
                #GenCall Score
                    #The first GenomeStudio parameter that should be optimized to obtain the highest genotyping accuracy is the GenCall score cutoff, or nocall threshold. The GenCall score is a quality metric calculated for each genotype (data point), and ranges from 0 to 1. GenCall scores generally decrease in value the further a sample is from the center of the cluster to which the data point is associated. The nocall threshold is the lower bound for calling genotypes relative to its associated cluster. Illumina FastTrack Genotyping Project Managers typically use a nocall threshold of 0.15 with Infinium data. This means that genotypes with a GenCall score less than 0.15 are not assigned genotypes because they are considered too far from the cluster centroid to make reliable genotype calls. They are instead assigned a “no call” for that locus. No calls on successful DNAs at successful loci contribute to lowering the call rate for the overall project. The standard 0.15 no-call threshold for Infinium data was determined using projects with trio and replicate information. A range of no-call thresholds were evaluated to optimize the call rate without compromising reproducibility or Mendelian consistency. The nocall threshold in GenomeStudio can be changed by selecting Tools | Options | Project | Nocall Threshold. Sample and SNP statistics must be recalculated after an adjustment to the nocall threshold. 
                    #This score is a quality metric that indicates the reliability of the genotypes called. A GenCall score value is calculated for every genotype and can range from 0.0 to 1.0. GenCall scores are calculated using information from the sample clustering algorithm. Each SNP is evaluated based on the angle of the clusters, dispersion of the clusters, overlap between clusters, and intensity. Genotypes with lower GenCall scores are located furthest from the center of a cluster and have lower reliability. There is no global interpretation of a GenCall score as it depends on the clustering of samples at each SNP. Clustering is affected by many different variables, including the quality of the samples and loci.
        #Calls: 
            #The total number of genotypes in each sample with a GenCall score above the no-call threshold.
        #Note about Calls - No_Calls
            #The sum of Calls + No_Calls is always 650181 in both batches
            #However, we have 654027 variants in the SNP map.
            #I think that No_Calls is NOT counting missing cases, because I guess you need to have genotype to obtain a GenCall score. Indeed, the parameter "# LOCI" in the same DNA report is 654027, so we indeed have the correct number of loci, but less genotypes, probably because of missing cases.
        #Call_Freq: 
            #Call_Freq is equal to #Calls /(#No_Calls + #Calls). Call_Freq is equivalent to Call Rate in the GenomeStudio Samples Table.
        #A/A_Freq: 
            #For each sample, the number of AA genotype calls divided by #Calls.
        #A/B_Freq: 
            #For each sample, the number of AB genotype calls divided by #Calls.
        #B/B_Freq: 
            #For each sample, the number of BB genotype calls divided by #Calls
        #Minor_Freq: 
            #If the number of AA calls is less than the number of BB calls for a sample, the frequency for the minor allele A is: [(2*AA) + AB] divided by [2*(AA+AB+BB)] across all called loci for that sample.
                #AB cases count the same for both alleles, so makes sense to use AA and BB
        #50%_GC_Score:
            #For each sample, this represents the 50th percentile of the distribution of GenCall scores across all called genotypes. For SNPs across all samples, this is referred to as the 50%_GC_Score. For samples across all loci, it is referred to as p50GC in the Samples Table.
        #10%_GC_Score: 
            #For each sample, this represents the 10th percentile of the distribution of GenCall scores across all called genotypes. For SNPs across all samples, this is referred to as the 10%_GC_Score. For samples across all loci, it is referred to as p10GC in the Samples Table. 
                #Note: Call frequency, 50% GenCall score, and 10% GenCall score are useful metrics for evaluating the quality and performance of DNA samples in an experiment.
        #0/1: 
            #GenomeStudio calculates a threshold from the distribution of 10%_GC_Score values across all samples in the DNA report. A ‘1’ is assigned to samples whose 10%_GC_Score is at or above this threshold. A ‘0’ is assigned to samples whose 10%_GC_Score is below this threshold. The equation defining this threshold is 0.85*90th percentile of 10%_GC_Score values for all samples in DNA Report (i.e., 0.85*90th percentile of column K in the DNA report).
    #SampleSheet
        #We have all the information we need about the samples in the final reports and the sample map
        #Sample_ID
            #Sample identifier (used only for display in the table).
        #Sample_Plate
            #The barcode of the sample plate for this sample (used only for display in the table).
        #Sample_Well
            #The well within the sample plate for this sample (used only for display in the table).
        #Sex
            #Male, Female, or Unknown.
        #SentrixBarcode_A
            #The barcode of the Universal Array Product that this sample was hybridized to for Manifest A
        #SentrixBarcode_A
            #The position within the Universal Array Product this sample was hybridized to for Manifest A (and similarly for _B, _C, etc. depending on how many manifests are used with your project).
        #https://www.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2011-1/genomestudio-gt-module-v1-0-user-guide-11319113-a.pdf
    #beadpool manifest (bpm file)
        #Also referred to as a SNP Manifest, this is a file containing the SNP-to-beadtype mapping, as well as all SNP annotations. For the Infinium assay, this is a *.bpm file in binary format.
            #https://www.biostars.org/p/82576/
        #this is used as input on genome studio to do quality control of sequencing and generate the final reports. Already done by AGFR.
    #ClusterFile.egt
        #This is used by genome studio
        #The cluster file contains the mean (R) and standard deviation (theta) of the cluster positions, in normalized coordinates, for every genotype, for every SNP. The cluster file also includes cluster score information, as well as the allele frequencies from the training set used to generate the cluster file. A cluster file is required for KaryoStudio. Illumina provides a standard cluster file for each product. Alternatively, customers may generate their own cluster file.
            #https://www.biostars.org/p/82576/
        #why we need a cluster?
            #Each SNP is analyzed independently to identify genotypes. Genotypes are called by comparing customer-generated data with those in the supplied cluster file. Genotype calls are highly accurate and unambiguous for high-quality samples. Generally, high-quality data with 99.5% call rates can be expected. However, accuracy is highly sample dependent. When samples do not perform as expected, experimenters can choose to reprocess (requeue) these samples to confirm or potentially improve results. Poorly performing samples can be systematically excluded from the project before recalling genotypes. Additionally, some samples (e.g., samples that have been whole-genome amplified) may not fit standard cluster positions well. A custom-generated cluster file may improve the call rate and accuracy for these aberrant cases.
            #https://www.illumina.com/Documents/products/technotes/technote_infinium_genotyping_data_analysis.pdf 
    #Reproducibility and Heritability Report.csv
        #We only have replicates in the second batch, so we do not have this file for the FIRST batch, only for the SECOND.
        #Reproducibility part is where we have data
            #Rep1_DNA_Name
                #Sample name designated as replicate #1.
            #Rep2_DNA_Name
                #Sample name designated as replicate #2
            #Correct
                #Number of loci with consistent replicate genotype comparisons
            #Errors
                #Number of loci with inconsistent replicate genotype comparisons.
            #Total
                #Number of total genotype comparisons (1 genotype comparison per locus per replicate pair). The report does not include genotypes with intensities that fall below the no-call threshold.
            #Repro_Freq
                #Reproducibility frequency. The error rate does not include genotype calls that fall below the no-call threshold
            #https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2-0/genomestudio-genotyping-module-v2-user-guide-11319113-01.pdf
        #SampleMap has a column for replicates, but this column is empty for all the 6 samples with similar IDs. Therefore, we cannot know what sample is the replicate. Also, it is worriesome that for two of the samples we also have duplicate IDs in the phenotypic data but with different phenotypes! So I will continue with the plan and remove these samples.
            #The Sample_ID of a sample that is a replicate to this sample (used in reproducibility error calculations).
                #https://www.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2011-1/genomestudio-gt-module-v1-0-user-guide-11319113-a.pdf
        #mail to bishop
            #Apologies for writing again, but I have found something else I think we should also ask to AGRF. I have found a Reproducibility report among the files of the second batch. This file includes information of samples that were genotyped twice to test the reproducibility of the results. Here we have the three samples that are duplicated (e.g., 1100JHJM_1 and 1100JHJM_2). This is still confusing to me because 1100JHJM appears two times in the phenotype file with different phenotypic values, so what sample did they sequenced twice? Also, the illumina sample map says there is no replicate of any sample, which contradicts the reproducibility report. Finally and most important, the reproducibility rate between each pair of duplicates is around 0.86. In our past genomic studies, we had a concordance rate of 99% between replicates, which is expected when you sequence the same sample twice. Maybe they calculated the concordance between samples in a different way, but I am a bit surprised about the low reproducibility.
    #plink files
        #I guess they generated ped files with all the genotypes for each batch, also the map of the snps
    #we also have binary files like Heredity.bin, Duplicates.bin... I have seen other people got the same files from illumina. 
        #If you see the linked biostars thread, there is someone with the illumina data within genome studio (software property of illumina). He has the same bin file I have, but the solution they gave him is that he needs to create the FinalReports and the Maps using genome studio in order to load the data in R and process it. We already have FinalReports and the maps, so I think we have everything we need to analyze the data.
            #https://www.biostars.org/p/2240/
    #FinalReports and snp/sample maps
        #I have received a Illumina report with 3 files ([link](https://www.biostars.org/p/51928/)):
            #The "FinalReport.txt" for Illumina raw genotype data generated from Genome Bead Studio for 2.5M (GSGT Version  2.0.4). This includes a header with number of SNPs, samples.... and then the data with sample index, sample names, alleles... the first row includes the column names. 
                #columns
                    #Sample Index    
                    #Sample ID       
                    #Sample Name     
                    #SNP Index       
                    #SNP Name        
                    #Chr     
                    #Position        
                    #GT Score
                        #Not sure what is this score     
                    #GC Score
                        #The GenCall score is a quality metric calculated for each genotype (data point), and ranges from 0 to 1. GenCall scores generally decrease in value the further a sample is from the center of the cluster to which the data point is associated.
                        #see above for further details
                    #Allele1 - AB    
                    #Allele2 - AB
                        #I guess this is only telling what allele is considered A in the notation of illumina, as 


                                #strand
            #Forward/Reverse is same as the one used in dbSNP, check if dbSNP always uses forward as reference in the last years
            #https://www.biostars.org/p/4885/
            #https://www.illumina.com/documents/products/technotes/technote_topbot.pdf
  
                    #Allele1 - Top   
                    #Allele2 - Top   
                    #Allele1 - Forward       
                    #Allele2 - Forward       
                    #Allele1 - Design        
                    #Allele2 - Design        
                    #Theta   
                    #R   

                            #The cluster file contains the mean (R) and standard deviation (theta) of the cluster positions, in normalized coordinates, for every genotype, for every SNP
    
                    #X Raw   
                    #Y Raw   
                    #X       
                    #Y       
                    #B Allele Freq   
                    #Log R Ratio     
                    #SNP Aux
                    #SNP     
                    #ILMN Strand     
                    #Top Genomic Sequence    
                    #Customer Strand

                #From this file, we can obtain the a lgen file. It is just taking the first few columns of the FinalReport. Plink has an option to read the .lgen format and convert it to PED or BED formats (see below; [link](https://www.biostars.org/p/13302/))
            #A SNP map file with the physical positions of the snps.
            #A sample map file with info for each individual, like the sex, ID..



        #try to use bmp file from genomestudio?


    #the phenotype csv has 1463 samples, not having phenotype data for 41 of them. 1248+216=1464, so we also lack 1 sample in the csv file (this is the replicate of 7800AGSO)


#IMPORTANT INFO ABOUT COORDINATE SYSTEM (1- vs. 0- based) in final reports
    #our final reports are 1-based and hg38:
        #I have checked if the physical position of three SNPs in the SNP map correspond with that in GeneBank, and that is the case. The three SNPs have the position indicated for hg38 in GeneBank, which is 1-based.
            #https://www.ncbi.nlm.nih.gov/snp/?term=rs77927848
            #https://www.ncbi.nlm.nih.gov/snp/?term=rs77928150
            #https://www.ncbi.nlm.nih.gov/snp/?term=rs77928688
        #These are also the coordinates in ensmeble. Importantly, for each SNP it is shown the VCF representation and the coordinate in that representation is also the same. Remember that VCF files are 1-based.
            #https://useast.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=6:64108056-64109056;v=rs77927848;vdb=variation;vf=185603098
            #https://useast.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=9:134543926-134544926;v=rs77928150;vdb=variation;vf=735478081
            #https://useast.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=13:72551528-72552528;v=rs77928688;vdb=variation;vf=53001601
            #https://useast.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=4:87992286-87993286;v=rs56917667;vdb=variation;vf=99078765
            #https://useast.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=3:61685181-61686181;v=rs73099117;vdb=variation;vf=99493412
        #blog post with summary of coordinate systems for several pages and programs
            #https://tidyomics.com/blog/2018/12/09/2018-12-09-the-devil-0-and-1-coordinate-system-in-genomics/
        #Therefore, we can say that the coordinates of the final reports are 1-based and hg38.

#IMPORTANT INFORMATION about integrity of the files
    #I get a warning when unzipping the zip of the second batch (CAGRF20093767.zip) with unzip, and I cannot access the data
    #with 7z I can access the data but I get a warning "headers error".

    #I am checking the integrity of the files using checksums:
        #https://www.tecmint.com/generate-verify-check-files-md5-checksum-linux/
        #A checksum is a digit which serves as a sum of correct digits in data, which can be used later to detect errors in the data during storage or transmission. MD5 (Message Digest 5) sums can be used as a checksum to verify files or strings in a Linux file system.
        #MD5 Sums are 128-bit character strings (numerals and letters) resulting from running the MD5 algorithm against a specific file. The MD5 algorithm is a popular hash function that generates 128-bit message digest referred to as a hash value, and when you generate one for a particular file, it is precisely unchanged on any machine no matter the number of times it is generated.
        #It is normally very difficult to find two distinct files that results in same strings. Therefore, you can use md5sum to check digital data integrity by determining that a file or ISO you downloaded is a BIT-FOR-BIT COPY OF THE REMOTE FILE or ISO.
        #you can do "md5sum file.txt" and get the 128-bit message. If you modify just a line, anything, from the file, you get a new 128-bit message.
        #you can obtain the message just piping:
            #echo "eso" | md5sum
            #echo "eso1" | md5sum #they have different message
        #This checks the content, but not the file name. So this would get the same checksum even having different names
            #echo "eso" > file_1; md5sum file_1
            #echo "eso" > file_2; md5sum file_2
        #You can redirect the hash value(s) of a file(s) into a text file and store, share them with others. For the two files above, you can issues the command below to redirect generated hash values into a text file for later use:
            #md5sum file_1 file_2 > checksum.md5
            #cat checksum.md5 # you get two lines, one per file, with checksum and file name
        #To check that the files have not been modified since you created the checksum, run "md5sum --check checksum.md5". You should be able to view the name of each file along with “OK”.
            #md5sum --check checksum.md5
                #file_1: OK
                #file_2: OK
        #Remember that after creating the checksum, you can not rename the files or else you get a “No such file or directory” error, when you try to verify the files with new names.
            #mv file_1 file_3
            #md5sum --check checksum.md5
                #file_1: FAILED open or read
                #file_2: OK
        #This is exactly what the PDFs about data integrity in both batches says
        #Our results
            #first batch: 
                #in "CombatGenes" run "md5sum --check checksums.md5"
                #checksums.md5: FAILED
                #example_ILGSA24-17303_FinalReport1.txt: OK
                #example_Sample_Map.txt: OK
                #example_SNP_Map.txt: OK
                #ILGSA24-17303.zip: OK
                #md5sum: WARNING: 1 computed checksum did NOT match
            #second batch
                #in "CombatGenes/17873" run "md5sum --check checksums.md5"
                #205771890087.zip: OK
                #205771890120.zip: OK
                #205771890129.zip: OK
                #205771890173.zip: OK
                #205785500018.zip: OK
                #205785500033.zip: OK
                #205785500038.zip: OK
                #205785500075.zip: OK
                #205857150085.zip: OK
                #205857150090.zip: OK
                #205857150136.zip: OK
                #205857150149.zip: OK
                #205955840060.zip: OK
                #205955840105.zip: OK
                #205955840108.zip: OK
                #205955840137.zip: OK
                #205960020112.zip: OK
                #206023350028.zip: OK
                #206023350029.zip: OK
                #206023350043.zip: OK
                #206023350156.zip: OK
                #206036460037.zip: OK
                #206036460040.zip: OK
                #206036460041.zip: OK
                #206036460042.zip: OK
                #206036460091.zip: OK
                #206036460093.zip: OK
                #206036460121.zip: OK
                #206036460124.zip: OK
                #206036460127.zip: OK
                #206036460128.zip: OK
                #206036460129.zip: OK
                #206036460164.zip: OK
                #206036460165.zip: OK
                #206036460169.zip: OK
                #206053690035.zip: OK
                #206063100046.zip: OK
                #206063100048.zip: OK
                #206063100056.zip: OK
                #206063100059.zip: OK
                #206063100076.zip: OK
                #206063100077.zip: OK
                #206063100078.zip: OK
                #206063100079.zip: OK
                #206063100080.zip: OK
                #206063100081.zip: OK
                #206063100119.zip: OK
                #206063100120.zip: OK
                #206063100131.zip: OK
                #206063100132.zip: OK
                #206123430018.zip: OK
                #md5sum: 206123430033.zip: No such file or directory
                #206123430033.zip: FAILED open or read
                #CAGRF20093767_CNMetrics.csv: OK
                #CAGRF20093767_DNAReport.csv: OK
                #CAGRF20093767_Reproducibility and Heritability Report.csv: OK
                #CAGRF20093767_SampleSheet.csv: OK
                #CAGRF20093767.zip: OK
                #checksums.md5: FAILED
                #GSA-24v3-0_A1_ClusterFile.egt: OK
                #GSA-24v3-0_A2.bpm: OK
                #PLINK_030222_0457.zip: OK
                #md5sum: WARNING: 1 listed file could not be read
                #md5sum: WARNING: 1 computed checksum did NOT match
            #In both batches we get some errors for some files, but the important thing is that the two zips we are using are EXACTLY THE SAME than when they were created in first place by the sequencing center. Therefore, all our problem transfering the data has not affected the data.
                #in the first batch we get Fail for checksum file itself (checksums.md5), but the rest of files that the checksum file targets, which are the one we are interested, are OK.
                #in the second bath we have problems again with the checksum file but also with 1 compressed file: 206123430033.zip is not present.
                    #As explained at the beginning of this script, I think we do not need these type of files (I think they are IDAT files), but just the FinalReports, so I think we are good here.

    #the check with sumchecks is great because we now know that the data I am using is exactly the same as created, but still we get the header warning in the zip of the second batch. As I said, that problem was not created by me nor David because the sumchecks are the same, we did not change anything, but maybe there was a problem when compressing the file by the sequencing center.

    #Using "7z t CombatGenes.zip", I have checked again the integrity of the whole zip and the two zips of each batch, thus checking specifically each FinalReport file: 
        #WHOLE ZIP: The result is OK
            #Scanning the drive for archives:
            #1 file, 79228316574 bytes (74 GiB)
            #
            #Testing archive: CombatGenes.zip
            #--
            #Path = CombatGenes.zip
            #Type = zip
            #Physical Size = 79228316574
            #64-bit = +
            #
            #Everything is Ok                                 
            #
            #Folders: 2
            #Files: 71
            #Size:       79228299532
            #Compressed: 79228316574
        #FIRST BATCH: The result is OK
            #Scanning the drive for archives:
            #1 file, 5664833599 bytes (5403 MiB)
            #
            #Testing archive: ILGSA24-17303.zip
            #--
            #Path = ILGSA24-17303.zip
            #Type = zip
            #Physical Size = 5664833599
            #64-bit = +
            #
            #Everything is Ok                                         
            #
            #Folders: 1
            #Files: 218
            #Size:       19929504446
            #Compressed: 5664833599
        #SECOND BATCH: ONLY 1 ERROR WITH THE HEADER
            #Scanning the drive for archives:
            #1 file, 36627556056 bytes (35 GiB)
            #
            #Testing archive: CAGRF20093767.zip
            #             
            #ERRORS:
            #Headers Error
            #
            #--
            #Path = CAGRF20093767.zip
            #Type = zip
            #ERRORS:
            #Headers Error
            #Physical Size = 36627556056
            #64-bit = +
            #
            #Archives with Errors: 1
            #
            #Open Errors: 1

    #I have written a question to "bioinformatics@agrf.org.au" to check this is ok:
        #I would like to ask two questions regarding the project "CAGRF20093767":
            #I have been checking the integrity of the specific compressed file containing the "FinalReports". The checksums are OK, but I get a warning when using "7z t": "Headers Error". I understand this is related to the signature header in the zip file, but should we worry about this? All FinalReports within the compressed file are OK (both according to 7z and checksums) and I have full access to them. 
            #When looking at the checksums for the whole "CAGRF20093767" project, all files are OK with two exceptions: there is one file that is not found ("206123430033.zip"), while "checksums.md5" failed. I am not particularly worried about the missing zip as I think I have all the information I need in the FinalReports, but I would like to check that it is ok for "checksums.md5" to fail if the rest of the files are OK.



#######################################
# Passing arugments of python program #
#######################################

#define input arugments to be passed when running this script in bash
#we will use a bash script to run this python program two times instead of doing a function and parallelize that function. In this way, the same python script is run for each batch and if we want to run again one batch, we just need to modify the bash script, not the python program. This makes sense because we have only two batches, and the parallelization will occur inside each batch. Importantly, we can also have separated .out files, so we can look at the checks separately.
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--batch_name", type=str, default="ILGSA24-17873", help="Name of the batch used as input. Always string.")
parser.add_argument("--n_cores", type=int, default=4, help="Number of cores/threads requested. Integer always, None does not work!")
parser.add_argument("--n_samples", type=int, default=10, help="Number of samples to be analyzed. Integer always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
batch_name = args.batch_name
n_cores = args.n_cores
n_samples = args.n_samples

#stop if the number of samples is None
if n_samples == None:
    raise ValueError("ERROR! FALSE! --n_samples is None but this script is not ready to deal with None as n_samples")

#starting
print("#################################################################################################################################\n#################################################################################################################################")
print("############################################# STARTING BATCH NUMBER " + batch_name + " USING " + str(n_cores) + " CORES AND " + str(n_samples) + " SAMPLES ##############################################")
print("#################################################################################################################################\n#################################################################################################################################")



###############################################
#### Unzip files of the batch in temp dict ####
###############################################

#set the name of the zip file to be unzzipped for each batch
if batch_name=="ILGSA24-17303":
    zip_name="ILGSA24-17303"
elif batch_name=="ILGSA24-17873":
    zip_name="CAGRF20093767"

#create temporary folder to save
import tempfile
temp_dir = tempfile.TemporaryDirectory()
#print(temp_dir.name)

#to remove the dir
#temp_dir.cleanup()
    #https://stackoverflow.com/questions/3223604/how-to-create-a-temporary-directory-and-get-its-path-file-name

#read only the zip and get list of files in it
import numpy as np
import zipfile
zipdata = zipfile.ZipFile("data/genetic_data/illumina_batches/" + zip_name + ".zip")
zipinfos = zipdata.infolist()

#select only FinalReports
zipinfos_subset = [zipinfo for zipinfo in zipinfos if zipinfo.filename.startswith(zip_name+"/"+zip_name+"_FinalReport")]

#reorder the zipinfos by natural order, 1, 2, ...10, 11...
names_to_sort = [zipinfo.filename for zipinfo in zipinfos_subset]
    #get the names of each file inside zipinfos_subset
from natsort import index_natsorted, order_by_index
index_order = index_natsorted(names_to_sort)
    #get the indexes under natural order as done by natsorted
zipinfos_subset = order_by_index(zipinfos_subset, index_order)
    #use these indexes to order the list of zipinfos

#if we have selected only a subset of samples
if (n_samples != None) and (batch_name=="ILGSA24-17303" and n_samples<216) | (batch_name=="ILGSA24-17873" and n_samples<1248):

    #select a subset of the samples
    zipinfos_subset = zipinfos_subset[0:n_samples]

#once we have the selected samples, add the paths for SNP and sample maps obtained from the original zipinfos
zipinfos_subset.append(zipinfos[np.where([zipinfo.filename == zip_name + "/SNP_Map.txt" for zipinfo in zipinfos])[0][0]])
zipinfos_subset.append(zipinfos[np.where([zipinfo.filename == zip_name + "/Sample_Map.txt" for zipinfo in zipinfos])[0][0]])

#get the name of each zip file
names_files = [zipinfo.filename for zipinfo in zipinfos_subset]

#check
print("\n#######################################\n#######################################")
print("Do we have selected the correct number of files to be decompressed according to arg --n_sample?")
print("#######################################\n#######################################")
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(len(names_files) == 216 + 2)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(len(names_files) == 1248 + 2)
else:
    print(len(names_files) == n_samples + 2)

#iterate across files to decompress
print("\n#######################################\n#######################################")
print("Unzipping data: ")
print("#######################################\n#######################################")
#we are using a loop because it is a bit faster than zipdata.extractall. Parallelizing does not seems to be a good solution for unzipping and I got problems with pool. This takes a few minutes anyways.
    #parallel slower
        #https://stackoverflow.com/questions/43313666/python-parallel-processing-to-unzip-files
    #count time
        #https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution
import numpy as np
#zipinfo=zipinfos_subset[1]
#zipinfo=zipinfos_subset[np.where([element == zip_name + "/SNP_Map.txt" for element in names_files])[0][0]]
#zipinfo=zipinfos_subset[np.where([element == zip_name + "/Sample_Map.txt" for element in names_files])[0][0]]
for zipinfo in zipinfos_subset:

    #if the file name starts with the adequate zip (batch) name and FinalReport
    if (zipinfo.filename.startswith(zip_name + "/" + zip_name + "_FinalReport")) | (zipinfo.filename.startswith(zip_name + "/SNP_Map.txt")) | (zipinfo.filename.startswith(zip_name + "/Sample_Map.txt")):

        #rename by removing the name of the parent folder, so we only get the file not the parent folder
        zipinfo.filename = zipinfo.filename.split(zip_name+"/")[1]
        print(zipinfo.filename)

        #extract the file in the temp dict
        zipdata.extract(zipinfo, temp_dir.name)
            #https://stackoverflow.com/questions/44079913/renaming-the-extracted-file-from-zipfile/56362289#56362289
            #https://stackoverflow.com/questions/13765486/how-to-unzip-specific-folder-from-a-zip-with-python

        #if the selected file is a FinalReport
        if zipinfo.filename.startswith(zip_name + "_FinalReport"):

            #remove the first 10 lines of the file, which is the header
            run_bash(" \
                cd " + temp_dir.name + "; \
                tail -n +11 " + zipinfo.filename + " > tmp.txt && mv tmp.txt " + zipinfo.filename)
                #if you use tail with "+number_line" you can get all lines starting from the selected number of line
                    #tail is faster than sed
                #then save the result in a temporal file and overwrite the original file with this new file without the header. The && will make sure that the file doesn't get overwritten when there is a problem.
                    #https://stackoverflow.com/questions/339483/how-can-i-remove-the-first-line-of-a-text-file-using-bash-sed-script
                    #https://www.baeldung.com/linux/remove-first-line-text-file

#count number of files and get the value to do check
n_files = int(run_bash("ls " + temp_dir.name + " | wc -l", return_value=True))

#check
print("\n#######################################\n#######################################")
print("the number of files extracted is correct according to arg --n_sample?")
print("#######################################\n#######################################")
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(n_files == 216 + 2)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(n_files == 1248 + 2)
else:
    print(n_files == n_samples + 2)

#list files present in temp
import glob
list_files = glob.glob(temp_dir.name + "/*")

#natural sorting, 1,10, 21.. that works with numbers + strings like 1b, 5c...
from natsort import natsorted
list_files = natsorted(list_files)
    #https://github.com/SethMMorton/natsort#installation

#check
print("\n#######################################\n#######################################")
print("the number of listed files is correct according to arg --n_sample?: ")
print("#######################################\n#######################################")
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(len(list_files) == 216 + 2)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(len(list_files) == 1248 + 2)
else:
    print(len(list_files) == n_samples + 2)

#extract the names of the final reports only, which start with the zip name
list_files_samples = [file for file in list_files if file.startswith(temp_dir.name + "/" + zip_name + "_FinalReport")]



###################################################
#### read illumina reports and maps with spark ####
###################################################

print("\n#####################\n#####################")
print("read illumina reports and maps with spark")
print("#####################\n#####################")

#spark gives us the possibility to analyze big datasets with much less ram and relatively fast. You can work with multiple final reports, being all connected, making possible to make queries on them, but not being all loaded in memory

#the header of the final reports has been removed, so now the first row has the column names of the genotype data, but no illumina info about the report

#set the number of threads per core
if n_cores==None:
    n_threads="*" #this gets all logical cores available (see below)
else:
    n_threads=str(n_cores)

#open a Spark session
#sparkcontext, which is used by TDI, is going to be deprecated
from pyspark.sql import SparkSession
spark = SparkSession.builder \
    .appName("working with illumina reports") \
    .master("local[" + n_threads + "]") \
    .config("spark.driver.memory", "10g") \
    .getOrCreate()
    #master: 
        #master indicates if you are in local machine or in a cluster. If cluster, you have to indicate the address. We are going to use the singularity container as a local cluster both here and in the HPC.
        #local[*] indicates that all logical cores, i.e., threads will be used. If the cores have multithreading, you can have two threads per core, so 12 cores would mean 24 threads. 
        #.master("local[1]") indicates the number of threads/logical cores manually, 1 in this case
        #I am not sure about the number of threads per core in the HPC, so I am going to set this as *.
        #Note that the number of threads impacts the number of partitions in which the data is split
        #links
            #https://stackoverflow.com/questions/32356143/what-does-setmaster-local-mean-in-spark
            #https://sparkbyexamples.com/spark/what-does-setmaster-local-mean-in-spark/
            #https://stackoverflow.com/questions/24622108/apache-spark-the-number-of-cores-vs-the-number-of-executors
    #.config() to change configuration
        #https://stackoverflow.com/questions/41886346/spark-2-1-0-session-config-settings-pyspark 
    #.config("spark.driver.memory", "6g") \
        #increase the memory for driver, if not, we get refused connection error, probably because java is killed due to out of memory problems
        #if you get more problems about JAVA stops, you can try to increase the memory for the per executor, where the operations are done.
        #https://stackoverflow.com/questions/49995044/bizarre-exception-from-pyspark-of-local-mode-in-docker-container
        #https://stackoverflow.com/questions/26562033/how-to-set-apache-spark-executor-memory
    #in case "refused connection" error, you have to commenting the first two lines "/etc/hosts", the ones with 127.0.0.1 IP. Indeed sparks "prefers" not to use this IP
        #https://stackoverflow.com/questions/24881757/apache-spark-connection-refused-for-worker

#we are going to work with data frames
    #https://spark.apache.org/docs/2.2.0/sql-programming-guide.html#spark-sql-dataframes-and-datasets-guide

#we need the same column order in all illumina reports in order to use all of them with spark. For that, extract the column names for each report using pandas. Only header
# Get first the paths for each final report
list_reports_files_full_path = glob.glob(temp_dir.name+ "/" + zip_name + "_FinalReport*.txt") 
    #I prefer using the glob module, as it does pattern matching and expansion.
        #https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
#read header of each report
import pandas as pd
list_df_reports = [
    pd.read_csv(
        path, 
        nrows=0, #no rows
        header=0, #only header
        sep='\t',
        low_memory=False) #to avoid problems with mixed types
    for path in list_reports_files_full_path]
print("\n#####################\n#####################")
print("all the columns in each report are the same than in the first report?")
print("#####################\n#####################")
print(all([(np.array(report.columns) == np.array(list_df_reports[0].columns)).all() for report in list_df_reports]))

#define the scheme, i.e., the type of the data
#YOU NEED THE COLUMN NAMES IN THE SAME ORDER THAN IN THE FILES
from pyspark.sql.types import StructType, IntegerType, StringType, DoubleType, BooleanType
schema_reports = StructType() \
    .add("Sample Index",IntegerType(), nullable=True) \
    .add("Sample ID",StringType(), nullable=True) \
    .add("Sample Name",StringType(), nullable=True) \
    .add("SNP Index",IntegerType(), nullable=True) \
    .add("SNP Name",StringType(), nullable=True) \
    .add("Chr",StringType(), nullable=True) \
    .add("Position",IntegerType(), nullable=True) \
    .add("GT Score",DoubleType(), nullable=True) \
    .add("GC Score",DoubleType(), nullable=True) \
    .add("Allele1 - AB",StringType(), nullable=True) \
    .add("Allele2 - AB",StringType(), nullable=True) \
    .add("Allele1 - Top",StringType(), nullable=True) \
    .add("Allele2 - Top",StringType(), nullable=True) \
    .add("Allele1 - Forward",StringType(), nullable=True) \
    .add("Allele2 - Forward",StringType(), nullable=True) \
    .add("Allele1 - Design",StringType(), nullable=True) \
    .add("Allele2 - Design",StringType(), nullable=True) \
    .add("Theta",DoubleType(), nullable=True) \
    .add("R",DoubleType(), nullable=True) \
    .add("X Raw",IntegerType(), nullable=True) \
    .add("Y Raw",IntegerType(), nullable=True) \
    .add("X",DoubleType(), nullable=True) \
    .add("Y",DoubleType(), nullable=True) \
    .add("B Allele Freq",DoubleType(), nullable=True) \
    .add("Log R Ratio",DoubleType(), nullable=True) \
    .add("SNP Aux",DoubleType(), nullable=True) \
    .add("SNP",StringType(), nullable=True) \
    .add("ILMN Strand",StringType(), nullable=True) \
    .add("Top Genomic Sequence",StringType(), nullable=True) \
    .add("Customer Strand",StringType(),nullable=True)
        #Chr has to be string because in case we have Y, X or other stuff
        #nullable means you can have null values
        #https://spark.apache.org/docs/3.1.3/api/python/reference/api/pyspark.sql.types.StructType.html

#create a spark DF with all the samples
df_samples = spark \
    .read \
        .option("delimiter", "\t") \
        .option("header", True) \
        .schema(schema_reports) \
        .csv(temp_dir.name + "/" + zip_name + "_FinalReport*.txt")
    #you can also use a previous list with paths to read
        #list_files_samples
    #Note this DF is NOT in memory, but you can make queries to go through the whole data while using few RAM
    #read using tab, using first row as header and then go to specific folder and select all files with the batch name and being final reports
    #use the previously defined schema
    #IMPORTANT, all the files have to have the same schema to work in this way. We already checked that all final reports are the same.
        #https://stackoverflow.com/questions/69350586/pyspark-read-multiple-csv-files-at-once

#check
print("\n#####################\n#####################")
print("check the column names of the reports DF are correct")
print("#####################\n#####################")
print(all(df_samples.schema.names == list_df_reports[0].columns))
    #this is an important step because if you have a wrong order in the schema, the names of the columns are not going to follow the order in the data. For example, if Sample ID is not first in schema input, then the first column with the sample ids will have another name.

#see schema and column names separately
print("\n#####################\n#####################")
print("see schema")
print("#####################\n#####################")
df_samples.printSchema()
#print(df_samples.schema.names)
    #https://stackoverflow.com/questions/39746752/how-to-get-name-of-dataframe-column-in-pyspark

#you can do operations on the columns
#df_samples.select(df_samples["Position"]+1).show()

#select only columns of interest
print("\n#####################\n#####################")
print("select only columns of interest")
print("#####################\n#####################")
df_samples_subset_raw = df_samples.select(["Sample Index", "Sample ID", "SNP Index", "SNP Name", "Chr", "Position", "Allele1 - Forward", "Allele2 - Forward"])
df_samples_subset_raw.show()

#add a column with the input file name
print("\n#####################\n#####################")
print("add a column with the input file name")
print("#####################\n#####################")
from pyspark.sql.functions import input_file_name
df_samples_subset = df_samples_subset_raw.withColumn("input_file", input_file_name())
df_samples_subset.show()

#load the maps
print("\n#####################\n#####################")
print("load SNP and Sample maps")
print("#####################\n#####################")
snp_map = pd.read_csv(temp_dir.name+ "/" + "SNP_Map.txt",
    delimiter="\t",
    header=0,
    low_memory=False)
sample_map = pd.read_csv(temp_dir.name+ "/" + "Sample_Map.txt",
    delimiter="\t",
    header=0,
    low_memory=False)
print(snp_map)
print(sample_map)

#you can also do SQL queries
#SQL vs pandas in spark, not great differences
    #https://www.confessionsofadataguy.com/dataframes-vs-sparksql-which-one-should-you-choose/
# Register the DataFrame as a SQL temporary view
#df_samples_subset.createOrReplaceTempView("df_samples_subset_sql")
    #do the query, asking for rows with chromosome 1
#spark \
#    .sql("SELECT * FROM df_samples_subset_sql WHERE Chr=1") \
#    .show()

#reorder using chromosome inverse
#df_samples_subset.orderBy("Chr", ascending=False).show()
    #to use in specific order
        #https://stackoverflow.com/questions/54071665/pyspark-read-multiple-csv-files-into-a-dataframe-in-order

#calculate the number of reports, i.e., samples, as the number of unique input_file names
n_unique_samples = df_samples_subset \
    .select("input_file") \
    .distinct() \
    .count()

#checks
print("\n#####################\n#####################")
print("Do we have the expected number of samples in spark as indicated with arg --n_sample?")
print("#####################\n#####################")
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(n_unique_samples == 216)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(n_unique_samples == 1248)
else:
    print(n_unique_samples == n_samples)
    #If no sample number is selected, we should have all the samples in both batches
print("\n#####################\n#####################")
print("We have as many unique samples in illumina reports as samples in sample map?")
print("#####################\n#####################")
print(n_unique_samples == sample_map.shape[0])
print(len(list_files_samples) == sample_map.shape[0])

#count the number of genotypes, i.e., total number of rows
n_genotypes = df_samples_subset.count()

#check
print("\n#####################\n#####################")
print("We have as many genotypes as snps*samples in the maps")
print("#####################\n#####################")
print(n_genotypes == snp_map.shape[0]*sample_map.shape[0])
print("\n#####################\n#####################")
print("We have as many genotypes as snps*samples according to arg --n_samples")
print("#####################\n#####################")
if (n_samples == None) or (n_samples==216 and batch_name=="ILGSA24-17303") or (n_samples==1248 and batch_name=="ILGSA24-17873"):
    print(n_genotypes == snp_map.shape[0]*sample_map.shape[0])
else:
    print(n_genotypes == snp_map.shape[0]*n_samples)    

#checks
print("\n#####################\n#####################")
print("Multiple checks within each sample")
print("#####################\n#####################")

#check that the sample ID is the same than that showed in the input file name
import pyspark.sql.functions as F
#define the two columns to be compared
col_input_name_check = F.split(F.col("input_file"), "file:" + temp_dir.name + "/" + zip_name + "_").getItem(1)
    #prepare a column from splitting the input_file column using the whole name until FinalReport as delimiter. Therefore, the split will have two parts, being FinalReportXX.txt the second, i.e., "getItem(1)". The column "input_file" was previously created.
col_index_check = F.concat(F.lit("FinalReport"), F.col("Sample Index"), F.lit(".txt"))
    #prepare a column concatenating FinalReport and .txt as literal values to the sample index, so we have the same format than in the input file name
    #in this way, we can check whether the sample index of each row corresponds with the sample index the input file.
df_samples_subset \
    .withColumn("check_1", col_index_check == col_input_name_check) \
    .select("check_1") \
    .distinct() \
    .show()
    #to the orignal data.frame, add a new column comparing the two previously defined columns, select the column with the check, get the distinct values and show.
    #Only true should be present.

#check we have the correct number of samples
n_distinct_samples = df_samples_subset \
    .select(df_samples_subset["Sample ID"]) \
    .distinct() \
    .count()
print("CHECK 2:") 
print(n_distinct_samples == sample_map.shape[0])
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(n_distinct_samples == 216)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(n_distinct_samples == 1248)
else:
    print(n_distinct_samples == n_samples)

#check we do not have duplicated SNP names
n_distinct_snp_names = df_samples_subset \
    .select(df_samples_subset["SNP Name"]) \
    .distinct() \
    .count()
print("CHECK 2b:") 
print(n_distinct_snp_names == snp_map.shape[0])

#check we do not have duplicated SNP positions within each chromosome
from pyspark.sql.functions import countDistinct
n_distinct_snp_positions = df_samples_subset \
    .groupBy("Chr") \
    .agg( \
        countDistinct("Position")) \
    .select(F.sum("count(Position)")) \
    .collect()
    #group by chromosome, then count the unique cases in each chromosome and sum all the distinct cases across chromosomes to get the total number of distinct SNPs by position
        #https://sparkbyexamples.com/pyspark/pyspark-groupby-count-distinct/
        #https://sparkbyexamples.com/pyspark/pyspark-count-distinct-from-dataframe/
        #https://stackoverflow.com/questions/48568214/how-to-sum-the-values-of-a-column-in-pyspark-dataframe
#calculate again the number of unique positions within each chromosome using pandas this time in the SNP map
unique_pos_snp_map = sum([len(snp_map.loc[snp_map["Chromosome"]==chromo, "Position"].unique()) for chromo in np.unique(snp_map["Chromosome"])])
    #for each chromosome, get all the positions of that chromosome in the snp map, select the distinct cases and calculate the number of these distinct cases
print("CHECK 2c:")
diff_n_snps_distinct_position = snp_map.shape[0] - np.array(n_distinct_snp_positions)[0][0]
if diff_n_snps_distinct_position != 0:
    print("IMPORTANT: WE DO HAVE DUPLICATED SNPS CONSIDERING POSITION WITHIN CHROMOSOMES")
    print(f"The difference between the total number of SNPs in the map and the distinct SNP positions is {diff_n_snps_distinct_position}")
    print("The number of distinct positions per chromosome is the same in spark and in pandas?")
    print(unique_pos_snp_map == np.array(n_distinct_snp_positions)[0][0])
    print("SEE BELOW for filtering of duplicated snps considering POSITION and ALLELES with plink")

#check that each sample has the same number of rows, i.e., SNPs
df_samples_subset \
    .groupBy("Sample Index") \
    .count() \
    .toDF("Sample Index", "count") \
    .withColumn("check_2d", F.col("count") == snp_map.shape[0]) \
    .select("check_2d") \
    .distinct() \
    .show()
    #group rows per sample, count the number of rows per sample, convert to DF, then create a new column checking whether count is equal to the number of snps in the map, select that column and see if only true

#do more checks
#define function to do the checks. this will work on each final report as a pandas DF
#pdf = pd.read_csv(temp_dir.name+ "/" + zip_name + "_FinalReport1.txt", delimiter="\t", header=0, low_memory=False)
def calc_checks(pdf):

    #reorder rows based on SNP name in both pdf and the snp map to avoid errors
    pdf = pdf.sort_values(by="SNP Name")
    snp_map_ordered = snp_map.sort_values(by="Name")

    #create a pandas DF with the column names of the checks
    df = pd.DataFrame(columns=["index", "check_3", "check_4", "check_5", "check_6", "check_7", "check_8", "check_9"])

    #fill the different columns
    df["index"] = pdf["Sample Index"]
        #get the index of the sample
    df["check_3"] = len(pdf["SNP Index"].unique()) == snp_map_ordered.shape[0]
        #the number of unique SNP indexes is exactly the number of SNPs in the map?
    df["check_4"] = pdf["SNP Index"].values == snp_map_ordered["Index"].values
    df["check_5"] = pdf["SNP Name"].values == snp_map_ordered["Name"].values
    df["check_6"] = pdf["Chr"].values == snp_map_ordered["Chromosome"].values
    df["check_7"] = pdf["Position"].values == snp_map_ordered["Position"].values
        #each row (i.e., genotype) has exactly the same data regarding SNP, position... than int eh snp map?
    df["check_8"] = np.isin(pdf["Sample ID"].values[0], sample_map["ID"].values)
    df["check_9"] = np.isin(pdf["Sample Index"].values[0], sample_map["Index"].values)
        #the sample ID (all rows of a report belongs to the same sample) is included in the sample map?
    return df
#define schema of resulting DF
schema_checks = StructType() \
    .add("index",IntegerType(), nullable=True) \
    .add("check_3", BooleanType(), nullable=True) \
    .add("check_4", BooleanType(), nullable=True) \
    .add("check_5", BooleanType(), nullable=True) \
    .add("check_6", BooleanType(), nullable=True) \
    .add("check_7", BooleanType(), nullable=True) \
    .add("check_8", BooleanType(), nullable=True) \
    .add("check_9", BooleanType(), nullable=True)
#define a function to calculate conditional count only when condition is satisfied
cnt_cond = lambda cond: F.sum(F.when(cond, 1).otherwise(0))
#run the checks
checks_raw = df_samples_subset \
    .groupby("Sample Index") \
    .applyInPandas(calc_checks, schema_checks) \
    .agg(
        cnt_cond(F.col("check_3") == True).alias('check_3_cnt'), 
        cnt_cond(F.col("check_4") == True).alias('check_4_cnt'), 
        cnt_cond(F.col("check_5") == True).alias('check_5_cnt'), 
        cnt_cond(F.col("check_6") == True).alias('check_6_cnt'), 
        cnt_cond(F.col("check_7") == True).alias('check_7_cnt'), 
        cnt_cond(F.col("check_8") == True).alias('check_8_cnt'), 
        cnt_cond(F.col("check_9") == True).alias('check_9_cnt')) \
    .collect()
        #group by sample, and within each sample calculate the checks and generate a DF. This DF has a row per genotype, so all genotypes of the same sample are repeated, because of this we aggreagte
        #count only those cases that are true for each check, so we get just 1 row and 5 columns with the total number of Trues
            #https://stackoverflow.com/questions/49021972/pyspark-count-rows-on-condition

#The total number of trues in each checks should be equal to the number of genotypes, although we did operations within sample
for index, check in enumerate(checks_raw[0]):
    print("CHECK " + str(index+3) + ": " + str(check == n_genotypes))

#check that the distinct sample IDs are the same than the IDs in the sample map
distinct_sample_ids = df_samples_subset \
    .select(F.col("Sample Index"), F.col("Sample ID")) \
    .distinct() \
    .collect()
    #get the distinct cases of sample ID and Index, which are the same because each sample ID is associated with a sample index
sample_index_check = [row["Sample Index"] for row in distinct_sample_ids]
    #get each sample index
from natsort import index_natsorted, order_by_index
index_order = index_natsorted(sample_index_check)
    #get the indexes under natural order as done by natsorted
distinct_sample_ids_ordered = order_by_index(distinct_sample_ids, index_order)
    #use these indexes to order the list of with the sample IDs, so we get first Sample ID of Sample 1, then sample 2...
    #I have found difficulties to load data in spark by Sample ID order. I can use orderBy by Sample ID and SNP Index, but this increases the computation time a bit for a few samples, so not sure if with hundreds of samples will be much worse. It seems that groupBy is called every time we do an operation with the new DF. So I am just ordering for this check.
ordered_sample_ids = [row["Sample ID"] for row in distinct_sample_ids_ordered]
    #extract the Sample IDs once we have the rows ordered
print("CHECK 10")
#if n_samples is None, we can just compare with the whole list of samples in the sample map. If n_sample is a number, then we need to compare with the corresponding number of selected samples
if (n_samples == None) or (n_samples==216 and batch_name=="ILGSA24-17303") or (n_samples==1248 and batch_name=="ILGSA24-17873"):
    print(all(np.equal(ordered_sample_ids, sample_map["ID"].values)))
else:
    print(all(np.equal(ordered_sample_ids, sample_map["ID"].values[0:n_samples])))

#see unique alleles
print("\n#####################\n#####################")
print("see unique alleles")
print("#####################\n#####################")
df_samples_subset \
    .select(df_samples_subset["Allele1 - Forward"]) \
    .distinct() \
    .show()
df_samples_subset \
    .select(df_samples_subset["Allele2 - Forward"]) \
    .distinct() \
    .show()



############################
#### Plink Installation ####
############################

#I have downloaded Plink ([link](https://www.cog-genomics.org/plink/)) and copied the executable ("plink") to `bin`, so we can use Plink from any place just typing `plink`. We are using Plink version 1.9 (see singularity recipe for the version used).
print("\n#####################\n#####################")
print("see plink version")
print("#####################\n#####################")
run_bash("plink --version")

#Note that there is Plink 1.9 ([link](https://www.cog-genomics.org/plink/1.9/)) and Plink 2.0 ([link](https://www.cog-genomics.org/plink/2.0/)), these are not connected but different programs. 
    #This [threat](https://www.biostars.org/p/299855/#:~:text=The%20main%20difference%20is%20that,for%20a%20while%20to%20come.) of biostars explains the differences:
        # "The main difference is that plink 1.9 is essentially finished, while plink 2.0 is an alpha-stage program which will have significant unfinished components for a while to come. As a consequence, current development priorities for plink 2.0 are centered around things which are impossible to do with plink 1.9, such as handling multiallelic/phased variants and dosage data and reliably tracking REF/ALT alleles; while things that plink 1.9 already handles perfectly well, such as working with .ped/.map file pairs, have been deliberately deprioritized for now.
        # So, **you should stick to 1.9 as long as it's good enough for the jobs you need to perform. But once you need to do something outside 1.9's scope, you're likely to find that 2.0 already has the additional feature you need** (or it'll be added quickly after you ask for it)"

#see the dev version of this script for details about the final report like allele strands, plink usage, etc....



####################################################
#### Conversion of final report to Plink inputs ####
####################################################

# We are going to convert the illumina report (BPM file) to lgen file that can be used as input in Plink. This seems to be trivial ([link](https://www.biostars.org/p/51928/)), but we are going to use a tutorial just in case ([link](https://www.youtube.com/watch?v=5_pNby7l2dE&t=1s&ab_channel=GenomicsBootCamp), [link](https://pastebin.com/pzgw7JVp)).

# The final report has 1 line for the SNP of one individual. We have separated final reports for each individual. This is close to the lgen format of plink.

# Therefore, our goal is to create lgen, also the fam and map files required to load the lgen in Plink. Some information might be missing in the final report, so you need to replace them.


#############
# lgen file #
#############

print("\n#####################\n#####################")
print("Preparing lgen file")
print("#####################\n#####################")

# **lgen file [plink info](https://www.cog-genomics.org/plink/1.9/formats#lgen)**
    # A text file with no header line, and one line per genotype call (or just not-homozygous-major calls if 'lgen-ref' was invoked) usually with the following five fields:
        #Family ID
        #Within-family ID
        #Variant identifier
        #Allele call 1 ('0' for missing)
        #Allele call 2

#subset all_reports selecting some columns. The non-selected columns were used for checks
lgen_file = df_samples_subset.select(["Sample ID", "SNP Name", "Allele1 - Forward", "Allele2 - Forward"])
lgen_file.show()

# Change to "0" those genotype calls with "--" to match the format of Plink
#allele 1
lgen_file = lgen_file \
    .withColumn("Allele1 - Forward_final", \
        F.when( \
            (lgen_file["Allele1 - Forward"] == "-") | (lgen_file["Allele1 - Forward"] == "--"), "0") \
            .otherwise(lgen_file["Allele1 - Forward"])) \
    .withColumn("Allele2 - Forward_final", \
        F.when( \
            (lgen_file["Allele2 - Forward"] == "-") | (lgen_file["Allele2 - Forward"] == "--"), "0") \
            .otherwise(lgen_file["Allele2 - Forward"])) \
    .drop("Allele1 - Forward") \
    .drop("Allele2 - Forward")
        #https://stackoverflow.com/questions/44773758/how-to-conditionally-replace-value-in-a-column-based-on-evaluation-of-expression
        #https://stackoverflow.com/questions/29600673/how-to-delete-columns-in-pyspark-dataframe
lgen_file.show()

#check
print("\n#####################\n#####################")
print("Check that all SNPs do NOT have '-' or '--' for allele 1 and 2")
print("#####################\n#####################")
(lgen_file \
    .filter( \
        (lgen_file["Allele1 - Forward_final"].isin(["-", "--"])) | (lgen_file["Allele2 - Forward_final"].isin(["-", "--"]))) \
    .count()) == 0
        #https://stackoverflow.com/questions/63330350/pyspark-dataframe-filter-column-contains-multiple-value

#Add additional columns that are required for lgen files
lgen_file = lgen_file \
    .withColumn("FID", F.lit("combat_" + batch_name))
    #add a column using combat+batch name as literal value

# Reorder the columns
lgen_file = lgen_file \
    .select(["FID", "Sample ID", "SNP Name", "Allele1 - Forward_final", "Allele2 - Forward_final"])
        #https://stackoverflow.com/questions/42912156/python-pyspark-data-frame-rearrange-columns

#change the column of sample ID
lgen_file = lgen_file \
    .withColumnRenamed("Sample ID", "sample_id")
        #this column name will be used when writing the data per sample, so it is better to avoid the spaces

#duplicate sample_id column because when you partition by that column during writing, it gets removed
lgen_file = lgen_file \
    .withColumn("sample", F.col("sample_id"))
    #https://stackoverflow.com/questions/44575911/spark-how-to-prevent-dataframewriter-from-dropping-the-partitioning-columns-on

#look
print("\n#####################\n#####################")
print("see lgen file before writing. Note that the second 'sample' column will be removed after using it to save with partitionBy ")
print("#####################\n#####################")
lgen_file.show()

#save the data

#we can save data in binary plink format split across many files, this can be then merged, so we can just save the results of spark in different csv file
    #book plink, tutorial chirstofer
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3

#prepare folders to save the lgen files
run_bash("rm -rf data/genetic_data/plink_inputs/" + batch_name + "; mkdir -p data/genetic_data/plink_inputs/" + batch_name)

#you can save each final report separately by ID
lgen_file \
    .write \
        .partitionBy("sample") \
        .option("header", False) \
        .option("delimiter", "\t") \
        .option("compression", "gzip") \
        .csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/")
        #https://stackoverflow.com/questions/47780397/saving-dataframe-records-in-a-tab-delimited-file
        #https://stackoverflow.com/questions/40163996/how-to-save-a-dataframe-as-compressed-gzipped-csv
        #https://sparkbyexamples.com/pyspark/pyspark-partitionby-example/

#NOTE about cores and partitions
    #Note that the number of files generated depends on the size of each file and the number of cores. So if you have 5 cores, i.e., local[5], at least you will get 5 partions of the data.
    #you can check the number of partitions with "lgen_file.rdd.getNumPartitions()".
    #I have been doing analyses with the partitions that my local/cluster let me, i.e., 5 in my laptop and tens in the HPC.
    #Then, when writing, I specify that I want the data of each sample separated in different folders with "partitionBy". I do not care if several files are saved in the same individual folder if all belong to the same Sample ID, later i will merge the data
        #https://sparkbyexamples.com/pyspark/pyspark-partitionby-example/
        #https://stackoverflow.com/questions/48143159/spark-write-to-disk-with-n-files-less-than-n-partitions

#or a single file you can use PARTITION OR COALESCEN
#lgen_file \
#    .coalesce(1) \
#    .write \
#        .option("header", True) \
#        .option("delimiter", "\t") \
#        .option("compression", "gzip") \
#        .csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/")
    #https://sparkbyexamples.com/spark/spark-write-dataframe-single-csv-file/
    #https://sparkbyexamples.com/spark/spark-repartition-vs-coalesce/#dataframe-%20coalesce

#remove the duplicated sample column
print("\n#####################\n#####################")
print("see lgen file after removing the duplicated sample column")
print("#####################\n#####################")
lgen_file = lgen_file \
    .drop("sample")
lgen_file.show()



############
# fam file #
############

print("\n#####################\n#####################")
print("Preparing fam file")
print("#####################\n#####################")

# **Fam file ([PLINK sample information file](https://www.cog-genomics.org/plink/1.9/formats#fam))**
    # A text file with no header line, and one line per sample with the following six fields:
        # - Family ID ('FID')
        # - Within-family ID ('IID'; cannot be '0')
        # - Within-family ID of father ('0' if father isn't in dataset)
        # - Within-family ID of mother ('0' if mother isn't in dataset)
        # - Sex code ('1' = male, '2' = female, '0' = unknown)
        # - Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

##check sex variable
#before creating the fam file, we need to see what is going on with the sex, because there are differences between phenotype dataset and the sample map of the first batch

print("\n#####################\n#####################")
print("Check sex in sample map and phenotype data")
print("#####################\n#####################")

#load pheno data, this include reported sex and VO2 max data. I have checked that the data is the same directly reading from excel than converting to csv
pheno_data = pd.read_excel(
    "data/pheno_data/combact gene DNA GWAS 23062022.xlsx",
    header=0,
    sheet_name="All DNA samples")
print(pheno_data)

#there are several phenotypes, see dev version for details about possible errors.
#here we are only going to modify an entry that is clearly wrong and do some checks with the sample map. Also, we are going to use the sex indicated in pheno_data because there are some samples with differences between illumina estimated sex and the one showed in the pheno_data

#beep test 8 has an error, "o" letter instead "0" number in one sample
print("\n#####################\n#####################")
print("Week 8 beep test for a sample is 11.1O, i.e., letter O instead number 0")
print("#####################\n#####################")
index_problematic_sample = np.where(pheno_data["Week 8 beep test"] == "11.1O")[0][0]
index_problematic_column = np.where(pheno_data.columns == "Week 8 beep test")[0][0]
print(pheno_data.iloc[index_problematic_sample, index_problematic_column])
print(pheno_data.iloc[index_problematic_sample,:])

#change 11.1O for 11.10
print("\n#####################\n#####################")
print("error solved")
print("#####################\n#####################")
pheno_data.iloc[index_problematic_sample, index_problematic_column] = 11.1
print(pheno_data.iloc[index_problematic_sample,:])

#check
print("\n#####################\n#####################")
print("All samples with genetic data are included in pheno_data?")
print("#####################\n#####################")
map_samples_in_pheno = sum(sample_map["ID"].isin(pheno_data["AGRF code"]))
    #count the number of sample IDs from the map that are in the phenotype data
if map_samples_in_pheno != sample_map.shape[0]:
    print("IMPORTANT, WE HAVE SAMPLES WITH GENETIC DATA NOT INCLUDED IN PHENOTYPE DATA")
    print(f'We have {sample_map.shape[0] - map_samples_in_pheno} samples out {sample_map.shape[0]} with genetic data but not included in phenotype data')
    #samples with genetic but no phenotipic data and viceversa will be removed when perform association analyses.
    #I think it can be useful to have these individual without pheno data even if they are not included in the analyses because it can help us to define better the existence of substructure. It is like we analyze only samples with adiposity in HELENA, but we use the PCAs considering all individuals genotyped.

#get sample IDs and gender from sample map, then modify using pheno_data
fam_file = sample_map.loc[:, ["ID", "Gender"]]

#codify the sex variable following plink notation and considering the sex in pheno_data
fam_file.loc[
    (fam_file["ID"].isin(pheno_data.loc[pheno_data["Gender"].isna(), "AGRF code"])) | 
    (~fam_file["ID"].isin(pheno_data["AGRF code"])), 
    "Gender"] = "0"
    #set as 0 the sex of those samples whose
        #ID is associated to NA gender in pheno_data
        #OR
        #ID is NOT included in pheno_data
fam_file.loc[
    fam_file["ID"].isin(pheno_data.loc[pheno_data["Gender"] == "M", "AGRF code"]), 
    "Gender"] = "1"
    #set as 1 the sex of those samples whose ID is associated with Gender=M in pheno_data
fam_file.loc[
    fam_file["ID"].isin(pheno_data.loc[pheno_data["Gender"] == "F", "AGRF code"]), 
    "Gender"] = "2"

#check
print("\n#####################\n#####################")
print("The new sex variable is correctly coded?")
print("#####################\n#####################")
print((len(fam_file["Gender"].unique()) == 3) & (np.isin(fam_file["Gender"].unique(), test_elements=["1", "2", "0"])).all())
    #only three possible sexes, being 1, 2 and 0.
print(all(fam_file.loc[fam_file["Gender"] == "1", "ID"].isin(pheno_data.loc[pheno_data["Gender"] == "M", "AGRF code"])))
print(all(fam_file.loc[fam_file["Gender"] == "2", "ID"].isin(pheno_data.loc[pheno_data["Gender"] == "F", "AGRF code"])))
print(all(~fam_file.loc[fam_file["Gender"] == "0", "ID"].isin(pheno_data.loc[pheno_data["Gender"].isin(["F", "M"]), "AGRF code"])))
    #In the new variable,
        #1 includes only IDs that are "M" in pheno data
        #2 includes only IDs that are "F" in pheno data
        #0 NOT includes IDs that are "M" or "F" in pheno data

# Add the family variables and the phenotype (not added for now)
fam_file["FID"] = "combat_" + batch_name #ID for the whole study
fam_file["IID_father"] = "0" #'0' if father isn't in dataset
fam_file["IID_mother"] = "0" 
fam_file["phenotype"] = -9 #this is no data for phenotype

# Reorder
fam_file = fam_file[["FID", "ID", "IID_father", "IID_mother", "Gender", "phenotype"]]

print("\n#####################\n#####################")
print("See fam file:")
print("#####################\n#####################")
print(fam_file)

#save without header:
fam_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + ".fam.gz",
    sep="\t",
    header=None,
    compression='gzip',
    index=False)



############
# map file #
############

print("\n#####################\n#####################")
print("Preparing map file")
print("#####################\n#####################")

# **Map file [Plink info](https://www.cog-genomics.org/plink/1.9/formats#map)**

# A text file with no header line, and one line per variant with the following 3-4 fields:
    # - Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
    # - Variant identifier
    # - Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
    # - Base-pair coordinate

# We can get this information from the SNP map.
# Make a copy with copy(), so modifications to the data or indices of the copy will not be reflected in the original object
map_file = snp_map.copy()
print(map_file)

#set genetic position as dummy
map_file["centimorgans"] = 0 #dummy value of zero

#select columns
map_file = map_file.loc[:, ["Chromosome", "Name", "centimorgans", "Position"]]
print(map_file)

# Save
map_file.to_csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + ".map.gz",
    sep="\t",
    header=None,
    compression='gzip',
    index=False)



##########################################
#### create plink binaries per sample ####
##########################################

print("\n#####################\n#####################")
print("create plink binaries per sample")
print("#####################\n#####################")

#list the of folders with the lgen data per sample
import glob
lgen_full_paths = glob.glob("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/sample=*", recursive=False)
#check
print("\n#####################\n#####################")
print("We have the correct number of lgen folders according to fam_file and --n_samples?")
print("#####################\n#####################")
print(len(lgen_full_paths) == fam_file.shape[0])
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(len(lgen_full_paths) == 216)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(len(lgen_full_paths) == 1248)
else:
    print(len(lgen_full_paths) == n_samples)

#get the extensions of the files in these folders
list_files = glob.glob("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/**/*.csv.gz", recursive=True)
list_files = [path.split("/")[-1] for path in list_files]
list_extensions = ["-".join(path.split("-")[2:]) for path in list_files]
    #load paths of all compressed partitions of the lgen files for all samples
    #split each one by "/" and select the last one to get the file name
    #split each one by "-" and select all except the first two, then join with "-" so we get only the extension avoiding the partition name
    
#then check all files have the same extension. They should even if the belong to different samples
print("\n#####################\n#####################")
print("lgen files have the correct extension?")
print("#####################\n#####################")
if all([extension == list_extensions[0] for extension in list_extensions]):

    #perfect
    print("TRUE")

    #extract the extension without file type
    lgen_extension = list_extensions[0]
    lgen_extension = lgen_extension.split(".csv.gz")[0]
else:
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM WITH THE EXTENSION OF THE LGEN FILES")

#natural sorting, 1,10, 21.. that works with numbers + strings like 1b, 5c...
from natsort import natsorted
lgen_full_paths = natsorted(lgen_full_paths)
    #https://github.com/SethMMorton/natsort#installation

#create a folder to save the bed files for the selected batch
run_bash(
    "rm -rf data/genetic_data/plink_bed_files/" + batch_name + "; \
    mkdir -p data/genetic_data/plink_bed_files/" + batch_name)
    #mkdir
        #p flag: A flag which enables the command to create parent directories as necessary. If the directories exist, no error is specified

#define function to calculate the bed file for each lgen file
#lgen_full_path = lgen_full_paths[0]
def plink_inputs_prep(lgen_full_path):

    #extract the ID of the sample from the path of the folder
    sample_id_file = lgen_full_path.split("/")[-1]

    #decompress the different partitions of the lgen file but remove just in case before these files in case they are already decompressed
    run_bash(" \
        cd ./data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name +  "_lgen_files/" + sample_id_file + "; \
        gunzip -kf *" + lgen_extension + ".csv.gz")
        #-k: keep original compressed file
        #-f: force compression or decompression even if the file has multiple links or the corresponding file already exists

    #get the path for the files decompressed
    files_selected_sample = glob.glob("./data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name +  "_lgen_files/" + sample_id_file + "/*" + lgen_extension + ".csv", recursive=False)
    files_selected_sample = natsorted(files_selected_sample)
        #to have the same order of snps than in the snp map when we concatenate all these files
    
    #load all the files as DFs into a list
    list_dfs = [
        pd.read_csv(
            path, 
            header=None,
            delimiter='\t',
            low_memory=False) #to avoid problems with mixed types
        for path in files_selected_sample]
    
    #check
    print("\n#####################\n#####################")
    print(f" We loaded all the partitions of {sample_id_file}?: {len(list_dfs) == len(files_selected_sample)}")
    print("#####################\n#####################")

    # Concatenate, because we have one row per SNP 
    lgen_file_selected_sample = pd.concat(
        objs=list_dfs, #list DFs
        axis=0, #concatenate along rows
        ignore_index=True) #clear the existing index and reset it
    
    #check
    print("\n#####################\n#####################")
    print(f' Do we have as many rows as the total sum of rows across the list of DFs for {sample_id_file}? {lgen_file_selected_sample.shape[0] == sum([element.shape[0] for element in list_dfs])}')
    print("#####################\n#####################")

    #get the unique cases for the sample ID
    selected_sample_id = np.unique(lgen_file_selected_sample[1])

    #get the unique cases for the batch name
    selected_sample_batch = np.unique(lgen_file_selected_sample[0])[0].split("_")[1]

    #check whether
        #the batch name within the file correspond with the correct batch name
        #the number of unique Sample IDs is equal to 1 AND
        #the ID in the file is the same than the ID in the name of the folder
        #the number of unique SNP names is not equal than the number of SNPs in the snp map OR
        #the SNP names are not exactly the same than in the snp map
    sample_checks = \
        (selected_sample_batch == batch_name) and \
        (len(selected_sample_id) == 1) and \
        (selected_sample_id == sample_id_file.split("=")[1]) and \
        (len(np.unique(lgen_file_selected_sample[2])) == snp_map.shape[0]) and \
        all(lgen_file_selected_sample[2] == snp_map["Name"])
    if sample_checks: #if all True
        print("WE PASS SAMPLE CHECKS FOR LGEN FILE OF SAMPLE " + selected_sample_id)
    else:
        raise ValueError("ERROR! FALSE! WE HAVE IMPORTANT ERRORS IN THE LGEN FILE OF SAMPLE")

    #save 
    lgen_file_selected_sample.to_csv("./data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/" + sample_id_file + "/" + sample_id_file + ".lgen", 
        sep='\t',
        header=False, 
        index=False)

    #get the row of the fam file corresponding to the selected sample
    #if we use the whole fam file, we would get a ped file with as many rows as sample, but as we are only using the final report of one sample, all the rest of rows would be zero.
    #we need a a fma file with one row, so we can generate a ped file with one row for the sample selected.
    fam_file_to_subset = pd.read_csv("data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + ".fam.gz", 
        delimiter="\t", 
        header=None, 
        low_memory=False)
    select_fam = fam_file_to_subset.loc[fam_file_to_subset[1] == selected_sample_id[0], :]
    
    #additional check about the number of samples in the lgen file
    if select_fam.shape[0] != 1:
        raise ValueError("ERROR! FALSE! WE DO NOT HAVE 1 SAMPLE IN THE FAM FILE")

    #save that as the fam file for this sample
    select_fam.to_csv("./data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/" + sample_id_file + "/" + sample_id_file + ".fam", 
        sep='\t',
        header=False, 
        index=False)

    #remove the lgen file
    del(lgen_file_selected_sample)
    del(list_dfs)

    #copy the whole map file, in this case, we have the same snps, i.e., the same map in all samples. Indeed, illumina gives one map for the whole batch. Save the file in the folder of the sample
    run_bash(" \
        cd data/genetic_data/plink_inputs/" + batch_name + "; \
        gunzip -c " + batch_name + ".map.gz > " + batch_name + "_lgen_files/" + sample_id_file + "/" + sample_id_file + ".map")
        #c flag: writes the output stream to stdout and then you can use > to redirect to another folder. This will leave the compressed file untouched.
            #https://superuser.com/questions/45650/how-do-you-gunzip-a-file-and-keep-the-gz-file

    #create a folder to save the bed file of the sample
    run_bash(" \
        rm -rf data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file + "; \
        mkdir -p data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file)

    #create ped files file using the lgen files
    print("\n#####################\n#####################")
    print("create ped files")
    print("#####################\n#####################")
    run_bash(" \
        cd ./data/genetic_data/; \
        plink \
            --lfile ./plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/" + sample_id_file + "/" + sample_id_file + " \
            --out ./plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file + "/" + batch_name + "_" + sample_id_file + " \
            --recode")
        #go to the folder with plink inputs
        #--lfile for loading the lgen file, which should be accompanied by a .fam and .map files having the same name except the extension.
            #https://www.cog-genomics.org/plink/1.9/formats
        #--recode creates a new text fileset, after applying sample/variant filters and other operations. By default, the fileset includes a .ped and a .map file, readable with --file.
            #https://www.cog-genomics.org/plink/1.9/data#recode
        #--out for the output name

    #the result is ped and map files 
    #.ped
        #Original standard text format for sample pedigree information and genotype calls. Normally must be accompanied by a .map file; Haploview requires an accompanying .info file instead. Loaded with --file, and produced by --recode.
        #Contains no header line, and one line per sample with 2V+6 fields (columns) where V is the number of variants (as many columns as variants). The first six fields are the same as those in a .fam file. The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.
        #If all alleles are single-character, PLINK 1.9 will correctly parse the more compact "compound genotype" variant of this format, where each genotype call is represented as a single two-character string. This does not require the use of an additional loading flag. You can produce such a file with "--recode compound-genotypes".    
    #map.
        #Variant information file accompanying a .ped text pedigree + genotype table. Also generated by "--recode rlist".
        #A text file with no header line, and one line per variant with the following 3-4 fields:
            #Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
            #Variant identifier
            #Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
            #Base-pair coordinate
        #All lines must have the same number of columns (so either no lines contain the morgans/centimorgans column, or all of them do)

    #compress the lgen and map file
    run_bash(" \
        cd ./data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/" + sample_id_file + "; \
        gzip -f " + sample_id_file + ".lgen; \
        gzip -f " + sample_id_file + ".map; \
        gzip -f " + sample_id_file + ".fam")
        #-f: force compression even a compressed file with the same name already exists

    #convert the ped file to bed
    print("\n#####################\n#####################")
    print("convert ped to bed")
    print("#####################\n#####################")
    run_bash(" \
        cd ./data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file + "; \
        plink \
            --file ./" + batch_name + "_" + sample_id_file + " \
            --out ./" + batch_name + "_" + sample_id_file + " \
            --make-bed")
        #--file let you load ped files
            #https://www.cog-genomics.org/plink/1.9/formats#ped
        #--make-bed creates bed file
            #https://www.cog-genomics.org/plink/1.9/formats#bed
            #https://zzz.bwh.harvard.edu/plink/binary.shtml
        #you could just do these two step in one loading the lgen file and using make-bed
            #https://icb.med.cornell.edu/wiki/index.php/Plink/howto
        #bed file is NOT the bed format of USCS!! It is a plink format more compacted and it seems is faster to work with. It is in hexadecimal, meaning that the genotypes are not stored as number but in just two letters. For example, in 0xdc, the two last letters are storing the genotype of several individuals. dc in binary is 11011100. This means that the first sample is is 00, i.e., homozygous for the first allele of this snp in the bim file, the second sample is mozygous for the second allele of this snp in the bim file, and so on... If you more samples, another hexa code is added until all samples are coded and then it moves to the next snp...
            #https://www.cog-genomics.org/plink/1.9/formats#bed
            #https://www.biostars.org/p/113166/
            #https://coolconversion.com/math/binary-octal-hexa-decimal/Convert_binary_number_11011100_in_decimal_

    #remove the files we are not interested in
    run_bash(" \
        cd ./data/genetic_data/plink_bed_files/" + batch_name + "/01_bed_per_sample/" + sample_id_file + "; \
        rm ./" + batch_name + "_" + sample_id_file + ".ped; \
        rm ./" + batch_name + "_" + sample_id_file + ".map")
        #the hh file says there is a genotype that is diplodi but should be haplo, like Y in male
        #it is created when doing other operations like calculate the frequency

    #remove the decompressed files of each data partition of the lgen file for the selected sample
    run_bash(" \
        cd ./data/genetic_data/plink_inputs/" + batch_name + "/" + batch_name + "_lgen_files/" + sample_id_file + "; \
        rm *" + lgen_extension + ".csv")

#open a pool
if n_cores==None:
    pool_processors = None
        #This uses all cores available
else:
    pool_processors = int(n_cores)
import multiprocessing as mp
pool = mp.Pool(processes=pool_processors)

#apply the function across the pool
results_pararllelize = [pool.apply_async(plink_inputs_prep, args=(lgen_file,)) for lgen_file in lgen_full_paths]

#close the pool
pool.close()

#wait for the completion of all scheduled jobs
pool.join() #without this, the script is finished without that all samples have been finished.
    #https://stackoverflow.com/questions/44896547/python-apply-async-doesnt-execute-function

#NOTE: 
    #when we have one sample per bed file, the bim file, which is a map file with also the alele options, will have partial information. for example, if the unique sample is homozigous for a SNP, e.g., AA, the map in that position will have A and 0, because there are no more alleles.
        #no problem, when a sample with the other allele is included, the bim file adds the other allele. I have checked several cases before and after mergin all bed/bim files into one file.



#######################################
#### merge binaries of all samples ####
#######################################

print("\n#####################\n#####################")
print("merge binaries of all samples")
print("#####################\n#####################")

#create a folder to save the binary files of all samples and then merge
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "; \
    rm -rf ./02_data_to_merge; \
    mkdir -p ./02_data_to_merge")

#copy the bed/bim/fam files of all samples to the same folder
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "; \
    find ./01_bed_per_sample/  | \
        egrep '" + batch_name + ".*(bed|bim|fam)' | \
        while read file_path; \
            do cp $file_path ./02_data_to_merge/; \
        done")
    #cd to the folder with bed files for the selected batch
    #find in a subfolder with bed files per sample
        #any file starting with batch name and ending with bed/bim/fam extensions
        #read the paths of these files
            #copy the file to the folder for merging
        #https://stackoverflow.com/questions/15617016/copy-all-files-with-a-certain-extension-from-all-subdirectories
        #https://unix.stackexchange.com/questions/189210/find-command-multiple-conditions

#check
print("\n#####################\n#####################")
print("we are going to merge the correct number of plink files")
print("#####################\n#####################")
if (n_samples==None) and (batch_name == "ILGSA24-17303"):
    print(len(glob.glob("data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/" + batch_name + "*", recursive=False)) == 216*3)
elif (n_samples==None) and (batch_name == "ILGSA24-17873"):
    print(len(glob.glob("data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/" + batch_name + "*", recursive=False)) == 1248*3)
else:
    print(len(glob.glob("data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/" + batch_name + "*", recursive=False)) == n_samples*3)
    #the total number of files in the folder for merging should be 3 times the number of samples


##error in this check, i.e., the number of bed files generated in the second batch.
if n_samples==1248 and batch_name=="ILGSA24-17873":
    #if we are working with the problematic batch and we have all the samples, then we can take a look to the problem

    print("\n#####################\n#####################")
    print("ERROR/FALSE in the number of bed files generated in the second batch. It is not 1248 but it is ok if it is 1242, because the 6 duplicated samples are left of due to their SNPs are the same than in SNP map but in a different order. It is ok because these samples are going out anyway")
    print("#####################\n#####################")

    #the number of bim, bed and fam files across the whole batch is 1242, i.e., 1242*3, when it should be 1248 as we have 1248 final reports. Therefore, we are lacking data from 6 final reports. The problem comes from 6 final reports that seems to be duplicated.
    print("Number of bed, bim and fam files generated")
    run_bash("ls data/genetic_data/plink_bed_files/ILGSA24-17873/02_data_to_merge/*bed | wc -l")
    run_bash("ls data/genetic_data/plink_bed_files/ILGSA24-17873/02_data_to_merge/*bim | wc -l")
    run_bash("ls data/genetic_data/plink_bed_files/ILGSA24-17873/02_data_to_merge/*fam | wc -l")
    
    #There are three samples (7800AGSO , 1100JHJM  and 1200JPJM) that have two different final reports, with the extension _1 and _2.
    print("see sample map and pheno data for problematic samples")
    print(sample_map.loc[sample_map["ID"].isin(["1100JHJM_1", "1100JHJM_2", "7800AGSO_1", "7800AGSO_2", "1200JPJM_1", "1200JPJM_2"]), :])
    print(pheno_data.loc[pheno_data["AGRF code"].isin(["1100JHJM", "7800AGSO", "1200JPJM"]), :])
        #index Name ID
        # 30   NaN  7800AGSO_1 
        # 59   NaN  7800AGSO_2 
        #409   NaN  1100JHJM_1 
        #413   NaN  1100JHJM_2 
        #576   NaN  1200JPJM_1 
        #585   NaN  1200JPJM_2
    
    #to see the final reports of these samples
        #you can them directly looking at the FinalReportX where X is the index.
        #[30, 59, 409, 413, 576, 585] but subtracting 1
        #zipinfos_subset=[zipinfos_subset[i] for i in [29, 58, 408, 412, 575, 584]]
    
    #1200JPJM and 1100JHJM are indeed duplicated in pheno data! For these two cases seems to exist an error, a duplication in the Id of the original database, so in each case, there are two individuals with the same AGRF code but different age, metrics... so it seems an error. Indeed, they have different genotypes in some snps, indicating that they are actually two different individuals. I think we cannot know if 1200JPJM_1 is referring to the first individual with 1200JPJM in the pheno_data or the other way around, so I would remove these two samples.
    print("see genotypes of the two samples duplicated in pheno data")
    print(pheno_data.loc[pheno_data["AGRF code"].duplicated(keep=False), :])
    df_samples_subset.filter(F.col("Sample ID") == "1100JHJM_1").show()
    df_samples_subset.filter(F.col("Sample ID") == "1100JHJM_2").show()
    df_samples_subset.filter(F.col("Sample ID") == "1200JPJM_1").show()
    df_samples_subset.filter(F.col("Sample ID") == "1200JPJM_2").show()
    
    #in the case of 7800AGSO, there is no duplication of AGRF code (see above)
    print("see genotypes of the sample with with two final reports but that is not duplicated in pheno data")
    df_samples_subset.filter(F.col("Sample ID") == "7800AGSO_1").show()
    df_samples_subset.filter(F.col("Sample ID") == "7800AGSO_2").show()
        #I cannot see differences in the 20 first snps, but we have to check the whole list of snps


    ##explore more 7800AGSO
    print("explore more the sample that is not duplicated in pheno data")
    #select the allele 1 and 2 of both 7800AGSO_1 and 7800AGSO_2
    allele_1_7800AGSO_1 = df_samples_subset.filter(F.col("Sample ID") == "7800AGSO_1").select("Allele1 - Forward").collect()
    allele_2_7800AGSO_1 = df_samples_subset.filter(F.col("Sample ID") == "7800AGSO_1").select("Allele2 - Forward").collect()
    allele_1_7800AGSO_2 = df_samples_subset.filter(F.col("Sample ID") == "7800AGSO_2").select("Allele1 - Forward").collect()
    allele_2_7800AGSO_2 = df_samples_subset.filter(F.col("Sample ID") == "7800AGSO_2").select("Allele2 - Forward").collect()
    
    #get the indexes of those SNPs for which 7800AGSO_1 and 7800AGSO_2 differ
    #as we are filtering on all rows of a sample, then we should have the same number of rows (SNPs) and order than in the snp map, so we can use this map to extract index
    snp_index_differ_allele_1 = [index+1 for index in range(0, snp_map.shape[0], 1) if allele_1_7800AGSO_1[index] != allele_1_7800AGSO_2[index]]
    snp_index_differ_allele_2 = [index+1 for index in range(0, snp_map.shape[0], 1) if allele_2_7800AGSO_1[index] != allele_2_7800AGSO_2[index]]
        #for each SNP, if the allele 1 or 2 is not the same between the two samples, get the index of the SNP
        #we have to sum 1 because we are going to select these SNPs using SNP Index column, which starts in 1 not in 0, in contrast with python indexing.
    
    #show the first ten of these snps for each sample
    df_samples_subset \
        .filter( \
            (F.col("Sample ID") == "7800AGSO_1") & \
            (F.col("SNP Index").isin(snp_index_differ_allele_1[0:10]))) \
            .show()
    df_samples_subset \
        .filter( \
            (F.col("Sample ID") == "7800AGSO_2") & \
            (F.col("SNP Index").isin(snp_index_differ_allele_1[0:10]))) \
        .show()
    df_samples_subset \
        .filter( \
            (F.col("Sample ID") == "7800AGSO_1") & \
            (F.col("SNP Index").isin(snp_index_differ_allele_2[0:10]))) \
            .show()
    df_samples_subset \
        .filter( \
            (F.col("Sample ID") == "7800AGSO_2") & \
            (F.col("SNP Index").isin(snp_index_differ_allele_2[0:10]))) \
        .show()
        #we can clearly see that the genotypes for these snps are different between 7800AGSO_1 and 7800AGSO_2. There is no duplication of this AGRF code, but still we have it two times in the final reports with different genotypes.
    
    #from the samples IDs that are NOT in pheno data
    sample_map.loc[~sample_map["ID"].isin(pheno_data["AGRF code"]),:]
        #we can see that we have the 1200JPJM_1/2 and 1100JHJM_1/2 that are actually in pheno data but with the same ID, and then we have 7800AGSO_1/2 that is not duplicated in pheno data. I think this is the missing sample in pheno data.
    
    #The total number of samples in the phenotype data is 1463. In contrast, we have in the two batches 1248+216=1464 samples. Therefore, there is a missing sample in the phenotype data that could be 7800AGSO_2, because that code is not duplicated in the pheno_data.


#get all bed files paths and names
bed_full_paths = glob.glob("data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/" + batch_name + "*.bed", recursive=False)
    #https://www.geeksforgeeks.org/how-to-use-glob-function-to-find-files-recursively-in-python/
bed_file_names = [path.split("/")[-1].split(".")[0] for path in bed_full_paths]

#save the names in a txt file
with open(r"./data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/list_to_merge.txt", "w") as fp:
    fp.write("\n".join(bed_file_names))
        #each name in a different line so we have to add "\n" to the name
        #https://pynative.com/python-write-list-to-file/

#create folder to save merged data
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "; \
    rm -rf 03_merged_data; \
    mkdir -p 03_merged_data")

#merge with plink
print("\n#####################\n#####################")
print("merge with plink")
print("#####################\n#####################")
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "/02_data_to_merge/; \
    plink \
        --merge-list ./list_to_merge.txt \
        --out ../03_merged_data/" + batch_name + "_merged_data")
    #run plink on 02_data_to_merge where all the bed files are located
    #then save the output in a folder outside this folder ("../")    
    #https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3



#############################
#### look for duplicates ####
#############################

print("\n#####################\n#####################")
print("look for duplicates")
print("#####################\n#####################")

#create general folder to do operations related to duplicates
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    rm -rf ./04_inspect_snp_dup; \
    mkdir -p ./04_inspect_snp_dup")

#there are duplicated positions in the original snp map, because of this we see a warning when merging. This is common in illumina data
    #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
    #https://www.biostars.org/p/281276/

#create a folder to list duplicates
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/; \
    rm -rf ./00_list_dup; \
    mkdir -p ./00_list_dup")

#list duplicates
run_bash(" \
    cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    plink \
        --bfile ./03_merged_data/" + batch_name + "_merged_data \
        --list-duplicate-vars suppress-first ids-only\
        --out  ./04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates")
        #list-duplicate-vars to list duplicates by POSITION
            #suppress-first prevents the first variant in each group from being reported (since, if you're removing duplicates, you probably want to keep one member of each group).
            #ids-only modifier removes the header and the position/allele columns, so the generated list can be used as input for --exclude
                #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
        #in plink 2 you can also look for duplicates in ID, but we already check SNP names duplicates with spark and it is seems is ok

#load the duplicate list
duplicate_cases = pd.read_csv(
    "./data/genetic_data/plink_bed_files/" + batch_name + "/04_inspect_snp_dup/00_list_dup/" + batch_name + "_duplicates.dupvar", 
    sep="\t", 
    header=None,
    low_memory=False)

#the file generated has in each row ALL snps that have the same position and alleles, except the first one, thus, the number of rows it is not the total number of duplicates, we need also to know the number of snps in each row.
n_duplicates_plink = sum([len(row.split(" ")) if " " in row else 1 for row in duplicate_cases[0]])
    #select each row in the first and only column of duplicate cases (only one columne due to "ids-only" flag in plink)
    #split the row by space and count the number of pieces, i.e., snps, if there are spaces in the row
    #if not, just count 1
    #then sum
    #https://stackoverflow.com/questions/4406389/if-else-in-a-list-comprehension

#see number of plink duplicates
print("\n#####################\n#####################")
print("see number of plink duplicates")
print("#####################\n#####################")
print(n_duplicates_plink)
print(f"This is a {(n_duplicates_plink/snp_map.shape[0])*100} percent respect to the total number of SNPs")

#check
print("\n#####################\n#####################")
print("the number of duplicates based only on POSITION should be equal or higher than the duplicates based on POSITION AND ALLELES, i.e., plink duplicates. Duplicates in the latter have the same position and alleles, thus they should be included in the former list. Of course, the former list could have more snps that share position but not alleles, thus fulfilling criteria of the former but not the latter list")
print("#####################\n#####################")
print(diff_n_snps_distinct_position >= n_duplicates_plink)
    #There is no correspondence between the duplicated snps calculated with numpy and those of plink because plink uses both position and alleles. So it only removes snps with the same position AND the same alleles, irrespectively of the strand. For example, two snps with the same position and having A/T and T/A for reference and derived, will be considered duplicated by plink, while two snps with the same position but different alleles not. This explains why plink detects less duplicates than our approach with numpy/spark. It also explains why, if I select distinct snps according to the concatenation of position and allele colums (e.g., 3243241_A/T), I get less duplicates than with plink, because plink considers AT and TA, not only AT. In summary, I cannot replicate the filtering done by plink and it is ok to have different number of duplicates.
        #https://www.cog-genomics.org/plink/1.9/data#list_duplicate_vars
    #Also note that we do not have duplicated SNP names, so we do not have to apply filters for that, which would have to be applied using plink2.
        #https://www.biostars.org/p/360918/


##duplicates we will be removed in the next steps of the pipeline
#duplicated positions should be merged or removed. In our case, we are talking about 1% of the snps, so it should not be a problem. But we will do this latter, after merging the batches, when we apply multiple QC filters. I guess it is cleaner to change nothing on the data, see any batch effect, modify and then do filters. I guess that having as many markers as possible could help to detect batch effect.s

print("\n#####################\n#####################")
print("duplicates we will be removed in the next steps of the pipeline")
print("#####################\n#####################")

#compress the bed/bim/fam files
run_bash(
    "cd ./data/genetic_data/plink_bed_files/" + batch_name + "/; \
    gzip ./03_merged_data/" + batch_name + "_merged_data.bed; \
    gzip ./03_merged_data/" + batch_name + "_merged_data.bim; \
    gzip ./03_merged_data/" + batch_name + "_merged_data.fam")



######################
#### finishing up ####
######################

#remove the temp dir
temp_dir.cleanup()
    #https://stackoverflow.com/questions/3223604/how-to-create-a-temporary-directory-and-get-its-path-file-name

#stop spark env
spark.stop()
