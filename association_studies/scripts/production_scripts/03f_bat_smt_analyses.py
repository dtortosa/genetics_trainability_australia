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



###################################################################
######## ANALYZE ASSOCIATIONS WITHIN RELEVANT SET OF GENES ########
###################################################################

#We are going to analyze the associations within relevant set of genes like the BAT connectome.



#FOR BAT ANALYSES
    #this would an additional step in this project that would be outside of the paper
    #take the 1000kb gene windows for all coding genes, liftover to hg38. If the USCS tool accepts genomic ranges, just use them as input, if not, split in two datasets the start and the end of the gene windows
    #for each phenotype (VO2, beep....), calculate the average (better than median because want influence of outliers within gene like in iHS, if a SNPs is veery important in a gene that should influence the info about the whole gene) p-value for the association of SNPs inside each gene
    #then, calculate 1000 random sets of genes, within each set, calculate the median association of all genes inside the set and compare with the BAT set to obtain a distribution and empirical p-value (is association lower in BAT? LOOF BAT PAPER SCRIPTS FOR THIS). Here I want median because i do not want a gene outliser change things, I want the overall impact of BAT genes in general, not just a few genes.
