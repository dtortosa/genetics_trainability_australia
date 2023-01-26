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

#run the python program for each batch. The programa uses python3.9 as interpreter (indicated at the top of the script)
./01b_illumina_report_to_plink.py --batch_name="ILGSA24-17303" > 01b_illumina_report_to_plink_ILGSA24_17303.out 2>&1 &
sleep 1s
	#we need delay the next process in order to leave time just in case
	#https://stackoverflow.com/questions/49944364/how-to-run-two-commands-but-with-a-delay-on-the-second-command-without-stopping
./01b_illumina_report_to_plink.py --batch_name="ILGSA24-17873" > 01b_illumina_report_to_plink_ILGSA24_17873.out 2>&1 &