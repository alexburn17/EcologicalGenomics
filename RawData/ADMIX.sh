#!/bin/bash

# Run ADMIXTURE to determine the number of genetic clusters in the SNP data, 
# and the ancestry proportions of each individual

# Here's the utility of 'for loops'...

for K in {1..10}

do 

admixture -C 0.000001 --cv ./SSW_all_biallelic.MAF0.02.Miss1.0.recode.vcf.geno $K \
| tee log${K}.out

done

# After the for loop finishes, you can use 'grep' to grab the values of the CV from each separate log file and append them into a new summary text file.

grep CV log*.out >chooseK.txt



