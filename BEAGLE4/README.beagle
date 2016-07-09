#*****************************************************************************#
#                    Parameter file description for running                   #
#                                BEAGLE VERSION 4                             #
#*****************************************************************************#

# NOTE -- VERY IMPORTANT
There are two bash scripts
A.
*******  BGLminor.sh *******    
  This is to run minor imputation on a (one) dataset 
  with few markers missing for some individuals

B. 
******  BGLmajor.sh ******
  This is to run major imputation on two different SNP chips
  Eg. Impute 50k to HD or 7k to 50k etc

C.
To ask for how to run any of the script just type:
BGLminor.sh help or BGLmajor.sh
The parameters needed to run the script will be print out
** Look at how to run the example file check the file "run_example.txt"

***************************************************************************************
* running ./BGLminor.sh script to undertake MINOR imputation with BEAGLE version 4 *
5 Arguments are needed to run ./BGLminor.sh script

*** Arguments *****
1. Reference file (The file should be a PLINK binary file with alleles coded as 11-12-22 or AA-AB-BB)
2. output name of Reference file
3. Output name of final file after imputation
4. start of chromosome
5. End of chromosome

** The final output is a plink binary file with its prefix as argument and suffix as _imp.bed, _imp.bim and _imp.fam
** We generate a folder called interMS-summary... 
-- This contains imputation per chromosomes and outher important output of BEAGLE

#########################################################################
running the Example file
#########################################################################
#for one chromosomes
./BGLminor.sh Example/ex01_ref ex01ref resultsREF 1 1

#for a range of chromosomes
./BGLminor.sh Example/ex01_ref ex01ref resultsREF 25 29

#for all chromosomes
./BGLminor Example/ex01_ref ex01ref resultsREF 1 29

OUTPUT files will have these names -- resultsREF_imp.bim; resultsREF_imp.bed; resultsREF_imp.fam;
*************************************************************************************************

running ./BGLmajor.sh script to undertake MAJOR imputation with BEAGLE version 4 ***#
7 Arguments are needed to run ./BGLmajor.sh script

*** Arguments *****
1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)
2. output name of Reference file
3. The file to be imputed (The file should be a PLINK binary file with alleles coded as 11, 12, 22)
4. output name of file to be imputed
5. Output name of final file after imputation
6. start of chromosome
7. End of chromosome

*** The final output is a plink binary file with its prefix as argument and suffix as _imp.bed, _imp.bim and _imp.fam

##################################################################################
running the example file
##################################################################################
#for one chromosomes
./BGLmajor.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 25 25

#for a range of chromosomes
./BGLmajor.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 1 2

#for all chromosomes
./BGLmajor.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 1 29

OUTPUT files will have these names -- resultsldchip_imp.bim; resultsldchip_imp.bed; resultsldchip_imp.fam;
**********************************************************************************************************
