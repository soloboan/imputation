# Readme file for running BEAGLE VERSION 4.1 
  note: that version 4 and 4.1 are different. Read from https://faculty.washington.edu/browning/beagle/beagle.html

## There are two bash scripts  
These script could be used for minor (sporadic missing, ungenotyped parents) and major imputation (from low density SNP panel to a high density panel)  


### A.  *BGLminor4n1.sh*  
  This is to run minor imputation on a (one) dataset  
  with few markers missing for some few individuals  

### B. *BGLmajor4n1.sh*   
  This is to run major imputation on two different SNP chips  
  Eg. Impute 50k to HD or 7k to 50k etc  

Parameter file specification could be obtained from the commnad line:  
just type: 
BGLminor4n1.sh help or BGLmajor4n1.sh help
The parameters needed to run the script will be print out  


## Running *./BGLminor4n1.sh* - to undertake MINOR imputation with BEAGLE version 4.1 *
  6 Arguments are needed to run ./BGLminor.sh script

      ******    Arguments    ******** 
        1. Reference file (The file should be a PLINK binary file with alleles coded as 11-12-22 or AA-AB-BB)
        2. output name of Reference file
        3. Output name of final file after imputation
        4. start of chromosome
        5. End of chromosome
        6. Allelecodes - if your data is 11/12/22 then use 12; if AA/AB/BB use AB.
        
        ### The final output is a plink binary file with its prefix as argument and suffix as _imp.bed, _imp.bim and _imp.fam
        ** We generate a folder called interMS-summary... 
        -- This contains imputation per chromosomes and outher important output of BEAGLE

        ***************************************************************
                      running the Example file
        ***************************************************************
        #for one chromosome only 
          ./BGLminor4n1.sh Example/ex01_ref ex01ref resultsREF 1 1 12

        #for a range of chromosomes
          ./BGLminor4n1.sh Example/ex01_ref ex01ref resultsREF 25 29 12

        #for all chromosomes
          ./BGLminor4n1.sh Example/ex01_ref ex01ref resultsREF 1 29 12

        OUTPUT files will have these names -- resultsREF_imp.bim; resultsREF_imp.bed; resultsREF_imp.fam;
        *****************************************************************************************************

## Running *./BGLmajor4n1.sh* script to undertake MAJOR imputation with BEAGLE version 4 ***#
  8 Arguments are needed to run ./BGLmajor.sh script

      ******    Arguments    ******** 
        1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)
        2. output name of Reference file
        3. The file to be imputed (The file should be a PLINK binary file with alleles coded as 11, 12, 22)
        4. output name of file to be imputed
        5. Output name of final file after imputation
        6. start of chromosome
        7. End of chromosome
        8. Allelecodes - if your data is 11/12/22 then use 12; if AA/AB/BB use AB.

        *** The final output is a plink binary file with its prefix as argument and suffix as _imp.bed, _imp.bim and _imp.fam
              We also generate a folder called interMS-summary... 
              -- This contains imputation per chromosomes and outher important output of BEAGLE
              
        ***************************************************************
                    running the example file
        ***************************************************************
          #for one chromosome only
            ./BGLmajor4n1.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 25 25 12

          #for a range of chromosomes
            ./BGLmajor4n1.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 1 2 12

          #for all chromosomes
            ./BGLmajor4n1.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 1 29 12

          OUTPUT files will have these names -- resultsldchip_imp.bim; resultsldchip_imp.bed; resultsldchip_imp.fam;
        ***************************************************************************************************************



