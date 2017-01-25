# Readme file for running FImpute  

## There are two bash scripts  
These script could could be used for minor (sporadic missing, ungenotyped parents) and major imputation (from low density SNP panel to a high density panel)  

### A. FIminor.sh  
  This is to run minor imputation on a (one) dataset  
  with few markers missing for some few individuals  

### B. FImajor.sh  
  This is to run major imputation on two different SNP chips  
  Eg. Impute 50k to HD or 7k to 50k etc  

Parameter file specification could be obtained from the commnad line:  
just type:  
FIminor.sh help or FImajor.sh  
The parameters needed to run the script will be print out  

## Running ./FIminor.sh - This will undertake MINOR imputation  
6 Arguments are needed to run ./FIminor.sh script  

    ******    Arguments    ********  
      1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
      2. output name of Reference file  
      3. Output name of final file after imputation  
      4. Output genotype format (either plink or genotypes format)  
      5. Allelecodes - if your data is 11/12/22 then use 12; if AA/AB/BB use AB.  
      6. Pedigree information (Optional -- with Progeny, Sire, Dam, Sex [M,F])  

     ### The final out is a plink binary file with its prefix as argument and _imp   

      ************************************************************************************  
                          runing the Example file  
      ************************************************************************************  
        without Pedigree information  
          * ./FIminor.sh Example/ex01_ref ref resultsREF plink *  
        
        with Pedigree information  
          * ./FIminor.sh Example/ex01_ref ref resultsREF plink Example/ex01.dat *   
          
        OUTPUT file -- resultsREF_imp.bim; resultsREF_imp.bed; resultsREF_imp.fam;  
      *************************************************************************************  

## Running ./FImajor.sh script to undertake MAJOR imputation with FIMPUTE  
9 Arguments are needed to run ./FImajor.sh script  

    ******    Arguments    ********  
      1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
      2. output name of Reference file  
      3. The file to be imputed (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
      4. output name of file to be imputed  
      5. Output name of final file after imputation  
      6. where R is located  (eg. /usr/bin/Rscript)
      7. Output genotype format (either plink or genotype format)  
      8. Allelecodes - if your data is 11/12/22 then use 12; if AA/AB/BB use AB.  
      9. Pedigree information (is "Optional" -- if provided, it should have 4 columns -- Progeny, Sire, Dam, Sex [M,F])  
      
      ### The final out is a plink binary file with its prefix as argument and _imp  


        **********************************************************************************  
                              runing the example file  
        **********************************************************************************  
          without Pedigree information  
            *./FImajor.sh Example/ex01_ref REF Example/ex01_valldchip val resultsVAL /usr/bin/Rscript plink*  
            
          with Pedigree information   
            *./FImajor.sh Example/ex01_ref REF Example/ex01_valldchip val resultsVAL /usr/bin/Rscript plink Example/ex01.dat*  
            
          OUTPUT will have names -- imp7kto50k_imp.bim; imp7kto50k_imp.bed; imp7kto50k_imp.fam;  
        **************************************************************************************************
