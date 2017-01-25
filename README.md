# imputation
Genotype Imputation bash script for BEAGLE v3, v4 & v4.1 and FImpute software.  
The bash script uses PLINK format data and PLINK software itself to undertake most of the task.

## imputation with FImpute
**Parameter file description for runing FIMPUTE**  

There are two bash scripts  
### A. _FIminor.sh_   
This is to run minor imputation on a (one) dataset with few markers missing for some few individuals

### B. _FImajor.sh_  
This is to run major imputation on two different SNP chips (Eg. Impute 50k to HD or 7k to 50k etc)

Get help by runing the following: (The parameters needed to run the script will be printed out)  
./FIminor.sh help  
./FImajor.sh help  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
### Running the FIminor.sh script to undertake MINOR imputation with FIMPUTE 
* 6 Arguments are needed to run FIminor.sh script  
Arguments  
 1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
 2. output name of Reference file  
 3. Output name of final file after imputation  
 4. Output genotype format (either plink or genotypes format) 
 5. Allelecode
 6. Pedigree information (Optional -- with Progeny, Sire, Dam, Sex)  
The final out is a plink binary file with its prefix as argument and _imp  

### Running FImajor.sh script to undertake MAJOR imputation with FIMPUTE
 * 9 Arguments are needed to run FImajor.sh script  

        Arguments  
         1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
         2. output name of Reference file  
         3. The file to be imputed (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
         4. output name of file to be imputed  
         5. Output name of final file after imputation  
         6. where R is located  
         7. Output genotype format (either plink or genotype format)
         8. Allelecode
         9. Pedigree information (is "Optional" -- if provided, it should have 4 columns -- Progeny, Sire, Dam, Sex)  
         
         The final out is a plink binary file with its prefix as argument and _imp  


        Runing the Example files
        *************************************************************************************
           Examples for minor imputation (FIminor.sh)  
            **without Pedigree information**  
              ./FIminor.sh Example/ex01_ref ref resultsREF plink 12 

            **with Pedigree information**  
              ./FIminor.sh Example/ex01_ref ref resultsREF plink 12 Example/ex01.dat  

            **OUTPUT file-names**  
              resultsREF_imp.bim, resultsREF_imp.bed & resultsREF_imp.fam  

          Examples for minor imputation (FImajor.sh)
            **without Pedigree information**
              ./FImajor.sh Example/ex01_ref REF Example/ex01_valldchip val resultsVAL /usr/bin/Rscript plink 12

           **with Pedigree information** 
              ./FImajor.sh Example/ex01_ref REF Example/ex01_valldchip val resultsVAL /usr/bin/Rscript plink 12 Example/ex01.dat

           **OUTPUT file-names**  
             imp7kto50k_imp.bim , imp7kto50k_imp.bed & imp7kto50k_imp.fam  


## imputation with BEAGLEv4

There are two bash scripts  
### A. _BGLminor.sh_  
This is to run minor imputation on a (one) dataset with few markers missing for some individuals  

### B. _BGLmajor.sh_  
This is to run major imputation on two different SNP chips. (Eg. Impute 50k to HD or 7k to 50k etc)  

Get help by runing the following: (The parameters needed to run the script will be printed out)  
./BGLminor.sh help  
./BGLmajor.sh help  

* Running BGLminor.sh script to undertake MINOR imputation with BEAGLE version 4  
5 Arguments are needed to run BGLminor.sh script  
* Arguments  
 1. Reference file (The file should be a PLINK binary file with alleles coded as 11-12-22 or AA-AB-BB)  
 2. output name of Reference file  
 3. Output name of final file after imputation  
 4. start of chromosome  
 5. End of chromosome  
The final output is a plink binary file with its prefix as argument and suffix as _imp.bed, _imp.bim and _imp.fam.  
We generate a folder called interMS-summary and this contains imputation per chromosomes and outher important output of BEAGLE

* Running ./BGLmajor.sh script to undertake MAJOR imputation with BEAGLE version 4  
7 Arguments are needed to run BGLmajor.sh scrip  
*Arguments  
 1. Reference file (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
 2. output name of Reference file  
 3. The file to be imputed (The file should be a PLINK binary file with alleles coded as 11, 12, 22)  
 4. output name of file to be imputed  
 5. Output name of final file after imputation  
 6. start of chromosome  
 7. End of chromosome  
The final output is a plink binary file with its prefix as argument and suffix as _imp.bed, _imp.bim and _imp.fam.  
We generate a folder called interMS-summary and this contains imputation per chromosomes and outher important output of BEAGLE


#### Runing the Example files

##### Examples for minor imputation (BGLminor.sh)
###### for one chromosomes
./BGLminor.sh Example/ex01_ref ex01ref resultsREF 1 1

###### for a range of chromosomes
./BGLminor.sh Example/ex01_ref ex01ref resultsREF 25 29

###### for all chromosomes
./BGLminor Example/ex01_ref ex01ref resultsREF 1 29

**OUTPUT file-names**  
 * resultsREF_imp.bim  
 * resultsREF_imp.bed  
 * resultsREF_imp.fam  
In addtion a folder with prefix interMS-summary will be generated and addiotnal results files stored

##### Examples for minor imputation (BGLmajor.sh)
###### for one chromosomes
./BGLmajor.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 25 25

###### for a range of chromosomes
./BGLmajor.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 1 2

###### for all chromosomes
./BGLmajor.sh Example/ex01_ref ex01ref Example/ex01_valldchip ex01val resultsldchip 1 29

**OUTPUT file-names**  
 * resultsldchip_imp.bim  
 * resultsldchip_imp.bed  
 * resultsldchip_imp.fam  
In addtion a folder with prefix interMS-summary will be generated and addiotnal results files stored

