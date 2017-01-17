#!/bin/bash
###############################
ref=$1
outref=$2
finaloutfile=$3
chrst=$4
chrend=$5
Allelecode=$6
###############################
echo " "
echo " "
echo "@***********************************************************************@"
echo "@           Genotype imputation of sporadic missing alleles             @"
echo "@                      using BEAGLE version 4.0                         @"
echo "@-----------------------------------------------------------------------@"
echo "@                      bash script written by:                          @"
echo "@                Solomon A. Boison |soloboan@yahoo.com|                 @"
echo "@                                                   bash_version 2.0.0  @"
echo "@                                                        22/09/2015     @"
echo "@***********************************************************************@"
echo " "
echo " "

#### Parameter file information
if [ ${ref} = help ]; then
 echo " These are the parameters needed (the order is important)"
 echo " ==== see the README file for more details ===="
 echo "*** Arguments *****
       1. Reference file (should be a PLINK binary file with alleles coded as 11, 12, 22)
       2. output name of Reference file
       3. Output name of final file after imputation
       4. start of chromosome
       5. End of chromosome
       6. The allele codes in your data (12, AB) [ACGT is not accepted] "
 echo " "
 echo " "
 exit
fi

# Download beagle from the web if not available
if [ ! -f beagle4.jar ] || [ ! -f beagle2vcf.jar ] || [ ! -f vcf2beagle.jar ]  || [ ! -f vcf2gprobs.jar ] || [ ! -f gprobsmetrics.jar ]; then
 echo "Beagle files were not found in the current directory, Download it from https://faculty.washington.edu/browning/beagle/b4_0.html ................"
 echo " Please download beagle version 4 "
 echo " rename the jar file to beagle4.jar "
fi

###################################################################
# # Download PLINK from the web if not available
if [ ! -f plink2 ]; then
echo " Please download Plink version 2 (or 1.9) "
echo ' Please rename the plink file to plink2 '
fi

#### Checking if beagle 4 jar file is available
if [ ! -f beagle4.jar ]; then
 echo "  Beagle jar file is not available -- Place beagle4.jar in this folder  "
 echo " "
 echo " "
 exit
fi
if [ ! -f beagle2vcf.jar ]; then
 echo "  Beagle jar file is not available -- Place beagle2vcf.jar in this folder  "
 echo " "
 echo " "
 exit
fi
if [ ! -f vcf2beagle.jar ]; then
 echo "  Beagle jar file is not available -- Place vcf2beagle.jar in this folder  "
 echo " "
 echo " "
 exit
fi
if [ ! -f vcf2gprobs.jar ]; then
 echo "  Beagle jar file is not available -- Place vcf2gprobs.jar in this folder  "
 echo " "
 echo " "
 exit
fi
if [ ! -f gprobsmetrics.jar ]; then
 echo "  Beagle jar file is not available -- Place gprobsmetrics.jar in this folder  "
 echo " "
 echo " "
 exit
fi

if [ -d interMS-summary${finaloutfile} ]; then
 echo "  Output directory exist, please delete "
 echo " "
 echo " "
 exit
fi

################################################
# Create temporary folder for analysis
FOLDER=tmp$RANDOM
mkdir ${FOLDER}
cp *.jar plink2 ./${FOLDER}
cd ${FOLDER}
#################################################

#### Checking if input file is available
if [ ! -f ../${ref}.bed ] || [ ! -f ../${ref}.bim ] || [ ! -f ../${ref}.fam ]; then
 echo "*** File " ${ref}.* " representing the REFERENCE genotype file of PLINK binary format were not found ***"
 echo " "
 echo " "
 exit
fi
#####

## check the allele code
#Allelecode=$(awk '{print $6}' ../${ref}.bim | sort | uniq | awk '{if ($1==1) print "12"; else if ($1=="B") print "AB"; else if($1=="G" || $1=="T" || $1=="C") print "ACGT"}')
Allelecode=$(echo ${Allelecode})
#####################################################################
echo 'Data Preparation started for ........ BEAGLE Version 4'
if [ $Allelecode = 12 ]; then
 cat ../${ref}.bim | awk '{print $2,1,2,"A","G"}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${ref} --update-alleles alleleupdate.txt --make-bed --out ${outref}_allelicUpd
rm alleleupdate.txt
elif [ $Allelecode = AB ]; then
 cat ../${ref}.bim | awk '{print $2,"A","B","A","G"}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${ref} --update-alleles alleleupdate.txt --make-bed --out ${outref}_allelicUpd
rm alleleupdate.txt
elif [ $Allelecode = ACGT ]; then
 echo "ACGT format is not allowed -- Please use an AB coding or 12"
 echo "The 12 format can be obtained with PLINK using the --recode12"
 echo " "
 echo " "
 exit
fi

#### Creating BEAGLE files for testing set
for i in `seq $chrst $chrend`
do
./plink2 --silent --cow --nonfounders --bfile ${outref}_allelicUpd --chr $i --recode --transpose --out tmpref
echo " "
echo " "
echo "Number of markers and samples for chromosome ...... "$i
#### generating the beagle 3 format
cat tmpref.tfam | awk '{print $2,$2}' | tr '\n' ' ' | awk '{print "I","sampleID",$0}' > headers
cat tmpref.tped | cut -d' ' -f2,5- | awk '{print "M",$0}' > gtypes
cat headers gtypes > ${outref}_chr$i.bgl
rm headers gtypes
nmarkers=$(cat tmpref.tped | awk 'END{print NR}')
niid=$(cat tmpref.tfam | awk 'END{print NR}')
echo 'Num of SNP Markers ... '$nmarkers
echo 'Num of samples     ... '$niid
rm tmpref.*
#############################################
echo " "

#### refereence snplist
cat ${outref}_allelicUpd.bim | 
awk '{if($5==0 && $6=="A") $5="G"; else if($5==0 && $6=="G") $5="A"; print $1,$2,$4,$5,$6}' > Afreqchr.txt
cat Afreqchr.txt | awk -v i=$i '{if ($1==i) print $2,$3,$4,$5}' > ${outref}chr$i.snplist
rm Afreqchr.txt 

###################################################################################
echo '*******************************************************************'
echo '*******              beagle version 4.0 imputation          ********'
echo '*******              chromosome' $i '...started !!!         ********'
echo "*******                                                     ********"
echo "********************************************************************"

nloci=$(awk 'END {print NR}' ${outref}chr$i.snplist)
if [ $nloci -le 50000 ]; then
#**** REF
 java -jar beagle2vcf.jar $i ${outref}chr$i.snplist ${outref}_chr$i.bgl 0 > ${outref}chr$i.vcf
 java -jar beagle4.jar gt=${outref}chr$i.vcf nsamples=10 burnin-its=10 phase-its=30 impute-its=30 out=${finaloutfile}chr$i.phased
elif [ $nloci -gt 50000 ]; then
#**** REF
 java -jar beagle2vcf.jar $i ${outref}chr$i.snplist ${outref}_chr$i.bgl 0 > ${outref}chr$i.vcf
 java -jar beagle4.jar gt=${outref}chr$i.vcf window=$nloci nsamples=10 burnin-its=10 phase-its=30 impute-its=30 out=${finaloutfile}chr$i.phased
fi

if [ ! -f ${finaloutfile}chr$i.phased.vcf.gz ]; then
 echo " Imputation with BEAGLE v4 aborted  "
 echo " Analysis now at on chromosome "$i
 echo " errrors were detectd, check the log file "
 echo " "
 echo " "
 exit
fi

zcat ${finaloutfile}chr$i.phased.vcf.gz | java -jar vcf2beagle.jar 0 ${finaloutfile}chr$i
zcat ${finaloutfile}chr$i.phased.vcf.gz | java -jar vcf2gprobs.jar > ${finaloutfile}chr$i.gprobs
cat ${finaloutfile}chr$i.gprobs | java -jar gprobsmetrics.jar > ${finaloutfile}chr$i.summary
echo " "
echo " "
done
rm -r ${outref}_allelicUpd*

if [ ! -d interMS-summary${finaloutfile} ]; then
 mkdir interMS-summary${finaloutfile}
 mv *chr* interMS-summary${finaloutfile}/.
 cp plink2 interMS-summary${finaloutfile}/.
 cd interMS-summary${finaloutfile}
else
 mv *chr* interMS-summary${finaloutfile}/.
 cp plink2 interMS-summary${finaloutfile}/.
 cd interMS-summary${finaloutfile}
fi

################################################################################
### Finally, merge all chromosomes and recode in PLINK binary data 
### data is recoded into original format 1 and 2
for i in `seq $chrst $chrend` 
do
echo "data transformation -- chromosome "$i 
cat ${finaloutfile}chr$i.markers | awk -v chr=$i '{print chr,$1,0,$2}' > map
zcat ${finaloutfile}chr$i.bgl.gz | awk 'NR==1' | sed 's/ /\n/g' | 
awk 'NR>2' | awk 'NR%2==1' | awk '{print $1,$1,0,0,0,-9}' > data.tfam
zcat ${finaloutfile}chr$i.bgl.gz | awk 'NR>1' | cut -d' ' -f3- > gtypes
paste -d' ' map gtypes > data.tped
nmarkers=$(cat data.tped | awk 'END{print NR}')
niid=$(cat data.tfam | awk 'END{print NR}')
echo 'Num of SNP Markers ... '$nmarkers
echo 'Num of samples     ... '$niid

if [ $Allelecode = 12 ];then
 cat map | awk '{print $2,"A","G",1,2}' > updallele.txt
 ./plink2 --silent --cow --nonfounders --allow-no-sex --tfile data --update-alleles updallele.txt --make-bed --out ${finaloutfile}chr${i}_imp
rm updallele.txt
elif [ $Allelecode = AB ];then
 cat map | awk '{print $2,"A","G","A","B"}' > updallele.txt
 ./plink2 --silent --cow --nonfounders --allow-no-sex --tfile data --update-alleles updallele.txt --make-bed --out ${finaloutfile}chr${i}_imp
rm updallele.txt
fi
done
rm data.* map gtypes

###############################################
ls -1 ${finaloutfile}chr*_imp.bed > bed
ls -1 ${finaloutfile}chr*_imp.bim > bim
ls -1 ${finaloutfile}chr*_imp.fam > fam
paste -d' ' bed bim fam > list
awk 'NR>1' list > merglist.txt

if [ $chrend -eq $chrst ]; then
 ./plink2 --silent --cow --nonfounders --allow-no-sex --bfile ${finaloutfile}chr${chrst}_imp --make-bed --out ../${finaloutfile}_imp
elif [ ! $chrend -eq $chrst ]; then
 bed=$(awk 'NR<2 {print $1}' list)
 bim=$(awk 'NR<2 {print $2}' list)
 fam=$(awk 'NR<2 {print $3}' list)
 ./plink2 --silent --cow --nonfounders --allow-no-sex --bed $bed --bim $bim --fam $fam --merge-list merglist.txt --make-bed --out imp
 ./plink2 --silent --cow --nonfounders --allow-no-sex --bfile imp --make-bed --out ../${finaloutfile}_imp
 rm imp.*
fi

#####################################
rm -r list bed fam bim merglist.txt *.log *.nosex
cd ..
rm *.nosex *.log *.jar
cd ..
mv $FOLDER/* .
rm -r $FOLDER
########################################

echo " "
echo " Imputation finished - files per chromosome are stored in the folder
       interMS-summary${finaloutfile} "
echo " "
echo " Imputed genotypes for all chromsomes are merge and stored in the currect directory as "
echo " 
       ${finaloutfile}_imp.bed 
       ${finaloutfile}_imp.bim 
       ${finaloutfile}_imp.fam "
echo " "
echo " "
echo " "
echo "@***********************************************************************@"
echo "@                   Report bugs to: soloboan@yahoo.com                  @"
echo "@                                    "$(date)"      @"
echo "@***********************************************************************@"
