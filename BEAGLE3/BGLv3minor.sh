#!/bin/bash
###############################
ref=$1
outref=$2
finaloutfile=$3
chrst=$4
chrend=$5
###############################
echo " "
echo " "
echo "@***********************************************************************@"
echo "@           Genotype imputation of sporadic missing alleles             @"
echo "@                      using BEAGLE version 3.3                         @"
echo "@-----------------------------------------------------------------------@"
echo "@                      bash script written by:                          @"
echo "@              Solomon A. Boison |soloboan@yahoo.com|                   @"
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
       5. End of chromosome"
 echo " "
 echo " "
 exit
fi

################################################
# Download beagle if unavailable
if [ ! -f beagle3.jar ]; then
 echo "Beagle file were not found in the current directory, Downloading them ................"
 echo " "
 wget http://faculty.washington.edu/browning/beagle/beagle.jar
 mv beagle.jar beagle3.jar
fi

###################################################################
# # Download PLINK from the web if not available
if [ ! -f plink2 ]; then
wget https://www.dropbox.com/s/e3igtqgwpwmd0di/plink2?dl=0
mv plink2\?dl\=0 plink2
chmod +x plink2
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


#### Checking if beagle 3.3 jar file is available
if [ ! -f beagle3.jar ]; then
 echo "  Beagle jar file is not available -- Place beagle3.jar in this folder  "
 echo " "
 echo " "
 exit
fi

#### Checking if input file is available
if [ ! -f ../${ref}.bed ] || [ ! -f ../${ref}.bim ] || [ ! -f ../${ref}.fam ]; then
 echo "*** File " ${ref}.* " representing the REFERENCE genotype file of PLINK binary format were not found ***"
 echo " "
 echo " "
 exit
fi
#####

## checking allele codes
#Allelecode=$(awk '{print $6}' ../${ref}.bim | sort | uniq | awk '{if ($1==1) print "12"; else if ($1=="B") print "AB"; else if($1=="G" || $1=="T" || $1=="C") print "ACGT"}')
Allelecode=$(echo '12')
###################################################################################
echo 'Data Preparation started for ........ imputation with BEAGLE Version 3'
if [ $Allelecode = ACGT ]; then
 echo "ACGT format is not allowed -- Please use an AB coding or 12"
 echo "The 12 format can be obtained with PLINK using the --recode12"
 echo " "
 echo " "
 exit
fi

#### Creating BEAGLE files for testing set
for i in `seq ${chrst} ${chrend}`
do
./plink2 --silent --cow --bfile ../${ref} --chr $i --recode transpose --out tmpref
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
#############################################
echo " "

#### refereence snplist
if [ $Allelecode = 12 ]; then
 cat ../${ref}.bim | 
 awk '{if($5==0 && $6==1) $5=2; else if($5==0 && $6==2) $5=1;print $1,$2,$4,$5,$6}' > Afreqchr.txt
 cat Afreqchr.txt | awk -v i=$i '{if ($1==i) print $2,$3,$4,$5}' > ${outref}chr$i.snplist
 rm Afreqchr.txt
elif [ $Allelecode = AB ]; then
 cat ../${ref}.bim | 
 awk '{if($5==0 && $6=="A") $5="B"; else if($5==0 && $6=="B") $5="A";print $1,$2,$4,$5,$6}' > Afreqchr.txt
 cat Afreqchr.txt | awk -v i=$i '{if ($1==i) print $2,$3,$4,$5}' > ${outref}chr$i.snplist
 rm Afreqchr.txt
fi

###################################################################################
echo "********************************************************"
echo '******         beagle version 3.3 imputation      ******'
echo '******         chromosome' $i '...started !!        ******'
echo "******                                            ******"
echo "********************************************************"
java -jar beagle3.jar missing=0 unphased=${outref}_chr$i.bgl nsamples=10 niterations=30 out=${finaloutfile}
echo " "
echo " "
if [ ! -f ${finaloutfile}.${outref}_chr$i.bgl.phased.gz ]; then
 echo " Imputation with BEAGLE v3.3 aborted  "
 echo " Analysis now at chromosome "$i
 echo " errrors were detectd, check the log file "
 echo " "
 echo " "
 exit
fi
done

################################################################################
### Finally, merge all chromosomes and recode in PLINK binary data 
### data is recoded into original format 1 and 2
for i in `seq $chrst $chrend`
do
echo "data transformation -- chromosome "$i 
cat ${outref}chr$i.snplist | awk -v chr=$i '{print chr,$1,0,$2}' > map
zcat ${finaloutfile}.${outref}_chr$i.bgl.phased.gz | 
awk 'NR==1' | sed 's/ /\n/g' |
awk 'NR>2' | awk 'NR%2==1' | awk '{print $1,$1,0,0,0,-9}' > data.tfam
zcat ${finaloutfile}.${outref}_chr$i.bgl.phased.gz | 
awk 'NR>1' | cut -d' ' -f3- > gtypes
paste -d' ' map gtypes > data.tped
nmarkers=$(cat data.tped | awk 'END{print NR}')
niid=$(cat data.tfam | awk 'END{print NR}')
echo 'Num of SNP Markers ... '$nmarkers
echo 'Num of samples     ... '$niid

#transfrom to plink binary data
./plink2 --silent --cow --tfile data --make-bed --out ${finaloutfile}chr${i}_imp
rm data.tped data.tfam
done

###############################################
ls -1 ${finaloutfile}chr*_imp.bed > bed
ls -1 ${finaloutfile}chr*_imp.bim > bim
ls -1 ${finaloutfile}chr*_imp.fam > fam
paste -d' ' bed bim fam > list
awk 'NR>1' list > merglist.txt

if [ $chrend -eq $chrst ]; then
 ./plink2 --silent --cow --bfile ${finaloutfile}chr${chrst}_imp --make-bed --out ${finaloutfile}_imp
elif [ ! $chrend -eq $chrst ]; then
 bed=$(awk 'NR<2 {print $1}' list)
 bim=$(awk 'NR<2 {print $2}' list)
 fam=$(awk 'NR<2 {print $3}' list)
 ./plink2 --silent --cow --bed $bed --bim $bim --fam $fam --merge-list merglist.txt --make-bed --out ${finaloutfile}_imp
 ./plink2 --silent --cow --bfile ${finaloutfile}_imp --make-bed --out ${finaloutfile}_imp
fi

mkdir interMS-summary${finaloutfile}
mv *chr* interMS-summary${finaloutfile}/.

########################################
rm list tmp* map bed fam bim gtypes merglist.txt *.nosex *.log *.jar
rm interMS-summary${finaloutfile}/*.nosex
rm interMS-summary${finaloutfile}/*.log
cp -r * ../.
cd ..
rm -r ${FOLDER}
##############################

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
echo "@***********************************************************************@"
echo "@                   Report bugs to: soloboan@yahoo.com                  @"
echo "@                                     "$(date)"      @"
echo "@***********************************************************************@"
