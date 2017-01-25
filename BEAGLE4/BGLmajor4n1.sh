#!/bin/bash
###############################
ref=$1
outref=$2
val=$3
outval=$4
finaloutfile=$5
chrst=$6
chrend=$7
Allelecode=$8
###############################
echo " "
echo " "
echo "@***********************************************************************@"
echo "@                      Genotype imputation from                         @"
echo "@             lower density chip to Higher density SNP panel            @"
echo "@                      using BEAGLE version 4.0                         @"
echo "@-----------------------------------------------------------------------@"
echo "@                      bash script written by:                          @"
echo "@               Solomon A. Boison |soloboan@yahoo.com|                  @"
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
       3. The file to be imputed (should PLINK binary file with alleles coded as 11, 12, 22)
       4. output name of file to be imputed
       5. Output name of final file after imputation
       6. start of chromosome
       7. End of chromosome
       8. Allelecode either 12 or AB for data with 11/12/22 or AA/AB/BB "
echo " "
echo " "
 exit
fi

# Download beagle from the web if not available
if [ ! -f beagle4n1.jar ] || [ ! -f beagle2vcf.jar ] || [ ! -f vcf2beagle.jar ]  || [ ! -f vcf2gprobs.jar ] || [ ! -f gprobsmetrics.jar ] || [ ! -f conform-gt.jar ]; then
 #echo "Beagle files were not found in the current directory, Downloading them ................"
 echo " Please download beagle version 4.1 "
 exit
fi

###################################################################
# # Download PLINK from the web if not available
if [ ! -f plink2 ]; then
#wget https://www.dropbox.com/s/e3igtqgwpwmd0di/plink2?dl=0
#mv plink2\?dl\=0 plink2
#chmod +x plink2
 echo " Please download plink version 2 and save the binary file as 'plink2' "
 exit
fi

#### Checking if beagle 4 jar file is available
if [ ! -f beagle4n1.jar ]; then
 echo "  Beagle jar file is not available -- Place beagle4n1.jar in this folder  "
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
if [ ! -f conform-gt.jar ]; then
 echo "  Beagle jar file is not available -- Place conform-gt.jar in this folder  "
 echo " "
 echo " "
 exit
fi

#### Checking if input file is available
if [ ! -f ${ref}.bed ] || [ ! -f ${ref}.bim ] || [ ! -f ${ref}.fam ]; then
 echo "*** File " ${ref}.* " representing the REFERENCE genotype file of PLINK binary format were not found ***"
 echo " "
 echo " "
 exit
fi
if [ ! -f ${val}.bed ] || [ ! -f ${val}.bim ] || [ ! -f ${val}.fam ]; then
 echo "*** File " ${val}.* " representing the TESTING genotype file of PLINK binary format were not found ***"
 echo " "
 echo " "
 exit
fi
####
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


## checking the allele code
#Allelecode_ref=$(awk '{print $6}' ../${ref}.bim | sort | uniq | awk '{if ($1==1) print "12"; else if ($1=="B") print "AB"; else if($1=="G" || $1=="T" || $1=="C") print "ACGT"}')
#Allelecode_val=$(awk '{print $6}' ../${val}.bim | sort | uniq | awk '{if ($1==1) print "12"; else if ($1=="B") print "AB"; else if($1=="G" || $1=="T" || $1=="C") print "ACGT"}')
Allelecode_ref=$(echo $Allelecode)
Allelecode_val=$(echo $Allelecode)

##################################################################
echo 'Data Preparation started for ........ BEAGLE Version 4'
######### Reference
if [ $Allelecode_ref = 12 ]; then
 cat ../${ref}.bim | awk '{print $2,1,2,"A","G"}' > alleleupdate.txt
 ./plink2 --silent --cow --noweb --nonfounders --bfile ../${ref} --update-alleles alleleupdate.txt --make-bed --out ${outref}_allelicUpd
 rm alleleupdate.txt
elif [ $Allelecode_ref = AB ]; then
 cat ../${ref}.bim | awk '{print $2,"A","B","A","G"}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${ref} --update-alleles alleleupdate.txt --make-bed --out ${outref}_allelicUpd
 rm alleleupdate.txt
elif [ $Allelecode_ref = ACGT ]; then
 echo "ACGT format is not allowed -- Please use an AB coding or 12"
 echo "The 12 format can be obtained with PLINK using the --recode12"
 echo " "
 echo " "
 exit
fi

## To be imputed set
if [ $Allelecode_val = 12 ]; then
 cat ../${val}.bim | awk '{print $2,1,2,"A","G"}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${val} --update-alleles alleleupdate.txt --make-bed --out ${outval}_allelicUpd
 rm alleleupdate.txt
elif [ $Allelecode_val = AB ]; then
 cat ../${val}.bim | awk '{print $2,"A","B","A","G"}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${val} --update-alleles alleleupdate.txt --make-bed --out ${outval}_allelicUpd
 rm alleleupdate.txt
elif [ $Allelecode_val = ACGT ]; then
 echo "ACGT format is not allowed -- Please use an AB or 12 coding"
 echo "The 12 format can be obtained with PLINK using the --recode12"
 echo " "
 echo " "
 exit
fi

#### Creating BEAGLE files for testing set
for i in `seq ${chrst} ${chrend}`
do
######### Reference
./plink2 --silent --cow --nonfounders --bfile ${outref}_allelicUpd --chr $i --recode transpose --out tmpref
echo " "
echo " "
echo "Number of markers and samples for chromosome ...... "$i "Reference SNP panel"
#### generating the beagle 3 format
cat tmpref.tfam | awk '{print $2,$2}' | tr '\n' ' ' | awk '{print "I","sampleID",$0}' > headers
cat tmpref.tped | cut -d' ' -f2,5- | awk '{print "M",$0}' > gtypes
cat headers gtypes > ${outref}_chr$i.bgl
rm headers gtypes
nmarkers=$(cat tmpref.tped | awk 'END{print NR}')
niid=$(cat tmpref.tfam | awk 'END{print NR}')
echo 'Num of SNP Markers ... '$nmarkers
echo 'Num of samples     ... '$niid

######### to be Imputed
./plink2 --silent --cow --nonfounders --bfile ${outval}_allelicUpd --chr $i --recode transpose --out tmpval
echo " "
echo "Number of markers and samples for chromosome ...... "$i "lower density SNP panel"
#### generating the beagle 3 format
cat tmpval.tfam | awk '{print $2,$2}' | tr '\n' ' ' | awk '{print "I","sampleID",$0}' > headers
cat tmpval.tped | cut -d' ' -f2,5- | awk '{print "M",$0}' > gtypes
cat headers gtypes > ${outval}_chr$i.bgl
rm headers gtypes
nmarkers=$(cat tmpval.tped | awk 'END{print NR}')
niid=$(cat tmpval.tfam | awk 'END{print NR}')
echo 'Num of SNP Markers ... '$nmarkers
echo 'Num of samples     ... '$niid
echo " "

#### refereence snplist
cat ${outref}_allelicUpd.bim | 
awk '{if($5==0 && $6=="A") $5="G"; else if($5==0 && $6=="G") $5="A"; print $1,$2,$4,$5,$6}' > Afreqchr.txt
cat Afreqchr.txt | awk -v i=$i '{if ($1==i) print $2,$3,$4,$5}' > ${outref}chr$i.snplist
cat ${outval}_allelicUpd.bim | 
awk '{if($5==0 && $6=="A") $5="G"; else if($5==0 && $6=="G") $5="A"; print $1,$2,$4,$5,$6}' > Afreqchr.txt
cat Afreqchr.txt | awk -v i=$i '{if ($1==i) print $2,$3,$4,$5}' > ${outval}chr$i.snplist
rm Afreqchr.txt tmpval.* tmpref.*

###################################################################################

echo '************************************************************************'
echo '*********            beagle version 4.0 imputation              ********'
echo '*********            chromosome' $i '...started !!!             ********'
echo "*********                                                       ********"
echo "*********            Phasing Reference samples                  ********"
echo "************************************************************************"

nloci=$(awk 'END {print NR}' ${outref}chr$i.snplist)
#**** REF
 java -jar beagle2vcf.jar $i ${outref}chr$i.snplist ${outref}_chr$i.bgl 0 > ${outref}chr$i.vcf
 java -jar beagle4n1.jar gt=${outref}chr$i.vcf nthreads=10 niterations=20 window=200000 overlap=50000 ne=200 gprobs=true out=${outref}chr$i.phased
if [ ! -f ${outref}chr$i.phased.vcf.gz ]; then
 echo " Imputation with BEAGLE v4 aborted  "
 echo " Analysis now at on chromosome "$i
 echo " errrors were detectd, check the log file "
 echo " "
 echo " "
 exit
fi
 echo " "
 echo '************************************************************************'
 echo "*********         Imputing lower denisty SNP panel              ********"
 echo "*********            chromosome' $i '...started !!!             ********"
 echo "************************************************************************"
#**** VAL
 java -jar beagle2vcf.jar $i ${outval}chr$i.snplist ${outval}_chr$i.bgl 0 > ${outval}_chr$i.vcf
 java -jar conform-gt.jar ref=${outref}chr$i.phased.vcf.gz gt=${outval}_chr$i.vcf chrom=$i out=${outval}_chr$i.vcfmod
#**** Impute
 java -jar beagle4n1.jar ref=${outref}chr$i.phased.vcf.gz gt=${outval}_chr$i.vcfmod.vcf.gz nthreads=10 niterations=20 window=200000 overlap=50000 ne=200 gprobs=true out=${finaloutfile}chr$i.imputed
if [ ! -f ${finaloutfile}chr$i.imputed.vcf.gz ]; then
 echo " Imputation with BEAGLE v4 aborted  "
 echo " Analysis now at on chromosome "$i
 echo " errrors were detectd, check the log file "
 echo " "
 echo " "
 exit
fi

## summary outputs
zcat ${finaloutfile}chr$i.imputed.vcf.gz | java -jar vcf2beagle.jar 0 ${finaloutfile}chr$i
zcat ${finaloutfile}chr$i.imputed.vcf.gz | java -jar vcf2gprobs.jar > ${finaloutfile}chr$i.gprobs
cat ${finaloutfile}chr$i.gprobs | java -jar gprobsmetrics.jar > ${finaloutfile}chr$i.summary
echo " "
done
rm ${outref}_allelicUpd* 
rm ${outval}_allelicUpd*

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
### data is recoded into original format 12 or AB 
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
echo " "

cat ${finaloutfile}chr$i.markers | awk -v chr=$i '{print chr,$1,0,$2}' > map
if [ $Allelecode_val = 12 ];then
 cat map | awk '{print $2,"A","G",1,2}' > updallele.txt
 ./plink2 --silent --cow --nonfounders --allow-no-sex --tfile data --update-alleles updallele.txt --make-bed --out ${finaloutfile}chr${i}_imp
rm updallele.txt
elif [ $Allelecode_val = AB ];then
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
fi

########################################
rm list bed fam bim merglist.txt *.nosex *.log imp.*
cd ..
rm *.nosex plink2 *.log *.jar
cd ..
mv $FOLDER/* .
rm -r $FOLDER
#####################################


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
