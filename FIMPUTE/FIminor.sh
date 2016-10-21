#!/bin/bash
###############################
ref=$1
outref=$2
finaloutfile=$3
outformat=$4
Pedigree=$5
###############################
echo " "
echo " "
echo "@***********************************************************************@"
echo "@           Genotype imputation of sporadic missing alleles             @"
echo "@                      using FImpute version 2.0                        @"
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
       4. Format of output file (either in 'plink' or 'genoytpes')
       5. Pedigree data (IID Sire Dam Sex) (optional) "
echo " "
echo " "
  exit
fi

###################################################################
# # Download PLINK from the web if not available
if [ ! -f plink2 ]; then
wget https://www.dropbox.com/s/e3igtqgwpwmd0di/plink2?dl=0
mv plink2\?dl\=0 plink2
chmod +x plink2
fi

#if [ ! -f plink ]; then
# echo "plink was not found in the current directory, thus it been download"
# echo " "
# wget http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip
# unzip plink-1.07-x86_64.zip
# cp plink-1.07-x86_64/plink .
# ./plink --noweb --silent --file plink-1.07-x86_64/test 
# rm -r plink-1.07-x86_64*
# if [ ! -f plink.log ]; then
#  wget http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-i686.zip
#  unzip plink-1.07-i686.zip
#  cp plink-1.07-i686/plink .
#  rm -r plink-1.07-i686*
#  echo " "
# fi
#rm plink.log
#echo " "
#fi

# Download FImpute from the web if not available
if [ ! -f FImpute ]; then
 echo "FImpute was not found in the current directory, thus it been download"
 echo " "
 wget http://www.aps.uoguelph.ca/~msargol/fimpute/FImpute_Linux.zip
 unzip FImpute_Linux.zip
 cp FImpute_Linux/FImpute .
 rm -r FImpute_Linux*
fi
###################################################################

#### Checking if input file is available
if [ ! -f ${ref}.bed ] || [ ! -f ${ref}.bim ] || [ ! -f ${ref}.fam ]; then
 echo "*** File " ${ref}.* " representing the genotype file of PLINK binary format were not found ***"
 echo " "
 echo " "
 exit
fi

if [ -d $finaloutfile ]; then
 echo "*** Directory "$finaloutfile" already exist, delete it and re-run the script ***"
 echo " "
 echo " "
 exit
fi
###

if [ $outformat = plink ]; then 
 echo " "
elif [ $outformat = genotypes ]; then
 echo ""
else
 echo "Specify the correct output format***"
 echo " either as 'plink' or 'genotypes' "
 echo " "
 echo " "
 exit
fi


################################################
# Create temporary folder for analysis
FOLDER=tmp$RANDOM
mkdir ${FOLDER}
cp FImpute plink2 ${FOLDER}/.
cd ${FOLDER}
#################################################


echo " "
echo "****    Data processing for imputation started     ****"
#Allelecode=$(awk '{print $6}' ../${ref}.bim | sort | uniq | awk '{if ($1==1) print "12"; else if ($1=="B") print "AB"; else if($1=="G" || $1=="T" || $1=="C") print "ACGT"}')
Allelecode=$(echo '12')
#####  REFERENCE  #######
if [ $Allelecode = 12 ]; then
 cat ../${ref}.bim | awk '{print $2,2}' > recodeallele.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${ref} --make-bed --out ref_upd
 ./plink2 --silent --cow --noweb --nonfounders --bfile ref_upd --recode A --recode-allele recodeallele.txt --out geno
rm recodeallele.txt

elif [ $Allelecode = AB ]; then
 cat ../${ref}.bim | awk '{print $2,"A","B",1,2}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --bfile ../${ref} --update-alleles alleleupdate.txt --make-bed --out ref_upd
 cat ref_upd.bim | awk '{print $2,2}' > recodeallele.txt
 ./plink2 --silent --cow --noweb --nonfounders --bfile ref_upd --recode A --recode-allele recodeallele.txt --out geno
rm alleleupdate.txt recodeallele.txt

elif [ $Allelecode = ACGT ]; then
 echo "ACGT format is not allowed -- Please use an AB coding or 12"
 echo "The 12 format can be obtained with PLINK using the --recode12"
 echo " "
 echo " "
 exit
fi

awk 'NR>1 {print $2,1}' geno.raw > IDs_sons.txt
awk 'NR>1' geno.raw | cut -d' ' -f7- | awk '{gsub(/NA/,5); print}' | 
awk 'BEGIN {FS=" ";OFS=""} {$1=$1; print}' | 
paste -d' ' IDs_sons.txt - > Fimpsons.geno
echo 'IID Chip Call.........' > header
cat header Fimpsons.geno > $outref.geno
cat ref_upd.bim | awk '{print $2,$1,$4,NR}' > tmp
echo 'SNP_ID Chr Pos Chip1' > chipheader
cat chipheader tmp > $outref.snpinfo

rm geno.* IDs_sons.txt Fimpsons.geno chipheader header tmp
rm ref_upd.*


if [ -z "${Pedigree}" ]; then
echo "title=*population based imputation of ${finaloutfile}*;
genotype_file=*${outref}.geno*;
snp_info_file=*${outref}.snpinfo*;
output_folder=*$finaloutfile*;
save_genotype;
save_hap_lib /diplotype;
njob=4;" > $finaloutfile.ctr
sed -i 's/*/"/g' $finaloutfile.ctr
else
thresmm=$(awk 'END {print 250/NR}' ../${ref}.bim)
thresm=$(awk 'END {print 25/NR}' ../${ref}.bim)
echo "title=*population based imputation of ${finaloutfile}*;
genotype_file=*${outref}.geno*;
snp_info_file=*${outref}.snpinfo*;
ped_file=*../${Pedigree}*;
output_folder=*${finaloutfile}*;
parentage_test /ert_mm=$thresmm /ert_m=$thresm /find_match_cnflt /find_match_mp /find_match_ugp /remove_conflict;
add_ungen /min_fsize=5 /output_min_fsize=5 /output_min_call_rate=0.95 /save_sep; 
save_genotype;
save_hap_lib /diplotype;
njob=4;" > $finaloutfile.ctr
sed -i 's/*/"/g' $finaloutfile.ctr
fi

echo " "
echo "data processing eneded, FImpute will start soon ............"
echo " "

#***** run FImpute *****#
./FImpute $finaloutfile.ctr

if [ ! -d $finaloutfile ]; then
 echo "*** Imputation unsuccessful errors were detected ***"
 echo " "
 echo " "
 exit
fi

if [ ! -f ${finaloutfile}/genotypes_imp.txt ]; then
 echo "***** Imputation unsuccessful errors were detected *****"
 if [ -f ${finaloutfile}/report.txt ]; then
 cat ${finaloutfile}/report.txt
 fi
 echo " "
 echo " "
 exit
fi


echo "**********************************************************"
echo "******           Imputation finished             *********"
echo "******                                           *********"
echo "**********************************************************"

if [ $outformat = plink ]; then
echo "**  Preparing imputed files into PLINK binary data format *"

nloci=$(awk 'END {print NR}' ../${ref}.bim)
nanim=$(awk 'END {print NR}' ../${ref}.fam)
if [ $nloci -gt 50000 ]; then
 echo "This process might take some time to complete depending on the "
 echo "1. number of markers and samples"
 echo "2. computer processor speed"
elif [ $nanim -gt 1000 ]; then
 echo "This process might take some time to complete depending on the "
 echo "1. number of markers and samples"
 echo "2. computer processor speed"
fi

#******** Extract the imputed data and make a PLINK file ***********#
cat ${finaloutfile}/genotypes_imp.txt | awk 'NR>1 {print $3}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' | 
awk '{for (i=1;i<=NF;i++) { if($i==0) $i="1 1"; else if($i==1) $i="1 2"; else if($i==2) $i="2 2"; else if($i==5) $i="0 0"}  print}' > geno
cat ${finaloutfile}/genotypes_imp.txt | awk 'NR>1 {print $1,$1,0,0,0,-9}' > ids.txt
paste -d' ' ids.txt geno > file.ped
cat $outref.snpinfo | awk 'NR>1 {print $2,$1,0,$3}' > file.map
rm ids.txt geno

if [ $Allelecode = 12 ]; then
 ./plink2 --silent --cow --nonfounders --file file --make-bed --out ${finaloutfile}_imp
 ./plink2 --silent --cow --nonfounders --bfile ${finaloutfile}_imp --make-bed --out ${finaloutfile}_imp
rm file.map file.ped
elif [ $Allelecode = AB ]; then
 cat ../${ref}.bim | awk '{print $2,1,2,"A","B"}' > alleleupdate.txt
 ./plink2 --silent --cow --nonfounders --file file --update-alleles alleleupdate.txt --make-bed --out ${finaloutfile}_imp
 ./plink2 --silent --cow --nonfounders --bfile ${finaloutfile}_imp --make-bed --out ${finaloutfile}_imp
rm alleleupdate.txt file.map file.ped
fi


#******** Extract the imputed ungenotyped animals and make a PLINK file ***********#
if [ -f ${finaloutfile}/genotypes_imp_chip0.txt ]; then
 cat ${finaloutfile}/genotypes_imp_chip0.txt | awk 'NR>1 {print $3}' | 
 awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' | 
 awk '{for (i=1;i<=NF;i++) { if($i==0) $i="1 1"; else if($i==1) $i="1 2"; else if($i==2) $i="2 2"; else if($i==5) $i="0 0"}  print}' > geno
 cat ${finaloutfile}/genotypes_imp_chip0.txt | awk 'NR>1 {print $1,$1,0,0,0,-9}' > ids.txt
 paste -d' ' ids.txt geno > file.ped
 cat $outref.snpinfo | awk 'NR>1 {print $2,$1,0,$3}' > file.map
 rm ids.txt geno
 if [ $Allelecode = 12 ]; then
  ./plink2 --silent --cow --nonfounders --file file --make-bed --out ${finaloutfile}_ungenoimp
  ./plink2 --silent --cow --nonfounders --bfile ${finaloutfile}_ungenoimp --make-bed --out ${finaloutfile}_ungenoimp
 rm file.map file.ped
 elif [ $Allelecode = AB ]; then
  cat ../${ref}.bim | awk '{print $2,1,2,"A","B"}' > alleleupdate.txt
  ./plink2 --silent --cow --nonfounders --file file --update-alleles alleleupdate.txt --make-bed --out ${finaloutfile}_ungenoimp
  ./plink2 --silent --cow --nonfounders --bfile ${finaloutfile}_ungenoimp --make-bed --out ${finaloutfile}_ungenoimp
 rm alleleupdate.txt file.map file.ped
fi
fi


elif [ $outformat = genotypes ]; then
echo "**  Preparing imputed files into genotype file format *"
echo ' '
nloci=$(awk 'END {print NR}' ../${ref}.bim)
nanim=$(awk 'END {print NR}' ../${ref}.fam)
if [ $nloci -gt 50000 ]; then
 echo "This process might take some time to complete depending on the "
 echo "1. number of markers and samples"
 echo "2. computer processor speed"
elif [ $nanim -gt 1000 ]; then
 echo "This process might take some time to complete depending on the "
 echo "1. number of markers and samples"
 echo "2. computer processor speed"
fi

#******** Extract the imputed data ***********#
cat ${finaloutfile}/genotypes_imp.txt | awk 'NR>1 {print $3}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' > geno
cat ${finaloutfile}/genotypes_imp.txt | awk 'NR>1 {print $1}' > ids.txt
paste -d' ' ids.txt geno > ${finaloutfile}_imp.genotype
cat $outref.snpinfo | awk 'NR>1 {print $2,$1,0,$3}' > ${finaloutfile}_imp.map
rm ids.txt geno

#******** Extract the imputed ungenotyped animals and make a PLINK file ***********#
if [ -f ${finaloutfile}/genotypes_imp_chip0.txt ]; then
 cat ${finaloutfile}/genotypes_imp_chip0.txt | awk 'NR>1 {print $3}' | 
 awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' > geno 
 cat ${finaloutfile}/genotypes_imp_chip0.txt | awk 'NR>1 {print $1}' > ids.txt
 paste -d' ' ids.txt geno > ${finaloutfile}_ungenoimp.genotype
 rm ids.txt geno
fi

fi
rm *.snpinfo *.ctr *.geno
rm plink2 FImpute

cp -r * ../.
cd ..
rm -r $FOLDER

echo " "
echo " Imputation finished - files per chromosome are stored in the folder
       interMS-summary${finaloutfile} "
echo " "
echo " Imputed genotypes for all chromsomes are merge and stored in the currect directory as "

if [ $outformat = plink ]; then
echo " 
       ${finaloutfile}_imp.bed 
       ${finaloutfile}_imp.bim 
       ${finaloutfile}_imp.fam "
elif [ $outformat = genotypes ]; then
echo " 
       ${finaloutfile}_imp.genotype 
       ${finaloutfile}_imp.map "
fi

if [ -f ${finaloutfile}/genotypes_imp_chip0.txt ]; then
echo " 
       ${finaloutfile}_ungenoimp"
fi

echo " "
echo " "
echo " "
echo "@***********************************************************************@"
echo "@                   Report bugs to: solobaon@yahoo.com                  @"
echo "@                                    "$(date)"      @"
echo "@***********************************************************************@"
