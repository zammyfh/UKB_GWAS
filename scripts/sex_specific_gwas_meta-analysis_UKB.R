#!/bin/sh
pheno=$1

echo "SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON" > scripts/temp/$pheno.metal.combine.sh

for sex in Females Males ; do

##unzipping BGEN files
gunzip output/hrc_only/$pheno"_"$sex.bgen.gz

##changing NA values
sed -i 's/\-nan/NA/g' output/hrc_only/$pheno"_"$sex.bgen

## rezip file
gzip output/hrc_only/$pheno"_"$sex.bgen

## start commands for metal put into the metal scripts
echo "MARKER SNP
ALLELE Allele1 Allele0
EFFECT BETA
FREQ A1FREQ
STDERR SE" >> scripts/temp/$pheno.metal.combine.sh

## -f gives true if the file specified exists
if [ -f ldsc_reg/intercept/$pheno.$sex.int ]; then

## this prints the contects of the intercept file into the meta-analysis script to tell metal the value for GC correction
awk '{print "GENOMICCONTROL", $1}' ldsc_reg/intercept/$pheno.$sex.int >> scripts/temp/$pheno.metal.combine.sh

## printing file names into metal scripts
echo "PROCESS "output/hrc_only/$pheno"_"$sex.bgen.gz >> scripts/temp/$pheno.metal.combine.sh 

else

## printing file names into metal scripts
echo "GENOMICCONTROL OFF
PROCESS "output/hrc_only/$pheno"_"$sex.bgen.gz >> scripts/temp/$pheno.metal.combine.sh 

fi

done

## adding further commands to metal shell
echo "OUTFILE output/hrc_only/sex_combined/"$pheno" .tbl" >> scripts/temp/$pheno.metal.combine.sh 

echo "ANALYZE HETEROGENEITY
QUIT" >> scripts/temp/$pheno.metal.combine.sh 


/well/got2d/zamfh/software/generic-metal/executables/metal < scripts/temp/$pheno.metal.combine.sh

