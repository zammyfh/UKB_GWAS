#### In R

setwd("/well/ukbbmcm/workDirectories/zammy")

withdrawl <- scan("../../phenotypes/w916_20170726_Withdrawals.csv")

non_eur <- read.table("/well/ukbbmcm/workDirectories/anubha/share/exclude.nonEuropeans.list", header=T)

rel <-read.table("../../phenotypes/ukb916_rel_s488374.dat", header = T)

twins <- rel[rel$Kinship>0.4,]

plink.fam <- read.table("/well/ukbbmcm/genotypes/directlyGenotyped/ukb_snp_chr1_v2.fam", header=F)

bgen.fam <- read.table("/well/ukbbmcm/phenotypes/ukb916_imp_chr1_v2_s487406.sample", header = T)

missing.fam <- plink.fam[,1][! plink.fam$V1 %in% bgen.fam$ID_1]

bolt_exclude <- unique(c(non_eur$FID, withdrawl, twins$ID1, missing.fam))

plink_exclude <- unique(c(non_eur$FID, withdrawl, rel$ID1))

write.table(cbind(bolt_exclude, bolt_exclude), file="/well/ukbbmcm/workDirectories/zammy/GWAS/id_files/bolt_individs.txt", quote=F, row.names=F, col.names=F)

write.table(cbind(plink_exclude, plink_exclude), file="/well/ukbbmcm/workDirectories/zammy/GWAS/id_files/plink_individs.txt", quote=F, row.names=F, col.names=F)



cd /well/ukbbmcm/workDirectories/zammy/GWAS


for i in {1..22} ; do

## removes the LD regions and plink undesirables from bfiles
/well/got2d/zamfh/software/plink --bfile ../../../genotypes/directlyGenotyped/ukb_snp_chr$i"_"v2 --remove ../herit/exclude_individs.txt --exclude ld_regions_rem.txt --range --make-bed --out bfiles/to_prune/ukb_chr$i"_"to_prune


## generates file of snps that can be pruned out based on unrelated european individuals
/well/got2d/zamfh/bolt_practice/Kuang/software/plink2 --bfile bfiles/to_prune/ukb_chr$i"_"to_prune --indep-pairwise 50 5 0.5 --out bfiles/model_snps/snps_1_chr$i


done

## appending all snps wanted for --model-snps 
cat bfiles/model_snps/*.in >> bfiles/model_snps/model_snps.txt


## creating corrects fam file

awk '{print $1,$2,$3,$4,$5,$5}' /well/ukbbmcm/genotypes/directlyGenotyped/ukb_snp_chr1_v2.fam > /well/ukbbmcm/workDirectories/zammy/GWAS/id_files/ukb_snp_chr1_correct_col6.fam


## extract chromo num and position from HRC data

awk '{print $1 ":" $2}' /well/ukbbmcm/workDirectories/wei/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | tail -n +2 >> hrc/hrc_snps.txt

cd hrc

awk -F\: '{print>$1 ".txt"}' hrc_snps.txt

## qctool to keep only the hrc imputed snps in the data file

for i in {1..22} ; do

../extractor/qctool -g /well/ukbbmcm/genotypes/imputed/ukb_imp_chr$i"_"v2.bgen -og impfiles/ukb_imp_chr$i"_".bgen -incl-positions hrc/$i.txt

done


## generating exclude file

cat /well/ukbbmcm/genotypes/imputed/ukb_mfi_chr*_v2.txt | awk '{print $1}' >> impfiles/bgen.snps

awk '{print $3}' /well/ukbbmcm/workDirectories/wei/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | tail -n +2 >> hrc/hrc_snps_rs.txt

comm -23 <(sort impfiles/bgen.snps) <(sort hrc/hrc_snps_rs.txt) > impfiles/hg_imp.snps 

grep -F -x -v -f hrc/hrc_snps_rs.txt impfiles/bgen.snps > impfiles/hg_2.snps


setwd("/well/ukbbmcm/workDirectories/zammy")

library(GenABEL)

library(genetics)

##install.packages("plyr")
library(plyr)

library(readr)

library(stringr)


phenotypes_ukbb <- read.table("/well/got2d/ukbbMcM/workDirectories/zammy/phenotypes/adipose_pheno.txt",header=T,sep="\t",stringsAsFactors=F)

## giving better names to athro measures
phenotypes_ukbb <- rename(phenotypes_ukbb, c("Hip.circumference" = "HIP", "Waist.circumference" = "WC", "Body.mass.index..BMI." = "BMI", "Age.when.attended.assessment.centre" = "age", "Sex" = "sex"))

## generating age squared variable
phenotypes_ukbb$age2<-phenotypes_ukbb$age^2

## generating whr
phenotypes_ukbb$WHR <- phenotypes_ukbb$WC / phenotypes_ukbb$HIP

phenotypes_ukbb$UK.Biobank.assessment.centre <- as.factor(phenotypes_ukbb$UK.Biobank.assessment.centre)

## adding in genetic stuff

pcs <- read.table("scores_analysis/input_files/pcs.txt",header=T, sep="", stringsAsFactors=F)

pcs_ids <-read.table("scores_analysis/input_files/pc_id.txt", header=F, sep="", stringsAsFactors=F)

pcs_tot <- cbind.data.frame(pcs_ids, pcs)

array <- read.table("../../qc/ukb_sqc_v2.txt", header = T)

array <- array[,1]

array  <- cbind.data.frame(pcs_ids, array)

names(phenotypes_ukbb)[1] <-"pt_id"

names(pcs_tot)[1] <-"pt_id"

names(array)[1] <- "pt_id"

phenotypes_ukbb <- merge(phenotypes_ukbb, pcs_tot, by="pt_id")

phenotypes_ukbb <- merge(phenotypes_ukbb, array, by="pt_id")

## need to add in excluding whites and the withdrawl people - generated earlier
exclude <- read.table("/well/ukbbmcm/workDirectories/zammy/GWAS/id_files/bolt_individs.txt")

phenotypes_ukbb <- phenotypes_ukbb[! phenotypes_ukbb$pt_id %in% exclude[,1], ]

## transforming anth traits
trait_list<-list("BMI", "WHR", "WC", "HIP")

phenotypes_ukbb <- phenotypes_ukbb[order(phenotypes_ukbb$sex),]

dataset_list <- split(phenotypes_ukbb, list(phenotypes_ukbb$sex)) # c, paste

####couldn't work out how to make this work with a loop....
for (i in 1:length(dataset_list)){
  
  for (trait in trait_list){
    
    dataset_list[[i]][, paste0(trait, "_transformed")] <-NA
    
    res <- resid(lm(dataset_list[[i]][, trait] ~ dataset_list[[i]]$age + dataset_list[[i]]$age2 + dataset_list[[i]]$UK.Biobank.assessment.centre + ., data=dataset_list[[i]][, paste0("PC", 1:6)]))
    
    dataset_list[[i]][complete.cases(dataset_list[[i]][,trait]), paste0(trait, "_transformed")] <- qnorm((rank(res, na.last="keep")-0.5)/sum(!is.na(res)))
  }
}

phenotypes_ukbb<- data.frame(do.call(rbind, dataset_list))

sex <- c("Females", "Males")

names(dataset_list) <- sex

output <- lapply(sex, function(sex)  dataset_list[[sex]][ ,c(1, 1, 7, 18:28)])
  
output <- lapply(output, function(tab) rename(tab, c("pt_id" = "FID", "pt_id.1" = "IID")))

names(output) <- sex

lapply(sex, function(sex) write.table(output[[sex]], paste0("GWAS/phenotype/", sex, "_pheno_file.txt"), quote = F, row.names = F))



## GWAS

#########pooled.ebmd.impute.geno.sh#########
#!/bin/bash
#file name pooled.ebmd.impute.geno.sh
#$ -cwd -V -N pooled.bmd
#$ -q long.qc
#$ -P mccarthy.prjc
#$ -j y

cd /well/ukbbmcm/workDirectories/zammy/GWAS


/well/got2d/zamfh/software/BOLT-LMM_v2.3/bolt \
--fam=/well/ukbbmcm/workDirectories/zammy/GWAS/id_files/ukb_snp_chr1_correct_col6.fam \
--bed=/well/ukbbmcm/genotypes/directlyGenotyped/ukb_snp_chr{1:22}_v2.bed \
--bim=/well/ukbbmcm/genotypes/directlyGenotyped/ukb_snp_chr{1:22}_v2.bim \
--remove=/well/ukbbmcm/workDirectories/zammy/GWAS/id_files/bolt_individs.txt \
--numThreads=10 \
--LDscoresFile=/well/got2d/weigan/tools/BOLT-LMM_v2.3/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/well/got2d/weigan/tools/BOLT-LMM_v2.3/tables/genetic_map_hg19_withX.txt.gz \
--lmmForceNonInf \
--modelSnps=/well/ukbbmcm/workDirectories/zammy/GWAS/bfiles/model_snps/model_snps.txt \
--phenoFile=/well/ukbbmcm/workDirectories/zammy/GWAS/phenotype/Females_pheno_file.txt \
--bgenFile=/well/ukbbmcm/genotypes/imputed/ukb_imp_chr{1:22}_v2.bgen \
--sampleFile=/well/ukbbmcm/phenotypes/ukb916_imp_chr1_v2_s487406.sample \
--bgenMinMAF=1e-3 \
--bgenMinINFO=0.3 \
--covarFile=/well/ukbbmcm/workDirectories/zammy/GWAS/phenotype/Females_pheno_file.txt \
--covarCol=array \
--qCovarCol=PC{1:6} \
--lmm \
--phenoCol=BMI_transformed \
--statsFile=/well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_BMI_Female_stats.gz \
--statsFileBgenSnps=/well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_BMI_Female_bgen_stats.gz \
--exclude=/well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/hg_2.snps \
--verboseStats
2>&1 | tee /well/ukbbmcm/workDirectories/zammy/GWAS/output/BMI_female.log



!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N females_bmi_bolt
#$ -P mccarthy.prjc -q long.qc
#$ -o /gpfs2/well/ukbbmcm/workDirectories/zammy/GWAS/log.out
#$ -e /gpfs2/well/ukbbmcm/workDirectories/zammy/GWAS/log.err
#$ -j y
#$ -pe shmem 10

/well/ukbbmcm/workDirectories/zammy/GWAS/scripts/bolt_bmi_females.sh



## R to generate a merged bfile file for PCA for bolt

setwd("/well/ukbbmcm/workDirectories/zammy/GWAS")

files <- sapply(2:22, function(i) {
  
  c(paste0("bfiles/region_LD_removed/t_ukb_chr", i, "_reg_ld.bed"), paste0("bfiles/region_LD_removed/t_ukb_chr", i, "_reg_ld.bim"), paste0("bfiles/region_LD_removed/t_ukb_chr", i, "_reg_ld.fam"))})


files <- t(files)

write.table(files, file="all_pruned_files.txt", quote=F, row.names=F, col.names=F)


## Bash to generate PCA files using plink 

## plink to merge the bfiles into one master bfile
/well/got2d/zamfh/software/plink --bfile bfiles/region_LD_removed/t_ukb_chr1_reg_ld --merge-list all_pruned_files.txt --make-bed --out bfiles/region_LD_removed/all_pruned_ukb

## converting bfiles to ped for EIGENSOFT pca caluclation
/well/got2d/zamfh/software/plink --bfile bfiles/region_LD_removed/all_pruned_ukb --recode --out all_pruned_ukb


## plink to generate pcas for the pruned snps - doesn't work because data too big and would take too long
/well/got2d/zamfh/software/plink --bfile bfiles/region_LD_removed/all_pruned_ukb --within unrelated_cluster.txt --pca --pca-cluster-names A --out ukb_bolt_pcs


