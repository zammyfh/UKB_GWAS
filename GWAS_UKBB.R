#### In R

setwd("/well/ukbbmcm/workDirectories/zammy")

## extracting individuals IDs to be excluded from analysis
phenotypes_ukbb <- read.table("/well/got2d/ukbbMcM/workDirectories/zammy/phenotypes/adipose_pheno.txt",header=T,sep="\t",stringsAsFactors=F)

exclude <- phenotypes_ukbb[phenotypes_ukbb$Ethnic.background != 1001 | is.na(phenotypes_ukbb$Ethnic.background), "Encoded.anonymised.participant.ID" ]

withdrawl <- scan("../../phenotypes/w916_20170726_Withdrawals.csv")

total <-c(exclude, withdrawl)

total <- unique(total)

write.table(cbind(total, total), file="/well/ukbbmcm/workDirectories/zammy/GWAS/exclude_individs.txt", quote=F, row.names=F, col.names=F)

rel <-read.table("../../phenotypes/ukb916_rel_s488374.dat", header = T)

rel <- rel[,"ID1"]

unrel <- phenotypes_ukbb$Encoded.anonymised.participant.ID[! phenotypes_ukbb$Encoded.anonymised.participant.ID %in% rel]

write.table(cbind(unrel, unrel, "A"), file="/well/ukbbmcm/workDirectories/zammy/GWAS/unrelated_cluster.txt", quote=F, row.names=F, col.names=F)



cd /well/ukbbmcm/workDirectories/zammy/GWAS


for i in {1..22} ; do

## removes the LD regions from bfiles
/well/got2d/zamfh/software/plink --bfile ../../../genotypes/directlyGenotyped/ukb_snp_chr$i"_"v2 --exclude ld_regions_rem.txt --range --make-bed --out bfiles/region_removed/ukb_chr$i"_"reg


## removes related people and non-europeans
/well/got2d/zamfh/bolt_practice/Kuang/software/plink2 --bfile bfiles/region_removed/ukb_chr$i"_"reg --remove ../herit/exclude_individs.txt --make-bed --out bfiles/region_related_removed/ukb_snp_chr$i"_"red_unrel


## generates file of snps that can be pruned out based on unrelated european individuals
/well/got2d/zamfh/bolt_practice/Kuang/software/plink2 --bfile bfiles/region_related_removed/ukb_snp_chr$i"_"red_unrel --indep-pairwise 50 5 0.5 --out bfiles/prune_snps/snps_1_chr$i


## prunes out snps identified in unrelated individuals and removes non-europeans
/well/got2d/zamfh/bolt_practice/Kuang/software/plink2 --bfile bfiles/region_removed/ukb_chr$i"_"reg --exclude bfiles/prune_snps/snps_1_chr$i.prune.out --remove exclude_noneur_individs.txt --make-bed --out bfiles/region_LD_removed/t_ukb_chr$i"_"reg_ld


done


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


