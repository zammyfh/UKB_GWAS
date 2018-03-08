setwd("/well/ukbbmcm/workDirectories/zammy/GWAS")

print("Hello")

library(data.table)

library(plyr)

library(qqman)

library(stringr)

## setting arguments for script
args <- commandArgs(TRUE)
trait <- args[1] ## summary stat file


### read in the info and chromo deets 
info <- fread(paste0('zcat output/hrc_only/BMI_Females.bgen.gz'), h=T, stringsAsFactors=F, select=c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO"))


### read in rc meta-analysed results   
metal <- fread(paste0("output/hrc_only/sex_combined/", trait, "1.tbl"), h=T, stringsAsFactors = F, fill=T)

names(metal)[1] <- "SNP"

## merge the meta-analysis results and chromo deets
full <- merge(metal[, c("SNP", "Effect", "StdErr", "P-value")], info, by="SNP")

## rename the columns for ldsc
full <- rename(full, c("Effect"="Beta", "StdErr"="SE", "P-value"="P_value", "ALLELE1"="A1", "ALLELE0"="A2"))

fullT <- full[full$INFO>=0.3 & full$A1FREQ>5e-3 & full$SE<10 & full$A1FREQ<0.995, ]

########### Generating files to run in LD score

## select out rows containing rs numbers using grepl
full_ld <- fullT[grepl("rs", fullT$SNP), ]

## fast write to a space seperated file
fwrite(full_ld, file=paste0("ldsc_reg/input/sex_combined/", trait, ".info.tbl"), sep=" ", quote=F, na="NA")

##### running LD score from within R

## setting sample size for relevant sex
n <- 444478

## creating command for ldsc to modify file for running in ldscore 
ldsc <- paste0("module load python
                /well/got2d/zamfh/software/LDSC/ldsc/munge_sumstats.py \ --sumstats ldsc_reg/input/sex_combined/", trait, ".info.tbl \ --out ldsc_reg/processed_input/sex_combined/", trait, " \ --N ", n, " \ --merge-alleles /well/got2d/zamfh/software/LDSC/ldsc/w_hm3.snplist")
system(ldsc)

## creating command to run ld score with the modifiled file generated from above command
ldsc <- paste0("module load python
               /well/got2d/zamfh/software/LDSC/ldsc/ldsc.py \ --h2 ldsc_reg/processed_input/sex_combined/", trait, ".sumstats.gz \ --ref-ld-chr /well/got2d/zamfh/software/LDSC/ldsc/eur_w_ld_chr/ \ --w-ld-chr /well/got2d/zamfh/software/LDSC/ldsc/eur_w_ld_chr/ \ --out ldsc_reg/output/sex_combined/", trait, "_ld_reg")

## this runs the ldsc command and saves the stdout into ldsc_output
ldsc_output <- system(ldsc, intern=T)


## delete the files used fo rldscore regression because otherwise will clog up the computer - just duplicates of each other
rm <- paste0("rm ldsc_reg/input/cohort_wide_sex_combined/", trait, ".info.tbl
              rm ldsc_reg/processed_input/cohort_wide_sex_combined/", trait, ".sumstats.gz")

system(rm)

## this selects the numeric values out from the ldsc_output for the ldscore intercept which is mix of characters and numerics
intercept <- as.numeric(str_extract_all(ldsc_output[29], "[0-9.]+")[[1]])

## correcting the P-values if the ldsc regression os significantly different to 1
if ((intercept[1] - intercept[2]*1.96) > 1) {
  
  print(TRUE)
  
  fullT$P_value <- fullT$P_value*intercept[1]
  
  write.table(intercept[1], file=paste0("ldsc_reg/intercept/", trait, ".int"), row.names = F, col.names = F)
}

########### Counting the number of 'significant' hits

top <- fullT[fullT$P_value<=5e-08, ]

top <- top[order(top$P_value), ]

snps_signif <- c()

while(nrow(top)>0){
  
  snp <- top$SNP[1]
  
  snps_signif <- c(snp, snps_signif)
  
  snp_bp <-top$BP[top$SNP==snp] 
  
  chr <- top$CHR[top$SNP==snp]
  
  range_min <- snp_bp-500000
  
  range_max <- snp_bp+500000
  
  window <- top[top$BP<range_max & top$BP>range_min, ]
  
  window_snp <- window$SNP[window$CHR==chr]
  
  top <- top[! top$SNP %in% window_snp, ]
  
}

##snps_by_trait <- c(length(snps_signif), snps_by_trait)

##names(snps_by_trait) <- c(paste(trait, sex, sep="."), names(snps_by_trait))

## Exporting the details of the sex specific significant hits
write.csv(fullT[fullT$SNP %in% snps_signif], file=paste0("GWAS_output/sex_combined/", trait, ".gwas.hits"), quote=F, row.names=F)

## checking to see if any of the SNPs are heterogenous by sex

hits <- metal[metal$SNP %in% snps_signif, ]

hits <- metal[metal$Het_PVal < 0.0005, ]

write.csv(hits, file=paste0("GWAS_output/sex_combined/", trait, "sex.heterogeneity"), quote=F, row.names = F)


## Preparing the data for trans-ethnic meta-analysis

full$SNP_safe <- full$SNP


######## Generating QC plots for the sex combined data

png(paste0("QC/QQ/", trait, "_QQ.png"), height=5, width=6, units="in", res=300)

qq(fullT$P_value, main = paste0("QQ plot of ", trait, " GWAS p-values"))

lmbx <- qchisq(median(fullT$P_value, na.rm=T),1,lower.tail=FALSE)/qchisq(0.50,1,lower.tail=FALSE)

legend('topleft',legend=bquote("Median"~lambda == .(round(lmbx,2))),bg='white',bty='n',cex=1)

dev.off()

CHR <- as.numeric(fullT$CHR)
BP <- as.numeric(fullT$BP)
P <- fullT$P_value
P[is.na(P)] <- 1

gwas <- cbind.data.frame(CHR,BP,P)


png(paste0("QC/manhat/", trait, "manhat.png"), height=6, width=10, units="in", res=500)

manhattan(gwas)

legend('topleft',legend=trait, bg='white',bty='n',cex=1)


dev.off()


#### end of R script 


