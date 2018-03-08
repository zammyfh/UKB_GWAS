setwd("/well/ukbbmcm/workDirectories/zammy/GWAS")

library(data.table)

library(plyr)

library(qqman)

library(stringr)

## setting arguments for script
args <- commandArgs(TRUE)
trait <- args[1] ## summary stat file
sex <- args[2]



## reads in the gzip files that have been produced by BOLT
full <- fread(paste0('zcat output/hrc_only/', trait, '_', sex, '.bgen.gz'), h=T, stringsAsFactors=F)


## rename the columns for ldsc
full <- rename(full, c("P_BOLT_LMM"="P_value", "A1FREQ"="A1_Freq", "ALLELE0"="ALLELE2"))

fullT <- full[full$INFO>=0.3 & full$A1_Freq>5e-3 & full$SE<10 & full$A1_Freq<0.955, ]


########### Generating files to run in LD score

## select out rows containing rs numbers using grepl
full_ld <- fullT[grepl("rs", fullT$SNP), ]

## fast write to a space seperated file
fwrite(full_ld, file=paste0("ldsc_reg/input/sex_specific/", trait, ".", sex, ".info.tbl"), sep=" ", quote=F, na="NA")

##### running LD score from within R

## providing info on the sample size for men and women
sample <- c(240745, 203733)

names(sample) <- c("Females", "Males")

## setting sample size for relevant sex
n <- sample[sex]

## creating command for ldsc to modify file for running in ldscore 
ldsc <- paste0("module load python
                /well/got2d/zamfh/software/LDSC/ldsc/munge_sumstats.py \ --sumstats ldsc_reg/input/sex_specific/", trait, ".", sex, ".info.tbl \ --out ldsc_reg/processed_input/sex_specific/", trait, ".", sex, " \ --N ", n, " \ --merge-alleles /well/got2d/zamfh/software/LDSC/ldsc/w_hm3.snplist")

system(ldsc)

## creating command to run ld score with the modifiled file generated from above command
ldsc <- paste0("module load python
               /well/got2d/zamfh/software/LDSC/ldsc/ldsc.py \ --h2 ldsc_reg/processed_input/sex_specific/", trait, ".", sex, ".sumstats.gz \ --ref-ld-chr /well/got2d/zamfh/software/LDSC/ldsc/eur_w_ld_chr/ \ --w-ld-chr /well/got2d/zamfh/software/LDSC/ldsc/eur_w_ld_chr/ \ --out ldsc_reg/output/sex_specific/", trait, ".", sex, "_ld_reg"
)

## this runs the ldsc command and saves the stdout into ldsc_output
ldsc_output <- system(ldsc, intern=T)

## this selects the numeric values out from the ldsc_output for the ldscore intercept which is mix of characters and numerics
intercept <- as.numeric(str_extract_all(ldsc_output[29], "[0-9.]+")[[1]])

## correcting the P-values if the ldsc regression os significantly different to 1
if ((intercept[1] - intercept[2]*1.96) > 1) {
  
  print(TRUE)
  
  fullT$P_value <- fullT$P_value*intercept[1]
  
  write.table(intercept[1], file=paste0("ldsc_reg/intercept/", trait, ".", sex, ".int"), row.names = F, col.names = F)
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
write.csv(fullT[fullT$SNP %in% snps_signif], file=paste0("GWAS_output/sex_specific/", trait, ".", sex, "gwas.hits"), quote=F, row.names=F)


######## Generating QC plots for the sex specific data

png(paste0("QC/QQ/", trait, "_", sex, "_QQ.png"), height=5, width=6, units="in", res=300)

qq(fullT$P_value, main = paste0("QQ plot of ", trait, " GWAS p-values in ", sex))

lmbx <- qchisq(median(fullT$P_value, na.rm=T),1,lower.tail=FALSE)/qchisq(0.50,1,lower.tail=FALSE)

legend('topleft',legend=bquote("Median"~lambda == .(round(lmbx,2))),bg='white',bty='n',cex=1)

dev.off()

CHR <- as.numeric(fullT$CHR)
BP <- as.numeric(fullT$BP)
P <- fullT$P_value
P[is.na(P)] <- 1

gwas <- cbind.data.frame(CHR,BP,P)


png(paste0("QC/manhat/", trait, "_", sex, "manhat.png"), height=6, width=10, units="in", res=500)

manhattan(gwas)

legend('topleft',legend=paste(trait, sex, sep=" "), bg='white',bty='n',cex=1)


dev.off()

#### end of R script 


