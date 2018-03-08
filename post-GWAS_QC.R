

### Removing the incorrectly imputed genotype data

grep -Fwf <(sort /well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/hg_2.snps) <(zless /well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_BMI_Females_bgen_stats.gz|awk '{print $1}' |sort) \
> /well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/bgen_gwas.excl



grep -Fwf /well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/bgen_gwas.excl <(zcat /well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_BMI_Females_bgen_stats.gz) \ 
>/well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/bgen_gwas_info.excl

comm -23 <(zcat /well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_BMI_Females_bgen_stats.gz|sort) \
<(sort /well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/bgen_gwas_info.excl) \
> /well/ukbbmcm/workDirectories/zammy/GWAS/test_wei.bgen



for sex in Females Males; do

for trait in BMI HIP WHR WCadjBMI HIPadjBMI WHRadjBMI; do


gunzip -c /well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_$trait"_"$sex"_"bgen_stats.gz | \
grep -F -v -w -f /well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/bgen_gwas.excl \
>/well/ukbbmcm/workDirectories/zammy/GWAS/output/hrc_only/$trait"_"$sex.bgen


done
done


gunzip -c /well/ukbbmcm/workDirectories/zammy/GWAS/output/pooled_BMI_Females_bgen_stats.gz | \
grep -F -v -w -f /well/ukbbmcm/workDirectories/zammy/GWAS/impfiles/hg_2.snps \
>/well/ukbbmcm/workDirectories/zammy/GWAS/test_try2.bgen

vi scripts/QC_plot.r

library(data.table)

setwd("/well/ukbbmcm/workDirectories/zammy/GWAS/")

args <- commandArgs(TRUE)
stats <- args[1] ## summary stat file


##reading in gwas data frame
gwas <- fread(paste0("/well/ukbbmcm/workDirectories/zammy/GWAS/output/hrc_only/", stats, ".bgen"), h=T, stringsAsFactors=F)


qq = function(pvector, ...) {
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( ppoints(length(pvector) ))
  plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
  abline(0,1,col="red")
}


manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  
  if (numchroms==1) {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }	else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...)) 
  }
  
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
}


png(paste0("QC/plots/QQ/", stats, "_QQ.png"), height=600, width=600, res=100)

qq(gwas$P_BOLT_LMM)

lmbx <- qchisq(median(gwas$P_BOLT_LMM,na.rm=T),1,lower.tail=FALSE)/qchisq(0.50,1,lower.tail=FALSE)

legend('topleft',legend=bquote("Median"~lambda == .(round(lmbx,2))),bg='white',bty='n',cex=1)

dev.off()

## subset to significant hits
gwas2 <- subset(gwas,gwas$P_BOLT_LMM < 0.001)

## Just keep CHR, BP, and P
CHR <- as.numeric(gwas2$CHR)
BP <- as.numeric(gwas2$BP)
P <- gwas2$P_BOLT_LMM
P[is.na(P)] <- 1

gwas3 <- cbind.data.frame(CHR,BP,P)

png(paste0("QC/plots/manhat/", stats, "_Manhattan.png"), height=600, width=1200, res=100, pointsize=3)

manhattan(gwas3,pch=16,cex=3,colors=c("lightskyblue2","midnightblue"),suggestiveline=FALSE,cex.axis=3,main='')

dev.off()


for sex in Males Females ; do

for trait in BMI HIP WHR WCadjBMI HIPadjBMI WHRadjBMI ; do

Rscript scripts/QC_plot.r $trait"_"$sex

done
done




vi scripts/QC_plot.r

library(data.table)
library(qqman)

setwd("/well/ukbbmcm/workDirectories/zammy/GWAS/")

args <- commandArgs(TRUE)
trait <- args[1] ## summary stat file
sex <- args[2]

##reading in gwas data frame
gwas <- fread(paste0("/well/ukbbmcm/workDirectories/zammy/GWAS/output/hrc_only/", trait, "_", sex, ".bgen"), h=T, stringsAsFactors=F)

gwas <- subset(gwas, gwas$SE<10)

png(paste0(trait, "_", sex, "_QQ.png"), height=5, width=6, units="in", res=300)

qq(gwas$P_BOLT_LMM, main = paste0("QQ plot of ", trait, " GWAS p-values in ", sex))

lmbx <- qchisq(median(gwas$P_BOLT_LMM,na.rm=T),1,lower.tail=FALSE)/qchisq(0.50,1,lower.tail=FALSE)

legend('topleft',legend=bquote("Median"~lambda == .(round(lmbx,2))),bg='white',bty='n',cex=1)

dev.off()
                                                                              
## subset to significant hits
gwas2 <- subset(gwas,gwas$P_BOLT_LMM < 0.001)

## Just keep CHR, BP, and P
CHR <- as.numeric(gwas2$CHR)
BP <- as.numeric(gwas2$BP)
P <- gwas2$P_BOLT_LMM
P[is.na(P)] <- 1

gwas3 <- cbind.data.frame(CHR,BP,P)


png(paste0("QC/plots/manhat/", trait, "_", sex, "_manhat.png"), height=4.5, width=7, units="in", res=300)

manhattan(gwas3)

dev.off()
 
for sex in Males Females ; do

for trait in HIP WHR WCadjBMI HIPadjBMI WHRadjBMI ; do

Rscript scripts/qqman_plot.r $trait $sex

done
done


