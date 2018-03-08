cd /well/ukbbmcm/workDirectories/zammy/GWAS

## running script that cleans the meta-analysis results slighlty, runs ldscore regression and plots QC graphs
## --- need to chat to Robin about whether to remove SNPS when there is heterogeneity between RCS in SNP effect 

for pheno in BMI HIP WHR HIPadjBMI WHRadjBMI WCadjBMI ; do

for sex in Females Males ; do

Rscript scripts/sex_specific_GWAS_results_UKB.R $pheno $sex &

done

done


## running script that meta-analyses the GWAS data to produce sex-combined results

for pheno in BMI HIP WHR HIPadjBMI WHRadjBMI WCadjBMI ; do

scripts/sex_specific_GWAS_meta-analysis_UKB $pheno

done



## running script to clean the sex-combined meta-analysed data - caluclates ldsc regression, top snps, sex heterogenous snps and does some plots

for pheno in BMI HIPadjBMI WHRadjBMI fat_per ; do

Rscript scripts/sex_combined_GWAS_results.R $pheno &
  
done




## running script to clean the sex-combined meta-analysed data - caluclates ldsc regression, top snps, sex heterogenous snps and does some plots

for pheno in BMI HIPadjBMI WHRadjBMI fat_per ; do

Rscript scripts/locus_zoom_plotting.R $pheno &
  
done
