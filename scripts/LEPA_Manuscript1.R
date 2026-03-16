#03/15/26
#data processing and analyses for Manuscript 1: wild ocelot genomics
#data aligned to the new ocelot genome
#vcf from the lab had minimal filtering, only basic quality control, and no imputation

####setting up####
#setting working directories
setwd("/Volumes/Expansion/2_TB_Working_Files") #WD if working out of the hard drive
setwd("~/Documents/Masters_Work/Analyses/1_Data/1_Working_Files") #WD if working off of the laptop
#installing and libraring required packages
library(ggplot2)
library(dplyr)
library(HardyWeinberg)
library(snpStats)
library(adegenet)
library(poppr)
library(vcfR)
library(promises)
library(adegenet)
library(vcfR)
library(poppr)
library(SNPRelate)
library(dartR)
##for karyotype plots
library(karyoploteR)
library(regioneR)
library(GenomicRanges)
library(data.table)
library(IRanges)
library(GenomicAlignments)
##for admixture plots
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)
library(tidyr)
library(stringr)
library(dplyr)

#troubleshooting code for mac
#to get plink to run on the mac, needed to delete the mac quarantine by using the following code in terminal:
    #xattr -d com.apple.quarantine ~/Documents/Masters_Work/Analyses/PLINK_Files/plink2_mac_20260228/plink2

####data pre processing in bcftools in terminal####
#setting working directory to use bcftools
  #cd /Users/tylerbostwick/bcftools
#creating an index file for the vcf, which is needed to run the subsetting
  #./bcftools index /Volumes/Expansion/2_TB_Working_Files/joint_call.LAO03M.20251202.vcf.gz
#subsetting the vcf to just contain the autosomes, the original file also has the sex chromosomes, mtDNA, and unplaced fragments
  #bcftools % ./bcftools view -r A1_RagTag,A2_RagTag,A3_RagTag,B1_RagTag,B2_RagTag,B3_RagTag,B4_RagTag,C1_RagTag,C2_RagTag,C3_RagTag,D1_RagTag,D2_RagTag,D3_RagTag,D4_RagTag,E1_RagTag,E2_RagTag,E3_RagTag \
  #/Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/joint_call.LAO03M.20251202.vcf.gz \
  #-O z -o /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/joint_call_autosomes.vcf.gz

####reading creating plink files from the joint_call VCF file####
system("./plink2 --vcf joint_call_autosomes.vcf.gz --keep-allele-order --allow-extra-chr --vcf-min-dp 10 --max-alleles 2 --chr-set 18 --make-bed --out SNP_AllChrom_AllInd")
      #89 individuals 107697881 variants remain after filter for depth and biallelic

####creating base unfiltered plink files####
#renaming chromosomes to a standard number based, which is what plink is expecting
bim <- read.table("SNP_AllChrom_AllInd.bim", stringsAsFactors = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")
#get the unique chr names
unique_chrs <- unique(bim$chr)
print(unique_chrs)
chr_map <- data.frame(
  old_chr = unique_chrs,
  new_chr = 1:length(unique_chrs),
  stringsAsFactors = FALSE
)
print(chr_map)
#apply that to the bim file directly
for (i in 1:nrow(chr_map)) {
  bim$chr[bim$chr == chr_map$old_chr[i]] <- chr_map$new_chr[i]
}
#write updated bim file
write.table(bim, "SNP_AllInd_unfilt_chrfix.bim", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)
#have plink write new bim bam bed files with the new bim file we created
system("./plink --bed SNP_AllChrom_AllInd.bed --bim SNP_AllInd_unfilt_chrfix.bim --fam SNP_AllChrom_AllInd.fam --make-bed --out SNP_AllInd_AllChrom_biallelic_chrfix")

###give snp ids -- labels every snps with a unique code
system("./plink2 --bfile SNP_AllInd_AllChrom_biallelic_chrfix --chr-set 17 --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out SNP_AllInd_AllChrom_biallelic_uniqueID")

###subset data -- remove mountain lion; all populations and origins get filtered together
system("./plink --bfile SNP_AllInd_AllChrom_biallelic_uniqueID --keep pop_subset_ocelot.txt --chr-set 17 --make-bed --out LEPA_biallelic_dp10")
    #Total genotyping rate in remaining samples is 0.694948; 107697881 variants and 86 samples pass filters and QC.
    #at this point, the only filters applied are coverage depth and biallelic

####SNP filtering to create standard data set for all LEPA####
#apply filters -- maf, miss, hwe> base, other adjustments can follow
###LEPA
system("./plink --bfile LEPA_biallelic_dp10 --chr-set 17 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_standard_final")
      #57793 variants and 86 samples pass filters and QC.

####Export full lepa standard vcf
system("./plink2 --bfile LEPA_standard_final --chr-set 17 --export vcf-4.2 bgz --out LEPA_standard_final")

####Subsetting the groups from the whole LEPA standard final after filtering together####
system("./plink --bfile LEPA_standard_final --keep wild_subset.txt --chr-set 17 --make-bed --out Wild_Standard_postfilt_subset") #wild populations
system("./plink --bfile LEPA_standard_final --keep zoo_subset.txt --chr-set 17 --make-bed --out Zoo_Standard_postfilt_subset") #zoo populations

####additional filtering for select analyses####
#LD pruning
#LEPA
system("./plink --bfile LEPA_standard_final --chr-set 17 --keep-allele-order --indep 50 5 2 --out LEPA_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("./plink --bfile LEPA_standard_final --extract LEPA_LDpruned_0.5_out.prune.in --chr-set 17 --make-bed --out LEPA_LDpruned_05") #extract SNPs and create new files
    #Total genotyping rate is 0.925572; 13161 variants and 86 samples pass filters and QC.
#write vcf
system("./plink2 --bfile LEPA_LDpruned_05 --chr-set 17 --export vcf-4.2 bgz --out LEPA_LDpruned_05")

#ROH filters; no MAF, miss 90, biallelic, coverage depth
#LEPA
system("./plink --bfile LEPA_biallelic_dp10 --chr-set 17 --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_roh_filter")
    #292296 variants and 86 samples pass filters and QC.


##############################################################
####Data Analysis
####PCA Analysis####
