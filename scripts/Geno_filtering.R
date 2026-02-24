###6/19/25
####filtering script -- taking unfiltered data and applying the necessary filters for each analysis type

#directory and package set up
###working directory
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/SNP_Manu_filters")

#####################Creating the base################################################

###create bim bed fam files from original vcf
system("plink --vcf 24041DeY-snp_filter_dp_gq__AllChrom_.vcf.gz --keep-allele-order --allow-extra-chr --chr-set 18 --make-bed --out SNP_dp_gq_AllChrom_AllInd")
###fix chromosome names
#read in bim file
bim <- read.table("SNP_dp_gq_AllChrom_AllInd.bim", stringsAsFactors = FALSE)
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
system("plink --bed SNP_dp_gq_AllChrom_AllInd.bed --bim SNP_AllInd_unfilt_chrfix.bim --fam SNP_dp_gq_AllChrom_AllInd.fam --make-bed --out SNP_AllInd_unfilt_chrfix")

###give snp ids
system("plink2 --bfile SNP_AllInd_unfilt_chrfix --chr-set 18 --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out SNP_AllInd_uniqueID")

###subset populations
#for F, Fis, nucleotide, and He/o, populations will be filtered completely separate
#for ROH, admixture, PCA, and outlier, populations will be filtered completely together
system("plink --bfile SNP_AllInd_uniqueID --keep ranch_subset.txt --chr-set 18 --make-bed --out Ranch_unfilt")
system("plink --bfile SNP_AllInd_uniqueID --keep refuge_subset.txt --chr-set 18 --make-bed --out Refuge_unfilt")
system("plink --bfile SNP_AllInd_uniqueID --keep generic_subset.txt --chr-set 18 --make-bed --out Generic_unfilt")
system("plink --bfile SNP_AllInd_uniqueID --keep brazilian_subset.txt --chr-set 18 --make-bed --out Brazilian_unfilt")
system("plink --bfile SNP_AllInd_uniqueID --keep pop_subset_ocelot.txt --chr-set 18 --make-bed --out LEPA_unfilt") ##mountain lions are removed at this point

####################standard filters for most analyses##############################
####apply filters -- maf, miss, hwe, biallelic> base, other adjustments can follow
###subsetted populations
##ranch
system("plink --bfile Ranch_unfilt --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out Ranch_maf05_miss90_hwe")
system("plink --bfile Ranch_maf05_miss90_hwe --chr-set 18 --keep-allele-order --biallelic-only --make-bed --out Ranch_standard_final")

##refuge
system("plink --bfile Refuge_unfilt --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out Refuge_maf05_miss90_hwe")
system("plink --bfile Refuge_maf05_miss90_hwe --chr-set 18 --keep-allele-order --biallelic-only --make-bed --out Refuge_standard_final")

##generic
system("plink --bfile Generic_unfilt --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out Generic_maf05_miss90_hwe")
system("plink --bfile Generic_maf05_miss90_hwe --chr-set 18 --keep-allele-order --biallelic-only --make-bed --out Generic_standard_final")

##brazilian
system("plink --bfile Brazilian_unfilt --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out Brazilian_maf05_miss90_hwe")
system("plink --bfile Brazilian_maf05_miss90_hwe --chr-set 18 --keep-allele-order --biallelic-only --make-bed --out Brazilian_standard_final")

###LEPA
system("plink --bfile LEPA_unfilt --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_maf05_miss90_hwe")
system("plink --bfile LEPA_maf05_miss90_hwe --chr-set 18 --keep-allele-order --biallelic-only --make-bed --out LEPA_standard_final")

####Export vcf's
system("plink2 --bfile Ranch_standard_final --chr-set 18 --export vcf-4.2 bgz --out Ranch_standard_final")
system("plink2 --bfile Refuge_standard_final --chr-set 18 --export vcf-4.2 bgz --out Refuge_standard_final")
system("plink2 --bfile Generic_standard_final --chr-set 18 --export vcf-4.2 bgz --out Generic_standard_final")
system("plink2 --bfile Brazilian_standard_final --chr-set 18 --export vcf-4.2 bgz --out Brazilian_standard_final")
system("plink2 --bfile LEPA_standard_final --chr-set 18 --export vcf-4.2 bgz --out LEPA_standard_final")

###################standard filters + LD pruning##################################
##Ranch
system("plink --bfile Ranch_standard_final --chr-set 18 --keep-allele-order --indep 50 5 2 --out Ranch_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Ranch_standard_final --extract Ranch_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out Ranch_LDpruned_05") #extract SNPs and create new files
#write vcf
system("plink2 --bfile Ranch_LDpruned_05 --chr-set 18 --export vcf-4.2 bgz --out Ranch_LDpruned_05")

##Refuge
system("plink --bfile Refuge_standard_final --chr-set 18 --keep-allele-order --indep 50 5 2 --out Refuge_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Refuge_standard_final --extract Refuge_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out Refuge_LDpruned_05") #extract SNPs and create new files
#write vcf
system("plink2 --bfile Refuge_LDpruned_05 --chr-set 18 --export vcf-4.2 bgz --out Refuge_LDpruned_05")

##Generic
system("plink --bfile Generic_standard_final --chr-set 18 --keep-allele-order --indep 50 5 2 --out Generic_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Generic_standard_final --extract Generic_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out Generic_LDpruned_05") #extract SNPs and create new files
#write vcf
system("plink2 --bfile Generic_LDpruned_05 --chr-set 18 --export vcf-4.2 bgz --out Generic_LDpruned_05")

##Brazilian
system("plink --bfile Brazilian_standard_final --chr-set 18 --keep-allele-order --indep 50 5 2 --out Brazilian_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Brazilian_standard_final --extract Brazilian_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out Brazilian_LDpruned_05") #extract SNPs and create new files
#write vcf
system("plink2 --bfile Brazilian_LDpruned_05 --chr-set 18 --export vcf-4.2 bgz --out Brazilian_LDpruned_05")

##LEPA
system("plink --bfile LEPA_standard_final --chr-set 18 --keep-allele-order --indep 50 5 2 --out LEPA_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile LEPA_standard_final --extract LEPA_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out LEPA_LDpruned_05") #extract SNPs and create new files
#write vcf
system("plink2 --bfile LEPA_LDpruned_05 --chr-set 18 --export vcf-4.2 bgz --out LEPA_LDpruned_05")

###################standard filters and thinning###################################################
#thinning final vcfs
system("plink --vcf Ranch_standard_final.vcf.gz --chr-set 18 --bp-space 1000 --export vcf bgz --out ranch_final_thin")
system("plink --vcf Refuge_standard_final.vcf.gz --chr-set 18 --bp-space 1000 --export vcf bgz --out refuge_final_thin")
system("plink --vcf Generic_standard_final.vcf.gz --chr-set 18 --bp-space 1000 --export vcf bgz --out generic_final_thin")
system("plink --vcf Brazilian_standard_final.vcf.gz --chr-set 18 --bp-space 1000 --export vcf bgz --out brazilian_final_thin")

######################SAMBA filters#################################################################

####unfiltered and thinned for initial use in Samba####
#thinning and writing vcf for samba -- unfiltered dataset
system("plink --bfile LEPA_unfilt --chr-set 18 --bp-space 1000 --make-bed --out lepa_unfilt_thin")
#export as vcf
system("plink2 --bfile lepa_unfilt_thin --chr-set 18 --export vcf bgz --out lepa_samba_nofilt_thin_final")

####filtering for outlier loci -- starting filtering from the standard unfiltered step at line 45####
#these analyses need the QC filters of maf, missingness, and biallelic
#no HWE, LD pruned, and thinned if necessary after LD pruning
#after population all filtered together, then subsetted for comparison between all groups

#applying intial QC filters
system("plink --bfile LEPA_unfilt --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --make-bed --out LEPA_maf05_miss90_OT")
system("plink --bfile LEPA_maf05_miss90_OT --chr-set 18 --keep-allele-order --biallelic-only --make-bed --out LEPA_OT_QC_filters")

#LD pruning
system("plink --bfile LEPA_OT_QC_filters --chr-set 18 --keep-allele-order --indep 50 5 2 --out LEPA_OT_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile LEPA_standard_final --extract LEPA_OT_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out LEPA_LDpruned_05_OT") #extract SNPs and create new files

#no thinning needed? separate files into pairs for comparison and export as vcf
system("plink --bfile LEPA_LDpruned_05_OT --chr-set 18 --export vcf bgz --out lepa_OT_filt") #lepa doesnt need subsetting, just export
system("plink --bfile LEPA_LDpruned_05_OT --keep gen_braz_subset.txt --chr-set 18 --export vcf bgz --out gen_braz_OT_filt") #generic and brazilian
system("plink --bfile LEPA_LDpruned_05_OT --keep ref_ran_subset.txt --chr-set 18 --export vcf bgz --out ref_ran_OT_filt") #refuge and ranch
system("plink --bfile LEPA_LDpruned_05_OT --keep gen_ref_subset.txt --chr-set 18 --export vcf bgz --out gen_ref_OT_filt") #generics and refuge
system("plink --bfile LEPA_LDpruned_05_OT --keep gen_ran_subset.txt --chr-set 18 --export vcf bgz --out gen_ran_OT_filt") #generics and ranch
system("plink --bfile LEPA_LDpruned_05_OT --keep braz_ref_subset.txt --chr-set 18 --export vcf bgz --out braz_ref_OT_filt") #brazilian and refuge
system("plink --bfile LEPA_LDpruned_05_OT --keep braz_ran_subset.txt --chr-set 18 --export vcf bgz --out braz_ran_OT_filt") #brazilian and ranch


####original standard filters that were thinned for dapc in samba -- not in use anymore####
#thinning for dapc
system("plink --bfile Wild_refiltered --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out wild_dapc_thin")
system("plink --bfile Zoo_refiltered --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out zoo_dapc_thin")
system("plink --bfile LEPA_refiltered --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out lepa_dapc_thin")

system("plink2 --bfile wild_dapc_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out wild_dapc_thin")
system("plink2 --bfile zoo_dapc_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out zoo_dapc_thin")
system("plink2 --bfile lepa_dapc_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out lepa_dapc_thin")

##############After this point, filtering needs to be updated with the new standard files###########
##########################ROH filters###############################################################
##ROH filters no MAF, miss 90, biallelic##
#filter missingness and hwe
system("plink --bfile Zoo_BaseFilter_Allchrom --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Zoo_miss90_hwe")
#filter biallelic
system("plink --bfile Zoo_miss90_hwe --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Zoo_roh_filter")

##Wild filtering
#filter MAF and missingness
system("plink --bfile Wild_BaseFilter_AllChrom --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Wild_miss90_hwe")
#filter biallelic
system("plink --bfile Wild_miss90_hwe --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Wild_roh_filter")


##thinning##
#--bp-space thins data by spacing out snps by 1000bps apart
system("plink --bfile Wild_roh_filter --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out wild_roh_thin")
#export as vcf
system("plink2 --bfile wild_roh_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out wild_thinned_roh")
