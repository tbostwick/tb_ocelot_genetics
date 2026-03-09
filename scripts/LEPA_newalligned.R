setwd("/Volumes/Expansion/TB_Working_Files")
##looking at the new genome files for with the ocelot alignment
##reading in and performing basic tasks with the file to see how the data is organized
##creating filtered file to use on the rest of the analyses -- so can run from the computer rather than the hard drive
##will need to update main code with the new file names reflected -- so that it is contiguous, even if the actual code is ran from here


#to get plink to run on the mac, needed to delete the mac quarantine by using the following code in terminal:
    #xattr -d com.apple.quarantine ~/Documents/Masters_Work/Analyses/PLINK_Files/plink2_mac_20260228/plink2


####reading creating plink files from the joint_call VCF file
system("./plink2 --vcf joint_call.LAO03M.20251202.vcf.gz --keep-allele-order --allow-extra-chr --vcf-min-dp 10 --max-alleles 2 --chr-set 18 --make-bed --out SNP_AllChrom_AllInd")
      #89 individuals 112237606 variants remain after filter for depth and biallelic

####running initial file renaming to create the base file on which the filters will be applied
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
print(chr_map) ###weird ass chromosome names..... dunno what to do with that tbh, may need to remove this step and just name the variants
#made it to here with the new data


#apply that to the bim file directly
for (i in 1:nrow(chr_map)) {
  bim$chr[bim$chr == chr_map$old_chr[i]] <- chr_map$new_chr[i]
}
#write updated bim file
write.table(bim, "SNP_AllInd_unfilt_chrfix.bim", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)
#have plink write new bim bam bed files with the new bim file we created
system("plink --bed SNP_dp_gq_biallelic_AllChrom_AllInd.bed --bim SNP_AllInd_unfilt_chrfix.bim --fam SNP_dp_gq_biallelic_AllChrom_AllInd.fam --make-bed --out SNP_AllInd_cd_biallelic_chrfix")

###give snp ids -- labels every snps with a unique code
system("plink2 --bfile SNP_AllInd_cd_biallelic_chrfix --chr-set 18 --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out SNP_AllInd_cd_biallelic_uniqueID")

###subset data -- remove mountain lion; all populations and origins get filtered together
system("plink --bfile SNP_AllInd_cd_biallelic_uniqueID --keep pop_subset_ocelot.txt --chr-set 18 --make-bed --out LEPA_cd_biallelic")