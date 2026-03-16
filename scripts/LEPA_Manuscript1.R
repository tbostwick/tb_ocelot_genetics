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

####creating plink files from the joint_call VCF file####
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
    #Total genotyping rate in remaining samples is 0.694513; 107697881 variants and 85 samples pass filters and QC.
    #at this point, the only filters applied are coverage depth and biallelic

####SNP filtering to create standard data set for all LEPA####
#apply filters -- maf, miss, hwe> base, other adjustments can follow
###LEPA
system("./plink --bfile LEPA_biallelic_dp10 --chr-set 17 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_standard_final")
      #61480 variants and 85 samples pass filters and QC

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
    #Total genotyping rate is 0.925572; 13965 variants and 85 samples pass filters and QC.
#write vcf
system("./plink2 --bfile LEPA_LDpruned_05 --chr-set 17 --export vcf-4.2 bgz --out LEPA_LDpruned_05")

#ROH filters; no MAF, miss 90, biallelic, coverage depth
#LEPA
system("./plink --bfile LEPA_biallelic_dp10 --chr-set 17 --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_roh_filter")
    #309052 variants and 85 samples pass filters and QC.


##############################################################
####Data Analysis
####PCA Analysis -- Done!####
##PCA on Wild individuals --- done!
#reading in data
system("./plink --bfile Wild_Standard_postfilt_subset --pca --chr-set 17 --out pca_wild") #run the pca code, specified the number of chromosomes
pca.data.wild <- read.table("pca_wild.eigenvec", header=FALSE) #load pca results from the eigenvec file
colnames(pca.data.wild)[2] <- "ID" #rename second column to allow for appending pop assignment later
eigenvalues <- read.table("pca_wild.eigenval", header=FALSE)$V1 #read in eigenvalues to calc variance
#calc variance
total_variance <- sum(eigenvalues)
pc1_variance <- round((eigenvalues[1] / total_variance) * 100, 2)
pc2_variance <- round((eigenvalues[2] / total_variance) * 100, 2)
#adding population information
wild_origins <- read.csv("wild_origins.csv", header = TRUE) #read in origins file
pca.data.wild.origins <- left_join(pca.data.wild, wild_origins, by = "ID") #append origin data to pca data
pca.data.wild.origins <- pca.data.wild.origins %>%
  mutate(ID = gsub("-.*", "", ID))
#plotting
wild_pca <- ggplot(pca.data.wild.origins, aes(x=V3,y=V4)) +  #plot with individual ID's and by origin
  geom_point(aes(shape = Pop, color = Pop), size = 3) +
  geom_rect(data = subset(pca.data.wild.origins, ID %in% c("E35M", "LO01F")),
            aes(xmin = min(V3) - 0.02, xmax = max(V3) + 0.04,
                ymin = min(V4) - 0.02, ymax = max(V4) + 0.04),
            fill = NA, color = "black", linewidth = 1, linetype = "dashed") +
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("LO03M")), #in correct position
            aes(label=ID), vjust=-0.6, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("OM331")), #in correct position
            aes(label=ID), vjust=-0.6, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("E35M")), #in correct position
            aes(label=ID), vjust=-0.6, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("LO01F")), #in correct position
            aes(label=ID), vjust=-1, hjust=0.02, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("E29M")), #in correct position
            aes(label=ID), vjust=-0.5, hjust=-0.2, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("E32M")),  #in correct position
            aes(label=ID), vjust=-0.4, hjust=-0.2, size=4, color = "black") + 
  scale_color_manual(name = "Origin", values = c("Refuge" = "#01004c", "Ranch" = "orchid")) +
  scale_shape_manual(name = "Origin", values = c("Refuge" = 19, "Ranch" = 17)) +
  labs(x = paste0("PC1 (", pc1_variance, "%)"),
       y = paste0("PC2 (", pc2_variance, "%)")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
wild_pca
ggsave("wild_pca.png", dpi = 300)

####Kinship -- mostly done, needs accuracy check####
##king-robust kingship estimator for WILD individuals -- come back and recheck to accuracy? why all related?
##kin values definitely over-representing the relationships, will need to address
system("./plink2 --bfile Wild_Standard_postfilt_subset --make-king-table --allow-extra-chr --chr-set 17 --out wild_KING_manu")
king.wild.matrix <- read.table("wild_KING_manu.kin0", header = FALSE) #make king table into object
colnames(king.wild.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
king.wild.matrix <- king.wild.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(king.wild.matrix, "wild_pairwise_kinship_manu.csv")

#heat map of kinship for wild individuals -- binned
ggplot(data = king.wild.matrix, aes(x=IID1, y=IID2, fill = KINSHIP)) +
  geom_tile(color="white") +
  scale_fill_stepsn(name = "Kinship", breaks = c(0, 0.04, 0.1, 0.2),
                    limit = c(0, 0.3),
                    labels = c("Unrelated", "3rd Degree", "2nd Degree", "1st Degree"),
                    aesthetics = "fill",
                    space = "Lab",
                    colors = c("grey100", "gold1", "orangered", "firebrick"))+
  labs(x="Individuals", y="Individuals")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


####Admixture -- not done####
##admixture plotting
#.q files from ADMIXTURE performed on a mac computer; code as follows:





###wild admixture plot
#reading in and preparing data
w_k2_table <- read.table("Wild_rf_standard.2.q")
w_ids <- read.table("wild_pop_samp_id.txt")
w_ids = rename(w_ids, sample_id = V1, pop_id = V2)
w_k2 <- cbind(w_ids, w_k2_table)
w_k2_long <- pivot_longer(
  w_k2, cols = c(V1, V2),
  names_to = "ancestry",
  values_to = "proportion"
) #pivots the table to the proportions are able to be plotted as stacked bars
w_k2_long$pop_id <- str_to_upper(w_k2_long$pop_id)
#plot code
w_k2plot <-
  ggplot(w_k2_long, aes(x = factor(sample_id), y = proportion, fill = ancestry)) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(pop_id), switch = "both", scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 12),
        strip.text.x = element_text(face = "bold", size = 15),
        strip.placement = "outside",
        axis.text.y = element_text(size = 15),
  ) +
  scale_fill_manual(values = c("V1" = "#01004c", "V2" = "#ffb2b0")) +
  guides(fill = "none")
w_k2plot
ggsave("wild_k2_admixture_Oct2025.png", w_k2plot, width = 15, height = 8, bg = "white")
