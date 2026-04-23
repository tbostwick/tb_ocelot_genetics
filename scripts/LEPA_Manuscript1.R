#04/22/26
#data processing and analyses for Manuscript 1: wild ocelot genomics
#data aligned to the new ocelot genome
#vcf from the lab had minimal filtering, only basic quality control, and no imputation

###~~~ Filtering steps~~~
#First in BCFTools:
    #select only autosomes
    #remove indels and non biallelic snps
    #filter for a genotype quality score above 9
#next, in plink for processesing ease:
    #rename chromosomes
    #give snps unique id's
    #for manuscript 1 -- subset just the wild
#then, filtering in plink for each analysis
    #depth of 7
    #generally, maf, miss, and hwe filters applied here

####setting up####
#setting working directories
setwd("/Volumes/Expansion/2_TB_Working_Files/Plink_files") #WD if working out of the hard drive
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


################################################################################
####data pre processing in bcftools in terminal####
#setting working directory to use bcftools
  #cd /Users/tylerbostwick/bcftools
#creating an index file for the vcf, which is needed to run the subsetting
  #./bcftools index /Volumes/Expansion/2_TB_Working_Files/joint_call.LAO03M.20251202.vcf.gz
#subsetting the vcf to just contain the autosomes, the original file also has the sex chromosomes, mtDNA, and unplaced fragments
  #bcftools % ./bcftools view -r A1_RagTag,A2_RagTag,A3_RagTag,B1_RagTag,B2_RagTag,B3_RagTag,B4_RagTag,C1_RagTag,C2_RagTag,C3_RagTag,D1_RagTag,D2_RagTag,D3_RagTag,D4_RagTag,E1_RagTag,E2_RagTag,E3_RagTag \
  #/Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/joint_call.LAO03M.20251202.vcf.gz \
  #-O z -o /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/joint_call_autosomes.vcf.gz

#removing indels and non biallelic snps from the vcf files in bcftools
  #./bcftools norm -m -any \ /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/joint_call_autosomes.vcf.gz \ | ./bcftools view -v snps -m2 -M2 \
  #-O z -o /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz
          #this code first normalizes the data set, separating the individual alleles. This allows the code to then remove indels (as some indels can be an alt allele to a snp),
          #then it removes any sites that have only 1 allele, and any that the have more than two -- essentially the biallelic filter
          #in total, this allows the indel filter to catch all instances of indels, then remove uninformative sites and multiallelic sites leaving only biallelic snps

#filtering for genotype quality scores in bcftools -- score of 9 used as min
#./bcftools filter -e 'FMT/GQ < 9' \ /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz \ -O z -o allInd_SNPs_autosomes_bi_gq9.vcf.gz
      #code tells bcftools to exclude (-e) and snps that have a genotype score or less than 9, then reads in the file, and outputs the new filtered file

####creating plink files from the base file created from BCFtools####
system("./plink2 --vcf allInd_SNPs_autosomes.vcf.gz --keep-allele-order --allow-extra-chr --vcf-min-dp 7 --vcf-max-dp 22 --max-alleles 2 --chr-set 17 --make-bed --out SNP_AllChrom_AllInd_dp7_gq9_bi")
      #89 individuals 103792041 variants remain after filter for depth and biallelic

####creating base plink files pre-standard filtering####
##renaming chromosomes to a standard number based, which is what plink is expecting
bim <- read.table("SNP_AllChrom_AllInd_dp7_gq9_bi.bim", stringsAsFactors = FALSE)
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
system("./plink --bed SNP_AllChrom_AllInd_dp7_gq9_bi.bed --bim SNP_AllInd_unfilt_chrfix.bim --fam SNP_AllChrom_AllInd_dp7_gq9_bi.fam --make-bed --out SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix")

###give snp ids -- labels every snps with a unique code
system("./plink2 --bfile SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix --chr-set 17 --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix_uniqueID")

###subset data -- remove mountain lion and duplicates, for manuscript 1 keeping only wild populations
system("./plink --bfile SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix_uniqueID --keep wild_subset.txt --chr-set 17 --make-bed --out SNP_wild_dp7_gq9_bi")
    #Total genotyping rate in remaining samples is 0.861001; 103792041 variants and 44 samples pass filters and QC.
    #at this point, coverage depth of 7, genotype quality of 9, and biallelic filters applied to create base file

####SNP filtering to create standard data set for all wild LEPA####
#apply filters -- maf, miss, hwe> base, other adjustments can follow
###Wild
system("./plink --bfile SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out wild_standard_final")
      #66269934 variants removed due to missing genotype data (--geno)
      #3576 variants removed due to Hardy-Weinberg exact test
      #33435588 variants removed due to minor allele threshold(s)
      #4082943 variants and 44 samples pass filters and QC.

####Export full wild standard vcf
system("./plink2 --bfile wild_standard_final --chr-set 17 --export vcf-4.2 bgz --out wild_standard_final")

####additional filtering for select analyses####
#LD pruning
system("./plink --bfile wild_standard_final --chr-set 17 --keep-allele-order --indep-pairwise 50 5 0.5 --out wild_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
      #change to indep-pairwise instead of indep
system("./plink --bfile wild_standard_final --extract wild_LDpruned_0.5_out.prune.in --chr-set 17 --make-bed --out wild_LDpruned_05") #extract SNPs and create new files
    #Total genotyping rate is 0.925318; 273852 variants and 44 samples pass filters and QC..
#write vcf
system("./plink2 --bfile wild_LDpruned_05 --chr-set 17 --export vcf-4.2 bgz --out wild_LDpruned_05")


#ROH and kinship filters; no MAF, miss 90, biallelic, coverage depth 7, genotype quality 9
#also used in kinship analyses; as is recommended not to filter for MAF or LD prune
system("./plink --bfile SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out wild_kin_roh_filter")
    #37518531 variants and 44 samples pass filters and QC.


##############################################################
####Data Analysis
####PCA Analysis -- Done!####
##PCA on Wild individuals --- done!
#reading in data
system("./plink --bfile wild_standard_final --pca --chr-set 17 --out pca_wild") #run the pca code, specified the number of chromosomes
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
            aes(label=ID), vjust=1.5, hjust=-0.1, size=4, color = "black") + 
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
##king-robust kingship estimator for WILD individuals 
system("./plink2 --bfile wild_kin_roh_filter --make-king-table --allow-extra-chr --chr-set 17 --out wild_KING_manu")
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
#using LD pruned dataset
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
####ROH -- not dome####
####Diversity stats -- not done####


####coverage stats from bcftools - can delete later to clean####
###THIRD iteration, running diagnostics in bcftools

#getting distribution of coverage depths in bcftools:
  #cd /Users/tylerbostwick/bcftools #setting working directory for bcftools
  # ./bcftools stats /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/joint_call_autosomes.vcf.gz | grep "^SN"
        #gives summary stats on the vcf file
        #output:
        #SN	0	number of samples:	89
        #SN	0	number of records:	120040661 -- total number of sites
        #SN	0	number of no-ALTs:	0 -- rows with no alt allele, essentially uninformative rows
        #SN	0	number of SNPs:	101914847 -- total number of snps
        #SN	0	number of MNPs:	0 -- total of multi-nucleotide polys
        #SN	0	number of indels:	20246516 -- total number of indels (need to remove)
        #SN	0	number of others:	0 -- anything else
        #SN	0	number of multiallelic sites:	12342780 -- sites with more than two alleles
        #SN	0	number of multiallelic SNP sites:	1673221 -- number of multiallelic sites that are snps
#after removing all indels and biallelic sites -- total snps left is 103792041
#computing mean coverage depth across sites, code queries the vcf for the depth at each site, the performs a pipe that calculates mean
#./bcftools query -f '[%DP\n]' /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz | \ awk '{sum+=$1; n++} END {print "Mean per-individual depth:", sum/n}'
      #output: Mean per-individual depth: 12.0631
#computing mean and standard deviation bounds for use in determining the depth filter:
#./bcftools query -f '[%DP\n]' /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz | \ awk '$1 > 0 {sum+=$1; sumsq+=$1*$1; n++} END {mean=sum/n; sd=sqrt(sumsq/n - mean^2); print "Mean:", mean; print "SD:", sd; print "Mean - 2SD:", mean-2*sd; print "Mean + 2SD:", mean+2*sd}'
      #output: Mean: 12.1187
              #SD: 5.08218
              #Mean - 2SD: 1.95431
              #Mean + 2SD: 22.2831






