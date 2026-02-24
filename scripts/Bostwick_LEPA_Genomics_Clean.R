##LEPA Manuscript Code Cleaned
#8/20/25
#Cleaned code for the analyses completed in Bostwick 2025
#analyses run with the following settings:
#MAF: 0.05, MISS: 90%, Biallelic, LD: 50 5 2 (corresponds to r2 of 0.5), HWE, read depth 10
#contains both the code for PLINK, adegenet DAPC, and Samba
##########################################################
####package install and working directory set up####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/SNP_LEPA_Manu_FelisFelisAligned")
##install packages if needed, below are packages that require specific code to install
BiocManager::install("snpStats")
install.packages("promises", version = ">= 1.3.2") #got error with incompatible versions, updating this
BiocManager::install("SNPRelate")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("karyoploteR", quietly = TRUE))
  BiocManager::install("karyoploteR")
if (!requireNamespace("regioneR", quietly = TRUE))
  BiocManager::install("regioneR")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE))
  BiocManager::install("IRanges")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
##library packages -- general packages needed
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

##########################################################
#Data Prep, filtering and LD Pruning
####Make bed bim fam files from vcf####
#includes read depth filter and biallelic filter, to remove snps with reads below 10x
system("plink2 --vcf 24041DeY-snp_filter_dpgq__AllChrom.vcf.gz --keep-allele-order --allow-extra-chr --vcf-min-dp 10 --max-alleles 2 --chr-set 18 --make-bed --out SNP_dp_gq_biallelic_AllChrom_AllInd")
####FINAL FILTERING FOR FUTURE ANALYSES, SUBSETTING-FILTERING-EXPORTING####
###fix chromosome names, changes chromosome names from the FCA naming to standard
#read in bim file
bim <- read.table("SNP_dp_gq_biallelic_AllChrom_AllInd.bim", stringsAsFactors = FALSE)
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
system("plink --bed SNP_dp_gq_biallelic_AllChrom_AllInd.bed --bim SNP_AllInd_unfilt_chrfix.bim --fam SNP_dp_gq_biallelic_AllChrom_AllInd.fam --make-bed --out SNP_AllInd_cd_biallelic_chrfix")

###give snp ids -- labels every snps with a unique code
system("plink2 --bfile SNP_AllInd_cd_biallelic_chrfix --chr-set 18 --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out SNP_AllInd_cd_biallelic_uniqueID")

###subset data -- remove mountain lion; all populations and origins get filtered together
system("plink --bfile SNP_AllInd_cd_biallelic_uniqueID --keep pop_subset_ocelot.txt --chr-set 18 --make-bed --out LEPA_cd_biallelic")

####apply filters -- maf, miss, hwe> base, other adjustments can follow

###LEPA
system("plink --bfile LEPA_cd_biallelic --chr-set 18 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_standard_final")

####Export vcf's
system("plink2 --bfile LEPA_standard_final --chr-set 18 --export vcf-4.2 bgz --out LEPA_standard_final")

####Subsetting the groups from the whole LEPA standard final after filtering together####
system("plink --bfile LEPA_standard_final --keep wild_subset.txt --chr-set 18 --make-bed --out Wild_Standard_postsubset") #wild populations
system("plink --bfile LEPA_standard_final --keep gen_braz_subset.txt --chr-set 18 --make-bed --out Zoo_Standard_postsubset") #zoo populations


####FINAL FILTERING -- LINKAGE DISEQUILIBRIUM PRUNING####
##LEPA
system("plink --bfile LEPA_standard_final --chr-set 18 --keep-allele-order --indep 50 5 2 --out LEPA_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile LEPA_standard_final --extract LEPA_LDpruned_0.5_out.prune.in --chr-set 18 --make-bed --out LEPA_LDpruned_05") #extract SNPs and create new files
#write vcf
system("plink2 --bfile LEPA_LDpruned_05 --chr-set 18 --export vcf-4.2 bgz --out LEPA_LDpruned_05")

####Filtering -- ROH filters no MAF, miss 90, biallelic, coverage depth####
#LEPA
system("plink --bfile LEPA_cd_biallelic --chr-set 18 --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_roh_filter")


####thinning#### -- might not need any more with the new coverage depth filter
#--bp-space thins data by spacing out snps by 1000bps apart
#for ROH -- using ROH filters
system("plink --bfile LEPA_roh_filter --chr-set 18 --bp-space 1000 --make-bed --out LEPA_thin")
#export as vcf
system("plink2 --bfile LEPA_roh_thin --chr-set 18  --export vcf bgz --out LEPA_thinned")

#thinning and writing vcf for samba -- unfiltered dataset -- check to make sure the biallelic filter works with samba
system("plink --bfile LEPA_cd_biallelic --chr-set 18 --bp-space 1000 --make-bed --out lepa_unfilt_thin")
#export as vcf
system("plink2 --bfile lepa_unfilt_thin --chr-set 18 --export vcf bgz --out lepa_samba_nofilt_thin_final")

#########################################################
#Data Analysis
####PCA####
#beyond just making a pca, this is good to check if the data was input correctly
#PCA on zoo individuals
system("plink --bfile Zoo_Standard_postsubset --pca --chr-set 18 --out pca_zoo") #run the pca code, specified the number of chromosomes
pca.data.zoo <- read.table("pca_zoo.eigenvec", header=FALSE) #load pca results from the eigenvec file
zoo_origins <- read.csv("zoo_origins.csv", header = TRUE) #read in origin metadata
pca.data.zoo.origins <- left_join(pca.data.zoo, zoo_origins, by = "V2") #append origin data to pca data
ggplot(pca.data.zoo.origins, aes(x=V3,y=V4, color = Cat.Group)) +  #plot with individual ID's and by origin
  geom_point() +
  geom_text(aes(label=V2), vjust=1, hjust=1, size=2) + #V2 is ID
  scale_color_manual(name = "Origin", values = c("Brazilian" = "dodgerblue3", "Generic" = "coral2")) +
  labs(x = "PC1", y = "PC2", title = "Zoo PCA by Individual") +
  theme_minimal()
#variance explained bar plot
pca.variance <- read.table("pca_zoo.eigenval", header = FALSE) #load pca variance results
hist(pca.variance$V1) ##not working right

##PCA on Wild individuals - updated for thesis
#reading in data
system("plink --bfile Wild_Standard_postsubset --pca --chr-set 18 --out pca_wild") #run the pca code, specified the number of chromosomes
pca.data.wild <- read.table("pca_wild.eigenvec", header=FALSE) #load pca results from the eigenvec file
eigenvalues <- read.table("pca_wild.eigenval", header=FALSE)$V1 #read in eigenvalues to calc variance
#calc variance
total_variance <- sum(eigenvalues)
pc1_variance <- round((eigenvalues[1] / total_variance) * 100, 2)
pc2_variance <- round((eigenvalues[2] / total_variance) * 100, 2)
#adding population information
wild_origins <- read.csv("wild_origins.csv", header = TRUE) #read in origins file
pca.data.wild.origins <- left_join(pca.data.wild, wild_origins, by = "V2") #append origin data to pca data
pca.data.wild.origins <- pca.data.wild.origins %>%
  mutate(V2 = gsub("-.*", "", V2))
#plotting
wild_pca <- ggplot(pca.data.wild.origins, aes(x=V3,y=V4)) +  #plot with individual ID's and by origin
  geom_point(aes(shape = Cat.Group, color = Cat.Group), size = 3) +
  geom_rect(data = subset(pca.data.wild.origins, V2 %in% c("E35M", "LO01F")),
            aes(xmin = min(V3) - 0.02, xmax = max(V3) + 0.04,
                ymin = min(V4) - 0.02, ymax = max(V4) + 0.04),
            fill = NA, color = "black", linewidth = 1, linetype = "dashed") +
  geom_text(data = subset(pca.data.wild.origins, V2 %in% c("E35M", "LO01F", "OM331", "LO03M")),
            aes(label=V2), vjust=-1, hjust=-0.2, size=4, color = "black") + #V2 is ID\
  geom_text(data = subset(pca.data.wild.origins, V2 %in% c("E29M")),
            aes(label=V2), vjust=1.3, hjust=-0.2, size=4, color = "black") + #V2 is ID
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

##PCA on all ocelots -- updated for thesis
#read in data
system("plink --bfile LEPA_standard_final --pca --chr-set 18 --allow-extra-chr --out pca_LP_final") #run the pca code, specified the number of chromosomes
pca.data.LP <- read.table("pca_LP_final.eigenvec", header=FALSE) #load pca results from the eigenvec file
eigenvalues <- read.table("pca_LP_final.eigenval", header=FALSE)$V1 #read in eigenvalues to calc variance
#calc variance
total_variance <- sum(eigenvalues)
pc1_variance <- round((eigenvalues[1] / total_variance) * 100, 2)
pc2_variance <- round((eigenvalues[2] / total_variance) * 100, 2)
#adding pop info
ocelot.origins <- read.csv("LP_origins.csv", header = TRUE) #read in origins file
pca.data.ocelot.origins <- left_join(pca.data.LP, ocelot.origins, by = "V2") #append origin data to pca data
pca.data.ocelot.origins <- pca.data.ocelot.origins %>%
  mutate(V2 = gsub("-.*", "", V2))
#plot
lepa_pca <- ggplot(pca.data.ocelot.origins, aes(x=V3,y=V4, color = Cat.Group)) +  #plot with individual ID's and by origin
  geom_point(aes(shape = Cat.Group, color = Cat.Group), size = 3) +
  scale_color_manual(name = "Origin", values = c("Refuge" = "purple4", "Ranch" = "orchid", 
                                                 "Brazilian" = "#3B967f", "Generic" = "#D66857")) +
  scale_shape_manual(name = "Origin", values = c("Refuge" = 19, "Ranch" = 19,
                                                 "Brazilian" = 17, "Generic" = 17)) +
  labs(x = paste0("PC1 (", pc1_variance, "%)"),
       y = paste0("PC2 (", pc2_variance, "%)")) +
  theme_minimal()+
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
lepa_pca


####KING Values####
#king-robust kingship estimator for zoo individuals -- UPDATED FOR THESIS
system("plink2 --bfile Zoo_Standard_postsubset --make-king-table --allow-extra-chr --chr-set 18 --out zoo_KING_manu")
zoo.king.matrix <- read.table("zoo_KING_manu.kin0", header = FALSE) #make king table into object
colnames(zoo.king.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
zoo.king.matrix <- zoo.king.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(zoo.king.matrix, "zoo_pairwise_kinship_manu.csv")

#heat map of kinship for zoo individuals -- binned
zoo.king.matrix <- read.csv("zoo_pairwise_kinship_manu.csv", header = TRUE) #read in data from previous step
zoo.king.matrix$IID1 <- as.factor(zoo.king.matrix$IID1) #changing the id names to factors to plot correctly
zoo.king.matrix$IID2 <- as.factor(zoo.king.matrix$IID2)
ggplot(data = zoo.king.matrix, aes(x=IID1, y=IID2, fill = KINSHIP)) +
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

#value distribution plot
hist(zoo.king.matrix$KINSHIP)

#king-robust kingship estimator for wild individuals
system("plink2 --bfile Wild_Standard_postsubset --make-king-table --allow-extra-chr --chr-set 18 --out wild_KING_manu")
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

#king-robust kingship estimator for ranch individuals -- dispersers from the ranch included
system("plink2 --bfile Ranch_Refiltered --allow-extra-chr --make-king-table --chr-set 18 --out ranch_KING")
king.ranch.matrix <- read.table("ranch_KING.kin0", header = FALSE) #make king table into object
colnames(king.ranch.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
king.ranch.matrix <- king.ranch.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(king.ranch.matrix, "ranch_pairwise_kinship_refiltered.csv")

#heat map of kinship for ranch individuals -- binned
melted_kin_ranch <- melt(king.ranch.matrix, measure.vars = "KINSHIP", variable.name = "KINSHIP")
ggplot(data = melted_kin_ranch, aes(x=IID1, y=IID2, fill = value)) +
  geom_tile(color="white") +
  scale_fill_stepsn(name = "Kinship", breaks = c(0, 0.04, 0.1, 0.2),
                    limit = c(0, 0.3),
                    labels = c("Unrelated", "First Cousins", "Half Sibs", "Full Sib/PO"),
                    aesthetics = "fill",
                    space = "Lab",
                    colors = c("grey100", "gold1", "orangered", "firebrick"))+
  labs(x="Individuals", y="Individuals")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#king-robust kingship estimator for refuge individuals -- dispersers from the ranch included
system("plink2 --bfile Refuge_Refiltered --allow-extra-chr --make-king-table --chr-set 18 --out refuge_KING")
king.refuge.matrix <- read.table("refuge_KING.kin0", header = FALSE) #make king table into object
colnames(king.refuge.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
king.refuge.matrix <- king.refuge.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(king.refuge.matrix, "refuge_pairwise_kinship_refiltered.csv")

#heat map of kinship for refuge individuals -- binned
melted_kin_refuge <- melt(king.refuge.matrix, measure.vars = "KINSHIP", variable.name = "KINSHIP")
ggplot(data = melted_kin_refuge, aes(x=IID1, y=IID2, fill = value)) +
  geom_tile(color="white") +
  scale_fill_stepsn(name = "Kinship", breaks = c(0, 0.04, 0.1, 0.2),
                    limit = c(0, 0.3),
                    labels = c("Unrelated", "First Cousins", "Half Sibs", "Full Sib/PO"),
                    aesthetics = "fill",
                    space = "Lab",
                    colors = c("grey100", "gold1", "orangered", "firebrick"))+
  labs(x="Individuals", y="Individuals")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

####ROH Selection -- PLINK####
#wild --thinned, third try of stricter saremi params for thinned data ---BEST SO FAR
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-het 20 --homozyg-window-missing 20 --homozyg-window-snp 100
       --homozyg-window-threshold 0.02 --out wild_roh_thin_test2")
wild_roh_thin_test2_results <- read.table("wild_roh_thin_test2.hom.indiv", header = T)

#zoo --thinned, third try of stricter saremi params for thinned data ---BEST SO FAR
system("plink --bfile zoo_dapc_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-het 20 --homozyg-window-missing 20 --homozyg-window-snp 100
       --homozyg-window-threshold 0.02 --out zoo_roh_thin_test2")
zoo_roh_thin_test2_results <- read.table("zoo_roh_thin_test2.hom.indiv", header = T)

####ROH assessment####
#average % genome in roh by population
population_mean_roh <- roh.df %>%
  group_by(pop_id) %>%
  summarise(Average = mean(Thin_test2, na.rm = TRUE),
            Count = n(),
            StdDev = sd(Thin_test2, na.rm = TRUE),
            StdError = (sd(Thin_test2, na.rm = TRUE)/sqrt(n())))
#plot1
ggplot() +
  geom_boxplot(data = roh.df,
               aes(x = pop_id, y = Thin_test2),
               width = 0.5, alpha = 0.7, outlier.shape = 1) +
  geom_point(data = population_mean_roh,
             aes(x = pop_id, y = Average),
             color = "red", size = 3) +
  geom_errorbar(data = population_mean_roh,
                aes(x = pop_id, ymin = Average - StdError, ymax = Average + StdError),
                width = 0.2, color = "red", linewidth = 1) +
  labs(title = "Distibution of ROH % by Population",
       x = "Population",
       y = "Percent genome in ROH") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))

#violin plot of wild FROH
l_pop <- read.csv("lepa_origins.csv")
l_pop <- rename(l_pop, ID = individual)
roh.df <- read.csv("wild_roh_percentgenome.csv", header = TRUE) #read in data
population_mean_roh <- read.csv("pop_mean_froh_merged.csv", header = TRUE) #read in data
roh.df <- rename(roh.df, ID = roh.df)
roh.df <- left_join(roh.df, l_pop, by = "ID")
ggplot() +
  # violin plot
  geom_violin(data = roh.df, aes(x = Pop, y = Thin_test2, fill = Pop), 
              alpha = 0.7) +
  # Add individual points
  geom_jitter(data = roh.df, aes(x = Pop, y = Thin_test2), 
              width = 0.1, alpha = 0.4, size = 1) +
  # Add population means with error bars
  geom_point(data = population_mean_roh, aes(x = pop_id, y = Average), 
             color = "red", size = 4, shape = 18) +
  geom_errorbar(data = population_mean_roh, 
                aes(x = pop_id, y = Average, 
                    ymin = Average - 2*StdDev, ymax = Average + 2*StdDev),
                color = "red", width = 0.2, size = 1) +
  scale_fill_manual(values = c("Ranch" = "#ffb2b0", "Refuge" = "#01004c")) +
  # Labels and theme
  labs(x = "Population", y = expression(F[ROH])) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")

#zoo -- populate dataframe with ID's and percent genome for each parameter test
zoo_roh <- read.table("zoo_roh_thin_test2.hom.indiv", header = TRUE)
roh.df.z <-zoo_roh$IID #populate ID's
roh.df.z <- as.data.frame(roh.df.z) #create the data frame
roh.df.z$Thin_test2 <- ((zoo_roh$KB*1000)/2425730029)*100
write.csv(roh.df.z, "zoo_roh_percentgenome.csv")
z_pop_id <-read.csv("zoo_origins.csv")
colnames(z_pop_id)[1] <- "IID"
colnames(roh.df.z)[1] <- "IID"
roh_zoo <- merge(roh.df.z, z_pop_id, by = "IID", all = TRUE)

origin_mean_roh <- roh_zoo %>%
  group_by(Cat.Group) %>%
  summarise(Average = mean(Thin_test2, na.rm = TRUE),
            Count = n(),
            StdDev = sd(Thin_test2, na.rm = TRUE),
            StdError = (sd(Thin_test2, na.rm = TRUE)/sqrt(n())))
#checking outliers
merged_df <- merge(roh_zoo, origin_mean_roh, by = "Cat.Group")
merged_df$z_score <- abs((merged_df$Thin_test2 - merged_df$Average) / merged_df$StdDev)
outliers_1sd <- sum(merged_df$z_score > 1, na.rm = TRUE)
outliers_2sd <- sum(merged_df$z_score > 2, na.rm = TRUE)
outliers_3sd <- sum(merged_df$z_score > 3, na.rm = TRUE)

outlier_summary <- data.frame(
  standard_deviations = c(1, 2, 3),
  count_outside = c(outliers_1sd, outliers_2sd, outliers_3sd),
  percent_outside = c(outliers_1sd/nrow(merged_df) * 100,
                      outliers_2sd/nrow(merged_df) * 100,
                      outliers_3sd/nrow(merged_df) * 100)
)

outliers_by_group <- merged_df %>%
  group_by(Cat.Group) %>%
  summarise(
    total_individuals = n(),
    outside_1sd = sum(z_score > 1, na.rm = TRUE),
    outside_2sd = sum(z_score > 2, na.rm = TRUE),
    outside_3sd = sum(z_score > 3, na.rm = TRUE),
    pct_outside_1sd = outside_1sd / total_individuals * 100,
    pct_outside_2sd = outside_2sd / total_individuals * 100,
    pct_outside_3sd = outside_3sd / total_individuals * 100
  )

outliers_1sd_individuals <- merged_df[merged_df$z_score > 1, ]
outliers_2sd_individuals <- merged_df[merged_df$z_score > 2, ]
outliers_3sd_individuals <- merged_df[merged_df$z_score > 3, ]

#adding wild individuals to the plot
roh.df.merge1 <- ind_merged$IID #populate ID's
colnames(roh.df.merge1)[1] <- "IID"
roh.df.merge1 <- as.data.frame(roh.df.merge1) #create the data frame
roh.df.merge1$thin2_merged <- ((ind_merged$KB*1000)/2425730029)*100
write.csv(roh.df.merge, "merged_roh_thin2.csv")
summary(thin2_merged$KB)
roh.df.merge1 <- merge(roh.df.merge1, w_pop_id_df, by = "IID", all = TRUE)
colnames(roh.df.merge1)[2] <- "Thin_test2"
colnames(roh.df.merge1)[3] <- "Cat.Group"
lepa_roh_df <- bind_rows(
  roh.df.merge1 %>% select(IID, Thin_test2, Cat.Group),  # adjust column names as needed
  roh_zoo %>% select(IID, Thin_test2, Cat.Group)
)
write.csv(lepa_roh_df, "supplemental_table_roh.csv")
#FROH violin plot with both wild and zoo present
lepa_roh_df <- read.csv("supplemental_table_roh.csv")
ggplot() +
  # violin plot
  geom_violin(data = lepa_roh_df, aes(x = Cat.Group, y = Thin_test2, fill = Cat.Group), 
              alpha = 0.7) +
  # Add individual points
  geom_jitter(data = lepa_roh_df, aes(x = Cat.Group, y = Thin_test2), 
              width = 0.1, alpha = 0.4, size = 1) +
  # Add population means with error bars
  geom_point(data = origin_mean_roh, aes(x = Cat.Group, y = Average), 
             color = "black", size = 4, shape = 18) +
  geom_errorbar(data = origin_mean_roh, 
                aes(x = Cat.Group, y = Average, 
                    ymin = pmax(0, Average - StdDev), 
                    ymax = Average + StdDev),
                color = "black", width = 0.2, size = 1) +
  scale_fill_manual(values = c("Ranch" = "#ffb2b0", "Refuge" = "#01004c",
                               "Generic" = "#D66857", "Brazilian" = "#3B967f")) +
  # Labels and theme
  labs(x = "Origin", y = expression(F[ROH])) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")
#FROH violin plot with just zoo
ggplot() +
  # violin plot
  geom_violin(data = roh_zoo, aes(x = Cat.Group, y = Thin_test2, fill = Cat.Group), alpha = 0.7) +
  # Add individual points
  geom_jitter(data = roh_zoo, aes(x = Cat.Group, y = Thin_test2),
              width = 0.1, alpha = 0.4, size = 1) +
  # Add population means with error bars
  geom_point(data = origin_mean_roh, aes(x = Cat.Group, y = Average),
             color = "black", size = 4, shape = 18) +
  geom_errorbar(data = origin_mean_roh, 
                aes(x = Cat.Group, y = Average, 
                    ymin = pmax(0, Average - 2*StdDev), ymax = Average + 2*StdDev),
                color = "black", width = 0.2, size = 1) +
  scale_fill_manual(values = c("Generic" = "#D66857", "Brazilian" = "#3B967f")) +
  # Labels and theme
  labs(x = "Origin", y = expression(F[ROH])) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")
#zoo descriptive stats
z_roh_t2 <- read.table("zoo_roh_thin_test2.hom", header = T)
z_roh_descript <- merge(z_roh_t2, z_pop_id, by = "IID", all = TRUE)
z_roh_descript$Length_MB <- z_roh_descript$KB/1000
z_roh_descript$Category <- cut(z_roh_descript$Length_MB,
                               breaks = c(0, 1, 2, 4, 6, 8, Inf),
                               labels = c("<1Mb", "1-2Mb", "2-4Mb", "4-6Mb", "6-8Mb", ">8Mb"),
                               include.lowest = TRUE)
z_population_category_counts <- z_roh_descript %>%
  group_by(Cat.Group, Category)  %>%
  summarise(Count = n(), Total_Length_MB = sum(Length_MB), .groups = "drop")
write.csv(z_population_category_counts, "zoo_roh_categories.csv")
write.csv(origin_mean_roh, "zoo_mean_roh.csv")
summary(z_roh_descript$Length_MB)

####proportion of ROH lengths
#read in .hom files
w_roh_t2 <- read.table("wild_roh_thin_test2.hom", header = T) #thin test 2
w_roh_ct <- read.table("wild_thin_composite1.hom", header = T) #composite test 1
#adding population information
colnames(w_pop_id_df)[1] <- "IID"
w_roh_t2 <- merge(w_roh_t2, w_pop_id_df, by = "IID", all = TRUE)
#calculate length in Mb -- kb to mb conversion
w_roh_t2$Length_MB <- w_roh_t2$KB/1000
w_roh_ct$Length_MB <- w_roh_ct$KB/1000
#look at resulting distribution
dev.off()
hist(w_roh_t2$Length_MB, main="Distribution of ROH lengths -- T2", xlab="Length (MB)")
hist(w_roh_ct$Length_MB, main = "Distribution of ROH lengths -- CT", xlab = "Length (MB)")
#define length categories
w_roh_t2$Category <- cut(w_roh_t2$Length_MB,
                         breaks = c(0, 1, 2, 4, 6, 8, Inf),
                         labels = c("<1Mb", "1-2Mb", "2-4Mb", "4-6Mb", "6-8Mb", ">8Mb"),
                         include.lowest = TRUE)
w_roh_ct$Category <- cut(w_roh_ct$Length_MB,
                         breaks = c(0, 1, 2, 4, 6, 8, Inf),
                         labels = c("<1Mb", "1-2Mb", "2-4Mb", "4-6Mb", "6-8Mb", ">8Mb"),
                         include.lowest = TRUE)
#get total proportion -- cumulative
w_total_length_all <- sum(w_roh_t2$Length_MB)
w_total_summary <- w_roh_t2 %>%
  group_by(Category) %>%
  summarize(
    Count = n(),
    Total_Length_MB = sum(Length_MB),
    Proportion_Length = sum(Length_MB)/w_total_length_all,
    Proportion_count = n()/nrow(w_roh_t2)
  )
w_total_length_all_ct <- sum(w_roh_ct$Length_MB)
w_total_summary_ct <- w_roh_ct %>%
  group_by(Category) %>%
  summarise(
    Count = n(),
    Total_Length_MB = (sum(Length_MB)),
    Proportion_Length = sum(Length_MB)/w_total_length_all_ct,
    Proportion_count = n()/nrow(w_roh_ct)
  )
#save results
write.csv(w_total_summary, "overall_roh_length_proportions_t2.csv")
write.csv(w_total_summary_ct, "overall_roh_length_proportions_ct1.csv")

#summary by individual

#summary by population 
population_category_counts <- w_roh_t2 %>%
  group_by(pop_id, Category)  %>%
  summarise(Count = n(), Total_Length_MB = sum(Length_MB), .groups = "drop")
#proportion by population
population_category_proportions <- population_category_counts %>%
  group_by(pop_id) %>%
  mutate(Proportion_Count = Count/sum(Count),
         Proportion_Length = Total_Length_MB/sum(Total_Length_MB))

#normalizing data for direct comparison
individuals_per_pop <- data.frame(
  pop_id = c("Ranch", "Refuge"),
  n_individuals = c(26, 18)
)

normalized_roh <- population_category_counts %>%
  left_join(individuals_per_pop, by = "pop_id") %>%
  group_by(pop_id) %>%
  mutate(
    # Within-population proportions
    Proportion_Count = Count / sum(Count),
    Proportion_Length = Total_Length_MB / sum(Total_Length_MB),
    
    # Per-individual normalized metrics
    Avg_ROH_Count_Per_Individual = Count / n_individuals,
    Avg_ROH_Length_MB_Per_Individual = Total_Length_MB / n_individuals
  ) %>%
  ungroup()
write.csv(normalized_roh, "pop_normalized_roh.csv")
write.csv(population_mean_roh, "pop_mean_froh.csv")
write.csv(w_total_summary, "roh_summary.csv")
####Visualizing ROH --karyotype plot -- wild####

#read in hom file
w_hom <- read.table("wild_thin_composite2.hom", header = TRUE)
#make data frame for plotting
w_roh_plot <- data.frame(
  chr = paste0(w_hom$CHR),  # Add 'chr' prefix if needed
  start = w_hom$POS1,
  end = w_hom$POS2,
  kb = w_hom$KB,           # ROH length in KB
  nsnp = w_hom$NSNP,       # Number of SNPs in ROH
  sample_id = w_hom$IID    # Individual ID
)
w_pop_id <- read.csv("wild_origins_2.csv") #adding population information
w_pop_id_df <-as.data.frame(w_pop_id)
w_roh_pop <- merge(w_roh_plot, w_pop_id_df, by = "sample_id", all = TRUE) #merging the population information with the data frame
#fixing chrome naming
chr_map <- c(
  "NC_058368.1" = "chr1",
  "NC_058369.1" = "chr2",
  "NC_058370.1" = "chr3",
  "NC_058371.1" = "chr4",
  "NC_058372.1" = "chr5",
  "NC_058373.1" = "chr6",
  "NC_058374.1" = "chr7",
  "NC_058375.1" = "chr8",
  "NC_058376.1" = "chr9",
  "NC_058377.1" = "chr10",
  "NC_058378.1" = "chr11",
  "NC_058379.1" = "chr12",
  "NC_058380.1" = "chr13",
  "NC_058381.1" = "chr14",
  "NC_058382.1" = "chr15",
  "NC_058383.1" = "chr16",
  "NC_058384.1" = "chr17",
  "NC_058385.1" = "chr18"
)
w_roh_plot$chr <- chr_map[w_roh_plot$chr]
#create custom feline genotype for karyoploteR, make a custom plot type for the 18 chr
feline_chr_sizes <- data.frame(
  chr = c(paste0("chr", 1:18)), 
  start = rep(1, 18),
  end = c(239367248,  # chr1
          169388855,  # chr2
          140443288,  # chr3
          205367284,  # chr4
          151959158,  # chr5
          148491486,  # chr6
          142168536,  # chr7
          221611373,  # chr8
          158578461,  # chr9
          115366950,  # chr10
          88083857,   # chr11
          94435393,   # chr12
          95154158,   # chr13
          61876196,   # chr14
          61988844,   # chr15
          41437797,   # chr16
          69239673,   # chr17
          83466477   # chr18
  )
)  
feline_genome_gr <- GRanges(seqnames = feline_chr_sizes$chr, ranges = IRanges(start = feline_chr_sizes$start, end = feline_chr_sizes$end)) #make into grange object
ocel_genome <- feline_genome_gr #making the custom plot type
seqlevels(ocel_genome) <- feline_chr_sizes$chr
seqlengths(ocel_genome) <- feline_chr_sizes$end

# roh convert to GRanges object
w_roh_gr <- GRanges(
  seqnames = w_roh_plot$chr,
  ranges = IRanges(start = w_roh_plot$start, end = w_roh_plot$end),
  kb = w_roh_plot$kb,
  nsnp = w_roh_plot$nsnp,
  sample_id = w_roh_plot$sample_id) #making hom file into grange

pop_gr <- GRanges(
  seqnames = w_roh_pop$chr,
  ranges = IRanges(start = w_roh_pop$start, end = w_roh_pop$end),
  kb = w_roh_pop$kb,
  nsnp = w_roh_pop$nsnp,
  sample_id = w_roh_pop$sample_id,
  population_id = w_roh_pop$pop_id)

#####creating function for plotting -- individual
plot_individual_roh <- function(sample_id, output_file = NULL) {
  # Filter ROH data for specific individual
  individual_roh <- w_roh_gr[mcols(w_roh_gr)$sample_id %in% sample_id]
  # Determine if output should go to a file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  #create the plot using the custom genome
  w_kp <- plotKaryotype(genome = ocel_genome, plot.type = 2, main = paste("ROH for", sample_id))
  #kpAddChromosomeNames(w_kp, srt = 45, cex = 0.8) #not using this line for now
  #plot individuals
  kpRect(w_kp, 
         chr = as.character(seqnames(individual_roh)), 
         x0 = start(individual_roh), 
         x1 = end(individual_roh),
         y0 = 0, 
         y1 = 1, 
         col = "#FF000080",  # Semi-transparent red
         border = "#FF0000",
         r0 = 0.5, r1 = 0.8)
  # Close file if opened
  if (!is.null(output_file)) {
    dev.off()
  }
  # Return the filtered data
  return(individual_roh)
}

####function for plotting -- population wide
plot_population_roh_smoothed <- function(population_id, output_file = NULL, window_size = 1e5, smooth = TRUE) {
  # Subset to individuals from the desired population
  pop_roh <- pop_gr[mcols(pop_gr)$population_id == population_id]
  
  if (length(pop_roh) == 0) {
    message("No ROH data found for population: ", population_id)
    return(NULL)
  }
  
  # Optional output to file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  
  # Plot karyotype
  pop_kp <- plotKaryotype(genome = ocel_genome, plot.type = 2, main = paste("ROH Density for", population_id))
  
  # Calculate raw coverage (how many ROHs overlap each base)
  roh_cov <- coverage(pop_roh)
  
  # Smooth and plot each chromosome
  for (chr in names(roh_cov)) {
    cov_vector <- as.numeric(roh_cov[[chr]])
    
    if (smooth) {
      # Apply rolling mean smoothing (simple moving average)
      kernel <- rep(1/window_size, window_size)
      smoothed <- stats::filter(cov_vector, kernel, sides = 2, circular = FALSE)
    } else {
      smoothed <- cov_vector
    }
    
    # Build GRanges object with smoothed values
    pos <- IRanges(start = seq_along(smoothed), width = 1)
    smoothed_gr <- GRanges(seqnames = chr, ranges = pos, score = smoothed)
    
    # Reduce resolution for plotting (optional)
    smoothed_gr <- smoothed_gr[start(smoothed_gr) %% 1000 == 0]  # sample every 1kb
    
    # Plot line track
    kpLines(pop_kp, chr = as.character(seqnames(smoothed_gr)), x = start(smoothed_gr), y = mcols(smoothed_gr)$score, col = "red", r0 = 0.5, r1 = 0.8)
  }
  
  if (!is.null(output_file)) {
    dev.off()
  }
  
  invisible(NULL)
}

##using the function -- plotting individual roh
unique_samples <- unique(w_roh_plot$sample_id)
print(paste("Found", length(unique_samples), "samples in the data"))

##using the function -- plotting population roh
unique_pops <- unique(mcols(pop_gr)$population_id)


# Sanitize function to make safe filenames
sanitize_filename <- function(name) {
  gsub("[^A-Za-z0-9_]", "_", name)  # Replace anything that's not a letter, number, or underscore
}

# Create directory for plots
dir.create("wild_roh_plots", showWarnings = FALSE)

dir.create("wild_roh_density_plots", showWarnings = FALSE)

# Loop through each sample and create sanitized output files
for (sample_id in unique_samples) {
  safe_id <- sanitize_filename(sample_id)
  output_file <- paste0("wild_roh_plots/", safe_id, "_roh_plot.pdf")
  plot_individual_roh(sample_id, output_file)
}

#loop through to create output files
for (pop in unique_pops) {
  safe_name <- sanitize_filename(pop)
  output_file <- paste0("wild_roh_density_plots/pop_", safe_name, "_roh_density.pdf")
  plot_population_roh_smoothed(population_id = pop, output_file = output_file)
}
dev.off()
###ROH overlap plot


####merged roh test -- post processing completed in command line####
##in command line, stitched roh that were separated by 1010 bases, resulting in:
#Original number of segments: 36348
#Merged number of segments: 33323
#Reduction: 3025 segments (8.32%)

#read in data
thin2_merged <- read.table("wild_thintest2_merged.hom", header = TRUE)
ind_merged <- read.table("wild_thintest2_merged.ind", header = TRUE)

#create data frame for analysis
roh.df.merge1 <- ind_merged$IID #populate ID's
colnames(roh.df.merge1)[1] <- "IID"
roh.df.merge1 <- as.data.frame(roh.df.merge1) #create the data frame
roh.df.merge1$thin2_merged <- ((ind_merged$KB*1000)/2425730029)*100
write.csv(roh.df.merge, "merged_roh_thin2.csv")
summary(thin2_merged$KB)
roh.df.merge1 <- merge(roh.df.merge1, w_pop_id_df, by = "IID", all = TRUE)

#average % genome in roh by population
colnames(w_pop_id_df)[1] <- "IID"
roh.df.merge.pop <- merge(roh.df.merge, w_pop_id_df, by = "IID", all = TRUE)
population_mean_roh <- roh.df.merge.pop %>%
  group_by(pop_id) %>%
  summarise(Average = mean(thin2_merged, na.rm = TRUE),
            Min = min(thin2_merged),
            Max = max(thin2_merged),
            Count = n(),
            StdDev = sd(thin2_merged, na.rm = TRUE),
            StdError = (sd(thin2_merged, na.rm = TRUE)/sqrt(n())))
write.csv(population_mean_roh, "pop_mean_froh.csv")

#calc bins
thin2_merged$Length_MB <- thin2_merged$KB/1000
thin2_merged$Category <- cut(thin2_merged$Length_MB,
                             breaks = c(0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, Inf),
                             labels = c("<1Mb", "1-2Mb", "2-4Mb", "4-6Mb", "6-8Mb", "8-10MB", "10-12Mb", "12-14Mb", 
                                        "14-16Mb", "16-18Mb", "18-20Mb", ">20Mb"),
                             include.lowest = TRUE)
#save new hom file
write.csv(thin2_merged, "wild_thintest2_merged_categories_hom.csv")
#summaries by individual
merged_ind_summary <- thin2_merged %>%
  group_by(IID) %>%
  summarize(
    Total_ROH_Count = n(),
    Total_ROH_Length_MB = sum(Length_MB),
    Mean_ROH_Length_MB = mean(Length_MB),
    Median_ROH_Length_MB = median(Length_MB),
    Min_ROH_Length_MB = min(Length_MB),
    Max_ROH_Length_MB = max(Length_MB)
  )
#proportions by individual
merged_ind_cat_summary <- thin2_merged %>%
  group_by(IID, Category) %>%
  summarize(
    Count = n(),
    Total_Length_MB = sum(Length_MB),
    .groups = "drop"
  ) %>%
  group_by(IID) %>%
  mutate(
    Individual_Total_Length = sum(Total_Length_MB),
    Proportion_Length = Total_Length_MB / Individual_Total_Length,
    Proportion_Count = Count / sum(Count)
  )
write.csv(merged_ind_cat_summary, "merged_individual_categories_roh.csv")
##proportions by populations
w_pop_id <- read.csv("wild_origins_2.csv") #adding population information
w_pop_id_df <-as.data.frame(w_pop_id) #making df
colnames(w_pop_id_df)[1] <- "IID"
thin2_merged <- merge(thin2_merged, w_pop_id_df, by = "IID", all = TRUE)

#pure counts, non normalized
population_category_counts <- thin2_merged %>%
  group_by(pop_id, Category)  %>%
  summarise(Count = n(), Total_Length_MB = sum(Length_MB), .groups = "drop")
#normalized for sample size differences
individuals_per_pop <- data.frame(
  pop_id = c("Ranch", "Refuge"),
  n_individuals = c(26, 18)
)

normalized_roh_merged <- population_category_counts %>%
  left_join(individuals_per_pop, by = "pop_id") %>%
  group_by(pop_id) %>%
  mutate(
    # Within-population proportions
    Proportion_Count = Count / sum(Count),
    Proportion_Length = Total_Length_MB / sum(Total_Length_MB),
    
    # Per-individual normalized metrics
    Avg_ROH_Count_Per_Individual = Count / n_individuals,
    Avg_ROH_Length_MB_Per_Individual = Total_Length_MB / n_individuals
  ) %>%
  ungroup()
write.csv(normalized_roh_merged, "roh_population_categories_merged.csv")

##plotting new merged roh bins
#read in hom file
merg_roh <- read.table("wild_thintest2_merged.hom", header = TRUE)
#make data frame for plotting
merged_plot <- data.frame(
  chr = paste0(merg_roh$CHR),  # Add 'chr' prefix if needed
  start = merg_roh$BP1,
  end = merg_roh$BP2,
  kb = merg_roh$KB,           # ROH length in KB
  nsnp = merg_roh$NSNP,       # Number of SNPs in ROH
  sample_id = merg_roh$IID    # Individual ID
)

#fixing chrome naming
chr_map <- c(
  "NC_058368.1" = "chr1",
  "NC_058369.1" = "chr2",
  "NC_058370.1" = "chr3",
  "NC_058371.1" = "chr4",
  "NC_058372.1" = "chr5",
  "NC_058373.1" = "chr6",
  "NC_058374.1" = "chr7",
  "NC_058375.1" = "chr8",
  "NC_058376.1" = "chr9",
  "NC_058377.1" = "chr10",
  "NC_058378.1" = "chr11",
  "NC_058379.1" = "chr12",
  "NC_058380.1" = "chr13",
  "NC_058381.1" = "chr14",
  "NC_058382.1" = "chr15",
  "NC_058383.1" = "chr16",
  "NC_058384.1" = "chr17",
  "NC_058385.1" = "chr18"
)
merged_plot$chr <- chr_map[merged_plot$chr]
#create custom feline genotype for karyoploteR, make a custom plot type for the 18 chr
feline_chr_sizes <- data.frame(
  chr = c(paste0("chr", 1:18)), 
  start = rep(1, 18),
  end = c(239367248,  # chr1
          169388855,  # chr2
          140443288,  # chr3
          205367284,  # chr4
          151959158,  # chr5
          148491486,  # chr6
          142168536,  # chr7
          221611373,  # chr8
          158578461,  # chr9
          115366950,  # chr10
          88083857,   # chr11
          94435393,   # chr12
          95154158,   # chr13
          61876196,   # chr14
          61988844,   # chr15
          41437797,   # chr16
          69239673,   # chr17
          83466477   # chr18
  )
)  
feline_genome_gr <- GRanges(seqnames = feline_chr_sizes$chr, ranges = IRanges(start = feline_chr_sizes$start, end = feline_chr_sizes$end)) #make into grange object
ocel_genome <- feline_genome_gr #making the custom plot type
seqlevels(ocel_genome) <- feline_chr_sizes$chr
seqlengths(ocel_genome) <- feline_chr_sizes$end

# roh convert to GRanges object
merg_gr <- GRanges(
  seqnames = merged_plot$chr,
  ranges = IRanges(start = merged_plot$start, end = merged_plot$end),
  kb = merged_plot$kb,
  nsnp = merged_plot$nsnp,
  sample_id = merged_plot$sample_id) #making hom file into grange


#####creating function for plotting -- individual
plot_individual_roh <- function(sample_id, output_file = NULL) {
  # Filter ROH data for specific individual
  individual_roh <- merg_gr[mcols(merg_gr)$sample_id %in% sample_id]
  # Determine if output should go to a file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  #create the plot using the custom genome
  w_kp <- plotKaryotype(genome = ocel_genome, plot.type = 2, main = paste("ROH for", sample_id))
  #kpAddChromosomeNames(w_kp, srt = 45, cex = 0.8) #not using this line for now
  #plot individuals
  kpRect(w_kp, 
         chr = as.character(seqnames(individual_roh)), 
         x0 = start(individual_roh), 
         x1 = end(individual_roh),
         y0 = 0, 
         y1 = 1, 
         col = "#FF000080",  # Semi-transparent red
         border = "#FF0000",
         r0 = 0.5, r1 = 0.8)
  # Close file if opened
  if (!is.null(output_file)) {
    dev.off()
  }
  # Return the filtered data
  return(individual_roh)
}

#using the function
unique_samples <- unique(merged_plot$sample_id)
sanitize_filename <- function(name) {
  gsub("[^A-Za-z0-9_]", "_", name)  # Replace anything that's not a letter, number, or underscore
}
dir.create("merg_roh_plot", showWarnings = FALSE)
for (sample_id in unique_samples) {
  safe_id <- sanitize_filename(sample_id)
  output_file <- paste0("merg_roh_plot/", safe_id, "_roh_plot.pdf")
  plot_individual_roh(sample_id, output_file)
}


###plotting chrom 10 only for three individuals
#read in hom file
merg_roh <- read.table("wild_thintest2_merged.hom", header = TRUE)
#make data frame for plotting
merged_plot <- data.frame(
  chr = paste0(merg_roh$CHR),  # Add 'chr' prefix if needed
  start = merg_roh$BP1,
  end = merg_roh$BP2,
  kb = merg_roh$KB,           # ROH length in KB
  nsnp = merg_roh$NSNP,       # Number of SNPs in ROH
  sample_id = merg_roh$IID    # Individual ID
)
#fixing chrome naming
chr_map <- c(
  "NC_058368.1" = "chr1",
  "NC_058369.1" = "chr2",
  "NC_058370.1" = "chr3",
  "NC_058371.1" = "chr4",
  "NC_058372.1" = "chr5",
  "NC_058373.1" = "chr6",
  "NC_058374.1" = "chr7",
  "NC_058375.1" = "chr8",
  "NC_058376.1" = "chr9",
  "NC_058377.1" = "chr10",
  "NC_058378.1" = "chr11",
  "NC_058379.1" = "chr12",
  "NC_058380.1" = "chr13",
  "NC_058381.1" = "chr14",
  "NC_058382.1" = "chr15",
  "NC_058383.1" = "chr16",
  "NC_058384.1" = "chr17",
  "NC_058385.1" = "chr18"
)
merged_plot$chr <- chr_map[merged_plot$chr]
#filtering dataframe for only chrome 10
merged_plot <- merged_plot[merged_plot$chr == "chr10", ]
#creating a custom genotype for karyoplotr for chrom 10 only
feline_chr_sizes <- data.frame(
  chr = "chr10", 
  start = 1,
  end = 115366950   # chr10 length only
) 
feline_genome_gr <- GRanges(seqnames = feline_chr_sizes$chr, ranges = IRanges(start = feline_chr_sizes$start, end = feline_chr_sizes$end)) #make into grange object
ocel_genome <- feline_genome_gr #making the custom plot type
seqlevels(ocel_genome) <- feline_chr_sizes$chr
seqlengths(ocel_genome) <- feline_chr_sizes$end
#making hom file into a granges file
merg_gr <- GRanges(
  seqnames = merged_plot$chr,
  ranges = IRanges(start = merged_plot$start, end = merged_plot$end),
  kb = merged_plot$kb,
  nsnp = merged_plot$nsnp,
  sample_id = merged_plot$sample_id)
#creating the function to plot
chr_10_plot <- function(sample_id, output_file = NULL) {
  # Filter ROH data for specific individual
  individual_roh <- merg_gr[mcols(merg_gr)$sample_id %in% sample_id]
  
  # Check if individual has any ROH on chr10
  if (length(individual_roh) == 0) {
    message(paste("No ROH found on chromosome 10 for", sample_id))
    return(NULL)
  }
  
  # Determine if output should go to a file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  #create the plot using the custom genome
  w_kp <- plotKaryotype(genome = ocel_genome, plot.type = 1, main = paste("ROH on Chromosome 10 for", sample_id))
  #kpAddChromosomeNames(w_kp, srt = 45, cex = 0.8) #not using this line for now
  #plot individuals
  kpRect(w_kp, 
         chr = as.character(seqnames(individual_roh)), 
         x0 = start(individual_roh), 
         x1 = end(individual_roh),
         y0 = 0, 
         y1 = 1, 
         col = "#FF000080",  # Semi-transparent red
         border = "#FF0000",
         r0 = 0.05, r1 = 0.95,
         data.panel = "ideogram") #should plot directly onto the ideogram
  # Close file if opened
  if (!is.null(output_file)) {
    dev.off()
  }
  # Return the filtered data
  return(individual_roh)
}
#using the function
unique_samples <- unique(merged_plot$sample_id)
sanitize_filename <- function(name) {
  gsub("[^A-Za-z0-9_]", "_", name)  # Replace anything that's not a letter, number, or underscore
}
dir.create("merg_roh_plot_chr10", showWarnings = FALSE)
for (sample_id in unique_samples) {
  safe_id <- sanitize_filename(sample_id)
  output_file <- paste0("merg_roh_plot_chr10/", safe_id, "_chr10_roh_plot.pdf")
  chr_10_plot(sample_id, output_file)
}


################################################################################
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Admixture_outputs")
##admixture plots
#.q files from ADMIXTURE performed on a mac computer
####creating admixture plots using ggplot -- wild####
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

####ranch only####
r_k2_table <- read.table("Ranch_standard_final.2.q")
r_ids <- read.table("ranch_subset.txt")
r_ids = subset(r_ids, select = -c(V1))
r_ids = rename(r_ids, sample_id = V2, pop_id = V3)
r_k2 <- cbind(r_ids, r_k2_table)
r_k2_long <- pivot_longer(
  r_k2, cols = c(V1, V2),
  names_to = "ancestry",
  values_to = "proportion"
) #pivots the table to the proportions are able to be plotted as stacked bars\
r_k2_long$pop_id <- str_to_upper(r_k2_long$pop_id)
#plot code
ranch_k2plot <-
  ggplot(r_k2_long, aes(x = factor(sample_id), y = proportion, fill = ancestry)) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(pop_id), switch = "both", scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 8),
        strip.text.x = element_text(face = "bold", size = 12),
        strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("V1" = "#3D1308", "V2" = "#ffb2b0")) +
  guides(fill = "none")
ranch_k2plot
ggsave("wild_k2_admixture_altcol.png", w_k2plot, width = 15, height = 8, bg = "white")

################################################################################
#fst prep for vcftools
###populations need to be filtered separately, and then merged
##No MAF and careful with HWE in filtering
#wild splitting and filtering
#dispersers are filtered with individuals at their current population, not natal population
system("plink --bfile Wild_BaseFilter_AllChrom --keep ranch_fst_subset.txt --allow-extra-chr --chr-set 18 --make-bed --out Ranch_base_Fst")
system("plink --bfile Wild_BaseFilter_AllChrom --keep refuge_fst_subset.txt --allow-extra-chr --chr-set 18 --make-bed --out Refuge_base_Fst")
system("plink --bfile Zoo_BaseFilter_Allchrom --keep generic_fst_subset.txt --allow-extra-chr --chr-set 18 --make-bed --out Generic_base_Fst")
system("plink --bfile Zoo_BaseFilter_Allchrom --keep brazilian_fst_subset.txt --allow-extra-chr --chr-set 18 --make-bed --out Brazilian_base_Fst")

system("plink --bfile Ranch_base_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Ranch_geno1_hwe_Fst")
system("plink --bfile Ranch_geno1_hwe_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Ranch_fst_filter")

system("plink --bfile Refuge_base_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Refuge_geno1_hwe_Fst")
system("plink --bfile Refuge_geno1_hwe_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Refuge_fst_filter")

system("plink --bfile Generic_base_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out generic_geno1_hwe_Fst")
system("plink --bfile generic_geno1_hwe_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Generic_fst_filter")

system("plink --bfile Brazilian_base_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Brazilian_geno1_hwe_Fst")
system("plink --bfile Brazilian_geno1_hwe_Fst --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Brazilian_fst_filter")
#zoo splitting and filtering
#exporting as vcf's
system("plink2 --bfile Ranch_fst_filter --chr-set 18 --allow-extra-chr --export vcf-4.2 bgz --out ranch_fst_in42")
system("plink2 --bfile Refuge_fst_filter --chr-set 18 --allow-extra-chr --export vcf-4.2 bgz --out refuge_fst_in42")
system("plink2 --bfile Generic_fst_filter --chr-set 18 --allow-extra-chr --export vcf-4.2 bgz --out generic_fst_in42")
system("plink2 --bfile Brazilian_fst_filter --chr-set 18 --allow-extra-chr --export vcf-4.2 bgz --out brazilian_fst_in42")
system("plink2 --bfile Zoo_refiltered --chr-set 18 --allow-extra-chr --export vcf-4.2 bgz --out zoo_refiltered_42")
system("plink2 --bfile Wild_refiltered --chr-set 18 --allow-extra-chr --export vcf-4.2 bgz --out wild_refiltered_42")

####Fst and nucleotide diversity####
summary_pi<- read.csv("summary_pi_per_individual.csv", header = TRUE)
#nucleotide diversitt summary -- reading in data
nucleotide_data_clean <- summary_pi %>%
  mutate(individual = gsub("-.*", "", individual),
         individual = gsub("standard_final_24040DeY_", "", individual))
origins <- as.data.frame(read.csv("lepa_origins.csv"))
as.data.frame(nucleotide_data_clean)
origins_clean <- origins %>%
  mutate(individual = gsub("-.*", "", individual))

nucleo_by_pop <- nucleotide_data_clean %>%
  left_join(origins_clean, by = "individual")

#making summary stats
n_stats <- nucleo_by_pop %>%
  group_by(Pop) %>%
  summarise(
    mean_pi_avg = mean(mean_pi, na.rm = TRUE),
    sd_pi_avg = sd(mean_pi, na.rm = TRUE),
    n_individuals = n(),
    .groups = 'drop'
  )
#plotting nucleotide diversity -- wild only
ggplot() +
  # violin plot
  geom_violin(data = subset(nucleo_by_pop, Pop %in% c("Ranch", "Refuge")), 
              aes(x = Pop, y = mean_pi, fill = Pop), 
              alpha = 0.7) +
  # Add individual points
  geom_jitter(data = subset(nucleo_by_pop, Pop %in% c("Ranch", "Refuge")), 
              aes(x = Pop, y = mean_pi), 
              width = 0.1, alpha = 0.4, size = 3) +
  # Add population means with error bars
  geom_point(data = subset(n_stats, Pop %in% c("Ranch", "Refuge")), aes(x = Pop, y = mean_pi_avg), 
             color = "black", size = 4, shape = 18) +
  geom_errorbar(data = subset(n_stats, Pop %in% c("Ranch", "Refuge")), 
                aes(x = Pop, y = mean_pi_avg, 
                    ymin = mean_pi_avg - 2*sd_pi_avg, ymax = mean_pi_avg + 2*sd_pi_avg),
                color = "black", width = 0.2, size = 1) +
  scale_fill_manual(values = c("Ranch" = "#ffb2b0", "Refuge" = "#01004c")) +
  # Labels and theme
  labs(x = "Population", y = "Nucleotide Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")

###nucleotide diversity all pops
ggplot() +
  # violin plot
  geom_violin(data = nucleo_by_pop, 
              aes(x = Pop, y = mean_pi, fill = Pop), 
              alpha = 0.7) +
  # Add individual points
  geom_jitter(data = nucleo_by_pop, aes(x = Pop, y = mean_pi), 
              width = 0.1, alpha = 0.6, size = 4) +
  # Add population means with error bars
  geom_point(data = n_stats, aes(x = Pop, y = mean_pi_avg), 
             color = "black", size = 4, shape = 18) +
  scale_fill_manual(values = c("Ranch" = "#ffb2b0", "Refuge" = "#01004c",
                               "Generic" = "#D66857", "Brazilian" = "#3B967f")) +
  # Labels and theme
  labs(x = "Population", y = "Nucleotide Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")


####expected and observed Heterozygosity####

# reading in data
Ra_vcf <- read.vcfR("ranch_final_thin.vcf.gz", verbose = TRUE)
ranch_genlight <- vcfR2genlight(Ra_vcf)

Re_vcf <- read.vcfR("refuge_final_thin.vcf.gz", verbose = TRUE)
Refuge_genlight <- vcfR2genlight(Re_vcf)

Ge_vcf <- read.vcfR("generic_final_thin.vcf.gz", verbose = TRUE)
Gen_genlight <- vcfR2genlight(Ge_vcf)

Br_vcf <- read.vcfR("brazilian_final_thin.vcf.gz", verbose = TRUE)
braz_genlight <- vcfR2genlight(Br_vcf)

##estimating values -- for loop
#create genlight list
gen_list <- list(
  "Ranch" = ranch_genlight,
  "Refuge" = Refuge_genlight,
  "Generic" = Gen_genlight,
  "Brazilian" = braz_genlight
)

#create dataframe for results
diversity_results <- data.frame(
  Population = character(),
  He = numeric(),
  He_SE = numeric(),
  Ho = numeric(),
  Ho_SE = numeric(),
  Fis = numeric(),
  Fis_SE = numeric(),
  stringsAsFactors = FALSE
)

#making for loop
for(i in 1:length(gen_list)) {
  pop_name <- names(gen_list)[i]
  gl_obj <- gen_list[[i]]
  #added to avoid error in the gl.report function
  gl_obj@other$loc.metrics.flags$monomorphs <- TRUE
  
  # Debug: Check the structure of het_stats
  # Uncomment the next line to see what's available:
  # print(str(gl.report.heterozygosity(gl_obj, method = "pop")))
  
  #calc heterozygosity stats - use method = "loc" to get per-locus values
  het_stats_pop <- gl.report.heterozygosity(gl_obj, method = "pop")
  het_stats_loc <- gl.report.heterozygosity(gl_obj, method = "loc")
  
  # Extract per-locus values for calculating standard error
  he_values <- het_stats_loc$He  # Expected heterozygosity per locus
  ho_values <- het_stats_loc$Ho  # Observed heterozygosity per locus
  
  # Calculate means
  he_mean <- mean(he_values, na.rm = TRUE)
  ho_mean <- mean(ho_values, na.rm = TRUE)
  
  # Calculate standard errors
  he_se <- sd(he_values, na.rm = TRUE) / sqrt(sum(!is.na(he_values)))
  ho_se <- sd(ho_values, na.rm = TRUE) / sqrt(sum(!is.na(ho_values)))
  
  # Calculate Fis per locus and its standard error
  fis_values <- 1 - (ho_values / he_values)
  fis_values <- fis_values[is.finite(fis_values)]  # Remove infinite/NaN values
  fis_mean <- mean(fis_values, na.rm = TRUE)
  fis_se <- sd(fis_values, na.rm = TRUE) / sqrt(sum(!is.na(fis_values)))
  
  #add to dataframe with standard errors
  diversity_results <- rbind(diversity_results, 
                             data.frame(Population = pop_name,
                                        He = he_mean,
                                        He_SE = he_se,
                                        Ho = ho_mean,
                                        Ho_SE = ho_se,
                                        Fis = fis_mean,
                                        Fis_SE = fis_se))
}
print(diversity_results)
write.csv(diversity_results, "diversity_statistics_by_pop.csv")

#diversity results by individuals
#creating data frame for results
individual_results <- data.frame(
  Individual = character(),
  Ho = numeric(),
  stringsAsFactors = FALSE
)

#making for loop
for(i in 1:length(gen_list)) {
  pop_name <- names(gen_list)[i]
  gl_obj <- gen_list[[i]]
  #added to avoid error in the gl.report function
  gl_obj@other$loc.metrics.flags$monomorphs <- TRUE
  #calc heterozygosity stats
  het_stats <- gl.report.heterozygosity(gl_obj, method = "ind")
  #extract values
  individual_ids <- het_stats$ind.name
  ho <- het_stats$Ho
  #add to dataframe
  individual_results <- rbind(individual_results, 
                              data.frame(Individual = individual_ids,
                                         Ho = ho))
}
print(individual_results)
write.csv(diversity_results, "diversity_statistics_by_pop.csv")