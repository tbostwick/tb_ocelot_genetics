##LEPA Analyses -- samba R package
#5/15/25
##Analyses of LEPA, Wild, and Zoo data sets using the samba R package with code provided from
  #the TWS2023 Genomic Analysis Workshop
##following the samba workshop, even if not all for outlier loci
##Analyses run on thinned data, with the following settings MAF: 0.05, MISS: 90%, Biallelic, HWE
##################################################################################################################
####Package set up and directory####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Samba_analyses")

## Load Packages ##
install.packages("devtools")
install.packages("tidyverse")
install.packages("hierfstat")
install.packages("gplots")
install.packages("pcadapt")
install.packages("Rtools")
BiocManager::install("LEA")
library(devtools)
install_github("drveera/ggman")
install_github("jdstorey/qvalue")
devtools::install_github("pievos101/PopGenome")
install.packages("VennDiagram")
library(LEA)
library(vcfR)
library(adegenet)
library(tidyverse)
library(hierfstat)
library(ggman)
library(ggplot2)
library(gplots)
library(pcadapt)
library(qvalue)
library(PopGenome)
library(VennDiagram)

#installing and loading samba
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt")
getpackages()
###############################################################################################################
#for these steps, samba requires you to filter in samba for the file structures to be correct
####loading in data as genlights -- unfiltered####
l_vcf <- read.vcfR("lepa_samba_nofilt_thin_final.vcf.gz", verbose = TRUE)
lepa_genlight <- vcfR2genlight(l_vcf)
lepa_genlight
##adding population assignments
l_pop <- read.csv("lepa_origins.csv")
l_df <- as.data.frame(l_pop) #creating dataframe from origins
l_df$ID <- paste0("24040DeY_", l_df$ID) #adding the prefix to the sample names

pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table
pop(lepa_genlight) <- pop_lookup[indNames(lepa_genlight)]# Directly assign populations based on IDs
pop(lepa_genlight) <- factor(pop(lepa_genlight), levels = levels(pop(lepa_genlight)))# Convert to factor with original levels
####initial look at data using samba####
#samba initialization and generation of quality control plots
genlight2sambar("lepa_genlight",do_confirm=TRUE)

#filtering
filterdata(snpmiss=0.1,min_mac=8,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=500, nchroms=18, silent=TRUE,maxprop_hefilter = 0.06)

#calculate kinship
calckinship()

#genetic diversity measurements
div_results <- calcdiversity(do_venn = FALSE) 
##################################################################################################################
#data can be pre-filtered outside of samba for the following steps
####loading packages and working directory####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Samba_analyses")
library(pcadapt)
library(vcfR)
library(qvalue)
library(VennDiagram)
library(grDevices)
library(dplyr)
####loading and prepping data####
#adding snp ids prior to analysis
system("plink2 --vcf lepa_dapc_thin.vcf.gz --chr-set 18 --allow-extra-chr --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out lepa_thinned_ids") #set variant ID -- shorter?
system("plink2 --bfile lepa_thinned_ids --chr-set 18 --allow-extra-chr --export vcf bgz --out lepa_thinned_ids")

system("plink2 --vcf zoo_dapc_thin.vcf.gz --chr-set 18 --allow-extra-chr --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out zoo_thinned_ids") #set variant ID -- shorter?
system("plink2 --bfile zoo_thinned_ids --chr-set 18 --allow-extra-chr --export vcf bgz --out zoo_thinned_ids")
#using the unique ID file generated earlier to avoid a merging error due to same variant names


z_vcf <- read.vcfR("zoo_thinned_ids.vcf.gz", verbose = TRUE)
zoo_genlight <- vcfR2genlight(z_vcf)
zoo_genlight

w_vcf <- read.vcfR("wild_dapc_thin.vcf.gz", verbose = TRUE)
wild_genlight <- vcfR2genlight(w_vcf)
wild_genlight

lepa_vcf <- read.vcfR("lepa_thinned_ids.vcf.gz", verbose = TRUE)
lepa_genlight <- vcfR2genlight(lepa_vcf)
lepa_genlight

##adding population assignments
w_pop <-read.csv("wild_origins.csv") #reading in origins
z_pop <- read.csv("zoo_origins.csv")
l_pop <- read.csv("lepa_origins.csv")
w_df <- as.data.frame(w_pop) #creating dataframe from origins
w_df$ID <- paste0("24040DeY_", w_df$ID) #adding the prefix to the sample names
z_df <-as.data.frame(z_pop)
z_df$ID <- paste0("24040DeY_", z_df$ID) #adding the prefix to the sample names
l_df <- as.data.frame(l_pop)
l_df$ID <- paste0("24040DeY_", l_df$ID) #adding the prefix to the sample names
##adding population id's to genlights
#wild
W_pop_lookup <- setNames(as.character(w_df$Pop), w_df$ID)# Create a lookup table
pop(wild_genlight) <- pop_lookup[indNames(wild_genlight)]# Directly assign populations based on IDs
pop(wild_genlight) <- factor(pop(wild_genlight), levels = levels(pop(wild_genlight)))# Convert to factor with original levels
#zoo
z_pop_lookup <- setNames(as.character(z_df$Pop), z_df$ID)# Create a lookup table
pop(zoo_genlight) <- pop_lookup[indNames(zoo_genlight)]# Directly assign populations based on IDs
pop(zoo_genlight) <- factor(pop(zoo_genlight), levels = levels(pop(zoo_genlight)))# Convert to factor with original levels
#lepa
l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table
pop(lepa_genlight) <- l_pop_lookup[indNames(lepa_genlight)]# Directly assign populations based on IDs
pop(lepa_genlight) <- factor(pop(lepa_genlight), levels = levels(pop(lepa_genlight)))# Convert to factor with original levels
####PCA -- wild####
w.pca.input <- read.pcadapt("wild_dapc_thin.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
wild_pca <- pcadapt(input = w.pca.input, K=20) #running pca with high k to check
plot(wild_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(wild_pca, option="screeplot", K=10, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(wild_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=wild_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(wild_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=wild_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(wild_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=wild_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(wild_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=wild_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(wild_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=wild_genlight$pop) #pc 5 and 6: structure yes?
PC67=plot(wild_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=wild_genlight$pop) #pc 6 and 7: structure, no not really
    #thus, retain 6 pc's

#calculate variance on each axis
EV <- (wild_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 18%, no sig decreases after 6 -- plateaus; so k=6

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
wild_k6 <- pcadapt(w.pca.input, K = 6)
summary(wild_k6)

####begin identifing outliers -- initial visualizations -- wild####

#basic manhattan plot -- just for visualizing the candidate loci
plot(wild_k6, option="manhattan")
#Q-Q Plots -- looking at how data is distributed
plot(wild_k6, option="qqplot")
#Histograms of the test statistic and of the p-values
hist(wild_k6$pvalues, xlab = "p-value",main = NULL, breaks = 50, col = "#2471A3")
    #plotting p-values and looking at significance, are majority of samples non-sig (.5 >)
    #samples with very low p-values are candidates for adaptive loci
plot(wild_k6, option = "stat.distribution")

####outlier tests -- wild####
## q-values  -- start looking for outlier loci
## for later analyses -- can exclude outlier to have a neutral and adaptive set of data for analyses
## do mult methods when doing outlier tests; several listed in samaba package
## can create a ven diagram to see where there is overlap; don't rely on just one method
##combined a better idea of what might be true outliers

#Choosing a cutoff for outlier detection
#To provide a list of outliers and choose a cutoff for outlier detection, there are several methods 
#that are listed below from the less conservative one to the more conservative one.
#q-value - least conservative
#bonferroni - most

#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs 
#for the geno3pops data, for an expected false discovery rate lower than 10%:
loci=wild_genlight$loc.names #calling locus names in genlight file
qval <- qvalue(wild_k6$pvalues)$qvalues #transformed p-values, making q-values column
alpha10 <- 0.1 #setting false discovery rate -- threshold; test at diff levels
outliers_qvalue10 <- which(qval < alpha10) #q values set at false discovery level of 10; looking for snps that have q value less than alpha
length(outliers_qvalue10) #How many snps detected?
# 69329 found
alpha5 <- 0.05 #5% more stat sig, but check at other levels
outliers_qvalue5.ints <- which(qval < alpha5)
length(outliers_qvalue5.ints)
# 46563 found;lost 20,000ish -- becoming more strict; can go low as want
outliers_qvalue5 <- as.character(loci)[outliers_qvalue5.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
outliers_qvalue5
write.csv(outliers_qvalue5, "adaptive_snps_q05.csv") #writing csv file
    ###snps found using q value method

#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
#code does the same thing, just different method
wild_BH <- p.adjust(wild_k6$pvalues,method="BH") #setting method
alpha10 <- 0.1
outliers_BH10 <- which(wild_BH < alpha10)
length(outliers_BH10) #same number as q values results -- 69329
#alpha 0.05
alpha5 <- 0.05
outliers_BH5.ints <- which(wild_BH < alpha5)
length(outliers_BH5.ints)
outliers_BH5 <- as.character(loci)[outliers_BH5.ints]
outliers_BH5 #same number as q value results -- 46463
write.csv(outliers_BH5, "adaptive_snps_BH5.csv") #writing csv file
#alpha .01
alpha01 <- 0.01
outliers_BH1.ints <- which(wild_BH < alpha01)
length(outliers_BH1.ints)
outliers_BH1 <- as.character(loci)[outliers_BH1.ints]
outliers_BH1 #24379 found

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
wild_BC <- p.adjust(wild_k6$pvalues,method="bonferroni")
alpha10 <- 0.1
outliers_BH10 <- which(wild_BC < alpha10)
length(outliers_BH10) #5411 resulting snps -- 40000 less than other methods

alpha5 <- 0.05
outliers_BC5.ints <- which(wild_BC < alpha5)
length(outliers_BC5.ints)
outliers_BC5 <- as.character(loci)[outliers_BC5.ints]
outliers_BC5 #5029 resulting snps
write.csv(outliers_BC5, "adaptive_snps_BC5.csv") #writing csv file

#If you want the pvalues for each SNP
snps_pvalues <- cbind(loci, wild_k6$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "All_Pvalues.txt", sep="\t", quote=FALSE)

####venn diagram to compare snps identified by all three methods####\
#reading in snps selected by each methods
qvalue_data <- read.csv("adaptive_snps_q05.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("adaptive_snps_BH5.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("adaptive_snps_BC5.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)
# Set up the directory for saving the Venn diagram
if (!dir.exists("venn_output")) {
  dir.create("venn_output")
}
#plotting venn diagram
venn.diagram(
  x = snp_lists,
  category.names = c("Q Value", "BH", "BC"),
  filename = "venn_output/snp_comparison_venn.png",
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3), alpha("#FDE725FF", 0.3)),
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  cat.fontface = "bold"
)

##################################################################################################################
#all lepa outlier loci tests
####PCA -- lepa####
l.pca.input <- read.pcadapt("lepa_thinned_ids.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
lepa_pca <- pcadapt(input = l.pca.input, K=30) #running pca with high k to check
plot(lepa_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(lepa_pca, option="screeplot", K=10, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(lepa_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(lepa_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(lepa_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(lepa_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 4 and 5: no structure?
PC56=plot(lepa_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 5 and 6: structure yes?
PC67=plot(lepa_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 6 and 7: structure, no not really
#thus, retain 5 pc's

#calculate variance on each axis
EV <- (lepa_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 18%, no sig decreases after 6 -- plateaus; so k=6

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
lepa_k5 <- pcadapt(l.pca.input, K = 5)
summary(lepa_k5)

####begin identifying outliers -- initial visualizations -- lepa####

#basic manhattan plot -- just for visualizing the candidate loci
plot(lepa_k5, option="manhattan")
#Q-Q Plots -- looking at how data is distributed
plot(lepa_k5, option="qqplot")
#Histograms of the test statistic and of the p-values
hist(lepa_k5$pvalues, xlab = "p-value",main = NULL, breaks = 50, col = "#2471A3")
#plotting p-values and looking at significance, are majority of samples non-sig (.5 >)
#samples with very low p-values are candidates for adaptive loci
plot(lepa_k5, option = "stat.distribution")

####outlier tests -- lepa####
## q-values  -- start looking for outlier loci
## for later analyses -- can exclude outlier to have a neutral and adaptive set of data for analyses
## do mult methods when doing outlier tests; several listed in samaba package
## can create a ven diagram to see where there is overlap; don't rely on just one method
##combined a better idea of what might be true outliers

#Choosing a cutoff for outlier detection
#To provide a list of outliers and choose a cutoff for outlier detection, there are several methods 
#that are listed below from the less conservative one to the more conservative one.
#q-value - least conservative
#bonferroni - most

#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs 
#for the geno3pops data, for an expected false discovery rate lower than 10%:
loci=lepa_genlight$loc.names #calling locus names in genlight file
qval <- qvalue(lepa_k5$pvalues)$qvalues #transformed p-values, making q-values column
alpha10 <- 0.1 #setting false discovery rate -- threshold; test at diff levels
outliers_qvalue10 <- which(qval < alpha10) #q values set at false discovery level of 10; looking for snps that have q value less than alpha
length(outliers_qvalue10) #How many snps detected?
# 457598 found
alpha5 <- 0.05 #5% more stat sig, but check at other levels
outliers_qvalue5.ints <- which(qval < alpha5)
length(outliers_qvalue5.ints)
# 410883 found;lost 40,000ish -- becoming more strict; can go low as want
outliers_qvalue5 <- as.character(loci)[outliers_qvalue5.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
outliers_qvalue5
write.csv(outliers_qvalue5, "lepa_adaptive_snps_q05.csv") #writing csv file
###snps found using q value method

#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
#code does the same thing, just different method
lepa_BH <- p.adjust(lepa_k5$pvalues,method="BH") #setting method
alpha10 <- 0.1
outliers_BH10 <- which(lepa_BH < alpha10)
length(outliers_BH10) #same number as q values results -- 457598
#alpha 0.05
alpha5 <- 0.05
outliers_BH5.ints <- which(lepa_BH < alpha5)
length(outliers_BH5.ints)
outliers_BH5 <- as.character(loci)[outliers_BH5.ints]
outliers_BH5 #same number as q value results -- 410883
write.csv(outliers_BH5, "lepa_adaptive_snps_BH5.csv") #writing csv file
#alpha .01
alpha01 <- 0.01
outliers_BH1.ints <- which(lepa_BH < alpha01)
length(outliers_BH1.ints)
outliers_BH1 <- as.character(loci)[outliers_BH1.ints]
outliers_BH1 #332181 found

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
lepa_BC <- p.adjust(lepa_k5$pvalues,method="bonferroni")
alpha10 <- 0.1
outliers_BH10 <- which(lepa_BC < alpha10)
length(outliers_BH10) #149169 resulting snps -- 200,000 less than other methods

alpha5 <- 0.05
outliers_BC5.ints <- which(lepa_BC < alpha5)
length(outliers_BC5.ints)
outliers_BC5 <- as.character(loci)[outliers_BC5.ints]
outliers_BC5 #143101 resulting snps
write.csv(outliers_BC5, "lepa_adaptive_snps_BC5.csv") #writing csv file

#If you want the pvalues for each SNP
snps_pvalues <- cbind(loci, lepa_k5$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "lepa_All_Pvalues.txt", sep="\t", quote=FALSE)

####venn diagram to compare snps identified by all three methods####
#reading in snps selected by each methods
qvalue_data <- read.csv("lepa_adaptive_snps_q05.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("lepa_adaptive_snps_BH5.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("lepa_adaptive_snps_BC5.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)
# Set up the directory for saving the Venn diagram
if (!dir.exists("venn_output")) {
  dir.create("venn_output")
}
#plotting venn diagram
venn.diagram(
  x = snp_lists,
  category.names = c("Q Value", "BH", "BC"),
  filename = "venn_output/lepa_snp_comparison_venn.png",
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3), alpha("#FDE725FF", 0.3)),
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  cat.fontface = "bold"
)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common)
write.table(all_common, "venn_output/lepa_adaptive_snps_allmethods.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

####extracting candidate loci and creating neutral and adaptive loci sets####
system("plink2 --vcf lepa_thinned_ids.vcf.gz --chr-set 18 --allow-extra-chr --extract lepa_adaptive_snps_allmethods.txt --export vcf bgz --out lepa_adaptive")
system("plink2 --vcf lepa_thinned_ids.vcf.gz --chr-set 18 --allow-extra-chr --exclude lepa_adaptive_snps_allmethods.txt --export vcf bgz --out lepa_neutral")


####Creating distance matrices for the trees####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/MEGA_adaptive_neutral")
install.packages("StAMPP")
install.packages("ape")
BiocManager::install("ggtree")
library(ggtree)
library(tidyverse)
library(StAMPP)
library(ape)
#load in vcfs
adaptive_snps <- read.vcfR("lepa_adaptive.vcf.gz")
neutral_snps <- read.vcfR("lepa_neutral.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations
l_pop <- read.csv("lepa_origins.csv")
l_df <- as.data.frame(l_pop) #creating dataframe from origins
l_df$ID <- paste0("24040DeY_", l_df$ID) #adding the prefix to the sample names

pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table
pop(adaptive_gl) <- pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))


#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Brazilian" = "#3B967f", "Generic" = "#D66857", 
                  "Refuge" = "#01004c", "Ranch" = "orchid")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")
###############################################################################################
#zoo adaptive
z.pca.input <- read.pcadapt("zoo_thinned_ids.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
zoo_pca <- pcadapt(input = z.pca.input, K=30) #running pca with high k to check
plot(zoo_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(zoo_pca, option="screeplot", K=15, plt.pkg = "ggplot") #zooming in on k15

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(zoo_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(zoo_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(zoo_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(zoo_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(zoo_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(zoo_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 6 and 7: structure, no not really
PC78=plot(zoo_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=zoo_genlight$pop) #pc 7 and 8: structure, no 
#thus, retain 6 pc's

#calculate variance on each axis
EV <- (zoo_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 18%, no sig decreases after 6 -- plateaus; so k=6

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
zoo_k6 <- pcadapt(z.pca.input, K = 6)
summary(zoo_k6)

####begin identifying outliers -- initial visualizations -- zoo####

#basic manhattan plot -- just for visualizing the candidate loci
plot(zoo_k6, option="manhattan")
#Q-Q Plots -- looking at how data is distributed
plot(zoo_k6, option="qqplot")
#Histograms of the test statistic and of the p-values
hist(zoo_k6$pvalues, xlab = "p-value",main = NULL, breaks = 50, col = "#2471A3")
#plotting p-values and looking at significance, are majority of samples non-sig (.5 >)
#samples with very low p-values are candidates for adaptive loci
plot(zoo_k6, option = "stat.distribution")

####outlier tests -- zoo####
#qvalue
loci=zoo_genlight$loc.names #calling locus names in genlight file
qval <- qvalue(zoo_k6$pvalues)$qvalues #transformed p-values, making q-values column
alpha5 <- 0.05 #5% more stat sig, but check at other levels
outliers_qvalue5.ints <- which(qval < alpha5)
length(outliers_qvalue5.ints)
outliers_qvalue5 <- as.character(loci)[outliers_qvalue5.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
outliers_qvalue5
write.csv(outliers_qvalue5, "zoo_adaptive_snps_q05.csv") #writing csv file

#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
#code does the same thing, just different method
zoo_BH <- p.adjust(zoo_k6$pvalues,method="BH") #setting method
#alpha 0.05
alpha5 <- 0.05
outliers_BH5.ints <- which(zoo_BH < alpha5)
length(outliers_BH5.ints)
outliers_BH5 <- as.character(loci)[outliers_BH5.ints]
outliers_BH5 #same number as q value results -- 410883
write.csv(outliers_BH5, "zoo_adaptive_snps_BH5.csv") #writing csv file


#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
zoo_BC <- p.adjust(zoo_k6$pvalues,method="bonferroni")

alpha5 <- 0.05
outliers_BC5.ints <- which(zoo_BC < alpha5)
length(outliers_BC5.ints)
outliers_BC5 <- as.character(loci)[outliers_BC5.ints]
outliers_BC5 #143101 resulting snps
write.csv(outliers_BC5, "zoo_adaptive_snps_BC5.csv") #writing csv file

#If you want the pvalues for each SNP
snps_pvalues <- cbind(loci, lepa_k5$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "lepa_All_Pvalues.txt", sep="\t", quote=FALSE)

####venn diagram to compare snps identified by all three methods####
#reading in snps selected by each methods
qvalue_data <- read.csv("zoo_adaptive_snps_q05.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("zoo_adaptive_snps_BH5.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("zoo_adaptive_snps_BC5.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)
#plotting venn diagram
venn.diagram(
  x = snp_lists,
  category.names = c("Q Value", "BH", "BC"),
  filename = "venn_output/lepa_snp_comparison_venn.png",
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3), alpha("#FDE725FF", 0.3)),
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  cat.fontface = "bold"
)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common)
write.table(all_common, "venn_output/zoo_adaptive_snps_allmethods.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

####extracting candidate loci and creating neutral and adaptive loci sets####
system("plink2 --vcf zoo_thinned_ids.vcf.gz --chr-set 18 --allow-extra-chr --extract zoo_adaptive_snps_allmethods.txt --export vcf bgz --out zoo_adaptive")
system("plink2 --vcf zoo_thinned_ids.vcf.gz --chr-set 18 --allow-extra-chr --exclude zoo_adaptive_snps_allmethods.txt --export vcf bgz --out zoo_neutral")

####Creating distance matrices for the trees -- zoo####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/MEGA_adaptive_neutral")
install.packages("StAMPP")
install.packages("ape")
BiocManager::install("ggtree")
library(ggtree)
library(tidyverse)
library(StAMPP)
library(ape)
#load in vcfs
adaptive_snps <- read.vcfR("lepa_adaptive.vcf.gz")
neutral_snps <- read.vcfR("lepa_neutral.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations
l_pop <- read.csv("lepa_origins.csv")
l_df <- as.data.frame(l_pop) #creating dataframe from origins
l_df$ID <- paste0("24040DeY_", l_df$ID) #adding the prefix to the sample names

pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table
pop(adaptive_gl) <- pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))


#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
adap_base +
  geom_tippoint(aes(color = population), size = 3) +
  geom_tiplab(aes(label = label), size = 3, hjust = -0.2) +
  theme_tree2() +
  theme(legend.position = "right")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  geom_tiplab(aes(label = label), size = 3, hjust = -0.2) +
  theme_tree2() +
  theme(legend.position = "right")

####create trees -- zoo####
#load in vcfs
adaptive_snps <- read.vcfR("zoo_adaptive.vcf.gz")
neutral_snps <- read.vcfR("zoo_neutral.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations
z_pop <- read.csv("zoo_origins.csv")
z_df <- as.data.frame(z_pop) #creating dataframe from origins
z_df$ID <- paste0("24040DeY_", z_df$ID) #adding the prefix to the sample names

pop_lookup <- setNames(as.character(z_df$Pop), z_df$ID)# Create a lookup table
pop(adaptive_gl) <- pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))


#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
z_pop_colors <- c("Brazilian" = "#3B967f", "Generic" = "#D66857")
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = z_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = z_pop_colors, name = "Population")


##########################manu edits#########################################################################
#contains code with the requested edits to the analysis during the thesis stage
#populations all filtered together, no HWE filter, and then separated into pairs for more in depth comparison
#data can be pre-filtered outside of samba for the following steps
####loading packages and working directory####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Samba/paired_outlier_test")
library(pcadapt)
library(vcfR)
library(qvalue)
library(VennDiagram)
library(grDevices)
library(dplyr)
library(adegenet)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(StAMPP)
library(ape)
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt")
getpackages()


####lepa first, redoing with the new filtering steps####

#creating a lookup data for population assignments -- can be used in future paired comparisons
l_pop <- read.csv("lepa_origins.csv")
l_df <- as.data.frame(l_pop)
l_df$ID <- paste0("24040DeY_", l_df$ID) #adding the prefix to the sample names
l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
lepa_vcf <- read.vcfR("lepa_OT_filt.vcf.gz", verbose = TRUE)
lepa_genlight <- vcfR2genlight(lepa_vcf)
lepa_genlight

pop(lepa_genlight) <- l_pop_lookup[indNames(lepa_genlight)]# Directly assign populations based on IDs
pop(lepa_genlight) <- factor(pop(lepa_genlight), levels = levels(pop(lepa_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
l.pca.input <- read.pcadapt("lepa_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
lepa_pca <- pcadapt(input = l.pca.input, K=30) #running pca with high k to check
plot(lepa_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(lepa_pca, option="screeplot", K=10, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(lepa_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(lepa_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(lepa_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(lepa_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(lepa_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(lepa_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(lepa_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(lepa_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=lepa_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (lepa_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 12%, no sig decreases after 7 -- plateaus; so k=7

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
lepa_k7 <- pcadapt(l.pca.input, K = 7)
summary(lepa_k7)

##initial visualizations
#basic manhattan plot -- just for visualizing the candidate loci
plot(lepa_k7, option="manhattan")
#Q-Q Plots -- looking at how data is distributed
plot(lepa_k7, option="qqplot")
#Histograms of the test statistic and of the p-values
hist(lepa_k7$pvalues, xlab = "p-value",main = NULL, breaks = 50, col = "#2471A3")
#plotting p-values and looking at significance, are majority of samples non-sig (.5 >)
#samples with very low p-values are candidates for adaptive loci
plot(lepa_k7, option = "stat.distribution")

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=lepa_genlight$loc.names #calling locus names in genlight file
qval <- qvalue(lepa_k7$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #316613 found

outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
outliers_qvalue01
write.csv(outliers_qvalue01, "lepa_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
lepa_BH <- p.adjust(lepa_k7$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(lepa_BH < alpha01)
length(outliers_BH01.ints)#316613 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "lepa_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
lepa_BC <- p.adjust(lepa_k7$pvalues,method="bonferroni")

alpha01 <- 0.01
outliers_BC01.ints <- which(lepa_BC < alpha01)
length(outliers_BC01.ints)#69243 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "lepa_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("lepa_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("lepa_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("lepa_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)
# Set up the directory for saving the Venn diagram
if (!dir.exists("venn_output")) {
  dir.create("venn_output")
}
#plotting venn diagram
venn.diagram(
  x = snp_lists,
  category.names = c("Q Value", "BH", "BC"),
  filename = "venn_output/lepa_snp_comparison_venn.png",
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3), alpha("#FDE725FF", 0.3)),
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.col = c("#440154FF", "#21908CFF", "#FDE725FF"),
  cat.fontface = "bold"
)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #69243
write.table(all_common, "lepa_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf lepa_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract lepa_adaptive_allmethods_manufilt.txt --export vcf bgz --out lepa_adaptive_new")
system("plink2 --vcf lepa_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude lepa_adaptive_allmethods_manufilt.txt --export vcf bgz --out lepa_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("lepa_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("lepa_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Brazilian" = "#3B967f", "Generic" = "#D66857", 
                  "Refuge" = "#01004c", "Ranch" = "orchid")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

####generic and brazilian####
#overwriting many r objects generated from the lepa section, make sure you have the outputs saved
#using the lookup table made in lepa section, found below, uncomment to use
#l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
GB_vcf <- read.vcfR("gen_braz_OT_filt.vcf.gz", verbose = TRUE)
GB_genlight <- vcfR2genlight(GB_vcf)

pop(GB_genlight) <- l_pop_lookup[indNames(GB_genlight)]# Directly assign populations based on IDs
pop(GB_genlight) <- factor(pop(GB_genlight), levels = levels(pop(GB_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
GB.pca.input <- read.pcadapt("gen_braz_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
GB_pca <- pcadapt(input = GB.pca.input, K=30) #running pca with high k to check
plot(GB_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(GB_pca, option="screeplot", K=20, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(GB_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(GB_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(GB_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(GB_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(GB_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(GB_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(GB_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(GB_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=GB_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (GB_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 10%, no sig decreases after 8 -- plateaus; so k=8

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
GB_k8 <- pcadapt(GB.pca.input, K = 8)
summary(GB_k8)

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=GB_genlight$loc.names #calling locus names in genlight file
GB_qval <- qvalue(GB_k8$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(GB_qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #99454 found
outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
write.csv(outliers_qvalue01, "GB_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
GB_BH <- p.adjust(GB_k8$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(GB_BH < alpha01)
length(outliers_BH01.ints)#99454 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "GB_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
GB_BC <- p.adjust(GB_k8$pvalues,method="bonferroni")
alpha01 <- 0.01
outliers_BC01.ints <- which(GB_BC < alpha01)
length(outliers_BC01.ints)#12731 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "GB_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("GB_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("GB_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("GB_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #7053
write.table(all_common, "GB_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf gen_braz_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract GB_adaptive_allmethods_manufilt.txt --export vcf bgz --out GB_adaptive_new")
system("plink2 --vcf gen_braz_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude GB_adaptive_allmethods_manufilt.txt --export vcf bgz --out GB_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("GB_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("GB_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Brazilian" = "#3B967f", "Generic" = "#D66857")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

####Ranch and refuge####
#overwriting many r objects generated from the lepa section, make sure you have the outputs saved
#using the lookup table made in lepa section, found below, uncomment to use
#l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
RR_vcf <- read.vcfR("ref_ran_OT_filt.vcf.gz", verbose = TRUE)
RR_genlight <- vcfR2genlight(RR_vcf)

pop(RR_genlight) <- l_pop_lookup[indNames(RR_genlight)]# Directly assign populations based on IDs
pop(RR_genlight) <- factor(pop(RR_genlight), levels = levels(pop(RR_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
RR.pca.input <- read.pcadapt("ref_ran_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
RR_pca <- pcadapt(input = RR.pca.input, K=30) #running pca with high k to check
plot(RR_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(RR_pca, option="screeplot", K=10, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(RR_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(RR_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(RR_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(RR_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(RR_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(RR_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(RR_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(RR_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=RR_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (RR_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 10%, no sig decreases after 8 -- plateaus; so k=8

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
RR_k6 <- pcadapt(RR.pca.input, K = 6)
summary(RR_k6)

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=RR_genlight$loc.names #calling locus names in genlight file
RR_qval <- qvalue(RR_k6$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(RR_qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #106771 found
outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
write.csv(outliers_qvalue01, "RR_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
RR_BH <- p.adjust(RR_k6$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(RR_BH < alpha01)
length(outliers_BH01.ints)#106771 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "RR_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
RR_BC <- p.adjust(RR_k6$pvalues,method="bonferroni")
alpha01 <- 0.01
outliers_BC01.ints <- which(RR_BC < alpha01)
length(outliers_BC01.ints)#13868 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "RR_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("RR_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("RR_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("RR_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #5506
write.table(all_common, "RR_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf ref_ran_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract RR_adaptive_allmethods_manufilt.txt --export vcf bgz --out RR_adaptive_new")
system("plink2 --vcf ref_ran_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude RR_adaptive_allmethods_manufilt.txt --export vcf bgz --out RR_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("RR_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("RR_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Refuge" = "#01004c", "Ranch" = "orchid")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

####generic vs refuge####
#overwriting many r objects generated from the lepa section, make sure you have the outputs saved
#using the lookup table made in lepa section, found below, uncomment to use
#l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
GF_vcf <- read.vcfR("gen_ref_OT_filt.vcf.gz", verbose = TRUE)
GF_genlight <- vcfR2genlight(GF_vcf)

pop(GF_genlight) <- l_pop_lookup[indNames(GF_genlight)]# Directly assign populations based on IDs
pop(GF_genlight) <- factor(pop(GF_genlight), levels = levels(pop(GF_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
GF.pca.input <- read.pcadapt("gen_ref_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
GF_pca <- pcadapt(input = GF.pca.input, K=30) #running pca with high k to check
plot(GF_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(GF_pca, option="screeplot", K=10, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(GF_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(GF_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(GF_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(GF_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(GF_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(GF_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(GF_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(GF_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=GF_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (GF_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 10%, no sig decreases after 8 -- plateaus; so k=8

#plot nice pca
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (18.5%)")+
  ylab("PC2 (7.3%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))

#set your K values for future tests
GF_k4 <- pcadapt(GF.pca.input, K = 4)
summary(GF_k4)

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=GF_genlight$loc.names #calling locus names in genlight file
GF_qval <- qvalue(GF_k4$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(GF_qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #216578 found
outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
write.csv(outliers_qvalue01, "GF_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
GF_BH <- p.adjust(GF_k4$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(GF_BH < alpha01)
length(outliers_BH01.ints)#216578 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "GF_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
GF_BC <- p.adjust(GF_k4$pvalues,method="bonferroni")
alpha01 <- 0.01
outliers_BC01.ints <- which(GF_BC < alpha01)
length(outliers_BC01.ints)#36255 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "GF_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("GF_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("GF_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("GF_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #4188
write.table(all_common, "GF_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf gen_ref_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract GF_adaptive_allmethods_manufilt.txt --export vcf bgz --out GF_adaptive_new")
system("plink2 --vcf gen_ref_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude GF_adaptive_allmethods_manufilt.txt --export vcf bgz --out GF_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("GF_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("GF_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Generic" = "#D66857", 
                  "Refuge" = "#01004c")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

####generic vs ranch####
#overwriting many r objects generated from the lepa section, make sure you have the outputs saved
#using the lookup table made in lepa section, found below, uncomment to use
#l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
GA_vcf <- read.vcfR("gen_ran_OT_filt.vcf.gz", verbose = TRUE)
GA_genlight <- vcfR2genlight(GA_vcf)

pop(GA_genlight) <- l_pop_lookup[indNames(GA_genlight)]# Directly assign populations based on IDs
pop(GA_genlight) <- factor(pop(GA_genlight), levels = levels(pop(GA_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
GA.pca.input <- read.pcadapt("gen_ran_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
GA_pca <- pcadapt(input = GA.pca.input, K=30) #running pca with high k to check
plot(GA_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(GA_pca, option="screeplot", K=10, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(GA_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(GA_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(GA_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(GA_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(GA_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(GA_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(GA_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(GA_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=GA_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (GA_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 10%, no sig decreases after 8 -- plateaus; so k=8

#set your K values for future tests
GA_k5 <- pcadapt(GA.pca.input, K = 5)
summary(GA_k5)

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=GA_genlight$loc.names #calling locus names in genlight file
GA_qval <- qvalue(GA_k5$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(GA_qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #394033 found
outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
write.csv(outliers_qvalue01, "GA_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
GA_BH <- p.adjust(GA_k5$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(GA_BH < alpha01)
length(outliers_BH01.ints)#394033 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "GA_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
GA_BC <- p.adjust(GA_k5$pvalues,method="bonferroni")
alpha01 <- 0.01
outliers_BC01.ints <- which(GA_BC < alpha01)
length(outliers_BC01.ints)#88608 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "GA_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("GA_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("GA_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("GA_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #8348
write.table(all_common, "GA_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf gen_ran_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract GA_adaptive_allmethods_manufilt.txt --export vcf bgz --out GA_adaptive_new")
system("plink2 --vcf gen_ran_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude GA_adaptive_allmethods_manufilt.txt --export vcf bgz --out GA_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("GA_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("GA_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Generic" = "#D66857", 
                  "Ranch" = "orchid")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

####brazilian vs refuge####
#overwriting many r objects generated from the lepa section, make sure you have the outputs saved
#using the lookup table made in lepa section, found below, uncomment to use
#l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
BF_vcf <- read.vcfR("braz_ref_OT_filt.vcf.gz", verbose = TRUE)
BF_genlight <- vcfR2genlight(BF_vcf)

pop(BF_genlight) <- l_pop_lookup[indNames(BF_genlight)]# Directly assign populations based on IDs
pop(BF_genlight) <- factor(pop(BF_genlight), levels = levels(pop(BF_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
BF.pca.input <- read.pcadapt("braz_ref_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
BF_pca <- pcadapt(input = BF.pca.input, K=30) #running pca with high k to check
plot(BF_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(BF_pca, option="screeplot", K=15, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(BF_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(BF_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(BF_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(BF_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(BF_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(BF_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(BF_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(BF_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=BF_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (BF_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 10%, no sig decreases after 8 -- plateaus; so k=8

#set your K values for future tests
BF_k7 <- pcadapt(BF.pca.input, K = 7)
summary(BF_k7)

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=BF_genlight$loc.names #calling locus names in genlight file
BF_qval <- qvalue(BF_k7$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(BF_qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #338551 found
outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
write.csv(outliers_qvalue01, "BF_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
BF_BH <- p.adjust(BF_k7$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(BF_BH < alpha01)
length(outliers_BH01.ints)#338551 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "BF_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
BF_BC <- p.adjust(BF_k7$pvalues,method="bonferroni")
alpha01 <- 0.01
outliers_BC01.ints <- which(BF_BC < alpha01)
length(outliers_BC01.ints)#54773 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "BF_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("BF_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("BF_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("BF_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #18786
write.table(all_common, "BF_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf braz_ref_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract BF_adaptive_allmethods_manufilt.txt --export vcf bgz --out BF_adaptive_new")
system("plink2 --vcf braz_ref_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude BF_adaptive_allmethods_manufilt.txt --export vcf bgz --out BF_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("BF_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("BF_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Brazilian" = "#3B967f", 
                  "Refuge" = "#01004c")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

####brazilian vs ranch####
#overwriting many r objects generated from the lepa section, make sure you have the outputs saved
#using the lookup table made in lepa section, found below, uncomment to use
#l_pop_lookup <- setNames(as.character(l_df$Pop), l_df$ID)# Create a lookup table

#reading in as genlight and adding pop id's
BA_vcf <- read.vcfR("braz_ran_OT_filt.vcf.gz", verbose = TRUE)
BA_genlight <- vcfR2genlight(BA_vcf)

pop(BA_genlight) <- l_pop_lookup[indNames(BA_genlight)]# Directly assign populations based on IDs
pop(BA_genlight) <- factor(pop(BA_genlight), levels = levels(pop(BA_genlight)))# Convert to factor with original levels

###PCA with the lepa data -- start of outlier process
BA.pca.input <- read.pcadapt("braz_ran_OT_filt.vcf.gz", type = "vcf") #making the vcf readable by pcaadapt
BA_pca <- pcadapt(input = BA.pca.input, K=30) #running pca with high k to check
plot(BA_pca, option = "screeplot", plt.pkg = "ggplot") #plotting pca screenplot
plot(BA_pca, option="screeplot", K=15, plt.pkg = "ggplot") #zooming in on k10

#making the pca plot and comparing pc's, looking to see where pc's no longer become informative -- or display no more structure
PC12=plot(BA_pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 1 and 2: structure yes
PC23=plot(BA_pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 2 and 3: strucutre yes
PC34=plot(BA_pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 3 and 4: structure yes
PC45=plot(BA_pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 4 and 5: structure yes
PC56=plot(BA_pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 5 and 6: structure yes
PC67=plot(BA_pca, option = "scores", i = 6, j = 7, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 6 and 7: structure yes
PC78=plot(BA_pca, option = "scores", i = 7, j = 8, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 7 and 8: structure yes
PC89=plot(BA_pca, option = "scores", i = 8, j = 9, plt.pkg = "ggplot", pop=BA_genlight$pop) #pc 8 and 9: structure not really


#calculate variance on each axis
EV <- (BA_pca$singular.values^2) #column created -- singular values squared to get variance
EV2 <-EV*100 #multiplied by 100 to get variance
EV2
#1 explains 10%, no sig decreases after 8 -- plateaus; so k=8

#set your K values for future tests
BA_k8 <- pcadapt(BA.pca.input, K = 8)
summary(BA_k8)

##outlier tests
#q-value - least conservative
#bonferroni - most
#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs.

#q value
loci=BA_genlight$loc.names #calling locus names in genlight file
BA_qval <- qvalue(BA_k8$pvalues)$qvalues #transformed p-values, making q-values column
alpha01 <- 0.01 #1% more stat sig, setting alpha level
outliers_qvalue01.ints <- which(BA_qval < alpha01) #identifing loci
length(outliers_qvalue01.ints) #551142 found
outliers_qvalue01 <- as.character(loci)[outliers_qvalue01.ints] #binding loci name to the info on the snps -- SNP name tied to ones it picked up with q values
write.csv(outliers_qvalue01, "BA_adaptive_snps_q01.csv") #writing csv file


#Benjamini-Hochberg Procedure -- another way to detect -- a bit more strict
BA_BH <- p.adjust(BA_k8$pvalues,method="BH") #setting method
#alpha .01
alpha01 <- 0.01
outliers_BH01.ints <- which(BA_BH < alpha01)
length(outliers_BH01.ints)#551142 found
outliers_BH01 <- as.character(loci)[outliers_BH1.ints]
write.csv(outliers_BH01, "BA_adaptive_snps_BH01.csv") #writing csv file

#Bonferroni correction -- accounting for smaller sample sizes; most conservative approach
BA_BC <- p.adjust(BA_k8$pvalues,method="bonferroni")
alpha01 <- 0.01
outliers_BC01.ints <- which(BA_BC < alpha01)
length(outliers_BC01.ints)#112935 resulting snps
outliers_BC01 <- as.character(loci)[outliers_BC01.ints]
write.csv(outliers_BC01, "BA_adaptive_snps_BC01.csv") #writing csv file

##venn diagram to compare snps identified by all three methods
#reading in snps selected by each methods
qvalue_data <- read.csv("BA_adaptive_snps_q01.csv", stringsAsFactors = FALSE)
names(qvalue_data)[2] <- "snp_ID" #changing column name
BH_data <- read.csv("BA_adaptive_snps_BH01.csv", stringsAsFactors = FALSE)
names(BH_data)[2] <- "snp_ID"
BC_data <- read.csv("BA_adaptive_snps_BC01.csv", stringsAsFactors = FALSE)
names(BC_data)[2] <- "snp_ID"
#extacting snps
qvalue_snps <- qvalue_data$snp_ID
BH_snps <- BH_data$snp_ID
BC_snps <- BC_data$snp_ID
#creating a list for the venn diagram
snp_lists <- list("Q value" = qvalue_snps,
                  "BH" = BH_snps,
                  "BC" = BC_snps)

#identifying and saving snps selected by all categories
all_common <- intersect(intersect(qvalue_snps, BH_snps), BC_snps)
n_all_common <- length(all_common) #29452
write.table(all_common, "BA_adaptive_allmethods_manufilt.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

##extracting values and creating new vcf
system("plink2 --vcf braz_ran_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --extract BA_adaptive_allmethods_manufilt.txt --export vcf bgz --out BA_adaptive_new")
system("plink2 --vcf braz_ran_OT_filt.vcf.gz --chr-set 18 --allow-extra-chr --exclude BA_adaptive_allmethods_manufilt.txt --export vcf bgz --out BA_neutral_new")

##create distance matrix for trees
#load in vcfs
adaptive_snps <- read.vcfR("BA_adaptive_new.vcf.gz")
neutral_snps <- read.vcfR("BA_neutral_new.vcf.gz")
#convert to genlight
adaptive_gl <- vcfR2genlight(adaptive_snps)
neutral_gl <- vcfR2genlight(neutral_snps)
#adding populations -- using lookup made in line 796
pop(adaptive_gl) <- l_pop_lookup[indNames(adaptive_gl)]# Directly assign populations based on IDs
pop(adaptive_gl) <- factor(pop(adaptive_gl), levels = levels(pop(adaptive_gl)))# Convert to factor with original levels
pop(neutral_gl) <- l_pop_lookup[indNames(neutral_gl)]# Directly assign populations based on IDs
pop(neutral_gl) <- factor(pop(neutral_gl), levels = levels(pop(neutral_gl)))

#calculate Nei's genetic distance
adaptive_dist <- stamppNeisD(adaptive_gl, pop = FALSE)
neutral_dist <- stamppNeisD(neutral_gl, pop = FALSE)
#use ape to make trees
adapt_dist_ob <- as.dist(adaptive_dist)
adaptive_nj_tree <- nj(adapt_dist_ob)
plot(adaptive_nj_tree)
neut_dist_ob <- as.dist(neutral_dist)
neut_nj_tree <- nj(neut_dist_ob)
plot(neut_nj_tree)

##make nice plots of the trees
#extract pop data
pop_data <- data.frame(sample = indNames(adaptive_gl),
                       population = as.character(pop(adaptive_gl)))

#cleaning names
pop_data_clean <- pop_data %>%
  mutate(sample = gsub("-.*", "", sample),
         sample = gsub("24040DeY_", "", sample))
adaptive_nj_tree$tip.label <- pop_data_clean$sample
neut_nj_tree$tip.label <- pop_data_clean$sample

#create tree
adap_base <- ggtree(adaptive_nj_tree, layout = "unrooted") %<+% pop_data_clean
l_pop_colors <- c("Brazilian" = "#3B967f", 
                  "Ranch" = "orchid")
adap_base +
  geom_tippoint(aes(color = population), size = 2) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")

neut_base <- ggtree(neut_nj_tree, layout = "unrooted") %<+% pop_data_clean
neut_base +
  geom_tippoint(aes(color = population), size = 3) +
  theme_tree2() +
  theme(legend.position = "right") +
  scale_color_manual(values = l_pop_colors, name = "Population")