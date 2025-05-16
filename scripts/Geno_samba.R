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
l_vcf <- read.vcfR("lepa_samba_nofilt_thin.vcf.gz", verbose = TRUE)
lepa_genlight <- vcfR2genlight(l_vcf)
lepa_genlight
##adding population assignments
l_pop <- read.csv("lepa_origins.csv")
l_df <- as.data.frame(l_pop) #creating dataframe from origins

strata(lepa_genlight) <- l_df #creating a strata for pop id's
setPop(lepa_genlight) <- ~Pop #setting the pop from the strata

####initial look at data using samba####
#samba initialization and generation of quality control plots
genlight2sambar("lepa_genlight",do_confirm=TRUE)

#filtering
filterdata(snpmiss=0.1,min_mac=8,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=500, nchroms=18, silent=TRUE,maxprop_hefilter = 0.06)

#calculate kinship
calckinship()

#genetic diversity measurements
calcdiversity(do_venn = FALSE) 
####################################################################################################################################
#data can be pre-filtered outside of samba for the following steps
####loading packages and working directory####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Samba_analyses")
library(pcadapt)
library(vcfR)
library(qvalue)
####loading and prepping data####
z_vcf <- read.vcfR("zoo_dapc_thin.vcf.gz", verbose = TRUE)
zoo_genlight <- vcfR2genlight(z_vcf)
zoo_genlight

w_vcf <- read.vcfR("wild_dapc_thin.vcf.gz", verbose = TRUE)
wild_genlight <- vcfR2genlight(w_vcf)
wild_genlight

lepa_vcf <- read.vcfR("lepa_dapc_thin.vcf.gz", verbose = TRUE)
lepa_genlight <- vcfR2genlight(lepa_vcf)
lepa_genlight

##adding population assignments
w_pop <-read.csv("wild_origins.csv") #reading in origins
z_pop <- read.csv("zoo_origins.csv")
l_pop <- read.csv("lepa_origins.csv")
w_df <- as.data.frame(w_pop) #creating dataframe from origins
z_df <-as.data.frame(z_pop)
l_df <- as.data.frame(l_pop)

strata(wild_genlight) <- w_df #creating a strata for pop id's
setPop(wild_genlight) <- ~Pop #setting the pop from the strata
strata(zoo_genlight) <- z_df
setPop(zoo_genlight) <- ~Pop
strata(lepa_genlight) <- l_df
setPop(lepa_genlight) <- ~Pop
####PCA####
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

####begin identifing outliers####
plot(wild_k6, option="manhattan")




















