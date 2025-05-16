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
install.packages("pcadapt")
install.packages("Rtools")
BiocManager::install("LEA")
library(devtools)
install_github("drveera/ggman")
install_github("jdstorey/qvalue")
devtools::install_github("pievos101/PopGenome")
library(LEA)
library(vcfR)
library(adegenet)
library(tidyverse)
library(hierfstat)
library(ggman)
library(ggplot2)
library(pcadapt)
library(qvalue)
library(PopGenome)

#installing and loading samba
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt")
getpackages()
###############################################################################################################
####loading in data as genlights####
Z_vcf <- read.vcfR("zoo_dapc_thin.vcf.gz", verbose = TRUE)
zoo_genlight <- vcfR2genlight(Z_vcf)
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

####initial look at data using samba####
#samba initialization and generation of quality control plots
genlight2sambar("wild_genlight",do_confirm=TRUE)
genlight2sambar("zoo_genlight",do_confirm=TRUE)
genlight2sambar("lepa_genlight",do_confirm=TRUE)

#genetic diversity measurements
inds$filter <- TRUE #next line is checking to make sure all data passes filtering, assigning passing status
snps$filter <- TRUE #mark all snps as passing the filter as well
calcdiversity(dohwe = FALSE, nsites = 1953946) 







