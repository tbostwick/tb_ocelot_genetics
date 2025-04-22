##LEPA Initial Analyses -- plink
#4/8/25
#Cleaned out unnecessary steps such as other LD pruning and refiltering settings
    #only kept code pertaining to the analyses going forward, removed test code
#genomic analyses and refiltering -- combine both the analysis and the refiltering code
#analyses run with the following settings:
    #MAF: 0.05, MISS: 90%, Biallelic, LD: 50 5 2 (corresponds to r2 of 0.5), HWE
    #contains both the code for PLINK analysis and adegenet DAPC
##########################################################
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/SNP_refiltered")
library(ggplot2)
library(dplyr)
install.packages("HardyWeinberg")
library(HardyWeinberg)

##########################################################
#Data Prep, filtering and LD Pruning
####Make bed bim fam files from vcf####
system("plink --vcf 24041DeY-snp_filter_dp_gq__AllChrom_.vcf.gz --keep-allele-order --allow-extra-chr --chr-set 18 --make-bed --out SNP_dp_gq_AllChrom_AllInd")

####Subsetting data by origin -- separate zoo and wild for filtering####
#invokes plink, calls the original file, "--keep" and then feed it a txt file with only the individuals you want
#tell it the number of chromosomes, make new bed bim fam files, and then what the output should be named
#create zoo file
system("plink --bfile SNP_dp_gq_AllChrom_AllInd --keep pop_subset_zoo.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Zoo_BaseFilter_Allchrom")
#create wild ocelots, no duplicates file
system("plink --bfile SNP_dp_gq_AllChrom_AllInd --keep pop_subset_wild_remove_dups.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Wild_BaseFilter_AllChrom")

####Filtering -- MAF 0.05 MISS 90 Bialleic####
##Zoo filtering
#filter MAF and missingness and hwe
system("plink --bfile Zoo_BaseFilter_Allchrom --chr-set 18 --allow-extra-chr --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out Zoo_maf05_miss90_hwe")
#filter biallelic
system("plink --bfile Zoo_maf05_miss90_hwe --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Zoo_refiltered")

##Wild filtering
#filter MAF and missingness
system("plink --bfile Wild_BaseFilter_AllChrom --chr-set 18 --allow-extra-chr --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out Wild_maf05_miss90_hwe")
#filter biallelic
system("plink --bfile Wild_maf05_miss90 --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Wild_refiltered")

####Filtering -- ROH filters no MAF, miss 90, biallelic####
#filter MAF and missingness and hwe
system("plink --bfile Zoo_BaseFilter_Allchrom --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Zoo_miss90_hwe")
#filter biallelic
system("plink --bfile Zoo_miss90_hwe --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Zoo_roh_filter")

##Wild filtering
#filter MAF and missingness
system("plink --bfile Wild_BaseFilter_AllChrom --chr-set 18 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Wild_miss90_hwe")
#filter biallelic
system("plink --bfile Wild_miss90_hwe --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Wild_roh_filter")


####Subsetting filtered data by population####
#create ranch file with dispersers included
system("plink --bfile Wild_refiltered --keep pop_subset_wild_ranch_wdispers.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Ranch_Refiltered")
#create refuge file with dispersers included
system("plink --bfile Wild_refiltered --keep pop_subset_wild_refuge_wdispers.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Refuge_Refiltered")
#create a brazilian file
system("plink --bfile Zoo_refiltered --keep pop_subset_brazilian.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Brazilain_refiltered")
#create a generic file
system("plink --bfile Zoo_refiltered --keep pop_subset_generic.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Generic_refiltered")

####Linkage Disequilibrium####
##Zoo LD Pruning
  #LD pruning spits out an error about duplicate ID names, in order for LD pruning to run smoothly the next line of code tells plink
  #to name variants under that specific naming scheme -- this doesn't affect nearly any other analysis, except for LD pruning
  #how the naming works is beyond me, but it sets the naming scheme so that all variants have a unique name which allows the 
  #pruning code to make it's comparisons in the second line and remove the proper variants in the last line
system("plink2 --bfile Zoo_refiltered --chr-set 18 --keep-allele-order --allow-extra-chr --set-all-var-ids @:#$r,$a --new-id-max-allele-len 323 --make-bed --out Zoo_refiltered_uniqueID") #set variant ID
system("plink --bfile Zoo_refiltered_uniqueID --chr-set 18 --keep-allele-order --allow-extra-chr --memory 24000 --indep 50 5 2 --out Zoo_RF_LDpruned_0.5") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Zoo_refiltered_uniqueID --extract Zoo_RF_LDpruned_0.5.prune.in --chr-set 18 --allow-extra-chr --memory 24000 --make-bed --out Zoo_LDpruned_0.5_RF") #extract SNPs and create new files
      #reduces SNPS from 27003602 to 4,475,851
#test effects with PCA
system("plink --bfile Zoo_LDpruned_0.5_RF --allow-extra-chr --pca --chr-set 18 --out pca_zoo_LD_0.5") #run the pca code, specified the number of chromosomes
pca.zoo_ld0.5 <- read.table("pca_zoo_LD_0.5.eigenvec", header=FALSE) #load pca results from the eigenvec file
zoo_origins <- read.csv("zoo_origins.csv", header = TRUE) #read in origin metadata
pca.zoo_ld0.5.origins <- left_join(pca.zoo_ld0.5, zoo_origins, by = "V2") #append origin data to pca data
ggplot(pca.zoo_ld0.5.origins, aes(x=V3,y=V4, color = Cat.Group)) +  #plot with individual ID's and by origin
  geom_point() +
  geom_text(aes(label=V2), vjust=1, hjust=1, size=2) + #V2 is ID
  scale_color_manual(name = "Origin", values = c("Generic" = "dodgerblue3", "Brazilian" = "coral2")) +
  labs(x = "PC1", y = "PC2", title = "Zoo PCA by Individual") +
  theme_minimal()
#Seems to invert PC1 and PC2

##Wild LD Pruning
system("plink2 --bfile Wild_refiltered --chr-set 18 --keep-allele-order --allow-extra-chr --set-all-var-ids @:#$r,$a --new-id-max-allele-len 323 --make-bed --out Wild_refiltered_uniqueID") #set variant ID
system("plink --bfile Wild_refiltered_uniqueID --chr-set 18 --keep-allele-order --allow-extra-chr --memory 24000 --indep 50 5 2 --out Wild_RF_LDpruned_0.5") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Wild_refiltered_uniqueID --extract Wild_RF_LDpruned_0.5.prune.in --chr-set 18 --allow-extra-chr --memory 24000 --make-bed --out Wild_LDpruned_0.5_RF") #extract SNPs and create new files
    #reduces SNPS from 12146279 to 1422635
#test effects with PCA
system("plink --bfile Wild_LDpruned_0.5_RF --allow-extra-chr --pca --chr-set 18 --out pca_wild_LD_0.5") #run the pca code, specified the number of chromosomes
pca.wild_ld0.5 <- read.table("pca_wild_LD_0.5.eigenvec", header=FALSE) #load pca results from the eigenvec file
wild_origins <- read.csv("wild_origins.csv", header = TRUE) #read in origins file
pca.wild.origins_ld0.5 <- left_join(pca.wild_ld0.5, wild_origins, by = "V2") #append origin data to pca data
ggplot(pca.wild.origins_ld0.5, aes(x=V3,y=V4, color = Cat.Group)) +  #plot with individual ID's and by origin
  geom_point() +
  geom_text(aes(label=V2), vjust=1, hjust=1, size=2) + #V2 is ID
  scale_color_manual(name = "Origin", values = c("Refuge" = "dodgerblue3", "Ranch" = "coral2")) +
  labs(x = "PC1", y = "PC2", title = "Wild PCA by Individual") +
  theme_minimal()
#seems to invert PC1 and PC2

####subset data after pruning into populations and origins
#create ranch file with dispersers included
system("plink --bfile Wild_LDpruned_0.5_RF --keep pop_subset_wild_ranch_wdispers.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Ranch_LDpruned_0.5_RF")
#create refuge file with dispersers included
system("plink --bfile Wild_LDpruned_0.5_RF --keep pop_subset_wild_refuge_wdispers.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Refuge_LDpruned_0.5_RF")
#create a brazilian file
system("plink --bfile Zoo_LDpruned_0.5_RF --keep pop_subset_brazilian.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Brazilian_LDpruned_0.5_RF")
#create a generic file
system("plink --bfile Zoo_LDpruned_0.5_RF --keep pop_subset_generic.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Generic_LDpruned_0.5_RF")

##test -- pruning populations separate from each other, then running inbreeding
#ranch
system("plink2 --bfile Ranch_Refiltered --chr-set 18 --keep-allele-order --allow-extra-chr --set-all-var-ids @:#$r,$a --new-id-max-allele-len 323 --make-bed --out Ranch_refiltered_uniqueID") #set variant ID
system("plink --bfile Ranch_refiltered_uniqueID --chr-set 18 --keep-allele-order --allow-extra-chr --memory 24000 --indep 50 5 2 --out Ranch_isolatetest_LDpruned_0.5") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Ranch_refiltered_uniqueID --extract Ranch_isolatetest_LDpruned_0.5.prune.in --chr-set 18 --allow-extra-chr --memory 24000 --make-bed --out Ranch_isolatetest_LDpruned_0.5") #extract SNPs and create new files
#reduces SNPS from 12146279 to 1255109
#refuge
system("plink2 --bfile Refuge_Refiltered --chr-set 18 --keep-allele-order --allow-extra-chr --set-all-var-ids @:#$r,$a --new-id-max-allele-len 323 --make-bed --out Refuge_refiltered_uniqueID") #set variant ID
system("plink --bfile Refuge_refiltered_uniqueID --chr-set 18 --keep-allele-order --allow-extra-chr --memory 24000 --indep 50 5 2 --out Refuge_isolatetest_LDpruned_0.5") #makes an out and in files of SNps to keep and SNPs to remove
system("plink --bfile Refuge_refiltered_uniqueID --extract Refuge_isolatetest_LDpruned_0.5.prune.in --chr-set 18 --allow-extra-chr --memory 24000 --make-bed --out Refuge_isolatetest_LDpruned_0.5") #extract SNPs and create new files
#reduces SNPS from 12146279 to 1135734
####thinning####
#--bp-space thins data by spacing out snps by 1000bps apart
system("plink --bfile Wild_roh_filter --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out wild_roh_thin")

#########################################################
#Data Analysis
####Data exploration and quality checks####
###hardy weinberg equilibrium testing
system("plink --bfile Zoo_Biallelic_Allchrom --chr-set 18 --hardy --out zoo_hwe") #zoo hwe
system("plink --bfile Wild_Biallelic_AllChrom --chr-set 18 --hardy --out wild_hwe") #wild hwe

##make ternary plot
#data wrangling to fit format for plotting
hardy_test <- read.table("zoo_hwe.hwe", header = TRUE) #read in data from the hwe test
genotype_counts <- do.call(rbind, strsplit(as.character(hardy_test$GENO), "/")) #split the geno column into three -- in order to plot the values
genotype_counts <- as.data.frame(genotype_counts) #make into a data from
colnames(genotype_counts) <- c("AA", "AB", "BB") #label columns
genotype_counts$AA <- as.numeric(genotype_counts$AA) #make into numeric format
genotype_counts$AB <- as.numeric(genotype_counts$AB) #make into numeric format
genotype_counts$BB <- as.numeric(genotype_counts$BB) #make into numeric format
head(genotype_counts)
#plot
HWTernaryPlot(genotype_counts, markercol = rgb(0,0,1,0.03),
              vertexlab = c("AA", "AB", "BB")) #so many snps it takes forever
#HWE filtering -- removing snps with extreme p-values
system("plink --bfile Zoo_Biallelic_Allchrom --chr-set 18 --hwe 1e-30 --make-bed --out Zoo_Biallelic_HWE_AllChrom")
    #did not remove any variants


####PCA####
#beyond just making a pca, this is good to check if the data was input correctly

#PCA on zoo individuals
system("plink --bfile Zoo_refiltered --allow-extra-chr --pca --chr-set 18 --out pca_zoo") #run the pca code, specified the number of chromosomes
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

#PCA on Wild individuals
system("plink --bfile Wild_refiltered --allow-extra-chr --pca --chr-set 18 --out pca_wild") #run the pca code, specified the number of chromosomes
pca.data.wild <- read.table("pca_wild.eigenvec", header=FALSE) #load pca results from the eigenvec file
wild_origins <- read.csv("wild_origins.csv", header = TRUE) #read in origins file
pca.data.wild.origins <- left_join(pca.data.wild, wild_origins, by = "V2") #append origin data to pca data

ggplot(pca.data.wild.origins, aes(x=V3,y=V4, color = Cat.Group)) +  #plot with individual ID's and by origin
  geom_point() +
  geom_text(aes(label=V2), vjust=1, hjust=1, size=2) + #V2 is ID
  scale_color_manual(name = "Origin", values = c("Refuge" = "dodgerblue3", "Ranch" = "coral2")) +
  labs(x = "PC1", y = "PC2", title = "Wild PCA by Individual") +
  theme_minimal()

#PCA on all ocelots -- need to make new subset?
system("plink --bfile LP_Biallelic_AllChrom --pca --chr-set 18 --out pca_LP") #run the pca code, specified the number of chromosomes
pca.data.LP <- read.table("pca_LP.eigenvec", header=FALSE) #load pca results from the eigenvec file
ocelot.origins <- read.csv("LP_origins.csv", header = TRUE) #read in origins file
pca.data.ocelot.origins <- left_join(pca.data.LP, ocelot.origins, by = "V2") #append origin data to pca data

ggplot(pca.data.ocelot.origins, aes(x=V3,y=V4, color = Cat.Group)) +  #plot with individual ID's and by origin
  geom_point() +
  geom_text(aes(label=V2), vjust=1, hjust=1, size=2) + #V2 is ID
  scale_color_manual(name = "Origin", values = c("Refuge" = "purple4", "Ranch" = "orchid1", 
                                                 "Brazilian" = "dodgerblue3", "Generic" = "coral2")) +
  labs(x = "PC1", y = "PC2", title = "Ocelot PCA by Individual") +
  theme_minimal()

####KING Values####
#king-robust kingship estimator for zoo individuals
system("plink2 --bfile Zoo_refiltered --make-king-table --allow-extra-chr --chr-set 18 --out zoo_KING")
zoo.king.matrix <- read.table("zoo_KING.kin0", header = FALSE) #make king table into object
colnames(zoo.king.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
zoo.king.matrix <- zoo.king.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(zoo.king.matrix, "zoo_pairwise_kinship_refiltered.csv")

#heat map of kinship for zoo individuals -- gradient
library(reshape2)
melted_kin_zoo <- melt(zoo.king.matrix, measure.vars = "KINSHIP", variable.name = "KINSHIP")
head(melted_kin_zoo)
ggplot(data = melted_kin_zoo, aes(x=IID1, y=IID2, fill = value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.4, 0.3), space = "Lab",
                       name = "Kinship") +
  theme_minimal()
#heat map of kinship for zoo individuals -- binned
ggplot(data = melted_kin_zoo, aes(x=IID1, y=IID2, fill = value)) +
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

#value distribution plot
hist(zoo.king.matrix$KINSHIP)

#king-robust kingship estimator for wild individuals
system("plink2 --bfile Wild_refiltered --make-king-table --allow-extra-chr --chr-set 18 --out wild_KING")
king.wild.matrix <- read.table("wild_KING.kin0", header = FALSE) #make king table into object
colnames(king.wild.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
king.wild.matrix <- king.wild.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(king.wild.matrix, "wild_pairwise_kinship_refiltered.csv")

#heatmap for kinship of wild individuals -- gradient
melted_kin_wild <- melt(king.wild.matrix, measure.vars = "KINSHIP", variable.name = "KINSHIP")
head(melted_kin_wild)
ggplot(data = melted_kin_wild, aes(x=IID1, y=IID2, fill = value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                      midpoint = 0, limit = c(-0.4, 0.3), space = "Lab",
                      name = "Kinship") +
  theme_minimal()
#heat map of kinship for wild individuals -- binned
ggplot(data = melted_kin_wild, aes(x=IID1, y=IID2, fill = value)) +
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

#heat map of kinship for ranch individuals -- binned
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
####Inbreeding -- zoo, lumped####
#Inbreeding analyses need to be performed on a dataset that has already been pruned for LD
#Data pruned at an r2 of 0.5 -- moderate pruning that retained more snps
  #and changed nothing from the more stringent pruning

#run expected vs observed heterozygosity of the alleles
system("plink --bfile Zoo_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --het --out zoo_inbreeding_results")
zoo_inbreeding <- read.table("zoo_inbreeding_results.het", header = TRUE) #read in het file

#basic stats
summary(zoo_inbreeding)

#plot histogram of F values
hist(zoo_inbreeding$F,
     main = "Distribution of Inbreeding Coefficients in Zoo Individuals",
     xlab = "F Coefficient",
     col = "lightblue",
     breaks = 30)

#identify outliers -- potentially inbred individuals
    #common threshold is +/- 3 standard deviations
mean_f_zoo <-mean(zoo_inbreeding$F)
sd_f_zoo <- sd(zoo_inbreeding$F)
outliers_zoo <- subset(zoo_inbreeding, F > mean_f_zoo + 3*sd_f_zoo | F < mean_f_zoo - 3*sd_f_zoo)
print(outliers_zoo) #identify individuals with high inbreeding values

#write table with individual inbreeding values
inbreeding_table_zoo <- data.frame(
  ID = zoo_inbreeding$IID,
  Observed_Homozygosity = zoo_inbreeding$O.HOM.,
  Expected_Homozygosity = zoo_inbreeding$E.HOM.,
  Num_Genotyped_SNPs = zoo_inbreeding$N.NM.,
  Inbreeding_Coefficient_F = zoo_inbreeding$F
)
write.csv(inbreeding_table_zoo, "Zoo_Individual_Inbreeding.csv", row.names = FALSE)
####Inbreeding -- zoo, stock separate -- brazilian####
#run expected vs observed heterozygosity of the alleles
system("plink --bfile Brazilian_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --het --out brazilian_inbreeding_results")
brazilian_inbreeding <- read.table("brazilian_inbreeding_results.het", header = TRUE) #read in het file
#basic stats
summary(brazilian_inbreeding)
#plot histogram of F values
hist(brazilian_inbreeding$F,
     main = "Distribution of Inbreeding Coefficients in Brazilian Individuals",
     xlab = "F Coefficient",
     col = "lightblue",
     breaks = 30)
#identify outliers -- potentially inbred individuals
#common threshold is +/- 3 standard deviations
mean_f_braz <-mean(brazilian_inbreeding$F)
sd_f_braz <- sd(brazilian_inbreeding$F)
outliers_braz <- subset(brazilian_inbreeding, F > mean_f_braz + 3*sd_f_braz | F < mean_f_braz - 3*sd_f_braz)
print(outliers_braz) #identify individuals with high inbreeding values
#write table with individual inbreeding values
inbreeding_table_braz <- data.frame(
  ID = brazilian_inbreeding$IID,
  Observed_Homozygosity = brazilian_inbreeding$O.HOM.,
  Expected_Homozygosity = brazilian_inbreeding$E.HOM.,
  Num_Genotyped_SNPs = brazilian_inbreeding$N.NM.,
  Inbreeding_Coefficient_F = brazilian_inbreeding$F
)
write.csv(inbreeding_table_braz, "Brazilian_Individual_Inbreeding.csv", row.names = FALSE)

####Inbreeding -- zoo, stock separate -- generic####
#run expected vs observed heterozygosity of the alleles
system("plink --bfile Generic_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --het --out generic_inbreeding_results")
gen_inbreeding <- read.table("generic_inbreeding_results.het", header = TRUE) #read in het file
#basic stats
summary(gen_inbreeding)
#plot histogram of F values
hist(gen_inbreeding$F,
     main = "Distribution of Inbreeding Coefficients in Generic Individuals",
     xlab = "F Coefficient",
     col = "lightblue",
     breaks = 30)
#identify outliers -- potentially inbred individuals
mean_f_gen <-mean(gen_inbreeding$F)
sd_f_gen <- sd(gen_inbreeding$F)
outliers_gen <- subset(gen_inbreeding, F > mean_f_gen + 3*sd_f_gen | F < mean_f_gen - 3*sd_f_gen)
print(outliers_gen) #identify individuals with high inbreeding values
#write table with individual inbreeding values
inbreeding_table_gen <- data.frame(
  ID = gen_inbreeding$IID,
  Observed_Homozygosity = gen_inbreeding$O.HOM.,
  Expected_Homozygosity = gen_inbreeding$E.HOM.,
  Num_Genotyped_SNPs = gen_inbreeding$N.NM.,
  Inbreeding_Coefficient_F = gen_inbreeding$F
)
write.csv(inbreeding_table_gen, "Generic_Individual_Inbreeding.csv", row.names = FALSE)

####Inbreeding -- wild, pops lumped####
#run expected vs observed heterozygosity of the alleles
system("plink --bfile Wild_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --het --out wild_inbreeding_results")
wild_inbreeding <- read.table("wild_inbreeding_results.het", header = TRUE) #read in het file

#basic stats
summary(wild_inbreeding)

#plot histogram of F values
hist(wild_inbreeding$F,
     main = "Distribution of Inbreeding Coefficients in Wild Individuals",
     xlab = "F Coefficient",
     col = "lightblue",
     breaks = 30)

#identify outliers -- potentially inbred individuals
#common threshold is +/- 3 standard deviations
mean_f_wild <-mean(wild_inbreeding$F)
sd_f_wild <- sd(wild_inbreeding$F)
outliers_wild <- subset(wild_inbreeding, F > mean_f_wild + 3*sd_f_wild | F < mean_f_wild - 3*sd_f_wild)
print(outliers_wild) #identify individuals with high inbreeding values

#write table with individual inbreeding values
inbreeding_table_wild <- data.frame(
  ID = wild_inbreeding$IID,
  Observed_Homozygosity = wild_inbreeding$O.HOM.,
  Expected_Homozygosity = wild_inbreeding$E.HOM.,
  Num_Genotyped_SNPs = wild_inbreeding$N.NM.,
  Inbreeding_Coefficient_F = wild_inbreeding$F
)
write.csv(inbreeding_table_wild, "Wild_Individual_Inbreeding.csv", row.names = FALSE)






####Inbreeding -- wild, pops separate, ranch####
#run expected vs observed heterozygosity of the 
#currently this code is rewritten for the test of pruning the populations seperately rather than subsetting after pruning
system("plink --bfile Ranch_isolatetest_LDpruned_0.5 --allow-extra-chr --chr-set 18 --het --out ranch_inbreeding_results_isolatetest")
ra_inbreeding <- read.table("ranch_inbreeding_results_isolatetest.het", header = TRUE) #read in het file
#basic stats
summary(ra_inbreeding)
#plot histogram of F values
hist(ra_inbreeding$F,
     main = "Distribution of Inbreeding Coefficients in Ranch Individuals",
     xlab = "F Coefficient",
     col = "lightblue",
     breaks = 30)
#identify outliers -- potentially inbred individuals
mean_f_ra <-mean(ra_inbreeding$F)
sd_f_ra <- sd(ra_inbreeding$F)
outliers_ra <- subset(ra_inbreeding, F > mean_f_ra + 3*sd_f_ra | F < mean_f_ra - 3*sd_f_ra)
print(outliers_ra) #identify individuals with high inbreeding values
#write table with individual inbreeding values
inbreeding_table_ra <- data.frame(
  ID = ra_inbreeding$IID,
  Observed_Homozygosity = ra_inbreeding$O.HOM.,
  Expected_Homozygosity = ra_inbreeding$E.HOM.,
  Num_Genotyped_SNPs = ra_inbreeding$N.NM.,
  Inbreeding_Coefficient_F = ra_inbreeding$F
)
write.csv(inbreeding_table_ra, "Ranch_Individual_Inbreeding_isolatetest.csv", row.names = FALSE)

####Inbreeding -- wild, pops separate, refuge####
#run expected vs observed heterozygosity of the alleles
system("plink --bfile Refuge_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --het --out refuge_inbreeding_results")
re_inbreeding <- read.table("refuge_inbreeding_results.het", header = TRUE) #read in het file
#basic stats
summary(re_inbreeding)
#plot histogram of F values
hist(re_inbreeding$F,
     main = "Distribution of Inbreeding Coefficients in Refuge Individuals",
     xlab = "F Coefficient",
     col = "lightblue",
     breaks = 30)
#identify outliers -- potentially inbred individuals
mean_f_re <-mean(re_inbreeding$F)
sd_f_re <- sd(re_inbreeding$F)
outliers_re <- subset(re_inbreeding, F > mean_f_re + 3*sd_f_re | F < mean_f_re - 3*sd_f_re)
print(outliers_re) #identify individuals with high inbreeding values
#write table with individual inbreeding values
inbreeding_table_re <- data.frame(
  ID = re_inbreeding$IID,
  Observed_Homozygosity = re_inbreeding$O.HOM.,
  Expected_Homozygosity = re_inbreeding$E.HOM.,
  Num_Genotyped_SNPs = re_inbreeding$N.NM.,
  Inbreeding_Coefficient_F = re_inbreeding$F
)
write.csv(inbreeding_table_re, "Refuge_Individual_Inbreeding.csv", row.names = FALSE)


####ROH -- not working####
##zoo -- all defaults in plink
system("plink --bfile Zoo_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --out zoo_roh_base")
zoo_base_roh_results <- read.table("zoo_roh_base.hom.indiv", header = T)
    #located 2 roh in 2 individuals
##zoo --less strigent search parameters, allowing for smaller ROH
system("plink --bfile Zoo_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --homozyg-het 0 --homozyg-snp 50 --out zoo_roh_t1")
zoo_t1_roh_results <- read.table("zoo_roh_t1.hom.indiv", header = T)
    #no roh found
##zoo -- following protocol of alemu et al
system("plink --bfile Zoo_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --homozyg-het 0 --homozyg-snp 50 --homozyg-kb 2000 --homozyg-density 100 --homozyg-gap 500 --out zoo_roh_t2")
zoo_t2_roh_results <- read.table("zoo_roh_t2.hom.indiv", header = T)
    #no roh found
##zoo -- no maf filtering, default plink
system("plink --bfile Zoo_roh_filter --allow-extra-chr --chr-set 18 --homozyg --out zoo_roh_nomaf_base")
zoo_nomaf_base_roh_results <- read.table("zoo_roh_nomaf_base.hom.indiv", header = T)

##wild -- all defaults in plink
system("plink --bfile Wild_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --out wild_roh_base")
wild_base_roh_results <- read.table("wild_roh_base.hom.indiv", header = T)
    #no roh found
##wild --less strigent search parameters, allowing for smaller ROH
system("plink --bfile Wild_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --homozyg-het 0 --homozyg-snp 50 --out wild_roh_t1")
wild_t1_roh_results <- read.table("wild_roh_t1.hom.indiv", header = T)
    #no roh found
##wild -- following protocol of alemu et al
system("plink --bfile Wild_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --homozyg-het 0 --homozyg-snp 50 --homozyg-kb 2000 --homozyg-density 100 --homozyg-gap 500 --out wild_roh_t2")
wild_t2_roh_results <- read.table("wild_roh_t2.hom.indiv", header = T)
    #no roh found
##wild -- loosest parameters
system("plink --bfile Wild_LDpruned_0.5_RF --allow-extra-chr --chr-set 18 --homozyg --homozyg-het 1 --homozyg-snp 25 --homozyg-kb 1000 --homozyg-density 200 --homozyg-gap 1000 --homozyg-window-het 1 --out wild_roh_t3")
wild_t3_roh_results <- read.table("wild_roh_t3.hom.indiv", header = T)
##wild -- no maf filtering, no ld pruning default plink
system("plink --bfile Wild_roh_filter --allow-extra-chr --chr-set 18 --homozyg --out wild_roh_nomaf_base")
wild_nomaf_base_roh_results <- read.table("wild_roh_nomaf_base.hom.indiv", header = T)
##wild -- no maf no ld pruning, loosest parameters
system("plink --bfile Wild_roh_filter --allow-extra-chr --chr-set 18 --homozyg --homozyg-het 1 --homozyg-snp 25 --homozyg-kb 1000 --homozyg-density 200 --homozyg-gap 1000 --homozyg-window-het 1 --out wild_roh_nomaf_t2")
wild_nomaf_t2_roh_results <- read.table("wild_roh_nomaf_t2.hom.indiv", header = T)
##wild -- meyermans et al parameters, no maf no ld
system("plink --bfile Wild_roh_filter --allow-extra-chr --chr-set 18 --homozyg --homozyg-kb 500 --homozyg-het 1 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-density 40 --out wild_meyer_roh")
wild_meyer_roh_results <- read.table("wild_meyer_roh.hom.indiv", header = T)


##########ROH WORKING!!!!!!#####
##wild -- M.Smith params
system("plink --bfile Wild_roh_filter --allow-extra-chr --chr-set 18 --homozyg --homozyg-density 50 --homozyg-gap 1000
       --homozyg-kb 300 --homozyg-snp 50 --homozyg-window-het 5 --homozyg-window-missing 5 --homozyg-window-snp 50
       --homozyg-window-threshold 0.05 --out wild_roh_smith")
wild_smith_results <- read.table("wild_roh_smith.hom.indiv", header = T)
##wild -- Saremi params -- more relaxed
system("plink --bfile Wild_roh_filter --allow-extra-chr --chr-set 18 --homozyg --homozyg-density 50 --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 20 --homozyg-window-missing 20 --homozyg-window-snp 50
       --homozyg-window-threshold 0.02 --homozyg-het 750 --out wild_roh_saremi")
wild_saremi_results <- read.table("wild_roh_saremi.hom.indiv", header = T) #more roh found, but increased average length while reducing number of segments
#wild -- thinned data saremi params  -- not working as well, 90% genome in ROH?
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-density 50 --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 20 --homozyg-window-missing 20 --homozyg-window-snp 50
       --homozyg-window-threshold 0.02 --homozyg-het 750 --out wild_roh_saremi_thin")
wild_saremi_thin_results <- read.table("wild_roh_saremi_thin.hom.indiv", header = T)
#wild -- thinned meyermans params -- best of the thinned, but not great
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-kb 500 --homozyg-het 1 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-density 40 --out wild_meyer_roh_thin")
wild_meyer_roh_thin_results <- read.table("wild_meyer_roh_thin.hom.indiv", header = T)
#wild -- thinned, stricter saremi params -- not working as well, 90% genome in ROH?
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-density 50 --homozyg-gap 1000 --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 10 --homozyg-window-missing 20 --homozyg-window-snp 50
       --homozyg-window-threshold 0.05 --homozyg-het 750 --out wild_roh_thin_strict1")
wild_roh_thin_strict1_results <- read.table("wild_roh_thin_strict1.hom.indiv", header = T)


##zoo -- M.Smith params
system("plink --bfile Zoo_roh_filter --allow-extra-chr --chr-set 18 --homozyg --homozyg-density 50 --homozyg-gap 1000
       --homozyg-kb 300 --homozyg-snp 50 --homozyg-window-het 5 --homozyg-window-missing 5 --homozyg-window-snp 50
       --homozyg-window-threshold 0.05 --out zoo_roh_smith")
zoo_smith_results <- read.table("zoo_roh_smith.hom.indiv", header = T)
#zoo --saremi params
system("plink --bfile zoo_roh_filter --allow-extra-chr --chr-set 18 --homozyg --homozyg-density 50 --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 20 --homozyg-window-missing 20 --homozyg-window-snp 50
       --homozyg-window-threshold 0.02 --homozyg-het 750 --out zoo_roh_saremi")
zoo_saremi_results <- read.table("zoo_roh_saremi.hom.indiv", header = T)


####ROH -- single chromosome test, filtering -- not working####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/SNP_chrom68")
###input chrom file
system("plink --vcf 24041DeY-snp_filter_dp_gq_NC_058368.1.vcf.gz --keep-allele-order --allow-extra-chr --make-bed --out SNP_dp_gq_chrom1")
###subset individuals
#create zoo file
system("plink --bfile SNP_dp_gq_chrom1 --keep pop_subset_zoo.txt --keep-allele-order --allow-extra-chr --make-bed --out Zoo_BaseFilter_chrom1")
#create wild ocelots, no duplicates file
system("plink --bfile SNP_dp_gq_chrom1 --keep pop_subset_wild_remove_dups.txt --keep-allele-order --allow-extra-chr --make-bed --out Wild_BaseFilter_chrom1")
###filter for miss, biallelic, and hwe
#zoo
#filter missingness and hwe
system("plink --bfile Zoo_BaseFilter_chrom1 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Zoo_miss90_hwe_chrom1")
#filter biallelic
system("plink --bfile Zoo_miss90_hwe_chrom1 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Zoo_roh_chrom1")

#Wild
#filter missingness and hwe
system("plink --bfile Wild_BaseFilter_chrom1 --allow-extra-chr --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out Wild_miss90_hwe_chrom1")
#filter biallelic
system("plink --bfile Wild_miss90_hwe_chrom1 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out Wild_roh_chrom1")


###ROH
#default plink zoo
system("plink --bfile Zoo_roh_chrom1 --allow-extra-chr --homozyg --out zoo_roh_chr1_default")
zoo_chr1_default_results <- read.table("zoo_roh_chr1_default.hom.indiv", header = T)
#meyermans et al zoo
system("plink --bfile Zoo_roh_chrom1 --allow-extra-chr --homozyg --homozyg-kb 500 --homozyg-het 1 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-density 40 --out zoo_meyer_chr1")
zoo_meyer_chr1_results <- read.table("zoo_meyer_chr1.hom.indiv", header = T)
#default plink wild
system("plink --bfile Wild_roh_chrom1 --allow-extra-chr --homozyg --out wild_roh_chr1_default")
wild_chr1_default_results <- read.table("wild_roh_chr1_default.hom.indiv", header = T)
#meyermans et al wild
system("plink --bfile Wild_roh_chrom1 --allow-extra-chr --homozyg --homozyg-kb 500 --homozyg-het 1 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-density 40 --out wild_meyer_chr1")
wild_meyer_chr1_results <- read.table("wild_meyer_chr1.hom.indiv", header = T)
################################################################################
##DAPC analysis using adegenet package -- 4/17/25
install.packages("adegenet")
library(adegenet)
library(poppr)
install.packages("vcfR")
library(vcfR)
install.packages("promises", version = ">= 1.3.2") #got error with incompatible versions, updating this
library(promises)
####dataset merging####
#for the species wide comparison, I am combining the separately filtered zoo and wild files
  #rather than performing the filtering steps on a combined file. This is to retain informative
  #snps that may be removed during the maf filtering steps as the populations are so different
  #should be noted that the zoo individuals have twice as many variants -- should consider for bias

#giving snps unique ID to facilitate merging and avoid duplicate errors 
#not working -- fix for the combined analysis
system("plink2 --bfile Zoo_Refiltered --chr-set 18 --allow-extra-chr --set-all-var-ids @:#:$1,$2 --new-id-max-allele-len 20 --make-bed --out Zoo_refiltered_uniqueID_short") #set variant ID -- shorter?
system("plink --bfile Zoo_refiltered_uniqueID --bmerge Wild_refiltered_uniqueID --chr-set 18 --allow-extra-chr --memory 25000 --make-bed --out LEPA_merged") #merging the datasets
    #using the unique ID file generated earlier to avoid a merging error due to same variant names

#making vcf's out of plink files
system("plink2 --bfile Zoo_refiltered --chr-set 18 --allow-extra-chr --export vcf bgz --out Zoo_VCF_DAPC")
system("plink2 --bfile Wild_refiltered --chr-set 18 --allow-extra-chr --export vcf bgz --out Wild_VCF_DAPC")
##switching working directories
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Adegenet_analyses")
####vcf to genlight objects -- both pops separate and merged####
vcf <- read.vcfR("Zoo_VCF_DAPC.vcf.gz", verbose = TRUE)
zoo_genlight <- vcfR2genlight(vcf)
zoo_genlight

w_vcf <- read.vcfR("Wild_VCF_DAPC.vcf.gz", verbose = TRUE)
wild_genlight <- vcfR2genlight(w_vcf)
wild_genlight

####finding clusters####
zoo_clust <-find.clusters(zoo_genlight, max.n.clust = 40)
