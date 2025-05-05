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
#create full LEPA file
system("plink --bfile SNP_dp_gq_AllChrom_AllInd --keep pop_subset_ocelot.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out LEPA_BaseFilter_Allchrom")
#create zoo file
system("plink --bfile SNP_dp_gq_AllChrom_AllInd --keep pop_subset_zoo.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Zoo_BaseFilter_Allchrom")
#create wild ocelots, no duplicates file
system("plink --bfile SNP_dp_gq_AllChrom_AllInd --keep pop_subset_wild_remove_dups.txt --chr-set 18 --keep-allele-order --allow-extra-chr --make-bed --out Wild_BaseFilter_AllChrom")

####Filtering -- MAF 0.05 MISS 90 Bialleic####
##LEPA filtering -- wild and zoo together
#filter MAF and missingness and hwe
system("plink --bfile LEPA_BaseFilter_Allchrom --chr-set 18 --allow-extra-chr --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out LEPA_maf05_miss90_hwe")
#filter biallelic
system("plink --bfile LEPA_maf05_miss90_hwe --chr-set 18 --allow-extra-chr --keep-allele-order --biallelic-only --make-bed --out LEPA_refiltered")

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
#export as vcf
system("plink2 --bfile wild_roh_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out wild_thinned_roh")

#thinning for dapc
system("plink --bfile Wild_refiltered --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out wild_dapc_thin")
system("plink --bfile Zoo_refiltered --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out zoo_dapc_thin")
system("plink --bfile LEPA_refiltered --chr-set 18 --allow-extra-chr --bp-space 1000 --make-bed --out lepa_dapc_thin")

system("plink2 --bfile wild_dapc_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out wild_dapc_thin")
system("plink2 --bfile zoo_dapc_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out zoo_dapc_thin")
system("plink2 --bfile lepa_dapc_thin --chr-set 18 --allow-extra-chr --export vcf bgz --out lepa_dapc_thin")

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

#wild --thinned, second try of stricter saremi params for thinned data -- still too permissive
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 300 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-window-missing 20 --homozyg-window-snp 50
       --homozyg-window-threshold 0.02 --out wild_roh_thin_test1")
wild_roh_thin_test1_results <- read.table("wild_roh_thin_test1.hom.indiv", header = T)

#wild --thinned, third try of stricter saremi params for thinned data ---BEST SO FAR
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-het 20 --homozyg-window-missing 20 --homozyg-window-snp 100
       --homozyg-window-threshold 0.02 --out wild_roh_thin_test2")
wild_roh_thin_test2_results <- read.table("wild_roh_thin_test2.hom.indiv", header = T)

#wild --thinned, third adjustment test -- no more than 20% variation across individuals
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 300 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-het 20 --homozyg-window-missing 10 --homozyg-window-snp 100
       --homozyg-window-threshold 0.05 --out wild_roh_thin_test3")
wild_roh_thin_test3_results <- read.table("wild_roh_thin_test3.hom.indiv", header = T)

#wild --thinned -- composite params test -- combined saremi, test3, and meyermans
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg  --homozyg-kb 300  --homozyg-snp 50 --homozyg-density 40 
       --homozyg-gap 500 --homozyg-het 3 --homozyg-window-snp 100 --homozyg-window-het 3  --homozyg-window-missing 5 
       --homozyg-window-threshold 0.05 --out wild_thin_composite1")
wild_roh_comp_results <- read.table("wild_thin_composite1.hom.indiv", header = T)

#wild --thinned -- composite params test 2 -- combined saremi, test3, and meyermans --adjusted homozyg-het and gap to potentially reduce breaks
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg  --homozyg-kb 300  --homozyg-snp 50 --homozyg-density 40 --homozyg-gap 1000 --homozyg-het 10 --homozyg-window-snp 100 --homozyg-window-het 3  --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out wild_thin_composite2")
wild_roh_comp2_results <- read.table("wild_thin_composite2.hom.indiv", header = T)

#wild --thinned, fourth adjustment test -- reduced homozyg-het flag
system("plink --bfile wild_roh_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 300 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-het 3 --homozyg-window-missing 10 --homozyg-window-snp 100
       --homozyg-window-threshold 0.05 --out wild_roh_thin_test4")
wild_roh_thin_test4_results <- read.table("wild_roh_thin_test4.hom.indiv", header = T)



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

#zoo --thinned, third try of stricter saremi params for thinned data ---BEST SO FAR
system("plink --bfile zoo_dapc_thin --allow-extra-chr --chr-set 18 --homozyg --homozyg-gap 1000
       --homozyg-kb 500 --homozyg-snp 50 --homozyg-window-het 3 --homozyg-het 20 --homozyg-window-missing 20 --homozyg-window-snp 100
       --homozyg-window-threshold 0.02 --out zoo_roh_thin_test2")
zoo_roh_thin_test2_results <- read.table("zoo_roh_thin_test2.hom.indiv", header = T)


####Comparing parameter effects of ROH####
#wild -- populate dataframe with ID's and percent genome in ROH for each parameter run
roh.df <- wild_smith_results$IID #populate ID's
roh.df <- as.data.frame(roh.df) #create the data frame
roh.df$Smith_full <- ((wild_smith_results$KB*1000)/2425730029)*100 #unthinned smith param
roh.df$Saremi_full <- ((wild_saremi_results$KB*1000)/2425730029)*100 #unthinned saremi param
roh.df$Saremi_thin <- ((wild_saremi_thin_results$KB*1000)/2425730029)*100 #thinned saremi param
roh.df$Meyer_thin <- ((wild_meyer_roh_thin_results$KB*1000)/2425730029)*100 #thinned meyermans param
roh.df$Strict_thin <- ((wild_roh_thin_strict1_results$KB*1000)/2425730029)*100
roh.df$Thin_test1 <- ((wild_roh_thin_test1_results$KB*1000)/2425730029)*100
roh.df$Thin_test2 <- ((wild_roh_thin_test2_results$KB*1000)/2425730029)*100
roh.df$Thin_test3 <- ((wild_roh_thin_test3_results$KB*1000)/2425730029)*100
roh.df$Thin_test4 <- ((wild_roh_thin_test4_results$KB*1000)/2425730029)*100
roh.df$Composite_thin <- ((wild_roh_comp_results$KB*1000)/2425730029)*100
roh.df$Composite_thin2 <- ((wild_roh_comp2_results$KB*1000)/2425730029)*100
##writing as csv
write.csv(roh.df, "wild_roh_percentgenome.csv")

#zoo -- populate dataframe with ID's and percent genome for each parameter test
roh.df.z <- zoo_saremi_results$IID #populate ID's
roh.df.z <- as.data.frame(roh.df.z) #create the data frame
roh.df.z$Saremi_full <- ((zoo_saremi_results$KB*1000)/2425730029)*100 #unthinned saremi param
roh.df.z$Thin_test2 <- ((zoo_roh_thin_test2_results$KB*1000)/2425730029)*100
write.csv(roh.df.z, "zoo_roh_percentgenome.csv")

####Visualizing ROH --karyotype plot -- wild####
#install packages
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
install.packages("GenomicAlignments")
library(karyoploteR)
library(regioneR)
library(GenomicRanges)
library(data.table)
library(IRanges)
library(GenomicAlignments)

#read in hom file
w_hom <- read.table("wild_thin_composite2.hom", header = TRUE)
#make data frame for plotting
w_roh_plot <- data.frame(
  chr = paste0("chr", w_hom$CHR),  # Add 'chr' prefix if needed
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
w_roh_plot[w_roh_plot == "chrNC_058368.1"] <- "chr1"
w_roh_plot[w_roh_plot == "chrNC_058369.1"] <- "chr2"
w_roh_plot[w_roh_plot == "chrNC_058370.1"] <- "chr3"
w_roh_plot[w_roh_plot == "chrNC_058371.1"] <- "chr4"
w_roh_plot[w_roh_plot == "chrNC_058372.1"] <- "chr5"
w_roh_plot[w_roh_plot == "chrNC_058373.1"] <- "chr6"
w_roh_plot[w_roh_plot == "chrNC_058374.1"] <- "chr7"
w_roh_plot[w_roh_plot == "chrNC_058375.1"] <- "chr8"
w_roh_plot[w_roh_plot == "chrNC_058376.1"] <- "chr9"
w_roh_plot[w_roh_plot == "chrNC_058377.1"] <- "chr10"
w_roh_plot[w_roh_plot == "chrNC_058378.1"] <- "chr11"
w_roh_plot[w_roh_plot == "chrNC_058379.1"] <- "chr12"
w_roh_plot[w_roh_plot == "chrNC_058380.1"] <- "chr13"
w_roh_plot[w_roh_plot == "chrNC_058381.1"] <- "chr14"
w_roh_plot[w_roh_plot == "chrNC_058382.1"] <- "chr15"
w_roh_plot[w_roh_plot == "chrNC_058383.1"] <- "chr16"
w_roh_plot[w_roh_plot == "chrNC_058384.1"] <- "chr17"
w_roh_plot[w_roh_plot == "chrNC_058385.1"] <- "chr18"
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
    kpPlotLines(w_kp, data = smoothed_gr, y = mcols(smoothed_gr)$score, col = "red", r0 = 0.5, r1 = 0.8)
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

###ROH overlap plot

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
####switching working directories####
setwd("C:/Users/kutab016/Documents/TB_Files/1_Thesis/3_Data/6_Cleaned Data/Genomics/1_WorkingGenomicsFiles/Adegenet_analyses")
####vcf to genlight objects -- both pops separate and merged####
#using thinned dataset as unthinned is too large
vcf <- read.vcfR("zoo_dapc_thin.vcf.gz", verbose = TRUE)
zoo_genlight <- vcfR2genlight(vcf)
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

####finding clusters####
zoo_clust <-find.clusters(zoo_genlight, max.n.clust = 40)
wild_clust <-find.clusters(wild_genlight, max.n.clust = 40)
lepa_clust <-find.clusters(lepa_genlight, max.n.clust = 40)

#looking at outputs
names(zoo_clust)
names(wild_clust)
zoo_clust$size
wild_clust$size
#looking at groups
table(pop(wild_genlight), wild_clust$grp)
table(pop(zoo_genlight), zoo_clust$grp)
table(pop(lepa_genlight), lepa_clust$grp) #gave weird groupings? with 4 clusters selected

####making DAPC####
w_dapc <- dapc(wild_genlight, wild_clust$grp)
w_dapc

scatter(w_dapc, posi.da = "bottomright", bg = "white")

l_dapc <- dapc(lepa_genlight, lepa_clust$grp)
scatter(l_dapc)
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(l_dapc, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")
