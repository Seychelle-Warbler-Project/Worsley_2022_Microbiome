##### AlphaDiversityFungal.R #####

# This R script contains the code necessary to 1) calculate the core fungal microbiome 
# 2) compare alpha diversity for bacteria and fungi in the same sample 
# 3) reproduce analyses investigating the association between environmental/host variables and gut microbiome alpha diversity, and 
# 4) investigate the association between fungal alpha diversity and host survival
# It uses the filtered ASV, taxonomy and metadata files generated in FungalSequences.R

#load packages
  library(ggplot2)
  library(ggpubr)
  library(phyloseq)
  library(ape)
  library(plyr)
  library(RColorBrewer)
  library(forcats)
  library(microbiome)
  library(arm)
  library(car)
  library(GGally)
  library(DHARMa)
  library(MuMIn)
  library(MASS)
  library(patchwork)


##################################
#### Generate Phyloseq object ####
##################################

#ASV table
  asv_table <- read.csv ("FilteredASVTable.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) 
  # should be 2540 features, 265 samples (unrarefied reads)
  # Note that all ASVs unassigned to phylum have been removed as well as potential contaminants
  # control samples have been filtered
  # faecal samples with <5000 reads have been removed
  # ASVs with <50 reads in total across all samples removed
  # duplicate extractions still present
  asv_table <- as.matrix (asv_table) #make into a matrix

#taxonomy
  taxonomy <- read.csv ("FilteredTaxonomyTable.csv", row.names=1)
  str (taxonomy) #should be 2540 observations, 7 taxonomic groupings, feature names as row names.
  taxonomy <- as.matrix (taxonomy)

#load filtered metadata with sample names as rownames- 265 observations
  metadata<-read.csv("FilteredMetadata.csv", row.names = 1)
  head(metadata)  
  str(metadata)

#import all as phyloseq objects
  ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
  TAX <- tax_table(taxonomy)
  META <- sample_data(metadata)
  head(META)

#check that the ASV and sample names are consistent across objects (e.g have dashes been replaced with dots?)
  str(taxa_names(TAX))
  str(taxa_names(ASV))

  str(sample_names(ASV))
  str(sample_names(META))

#### MERGE INTO PHYLOSEQ OBJECT ####
  physeq <- phyloseq(ASV, TAX, META)
  physeq #2540 taxa 265 samples
  
  sample_data(physeq)$SampleReads<- sample_sums(physeq)
  
######################
#### Rarefy Reads ####
######################
  
### rarefy to 5000 and set seed before rarefying (89754), so that results are reproducible  
  
  physeqRare<-rarefy_even_depth(physeq, 5000, rngseed = 89754)
  
  # 10 ASVs removed after subsampling- leaves 2530 ASVs across 265 samples
  physeqRare
  sample_sums(physeqRare)
  
  
### REMOVE EXTRACTION DUPLICATES ###
# retain the one with the highest read count but exclude both duplicates of 408 as bird sampled again in 2020 (tube 490)
  
  TubesToRemove<- c("FUN370R", "FUN387", "FUN408", "FUN408R", "FUN424", "FUN453")
  
  physeqRare2<- prune_samples(!sample_data(physeqRare)$TubeNoFungal %in% TubesToRemove, physeqRare)
  physeqRare2 # 259 samples remain after duplicates removed  
  

#########################################
##### PLOT RELATIVE ABUNDANCES ##########
#########################################
  
#check ASV number per sample in rarefied samples
  ASVno<-data.frame(estimate_richness(physeqRare2, split=TRUE, measures= c("Observed")))
  head(ASVno)
  summary(ASVno$Observed) #mean ASV no= 49.49
  sd(ASVno$Observed) #sd 26.02
  
#transform counts to relative abundances
  relabund <- transform_sample_counts(physeqRare2, function(x){(x / sum(x))*100})
  relabundTable<-otu_table(relabund)
  head(relabundTable)
  colSums(relabundTable) 
  
#### PHYLUM LEVEL ####
  
  ps.phylum<- tax_glom(relabund, taxrank="Phylum", NArm=FALSE)
  ps.phylum #detected 7 phyla
  plot_bar(ps.phylum, fill="Phylum")
  
  dat <- psmelt(ps.phylum) #create a dataframe of agglomerated data at phylum level.
  str(dat)
  
#barplot:relative abundance of each phylum per sample.
  #order the bars by the abundance of Ascomycota
  
  unique(dat$Phylum)
  dat$Phylum<-factor(dat$Phylum, c("Ascomycota", "Basidiomycota", "Mortierellomycota", "Mucoromycota",
                                   "Basidiobolomycota", "Rozellomycota", "Chytridiomycota"))
  dfAsco<-dat[dat$Phylum=="Ascomycota",]
  str(dfAsco)
  dfAsco<-dfAsco[order(-dfAsco$Abundance),]
  sampleOrder<-dfAsco$Sample
  head(sampleOrder)
  sampleOrder[1:37]
  
  dat$Sample<-factor(dat$Sample, sampleOrder)
  str(dat)
  pal.phy<- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F")

  
  phylumPlot <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Phylum))
  phylumPlot<- phylumPlot + geom_bar(aes(), stat="identity", position="stack") + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values=pal.phy)  +
    guides(fill=guide_legend(nrow=7)) + ylab("Relative abundance (%)")+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))+ ylab("Relative Abundance (%)")+
    theme(axis.title = element_text(size=18), axis.text = element_text(size=16),
          legend.text = element_text(size=12, face="italic"), legend.title = element_text(size=12))
  


#### FAMILY LEVEL ####
  
  ps.fam<- tax_glom(relabund, taxrank="Family", NArm=FALSE)
  
  datFam <- psmelt(ps.fam) #create a dataframe of agglomerated data at phylum level.
  str(datFam)
  mediansFam <- ddply(datFam, ~Family, function(x) c(median=median(x$Abundance))) #find the median count per phylum
  mediansFam[order(mediansFam$median),]
  
  
  # find families whose median relabundance is less than 1%
  other <- mediansFam[mediansFam$median <= 0.01,]$Family
  str(other)
  
  # change their name to "Other" to make plot less complex
  datFam[datFam$Family %in% other,]$Family <- 'Other'

  
#barplot:relative abundance of each family per sample.
  #order the bars by the abundance of Cladosporiaceae

  unique(datFam$Family)
  datFam$Family<-factor(datFam$Family, c( "Cladosporiaceae","Dothideales_fam_Incertae_sedis", "Aspergillaceae", "Microstromatales_fam_Incertae_sedis", "Bulleribasidiaceae",
                                          "Bionectriaceae", "Mycosphaerellaceae", "Hypocreales_fam_Incertae_sedis", "Didymellaceae", "Nectriaceae",
                                          "Neodevriesiaceae", "Phaeosphaeriaceae", "Teratosphaeriaceae", "Other"))
  dfClado<-datFam[datFam$Family=="Cladosporiaceae",]
  dfClado<-dfClado[order(-dfClado$Abundance),]
  sampleOrderFam<-dfClado$Sample
  head(sampleOrderFam)
  
  datFam$Sample<-factor(datFam$Sample, sampleOrderFam)
  str(datFam)
  
  #make palette
  brewer.pal(n = 12, name = "Paired")
  
  pal.fam<- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
           "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "grey75", "grey44")

  
  plotFam <- ggplot(data=datFam, aes(x=Sample, y=Abundance, fill=Family))
  plotFam<- plotFam + geom_bar(aes(), stat="identity", position="stack")  + 
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    guides(fill=guide_legend(nrow=14)) + ylab("Relative abundance (%)")+
    scale_fill_manual(values=pal.fam)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))+
    theme(axis.title = element_text(size=18), axis.text = element_text(size=16),
          legend.text = element_text(size=12, face="italic"), legend.title = element_text(size=12))
  
  
#plot
  
  (phylumPlot/plotFam)
  

#summarise the numbers of phyla, class etc... identified
  ps.phylum #7 
  ps.Class<- tax_glom(relabund, taxrank="Class", NArm=FALSE)
  ps.Class # 39
  ps.fam # 249
  ps.Gen<- tax_glom(relabund, taxrank="Genus", NArm=FALSE)
  ps.Gen #475


########################  
#### CORE MYCOBIOME ####
########################
  
#https://web.stanford.edu/class/bios221/Pune/Labs/Lab_microbiome/Lab_microbiome.html 
  
# Transform to compositional abundances
  pseq.rel <- microbiome::transform(physeqRare2, "compositional")
  otu_table(pseq.rel)
  sample_sums(pseq.rel)
  pseq.rel #2530 ASVs, 259 samples
  
# Calculate prevalences for all taxonomic groups
  head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))
  
# Define core ASVs (>0.1% relative abundance in >50% of the samples)
  pseq.core <- core(pseq.rel, detection = 0.1/100, prevalence = 50/100)
  taxonomyCoreAsv<- data.frame(tax_table(pseq.core))
  taxonomyCoreAsv$ASV<- row.names(taxonomyCoreAsv)

  abundCore<- data.frame(otu_table(pseq.core)) #How abundant are they across samples?
  str(abundCore)
  head(abundCore)
  MeanAbundance<- data.frame(100*(apply(abundCore[ ,1:259],1,mean)))#mean abundance
  colnames(MeanAbundance)<- "MeanAbundance"
  SDAbundances<- data.frame(100*(apply(abundCore[ ,1:259],1,sd))) #standard deviation of abundance across samples
  colnames(SDAbundances)<- "SD"
  PrevCoreASV<- data.frame(prevalence(pseq.core, detection = 0.1/100, sort = TRUE)) #prevelance of each core ASV- proportion of samples they are in (out of 1).
  colnames(PrevCoreASV)<-"Prevalence"
  
  taxonomyCoreAsv$MeanAbundance<- MeanAbundance$MeanAbundance  
  taxonomyCoreAsv$SDAbundance<- SDAbundances$SD
  taxonomyCoreAsv$Prevalence<- PrevCoreASV$Prevalence
  
  #write.csv(taxonomyCoreAsv, "CoreASVsMycobiome.csv")
  
# phylum level
  pseq.relphylum <- microbiome::transform(ps.phylum, "compositional")
  otu_table(pseq.relphylum)
  sample_sums(pseq.relphylum)
  # Pick the core (>0.1% relative abundance in >50% of the samples)
  pseq.corePhylum <- core(pseq.relphylum, detection = 0.1/100, prevalence = 50/100)
  tax_table(pseq.corePhylum) # 2 core phyla- Ascomycota and Basidiomycota
  
  abundCorephylum<- data.frame(otu_table(pseq.corePhylum)) #How abundant are they across samples?
  head(abundCorephylum)
  MeanAbundanceP<- 100*(apply(abundCorephylum[ ,1:259],1,mean)) #mean abundance 
  SDAbundancesP<- 100*(apply(abundCorephylum[ ,1:259],1,sd)) #standard deviation of abundance across samples
  prevalence(pseq.corePhylum, detection = 0.1/100, sort = TRUE) #prevelance of each core ASV- proportion samples they are in.
  
# Family level
  pseq.relFam<- microbiome::transform(ps.fam, "compositional")
  otu_table(pseq.relFam)
  sample_sums(pseq.relFam)
  # Pick the core (>0.1% relative abundance in >50% of the samples)
  pseq.coreFam <- core(pseq.relFam, detection = 0.1/100, prevalence = 50/100)
  tax_table(pseq.coreFam) # 13 core families (note 2 are unidentified at family level)
  taxonomyCoreFam<- data.frame(tax_table(pseq.coreFam))
  
  abundCoreFam<- data.frame(otu_table(pseq.coreFam)) #How abundant are they across samples?
  head(abundCoreFam)
  MeanAbundanceF<- data.frame(100*(apply(abundCoreFam[ ,1:259],1,mean))) #mean abundance in captive samples
  colnames(MeanAbundanceF)<- "MeanAbundance"
  SDAbundanceF<- data.frame(100*(apply(abundCoreFam[ ,1:259],1,sd))) #standard deviation of abundance across samples
  colnames(SDAbundanceF)<-"SDAbundance"
  PrevF<- data.frame(prevalence(pseq.coreFam, detection = 0.1/100, sort = TRUE)) #prevelance of each core ASV- proportion samples they are in.
  colnames(PrevF)<-"Prevalence"
  
  taxonomyCoreFam$MeanAbundance<- MeanAbundanceF$MeanAbundance
  taxonomyCoreFam$SDAbundance<- SDAbundanceF$SDAbundance
  taxonomyCoreFam$Prevalence<- PrevF$Prevalence
  
  #write.csv(taxonomyCoreFam, "CoreFamiliesFungal.csv")
  
# Genus level
  ps.Gen #detected 475
  pseq.relGen<- microbiome::transform(ps.Gen, "compositional")
  otu_table(pseq.relGen)
  sample_sums(pseq.relGen)
  # Pick the core (>0.1% relative abundance in >50% of the samples)
  pseq.coreGen <- core(pseq.relGen, detection = 0.1/100, prevalence = 50/100)
  tax_table(pseq.coreGen) # 13 core families (note 2 are unidentified at family level)
  taxonomyCoreGen<- data.frame(tax_table(pseq.coreGen))
  
  abundCoreGen<- data.frame(otu_table(pseq.coreGen)) #How abundant are they across samples?
  head(abundCoreGen)
  MeanAbundanceG<- data.frame(100*(apply(abundCoreGen[ ,1:259],1,mean))) #mean abundance in captive samples
  colnames(MeanAbundanceG)<- "MeanAbundance"
  SDAbundanceG<- data.frame(100*(apply(abundCoreGen[ ,1:259],1,sd))) #standard deviation of abundance across samples
  colnames(SDAbundanceG)<-"SDAbundance"
  PrevG<- data.frame(prevalence(pseq.coreGen, detection = 0.1/100, sort = TRUE)) #prevelance of each core ASV- proportion samples they are in.
  colnames(PrevG)<-"Prevalence"
  
  taxonomyCoreGen$MeanAbundance<- MeanAbundanceG$MeanAbundance
  taxonomyCoreGen$SDAbundance<- SDAbundanceG$SDAbundance
  taxonomyCoreGen$Prevalence<- PrevG$Prevalence
  
  #write.csv(taxonomyCoreGen, "CoreGenusFungal.csv")
  

  
#################################  
#### CALCULATE ALPHA DIVERSITY ##
#################################
  
# calculate diversity metrics using estimate_richness()
  richnessEstRare<-estimate_richness(physeqRare2, split=TRUE, measures= c("Chao1", "Shannon", "observed"))
  head(richnessEstRare)
  str(richnessEstRare)
  
# add alpha diversity metrics to metadata
  physeqRareMeta <- as.data.frame(sample_data(physeqRare2))
  head(physeqRareMeta)
  str(physeqRareMeta)
  physeqRareMeta$Chao1Fungi <- richnessEstRare$Chao1
  physeqRareMeta$ShannonFungi <- richnessEstRare$Shannon
  physeqRareMeta$ObservedFungi <- richnessEstRare$Observed
  
  #write.csv(physeqRareMeta, "FungalAlphaDiversity.csv")
  
  
################################################################  
######### 3. alpha diversity- bring in bacterial sequences #####
################################################################  

#### Generate Phyloseq object ####

  #ASV table
  asv_tableB <- read.csv ("FilteredASVTableBacteria.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_tableB) 
  # should be 18878 taxa and 586 samples (unrarefied reads)
  # Note that all ASVs unassigned to phylum,archaea etc... have been removed as well as potential contaminants (filtered using decontam)
  # control samples have been filtered
  # faecal samples with <10000 reads have been removed
  # ASVs with <50 reads across all samples have been removed
  # duplicate extractions still present
  asv_tableB <- as.matrix (asv_tableB) #make into a matrix
  
  #taxonomy
  taxonomyB <- read.csv ("FilteredTaxonomyTableBacteria.csv", row.names=1)
  str (taxonomyB) #should be 18878 observations, 7 taxonomic groupings, feature names as row names.
  taxonomyB <- as.matrix (taxonomyB)
  
  #load filtered metadata with sample names as rownames- 586 observations
  metadataB<-read.csv("FilteredMetadataBacteria.csv", row.names = 1)
  head(metadataB)  
  str(metadataB)
  
  #import all as phyloseq objects
  ASVbact <- otu_table(asv_tableB, taxa_are_rows = TRUE)
  TAXbact <- tax_table(taxonomyB)
  METAbact <- sample_data(metadataB)
  head(METAbact)
  
  #check that the ASV and sample names are consistent across objects (e.g have dashes been replaced with dots?)
  str(taxa_names(TAXbact))
  str(taxa_names(ASVbact))
  
  str(sample_names(ASVbact))
  str(sample_names(METAbact))
  
  #### MERGE INTO PHYLOSEQ OBJECT ####
  physeqBacteria <- phyloseq(ASVbact, TAXbact, METAbact)
  physeqBacteria 
  

##### RAREFY READS TO MIN SAMPLING DEPTH ######
  
  #rarefy to 10000 and set seed before rarefying, so that results are reproducible  
  physeqRareBacteria<-rarefy_even_depth(physeqBacteria, 10000, rngseed = 287)
  
  # 260 ASVs removed after subsampling- leaves 18852 across 586 samples
  physeqRareBacteria
  sample_sums(physeqRareBacteria)

  
##### CALCULATE ALPHA DIVERSITY METRICS ####
  
  # calculate observed, chao1, shannon using estimate_richness()
  richnessEstRareBacteria<-estimate_richness(physeqRareBacteria, split=TRUE, measures= c("Chao1", "Shannon", "Observed"))
  head(richnessEstRareBacteria)
  str(richnessEstRareBacteria)
  
  physeqBacteriaMeta <- as.data.frame(sample_data(physeqRareBacteria))
  physeqBacteriaMeta$Chao1 <- richnessEstRareBacteria$Chao1
  physeqBacteriaMeta$Shannon <- richnessEstRareBacteria$Shannon
  physeqBacteriaMeta$Observed <- richnessEstRareBacteria$Observed
  head(physeqBacteriaMeta)
  str(physeqBacteriaMeta)
  #write.csv(physeqBacteriaMeta, "BacterialAlphaDiversity.csv")
  
  
################################
## Merge fungal and bacterial ##
################################

#fungal alpha diversity data
  fungiAlpha<- read.csv("FungalAlphaDiversity.csv", row.names = 1)
  str(fungiAlpha) #259 samples
  length(unique(fungiAlpha$BirdID)) #259 birds
  
#bacterial alpha diversity data
  bacteriaAlpha<- read.csv("BacterialAlphaDiversity.csv")
  bacteriaAlpha2<- bacteriaAlpha[,c("Sample.ID", "BirdID", "TubeNumber","Chao1", "Shannon", "Observed")]  
  colnames(bacteriaAlpha2)<- c("BacterialSequencingSample.ID", "BirdIDbacteria", "TubeNumberBacteria","Chao1Bacteria", "ShannonBacteria", "ObservedBacteria")  

#merge based on bacterial sequencing id
  
  FungiAndBacteria<- merge(fungiAlpha, bacteriaAlpha2, by="BacterialSequencingSample.ID", all.x = T)
  str(FungiAndBacteria) #259 samples
  FungiAndBacteria$ObservedBacteria
  FungiAndBacteria$ShannonBacteria
  #note that 2 of the bacterial samples didn't reach threshold of 10'000 reads- hence NA
  
  NoASVsBacteria<- FungiAndBacteria[!is.na(FungiAndBacteria$ObservedBacteria),]
  NoASVsBacteria<-as.vector(NoASVsBacteria$ObservedBacteria)
  mean(NoASVsBacteria) #278.53
  sd(NoASVsBacteria) #167.77
  
  ShannBacteria<- FungiAndBacteria[!is.na(FungiAndBacteria$ShannonBacteria),]
  ShannBacteria<-as.vector(ShannBacteria$ShannonBacteria)
  mean(ShannBacteria) #3.86
  sd(ShannBacteria) #1.15
  
  NoASVsFungi<- as.vector(FungiAndBacteria$ObservedFungi)
  mean(NoASVsFungi) #49.49- for the 259 samples
  sd(NoASVsFungi) #26.02
  
  #write.csv(FungiAndBacteria, "FungiAndBacteria.csv")

  
# plot of fungi vs bacterial diversity & correlation analysis
  
  FungiAndBacteria<- read.csv("FungiAndBacteria.csv")
  head(FungiAndBacteria)
  
  cor.test(FungiAndBacteria$ObservedBacteria, FungiAndBacteria$ObservedFungi)
  cor.test(FungiAndBacteria$ShannonBacteria, FungiAndBacteria$ShannonFungi)
  
  ObservedCor<- ggscatter(FungiAndBacteria, x = "ObservedFungi", y = "ObservedBacteria", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Observed Fungal ASVs", ylab = "Observed Bacterial ASVs")+
                  font("xlab", size = 18)+
                  font("ylab", size = 18)+
                  font("xy.text", size = 14)
  ObservedCor
  
  ShannonCorr<- ggscatter(FungiAndBacteria, x = "ShannonFungi", y = "ShannonBacteria", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Fungal Shannon diversity", ylab = "Bacterial Shannon diversity")+
                  font("xlab", size = 18)+
                  font("ylab", size = 18)+
                  font("xy.text", size = 14)
  ShannonCorr
  
  
################################
#### Alpha diversity Models ####
################################
  
  AlphaMeta<- read.csv("FungiAndBacteria.csv")
  str(AlphaMeta)  #259 samples
  mean(AlphaMeta$ShannonFungi)
  sd(AlphaMeta$ShannonFungi)
  mean(AlphaMeta$ObservedFungi)
  sd(AlphaMeta$ObservedFungi)
  

  
#######################  
#### (1) All data #####
#######################
  
  AlphaMeta$SexEstimate<- as.factor(AlphaMeta$SexEstimate)
  AlphaMeta$FieldPeriodIDBinary<- as.factor(AlphaMeta$FieldPeriodIDBinary)
  AlphaMeta$TimeOfDayBinary<- as.factor(ifelse(AlphaMeta$MinutesSinceSunrise<360, "AM", "PM"))
  AlphaMeta$AgeYears<- (AlphaMeta$AgeDays)/365.25
  AlphaMeta$HS_2<- (AlphaMeta$Hs_obs^2) # heterozygosity squared for models
  AlphaMeta$AgeYears_2<- (AlphaMeta$AgeYears^2) # Age squared for models
  AlphaMeta$SampleYear<- as.factor(AlphaMeta$SampleYear)
  

# check for missing TQ values
  AlphaMeta2<- AlphaMeta[!is.na(AlphaMeta$TQcorrectedWithAvBA),] #4 floaters removed- now 255 samples
  str(AlphaMeta2)
  
# subset to variables to use in model
  AlphaMetaModel<-  AlphaMeta2[,c("BirdID","sample.id","TubeNoFungal","ShannonFungi","ObservedFungi","ShannonBacteria",
                                  "Chao1Fungi", "FieldPeriodIDBinary", "SexEstimate","AgeYears", "SampleYear",
                                  "TQcorrectedWithAvBA","MinutesSinceSunrise",
                                  "TimeOfDayBinary", "Hs_obs", "HS_2", "AgeYears_2", 
                                  "DateToFreeze")]
  
#Use rescale to standardise continuous predictors (centre and divide by 2 standard deviations
  # values have a mean of 0, sd of 0.5)
  
  AlphaMetaModel$AgeYears_scaled<-rescale(AlphaMetaModel$AgeYears, binary.inputs = "full")
  AlphaMetaModel$AgeYears_2_scaled<-rescale(AlphaMetaModel$AgeYears_2, binary.inputs = "full")
  AlphaMetaModel$HS_scaled<-rescale(AlphaMetaModel$Hs_obs, binary.inputs = "full") 
  AlphaMetaModel$HS_2_scaled<-rescale(AlphaMetaModel$HS_2, binary.inputs = "full") 
  AlphaMetaModel$TQ_scaled<- rescale(AlphaMetaModel$TQcorrectedWithAvBA, binary.inputs = "full") 
  AlphaMetaModel$TimeToFreeze_scaled<- rescale(AlphaMetaModel$DateToFreeze, binary.inputs = "full") 
  
# check for outliers in response
  
  dotchart(AlphaMetaModel$ShannonFungi) #look good
  dotchart(AlphaMetaModel$ObservedFungi) #look good
  
# histograms of response
  
  hist(AlphaMetaModel$ShannonFungi) #normal
  hist(AlphaMetaModel$ObservedFungi) #negative binomial model needed as overdispersed count data
  mean(AlphaMetaModel$ObservedFungi)
  var(AlphaMetaModel$ObservedFungi)
  
### All data Shannon diversity model ###
  
#check VIFs
  DivShanVif<- lm(ShannonFungi ~ SexEstimate + AgeYears_scaled + HS_scaled + 
                    FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled, data=AlphaMetaModel)
  
  vif(DivShanVif) # all look good- no correlations
  
# Shannon model 1
  DiversityShannon<- lm(ShannonFungi ~ SexEstimate + AgeYears_scaled + AgeYears_2_scaled + HS_scaled + 
                          HS_2_scaled + FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled+ 
                          AgeYears_scaled*SexEstimate,
                        data=AlphaMetaModel)
  
  summary(DiversityShannon)
  Anova(DiversityShannon)
  
  ggplot(AlphaMetaModel, aes(SampleYear, ShannonFungi)) + geom_boxplot()

# Drop the interaction term as not significant - model av more tricky with interaction terms
  DiversityShannon2<- lm(ShannonFungi ~ SexEstimate + AgeYears_scaled + AgeYears_2_scaled + HS_scaled + 
                          HS_2_scaled + FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                        data=AlphaMetaModel)
  
  summary(DiversityShannon2)
  Anova(DiversityShannon2)
  
  #Residuals
  simulationOutput <- simulateResiduals(fittedModel = DiversityShannon2, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form =AlphaMetaModel$AgeYears_scaled)
  testDispersion(simulationOutput)
  
  
# Full model and model averaging
  
  DiversityShannonFull<-  lm(ShannonFungi ~ SexEstimate + AgeYears_scaled + AgeYears_2_scaled + HS_scaled + 
                              HS_2_scaled + FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                            data=AlphaMetaModel, na.action = "na.fail")
  
#use MUMIN::DREDGE to test all models options and find best one and average 
  #(make sure squared terms only occur with lower order term and not without)
  msubset <- expression(dc(AgeYears_scaled, AgeYears_2_scaled) &
                          dc(HS_scaled, HS_2_scaled))
  ShannonDiv_dredge<-dredge(DiversityShannonFull, subset=msubset)
  
#average the model set with delta 5AIC cutoff
  TopModels5AIC<- get.models(ShannonDiv_dredge, subset = delta < 5)
  TopModels5AIC
  averaged_shannon5AIC<- model.avg(TopModels5AIC)
  summary(averaged_shannon5AIC)
  
#average the model set with 7AIC cutoff- consistent results
  TopModels7AIC<- get.models(ShannonDiv_dredge, subset = delta < 7)
  summary(TopModels7AIC)
  averaged_shannon7AIC<- model.avg(TopModels7AIC)
  Top7div_Shannon<- summary(averaged_shannon7AIC)
  CoeffDivShannon<- Top7div_Shannon$coefmat.subset
  #write.csv(CoeffDivShannon, "CoeffShannonAllDataModel.csv")
  
#drop the quadratic terms for age and for heterozygosity as not significant
  # enable interpretation of higher order effects
  
# Full model and model averaging without the squared terms
  
  DiversityShannonFull_2<-  lm(ShannonFungi ~ SexEstimate + AgeYears_scaled + HS_scaled + 
                              FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                             data=AlphaMetaModel, na.action = "na.fail")
  
  #use MUMIN::DREDGE to test all models options and find best one and average
  ShannonDiv_dredge_2<-dredge(DiversityShannonFull_2)
  
  #average the model set with 7AIC cutoff- consistent results
  TopModels7AIC_2<- get.models(ShannonDiv_dredge_2, subset = delta < 7)
  summary(TopModels7AIC_2)
  averaged_shannon7AIC_2<- model.avg(TopModels7AIC_2)
  Top7div_Shannon_2<- summary(averaged_shannon7AIC_2)
  CoeffDivShannon_2<- Top7div_Shannon_2$coefmat.subset
  #write.csv(CoeffDivShannon_2, "CoeffShannonAllDataModel_noQuadratics.csv")
  
  
  
#Plot of Field Periods: only significant predictor
  
  SampleSize<- data.frame(table(AlphaMetaModel$FieldPeriodIDBinary))
  head(SampleSize)
  colnames(SampleSize)<-c("FieldPeriodIDBinary", "Size")
  
  ggplot(AlphaMetaModel, aes(FieldPeriodIDBinary, ShannonFungi)) + geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitter(width = 0.1)) +theme_bw() + theme(panel.grid = element_blank())+
    theme(axis.text = element_text(size=18), axis.title = element_text(size = 20))+
    xlab("\nSeason") + ylab("Fungal Shannon diversity\n")+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    geom_text(data=SampleSize,aes(x=FieldPeriodIDBinary, y=4.3, label=Size), size=6)
  
 
  

### All data Observed ASVs model ###
  
  #check VIFs
  DivObsVif<- glm.nb(ObservedFungi ~ SexEstimate + AgeYears_scaled + HS_scaled + SampleYear + 
                    FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled, data=AlphaMetaModel)
  
  vif(DivObsVif) # all look good- no correlations
  
  # Observed ASVs model 1
  DiversityObserved<- glm.nb(ObservedFungi ~ SexEstimate + AgeYears_scaled + AgeYears_2_scaled + HS_scaled + SampleYear+
                          HS_2_scaled + FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled+
                          AgeYears_scaled*SexEstimate,
                        data=AlphaMetaModel)
  
  summary(DiversityObserved)
  Anova(DiversityObserved)
  
  
  # Drop the interaction term as not significant 
  DiversityObserved2<- lm(ObservedFungi ~ SexEstimate + AgeYears_scaled + AgeYears_2_scaled + HS_scaled + SampleYear + 
                           HS_2_scaled + FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                         data=AlphaMetaModel)
  
  summary(DiversityObserved2)
  Anova(DiversityObserved2)
  
  #Residuals
  simulationOutput <- simulateResiduals(fittedModel = DiversityObserved2, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form =AlphaMetaModel$AgeYears_scaled)
  testDispersion(simulationOutput)
  
  
  # Full model and model averaging
  
  DiversityObservedFull<-  glm.nb(ObservedFungi ~ SexEstimate + AgeYears_scaled + AgeYears_2_scaled + HS_scaled + 
                               HS_2_scaled + FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                             data=AlphaMetaModel, na.action = "na.fail")
  
  #use MUMIN::DREDGE to test all models options and find best one and average
  msubset <- expression(dc(AgeYears_scaled, AgeYears_2_scaled) &
                          dc(HS_scaled, HS_2_scaled))
  ObservedDiv_dredge<-dredge(DiversityObservedFull, subset=msubset)
  
  #average the model set with delta 5AIC cutoff
  TopModels5AIC<- get.models(ObservedDiv_dredge, subset = delta < 5)
  TopModels5AIC
  averaged_Observed5AIC<- model.avg(TopModels5AIC)
  summary(averaged_Observed5AIC)
  
  #average the model set with 7AIC cutoff- consistent results
  TopModels7AIC<- get.models(ObservedDiv_dredge, subset = delta < 7)
  summary(TopModels7AIC)
  averaged_Observed7AIC<- model.avg(TopModels7AIC)
  Top7div_Observed<- summary(averaged_Observed7AIC)
  CoeffDivObserved<- Top7div_Observed$coefmat.subset
  #write.csv(CoeffDivObserved, "CoeffObservedAllDataModel.csv")
  
# remove the squared terms for age and hs as not sig- allows interpretation of main effects
  
  
# Full model and model averaging without squared terms
  
  DiversityObservedFull_2<-  glm.nb(ObservedFungi ~ SexEstimate + AgeYears_scaled + + HS_scaled + 
                                    FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                                  data=AlphaMetaModel, na.action = "na.fail")
  
  #use MUMIN::DREDGE to test all models options and find best one and average
  ObservedDiv_dredge_2<-dredge(DiversityObservedFull_2)
  
  #average the model set with 7AIC cutoff- consistent results
  TopModels7AIC_2<- get.models(ObservedDiv_dredge_2, subset = delta < 7)
  summary(TopModels7AIC_2)
  averaged_Observed7AIC_2<- model.avg(TopModels7AIC_2)
  Top7div_Observed_2<- summary(averaged_Observed7AIC_2)
  CoeffDivObserved_2<- Top7div_Observed_2$coefmat.subset
  write.csv(CoeffDivObserved_2, "CoeffObservedAllDataModel_noSquaredTerms.csv")
  
  
  
  #Plot of Field Periods: significant predictor
  ggplot(AlphaMetaModel, aes(FieldPeriodIDBinary, ObservedFungi)) + geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitter(width = 0.1)) +theme_bw() + theme(panel.grid = element_blank())+
    theme(axis.text = element_text(size=18), axis.title = element_text(size = 20))+
    xlab("\nSeason") + ylab("Observed fungal ASV richness\n")+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") +
    geom_text(data=SampleSize,aes(x=FieldPeriodIDBinary, y=158, label=Size), size=6)
  
  table(AlphaMetaModel$FieldPeriodIDBinary)

  #plot of time to freezing
  ggplot(AlphaMetaModel, aes(DateToFreeze, ObservedFungi)) + geom_point()
  
  
  
##################################    
#### (2) smaller MHC dataset #####
##################################    
  
  MHCAlphameta<- read.csv("FungiAndBacteria.csv")
  table(MHCAlphameta$Both_MHC_present)  #192 birds have complete MHC data (17 don't have both and 50 don't have any data)
  
# subset to complete cases for MHC data
  
  MHCDataAvailable<- MHCAlphameta[!is.na(MHCAlphameta$Both_MHC_present),] #209
  MHCDataBothPresent <- MHCDataAvailable[MHCDataAvailable$Both_MHC_present == "Y",]
  str(MHCDataBothPresent) #192 samples/birds
  
#check for missing TQ values
  MHCDataBothPresent[is.na(MHCDataBothPresent$TQcorrectedWithAvBA),] #4 floaters
  
  
#subset to variables used in model
  
  MHCAlphaModelsSubset<- MHCDataBothPresent[,c("BirdID","sample.id","TubeNoFungal","ShannonFungi","ObservedFungi",
                                         "Chao1Fungi", "MHC1_Diversity", "MHC2_Diversity",
                                         "MHC1_Diversity_2", "MHC2_Diversity_2", "Ase.ua5", "Ase.ua3","Ase.ua6",
                                         "Ase.ua7","Ase.ua9","Ase.ua11","Ase.ua8","Ase.ua4","Ase.ua1","Ase.ua10",
                                         "Ase.dab1", "Ase.dab2",  "Ase.dab3", "Ase.dab4",  "Ase.dab5",
                                         "TLR3","FieldPeriodIDBinary", "SexEstimate","AgeDays","Ageclass",
                                         "TQcorrectedWithAvBA","MinutesSinceSunrise","Hs_obs",
                                         "DateToFreeze")]
  
   
  #write.csv(MHCAlphaModelsSubset, "MHCAlphaModels.csv", row.names = F)
  

  
##### MHC diversity ######
  
  MHCAlphaModels<- read.csv("MHCAlphaModels.csv")
  str(MHCAlphaModels) #192 samples/birds
  
#remove samples with no territory quality values- from three floaters with no assigned territory
  
  MHCAlphaModels<-MHCAlphaModels[!is.na(MHCAlphaModels$TQcorrectedWithAvBA),]
  str(MHCAlphaModels) #189 samples remain
  
  
#check for correlations between variables in MHC diversity model

  cor_mhcDIV<- MHCAlphaModels[,c("MHC1_Diversity", "MHC2_Diversity","Hs_obs", "TQcorrectedWithAvBA")]
  MHC_Div_CORREL<-ggcorr(cor_mhcDIV,label=TRUE)
  MHC_Div_CORREL # no strong correlation between continuous variables
  

#format variables correctly  
  MHCAlphaModels$BirdID<- as.factor(MHCAlphaModels$BirdID)
  MHCAlphaModels$sample.id<- as.factor(MHCAlphaModels$sample.id)
  MHCAlphaModels$TubeNoFungal<- as.factor(MHCAlphaModels$TubeNoFungal)
  MHCAlphaModels$SexEstimate<- as.factor(MHCAlphaModels$SexEstimate)
  MHCAlphaModels$FieldPeriodIDBinary<- as.factor(MHCAlphaModels$FieldPeriodIDBinary)
  MHCAlphaModels$Ageclass<- as.factor(MHCAlphaModels$Ageclass)
  MHCAlphaModels$TLR3<- as.factor(MHCAlphaModels$TLR3)
  MHCAlphaModels$Ase.ua10<- as.factor(MHCAlphaModels$Ase.ua10) #P=39,A=150 -20.63%
  MHCAlphaModels$Ase.ua1<- as.factor(MHCAlphaModels$Ase.ua1) #P=39,A=150 - 20.63%
  MHCAlphaModels$Ase.ua4<- as.factor(MHCAlphaModels$Ase.ua4) #P=35,A=154 - 18.52%
  MHCAlphaModels$Ase.ua8<- as.factor(MHCAlphaModels$Ase.ua8) #P=38,A=151 - 20.12 %
  MHCAlphaModels$Ase.ua11<- as.factor(MHCAlphaModels$Ase.ua11) #P=108,A=81 -57.14 %
  MHCAlphaModels$Ase.ua9<- as.factor(MHCAlphaModels$Ase.ua9) #P=125,A=64 - 66.14%
  MHCAlphaModels$Ase.ua7<- as.factor(MHCAlphaModels$Ase.ua7) #P=144,A=45 - 76.19%
  MHCAlphaModels$Ase.ua6<- as.factor(MHCAlphaModels$Ase.ua6) #P=116,A=73 - 61.38%
  MHCAlphaModels$Ase.ua3<- as.factor(MHCAlphaModels$Ase.ua3) #P=122,A=67 - 64.55%
  MHCAlphaModels$Ase.ua5<- as.factor(MHCAlphaModels$Ase.ua5) #P=132,A=57 - 69.84%
  MHCAlphaModels$Ase.dab1<- as.factor(MHCAlphaModels$Ase.dab1) # P=187,A=2 - exclude from analysis as present in almost all- 98.94%
  MHCAlphaModels$Ase.dab2<- as.factor(MHCAlphaModels$Ase.dab2) # P=184, A= 5 - exclude from analysis as present in almost all- 97.35%
  MHCAlphaModels$Ase.dab3<- as.factor(MHCAlphaModels$Ase.dab3) # P=51, A=138 - 26.98%
  MHCAlphaModels$Ase.dab4<- as.factor(MHCAlphaModels$Ase.dab4) # P=47, A=142 - 24.87%
  MHCAlphaModels$Ase.dab5<- as.factor(MHCAlphaModels$Ase.dab5) # P=41, A=148 - 21.69%
  MHCAlphaModels$AgeYears<- (MHCAlphaModels$AgeDays)/365.25
  MHCAlphaModels$TimeOfDayBinary<- as.factor(ifelse(MHCAlphaModels$MinutesSinceSunrise<360, "AM", "PM"))
  str(MHCAlphaModels)
  
### Plot percentage of birds that each allele is present in
  AlleleNames<- as.factor(colnames(MHCAlphaModels[,11:25]))
  PercentAlleles<- c(69.84, 64.55, 61.38, 76.19, 66.14, 57.14, 20.12, 18.52, 20.63, 20.63, 98.94, 97.35, 26.98, 24.87, 21.69)
  ColourAllele<- c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue",
                   "cornflowerblue","cornflowerblue","goldenrod", "goldenrod","cornflowerblue","cornflowerblue","cornflowerblue")
  
  AllelePrevalence<- data.frame(AlleleNames, PercentAlleles, ColourAllele)
  str(AllelePrevalence)
  par(mar=c(5, 12, 2, 2))
  barplot(AllelePrevalence$PercentAlleles,
          xlab = "Allele prevalence",
          names.arg = AlleleNames,
          col = AllelePrevalence$ColourAllele,
          horiz = TRUE, 
          las=1, 
          xlim=c(0,100),
          ylim=c(0,18))
  


### 1) Shannon Diversity and MHC diversity ###
  
#Use rescale to standardise continuous predictors (centre and divide by 2 standard deviations- values have a mean of 0, sd of 0.5)
  
  MHCAlphaModels$MHC1_Div_scaled<-rescale(MHCAlphaModels$MHC1_Diversity, binary.inputs = "full") #MHC1 diversity
  MHCAlphaModels$MHC2_Div_scaled<-rescale(MHCAlphaModels$MHC2_Diversity, binary.inputs = "full") #MHC2 diversity
  MHCAlphaModels$MHC1_Div2_scaled<-rescale(MHCAlphaModels$MHC1_Diversity_2, binary.inputs = "full") #MHC1 diversity squared
  MHCAlphaModels$MHC2_Div2_scaled<-rescale(MHCAlphaModels$MHC2_Diversity_2, binary.inputs = "full") #MHC2 diversity squared
  MHCAlphaModels$HS_scaled<-rescale(MHCAlphaModels$Hs_obs, binary.inputs = "full") #heterozygosity
  MHCAlphaModels$TQ_scaled<-rescale(MHCAlphaModels$TQcorrectedWithAvBA, binary.inputs = "full") #Territory quality
  MHCAlphaModels$Age_scaled<-rescale(MHCAlphaModels$AgeYears, binary.inputs = "full") # Age
  MHCAlphaModels$TimeToFreeze_scaled<-rescale(MHCAlphaModels$DateToFreeze, binary.inputs = "full") # storage at 4 degree

#VIFs
  VIFCheck<- lm(ShannonFungi ~ MHC1_Div_scaled + MHC2_Div_scaled + HS_scaled + SexEstimate + TLR3+
                  Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled, data=MHCAlphaModels)
  summary(VIFCheck)  
  Anova(VIFCheck)
  
  vif(VIFCheck) #all below 2
  
#Residuals
  simulationOutput <- simulateResiduals(fittedModel = VIFCheck, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form = MHCAlphaModels$MHC1_Diversity)
  testDispersion(simulationOutput)
  

# 1) Full model (includes squared terms for MHC1 and MHC2 diversity)
  
 ShannonMHCDiversity<- lm(ShannonFungi ~ MHC1_Div_scaled + MHC2_Div_scaled + MHC1_Div2_scaled + MHC2_Div2_scaled +
                            HS_scaled + SexEstimate + TLR3+
                            Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                          data=MHCAlphaModels, na.action = "na.fail")
  
 summary(ShannonMHCDiversity)  
 Anova(ShannonMHCDiversity)
 
#use MUMIN::DREDGE to test all models options and find best one and average
 msubsetMHC <- expression(dc(MHC1_Div_scaled, MHC1_Div2_scaled) &
                         dc(MHC2_Div_scaled, MHC2_Div2_scaled))
 ShannonMHCDiv_dredge<-dredge(ShannonMHCDiversity, subset=msubsetMHC)
 
#average the model set with delta 5AIC cutoff
 TopModels5AIC<- get.models(ShannonMHCDiv_dredge, subset = delta < 5)
 TopModels5AIC
 averaged_shannon5AIC<- model.avg(TopModels5AIC)
 summary(averaged_shannon5AIC)
 
#average the model set with 7AIC cutoff- consistent results
 TopModels7AIC<- get.models(ShannonMHCDiv_dredge, subset = delta < 7)
 summary(TopModels7AIC)
 averaged_shannon7AIC<- model.avg(TopModels7AIC)
 Top7MHCdiv_Shannon<- summary(averaged_shannon7AIC)
 CoeffMHCDivShannon<- Top7MHCdiv_Shannon$coefmat.subset
 
 #write.csv(CoeffMHCDivShannon, "CoeffMHCDivShannonModel.csv")

# remove the squared terms for MHC diversity as not sig.
 
 
# 2) Full model (remove squared terms for MHC1 and MHC2 diversity)
 
 ShannonMHCDiversity_2<- lm(ShannonFungi ~ MHC1_Div_scaled + MHC2_Div_scaled +
                            HS_scaled + SexEstimate + TLR3+
                            Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                          data=MHCAlphaModels, na.action = "na.fail")
 
 summary(ShannonMHCDiversity_2)  
 Anova(ShannonMHCDiversity_2)
 
 #use MUMIN::DREDGE to test all models options and find best one and average
 ShannonMHCDiv_dredge_2<-dredge(ShannonMHCDiversity_2)

 #average the model set with 7AIC cutoff- consistent results
 TopModels7AIC_2<- get.models(ShannonMHCDiv_dredge_2, subset = delta < 7)
 summary(TopModels7AIC_2)
 averaged_shannon7AIC_2<- model.avg(TopModels7AIC_2)
 Top7MHCdiv_Shannon_2<- summary(averaged_shannon7AIC_2)
 CoeffMHCDivShannon_2<- Top7MHCdiv_Shannon_2$coefmat.subset
 
 #write.csv(CoeffMHCDivShannon_2, "CoeffMHCDivShannonModel_noSquaredTerms.csv")
 

  
# Plot of MHC diversity and fungal shannon diversity - Not significant
  MHCAlphaModels$MHC1_DiversityFactor<- as.factor(MHCAlphaModels$MHC1_Diversity)
  MHCAlphaModels$MHC2_DiversityFactor<- as.factor(MHCAlphaModels$MHC2_Diversity)
  
  ggplot(MHCAlphaModels, aes(MHC1_Diversity, ShannonFungi)) + geom_point() 
  ggplot(MHCAlphaModels, aes(MHC2_Diversity, ShannonFungi)) + geom_point() 
 
  ggplot(MHCAlphaModels, aes(MHC1_DiversityFactor, ShannonFungi)) + geom_violin() + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") + theme_bw() + theme(panel.grid = element_blank())+
    theme(axis.text = element_text(size=14), axis.title = element_text(size = 14))+
    xlab("MHC-I diversity") + ylab("Fungal Shannon diversity")
  
  ggplot(MHCAlphaModels, aes(MHC2_DiversityFactor, ShannonFungi)) + geom_violin() + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") + theme_bw() + theme(panel.grid = element_blank())+
    theme(axis.text = element_text(size=14), axis.title = element_text(size = 14))+
    xlab("MHC-II diversity") + ylab("Fungal Shannon diversity")
  
  
### 2) Observed ASVs and MHC diversity ### 
  
# Observed ASVs are count data- probably overdispersed
  mean(MHCAlphaModels$ObservedFungi)
  var(MHCAlphaModels$ObservedFungi) # variance is much greater than the mean- overdispersed- need negative binomial model
  
  
#VIFs
  VIFCheckObs1<- glm.nb(ObservedFungi ~ MHC1_Div_scaled + MHC2_Div_scaled + HS_scaled + SexEstimate + TLR3+
                  Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary +TimeToFreeze_scaled, data=MHCAlphaModels)
  
  summary(VIFCheckObs1)
  Anova(VIFCheckObs1)
  
  vif(VIFCheckObs1) #all below 2
  
#Residuals
  simulationOutput <- simulateResiduals(fittedModel = VIFCheckObs1, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form = MHCAlphaModels$MHC1_Div_scaled)
  testDispersion(simulationOutput)
  
  
# 1) Full model (includes squared values for MHC1 and MHC2 diversity)
  
  ObservedMHCDiversity<- glm.nb(ObservedFungi ~ MHC1_Div_scaled + MHC2_Div_scaled + MHC1_Div2_scaled + MHC2_Div2_scaled +
                                  HS_scaled + SexEstimate + TLR3+
                                  Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                                data=MHCAlphaModels, na.action = "na.fail")
  
  summary(ObservedMHCDiversity)  
  Anova(ObservedMHCDiversity)
  
  #use MUMIN::DREDGE to test all models options and find best one and average
  msubsetMHC <- expression(dc(MHC1_Div_scaled, MHC1_Div2_scaled) &
                             dc(MHC2_Div_scaled, MHC2_Div2_scaled))
  ObservedMHCDiv_dredge<-dredge(ObservedMHCDiversity)
  
  #average the model set with 2AIC cutoff
  TopModels5AIC_ASVsObserved<- get.models(ObservedMHCDiv_dredge, subset = delta < 5)
  TopModels5AIC_ASVsObserved 
  averaged_Observed5AIC<- model.avg(TopModels5AIC_ASVsObserved)
  summary(averaged_Observed5AIC)
  
  #average the model set with 7AIC cutoff
  TopModels7AIC_ASVsObserved<- get.models(ObservedMHCDiv_dredge, subset = delta < 7)
  TopModels7AIC_ASVsObserved 
  averaged_Observed7AIC<- model.avg(TopModels7AIC_ASVsObserved)
  summary(averaged_Observed7AIC)
  Top7MHCdiv_Observed<- summary(averaged_Observed7AIC)
  CoeffMHCDivObserved<- Top7MHCdiv_Observed$coefmat.subset
  
  #write.csv(CoeffMHCDivObserved, "CoeffMHCDivObservedModel.csv")
  

# 2) Full model (remove squared values for MHC1 and MHC2 diversity)
  
  ObservedMHCDiversity_2<- glm.nb(ObservedFungi ~ MHC1_Div_scaled + MHC2_Div_scaled +
                                  HS_scaled + SexEstimate + TLR3+
                                  Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled,
                                data=MHCAlphaModels, na.action = "na.fail")
  
  summary(ObservedMHCDiversity_2)  
  Anova(ObservedMHCDiversity_2)
  
  #use MUMIN::DREDGE to test all models options and find best one and average
  ObservedMHCDiv_dredge_2<-dredge(ObservedMHCDiversity_2)
  
  #average the model set with 7AIC cutoff
  TopModels7AIC_ASVsObserved_2<- get.models(ObservedMHCDiv_dredge_2, subset = delta < 7)
  TopModels7AIC_ASVsObserved_2 
  averaged_Observed7AIC_2<- model.avg(TopModels7AIC_ASVsObserved_2)
  summary(averaged_Observed7AIC_2)
  Top7MHCdiv_Observed_2<- summary(averaged_Observed7AIC_2)
  CoeffMHCDivObserved_2<- Top7MHCdiv_Observed_2$coefmat.subset
  
  #write.csv(CoeffMHCDivObserved_2, "CoeffMHCDivObservedModel_noSquaredTerms.csv")
  
  
  
#Plot to check field period- not significant
  ggplot(MHCAlphaModels, aes(FieldPeriodIDBinary, ObservedFungi)) + geom_boxplot()+
    geom_point(position=position_jitter(width = 0.1)) +theme_bw() + theme(panel.grid = element_blank())+
    theme(axis.text = element_text(size=14), axis.title = element_text(size = 14))+
    xlab("Field Period") + ylab("Fungal Shannon diversity")+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") 
  table(MHCAlphaModels$FieldPeriodIDBinary)
  
#Plots of MHC diversity- not significant
  ggplot(MHCAlphaModels, aes(MHC1_DiversityFactor, ObservedFungi)) + geom_violin() + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red")
  
  ggplot(MHCAlphaModels, aes(MHC2_DiversityFactor, ObservedFungi)) + geom_violin() + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red")

  
################################
#### ALLELE PRESENCE/ABSENCE ###
################################
  
 
### check for correlations between all alleles (presence/absence)
 
 cor_mhc<- MHCAlphaModelsSubset[,c("Ase.ua5", "Ase.ua3","Ase.ua6",
                             "Ase.ua7","Ase.ua9","Ase.ua11","Ase.ua8","Ase.ua4","Ase.ua1","Ase.ua10",
                             "Ase.dab3", "Ase.dab4",  "Ase.dab5")]
  
  
 MHC_allele_CORREL<-ggcorr(cor_mhc,label=TRUE)
 MHC_allele_CORREL #some level of correlation present- need to check VIFs and P/A in more detail
 
 
 MHCAlphaModels1and10<- as.factor(MHCAlphaModels$Ase.ua1==MHCAlphaModels$Ase.ua10)
 table(MHCAlphaModels1and10) # Ase1 and Ase10 are the same (present or absent) 98% of the time
 #only include one of these in the model
 

 
#### CheckingVIFs #####
 
 ## all alleles in model 
 VIFCheckAllelesShannonAll<- glm(ShannonFungi ~ Ase.ua5 + Ase.ua3 + Ase.ua6 + Ase.ua7 + Ase.ua9 + 
                Ase.ua11 + Ase.ua8 + Ase.ua4 + Ase.ua1 + Ase.ua10 + Ase.dab3 + Ase.dab4 + TLR3+
                  Ase.dab5 + HS_scaled + SexEstimate +
                  Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary, data=MHCAlphaModels)

  summary(VIFCheckAllelesShannonAll)  
  vif(VIFCheckAllelesShannonAll) # Vifs very high (~9) for Ase.ua1 and 10- remove ua10 as the highest
  
  ## remove ua10
  VIFCheckAllelesShannon2<- glm(ShannonFungi ~ Ase.ua1 + Ase.ua5 + Ase.ua3 + Ase.ua6 + Ase.ua7 + Ase.ua9 + 
                                   Ase.ua11 + Ase.ua8 + Ase.ua4 + Ase.dab3 + Ase.dab4 +
                                   Ase.dab5 + HS_scaled + SexEstimate + TLR3 +
                                   Age_scaled +FieldPeriodIDBinary +TQ_scaled + TimeOfDayBinary, data=MHCAlphaModels)
  
  vif(VIFCheckAllelesShannon2) # Vifs are now ~5 or less for all other terms - remove this from the model
  summary(VIFCheckAllelesShannon2)
  

### check number of genotypes lost by removing variables ###
  
# MHC1 genotypes- all alleles present
  MHCAlphaModels$MHC1Genotype<- paste(MHCAlphaModels$Ase.ua1,MHCAlphaModels$Ase.ua3, MHCAlphaModels$Ase.ua4, 
                                      MHCAlphaModels$Ase.ua5, MHCAlphaModels$Ase.ua6, MHCAlphaModels$Ase.ua7, 
                                      MHCAlphaModels$Ase.ua8, MHCAlphaModels$Ase.ua9, MHCAlphaModels$Ase.ua11)
  length(unique(MHCAlphaModels$MHC1Genotype)) # 41 different combinations (as opposed to 44 whenAse10 also present)
  ggplot(MHCAlphaModels, aes(MHC1Genotype)) + geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #most common = 010111011 - ua3, 5, 6, 7, 9, 11 are all present
  #next most common = 010001011 - ua3, 7, 9, 11 are all present
  #next most common = 000110000 - ua5 and 6 present
  
  
  
#### 1) Shannon diversity and P/A of MHC alleles
  
#### Full model (no Ase.ua10) ###
  ShannonMHCAlleles<- lm(ShannonFungi ~ Ase.ua1 + Ase.ua3 + Ase.ua4 + Ase.ua5 +  Ase.ua6 + Ase.ua7 +
                           Ase.ua8 + Ase.ua9 + Ase.ua11+ Ase.dab3 + Ase.dab4 + Ase.dab5 +
                           HS_scaled + SexEstimate + TLR3 + Age_scaled +FieldPeriodIDBinary +
                           TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled, data=MHCAlphaModels, na.action="na.fail")
  
  
  summary(ShannonMHCAlleles)  
  Anova(ShannonMHCAlleles)
  
#Residuals
  simulationOutput <- simulateResiduals(fittedModel = ShannonMHCAlleles, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form = MHCAlphaModels$Ase.ua6)
  testDispersion(simulationOutput)
  
#use MUMIN::DREDGE to test all models options and find best one and average
  ShannonMHCAlleles_dredge<-dredge(ShannonMHCAlleles)
  
#average the model set with 5AIC cutoff
  TopModels5AIC_allelesS<- get.models(ShannonMHCAlleles_dredge, subset = delta < 5)
  TopModels5AIC_allelesS 
  averaged_shannon5AIC_Alleles<- model.avg(TopModels5AIC_allelesS)
  summary(averaged_shannon5AIC_Alleles)
  
  #average the model set with 7AIC cutoff
  TopModels7AIC_allelesS<- get.models(ShannonMHCAlleles_dredge, subset = delta < 7)
  TopModels7AIC_allelesS 
  averaged_shannon7AIC_Alleles<- model.avg(TopModels7AIC_allelesS)
  Top7ASVAlleles_Shannon<- summary(averaged_shannon7AIC_Alleles)
  CoeffMHCAllelesShannon<- Top7ASVAlleles_Shannon$coefmat.subset
  
  #write.csv(CoeffMHCAllelesShannon, "CoeffMHC_Alleles_ShannonModel.csv")

  
  ggplot(MHCAlphaModels, aes(Ase.ua11, ShannonFungi)) + geom_violin() + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red")
  
  ggplot(MHCAlphaModels, aes(Ase.ua11, ShannonFungi)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=4, fill="red") + 
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18))
  
  
#plot model outcomes
  
  ModelShannResults<- read.csv("CoeffMHC_Alleles_ShannonModel.csv") # file with variables in correct order
  ModelShannResults
  ModelShannResults<-ModelShannResults[c(-1,-8),] #remove intercept
  ModelShannResults
  str(ModelShannResults)
  
  ModelShannResults$SignificancePoints<- c(rep("Not significant", times=2), "Significant", rep("Not significant", times=13),"Significant", "Not significant", "Not significant", "Not significant")
  
  ModelShannResults$X<-factor(ModelShannResults$X, c("Ase.dab5","Ase.dab4", "Ase.dab3", "Ase.ua11", "Ase.ua9", "Ase.ua8", "Ase.ua7",
                                                     "Ase.ua6", "Ase.ua5","Ase.ua4", "Ase.ua3", "Ase.ua1","AC","CC","TimeFreeze",
                                                     "TimeDay","TerritoryQuality","HS", "Season","Sex","Age"))
  
  ShannModOut<- ggplot(ModelShannResults, aes(X,Estimate, fill=SignificancePoints, col=SignificancePoints)) + theme_bw()+
    geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE, col=SignificancePoints), width=0.2, size=1) +
    geom_point(size=4, shape=21) + 
    scale_fill_manual(values=c("grey70","goldenrod")) + 
    scale_colour_manual(values=c("grey70","goldenrod")) + 
    coord_flip() + xlab("") + ylab("\nEstimate") + geom_hline(aes(yintercept=0), linetype="dashed") +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  ShannModOut
  
  ShannModOutNoLables<- ggplot(ModelShannResults, aes(X,Estimate, fill=SignificancePoints, col=SignificancePoints)) + theme_bw()+
    geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE, col=SignificancePoints), width=0.2, size=1) +
    geom_point(size=4, shape=21) +
    scale_fill_manual(values=c("grey70","goldenrod")) + 
    scale_colour_manual(values=c("grey70","goldenrod")) + 
    coord_flip() + xlab("") + ylab("") + geom_hline(aes(yintercept=0), linetype="dashed") +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  ShannModOutNoLables
  
### 2) Observed ASVs and Allele presence/absence ### 

  
### Full model of allele p/a
  
  Alleles_observedASV<- glm.nb(ObservedFungi ~ Ase.ua1+ Ase.ua3 + Ase.ua4 + Ase.ua5 +  Ase.ua6 + Ase.ua7 +
                                 Ase.ua8 + Ase.ua9 + Ase.ua11+ Ase.dab3 + Ase.dab4 + Ase.dab5 +
                                 HS_scaled + SexEstimate + TLR3 + Age_scaled +FieldPeriodIDBinary +
                                 TQ_scaled + TimeOfDayBinary + TimeToFreeze_scaled, data=MHCAlphaModels, na.action="na.fail")
  
  
  vif(Alleles_observedASV) # all ~5 or less
  summary(Alleles_observedASV)  
  Anova(Alleles_observedASV)
  
  #use MUMIN::DREDGE to test all models options and find best one and average
  ObservedMHCAlleles_dredge<-dredge(Alleles_observedASV, trace=2)
  
  #average the model set with 5 AIC cutoff
  TopModels5AIC_ASVsObserved<- get.models(ObservedMHCAlleles_dredge, subset = delta < 5)
  TopModels5AIC_ASVsObserved #results in a set of 12 models
  averaged_Observed5AIC<- model.avg(TopModels5AIC_ASVsObserved)
  summary(averaged_Observed5AIC)
  
  #average the model set with 7AIC cutoff
  TopModels7AIC_ASVsObserved<- get.models(ObservedMHCAlleles_dredge, subset = delta < 7)
  TopModels7AIC_ASVsObserved 
  averaged_Observed7AIC<- model.avg(TopModels7AIC_ASVsObserved)
  Top7ASVAlleles_Observed<- summary(averaged_Observed7AIC)
  CoeffMHCAllelesObserved<- Top7ASVAlleles_Observed$coefmat.subset
  
  write.csv(CoeffMHCAllelesObserved, "CoeffMHC_Alleles_ObservedModel.csv")
  
  
 Asedab3<- ggplot(MHCAlphaModels, aes(Ase.dab3, ObservedFungi)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="red") + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.y=element_blank(), axis.text = element_text(size=16), axis.title = element_text(size=18)) +
    labs(x ="Ase-dab3")+
    theme(axis.title.x = element_text(face = "italic"))+
    scale_x_discrete(labels=c("0" = "absent", "1" = "present"))+
   theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))

  
  Aseua4<- ggplot(MHCAlphaModels, aes(Ase.ua4, ObservedFungi)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="red") + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18)) + ylab("Observed ASVs\n")+
    labs(x ="Ase-ua4")+
    theme(axis.title.x = element_text(face = "italic"))+
    scale_x_discrete(labels=c("0" = "absent", "1" = "present"))+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))
  
  Aseua7<- ggplot(MHCAlphaModels, aes(Ase.ua7, ObservedFungi)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="red") + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.y= element_blank(), axis.text = element_text(size=16), axis.title.x = element_text(size=18)) +
    labs(x ="Ase-ua7")+
    theme(axis.title.x = element_text(face = "italic"))+
    scale_x_discrete(labels=c("0" = "absent", "1" = "present"))+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))
  
  Aseua8<- ggplot(MHCAlphaModels, aes(Ase.ua8, ObservedFungi)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width=0.1))+
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="red") + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18)) + ylab("Observed ASVs")+
    labs(x ="Ase-ua8")+
    theme(axis.title.x = element_text(face = "italic"))+
    scale_x_discrete(labels=c("0" = "absent", "1" = "present"))+
    theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  

  plot((Aseua4|Aseua7)/(Aseua8|Asedab3))


#calculate means

  Ase.dab3mean<- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.dab3, mean) 
  Ase.dab3SD<- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.dab3, sd) 

  Ase.ua4mean<- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.ua4, mean) 
  Ase.ua4SD <- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.ua4, sd) 

  Ase.ua7mean<- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.ua7, mean) 
  Ase.ua7SD <- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.ua7, sd) 

  Ase.ua8mean<- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.ua8, mean) 
  Ase.ua8SD <- tapply(MHCAlphaModels$ObservedFungi,MHCAlphaModels$Ase.ua8, sd) 

  
#plot model outcomes
  
  ModelOResults<- read.csv("CoeffMHC_Alleles_ObservedModel.csv")
  ModelOResults
  ModelOResults<-ModelOResults[c(-1),]
  ModelOResults
  str(ModelOResults)
  
  ModelOResults$Significance<- c(rep("Not significant",times=11),"Significant","Not significant","Not significant", "Significant","Significant",
                                 rep("Not significant", times=2),"Significant","Not significant", "Not significant" )
  
  ModelOResults$X<-factor(ModelOResults$X, c("Ase.dab5","Ase.dab4", "Ase.dab3", "Ase.ua11", "Ase.ua9", "Ase.ua8", "Ase.ua7",
                                                     "Ase.ua6", "Ase.ua5","Ase.ua4", "Ase.ua3", "Ase.ua1","TLR3AC","TLR3CC","TimeStored",
                                                     "TimeDay","TQ", "Season","HS","Sex","Age"))
  
  
  ObsModOut<- ggplot(ModelOResults, aes(X,Estimate, fill=Significance, col=Significance)) + theme_bw()+
    geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE, col=Significance), width=0.2, size=1) +
    geom_point(size=4, shape=21) +
    scale_fill_manual(values=c("grey70","goldenrod")) + 
    scale_colour_manual(values=c("grey70","goldenrod")) + 
    coord_flip() + xlab("") + ylab("") + geom_hline(aes(yintercept=0), linetype="dashed") +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  ObsModOut
  
  ObsModOutNoLabels<- ggplot(ModelOResults, aes(X,Estimate, fill=Significance, col=Significance)) + theme_bw()+
    geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE, col=Significance), width=0.2, size=1) +
    geom_point(size=4, shape=21) +
    scale_fill_manual(values=c("grey70","goldenrod")) + 
    scale_colour_manual(values=c("grey70","goldenrod")) + 
    coord_flip() + xlab("") + ylab("\nEstimate") + geom_hline(aes(yintercept=0), linetype="dashed") +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ObsModOutNoLabels
  
  
  
  
  
################################
####### Survival Models ########
################################
  
  AlphaSurvival<- read.csv("FungiAndBacteria.csv")
  str(AlphaSurvival)
  AlphaSurvival$SurvivedNextSeasonFactor<- as.factor(AlphaSurvival$SurvivedNextSeasonFactor)
  AlphaSurvival$SampleYear<- as.factor(AlphaSurvival$SampleYear)
  AlphaSurvival$ShannonFungi_2<- AlphaSurvival$ShannonFungi^2
  AlphaSurvival$ObservedFungi_2<- AlphaSurvival$ObservedFungi^2
  AlphaSurvival$SexEstimate<- as.factor(AlphaSurvival$SexEstimate)
  AlphaSurvival$FieldPeriodIDBinary<- as.factor(AlphaSurvival$FieldPeriodIDBinary)
  AlphaSurvival$TimeOfDayBinary<- as.factor(ifelse(AlphaSurvival$MinutesSinceSunrise<360, "AM", "PM"))
  AlphaSurvival$AgeYears<- (AlphaSurvival$AgeDays)/365.25
  AlphaSurvival$Observer<- as.factor(AlphaSurvival$Observer)
  AlphaSurvival$AgeYears_2<- (AlphaSurvival$AgeYears^2) # Age squared for models
  
# check for missing TQ values
  AlphaSurvival2<- AlphaSurvival[!is.na(AlphaSurvival$TQcorrectedWithAvBA),] #4 floaters removed- now 255 samples
  str(AlphaSurvival2) #255 samples
  
# How many survived/died?
  table(AlphaSurvival2$SurvivedNextSeasonFactor)  #37 died, 218 survived by the next breeding season
  
# rescale predictors
  AlphaSurvival2$AgeYears_scaled<-rescale(AlphaSurvival2$AgeYears, binary.inputs = "full")
  AlphaSurvival2$AgeYears_2_scaled<-rescale(AlphaSurvival2$AgeYears_2, binary.inputs = "full")
  AlphaSurvival2$TQ_scaled<- rescale(AlphaSurvival2$TQcorrectedWithAvBA, binary.inputs = "full") 
  AlphaSurvival2$TimeOfDay_scaled <- rescale(AlphaSurvival2$MinutesSinceSunrise, binary.inputs = "full") 
  AlphaSurvival2$ShannonFungi_scaled<- rescale(AlphaSurvival2$ShannonFungi)
  AlphaSurvival2$ShannonFungi_2_scaled<- rescale(AlphaSurvival2$ShannonFungi_2)
  AlphaSurvival2$ObservedFungi_scaled<- rescale(AlphaSurvival2$ObservedFungi)
  AlphaSurvival2$ObservedFungi_2_scaled<- rescale(AlphaSurvival2$ObservedFungi_2)
  AlphaSurvival2$TerritoryID<- as.factor(AlphaSurvival2$TerritoryID)
  
  #write.csv(AlphaSurvival2,"survivalData.csv")
  
  
  
### 1) Shannon model ###
  
# Model vifs
  
  SurvivalShannonVif<- glm(SurvivedNextSeason~ ShannonFungi_scaled + AgeYears_scaled + SexEstimate + TQ_scaled+ FieldPeriodIDBinary +SampleYear, data= AlphaSurvival2, family=binomial(link=logit))
  vif(SurvivalShannonVif)
  
  
# Full model
  
  SurvivalShannonMod1<- glm(SurvivedNextSeason~ ShannonFungi_scaled + ShannonFungi_2_scaled+ AgeYears_scaled + AgeYears_2_scaled+
                              SexEstimate + TQ_scaled+ FieldPeriodIDBinary + SampleYear+ ShannonFungi_scaled*AgeYears_scaled,
                            data= AlphaSurvival2, family=binomial(link=logit), na.action = "na.fail")
  summary(SurvivalShannonMod1)
  
  
# remove interaction with age
  SurvivalShannonMod2<- glm(SurvivedNextSeason~ ShannonFungi_scaled + ShannonFungi_2_scaled+ AgeYears_scaled + AgeYears_2_scaled+
                              SexEstimate + TQ_scaled+ FieldPeriodIDBinary + SampleYear,
                            data= AlphaSurvival2, family=binomial(link=logit), na.action = "na.fail")
  summary(SurvivalShannonMod2)
  
  
#Residuals
  simulationOutput <- simulateResiduals(fittedModel =  SurvivalShannonMod2, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form = AlphaSurvival2$ShannonFungi_scaled)
  testDispersion(simulationOutput)
  
  
#use MUMIN::DREDGE to test all models options and find best one and average
  msubset <- expression(dc(AgeYears_scaled, AgeYears_2_scaled) &
                          dc(ShannonFungi_scaled, ShannonFungi_2_scaled))
  ShannonDivSurv_dredge<-dredge(SurvivalShannonMod2, subset=msubset)
  
  
#average the model set with 7AIC cutoff
  TopModels7AIC<- get.models(ShannonDivSurv_dredge, subset = delta < 7)
  summary(TopModels7AIC)
  averaged_shannon7AIC<- model.avg(TopModels7AIC)
  Top7div_Shannon<- summary(averaged_shannon7AIC)
  Top7div_Shannon
  CoeffDivShannon<- Top7div_Shannon$coefmat.subset
  #write.csv(CoeffDivShannon, "CoeffShannonSurvival.csv")
  
# remove quadratic terms for age and shannon diversity as not significant
  
# Full model no squared terms
  
  SurvivalShannonMod1_noSq<- glm(SurvivedNextSeason~ ShannonFungi_scaled + AgeYears_scaled +
                                   SexEstimate + TQ_scaled+ FieldPeriodIDBinary + SampleYear,
                                 data= AlphaSurvival2, family=binomial(link=logit), na.action = "na.fail")
  summary(SurvivalShannonMod1_noSq)
  
#use MUMIN::DREDGE to test all models options and find best one and average
  ShannonDivSurv_dredge_2<-dredge(SurvivalShannonMod1_noSq)
  
  
#average the model set with 7AIC cutoff
  TopModels7AIC_2<- get.models(ShannonDivSurv_dredge_2, subset = delta < 7)
  summary(TopModels7AIC_2)
  averaged_shannon7AIC_2<- model.avg(TopModels7AIC_2)
  Top7div_Shannon_2<- summary(averaged_shannon7AIC_2)
  Top7div_Shannon_2
  CoeffDivShannon_2<- Top7div_Shannon_2$coefmat.subset
  #write.csv(CoeffDivShannon_2, "CoeffShannonSurvivalNosq.csv")
  
  
  
### 2) Observed ASVs model ###
  
# Model vifs
  SurvivalObservedVif<- glm(SurvivedNextSeason~ ObservedFungi_scaled + AgeYears_scaled + SexEstimate + TQ_scaled+
                              FieldPeriodIDBinary + SampleYear, data= AlphaSurvival2, family=binomial(link=logit))
  vif(SurvivalObservedVif)
  
  
# Full model
  SurvivalObservedMod1<- glm(SurvivedNextSeason~ ObservedFungi_scaled + ObservedFungi_2_scaled+ AgeYears_scaled + AgeYears_2_scaled+
                               SexEstimate + TQ_scaled+ FieldPeriodIDBinary+ SampleYear + ObservedFungi_scaled*AgeYears_scaled,
                             data= AlphaSurvival2, family=binomial(link=logit), na.action = "na.fail")
  summary(SurvivalObservedMod1)
  
# remove interaction with age
  SurvivalObservedMod2<- glm(SurvivedNextSeason~ ObservedFungi_scaled + ObservedFungi_2_scaled + AgeYears_scaled + AgeYears_2_scaled+
                               SexEstimate + TQ_scaled+ FieldPeriodIDBinary+ SampleYear,
                             data= AlphaSurvival2, family=binomial(link=logit), na.action = "na.fail")
  summary(SurvivalShannonMod2)
  
  
#Residuals
  simulationOutput <- simulateResiduals(fittedModel =  SurvivalObservedMod2, plot = F)
  plot(simulationOutput) # all look good- no significant deviations or over dispersion detected
  plotResiduals(simulationOutput, form = AlphaSurvival2$ObservedFungi_scaled)
  testDispersion(simulationOutput)
  
  
#use MUMIN::DREDGE to test all models options and find best one and average
  msubsetMHC <- expression(dc(MHC1_Div_scaled, MHC1_Div2_scaled) &
                             dc(MHC2_Div_scaled, MHC2_Div2_scaled))
  ObservedDivSurv_dredge<-dredge(SurvivalObservedMod2, subset=msubset)
  
  
#average the model set with 7AIC cutoff
  TopModels7AICObs<- get.models(ObservedDivSurv_dredge, subset = delta < 7)
  summary(TopModels7AICObs)
  averaged_Obs7AIC<- model.avg(TopModels7AICObs)
  Top7div_Observed<- summary(averaged_Obs7AIC)
  Top7div_Observed
  CoeffDivObserved<- Top7div_Observed$coefmat.subset
  write.csv(CoeffDivObserved, "CoeffObservedSurvival.csv")
  
  
# Full model - remove the square terms as not sig
  
  SurvivalObservedMod1_2<- glm(SurvivedNextSeason~ ObservedFungi_scaled +  AgeYears_scaled +
                                 SexEstimate + TQ_scaled+ FieldPeriodIDBinary+ SampleYear,
                               data= AlphaSurvival2, family=binomial(link=logit), na.action = "na.fail")
  summary(SurvivalObservedMod1_2)
  
#use MUMIN::DREDGE to test all models options and find best one and average
  ObservedDivSurv_dredge_2<-dredge(SurvivalObservedMod1_2)
  
  
#average the model set with 7AIC cutoff
  TopModels7AICObs_2<- get.models(ObservedDivSurv_dredge_2, subset = delta < 7)
  summary(TopModels7AICObs_2)
  averaged_Obs7AIC_2<- model.avg(TopModels7AICObs_2)
  Top7div_Observed_2<- summary(averaged_Obs7AIC_2)
  Top7div_Observed_2
  CoeffDivObserved_2<- Top7div_Observed_2$coefmat.subset
  #write.csv(CoeffDivObserved_2, "CoeffObservedSurvivalNosq.csv")
  
  
  
  