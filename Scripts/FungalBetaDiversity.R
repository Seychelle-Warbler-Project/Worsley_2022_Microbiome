#### FungalBetaDiversity.R ####
#This R script contains the code necessary to reproduce the results of the analysis investigating the 
#association between environment/host variables and gut microbiome beta diversity. 
#It uses the filtered ASV, taxonomy and metadata files generated in FungalSequences.R


#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(plyr)
  library(RColorBrewer)
  library(forcats)
  library(microbiome)
  library(vegan)
  library(BiodiversityR)
  library(ggordiplots)
  library(data.table)
  library(dplyr)
  library(nloptr)
  library(ANCOMBC)
  library(corncob)
  library(stringr)

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

  
### REMOVE EXTRACTION DUPLICATES ###
# retain the one with the highest read count but exclude both duplicates of 408 as bird sampled again in 2020 (tube 490)
  
  TubesToRemove<- c("FUN370R", "FUN387", "FUN408", "FUN408R", "FUN424", "FUN453")
  
  physeq2<- prune_samples(!sample_data(physeq)$TubeNoFungal %in% TubesToRemove, physeq)
  physeq2 # 259 samples remain after duplicates removed  

  
### PRUNE RARE TAXA (LOW PREVALENCE) ###
  
# Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(physeq2),
                 MARGIN = ifelse(taxa_are_rows(physeq2), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq2),
                      tax_table(physeq2))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
  # Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq2, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq2),color=Phylum)) +
    # Include a guess for the threshold- here 1% = approx 3 samples
    geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 1% samples
  prevalenceThreshold<-1*(259/100)
  
# Execute the prevalence filter, using `prune_taxa()` function
  head(prevdf1)
  
  KeepTaxa1prev<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1prev)
  ASVNamesToKeep<- rownames(KeepTaxa1prev)
  
  physeq2
  physeqfiltered<- prune_taxa(ASVNamesToKeep, physeq2)
  physeqfiltered # 794 taxa, 259 samples  
  
####################################
### CLR TRANSFORM ASV ABUNDANCES ###
####################################  

  physeq_clr <- microbiome::transform(physeqfiltered, "clr")
  head(otu_table(physeq_clr))
  
  CLRSamples<- data.frame(sample_data(physeq_clr))
  str(CLRSamples) #259
  
  
### analyse all data #### 
# Remove the four samples with no TQ value (these are floaters)
  
  floatersRemoved<- CLRSamples[!is.na(CLRSamples$TQcorrectedWithAvBA),]
  str(floatersRemoved) #255 samples/birds
  
  SampleIDsToKeep<- floatersRemoved$sample.id
  
  physeq_clr2<- subset_samples(physeq_clr, sample.id %in% SampleIDsToKeep) #255 samples, 794 taxa
  
  
  #function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  
#Extract ASV Matrix and Sample Data  
  clr_MatrixAll<-vegan_otu(physeq_clr2)
  clr_SampleDataAll<-as(sample_data(physeq_clr2),"data.frame")
  str(clr_SampleDataAll) #255 samples
  

# subset to variables in model and correct format of sample data
  
  All_SampleData<- clr_SampleDataAll[,c("BirdID","sample.id","TubeNoFungal", "SampleYear",
                                     "FieldPeriodIDBinary", "SexEstimate","AgeDays","Ageclass",
                                     "TQcorrectedWithAvBA","MinutesSinceSunrise","Hs_obs", "SurvivedNextSeasonFactor", 
                                     "SurvivedNextSeason", "DateToFreeze", "FieldPeriodID")]
  
  str(All_SampleData)
  
  All_SampleData$SexEstimate <- as.factor(All_SampleData$SexEstimate)
  All_SampleData$FieldPeriodID <- as.factor(All_SampleData$FieldPeriodID)
  All_SampleData$SampleYear <- as.factor(All_SampleData$SampleYear)
  All_SampleData$FieldPeriodIDBinary <- as.factor(All_SampleData$FieldPeriodIDBinary)
  All_SampleData$TimeOfDayBinary <- as.factor(ifelse(All_SampleData$MinutesSinceSunrise < 360, "AM", "PM"))
  All_SampleData$AgeYears<- All_SampleData$AgeDays/365.25
  All_SampleData$SurvivedNextSeasonFactor<- as.factor(All_SampleData$SurvivedNextSeasonFactor)
  All_SampleData$Ageclass<- as.factor(All_SampleData$Ageclass)
  
  
ggplot(All_SampleData, aes(SampleYear, TQcorrectedWithAvBA)) + geom_boxplot()
ggplot(All_SampleData, aes(FieldPeriodID, TQcorrectedWithAvBA)) + geom_boxplot()
  
### PERMANOVA with larger dataset ###
  
  perm <- how(nperm = 9999)
  set.seed(3779)
  permanovaAll<- adonis2(clr_MatrixAll ~ AgeYears + SexEstimate + Hs_obs+ TQcorrectedWithAvBA + 
                                    FieldPeriodIDBinary + TimeOfDayBinary + SurvivedNextSeasonFactor + DateToFreeze,
                                  data=All_SampleData, permutations = perm, method = "euclidean", by= "margin")
  permanovaAll
  
  ResultsPermanovaAll<- data.frame(permanovaAll)
  #write.csv(ResultsPermanovaAll, "PermanovaLargerDataset.csv")

  
### betadisper tests to check homogeneity of variance ###
  
  distMatrixBetaAll <- vegdist(clr_MatrixAll, method="euclidean")
  
# field period betadisper 
  BDFieldPeriodAll<-betadisper(distMatrixBetaAll,All_SampleData$FieldPeriodIDBinary)
  set.seed(2699)
  permutest(BDFieldPeriodAll, permutations = 9999) 
  par(mar=c(6,6,4,4))  
  boxplot(BDFieldPeriodAll,ylab="Distance to centroid\n", xlab="\nFieldPeriod", cex.axis=1.5, cex.lab=1.5) #minor higher variance
  

# TQ betadisper  
  summary(All_SampleData$TQcorrectedWithAvBA) 
  hist(All_SampleData$TQcorrectedWithAvBA)
 
  All_SampleData$TQFactor<- as.factor(ifelse(All_SampleData$TQcorrectedWithAvBA < 12678, "low","high")) #above or below the median
                                          
  plot( All_SampleData$TQFactor)
  
  BDTQall<-betadisper(distMatrixBetaAll,All_SampleData$TQFactor)
  set.seed(269789)
  permutest(BDTQall, permutations = 9999) #no difference in variance between categories for TQ
  

# Time of day betadisper
  BDTimeAll<-betadisper(distMatrixBetaAll,All_SampleData$TimeOfDayBinary)
  set.seed(2959)
  permutest(BDTimeAll, permutations = 9999) #no significant difference in variability (0.0541)
  boxplot(BDTimeAll,ylab="Distance to centroid\n", xlab="\nTime of day", cex.axis=1.5, cex.lab=1.5)
  
# Survival betadisper
  BDSurvivalAll<-betadisper(distMatrixBetaAll,All_SampleData$SurvivedNextSeasonFactor)
  set.seed(2959)
  permutest(BDSurvivalAll, permutations = 9999) #no significant difference in variability (0.0541)
  boxplot(BDSurvivalAll,ylab="Distance to centroid\n", xlab="\nSurvival", cex.axis=1.5, cex.lab=1.5)
  
  
 # Freezing betadisper
  summary(All_SampleData$DateToFreeze)
  All_SampleData$DateToFreezeFactor<- ifelse(All_SampleData$DateToFreeze<39, "low", "high")
  BDDateToFreeze<-betadisper(distMatrixBetaAll,All_SampleData$DateToFreezeFactor)
  set.seed(2959)
  permutest(BDDateToFreeze, permutations = 9999) #no significant difference in variability



##############
### envfit ###
##############

  clr_pcaAll<-rda(clr_MatrixAll, scale=FALSE)
  
#uses linear model permutations to map variables onto an ordination#
  
  envfitDat<- All_SampleData[,c("AgeYears", "SexEstimate", "Hs_obs", "TQcorrectedWithAvBA",
                                "FieldPeriodIDBinary", "TimeOfDayBinary", "SurvivedNextSeasonFactor", "DateToFreeze")]
  set.seed(67689)
  envfitmod<- envfit(clr_pcaAll, envfitDat, permutations = 9999, choices=c(1,2))
  envfitmod #relationship with axes 1,2 - datetofreeze, field period and time of day (greatest R2) associated with shifts
  
  
  set.seed(67999)
  envfitmod3<- envfit(clr_pcaAll, envfitDat, permutations = 9999, choices=c(3,4))
  envfitmod3 # TQ associated with shifts along 3 and 4 (primarily 4), date to freeze, survival, field period also associated (primarily shifted along axis 3)
  
  
#plotting envfit arrows onto PCA for TQ#
  #Principal components 3 and 4 - territory quality is primarily associated with shift along PC4
  
  #plotting envfit arrows onto PCA#
  arrow_factor <- ordiArrowMul(envfitmod3, fill=0.5) #get proportions, longer arrows= stronger predictor
  spp.scrs <- as.data.frame(scores(envfitmod3, display = "vectors"))* arrow_factor
  spp.scrs$Species <- c("AgeYears", "Hs_obs", "TQcorrectedWithAvBA", "DateToFreeze")
  
  scorePCA<- data.frame(scores(clr_pcaAll, choices=c(3,4), display = "sites", scaling=1))
  summary(scorePCA)
  colnames(scorePCA)<- c("PC3", "PC4")
  scorePCA$FieldPeriod<- as.factor(All_SampleData$FieldPeriodIDBinary)
  scorePCA$TQFactor<- as.factor(All_SampleData$TQFactor)
  scorePCA$Survival<- as.factor(All_SampleData$SurvivedNextSeasonFactor)
  
  
  EnvPlot34<- ggplot(data = scorePCA, aes(x = PC3, y = PC4)) + theme_bw() + geom_point(size = 2, stroke=2, aes(col=FieldPeriod)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC3 2.41%") + ylab("PC4 2.09%\n") +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18),
          legend.text = element_text(size=18))+
    coord_fixed()+
    geom_segment(data = spp.scrs,
                 aes(x = 0, xend = PC3, y = 0, yend = PC4),
                 arrow = arrow(length = unit(0.25, "cm")), colour = "grey31", size=1.5)+
    geom_text(data = spp.scrs, aes(label= Species, x = PC3*1.3, y = PC4*1.3),
              size = 5)
  
  EnvPlot34
  
  
  
#####################################  
### Principal Components Analysis ###
#####################################
  
  
  clr_pcaAll<-rda(clr_MatrixAll)
  
  screeplot(clr_pcaAll) # PC1 gives most information, followed by 2.
  
  sig <- PCAsignificance (clr_pcaAll, axes = 8) 
  sig # PC1 = 4.314, PC2 = 3.542, PC3 = 2.481, PC4 = 2.094

  
## PCA of field period ##

  plotPCAFPeriodAll<-gg_ordiplot(clr_pcaAll,groups=All_SampleData$FieldPeriodIDBinary,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCAFPeriodAll
  
  plotPCAFPeriodAll2<-gg_ordiplot(clr_pcaAll,groups=All_SampleData$FieldPeriodIDBinary,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(3,4)) 
  plotPCAFPeriodAll2
  
  ord.data <- plotPCAFPeriodAll2$df_ord
  str(ord.data)
  ord.data$FieldPeriod<- All_SampleData$FieldPeriodIDBinary
  ord.data$Sample.ID<- All_SampleData$sample.id
  head(ord.data)
  
  #centroids
  pca_scores<-scores(clr_pcaAll, choices=c(1,2,3,4), scaling=1)
  sample_scores<- data.frame(pca_scores$sites)
  str(sample_scores)
  summary(sample_scores)
  sample_scores<-data.frame(setDT(sample_scores, keep.rownames = TRUE)[])
  str(sample_scores)
  colnames(sample_scores)<-c("Sample.ID", "PC1", "PC2","PC3", "PC4")
  sample_scores$FieldPeriod<- All_SampleData$FieldPeriodIDBinary
  
  centroidsPC1<-data.frame(ddply(sample_scores, .(FieldPeriod), summarize, mean=mean(PC1)))
  centroidsPC2<-data.frame(ddply(sample_scores, .(FieldPeriod), summarize, mean=mean(PC2)))
  centroidsPC3<-data.frame(ddply(sample_scores, .(FieldPeriod), summarize, mean=mean(PC3)))
  centroidsPC4<-data.frame(ddply(sample_scores, .(FieldPeriod), summarize, mean=mean(PC4)))
  centroids<- merge(centroidsPC3, centroidsPC4, by="FieldPeriod")
  colnames(centroids)<- c("FieldPeriod", "x", "y")
  
  PlotFPeriod<- ggplot(data = ord.data, aes(x = x, y = y, fill = FieldPeriod)) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC3 2.90%") + ylab("PC4 2.48%\n") + 
    scale_fill_manual(values=c("goldenrod2", "purple3")) +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE)) +
    geom_point(size = 4, stroke=1, shape=21, alpha=0.8, colour="grey55")+ geom_point(data=centroids, aes(x,y), shape=21, stroke=4, size=8, colour="black")
    
  
  PlotFPeriod
  
  
  
## PCA of TQ ##
  
  plot(All_SampleData$TQFactor)
  
  #PC1/PC2
  plotPCATQ<-gg_ordiplot(clr_pcaAll,groups=All_SampleData$TQFactor,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCATQ
  
  ord.dataTQ <- plotPCATQ$df_ord
  str(ord.dataTQ)
  ord.dataTQ$TQFactor<- All_SampleData$TQFactor
  ord.dataTQ$TQ<- All_SampleData$TQcorrectedWithAvBA
  head(ord.dataTQ)
  
  ggplot(ord.dataTQ, aes(TQ, x)) + geom_point() +ylab("PC1")
  ggplot(ord.dataTQ, aes(TQ, y)) + geom_point() +ylab("PC2")
  
  cor.test(ord.dataTQ$TQ, ord.dataTQ$x)
  cor.test(ord.dataTQ$TQ, ord.dataTQ$x)  
  
  #PC3/PC4
  plotPCATQ2<-gg_ordiplot(clr_pcaAll,groups=All_SampleData$TQFactor,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(3,4)) 
  plotPCATQ2
  
  ord.dataTQ2 <- plotPCATQ2$df_ord
  str(ord.dataTQ2)
  ord.dataTQ2$TQFactor<- All_SampleData$TQFactor
  ord.dataTQ2$TQ<- All_SampleData$TQcorrectedWithAvBA
  head(ord.dataTQ2)
  
  ggplot(ord.dataTQ2, aes(TQ, x)) + geom_point() +ylab("PC3")
  ggplot(ord.dataTQ2, aes(TQ, y)) + geom_point() +ylab("PC4")
  
  cor.test(ord.dataTQ2$TQ, ord.dataTQ2$x)
  cor.test(ord.dataTQ2$TQ, ord.dataTQ2$y) #correlated with PC4
  
  #centroids for PC3/4
  pca_scoresTQ<-scores(clr_pcaAll, choices=c(1,2,3,4), scaling=1)
  sample_scoresTQ<- data.frame(pca_scoresTQ$sites)
  str(sample_scoresTQ)
  summary(sample_scoresTQ)
  sample_scoresTQ<-data.frame(setDT(sample_scoresTQ, keep.rownames = TRUE)[])
  str(sample_scoresTQ)
  colnames(sample_scoresTQ)<-c("Sample.ID", "PC1", "PC2","PC3", "PC4")
  sample_scoresTQ$TQFactor<- All_SampleData$TQFactor
  
  centroidsPC3TQ<-data.frame(ddply(sample_scoresTQ, .(TQFactor), summarize, mean=mean(PC3)))
  centroidsPC4TQ<-data.frame(ddply(sample_scoresTQ, .(TQFactor), summarize, mean=mean(PC4)))
  centroidsTQ<- merge(centroidsPC3TQ, centroidsPC4TQ, by="TQFactor")
  colnames(centroidsTQ)<- c("TQFactor", "x", "y")
  
  PlotTQ<- ggplot(data = ord.dataTQ2, aes(x = x, y = y, fill = TQFactor)) + 
    theme_bw() + geom_point(size = 4, stroke=1, shape=21, alpha=0.8, colour="grey55")+
    geom_point(data=centroidsTQ, aes(x,y), shape=21, stroke=4, size=8, alpha=0.9, colour="black")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC3 2.90%") + ylab("PC4 2.48%\n") + 
    scale_fill_manual(values=c("purple3", "goldenrod2")) +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotTQ
  
  
## PCA of Survival ##
  
  plotPCASurvivalall<-gg_ordiplot(clr_pcaAll,groups=All_SampleData$SurvivedNextSeasonFactor,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(3,4)) 
  plotPCASurvivalall  
  
  ord.dataSurvival <- plotPCASurvivalall$df_ord
  str(ord.dataSurvival)
  ord.dataSurvival$Survival<- All_SampleData$SurvivedNextSeasonFactor
  ord.dataSurvival$TQFactor<- All_SampleData$TQFactor
  head(ord.dataSurvival)
  
  ggplot(ord.dataSurvival, aes(Survival, x)) + geom_boxplot() +ylab("PC3")
  ggplot(ord.dataSurvival, aes(Survival, y)) + geom_boxplot() +ylab("PC4")
  
  #centroids for PC3,4
  pca_scoresSurvival<-scores(clr_pcaAll, choices=c(1,2,3,4), scaling=1)
  sample_scoresSurvival<- data.frame(pca_scoresSurvival$sites)
  str(sample_scoresSurvival)
  summary(sample_scoresSurvival)
  sample_scoresSurvival<-data.frame(setDT(sample_scoresSurvival, keep.rownames = TRUE)[])
  str(sample_scoresSurvival)
  colnames(sample_scoresSurvival)<-c("Sample.ID", "PC1", "PC2","PC3", "PC4")
  sample_scoresSurvival$Survival<- All_SampleData$SurvivedNextSeasonFactor
  
  centroidsPC3Survival<-data.frame(ddply(sample_scoresSurvival, .(Survival), summarize, mean=mean(PC3)))
  centroidsPC4Survival<-data.frame(ddply(sample_scoresSurvival, .(Survival), summarize, mean=mean(PC4)))
  centroidsSurvival<- merge(centroidsPC3Survival, centroidsPC4Survival, by="Survival")
  colnames(centroidsSurvival)<- c("Survival", "x", "y")
  
  PlotSurvival<- ggplot(data = ord.dataSurvival, aes(x = x, y = y, fill = Survival)) + 
    theme_bw() + geom_point(size = 4, stroke=1, alpha=0.8, colour="grey55", shape=21)+
    geom_point(data=centroidsSurvival, aes(x,y), shape=21, stroke=4, size=8, alpha=0.9, colour="black")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC3 2.90%") + ylab("PC4 2.48%\n") + 
    scale_fill_manual(values=c("violetred4", "lightsteelblue1"), labels=c("died", "survived")) +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), 
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotSurvival
  
## PCA of Time of day ##
  
  plotPCATimeall<-gg_ordiplot(clr_pcaAll,groups=All_SampleData$TimeOfDayBinary,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCATimeall  
  
  ord.dataTime <- plotPCATimeall$df_ord
  str(ord.dataTime)
  ord.dataTime$Time<- All_SampleData$TimeOfDayBinary
  head(ord.dataTime)
  
  ggplot(ord.dataTime, aes(Time, x)) + geom_boxplot() +ylab("PC1")
  ggplot(ord.dataTime, aes(Time, y)) + geom_boxplot() +ylab("PC2")
  
  #centroids for PC1,2
  pca_scoresTime<-scores(clr_pcaAll, choices=c(1,2,3,4), scaling=1)
  sample_scoresTime<- data.frame(pca_scoresTime$sites)
  str(sample_scoresTime)
  summary(sample_scoresTime)
  sample_scoresTime<-data.frame(setDT(sample_scoresTime, keep.rownames = TRUE)[])
  str(sample_scoresTime)
  colnames(sample_scoresTime)<-c("Sample.ID", "PC1", "PC2","PC3", "PC4")
  sample_scoresTime$Time<- All_SampleData$TimeOfDayBinary
  
  centroidsPC1Time<-data.frame(ddply(sample_scoresTime, .(Time), summarize, mean=mean(PC1)))
  centroidsPC2Time<-data.frame(ddply(sample_scoresTime, .(Time), summarize, mean=mean(PC2)))
  centroidsTime<- merge(centroidsPC1Time, centroidsPC2Time, by="Time")
  colnames(centroidsTime)<- c("Time", "x", "y")
  
  PlotTime<- ggplot(data = ord.dataTime, aes(x = x, y = y, fill = Time)) + 
    theme_bw() + geom_point(size = 4, stroke=1, shape=21, alpha=0.8, colour="grey55")+
    geom_point(data=centroidsTime, aes(x,y), shape=21, stroke=4, size=8, alpha=0.9, colour="black")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 4.31%") + ylab("PC2 3.54%\n") + 
    scale_fill_manual(values=c("goldenrod2", "purple3")) +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), 
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotTime
  
  
  
  
### Differential abundance testing ###
#### Ancom BC https://www.nature.com/articles/s41467-020-17041-7
  
  #install latest release of AncomBC
  
  #install.packages("remotes")
  #remotes::install_github("FrederickHuangLin/ANCOMBC")
  
  library(ANCOMBC)
  
  physeqfiltered # unrarefied reads, but filtered to remove rare taxa, not CLR transformed
  phyfilteredDat<- data.frame(sample_data(physeqfiltered))
  
  floatersRemoved<- phyfilteredDat[!is.na(phyfilteredDat$TQcorrectedWithAvBA),]
  str(floatersRemoved) #255 samples/birds
  
  SampleIDsToKeep<- floatersRemoved$sample.id
  
  physeqfiltered2<- subset_samples(physeqfiltered, sample.id %in% SampleIDsToKeep) #255 samples, 794 taxa
  
  str(sample_data(physeqfiltered2))
  sample_data(physeqfiltered2)$SexEstimate <- as.factor(sample_data(physeqfiltered2)$SexEstimate)
  sample_data(physeqfiltered2)$FieldPeriodIDBinary <- as.factor(sample_data(physeqfiltered2)$FieldPeriodIDBinary)
  sample_data(physeqfiltered2)$SurvivedNextSeasonFactor <- as.factor(sample_data(physeqfiltered2)$SurvivedNextSeasonFactor)
  sample_data(physeqfiltered2)$TimeOfDayBinary<- as.factor(ifelse(sample_data(physeqfiltered2)$MinutesSinceSunrise< 360, "AM", "PM"))
  sample_data(physeqfiltered2)$AgeYears<- sample_data(physeqfiltered2)$AgeDays/365.25
  sample_data(physeqfiltered2)$TQFactor<- as.factor(ifelse(sample_data(physeqfiltered2)$TQcorrectedWithAvBA < 12678, "low", "high"))
   
  

# running Ancom BC
  
  AncomDA<-ancombc(physeqfiltered2, formula= "AgeYears + TQFactor + FieldPeriodIDBinary + SexEstimate + Hs_obs +
                   TimeOfDayBinary + SurvivedNextSeasonFactor+ DateToFreeze", p_adj_method = "BH",
                   lib_cut = 0, alpha = 0.05, neg_lb = TRUE, group = "FieldPeriodIDBinary", struc_zero = TRUE,global = TRUE)
  str(AncomDA)
  
  results<-AncomDA$res
  str(results)
  head(results)
  
  diffAbund<-as.data.frame(results$diff_abn)
  head(diffAbund)
  str(diffAbund)
  
### FieldPeriod - 31 DA taxa
  DiffMinvMaj<-diffAbund[diffAbund$FieldPeriodIDBinaryMinor=="TRUE",]
  str(DiffMinvMaj)
  rowsFP<-row.names(DiffMinvMaj)
  
  results2<-as.data.frame(results)
  head(results2)
  significantMinVMaj<-subset(results2, rownames(results2) %in% rowsFP)
  head(significantMinVMaj)
  str(significantMinVMaj) # 31 observations 
  
  library(corncob)
  taxonMinVMaj<-otu_to_taxonomy(OTU = row.names(DiffMinvMaj), data = physeqfiltered2)
  taxonFP<-as.data.frame(taxonMinVMaj)
  significantMinVMaj$taxonomy<-taxonFP$taxon
  #View(significantMinVMaj)
  
  #write.csv(significantMinVMaj, "MajorVMinorAncom.csv")
  #separate out the taxonomy fields in the excel
  
 
## make a plot for major/minor
  
  sigSeason<- read.csv("MajorVMinorAncom.csv")
  str(sigSeason)
  head(sigSeason)
  
  significantSeason<-data.frame(sigSeason$beta.FieldPeriodIDBinaryMinor,sigSeason$Phylum, sigSeason$Class, sigSeason$Order, sigSeason$Family, sigSeason$se.FieldPeriodIDBinaryMinor)
  str(significantSeason)
  colnames(significantSeason)<-c("Beta",  "Phylum","Class", "Order", "Family", "se")
  significantSeason$Order<-as.factor(significantSeason$Order)
  significantSeason$Class<-as.factor(significantSeason$Class)
  significantSeason$Phylum<-as.factor(significantSeason$Phylum)
  significantSeason$Family<-as.factor(significantSeason$Family)
  significantSeason$Species<-as.factor(seq(from=1, to=31, by=1))
  
  #colour palette
  brewer.pal(n = 12, name = "Paired")
  
  ggplot(significantSeason, aes(x= Beta, y=fct_inorder(Family), colour=str_wrap(Order,15), group=Species)) + 
    geom_point(position=position_dodge(width = 0.8), size=2) +
    geom_errorbar(aes(xmin=Beta -se, xmax=Beta + se), width=0,lwd=1, position= position_dodge(width=0.8)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(name="\nLog fold change", breaks=c(-3,-2, -1,0, 1, 2)) + 
    ylab("Family\n")  + geom_vline(xintercept = 0, linetype="dotted", size=1.5) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=22),
          legend.text = element_text(size=20), legend.title = element_text(size=22)) +
    scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                                "#FFFF99", "#B15928","grey10", "grey55" ))+
    labs(colour="Order")
  
  
  
### TQ
    
  DiffTQHvL<-diffAbund[diffAbund$TQFactor=="TRUE",]
  str(DiffTQHvL)
  rowsTQ<-row.names(DiffTQHvL)
  
  results2<-as.data.frame(results)
  head(results2)
  significantTQ<-subset(results2, rownames(results2) %in% rowsTQ)
  head(significantTQ)
  str(significantTQ) # 7 observations 
  
  library(corncob)
  taxonTQ<-otu_to_taxonomy(OTU = row.names(DiffTQHvL), data = physeqfiltered2)
  taxonTQ<-as.data.frame(taxonTQ)
  significantTQ$taxonomy<-taxonTQ$taxon
  #View(significantTQ)
  
  #write.csv(significantTQ, "TQLowVhighAncom.csv")
  #separate out the taxonomy fields in the excel

  
## make a plot for high/low
  
  sigTQ<- read.csv("TQLowVhighAncom.csv")
  str(sigTQ)
  
  significantTerr<-data.frame(sigTQ$beta.TQFactorlow,sigTQ$Phylum, sigTQ$Class, sigTQ$Order, sigTQ$Family, sigTQ$se.TQFactorlow)
  str(significantTerr)
  colnames(significantTerr)<-c("Beta",  "Phylum","Class", "Order", "Family", "se")
  significantTerr$Order<-as.factor(significantTerr$Order)
  significantTerr$Class<-as.factor(significantTerr$Class)
  significantTerr$Phylum<-as.factor(significantTerr$Phylum)
  significantTerr$Family<-as.factor(significantTerr$Family)
  significantTerr$Species<-as.factor(seq(from=1, to=7, by=1))
  
  ## plot of all differentially abundant taxa
  
  ggplot(significantTerr, aes(x= Beta, y=Family, colour=Order, group=Species)) + 
    geom_point(position=position_dodge(width = 0.8), size=2) +
    geom_errorbar(aes(xmin=Beta -se, xmax=Beta + se), width=0,lwd=1, position= position_dodge(width=0.8)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(name="\nLog fold change", breaks=c(-3,-2, -1,0, 1, 2)) + 
    ylab("Family\n")  + geom_vline(xintercept = 0, linetype="dotted", size=1.5) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=22),
          legend.text = element_text(size=20), legend.title = element_text(size=22)) 
  
  
  
### Survival
  
  DiffSurvival<-diffAbund[diffAbund$SurvivedNextSeasonFactorTRUE=="TRUE",]
  str(DiffSurvival)
  Survrows<-row.names(DiffSurvival)
  
  results2<-as.data.frame(results)
  head(results2)
  significantSurv<-subset(results2, rownames(results2) %in% Survrows)
  head(significantSurv)
  str(significantSurv) # 6 observations 
  
  library(corncob)
  taxonSurv<-otu_to_taxonomy(OTU = row.names(DiffSurvival), data = physeqfiltered2)
  taxonSurv<-as.data.frame(taxonSurv)
  significantSurv$taxonomy<-taxonSurv$taxon
  #View(significantSurv)
  
  #write.csv(significantSurv, "SurvivalAncom.csv")
  #separate out the taxonomy fields in the excel
  
  
  ## make a plot for survival
  
  sigSurv<- read.csv("SurvivalAncom.csv")
  str(sigSurv)
  
  significantSurvival<-data.frame(sigSurv$beta.SurvivedNextSeasonFactorTRUE,sigSurv$Phylum, sigSurv$Class, sigSurv$Order, sigSurv$Family, sigSurv$se.SurvivedNextSeasonFactorTRUE)
  str(significantSurvival)
  colnames(significantSurvival)<-c("Beta",  "Phylum","Class", "Order", "Family", "se")
  significantSurvival$Order<-as.factor( significantSurvival$Order)
  significantSurvival$Class<-as.factor( significantSurvival$Class)
  significantSurvival$Phylum<-as.factor( significantSurvival$Phylum)
  significantSurvival$Family<-as.factor( significantSurvival$Family)
  significantSurvival$Species<-as.factor(seq(from=1, to=6, by=1))
  
  ## plot of all differentially abundant taxa
  
  ggplot(significantSurvival, aes(x= Beta, y=Family, colour=Order, group=Species)) + 
    geom_point(position=position_dodge(width = 0.8), size=2) +
    geom_errorbar(aes(xmin=Beta -se, xmax=Beta + se), width=0,lwd=1, position= position_dodge(width=0.8)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(name="\nLog fold change", breaks=c(-3,-2, -1,0, 1, 2)) + 
    ylab("Family\n")  + geom_vline(xintercept = 0, linetype="dotted", size=1.5) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=22),
          legend.text = element_text(size=20), legend.title = element_text(size=22)) 
  
  scale_colour_manual(values=c("goldenrod2","#CC79A7", "#56B4E9",  "aquamarine4", "darkgrey"))
  
  

  
### time of day
  DiffTime<-diffAbund[diffAbund$TimeOfDayBinaryPM=="TRUE",]
  str(DiffTime)
  rows<-row.names(DiffTime)
  
  results2<-as.data.frame(results)
  head(results2)
  significantTime<-subset(results2, rownames(results2) %in% rows)
  head(significantTime)
  str(significantTime) # 0 observations 


  
  
####################  
### MHC analysis ###
####################
  
# Subset to complete MHC cases
  
  clr2Data <- data.frame(sample_data(physeq_clr2))
  str(clr2Data) #255 samples
  
  MHCDataAvailable<- clr2Data[!is.na(clr2Data$Both_MHC_present),]
  MHCDataBothPresent <- MHCDataAvailable[MHCDataAvailable$Both_MHC_present == "Y",]
  str(MHCDataBothPresent) #189 samples/birds
  
  
  SampleIDsToKeep2<- MHCDataBothPresent$sample.id
  
  physeq_clr3<- subset_samples(physeq_clr2, sample.id %in% SampleIDsToKeep2) #189 samples, 794 taxa

  
#function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
#Extract ASV Matrix and Sample Data  
  clr_MatrixMHC<-vegan_otu(physeq_clr3)
  clr_SampleDataMHC<-as(sample_data(physeq_clr3),"data.frame")
  str(clr_SampleDataMHC)

# subset to variables in model and correct format of sample data
  
  MHC_SampleData<- clr_SampleDataMHC[,c("BirdID","sample.id","TubeNoFungal", "MHC1_Diversity", "MHC2_Diversity",
                                               "MHC1_Diversity_2", "MHC2_Diversity_2", "Ase.ua5", "Ase.ua3","Ase.ua6",
                                               "Ase.ua7","Ase.ua9","Ase.ua11","Ase.ua8","Ase.ua4","Ase.ua1","Ase.ua10",
                                               "Ase.dab1", "Ase.dab2",  "Ase.dab3", "Ase.dab4",  "Ase.dab5",
                                               "TLR3","FieldPeriodIDBinary", "SexEstimate","AgeDays","Ageclass",
                                               "TQcorrectedWithAvBA","MinutesSinceSunrise","Hs_obs",
                                                "SurvivedNextSeasonFactor", "SurvivedNextSeason","DateToFreeze")]
  
  str(MHC_SampleData)
  
  MHC_SampleData$SexEstimate <- as.factor(MHC_SampleData$SexEstimate)
  MHC_SampleData$TLR3 <- as.factor(MHC_SampleData$TLR3)
  MHC_SampleData$FieldPeriodIDBinary <- as.factor(MHC_SampleData$FieldPeriodIDBinary)
  MHC_SampleData$TimeOfDayBinary <- as.factor(ifelse(MHC_SampleData$MinutesSinceSunrise < 360, "AM", "PM"))
  MHC_SampleData$AgeYears<- MHC_SampleData$AgeDays/365.25
  MHC_SampleData$SurvivedNextSeasonFactor<- as.factor(MHC_SampleData$SurvivedNextSeasonFactor)
  table(MHC_SampleData$SurvivedNextSeasonFactor) #24 died, 165 survived

  
  summary(glm(SurvivedNextSeason~ MHC1_Diversity,
      data= MHC_SampleData, family=binomial(link=logit))) # MHC 1 diversity not predictive of survival
  
  summary(glm(SurvivedNextSeason~ MHC2_Diversity,
              data= MHC_SampleData, family=binomial(link=logit))) # 2 MHC diversity not predictive of survival

#######################    
### 1) MHC diversity ##
####################### 
  

### PERMANOVA ###
  
  perm <- how(nperm = 9999)
  set.seed(57789)
  
  permanovaMHCdiversity<- adonis2(clr_MatrixMHC ~ AgeYears + SexEstimate + Hs_obs+ TQcorrectedWithAvBA +
                              FieldPeriodIDBinary + TimeOfDayBinary + SurvivedNextSeasonFactor +
                                TLR3 + MHC1_Diversity + MHC2_Diversity +DateToFreeze,
                            data=MHC_SampleData, permutations = perm, method = "euclidean", by= "margin")
  permanovaMHCdiversity
  
  ResultsPermanovaMHCDiv<- data.frame(permanovaMHCdiversity)
  #write.csv(ResultsPermanovaMHCDiv, "PermanovaMHCDiversity.csv")
  
  clr_SampleDataMHC$AgeYears<- clr_SampleDataMHC$AgeDays/365.25
  summary(glm(SurvivedNextSeason ~ AgeYears + SexEstimate + Hs_obs + TQcorrectedWithAvBA + SampleYear+
                MHC1_Diversity, data= clr_SampleDataMHC, family = "binomial"))

  
#betadisper tests to check homogeneity of variance
  
  distMatrixBeta <- vegdist(clr_MatrixMHC, method="euclidean")
  
  BDFieldPeriod<-betadisper(distMatrixBeta,MHC_SampleData$FieldPeriodIDBinary)
  set.seed(2699)
  permutest(BDFieldPeriod, permutations = 9999) 
  boxplot(BDFieldPeriod,ylab="Distance to centroid\n", xlab="\nFieldPeriod", cex.axis=1.5, cex.lab=1.5)

  BDMHCdiv<-betadisper(distMatrixBeta,MHC_SampleData$MHC1_Diversity)
  set.seed(269789)
  permutest(BDMHCdiv, permutations = 9999) 
  boxplot(BDMHCdiv,ylab="Distance to centroid\n", xlab="\nMHCdiversity", cex.axis=1.5, cex.lab=1.5)
  
  BDSurvivalMHC<-betadisper(distMatrixBeta,MHC_SampleData$SurvivedNextSeasonFactor)
  set.seed(269789)
  permutest(BDSurvivalMHC, permutations = 9999) 
  
  
  
  
### envfit ###
  
  clr_pcaMHC<-rda(clr_MatrixMHC)
  
#uses linear model permutations to map variables onto an ordination###
  
  envfitDat<- MHC_SampleData[,c("AgeYears", "SexEstimate", "Hs_obs", "TQcorrectedWithAvBA",
                                "FieldPeriodIDBinary", "TimeOfDayBinary", "SurvivedNextSeasonFactor",
                                "TLR3", "MHC1_Diversity", "MHC2_Diversity", "DateToFreeze")]
  set.seed(67689)
  envfitmod<- envfit(clr_pcaMHC, envfitDat, permutations = 9999, choices=c(1,2))
  envfitmod #relationship with axes 1,2 - MHCI diversity, field period and survival
  
  
  set.seed(67999)
  envfitmod3<- envfit(clr_pcaMHC, envfitDat, permutations = 9999, choices=c(3,4))
  envfitmod3 # TQ associated with shifts along 3 and 4 (primarily 4), field period also associated (primarily shifted along axis 3)

  
  #plotting envfit arrows onto PCA#
  arrow_factor <- ordiArrowMul(envfitmod, fill=0.1) #get proportions, longer arrows= stronger predictor
  spp.scrs <- as.data.frame(scores(envfitmod, display = "vectors"))* arrow_factor
  spp.scrs$Species <- c("AgeYears", "Hs_obs", "TQcorrectedWithAvBA","MHC1_diversity", "MHC2_diversity", "DateToFreeze")
  
  scorePCA<- data.frame(scores(clr_pcaMHC, choices=c(1,2), display = "sites", scaling=1))
  summary(scorePCA)
  colnames(scorePCA)<- c("PC1", "PC2")
  scorePCA$FieldPeriod<- as.factor(MHC_SampleData$FieldPeriodIDBinary)
  scorePCA$MHCIdiversity<- as.factor(MHC_SampleData$MHC1_Diversity)
  scorePCA$Survival<- as.factor(MHC_SampleData$SurvivedNextSeasonFactor)
  
  centroidsPC1S<-data.frame(ddply(scorePCA, .(FieldPeriod), summarize, mean=mean(PC1)))
  centroidsPC2S<-data.frame(ddply(scorePCA, .(FieldPeriod), summarize, mean=mean(PC2)))
  centroidsMHCsurv<- merge(centroidsPC1S, centroidsPC2S, by="FieldPeriod")
  colnames(centroidsMHCsurv)<- c("Survival", "x", "y")
  
  
  EnvPlot12MHC<- ggplot(data = scorePCA, aes(x = PC1, y = PC2)) + theme_bw() + geom_point(size = 2, stroke=2, aes(col=MHCIdiversity, shape=FieldPeriod)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 4.32%") + ylab("PC2 3.66%\n") +
    geom_point(data=centroidsMHCsurv, aes(x,y), shape=21, stroke=4, size=4, alpha=0.9, colour="black")+
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18),
          legend.text = element_text(size=18))+
    coord_fixed()+
    geom_segment(data = spp.scrs,
                 aes(x = 0, xend = PC1, y = 0, yend = PC2),
                 arrow = arrow(length = unit(0.25, "cm")), colour = "grey31", size=1.5)+
    geom_text(data = spp.scrs, aes(label= Species, x = PC1*1.3, y = PC2*1.3),
              size = 5)
  
  EnvPlot12MHC
  
  
### Principal Components Analysis ###

  clr_pcaMHC<-rda(clr_MatrixMHC)
  
  sigMHC <- PCAsignificance (clr_pcaMHC, axes = 8) 
  sigMHC # PC1 = 4.319, PC2 = 3.655, PC3 = 2.796, PC4 = 2.618
  

## PCA of MHC ##
  
  plotPCAMHC<-gg_ordiplot(clr_pcaMHC,groups=MHC_SampleData$MHC1_Diversity,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCAMHC

  
  ord.dataMHC <- plotPCAMHC$df_ord
  ord.dataMHC$MHCIdiversity<- as.factor(MHC_SampleData$MHC1_Diversity)
  ord.dataMHC$MHCdiversity<- as.factor(ifelse(MHC_SampleData$MHC1_Diversity <5,"<5","5+"))
  ord.dataMHC$MHCdiversity2<- as.factor(ifelse(MHC_SampleData$MHC1_Diversity <4, "2or3", 
                                               ifelse(MHC_SampleData$MHC1_Diversity<6, "4or5","6or7")))
  head(ord.dataMHC)

  
  #centroids for PC1,2
  pca_scoresMHC<-scores(clr_pcaMHC, choices=c(1,2,3,4), scaling=1)
  sample_scoresMHC <- data.frame(pca_scoresMHC$sites)
  str(sample_scoresMHC)
  summary(sample_scoresMHC)
  sample_scoresMHC<-data.frame(setDT(sample_scoresMHC, keep.rownames = TRUE)[])
  str(sample_scoresMHC)
  colnames(sample_scoresMHC)<-c("Sample.ID", "PC1", "PC2","PC3", "PC4")
  sample_scoresMHC$MHC<- as.factor(MHC_SampleData$MHC1_Diversity)
  sample_scoresMHC$MHCdiversity<- as.factor(ifelse(MHC_SampleData$MHC1_Diversity <5,"<5","5+"))
 
  centroidsPC1MHC<-data.frame(ddply(sample_scoresMHC, .(MHC), summarize, mean=mean(PC1)))
  centroidsPC2MHC<-data.frame(ddply(sample_scoresMHC, .(MHC), summarize, mean=mean(PC2)))
  centroidsMHC<- merge(centroidsPC1MHC, centroidsPC2MHC, by="MHC")
  colnames(centroidsMHC)<- c("MHCIdiversity", "x", "y")
  
  PlotMHC<- ggplot(data = ord.dataMHC, aes(x = x, y = y, fill = MHCIdiversity)) + 
    theme_bw() + geom_point(size = 4, stroke=1, shape=21, alpha=0.8, colour="grey55")+
    geom_point(data=centroidsMHC, aes(x,y), shape=21, stroke=4, size=6, alpha=0.9, colour="black")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 4.32%") + ylab("PC2 3.66%\n") + 
    scale_fill_manual(values=c("darkseagreen1","darkolivegreen3","darkgreen", "skyblue1", "slateblue", "navyblue")) +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18), legend.text = element_text(size=18))
  
  PlotMHC
  
  
  centroidsPC1MHCb<-data.frame(ddply(sample_scoresMHC, .(MHCdiversity), summarize, mean=mean(PC1)))
  centroidsPC2MHCb<-data.frame(ddply(sample_scoresMHC, .(MHCdiversity), summarize, mean=mean(PC2)))
  centroidsMHCb<- merge(centroidsPC1MHCb, centroidsPC2MHCb, by="MHCdiversity")
  colnames(centroidsMHCb)<- c("MHCdiversity", "x", "y")
  
  PlotMHC<- ggplot(data = ord.dataMHC, aes(x = x, y = y, fill = MHCdiversity)) + 
    theme_bw() + geom_point(size = 4, stroke=1, shape=21, alpha=0.8, colour="grey55")+
    geom_point(data=centroidsMHCb, aes(x,y), shape=21, stroke=4, size=8, alpha=0.9, colour="black")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 4.32%") + ylab("PC2 3.66%\n") + 
    scale_fill_manual(values=c("forestgreen", "lemonchiffon2"), name="MHC diversity") +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22), legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.5, 'cm'))+ guides(fill = guide_legend(byrow = TRUE))
  
  PlotMHC
  
  
  
### Differential abundance testing ###
  #### Ancom BC
  
  
  library(ANCOMBC)
  
  physeqfiltered # unrarefied reads, but filtered to remove rare taxa, not CLR transformed
  phyfilteredDat<- data.frame(sample_data(physeqfiltered))
  
  floatersRemoved<- phyfilteredDat[!is.na(phyfilteredDat$TQcorrectedWithAvBA),]
  str(floatersRemoved) #255 samples/birds
  
  SampleIDsToKeep<- floatersRemoved$sample.id
  
  physeqfiltered2<- subset_samples(physeqfiltered, sample.id %in% SampleIDsToKeep) #255 samples, 794 taxa
  phyfilteredDat2<- data.frame(sample_data(physeqfiltered2))
  
  MHCDataAvailable<- phyfilteredDat2[!is.na(phyfilteredDat2$Both_MHC_present),]
  MHCDataBothPresent <- MHCDataAvailable[MHCDataAvailable$Both_MHC_present == "Y",]
  str(MHCDataBothPresent) #189 samples/birds
  
  SampleIDsToKeep2<- MHCDataBothPresent$sample.id
  
  physeqfiltered3<- subset_samples(physeqfiltered2, sample.id %in% SampleIDsToKeep2) #189 samples, 794 taxa
  
  
  str(sample_data(physeqfiltered3))
  sample_data(physeqfiltered3)$SexEstimate <- as.factor(sample_data(physeqfiltered3)$SexEstimate)
  sample_data(physeqfiltered3)$FieldPeriodIDBinary <- as.factor(sample_data(physeqfiltered3)$FieldPeriodIDBinary)
  sample_data(physeqfiltered3)$SurvivedNextSeasonFactor <- as.factor(sample_data(physeqfiltered3)$SurvivedNextSeasonFactor)
  sample_data(physeqfiltered3)$TimeOfDayBinary<- as.factor(ifelse(sample_data(physeqfiltered3)$MinutesSinceSunrise< 360, "AM", "PM"))
  sample_data(physeqfiltered3)$AgeYears<- sample_data(physeqfiltered3)$AgeDays/365.25
  sample_data(physeqfiltered3)$TQFactor<- as.factor(ifelse(sample_data(physeqfiltered3)$TQcorrectedWithAvBA < 12678, "low", "high"))
  sample_data(physeqfiltered3)$TQFactor2<- as.factor(ifelse(sample_data(physeqfiltered3)$TQcorrectedWithAvBA <7624, "1stQuartile",
                                                            ifelse(sample_data(physeqfiltered3)$TQcorrectedWithAvBA < 25866, "interquartile", "4thQuartile")))
  sample_data(physeqfiltered3)$MHC1diversityFactor<- as.factor(ifelse(sample_data(physeqfiltered3)$MHC1_Diversity <5, "2-4","5-7"))
  sample_data(physeqfiltered3)$MHC1diversity<- as.factor(sample_data(physeqfiltered3)$MHC1_Diversity)
  
  

# running Ancom BC
  
  AncomDA<-ancombc(physeqfiltered3, formula= "MHC1diversityFactor +AgeYears + TQFactor + FieldPeriodIDBinary +
                  SexEstimate + Hs_obs + TimeOfDayBinary + SurvivedNextSeasonFactor + DateToFreeze", p_adj_method = "BH",
                   lib_cut = 0, alpha = 0.05, neg_lb = TRUE, group = "MHC1diversityFactor", struc_zero = TRUE,global = FALSE)
  str(AncomDA)
  
  results<-AncomDA$res
  str(results)
  head(results)
  
  diffAbund<-as.data.frame(results$diff_abn)
  head(diffAbund)
  str(diffAbund)
  
  
  ### MHC diversity
  DiffMHC<-diffAbund[diffAbund$"MHC1diversityFactor5-7" =="TRUE",]
  str(DiffMHC)
  rowsMHC<-row.names(DiffMHC)
  
  results2<-as.data.frame(results)
  head(results2)
  significantMHC<-subset(results2, rownames(results2) %in% rowsMHC)
  head(significantMHC)
  str(significantMHC) 
  
  library(corncob)
  taxonMHC<-otu_to_taxonomy(OTU = row.names(DiffMHC), data = physeqfiltered3)
  taxonMHC<-as.data.frame(taxonMHC)
  str(significantMHC)
  significantMHC$taxonomy<-taxonMHC$taxon
  
  tax<-significantMHC[,c("beta.MHC1diversityFactor5.7","se.MHC1diversityFactor5.7", "taxonomy")]
  View(tax)
  
  table(sample_data(physeqfiltered3)$MHC1diversity)
  table(sample_data(physeqfiltered3)$MHC1diversity)
  table(sample_data(physeqfiltered3)$MHC1diversityFactor2)
  
  #write.csv(significantMHC, "MHCdiversityAncom.csv")
  #separate out the taxonomy fields in the excel
  
  
  
  
#######################    
### 2) MHC Alleles ####
####################### 
  
  str(MHC_SampleData)
  
  MHC_SampleData$Ase.dab3<- as.factor(MHC_SampleData$Ase.dab3)
  MHC_SampleData$Ase.dab4<- as.factor(MHC_SampleData$Ase.dab4)
  MHC_SampleData$Ase.dab5<- as.factor(MHC_SampleData$Ase.dab5)
  MHC_SampleData$Ase.ua1<- as.factor(MHC_SampleData$Ase.ua1)
  MHC_SampleData$Ase.ua3<- as.factor(MHC_SampleData$Ase.ua3)
  MHC_SampleData$Ase.ua4<- as.factor(MHC_SampleData$Ase.ua4)
  MHC_SampleData$Ase.ua5<- as.factor(MHC_SampleData$Ase.ua5)
  MHC_SampleData$Ase.ua6<- as.factor(MHC_SampleData$Ase.ua6)
  MHC_SampleData$Ase.ua7<- as.factor(MHC_SampleData$Ase.ua7)
  MHC_SampleData$Ase.ua8<- as.factor(MHC_SampleData$Ase.ua8)
  MHC_SampleData$Ase.ua9<- as.factor(MHC_SampleData$Ase.ua9)
  MHC_SampleData$Ase.ua10<- as.factor(MHC_SampleData$Ase.ua10)
  MHC_SampleData$Ase.ua11<- as.factor(MHC_SampleData$Ase.ua11)

  
  summary(glm(SurvivedNextSeason~ Ase.dab3 + Ase.dab4 + Ase.dab5 + Ase.ua1 + Ase.ua3 + 
                Ase.ua4 + Ase.ua5 + Ase.ua6 + Ase.ua7 + Ase.ua8 + Ase.ua9 + Ase.ua11,
              data= MHC_SampleData, family=binomial(link=logit))) # not predictive of survival
  
 
  
### PERMANOVA ###
  #remove Ase.ua10 as correlated heavily with Ase.ua1 

  perm <- how(nperm = 9999)
  set.seed(23966)
  permanovaMHCAlleles<- adonis2(clr_MatrixMHC ~ AgeYears + SexEstimate + Hs_obs+ TQcorrectedWithAvBA +
                                    FieldPeriodIDBinary + TimeOfDayBinary + DateToFreeze + SurvivedNextSeasonFactor + TLR3 +
                                    Ase.ua1+ Ase.ua3 + Ase.ua4 + Ase.ua5 + Ase.ua6 + Ase.ua7 + Ase.ua8 + 
                                    Ase.ua9 + Ase.ua11 + Ase.dab3 + Ase.dab4 + Ase.dab5, 
                                  data=MHC_SampleData, permutations = perm, method = "euclidean", by= "margin")
  permanovaMHCAlleles
  
  ResultsPermanovaMHCalleles<- data.frame(permanovaMHCAlleles)
  write.csv(ResultsPermanovaMHCalleles, "PermanovaMHCAlleles.csv")

  
###Principal Components Analysis###
  
  clr_pcaMHC<-rda(clr_MatrixMHC)
  
  screeplot(clr_pcaMHC) # PC1 gives most information, followed by 2.
  
  sigMHC <- PCAsignificance (clr_pcaMHC, axes = 8) 
  sigMHC # PC1 = 4.319, PC2 = 3.655, PC3 = 2.796, PC4 = 2.618

  plotPCA_aseua11<-gg_ordiplot(clr_pcaMHC,groups=clr_SampleDataMHC$Ase.ua11,ellipse = FALSE,plot=FALSE, pt.size=2, scaling=1, choices = c(2,3)) 
  plotPCA_aseua11

  

