#### FungalSequences.R ####
#This R script contains the code necessary to reproduce all initial pre-processing steps using outputs from QIIMEII
#(e.g. decontamination, sample filtering, extraction repeatability tests). 
#The outputs from this code are used in all downstream analysis.

#load packages
  library(ggplot2)
  library(phyloseq)
  library(ape)
  library(iNEXT)
  library(patchwork) 
  library(decontam)  
  library(dplyr)
  library(vegan)
  library(microbiome)

#ASV table
  asv_table <- read.csv ("feature_table.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) #should be 6942 features, 298 samples (2 faecal samplers were filtered out by QIIME due to low numbers of reads)
  head (asv_table)
  asv_table <- as.matrix (asv_table) #make into a matrix

#taxonomy
  taxonomy <- read.csv ("taxonomy.csv", row.names=1)
  str (taxonomy) #should be 6942 observations, 7 taxonomic groupings, feature names as row names.
  taxonomy <- as.matrix (taxonomy)

#load filtered metadata with sample names as rownames
  metadata<-read.csv("FungiMetadata.csv", row.names = 1)
  str(metadata)  
  
  metaF<- metadata[metadata$SampleType4way=="F",]
  table(metaF$SampleYear)
  table(metaF$FieldPeriodIDFactor)
  Div<- metaF[!is.na(metaF$MHC2_Diversity),]
  mean(Div$MHC2_Diversity)
  Div2<- metaF[!is.na(metaF$MHC1_Diversity),]
  mean(Div2$MHC1_Diversity)
  metaF$Yrs<- metaF$AgeDays/365.25
  range(metaF$Yrs)
  
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
  physeq #check there are the right numbers of samples/features/variables etc... 

  table(sample_data(physeq)$SampleType4way)

### 2 faecal samples were lost in qiime filtering steps- 280 F instead of 282, 298 samples instead of 300
### 7 Collection controls (4 hands, 2 bag, 1 tray)
### 4 PCR negatives
### 2 Positive controls - 1 zymo (2 species), 1 mock community (12 ASVs)

#check the number of reads per sample
  sample_reads<-data.frame(reads=sample_sums(physeq))
  head(sample_reads)
  sum(sample_reads$reads)
  mean(sample_reads$reads)
  sd(sample_reads$reads)
  
#Postive controls
  PosControls<- subset_samples(physeq, SampleType4way == "Positive control")
  PosControlsASVs<- data.frame(otu_table(PosControls))
  PosControlsASVs<- merge(PosControlsASVs, taxonomy, by="row.names")
  
  Sample_298<- PosControlsASVs[PosControlsASVs$Sample_298 > 0,] # Cryptococcus and Saccharomyces identified in Zymo mock community
  Sample_299<- PosControlsASVs[PosControlsASVs$Sample_299 > 0,] # 12 ASVs are abundant (read count >50) in synthetic positive control- this contained 12 fungal sequences
  

 
##################  
### FILTERING ####
##################

#### filter to remove instances where features are not assigned at phylum level #####

#summarise number of features at kingdom level- all Fungi
  table(tax_table(physeq)[,"Kingdom"]) 

#print features per phylum- can see that 3389 are phylum unassigned (have a blank name)
  table(tax_table(physeq)[,"Phylum"])
# removes all labelled as unassigned at phylum level
  physeq<-subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "unidentified"))
  table(tax_table(physeq)[,"Phylum"]) #check removed- gives the number of features per phylum

#check adjusted numbers - now 3552 taxa in 298 samples  
  physeq

#remove those with 0 reads following removal of ASVs unassigned to phylum
  sample_data(physeq)$SampleReads<- sample_sums(physeq)
  DataPhyseq<- data.frame(sample_data(physeq))
  ZeroCounts<- DataPhyseq[DataPhyseq$SampleReads==0,]
  ZeroCounts #6 samples now have zero reads (all had low original read counts- 3 PCR negatives, 1 extraction control, 1 positive control, 1 faecal sample)
  ZeroCountIDs<- row.names(ZeroCounts)
  
  physeq2<- subset_samples(physeq, !(sample.id %in% ZeroCountIDs)) # 292 samples, 3552 taxa
  table(sample_data(physeq2)$SampleType4way)
  
  
### DECONTAM ###
  
# package decontam to remove ASVs found in negative controls (lab or collection controls)

#plot read counts per sample
  df <- as.data.frame(sample_data(physeq2)) # Put sample_data into a data.frame
  df$LibrarySize <- sample_sums(physeq2)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType4way)) + geom_point(size=0.4)
  #extraction and PCR controls fall at lower end of library size scale
  #collection controls are distributed throughout
  
## Filter based on prevalence ##

# 1) Filter based on ASVs in extraction controls/PCR negatives- 5 remaining
  
  physeq2 #3552 ASVs
  sample_data(physeq2)$is.negLab <- ifelse(sample_data(physeq2)$SampleType4way == "Extraction control", TRUE, 
                                           ifelse(sample_data(physeq2)$SampleType4way == "PCRNegative", TRUE, FALSE))

  
  
  decontamdf.prev.Lab <- isContaminant(physeq2, method="prevalence", neg="is.negLab")
  table(decontamdf.prev.Lab$contaminant) # 7 contaminants identified
  ContamsLabControls<- decontamdf.prev.Lab[decontamdf.prev.Lab$contaminant=="TRUE",]
  
  ContamsLabControls<- merge(ContamsLabControls, taxonomy, by="row.names")
  str(ContamsLabControls)
  #write.csv(ContamsLabControls, "ContaminantsLabControls.csv") 
  # several Candida ASVs- can contaminate reagents
  

# Make phyloseq object of presence-absence in negative controls and true samples
  ps.pa <- transform_sample_counts(physeq2, function(abund) 1*(abund>0))
  ps.pa.Lab <- prune_samples(sample_data(ps.pa)$is.negLab== "TRUE", ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.negLab== "FALSE", ps.pa)
  
# Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.Lab),
                      contaminant=decontamdf.prev.Lab$contaminant)
  ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
#Remove contaminating taxa
  taxaNames<- taxa_names(physeq2)
  str(taxaNames)
  IDsLabContam<- as.vector(ContamsLabControls$Row.names)
  taxaNames2<- taxaNames[!taxaNames %in% IDsLabContam]
  str(taxaNames2)
  physeq3<-prune_taxa(taxaNames2, physeq2) 
  physeq3 #3545 ASVs
  
  sample_sums(physeq3) # note 1 PCR control now has 0 reads

  
  
# 2) Remove ASVs that are collection contaminants
  
  sample_data(physeq3)$is.negCollection <- sample_data(physeq3)$SampleType4way == "Collection control"
  
  decontamdf.prev.Collection <- isContaminant(physeq3, method="prevalence", neg="is.negCollection")
  table(decontamdf.prev.Collection$contaminant) #104 contaminants
  ContamsCollectionControls<- decontamdf.prev.Collection[decontamdf.prev.Collection$contaminant=="TRUE",]
  
  ContamsCollectionControls<- merge(ContamsCollectionControls, taxonomy, by="row.names")
  str(ContamsCollectionControls)
  #write.csv(ContamsCollectionControls, "ContaminantsCollectionControls.csv")
  
# Make phyloseq object of presence-absence in negative controls and true samples
  ps.pa.Collection1 <- transform_sample_counts(physeq3, function(abund) 1*(abund>0))
  ps.pa.Collection2 <- prune_samples(sample_data(ps.pa.Collection1)$is.negCollection== "TRUE", ps.pa.Collection1)
  ps.pa.posCollection <- prune_samples(sample_data(ps.pa.Collection1)$is.negLab== "FALSE", ps.pa.Collection1)
  
  # Make data.frame of prevalence in positive and negative samples
  df.paCollection <- data.frame(pa.pos=taxa_sums(ps.pa.posCollection), pa.neg=taxa_sums(ps.pa.Collection2),
                      contaminant=decontamdf.prev.Collection$contaminant)
  ggplot(data=df.paCollection, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
#Remove contaminating taxa as all seem reasonable (note the most prevalent one tends to be in samples taken from trousers/cotton bags/hands, or tiny samples)
  taxaNamesAll<- taxa_names(physeq3)
  str(taxaNamesAll)
  IDsCollectionContam<- as.vector(ContamsCollectionControls$Row.names)
  taxaNamesToKeep<- taxaNamesAll[!taxaNamesAll %in% IDsCollectionContam]
  str(taxaNamesToKeep)
  physeq4<-prune_taxa(taxaNamesToKeep, physeq3) 
  physeq4 #(3441 ASVs)
  
  sample_sums(physeq4)
  

#remove samples with 0 reads following decontam
  sample_data(physeq4)$SampleReads<- sample_sums(physeq4)
  DataPhyseq<- data.frame(sample_data(physeq4))
  ZeroCounts<- DataPhyseq[DataPhyseq$SampleReads==0,]
  ZeroCounts #1 sample now has zero reads (this is a PCR negative)
  ZeroCountIDs<- row.names(ZeroCounts)
  
  physeq5<- subset_samples(physeq4, !(sample.id %in% ZeroCountIDs)) # 291 samples, 3441 taxa
  table(sample_data(physeq5)$SampleType4way)
  sample_sums(physeq5)
  

#### remove control samples
  
  physeq6<- subset_samples(physeq5, SampleType == "F") 
  physeq6 # 3441 ASVs 279 faecal samples
  table(sample_data(physeq6)$SampleType4way)
  
  
#### SAMPLE COMPLETENESS/alpha rarefaction plotting (species accumulation) using iNEXT #######
  
  head(otu_table(physeq6))
  
  #make the otu abundance table into a matrix and check by printing first 2 rows/columns
  abund <- as(otu_table(physeq6), "matrix") 
  abund[1:2,1:2]
  
  #convert to a dataframe
  abund2 <- as.data.frame(abund)
  str(abund2)
  
  #iNEXT only takes numeric values, so change all values in the dataframe to numeric values instead of integers.
  df2 <- mutate_all(abund2, function(x) as.numeric(x)) 
  str(df2)
  
# Install and load inext
  #iNEXT= https://cran.r-project.org/web/packages/iNEXT/iNEXT.pdf
  
#Run inext
  #q=0 specifies that the function should use species richness (rather than shannon (1) or simpson (2) indices) for the rarefaction. 
  #Datatype= abundance because you have raw abundances. 
  #endpoint=20000 specifies the number of reads (sample size) that should be used as an endpoint for the rarefaction/extrapolation. 
  
  inext<-iNEXT(df2, q=0, datatype="abundance", endpoint=20000)
  
  #plot rarefaction curve
  rarefaction<- ggiNEXT(inext, type=1, se=TRUE, facet.var="none", grey= TRUE) + theme(legend.position = "none")+ xlab("Sequencing depth") + ylab("Observed ASVs")
  rarefaction
  
  #Plot sample completeness curve
  completeness<-ggiNEXT(inext, type=2, se=TRUE, facet.var="none", grey=TRUE)+scale_shape_manual(values=rep(20,45))+ theme(legend.position = "none") +xlab("Read count") +ylab("Sample completeness")
  completeness + geom_vline(xintercept=5000, alpha=0.5, linetype=2)
  
  
## Identify samples with <5000 reads
  SampleReadsphyseq6<- data.frame(sample_sums(physeq6))
  SampleReadsphyseq6$Sample.ID<- row.names(SampleReadsphyseq6) 
  colnames(SampleReadsphyseq6)<- c("reads","sample")
  SampleReadsphyseq6 <-  SampleReadsphyseq6[order(SampleReadsphyseq6$reads),]
  SampleReadsphyseq6[c(1:30),]
  samplesLowReads<- SampleReadsphyseq6[SampleReadsphyseq6$reads <5000,]
  samplesLowReadsID<- samplesLowReads$sample
  
  physeq6Data<- data.frame(sample_data(physeq6))
  physeq6DataLowReads<- physeq6Data[physeq6Data$sample.id %in% samplesLowReadsID,] # almost all samples had very little material- see extraction notes
  physeq6DataLowReads$SampleExtractionNotes
  physeq6DataLowReads$Qubit
  physeq6DataLowReads$SampleReads
  

### remove samples with less than 5000 reads 
  physeq7<-prune_samples(sample_sums(physeq6)>=5000, physeq6)
  physeq7 # removes14 samples - 265 samples, 3441 taxa
  
  
### Filter out ASVs with fewer than 50 reads in total across all samples (these might be spurious as present in positive controls)
  
# Compute prevalence of each feature (total number of samples in which a taxon appears at least once)
  prevdf = apply(X = otu_table(physeq7),
                 MARGIN = ifelse(taxa_are_rows(physeq7), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq7),
                      tax_table(physeq7))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
  # Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq7, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq7),color=Phylum)) +
    # Include threshold: total abundance across all samples set to 50
    geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
# Define abundance threshold as 50 total reads across samples
  
  abundanceThreshold<-50
  
# Execute the filter, using `prune_taxa()` function
  
  head(prevdf1)
  
  KeepTaxa1<-prevdf1[prevdf1$TotalAbundance>=abundanceThreshold,]
  str(KeepTaxa1)
  head(KeepTaxa1)
  KeepTaxaNames<- rownames(KeepTaxa1)
  
  physeq8<- prune_taxa(KeepTaxaNames, physeq7)
  physeq8 # 2540 taxa (down from 3441)
  
 
#Find mean number of ASVs across samples prior to rarefaction
  
  range(sample_sums(physeq8)) #range of read numbers across samples: 5013: 1204831
  observationThreshold <- 1
  sumASVs<- (apply(otu_table(physeq8) > observationThreshold, 2, sum))
  str(sumASVs)
  sumASVs <- data.frame(sumASVs)
  str(sumASVs)
  mean(sumASVs[,1]) #mean ASVs per sample = 51.69
  sd(sumASVs[,1]) #Standard deviation of ASVs per sample = 27.90
  range(sumASVs[,1]) # range 3-160
  
  ReadsAfterFiltering<- sample_sums(physeq8)
  head(ReadsAfterFiltering)
  length(ReadsAfterFiltering)
  range(ReadsAfterFiltering) #5013, 1204831
  mean(ReadsAfterFiltering) #66085.38
  sd(ReadsAfterFiltering) #77067.16

  
### extract the filtered taxonomy,asv and metadata files for other downstream analysis
  
  filteredMeta<- data.frame(sample_data(physeq8))
  #write.csv(filteredMeta, "FilteredMetadata.csv")
  
  filteredASV<- data.frame(otu_table(physeq8))
  #write.csv(filteredASV, "FilteredASVTable.csv")
  
  filteredTax<- data.frame(tax_table(physeq8))
  #write.csv(filteredTax, "FilteredTaxonomyTable.csv")

  
#################################################  
###Repeatability of alpha diversity measures ####
#################################################  
  
  
##### RAREFY READS TO MIN SAMPLING DEPTH: 5000 reads ######
  
  #rarefy to 5000 and set seed before rarefying (89754), so that results are reproducible  
  physeqRare<-rarefy_even_depth(physeq8, 5000, rngseed = 89754)
  
  # 10 ASVs removed after subsampling- leaves 2530 ASVs across 265 samples
  physeqRare
  sample_sums(physeqRare)
  

##### CALCULATE ALPHA DIVERSITY METRICS ####
  
  # calculate shannon diversity using estimate_richness()
  richnessEstRare<-estimate_richness(physeqRare, split=TRUE, measures= c("Chao1", "Shannon", "Observed"))
  head(richnessEstRare)
  str(richnessEstRare)
  
  
  #add alpha diversity metrics to metadata
  physeqRareMeta <- as.data.frame(sample_data(physeqRare))
  head(physeqRareMeta)
  physeqRareMeta$Chao1 <- richnessEstRare$Chao1
  physeqRareMeta$Shannon <- richnessEstRare$Shannon
  physeqRareMeta$Observed <- richnessEstRare$Observed
  
  #write.csv(physeqRareMeta, "MetaForMatrix.csv")

### IDENTIFY SAMPLES EXTRACTED TWICE ######
  
  MetaForMatrix<- read.csv("MetaForMatrix.csv", row.names = 1)
  head(MetaForMatrix)
  str(MetaForMatrix)

  #extract those extracted twice - add these as extraction duplicates to excel
  row_namesExDup <- as.vector(grep("R", MetaForMatrix$TubeNoFungal, value=TRUE))
  str(row_namesExDup) #5 extraction reps
  
  #Filter metadata to only include duplicated tube numbers  
  
  DupsMeta<- MetaForMatrix[MetaForMatrix$TubeNoFungal %in% row_namesExDup,]
  TubeNosDups<-as.vector(DupsMeta$TubeNumber) # get tube numbers of all duplicates
  
  TubesForMatrix<- MetaForMatrix[MetaForMatrix$TubeNumber %in% TubeNosDups,] #subset to just those samples that are duplicated
  str(TubesForMatrix) #10 tube numbers (5 samples each extracted twice)
  Tubenos<- TubesForMatrix[,c(2,5)]
  str(Tubenos)
  

###### calculate the repeatability of alpha diversity metrics #####
  
  #Make a distance matrix with euclidean distances of richness
  Observed <- TubesForMatrix$Observed
  samples <- TubesForMatrix$OriginalTubeNumberUnique
  names(Observed) <- samples

  distMatrix <- vegdist(Observed, method="euclidean")
  distMatrix <- as.matrix(distMatrix)
  distMatrix[1:3,1:3]
  
  
  #extract as a dataframe (just upper triangle- one way comparisons)
  pairwiseDist <- t(combn(colnames(distMatrix), 2))
  pairwiseDist<- data.frame(pairwiseDist, dist=distMatrix[pairwiseDist])
  head(pairwiseDist)
  colnames(pairwiseDist)<- c("ID1","ID2", "dist")
  str(pairwiseDist)
  
  #Copy the two Id columns so that we can check for matches: ID1extr, ID2extr
  pairwiseDist$ID1extr<- pairwiseDist$ID1
  pairwiseDist$ID2extr<- pairwiseDist$ID2
  str(pairwiseDist)
  
  #remove "R" from ID1extr and ID2extr to check for extraction replicates (note, the second extraction has R after the tube number)
  pairwiseDist$ID1extr<-gsub("R","",as.character(pairwiseDist$ID1extr))
  pairwiseDist$ID2extr<-gsub("R","",as.character(pairwiseDist$ID2extr))
  str(pairwiseDist)
  
  #name type of comparison
  pairwiseDist$Comparison<-as.factor(ifelse(pairwiseDist$ID1extr==pairwiseDist$ID2extr, "ExDuplicate", "Between"))
  str(pairwiseDist)
  
  #plot the distances
  samplePairs<- c("Different samples", "Extraction duplicates")
  AlphaPlot<- ggplot(pairwiseDist, aes(Comparison, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2, stroke=1)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4, stroke=2) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=samplePairs) + theme(element_blank())

# Kruskal-Wallis test
  
  kruskal.test(dist ~ Comparison, data = pairwiseDist)

  
  
#######################################
##### Beta diversity repeatability ####
#######################################
  
  physeq8 #use the unrarefied reads
  
# filter rare taxa from phyloseq object (those with low prevalence)
  
  # Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(physeq8),
                 MARGIN = ifelse(taxa_are_rows(physeq8), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq8),
                      tax_table(physeq8))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum

  # Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq8, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq8),color=Phylum)) +
    # Include a guess for the threshold- here 1% = approx 3 samples
    geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
# Define prevalence threshold as 1% samples and abundance as 50 total reads across samples
  prevalenceThreshold<-1*(265/100)
  
# Execute the prevalence filter, using `prune_taxa()` function
  head(prevdf1)
  
  KeepTaxa1prev<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1prev)
  ASVNamesToKeep<- rownames(KeepTaxa1prev)
  
  physeq8
  physeqfiltered<- prune_taxa(ASVNamesToKeep, physeq8)
  physeqfiltered # 809 taxa, 265 samples
  
  
#Filter the phyloseq object to contain just the sequencing and repeat extractions
  
  TubeNosDups #tube nos of seq dups and repeat extractions
  
  physeqFilteredDups <- prune_samples(sample_data(physeqfiltered)$TubeNumber %in% TubeNosDups, physeqfiltered)
  physeqFilteredDups #10 samples (5 extraction duplicates)
  
#CLR transformation of ASV abundances and distance matrix
  
  physeq_clrDups <- microbiome::transform(physeqFilteredDups, "clr")
  
  #function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  #Extract OTU Matrix and Sample Data  
  clr_v<-vegan_otu(physeq_clrDups)
  clr_s<-as(sample_data(physeq_clrDups),"data.frame")  
  
  distMatrixBeta <- vegdist(clr_v, method="euclidean")
  str(distMatrixBeta)
  distMatrixBeta <- as.matrix(distMatrixBeta)
  distMatrixBeta[1:3,1:3]
  
#extract as a dataframe (just upper triangle- one way comparisons)
  pairwiseDistBeta <- t(combn(colnames(distMatrixBeta), 2))
  head(pairwiseDistBeta)
  pairwiseDistBeta<- data.frame(pairwiseDistBeta, dist=distMatrixBeta[pairwiseDistBeta])
  head(pairwiseDistBeta)
  str(pairwiseDistBeta)
  
#match sample ids to unique tube IDs
  TubeMeta<- data.frame(sample_data(physeq_clrDups))
  str(TubeMeta)
  TubeMetaNos<- TubeMeta[,c("sample.id","OriginalTubeNumberUnique")]
  str(TubeMetaNos)
  
  colnames(pairwiseDistBeta)<- c("sample.id","ID2", "dist")
  pairwiseDistBeta2<- merge(pairwiseDistBeta, TubeMetaNos, by="sample.id")
  str(pairwiseDistBeta2)
  
  colnames(pairwiseDistBeta2)<- c("ID1","sample.id", "dist", "TubeID1")
  pairwiseDistBeta3<- merge(pairwiseDistBeta2, TubeMetaNos, by="sample.id")
  str(pairwiseDistBeta3)
  colnames(pairwiseDistBeta3)<- c("ID2","ID1", "dist", "TubeID1", "TubeID2")
  head(pairwiseDistBeta3)
  
#Copy the two Id columns so that we can check for matches between extraction dups:
  pairwiseDistBeta3$ID1extr<- pairwiseDistBeta3$TubeID1
  pairwiseDistBeta3$ID2extr<- pairwiseDistBeta3$TubeID2
  str(pairwiseDistBeta3)
  
  #remove "R" from ID1extr and ID2extr to check for extraction replicates (the second extraction has R after the tube number)
  pairwiseDistBeta3$ID1extr<-gsub("R","",as.character(pairwiseDistBeta3$ID1extr))
  pairwiseDistBeta3$ID2extr<-gsub("R","",as.character(pairwiseDistBeta3$ID2extr))
  head(pairwiseDistBeta3)
  
  #name type of comparison
  pairwiseDistBeta3$Comparison<-as.factor(ifelse(pairwiseDistBeta3$ID1extr==pairwiseDistBeta3$ID2extr, "ExDuplicate", "Between"))
  
#plot of within and between samples
  BetaPlot<-ggplot(pairwiseDistBeta3, aes(Comparison, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2, stroke=1)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4, stroke=2) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() +
    theme(axis.title.x = element_text(size=20), axis.text = element_text(size=16), axis.title.y = element_blank()) +
    theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=samplePairs) + theme(element_blank())

  
# Kruskal-Wallis test
  
  kruskal.test(dist ~ Comparison, data = pairwiseDistBeta3)


# combined plot  
  
  (AlphaPlot|BetaPlot)


