This folder contains data and R script needed to reproduce the results in the paper:
Worsley, S.F., Davies, C.S., Mannarelli, ME. et al. Assessing the causes and consequences of gut mycobiome variation in a wild population of the Seychelles warbler. Microbiome 10, 242 (2022). https://doi.org/10.1186/s40168-022-01432-7

All variable names in these data should be self-explanatory, but please don't hesitate to send an email to the corresponding author if anything is unclear.

Files included:

### Scripts ###

#FungalSequences.R
This R script contains the code necessary to reproduce all initial pre-processing steps using outputs from QIIMEII (e.g. decontamination, sample filtering, extraction repeatability tests). The outputs from this code are used in all downstream analysis.

#AlphaDiversityFungal.R
This R script contains the code necessary to 1) calculate the core fungal microbiome 2) compare alpha diversity for bacteria and fungi in the same sample 3) reproduce analyses investigating the association between environmental/host variables and gut microbiome alpha diversity, and 4) investigate the association between fungal alpha diversity and host survival. It uses the filtered ASV, taxonomy and metadata files generated in FungalSequences.R

#FungalBetaDiversity.R
This R script contains the code necessary to reproduce the results of the analysis investigating the association between environment/host variables and gut microbiome beta diversity. It uses the filtered ASV, taxonomy and metadata files generated in FungalSequences.R


### Data ###

# feature_table.csv
This file contains the raw abundances of fungal Amplicon Sequence Variants (ASV) inferred for each faecal sample using QIIMEII. This ASV table is used to generate phyloseq objects. The paired-end sequences for each sample used to generate the table have been deposited with the ENA under the accession number PRJEB54641. 

# taxonomy.csv
This file contains the taxonomy of fungal ASVs in feature_table.csv. It was generated as part of the QIIMEII pipeline. This file is used to generate phyloseq objects.

# FungiMetadata.csv
This file contains the sample metadata. Columns correspond to:

sample.id = Unique identifier for each faecal sample
UniqueID = NBAF sequencing ID combined with unique sample tube number for samples that underwent fungal sequencing
TubeNoFungal = Unique tube number for each faecal sample- prefix "FUN" added to differentiate samples that underwent fungal sequencing (suffix "R" means repeat extraction of the same sample)
OriginalTubeNumberUnique = Unique tube number for each faecal sample given at time of sampling (suffix "R" means repeat extraction of the same sample)
BacteriaOriginalSampleIDUnique = the corresponding unique tube number used in metadata for bacterial sequencing (suffix "R" means repeat extraction of the same sample)
BacterialSequencingSample.ID= NBAF sequencing ID combined with unique sample tube number for samples that underwent bacterial sequencing.
BacterialNBAFSAMPLEID= NBAF sequencing ID for samples that underwent bacterial sequencing.
BirdID = Unique identifier for each bird
BTO = Unique BTO ring number for each Bird
FieldRing = Unique colour ring combination for each Bird
TubeNumber = Tube number for each faecal sample (note there are duplicates if the sample has been sequenced and/or extracted twice- i.e. they are not unique).
SampleDate = date faecal sample taken
SampleType = Sample type (C= control or F= faecal)
SampleType4way = Sample type in further detail. Collection control; Extraction control (blank swab) or PCR negative (water control in PCR); Positive control (mock community); F (faecal sample).
SampleYear = Year faecal sample take (2017-2020)
FreezerDate = date samples were frozen at -80 degrees C.
DateToFreeze = number of days faecal sample stored at 4 degrees C prior to freezing.
ExtractionDate = date of DNA extraction from faecal sample
Qubit = Concentration (ng/μl) of DNA extracted from sample
CatchID = Unique identifier for catch
TerritoryIDSampled = Territory ID where bird was caught/sampled (note not necessarily breed group territory).
BacterialSequenceRunNo = sequence run sample was included in for bacterial sequencing (1-4)
Island = Island where bird was samples (CN: Cousin Island)
FieldPeriodID = Unique identifier for field period sample was taken in (164 = major17, 166= minor18, 167= major18, 171= minor19, 173= major19, 174= minor20)
FieldPeriodIDFactor = Field period sample was taken (Major17, Minor18, Major18, Minor19, Major19, Minor20)
FieldPeriodIDBinary = Breeding season sample was taken in (Minor or Major)
BirthDate = estimated birth date for each bird
SexEstimate = sex of bird (M= male, F= female)
AgeDays = Age in days
Ageclass = ageclass of bird (FL= fledgling, OFL= old fledgling, SA= Sub-adult, A= adults)
Eyecolour= eyecolour of bird (G= grey, LB= light brown, RB= red brown).
BodyMass = Body mass of bird (in grams)
RightTarsus = length of right tarsus (in mm)
EggPresent = an egg present in females (0= no, -1= yes)
BroodPatch = is a brood patch present- suggests incubating (0= no, -1= yes)
Status = status of bird in population census
TerritoryID = Territory ID of that bird's breed group
TQcorrected = quality of territory calculated using the "corrected" method in the database (note there are missing values for some seasons)
TQcorrectedWithAvBA = quality of territory (as above). For territories with missing qualities an average of the preceding and following seasons was taken.
SurvivedNextSeason = Did the bird survive until the following breeding season, binary (0= no/died, 1= yes/survived)
SurvivedNextSeasonFactor = Did the bird survive until the following breeding season (no= died or yes= survived)
CatchTime = time of day bird was caught/sampled
MinutesSinceSunrise = time bird was caught as minutes since sunrise (06:00 am)
Observer = Initials of observer who caught/sampled bird
Hs_obs = genome-wide heterozygosity
Both_MHC_present = data available for both MHC-I and MHC-II? (Y= yes, N = no only one class genotyped, NA = no data for either locus)
MHC1_Diversity = the number of MHC-I alleles carried by an individual
MHC1_Diversity_2 = the number of MHC-I alleles carried by an individual squared
Ase.ua1 to Ase.ua11 = presence/absence of specific MHC-I alleles in an individual (0= absent, 1= present)- only alleles present in > 5% but < 95% of individuals are displayed here.
MHC2_Diversity = the number of MHC-II alleles carried by an individual
MHC2_Diversity_2 = the number of MHC-II alleles carried by an individual squared
Ase.dab1 to Ase.dab5 = presence/absence of specific MHC-II alleles in an individual (0= absent, 1= present)- only alleles present in > 5% but < 95% of individuals are displayed here.
TLR3 = genotype at TLR3 locus (AA, AC, or CC)


# FilteredASVTableBacteria.csv
This file contains the raw abundances of bacterial Amplicon Sequence Variants (ASV) inferred for each faecal sample using QIIMEII analysis. The ASV table is used to generate phyloseq objects as part of the analyses in AlphaDiversityFungal.R

# FilteredTaxonomyTableBacteria.csv
This file contains the taxonomy of bacterial ASVs in FilteredASVTableBacteria.csv. This file is used to generate phyloseq objects as part of the analyses in AlphaDiversityFungal.R

# FilteredMetadataBacteria.csv
This file contains the metadata for samples that underwent 16S V4 amplicon sequencing as part of a previous paper (https://doi.org/10.1186/s42523-021-00149-6). This file is used to generate a bacterial phyloseq objects as part of the analyses in AlphaDiversityFungal.R. Column names can be interpreted in the same way as for FungiMetadata.csv (see above).