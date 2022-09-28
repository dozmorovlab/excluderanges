###
### =========================================================================
### excluderanges is an AnnotationHub package that stores genomic coordinates of
### problematic genomic regions as GRanges.
### -------------------------------------------------------------------------
###
# SECOND VERSION: AUTOMATIC DOWNLOAD, MORE OBJECTS, GENOME ASSEMBLIES, ORGANISMS.
# INCLUDES T2T AND MM39 EXCLUDERANGES GENERATED WITH BOYLE-LAB/BLACKLIST.
# USES STANDARD CHROMOSOMES (AUTOSOMAL, PRIMARY ASSEMBLY).
# OBJECTS COMMONLY GENERATED IN THE FIRST AND SECOND VERSIONS ARE IDENTICAL FOR
# STANDARD CHROMOSOMES.

# 02a_preprocessing.Rmd
# Creates common (non-gap) excluderanges (classic blacklist regions)
# 02b_preprocessing.Rmd
# Creates gap excluderanges (problematic genomic regions from the UCSC database)

library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(readr)
library(readxl)

# Project folder path
dir_data <- "/Users/mdozmorov/Documents/Work/GitHub/excluderanges.dev/test"
# Results folder, create if not exist
dir_results <- file.path(dir_data, "excludableSets_bed")
if (!dir.exists(dir_results)) dir.create(dir_results)
# Create subdirs for specific results
if (!dir.exists(file.path(dir_results, "bed"))) dir.create(file.path(dir_results, "bed"))
if (!dir.exists(file.path(dir_results, "rds"))) dir.create(file.path(dir_results, "rds"))
# Increase download file timeout 
# https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r
options(timeout=100)

# Download T2T and mm39 manually generated lists from Google Drive
download.file(url = "https://drive.google.com/uc?export=download&id=1xxQwn5meQ0U0AtatDmW11V98LqSVqrCS", destfile = file.path(dir_results, "T2T.excluderanges.bed.gz"))
download.file(url = "https://drive.google.com/uc?export=download&id=17LSNuRTwg5RQCyTLK3GEf8Am1mX_3nyB", destfile = file.path(dir_results, "mm39.excluderanges.bed.gz"))

# Manually created matrix with information for download
mtx <- rbind(
  ### HUMAN excludable sets
  c("T2T.excluderanges", "T2T", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2022", "excluderanges", "excluderanges"), # checked
  c("hg38.Kundaje.GRCh38_unified_Excludable", "hg38", "Defined as a combination of hg38.Lareau.hg38_peaks, hg38.Boyle.hg38-Excludable.v2, and hg38.Wimberley.peakPass60Perc_sorted, followed by manual curation. Supersedes hg38.Kundaje.GRCh38.Excludable.", "2020", "ENCODE", "ENCFF356LFX"), # checked, Kundaje-unified, ENCFF356LFX, added 2020-05-05. checked twice.
  c("hg38.Bernstein.Mint_Excludable_GRCh38", "hg38", "Defined from Mint-ChIP (low input, multiplexed ChIP-seq) data", "2019", "ENCODE", "ENCFF023CZC"), # checked -- Bernstein, ENCFF023CZC, added 2019-07-22. checked twice.
  c("hg38.Boyle.hg38-Excludable.v2", "hg38", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"), # checked -- Boyle, commit Nov 13, 2018. checked twice.
  c("hg38.Kundaje.GRCh38.Excludable", "hg38", "Defined by Anshul Kundaje as a part of ENCODE and modENCODE consortia", "2016", "ENCODE", "ENCFF419RSJ"), # checked, Kundaje, ENCFF419RSJ, added 2016-10-16. checked twice.
  c("hg38.Lareau.hg38.full.Excludable", "hg38", "ENCODE excludable regions combined with regions of high homology to mtDNA (NUMT regions)", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/combinedBlacklist/hg38.full.blacklist.bed"), # checked -- Lareau-full, commit May 29, 2017. checked twice.
  c("hg38.Reddy.wgEncodeDacMapabilityConsensusExcludable.hg38", "hg38", "Defined by the ENCODE consortium, includes satellite repeats (CATTC, GAATG, GAGTG, ACRO1), RepeatMasker repeats (ALR/Alpha, BSR/Beta), centromeric repeats, chrM, High/Low mappability islands. Has extra chromosomes, use keepStandardChromosomes() filtering", "2016", "ENCODE", "ENCFF220FIN"), # checked -- Reddy, ENCFF220FIN, added 2016-03-11. checked twice.
  c("hg38.Wimberley.peakPass60Perc_sorted", "hg38", "Defined by the ewimberley/peakPass software", "2021", "GitHub", "https://github.com/ewimberley/peakPass/raw/main/excludedlists/hg38/peakPass60Perc_sorted.bed"), # checked -- Wimberley, peakpass, hg38, commit Feb 11, 2021. checked twice
  c("hg38.Wold.hg38mitoExcludable", "hg38", "Definition method unknown", "2016", "ENCODE", "ENCFF940NTE"), # checked -- Wold, ENCFF940NTE, added 2016-04-21. checked twice.
  c("hg38.Yeo.eCLIP_Excludableregions.hg38liftover.bed.fixed", "hg38", "Defined from eCLIP data", "2019", "ENCODE", "ENCFF269URO"), # checked -- Yeo, ENCFF269URO, added 2019-11-08. checked twice.
  c("hg19.Boyle.hg19-Excludable.v2", "hg19", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz?raw=true"), # checked -- Boyle, commit Nov 13, 2018. checked twice.
  c("hg19.Bernstein.Mint_Excludable_hg19", "hg19", "Defined from Mint-ChIP (low input, multiplexed ChIP-seq) data", "2019", "ENCODE", "ENCFF200UUD"), # checked-- Bernstein, ENCFF200UUD, added 2019-07-22. checked twice.
  c("hg19.Birney.wgEncodeDacMapabilityConsensusExcludable", "hg19","Defined by the ENCODE consortium, includes satellite repeats (CATTC, GAATG, GAGTG, ACRO1), RepeatMasker repeats (ALR/Alpha, BSR/Beta), centromeric repeats, chrM, High/Low mappability islands", "2011", "ENCODE", "ENCFF001TDO"), # checked -- Birney, ENCFF001TDO, added 2011-05-04. checked twice.
  c("hg19.Crawford.wgEncodeDukeMapabilityRegionsExcludable", "hg19", "Defined by the ENCODE consortium, includes satellite repeats (CATTC, GAATG, GAGTG, ACRO1), RepeatMasker repeats (ALR/Alpha, BSR/Beta), human satellite repeat HSATII, chrM, ribosomal subunit consensus sequences LSU-rRNA_Hsa, SSU-rRNA_Hsa. Has extra chromosomes, use keepStandardChromosomes() filtering", "2011", "ENCODE", "ENCFF001THR"), # checked -- Crawford, ENCFF001THR, added 2011-03-28. checked twice.
  c("hg19.Lareau.hg19.full.Excludable", "hg19", "ENCODE excludable regions combined with regions of high homology to mtDNA (NUMT regions)", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/combinedBlacklist/hg19.full.blacklist.bed"), # checked -- Lareau-full, commit May 29, 2017. checked twice.
  c("hg19.Wold.hg19mitoExcludable", "hg19", "Definition method unknown", "2016", "ENCODE", "ENCFF055QTV"), # checked -- Wold, ENCFF055QTV, added 2016-04-21. checked twice.
  c("hg19.Yeo.eCLIP_Excludableregions.hg19", "hg19", "Defined from eCLIP data, includes skyscraper, rRNA pseudogene, unreliably mapped satellite repeat, and low complexity skyscraper peak regions", "2019", "ENCODE", "ENCFF039QTN"), # checked -- Yeo, ENCFF039QTN, added 2019-11-08. checked twice.
  
  ### MOUSE excludable sets
  c("mm39.excluderanges", "mm39", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2022", "excluderanges", "excluderanges"), # checked
  c("mm10.Boyle.mm10-Excludable.v2", "mm10", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz"), # checked -- Boyle, mm10, commit Nov 19, 2018. checked twice.
  c("mm10.Hardison.Excludable.full", "mm10", "Definition method unknown", "2016", "ENCODE", "ENCFF790DJT"), # checked -- Hardison, ENCFF790DJT, added 2016-03-28. checked twice.
  
  c("mm10.Hardison.psuExcludable.mm10", "mm10", "Definition method unknown", "2016", "ENCODE", "ENCFF226BDM"), # checked -- Hardison, ENCFF226BDM, added 2016-03-28. checked twice.
  c("mm10.Kundaje.anshul.Excludable.mm10", "mm10", "Defined by Anshul Kundaje as a part of ENCODE and modENCODE consortia", "2016", "ENCODE", "ENCFF999QPV"), # checked -- Kundaje, ENCFF999QPV, added 2016-03-28. checked twice.
  c("mm10.Kundaje.mm10.Excludable", "mm10", "Defined by Anshul Kundaje as a part of ENCODE and modENCODE consortia", "2016", "ENCODE", "ENCFF547MET"), # checked -- Kundaje, ENCFF547MET, added 2016-10-16. checked twice.
  c("mm10.Lareau.mm10.full.Excludable", "mm10", "ENCODE excludable regions combined with regions of high homology to mtDNA (NUMT regions)", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/combinedBlacklist/mm10.full.blacklist.bed"), # checked -- Lareau-full mm10, commit May 29, 2017. checked twice.
  c("mm10.Wold.mm10mitoExcludable", "mm10", "Definition method unknown", "2016", "ENCODE", "ENCFF759PJK"), # checked -- Wold, ENCFF759PJK, added 2016-04-21. checked twice.
  c("mm9.Lareau.mm9.full.Excludable", "mm9", "ENCODE excludable regions combined with regions of high homology to mtDNA (NUMT regions)", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/combinedBlacklist/mm9.full.blacklist.bed"), # checked -- Lareau-full mm9, commit May 29, 2017. checked twice.
  c("mm9.Wold.mm9mitoExcludable", "mm9", "Definition method unknown", "2016", "ENCODE", "ENCFF299EZH"), # checked -- Wold, ENCFF299EZH, add 2016-04-21. checked twice.
  
  ### WORM excludable sets
  c("ce11.Boyle.ce11-Excludable.v2", "ce11", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz"), # checked -- Boyle ce11, commit Nov 19, 2018. checked twice.
  c("ce10.Boyle.ce10-Excludable.v2", "ce10", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce10-blacklist.v2.bed.gz"), # checked -- Boyle ce10, commit Nov 19, 2018. checked twice.
  c("ce10.Kundaje.ce10-Excludable", "ce10", "Defined by Anshul Kundaje, superseded by ce10.Boyle.ce10-Excludable.v2", "2012", "Stanford.edu", "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz"), # checked, Kundaje ce10, last modified 2012-11-13. checked twice.
  
  ### ZEBRAFISH excludable sets
  c("danRer10.Domingues.Excludableed", "danRer10", "Defined manually using total RNA-seq.", "2020", "GitHub", "https://github.com/adomingues/redl_domingues_et_al_dev_2020/raw/main/blacklisted.bed"), # checked -- Domingues danRer10, commit Dec 27, 2020. checked twice.
  c("danRer10.Yang.Supplemental_Table_19.ChIP-seq_black_list_in_the_zebrafish_genome", "danRer10", "Defined via MACS2 peak calling using ChIP-seq (PMID: 33239788)", "2020", "Publication", "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8183574/bin/NIHMS1695157-supplement-Supplementary_Table_1-19.zip"), # checked -- published 2020 Nov 25. checked twice.
  
  ### FLY excludable sets
  c("dm6.Boyle.dm6-Excludable.v2", "dm6", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/dm6-blacklist.v2.bed.gz"), # checked -- Boyle dm6, commit Nov 19, 2018. checked twice.
  c("dm3.Boyle.dm3-Excludable.v2", "dm3", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions", "2018", "GitHub", "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/dm3-blacklist.v2.bed.gz"), # checked -- Boyle dm3, commit Nov 19, 2018. checked twice.
  c("dm3.Kundaje.dm3-Excludable", "dm3", "Defined by Anshul Kundaje. Contains heterochromatin chromosomes chr2LHet. Superseded by dm3.Boyle.dm3-Excludable.v2", "2012", "Stanford.edu", "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz"), # checked -- Kundaje, dm3, last modified 2012-11-20. checked twice.
  
  ### PLANT excludable sets
  c("TAIR10.Wimberley.predicted_excluded_list_sorted_0.6", "TAIR10", "Defined by the ewimberley/peakPass software", "2021", "GitHub", "https://github.com/ewimberley/peakPass/raw/main/excludedlists/tair10/predicted_excluded_list_sorted_0.6.bed"), # checked, Wimberley peakpass TAIR10, commit Feb 11, 2021. checked twice.
  c("TAIR10.Klasfeld.arabidopsis_Excludable_20inputs", "TAIR10", "Defined by the Boyle-Lab/Blacklist software, High Signal and Low Mappability regions (DOI: 10.1101/2022.02.27.482177)", "2021", "GitHub", "https://github.com/sklasfeld/GreenscreenProject/raw/main/data/arabidopsis_blacklist_20inputs.bed"), # checked, Klasfeld TAIR10 Blacklist, commit Apr 29, 2021. checked twice.
  c("TAIR10.Klasfeld.arabidopsis_greenscreen_20inputs", "TAIR10", "Defined by the green screen pipeline (DOI: 10.1101/2022.02.27.482177)", "2021", "GitHub", "https://github.com/sklasfeld/GreenscreenProject/raw/main/data/arabidopsis_greenscreen_20inputs.bed"), # checked, Klasfeld TAIR10 Greenscreen, commit Apr 29, 2021. checked twice.
  
  ### HUMAN NUMT sets
  c("T2T.Lareau.chm13v2.0_peaks", "T2T", "Regions of high homology to mtDNA (NUMT regions) defined by caleblareau/mitoblacklist", "2022", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/peaks/chm13v2.0_peaks.narrowPeak"), # checked -- Lareau-peaks, commit Jun 10 [2022]. checked twice.
  c("hg38.Lareau.hg38_peaks", "hg38", "Regions of high homology to mtDNA (NUMT regions) defined by caleblareau/mitoblacklist", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/peaks/hg38_peaks.narrowPeak"), # checked -- Lareau-peaks, commit Apr 24, 2017. checked twice.
  c("hg19.Lareau.hg19_peaks", "hg19", "Regions of high homology to mtDNA (NUMT regions) defined by caleblareau/mitoblacklist", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/peaks/hg19_peaks.narrowPeak"), # checked -- Lareau-peaks, commit Apr 24, 2017. checked twice.
  
  ### MOUSE NUMT sets
  c("mm10.Lareau.mm10_peaks", "mm10", "Regions of high homology to mtDNA (NUMT regions) defined by caleblareau/mitoblacklist", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/peaks/mm10_peaks.narrowPeak"), # checked -- Lareau-peaks, mm10, commit Apr 24, 2017. checked twice
  c("mm9.Lareau.mm9_peaks", "mm9", "Regions of high homology to mtDNA (NUMT regions) defined by caleblareau/mitoblacklist", "2017", "GitHub", "https://github.com/caleblareau/mitoblacklist/raw/master/peaks/mm9_peaks.narrowPeak") # checked -- Lareau-peaks, mm9, commit Apr 24, 2017. checked twice.
)
# Make data frame
mtx <- as.data.frame(mtx)
# Assign column names
colnames(mtx) <- c("Name", "Assembly", "Description", "Year last updated", "Source", "ID/URL")

# Excluderanges
for (i in 1:nrow(mtx)) {
  # Process each record
  print(paste("Processing", i, mtx$Name[i]))
  # If source is ENCODE, construct URL
  if (mtx$Source[i] == "ENCODE") {
    URL <- paste0("https://www.encodeproject.org/files/", mtx$`ID/URL`[i], "/@@download/", mtx$`ID/URL`[i], ".bed.gz")
  } else {
    URL <- mtx$`ID/URL`[i]
  }

  # Process only non-publication data. Publication data requires manual processing
  if (!(mtx$Source[i] %in% c("Publication", "excluderanges"))) {
    # Download output file name
    if (mtx$Source[i] == "ENCODE") {
      # ENCODE IDs need to have ".bed.gz" extension
      fileNameOut1 <- file.path(dir_results, paste0(mtx$Name[i], ".bed.gz"))
    } else {
      # Directly downloadable files have their own extension
      fileNameOut1 <- file.path(dir_results, paste0(mtx$Name[i], ".", tools::file_ext(mtx$`ID/URL`[i])))
    }
    # BED output file name
    fileNameOut2 <- file.path(dir_results, "bed", paste0(mtx$Name[i], ".bed"))
    # RDS output file name
    fileNameOut3 <- file.path(dir_results, "rds", paste0(mtx$Name[i], ".rds"))
    
    # Download file, if doesn't exist
    if (!file.exists(fileNameOut1)) download.file(URL, fileNameOut1)
    
    # Special provisions for selected files
    if (mtx$Name[i] == "T2T.Lareau.chm13v2.0_peaks") {
      # T2T.Lareau.chm13v2.0_peaks has spaces instead of tabs
      # Convert spaces to tabs in place
      system(paste0("awk -v OFS=\"\\t\" '$1=$1' ", fileNameOut1, " > tmp.bed"))
      system(paste0("mv tmp.bed ", fileNameOut1))
    }
    if (mtx$Name[i] == "hg38.Yeo.eCLIP_Excludableregions.hg38liftover.bed.fixed") {
      # hg38.Yeo.eCLIP_Excludableregions.hg38liftover.bed.fixed has "-" instead of scores
      excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
      excludeBED[, 5] <- 0 # Replace scores by 0
      write_tsv(excludeBED, str_remove(fileNameOut1, ".gz$"), col_names = FALSE)
      system(paste0("gzip -f ", str_remove(fileNameOut1, ".gz$")))
    }
    if (mtx$Assembly[i] == "TAIR10") {
      # TAIR10 chromosomes are capitalized, like "Chr1". Convert to lower case
      excludeBED <- read_tsv(file = fileNameOut1, col_names = FALSE)
      # Sometimes chromosomes are sjut numbers, append "chr"
      if (sum(grepl("chr", excludeBED$X1, ignore.case = TRUE)) == 0) {
        excludeBED$X1 <- paste0("chr", excludeBED$X1)
      }
      # Convert to lower case and save
      excludeBED$X1 <- tolower(excludeBED$X1)
      write_tsv(excludeBED, file = fileNameOut1, col_names = FALSE)
    }
    
    # Have a look at the original data
    excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
    # Assign column names depending on the number of columns
    all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
    # Convert to GRanges object
    denyGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
    
    # Genome abbreviation
    genome_id <- mtx$Assembly[i]
    if (genome_id == "T2T") {
      # Seqinfo for T2T genome
      chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = "T2T-CHM13v2.0", assembled.molecules.only = TRUE)
      chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule)) # GCA_009914755.4
      # Make the same format as UCSC chromosome data
      chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                               size = chrom_data$SequenceLength,
                               assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                               circular = chrom_data$circular)
      # Keep standard chromosomes
      chromosomes_standard <- chrom_data$chrom
      # Subset and match to chromosomes in the denyGR object
      # Common chromosomes
      chromosomes_common <- intersect(chrom_data$chrom, seqlevels(denyGR))
      # Subset denyGR
      denyGR <- keepSeqlevels(denyGR, chromosomes_common, pruning.mode = "tidy")      
      # Subset chrom_data
      chrom_data <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
      # Match objects
      chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
      # Check if chromosome order is the same
      if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
        print(paste("Chromosome order does not match for", genome_id, "genome."))
        break
      }
      # Assign seqinfo data
      seqlengths(denyGR) <- chrom_data$size
      isCircular(denyGR) <- ifelse(is.na(chrom_data$circular), FALSE, TRUE)
      genome(denyGR)     <- "T2T-CHM13v2.0"
    } else if (genome_id == "TAIR10") {
      # Seqinfo for TAIR10 genome
      chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = "TAIR10", assembled.molecules.only = TRUE)
      chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule)) # GCA_000001735.1
      # Make the same format as UCSC chromosome data
      chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                               size = chrom_data$SequenceLength,
                               assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                               circular = chrom_data$circular)
      # Keep standard chromosomes
      chromosomes_standard <- chrom_data$chrom
      # Subset and match to chromosomes in the denyGR object
      # Common chromosomes
      chromosomes_common <- intersect(chrom_data$chrom, seqlevels(denyGR))
      # Subset denyGR
      denyGR <- keepSeqlevels(denyGR, chromosomes_common, pruning.mode = "tidy")      
      # Subset chrom_data
      chrom_data <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
      # Match objects
      chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
      # Check if chromosome order is the same
      if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
        print(paste("Chromosome order does not match for", genome_id, "genome."))
        break
      }
      # Assign seqinfo data
      seqlengths(denyGR) <- chrom_data$size
      isCircular(denyGR) <- ifelse(is.na(chrom_data$circular), FALSE, TRUE)
      genome(denyGR)     <- "TAIR10" # "GCA_000001735.1"
    } else {
      # Get chromosome info
      chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
      # Exclude chrMT, only in hg19 genome
      chrom_data <- chrom_data[chrom_data$chrom != "chrMT", ]
      # Keep standard chromosomes
      chromosomes_standard <- chrom_data$chrom
      # Subset and match to chromosomes in the denyGR object
      # Common chromosomes
      chromosomes_common <- intersect(chrom_data$chrom, seqlevels(denyGR))
      # Subset denyGR
      denyGR <- keepSeqlevels(denyGR, chromosomes_common, pruning.mode = "tidy")      
      # Subset chrom_data
      chrom_data <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
      # Match objects
      chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
      # Check if chromosome order is the same
      if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
        print(paste("Chromosome order does not match for", genome_id, "genome."))
        break
      }
      # Assign seqinfo data
      seqlengths(denyGR) <- chrom_data$size
      isCircular(denyGR) <- chrom_data$circular
      genome(denyGR)     <- genome_id
    }
    
    # Keep autosomes/sex/M chromosomes
    denyGR <- keepStandardChromosomes(x = denyGR, pruning.mode = "tidy")
    # Sort the object
    denyGR <- sort(denyGR)
    
    # Proceed only if non-empty. Note empty ones
    if (length(denyGR) == 0) {
      print(paste(mtx$Name[i], "is empty")) # Do nothing, just print
    } else {
      # Process with saving objects
      # Save BED file
      write_tsv(as.data.frame(denyGR), file = fileNameOut2, col_names = FALSE)
      # Save as Rds object
      saveRDS(object = denyGR, file = fileNameOut3)
      # excludeGR <- readRDS(file = fileNameOut3)
    }
  } else {
    print(paste("Skipped:", mtx$Name[i]))
  }
}

# Custom processing 
## danRer10.Yang.Supplemental_Table_19.ChIP-seq_black_list_in_the_zebrafish_genome
i <- which(mtx$Name == "danRer10.Yang.Supplemental_Table_19.ChIP-seq_black_list_in_the_zebrafish_genome")
# Download URL
URL <- mtx$`ID/URL`[i]
# Directly downloadable files have their own extension
fileNameOut1 <- file.path(dir_results, paste0(mtx$Name[i], ".", tools::file_ext(mtx$`ID/URL`[i])))
# BED output file name
fileNameOut2 <- file.path(dir_results, "bed", paste0(mtx$Name[i], ".bed"))
# RDS output file name
fileNameOut3 <- file.path(dir_results, "rds", paste0(mtx$Name[i], ".rds"))

# Download file, if doesn't exist
if (!file.exists(fileNameOut1)) {
  download.file(URL, fileNameOut1)
  # Unzip
  unzip(fileNameOut1, files = c("Supplemental_Tables/Supplemental Table 19. ChIP-seq black list in the zebrafish genome.xlsx"), junkpaths = TRUE, exdir = dir_results)
}
# Have a look inside
excludeBED <- read_xlsx(file.path(dir_results, "Supplemental Table 19. ChIP-seq black list in the zebrafish genome.xlsx"), col_names = FALSE)
# Save as tab-separated for rtracklayer
fileNameOut1 <- file.path(dir_results, paste0(mtx$Name[i], ".bed"))
write_tsv(excludeBED, fileNameOut1, col_names = FALSE)
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
denyGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)

# Genome abbreviation
genome_id <- mtx$Assembly[i]
# Get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the denyGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(denyGR))
# Subset denyGR
denyGR <- keepSeqlevels(denyGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(denyGR) <- chrom_data$size
isCircular(denyGR) <- chrom_data$circular
genome(denyGR)     <- genome_id

# Keep autosomes/sex/M chromosomes
denyGR <- keepStandardChromosomes(x = denyGR, pruning.mode = "tidy")
# Sort the object
denyGR <- sort(denyGR)
# Process with saving objects
# Save BED file
write_tsv(as.data.frame(denyGR), file = fileNameOut2, col_names = FALSE)
# Save as Rds object
saveRDS(object = denyGR, file = fileNameOut3)
# excludeGR <- readRDS(file = fileNameOut3)

## T2T.excluderanges 
i <- which(mtx$Name == "T2T.excluderanges")
URL <- mtx$`ID/URL`[i]
# Excluderanges need to have ".bed.gz" extension
fileNameOut1 <- file.path(dir_results, paste0(mtx$Name[i], ".bed.gz"))
# BED output file name
fileNameOut2 <- file.path(dir_results, "bed", paste0(mtx$Name[i], ".bed"))
# RDS output file name
fileNameOut3 <- file.path(dir_results, "rds", paste0(mtx$Name[i], ".rds"))
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
denyGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)

# Genome abbreviation
genome_id <- mtx$Assembly[i]
if (genome_id == "T2T") {
  # Seqinfo for T2T genome
  chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = "T2T-CHM13v2.0", assembled.molecules.only = TRUE) # GCA_009914755.4
  chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
  # Make the same format as UCSC chromosome data
  chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                           size = chrom_data$SequenceLength,
                           assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                           circular = chrom_data$circular)
  # Keep standard chromosomes
  chromosomes_standard <- chrom_data$chrom
  # Subset and match to chromosomes in the denyGR object
  # Common chromosomes
  chromosomes_common <- intersect(chrom_data$chrom, seqlevels(denyGR))
  # Subset denyGR
  denyGR <- keepSeqlevels(denyGR, chromosomes_common, pruning.mode = "tidy")      
  # Subset chrom_data
  chrom_data <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
  # Match objects
  chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  # Assign seqinfo data
  seqlengths(denyGR) <- chrom_data$size
  isCircular(denyGR) <- ifelse(is.na(chrom_data$circular), FALSE, TRUE)
  genome(denyGR)     <- "T2T-CHM13v2.0" # "GCA_009914755.4"
} 

# Keep autosomes/sex/M chromosomes
denyGR <- keepStandardChromosomes(x = denyGR, pruning.mode = "tidy")
# Sort the object
denyGR <- sort(denyGR)
# Process with saving objects
# Save BED file
write_tsv(as.data.frame(denyGR), file = fileNameOut2, col_names = FALSE)
# Save as Rds object
saveRDS(object = denyGR, file = fileNameOut3)
# excludeGR <- readRDS(file = fileNameOut3)

## mm39.excluderanges
i <- which(mtx$Name == "mm39.excluderanges")
URL <- mtx$`ID/URL`[i]
# Excluderanges need to have ".bed.gz" extension
fileNameOut1 <- file.path(dir_results, paste0(mtx$Name[i], ".bed.gz"))
# BED output file name
fileNameOut2 <- file.path(dir_results, "bed", paste0(mtx$Name[i], ".bed"))
# RDS output file name
fileNameOut3 <- file.path(dir_results, "rds", paste0(mtx$Name[i], ".rds"))
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
denyGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)

# Genome abbreviation
genome_id <- mtx$Assembly[i]
# Get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the denyGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(denyGR))
# Subset denyGR
denyGR <- keepSeqlevels(denyGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(denyGR) <- chrom_data$size
isCircular(denyGR) <- chrom_data$circular
genome(denyGR)     <- genome_id

# Keep autosomes/sex/M chromosomes
denyGR <- keepStandardChromosomes(x = denyGR, pruning.mode = "tidy")
# Sort the object
denyGR <- sort(denyGR)
# Process with saving objects
# Save BED file
write_tsv(as.data.frame(denyGR), file = fileNameOut2, col_names = FALSE)
# Save as Rds object
saveRDS(object = denyGR, file = fileNameOut3)
# excludeGR <- readRDS(file = fileNameOut3)

# 02b_preprocessing.Rmd
# Creates gap excluderanges (problematic genomic regions from the UCSC database)

# Download bigBedtoBed utility
# Tested on Mac, needs to be changed depending on OS.
if (!file.exists(file.path(dir_results, 'bigBedtoBed'))) {
  # download bigBedtoBed converter for mac -- this needs to be adjusted to work with other OS'.
  # location of macOS bigBedtoBed
  bigBedtoBed.URL <- 'http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed'
  # out file name
  bigBedtoBed.out <- file.path(dir_results, 'bigBedtoBed')
  # download bigBedtoBed
  download.file(bigBedtoBed.URL, bigBedtoBed.out)
  # adjust permissions
  system(paste0('chmod 755 ',bigBedtoBed.out))
} else {
  bigBedtoBed.out <- file.path(dir_results, 'bigBedtoBed')
}

# Create a table
# Manually created matrix with information for download
mtx <- rbind(
  ### HUMAN NUMTs
  c("hg19.UCSC.numtS", "hg19", "Human NumtS mitochondrial sequence", "2011", "UCSC", "numtS"),
  ### MOUSE NUMTs
  c("mm9.UCSC.numtS", "mm9", "Mouse NumtS mitochondrial sequence", "2011", "UCSC", "numtS"),
  
  ### HUMAN gap sets
  ## T2T
  c("T2T.CHM13.chm13.draft_v2.0.cen_mask", "T2T", "Centromeric satellite masking bed file (v2.0)", "2022", "CHM13", "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.cen_mask.bed"),
  c("T2T.CHM13.chm13.draft_v1.1.telomere", "T2T", "Telomere identified by the VGP pipeline (v1.1)", "2022", "CHM13", "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v1.1.telomere.bed.gz"),
  c("T2T.UCSC.censat", "T2T", "T2T peri/centromeric satellite annotation (v2.0, 20220329, CHM13 v2.0)", "2022", "UCSChub", "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.censat/censat.bb"),
  c("T2T.UCSC.gap", "T2T", "Locations of assembly gaps, as determine by strings of 'N' characters (v1.0)", "2021", "UCSChub", "http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/gap.bigBed"),
  c("T2T.UCSC.hgUnique.hg38", "T2T", "Regions unique to the T2T-CHM13 v2.0 assembly compared to the GRCh38/hg38 and GRCh37/hg19 reference assemblies", "2022", "UCSChub", "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.hgUnique/hgUnique.hg38.bb"), # checked --T2T unique regions, updated 2022-04-09. 
  ## hg38
  c("hg38.UCSC.centromere", "hg38", "Gaps from centromeres", "2014", "UCSC", "centromeres"), # checked, hg38 centromeres, updated 2014-01-09
  c("hg38.UCSC.telomere", "hg38", "Gaps from telomeres", "2018", "UCSC", "gap"), # checked -- hg38 gaps table, updated 2018-08-06. Unsure of best link for source.
  c("hg38.UCSC.short_arm", "hg38", "Gaps on the short arm of the chromosome", "2018", "UCSC", "gap"), # checked -- hg38 gaps table, updated 2018-08-06. Unsure of best link for source.
  c("hg38.UCSC.heterochromatin", "hg38", "Gaps from large blocks of heterochromatin", "2018", "UCSC", "gap"), # checked -- hg38 gaps table, updated 2018-08-06. Unsure of best link for source.
  c("hg38.UCSC.contig", "hg38", "Gaps between contigs in scaffolds", "2018", "UCSC", "gap"), # checked -- hg38 gaps table, updated 2018-08-06. Unsure of best link for source.
  c("hg38.UCSC.scaffold", "hg38", "Gaps between scaffolds in chromosome assemblies. Has extra chromosomes, use keepStandardChromosomes() filtering", "2018", "UCSC", "gap"), # checked -- hg38 gaps table, updated 2018-08-06. Unsure of best link for source.
  ## hg19
  c("hg19.UCSC.centromere", "hg19", "Gaps from centromeres", "2020", "UCSC", "gap"), # checked -- hg19 gaps table, updated 2020-02-20. Unsure of best link for source. 
  c("hg19.UCSC.telomere", "hg19", "Gaps from telomeres", "2020", "UCSC", "gap"), # checked -- hg19 gaps table, updated 2020-02-20. Unsure of best link for source.
  c("hg19.UCSC.short_arm", "hg19", "Gaps on the short arm of the chromosome", "2020", "UCSC", "gap"), # checked -- hg19 gaps table, updated 2020-02-20. Unsure of best link for source. 
  c("hg19.UCSC.heterochromatin", "hg19", "Gaps from large blocks of heterochromatin", "2020", "UCSC", "gap"), # checked -- hg19 gaps table, updated 2020-02-20. Unsure of best link for source. 
  c("hg19.UCSC.clone", "hg19", "Gaps between clones in the same map contig. Has extra chromosomes, use keepStandardChromosomes() filtering", "2020", "UCSC", "gap"), # checked -- hg19 gaps table, updated 2020-02-20. Unsure of best link for source. 
  c("hg19.UCSC.contig", "hg19", "Gaps between contigs in scaffolds", "2020", "UCSC", "gap"), # checked -- hg19 gaps table, updated 2020-02-20. Unsure of best link for source. 
  c("hg19.UCSC.scaffold", "hg19", "Gaps between scaffolds in chromosome assemblies. Only non-autosomal chromosomes", "2020", "UCSC", "gap"), # checked -- hg19 scaffold NOT described on table schema, but IS in manual and automated gap table download.
  
  ### MOUSE gap sets
  ## mm39
  c("mm39.UCSC.centromere", "mm39", "Gaps from centromeres", "2020", "UCSC", "gap"), # checked -- mm39 gaps tables, updated 2020-07-27. Unsure of best link for source.
  c("mm39.UCSC.telomere", "mm39", "Gaps from telomeres", "2020", "UCSC", "gap"), # checked -- mm39 gaps tables, updated 2020-07-27. Unsure of best link for source.
  c("mm39.UCSC.short_arm", "mm39", "Gaps on the short arm of the chromosome" ,"2020", "UCSC", "gap"), # checked -- mm39 gaps tables, updated 2020-07-27. Unsure of best link for source.
  c("mm39.UCSC.contig", "mm39", "Gaps between contigs in scaffolds", "2020", "UCSC", "gap"), # checked -- mm39 gaps tables, updated 2020-07-27. Unsure of best link for source.
  c("mm39.UCSC.scaffold", "mm39", "Gaps between scaffolds in chromosome assemblies", "2020", "UCSC", "gap"), # checked -- mm39 gaps tables, updated 2020-07-27. Unsure of best link for source.
  ## mm10
  c("mm10.UCSC.centromere", "mm10", "Gaps from centromeres", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  c("mm10.UCSC.telomere",   "mm10", "Gaps from telomeres", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  c("mm10.UCSC.short_arm", "mm10", "Gaps on the short arm of the chromosome", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  c("mm10.UCSC.clone", "mm10", "Gaps between clones in the same map contig. Has extra chromosomes, use keepStandardChromosomes() filtering", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  c("mm10.UCSC.contig", "mm10", "Gaps between contigs in scaffolds", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  c("mm10.UCSC.scaffold", "mm10", "Gaps between scaffolds in chromosome assemblies", "2021", "UCSC", "gap"), # checked -- mm10 scaffold NOT described on table schema, but IS in manual and automated gap table download.
  c("mm10.UCSC.other", "mm10", "Sequence of Ns in the assembly that were not marked as gaps in the AGP (A Golden Path) assembly definition file. Has extra chromosomes, use keepStandardChromosomes() filtering", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  c("mm10.UCSC.fragment", "mm10", "A single gap of 31 bases in chrX_GL456233_random", "2021", "UCSC", "gap"), # checked -- mm10 gaps table, updated 2021-04-08. Unsure of best link for source.
  ## mm9
  c("mm9.UCSC.centromere", "mm9", "Gaps from centromeres", "2007", "UCSC", "gap"), # checked -- mm9 gaps table, updated 2007-07-23. Unsure of best link for source.
  c("mm9.UCSC.fragment", "mm9", "Gaps between the contigs of a draft clone. (In this context, a contig is a set of overlapping sequence reads). Has extra chromosomes, use keepStandardChromosomes() filtering", "2007", "UCSC", "gap"), # checked -- mm9 gaps table, updated 2007-07-23. Unsure of best link for source.
  c("mm9.UCSC.contig", "mm9", "Gaps between contigs in scaffolds. Has extra chromosomes, use keepStandardChromosomes() filtering", "2007", "UCSC", "gap"), # checked -- mm9 gaps table, updated 2007-07-23. Unsure of best link for source.
  
  ### ZEBRAFISH gap sets
  c("danRer10.UCSC.contig", "danRer10", "Gaps between contigs in scaffolds", "2015", "UCSC", "gap"), # checked -- danRer10 gap table, updated 2015-01-22
  c("danRer10.UCSC.scaffold", "danRer10", "Gaps between scaffolds in chromosome assemblies", "2015", "UCSC", "gap"), # checked -- danRer10 gap table, updated 2015-01-22
  
  ### FLY gap sets
  c("dm6.UCSC.other", "dm6", "Sequence of Ns in the assembly that were not marked as gaps in the AGP (A Golden Path) assembly definition file", "2014", "UCSC", "gap"), # checked -- dm6 gap table, updated 2014-08-28
  c("dm3.UCSC.contig", "dm3", "Gaps between contigs in scaffolds", "2006", "UCSC", "gap"), # checked -- dm3 gap table, updated 2006-07-10
  c("dm3.UCSC.scaffold", "dm3", "Gaps between scaffolds in chromosome assemblies", "2006", "UCSC", "gap"), # checked -- dm3 gap table, updated 2006-07-10
  
  ### Arabidopsis
  c("TAIR10.UCSC.araTha1.gap", "TAIR10", "Gaps in the May 2011 Arabidopsis thaliana genome assembly", "2013", "UCSChub", "https://genome-test.gi.ucsc.edu/~hiram/hubs/Plants/araTha1/bbi/araTha1.gap.bb")
)
# Make data frame
mtx <- as.data.frame(mtx)
# Assign column names
colnames(mtx) <- c("Name", "Assembly", "Description", "Year last updated", "Source", "ID/URL")

# Process each genome
for (genome_id in unique(mtx$Assembly) ) {
  # T2T and TAIR10 processed later
  if (genome_id %in% (c("T2T", "TAIR10"))) next
  # Get chromosome info
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
  # Keep standard chromosomes
  chromosomes_standard <- chrom_data$chrom
  # Get genome-specific gaps table
  mySession <- browserSession()
  genome(mySession) <- genome_id
  # gaps <- getTable(ucscTableQuery(mySession, track = "gap"))
  query <- ucscTableQuery(mySession, table = "gap")
  gaps <- getTable(query)
  # Look inside Original data
  print(head(gaps))
  # Summary of gap types
  gap_types <- table(gaps$type)
  print(gap_types)
  # Process each gap type
  for (gap_type in names(gap_types)) {
    gaps_selected <- gaps[gaps$type == gap_type, ]
    # UCSC coordinates are 0 based: 
    # "chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    #             ... For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100"
    gapsGR <- makeGRangesFromDataFrame(gaps_selected, keep.extra.columns = TRUE) # starts.in.df.are.0based=T, 
    # Sort GR object
    gapsGR <- sort(gapsGR)
    # Select seqinfo data for the gaps object
    # Subset and match to chromosomes in the gapsGR object
    # Common chromosomes
    chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
    # Subset gapsGR
    gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
    # Subset chrom_data
    chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
    # Match objects
    chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
    # Check if chromosome order is the same
    if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
      print(paste("Chromosome order does not match for", genome_id, "genome."))
      break
    }
    
    # Assign seqinfo data
    seqlengths(gapsGR) <- chrom_data_subset$size
    isCircular(gapsGR) <- chrom_data_subset$circular
    genome(gapsGR)     <- genome_id
        # Construct file name, e.g., hg19.UCSC.gap_centromere.bed
        fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".rds")
        # Save as Rds object
        saveRDS(object = gapsGR, file = file.path(dir_results, "rds", fileNameOut))
        # Save as bed file
        fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".bed")
        write_tsv(as.data.frame(gapsGR), file = file.path(dir_results, "bed", fileNameOut), col_names = FALSE)
  }
}

## NUMTs for hg19 and mm9
for (genome_id in c("hg19", "mm9")) {
  # Get chromosome info
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
  # Keep standard chromosomes
  chromosomes_standard <- chrom_data$chrom
  # Get numtS
  mySession <- browserSession()
  genome(mySession) <- genome_id
  query <- ucscTableQuery(mySession, table = "numtS")
  gaps <- getTable(query)
  # Look inside Original data
  print(head(gaps))
  # Process
  gap_type <- 'numtS'
  # mm9 has additional columns throwing an error when creating GRanges
  gapsGR <- makeGRangesFromDataFrame(gaps[, 1:7], keep.extra.columns = TRUE) # starts.in.df.are.0based=T, 
  # Sort GR object
  gapsGR <- sort(gapsGR)
  # Look inside 
  print(gapsGR)
  # Select seqinfo data for the centromeres object
  # Subset and match to chromosomes in the gapsGR object
  # Common chromosomes
  chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
  # Subset gapsGR
  gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
  # Subset chrom_data
  chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
  # Match objects
  chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  
  # Assign seqinfo data
  seqlengths(gapsGR) <- chrom_data_subset$size
  isCircular(gapsGR) <- chrom_data_subset$circular
  genome(gapsGR)     <- genome_id
    # Construct file name, e.g., hg19.UCSC.gap_centromere.rds
    fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".rds")
    # Save as Rds object
    saveRDS(object = gapsGR, file = file.path(dir_results, "rds", fileNameOut))
    # save as bed file
    fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".bed")
    write_tsv(as.data.frame(gapsGR), file = file.path(dir_results, "bed", fileNameOut), col_names = FALSE)
}

## Centromeres for hg38
# Get hg38 centromeres track, not part of gap table
genome_id <- 'hg38'
chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id, assembled.molecules.only = TRUE)
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Get hg38 centromere gaps
mySession <- browserSession()
genome(mySession) <- genome_id
query <- ucscTableQuery(mySession, table = "centromeres")
gaps <- getTable(query)
# Process
gap_type <- 'centromere'
gapsGR <- makeGRangesFromDataFrame(gaps, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# Look inside 
print(gapsGR)
# Select seqinfo data for the centromeres object
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}

# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- chrom_data_subset$circular
genome(gapsGR)     <- genome_id
  # Construct file name, e.g., hg19.UCSC.gap_centromere.rds
  fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".rds")
  # Save as Rds object
  saveRDS(object = gapsGR, file = file.path(dir_results, "rds", fileNameOut))
  # save as bed file
  fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".bed")
  write_tsv(as.data.frame(gapsGR), file = file.path(dir_results, "bed", fileNameOut), col_names = FALSE)

# T2T centromere
### Get T2T-unique track
genome_id <- "T2T-CHM13v2.0"
# T2T is NOT registered in GenomeInfoDb, see GenomeInfoDb::registered_UCSC_genomes()
# location of T2T unique regions
fileNameURL <- "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.cen_mask.bed"
# out file names
# BED
fileNameOut1 <- file.path(dir_results,'chm13.draft_v2.0.cen_mask.bed')
# bed
fileNameOut2 <- file.path(dir_results, "bed", 'T2T.CHM13.chm13.draft_v2.0.cen_mask.bed')
# rds
fileNameOut3 <- file.path(dir_results, "rds", 'T2T.CHM13.chm13.draft_v2.0.cen_mask.rds')
# download T2T unique regions
if (!file.exists(fileNameOut1)) download.file(fileNameURL, fileNameOut1)
# convert T2T unique to bed
# system(paste0(bigBedtoBed.out,' ',fileNameOut1,' ',fileNameOut2))
# read in T2T unique
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
gapsGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = genome_id, assembled.molecules.only = TRUE)
chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
# Make the same format as UCSC chromosome data
chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                         size = chrom_data$SequenceLength,
                         assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                         circular = chrom_data$circular)
# Rename chrMT to chrM, for T2T assembly
chrom_data$chrom[chrom_data$chrom == "chrMT"] <- "chrM"
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- ifelse(is.na(chrom_data_subset$circular), FALSE, TRUE)
genome(gapsGR)     <- genome_id
  # save as .rds
  saveRDS(object = gapsGR, file = fileNameOut3)
  # save as .bed
  write_tsv(as.data.frame(gapsGR), file = fileNameOut2, col_names = FALSE)

# T2T telomere
### Get T2T-unique track
genome_id <- "T2T-CHM13v1.1"
# T2T is NOT registered in GenomeInfoDb, see GenomeInfoDb::registered_UCSC_genomes()
# location of T2T unique regions
fileNameURL <- 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v1.1.telomere.bed.gz'
# out file names
# BED
fileNameOut1 <- file.path(dir_results,'chm13.draft_v1.1.telomere.bed.gz')
# bed
fileNameOut2 <- file.path(dir_results, "bed", 'T2T.CHM13.chm13.draft_v1.1.telomere.bed')
# rds
fileNameOut3 <- file.path(dir_results, "rds", 'T2T.CHM13.chm13.draft_v1.1.telomere.rds')
# download T2T unique regions
if (!file.exists(fileNameOut1)) download.file(fileNameURL, fileNameOut1)
# convert T2T unique to bed
# system(paste0(bigBedtoBed.out,' ',fileNameOut1,' ',fileNameOut2))
# read in T2T unique
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut1, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
gapsGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# Look inside 
print(gapsGR)

# get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = genome_id, assembled.molecules.only = TRUE)
chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
# Make the same format as UCSC chromosome data
chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                         size = chrom_data$SequenceLength,
                         assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                         circular = chrom_data$circular)
# Rename chrMT to chrM, for T2T assembly
chrom_data$chrom[chrom_data$chrom == "chrMT"] <- "chrM"
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- ifelse(is.na(chrom_data_subset$circular), FALSE, TRUE)
genome(gapsGR)     <- genome_id
  # save as .rds
  saveRDS(object = gapsGR, file = fileNameOut3)
  # save as .bed
  write_tsv(as.data.frame(gapsGR), file = fileNameOut2, col_names = FALSE)

# T2T censat
### Get T2T-unique track
genome_id <- "T2T-CHM13v2.0"
# T2T is NOT registered in GenomeInfoDb, see GenomeInfoDb::registered_UCSC_genomes()
# location of T2T unique regions
fileNameURL <- 'https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.censat/censat.bb'
# out file names
# BED
fileNameOut1 <- file.path(dir_results,'censat.bb')
# bed
fileNameOut2 <- file.path(dir_results, "bed", 'T2T.UCSC.censat.bed')
# rds
fileNameOut3 <- file.path(dir_results, "rds", 'T2T.UCSC.censat.rds')
# download T2T unique regions
if (!file.exists(fileNameOut1)) download.file(fileNameURL, fileNameOut1)
# convert T2T unique to bed
system(paste0(bigBedtoBed.out,' ',fileNameOut1,' ',fileNameOut2))
# read in T2T unique
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut2, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
gapsGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# Add censat type
gapsGR$type <- sapply(gapsGR$name, function(x) strsplit(x, "_")[[1]][1])
# get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = genome_id, assembled.molecules.only = TRUE)
chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
# Make the same format as UCSC chromosome data
chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                         size = chrom_data$SequenceLength,
                         assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                         circular = chrom_data$circular)
# Rename chrMT to chrM, for T2T assembly
chrom_data$chrom[chrom_data$chrom == "chrMT"] <- "chrM"
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- ifelse(is.na(chrom_data_subset$circular), FALSE, TRUE)
genome(gapsGR)     <- genome_id
  # save as .rds
  saveRDS(object = gapsGR, file = fileNameOut3)
  # save as .bed
  write_tsv(as.data.frame(gapsGR), file = fileNameOut2, col_names = FALSE)

# T2T gap
### Get T2T-unique track
genome_id <- "T2T-CHM13v1.0"
# T2T is NOT registered in GenomeInfoDb, see GenomeInfoDb::registered_UCSC_genomes()
# location of T2T unique regions
fileNameURL <- 'http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/gap.bigBed'
# out file names
# BED
fileNameOut1 <- file.path(dir_results,'gap.bigBed')
# bed
fileNameOut2 <- file.path(dir_results, "bed", 'T2T.UCSC.gap.bed')
# rds
fileNameOut3 <- file.path(dir_results, "rds", 'T2T.UCSC.gap.rds')
# download T2T unique regions
if (!file.exists(fileNameOut1)) download.file(fileNameURL, fileNameOut1)
# convert T2T unique to bed
system(paste0(bigBedtoBed.out,' ',fileNameOut1,' ',fileNameOut2))
# read in T2T unique
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut2, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
gapsGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = genome_id, assembled.molecules.only = TRUE)
chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
# Make the same format as UCSC chromosome data
chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                         size = chrom_data$SequenceLength,
                         assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                         circular = chrom_data$circular)
# Rename chrMT to chrM, for T2T assembly
chrom_data$chrom[chrom_data$chrom == "chrMT"] <- "chrM"
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- ifelse(is.na(chrom_data_subset$circular), FALSE, TRUE)
genome(gapsGR)     <- genome_id
  # save as .rds
  saveRDS(object = gapsGR, file = fileNameOut3)
  # save as .bed
  write_tsv(as.data.frame(gapsGR), file = fileNameOut2, col_names = FALSE)

# T2T hgUniqueHg38
### Get T2T-unique track
genome_id <- "T2T-CHM13v2.0"
# T2T is NOT registered in GenomeInfoDb, see GenomeInfoDb::registered_UCSC_genomes()
# location of T2T unique regions
fileNameURL <- 'https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.hgUnique/hgUnique.hg38.bb'
# out file names
# bigbed
fileNameOut1 <- file.path(dir_results,'T2T.UCSC.hgUnique.hg38.bb')
# bed
fileNameOut2 <- file.path(dir_results, "bed", 'T2T.UCSC.hgUnique.hg38.bed')
# rds
fileNameOut3 <- file.path(dir_results, "rds", 'T2T.UCSC.hgUnique.hg38.rds')
# download T2T unique regions
if (!file.exists(fileNameOut1)) download.file(fileNameURL, fileNameOut1)
# convert T2T unique to bed
system(paste0(bigBedtoBed.out,' ',fileNameOut1,' ',fileNameOut2))
# read in T2T unique
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut2, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
gapsGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = genome_id, assembled.molecules.only = TRUE)
chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
# Make the same format as UCSC chromosome data
chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                         size = chrom_data$SequenceLength,
                         assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                         circular = chrom_data$circular)
# Rename chrMT to chrM, for T2T assembly
chrom_data$chrom[chrom_data$chrom == "chrMT"] <- "chrM"
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- ifelse(is.na(chrom_data_subset$circular), FALSE, TRUE)
genome(gapsGR)     <- genome_id
  # save as .rds
  saveRDS(object = gapsGR, file = fileNameOut3)
  # save as .bed
  write_tsv(as.data.frame(gapsGR), file = fileNameOut2, col_names = FALSE)

# TAIR10
### Get T2T-unique track
genome_id <- "TAIR10" # "GCA_000001735.1"
# location of araTha1 gap regions
fileNameURL <- 'https://genome-test.gi.ucsc.edu/~hiram/hubs/Plants/araTha1/bbi/araTha1.gap.bb'
# out file names
# BED
fileNameOut1 <- file.path(dir_results,'araTha1.gap.bb')
# bed
fileNameOut2 <- file.path(dir_results,"bed", 'TAIR10.UCSC.araTha1.gap.bed')
# rds
fileNameOut3 <- file.path(dir_results,"rds", 'TAIR10.UCSC.araTha1.gap.rds')
# download araTha1 gap regions
if (!file.exists(fileNameOut1)) download.file(fileNameURL, fileNameOut1)
# convert araTha1 gap to bed
system(paste0(bigBedtoBed.out,' ',fileNameOut1,' ',fileNameOut2))
# read in araTha1 gap
# Have a look at the original data
excludeBED <- read_tsv(fileNameOut2, col_names = FALSE)
# Assign column names depending on the number of columns
all_columns <- c("chr", "start", "stop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
colnames(excludeBED) <- all_columns[1:ncol(excludeBED)]
# Convert to GRanges object
gapsGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
# Sort GR object
gapsGR <- sort(gapsGR)
# get chromosome info
chrom_data <- GenomeInfoDb::getChromInfoFromNCBI(assembly = genome_id, assembled.molecules.only = TRUE)
chrom_data$AssignedMolecule <- as.character(paste0("chr", chrom_data$AssignedMolecule))
# Make the same format as UCSC chromosome data
chrom_data <- data.frame(chrom = chrom_data$AssignedMolecule,
                         size = chrom_data$SequenceLength,
                         assembled = ifelse(chrom_data$AssemblyUnit == "Primary Assembly", TRUE, FALSE),
                         circular = chrom_data$circular)
# Keep standard chromosomes
chromosomes_standard <- chrom_data$chrom
# Subset and match to chromosomes in the gapsGR object
# Common chromosomes
chromosomes_common <- intersect(chrom_data$chrom, seqlevels(gapsGR))
# Subset gapsGR
gapsGR <- keepSeqlevels(gapsGR, chromosomes_common, pruning.mode = "tidy")      
# Subset chrom_data
chrom_data_subset <- chrom_data[chrom_data$chrom %in% chromosomes_common, ]
# Match objects
chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
# Check if chromosome order is the same
if (!all.equal(seqlevels(gapsGR), chrom_data_subset$chrom)) {
  print(paste("Chromosome order does not match for", genome_id, "genome."))
  break
}
# Assign seqinfo data
seqlengths(gapsGR) <- chrom_data_subset$size
isCircular(gapsGR) <- ifelse(is.na(chrom_data_subset$circular), FALSE, TRUE)
genome(gapsGR)     <- genome_id
  # save as .rds
  saveRDS(object = gapsGR, file = fileNameOut3)
  # save as .bed
  write_tsv(as.data.frame(gapsGR), file = fileNameOut2, col_names = FALSE)



# FIRST VERSION: MANUAL DOWNLOAD, LIMITED OBJECTS, HUMAN AND MOUSE ORGANISMS.
# INCLUDED NON-PRIMARY ASSEMBLY CHROMOSOMES. SUPERSEDED BY THE SECOND VERSION. 
# 
# excluderanges package currently contains 19 Rds objects. The original and 
# processed data are available at 
# https://drive.google.com/drive/folders/124DZtsU0YVWqkb7dgu8Nk6b3N8-ShVSC?usp=sharing

# The object names are structured as "<genome assembly>.<lab>.<original file name>", 
# e.g., "hg19.Birney.wgEncodeDacMapabilityConsensusExcludable".

# # ENCODE
# ## hg19
# 
# # Mint_Blacklist_hg19.bed.gz
# # https://www.encodeproject.org/files/ENCFF200UUD/, Bradley Bernstein, Broad
# wget https://www.encodeproject.org/files/ENCFF200UUD/@@download/ENCFF200UUD.bed.gz
# gzip --decompress --stdout ENCFF200UUD.bed.gz | bedtools sort -i - > hg19.Bernstein.Mint_Blacklist_hg19.bed
# 
# # wgEncodeDacMapabilityConsensusExcludable.bed.gz, Ewan Birney, EBI
# # https://www.encodeproject.org/files/ENCFF001TDO/, Ewan Birney, EBI
# wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
# gzip --decompress --stdout ENCFF001TDO.bed.gz | bedtools sort -i - > hg19.Birney.wgEncodeDacMapabilityConsensusExcludable.bed
# 
# # wgEncodeDukeMapabilityRegionsExcludable.bed.gz
# # https://www.encodeproject.org/files/ENCFF001THR/, Gregory Crawford, Duke
# wget https://www.encodeproject.org/files/ENCFF001THR/@@download/ENCFF001THR.bed.gz
# gzip --decompress --stdout ENCFF001THR.bed.gz | bedtools sort -i - > hg19.Crawford.wgEncodeDukeMapabilityRegionsExcludable.bed
# 
# # hg19mitoblack.bed.gz
# # https://www.encodeproject.org/files/ENCFF055QTV/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF055QTV/@@download/ENCFF055QTV.bed.gz
# gzip --decompress --stdout ENCFF055QTV.bed.gz | bedtools sort -i - > hg19.Wold.hg19mitoblack.bed
# 
# # eCLIP_blacklistregions.hg19.bed.gz
# # https://www.encodeproject.org/files/ENCFF039QTN/, Gene Yeo, UCSD
# wget https://www.encodeproject.org/files/ENCFF039QTN/@@download/ENCFF039QTN.bed.gz
# gzip --decompress --stdout ENCFF039QTN.bed.gz | bedtools sort -i - > hg19.Yeo.eCLIP_blacklistregions.hg19.bed
# 
# ## hg38
# 
# # Mint_Blacklist_GRCh38.bed.gz
# # https://www.encodeproject.org/files/ENCFF023CZC/, Bradley Bernstein, Broad
# wget https://www.encodeproject.org/files/ENCFF023CZC/@@download/ENCFF023CZC.bed.gz
# gzip --decompress --stdout ENCFF023CZC.bed.gz | bedtools sort -i - > hg38.Bernstein.Mint_Blacklist_GRCh38.bed
# 
# # GRCh38_unified_blacklist.bed.gz
# # https://www.encodeproject.org/files/ENCFF356LFX/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
# gzip --decompress --stdout ENCFF356LFX.bed.gz | bedtools sort -i - > hg38.Kundaje.GRCh38_unified_blacklist.bed
# 
# # GRCh38.blacklist.bed.gz
# # https://www.encodeproject.org/files/ENCFF419RSJ/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz
# gzip --decompress --stdout ENCFF419RSJ.bed.gz | bedtools sort -i - > hg38.Kundaje.GRCh38.blacklist.bed
# 
# # wgEncodeDacMapabilityConsensusExcludable.hg38.bed.gz
# # https://www.encodeproject.org/files/ENCFF220FIN/, Tim Reddy, Duke
# wget https://www.encodeproject.org/files/ENCFF220FIN/@@download/ENCFF220FIN.bed.gz
# gzip --decompress --stdout ENCFF220FIN.bed.gz | bedtools sort -i - > hg38.Reddy.wgEncodeDacMapabilityConsensusExcludable.hg38.bed
# 
# # hg38mitoblack.bed.gz
# # https://www.encodeproject.org/files/ENCFF940NTE/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF940NTE/@@download/ENCFF940NTE.bed.gz
# gzip --decompress --stdout ENCFF940NTE.bed.gz | bedtools sort -i - > hg38.Wold.hg38mitoblack.bed
# 
# # eCLIP_blacklistregions.hg38liftover.bed.fixed.bed.gz
# # https://www.encodeproject.org/files/ENCFF269URO/, Gene Yeo, UCSD
# wget https://www.encodeproject.org/files/ENCFF269URO/@@download/ENCFF269URO.bed.gz
# gzip --decompress --stdout ENCFF269URO.bed.gz | bedtools sort -i - > hg38.Yeo.eCLIP_blacklistregions.hg38liftover.bed.fixed.bed
# 
# ## mm10
# 
# # blacklist.full.bed.gz
# # https://www.encodeproject.org/files/ENCFF790DJT/, Ross Hardison, PennState
# wget https://www.encodeproject.org/files/ENCFF790DJT/@@download/ENCFF790DJT.bed.gz
# gzip --decompress --stdout ENCFF790DJT.bed.gz | bedtools sort -i - > mm10.Hardison.blacklist.full.bed
# 
# # psublacklist.mm10.bed.gz
# # https://www.encodeproject.org/files/ENCFF226BDM/, Ross Hardison, PennState
# wget https://www.encodeproject.org/files/ENCFF226BDM/@@download/ENCFF226BDM.bed.gz
# gzip --decompress --stdout ENCFF226BDM.bed.gz | bedtools sort -i - > mm10.Hardison.psublacklist.mm10.bed
# 
# # anshul.blacklist.mm10.bed.gz
# # https://www.encodeproject.org/files/ENCFF999QPV/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF999QPV/@@download/ENCFF999QPV.bed.gz
# gzip --decompress --stdout ENCFF999QPV.bed.gz | bedtools sort -i - > mm10.Kundaje.anshul.blacklist.mm10.bed
# 
# # mm10.blacklist.bed.gz
# # https://www.encodeproject.org/files/ENCFF547MET/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz
# gzip --decompress --stdout ENCFF547MET.bed.gz | bedtools sort -i - > mm10.Kundaje.mm10.blacklist.bed
# 
# # mm10mitoblack.bed.gz
# # https://www.encodeproject.org/files/ENCFF759PJK/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF759PJK/@@download/ENCFF759PJK.bed.gz
# gzip --decompress --stdout ENCFF759PJK.bed.gz | bedtools sort -i - > mm10.Wold.mm10mitoblack.bed
# 
# ## mm9 
# 
# # mm9mitoblack.bed.gz
# # https://www.encodeproject.org/files/ENCFF299EZH/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF299EZH/@@download/ENCFF299EZH.bed.gz
# gzip --decompress --stdout ENCFF299EZH.bed.gz | bedtools sort -i - > mm9.Wold.mm9mitoblack.bed
# 
# # http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/
# 
# # http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz
# wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz
# gzip --decompress --stdout dm3-blacklist.bed.gz | bedtools sort -i - > dm3.Kundaje.dm3-blacklist.bed
# 
# # http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz
# wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz
# gzip --decompress --stdout ce10-blacklist.bed.gz | bedtools sort -i - > ce10.Kundaje.ce10-blacklist.bed

# The following example demonstrate how the coordinates of assembly-specific 
# problematic regions were converted into Rds objects

# Download results data from https://drive.google.com/drive/folders/124DZtsU0YVWqkb7dgu8Nk6b3N8-ShVSC?usp=sharing
# library(GenomicRanges)
# library(GenomeInfoDb)
# library(rtracklayer)
# library(readr)
# # Folder with results
# dir_in <- "/Users/mdozmorov/Documents/Data/GoogleDrive/excludedata"
# 
# # All BED files
# files <- list.files(path = dir_in, pattern = "bed$")
# # In each subfolder
# for (file in files) {
#   # Read "fimo.bed" created by "fimo.qsub"
#   excludeBED <- read_tsv(file.path(dir_in, file), col_names = FALSE)
#   # Assign column names depending on the number of columns
#   if (ncol(excludeBED) == 3) { # Only 3 columns
#     colnames(excludeBED) <- c("chr", "start", "stop")
#   } else { # If more than 3 columns, consider the first 6
#     colnames(excludeBED) <- c("chr", "start", "stop", "name", "score", "strand")
#   }
#   # Convert to GRanges object
#   denyGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)
#   # Add seqinfo
#   # Parse out genome ID from the file name, to get hg19, hg38, mm9, mm10, etc.
#   genome_id <- strsplit(file, ".", fixed = TRUE)[[1]][1]
#   # Get chromosome info and match it to the chromosome order in excludeBED
#   chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
#   chrom_data <- chrom_data[chrom_data$chrom %in% seqlevels(denyGR), ]
#   chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
#   # Check if chromosome order is the same
#   if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
#     print(paste("Chromosome order does not match for", genome_id, "genome."))
#     break
#   }
#   # Assign seqinfo data
#   seqlengths(denyGR) <- chrom_data$size
#   isCircular(denyGR) <- chrom_data$circular
#   genome(denyGR)     <- genome_id
#   # Reformat output file name
#   fileNameOut <- sub("blacklist", "Excludable", file, ignore.case = TRUE)
#   fileNameOut <- sub("black", "Excludable", fileNameOut, ignore.case = TRUE)
#   
#   # Save as Rds object. file extension is changed to rds
#   saveRDS(object = denyGR, file = file.path(dir_in, sub("bed$", "rds", fileNameOut)))
#   # excludeGR <- readRDS(file = file.path(dir_in, sub("bed$", "rds", fileNameOut)))
# }
# 
# # The following example demonstrates how the UCSC gaps data were processed
# # All genomes
# genomes <- c("hg19", "hg38", "mm9", "mm10")
# # Process each genome
# for (genome_id in genomes) {
#   print(paste("Genome", genome_id))
#   # Get chromosome info
#   chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
#   # Get genome-specific gaps table
#   mySession <- browserSession()
#   genome(mySession) <- genome_id
#   # gaps <- getTable(ucscTableQuery(mySession, track = "gap"))
#   query <- ucscTableQuery(mySession, table = "gap")
#   gaps <- getTable(query)
#   # Process each gap type
#   for (gap_type in sort(unique(gaps$type))) {
#     # gap_type <- "heterochromatin"
#     gaps_selected <- gaps[gaps$type == gap_type, ]
#     gapsGR <- makeGRangesFromDataFrame(gaps_selected, keep.extra.columns = TRUE)
#     # Select seqinfo data for the gaps object
#     chrom_data_subset <- chrom_data[chrom_data$chrom %in% seqlevels(gapsGR), ]
#     chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
#     if (!all(seqlevels(gapsGR) == chrom_data_subset$chrom)) {
#       stop("Chromosome names do not match.")
#     }
#     # Assign seqinfo data
#     seqlengths(gapsGR) <- chrom_data_subset$size
#     isCircular(gapsGR) <- chrom_data_subset$circular
#     genome(gapsGR)     <- genome_id
#     # Construct file name, e.g., hg19.UCSC.gap_centromere.rds
#     fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".rds")
#     # Save as Rds object
#     saveRDS(object = gapsGR, file = file.path(dir_in, fileNameOut))
#     print(paste("Length", gap_type, length(gapsGR)))
#   }
# }
# 
# 
# # Get genome-specific gaps table
# mySession <- browserSession()
# genome(mySession) <- "hg19"
# query <- ucscTableQuery(mySession, table = "gap")
# gaps <- getTable(query)
# # Number of regions per gap type
# library(ggplot2)
# mtx_to_plot <- as.data.frame(table(gaps$type))
# colnames(mtx_to_plot) <- c("Type", "Number")
# mtx_to_plot <- mtx_to_plot[order(mtx_to_plot$Number), ]
# mtx_to_plot$Type <- factor(mtx_to_plot$Type, levels = mtx_to_plot$Type)
# ggplot(mtx_to_plot, aes(x = Number, y = Type, fill = Type)) +
#   geom_bar(stat="identity") +
#   theme_bw() + theme(legend.position = "none")
# ggsave("man/figures/excluderanges_hg19_gaps_number.png", width = 5, height = 2.5)



