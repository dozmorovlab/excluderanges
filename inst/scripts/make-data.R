### =========================================================================
### excluderanges is an AnnotationHub package that stores genomic coordinates of
### problematic genomic regions as GRanges.
### -------------------------------------------------------------------------
###

# [ --- ] excluderanges package currently contains 19 Rds objects. The original and 
# processed data are available at 
# [ --- ] https://drive.google.com/drive/folders/124DZtsU0YVWqkb7dgu8Nk6b3N8-ShVSC?usp=sharing

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
# sed 's/\ /_/g' hg19.Yeo.eCLIP_blacklistregions.hg19.bed > temp && mv temp hg19.Yeo.eCLIP_blacklistregions.hg19.bed
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
# # hg38.Kundaje.with_Boyle.v2.bed.gz
# https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
#
# # hg38.Lareau.MT_excludable_set.bed.gz
# https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/encodeBlacklist/hg38.encode.blacklist.bed
# 
# # wgEncodeDacMapabilityConsensusExcludable.hg38.bed.gz
# # https://www.encodeproject.org/files/ENCFF220FIN/, Tim Reddy, Duke
# wget https://www.encodeproject.org/files/ENCFF220FIN/@@download/ENCFF220FIN.bed.gz
# gzip --decompress --stdout ENCFF220FIN.bed.gz | bedtools sort -i - > hg38.Reddy.wgEncodeDacMapabilityConsensusExcludable.hg38.bed
# 
# # hg38.Wimberley.peakpass_excludable_set.bed.gz
# https://raw.githubusercontent.com/ewimberley/peakPass/main/excludedlists/hg38/peakPass60Perc_sorted.bed
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
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(readr)
# Folder with results
# dir_in <- "/Users/mdozmorov/Documents/Data/GoogleDrive/excludedata"
dir_in <- "/Users/jogata/Documents/decluderanges_data/package_data/denydata"

# directory where output tables are saved
# save_dir <- file.path('~', 'Documents', 'Github', 'decluderanges.dev', 'manuscript')

exclude_information <- matrix(NA, nrow = 1, ncol = 5)
colnames(exclude_information) <- c('Object', 
                                   'Assembly', 
                                   'Number of regions', 
                                   'Min : median : max of set', 
                                   'Chromosomes with no regions'
                                   )


# All BED files
files <- list.files(path = dir_in, pattern = "bed$", ignore.case=T)
# In each subfolder
for (file in files) {
  
  # Read "fimo.bed" created by "fimo.qsub"
  excludeBED <- read.table(file.path(dir_in, file))
  
  # Assign column names depending on the number of columns
  if (ncol(excludeBED) == 3) { # Only 3 columns
    colnames(excludeBED) <- c("chr", "start", "stop")
  } else { # If more than 3 columns, consider the first 6
    colnames(excludeBED) <- c("chr", "start", "stop", "name", "score", "strand")
  }

  # Convert to GRanges object
  denyGR <- GenomicRanges::makeGRangesFromDataFrame(excludeBED, keep.extra.columns = TRUE)

  # Parse out genome ID from the file name, to get hg19, hg38, mm9, mm10, etc.
  genome_id <- strsplit(file, ".", fixed = TRUE)[[1]][1]

  # Get chromosome info and match it to the chromosome order in excludeBED
  if (genome_id=="TAIR10"){
    
    chrom_data <- GenomeInfoDb::getChromInfoFromNCBI("TAIR10") # TAIR10 not on UCSC
    main_chroms <- chrom_data[!grepl("_", chrom_data$SequenceName),] 
    chrom_data <- chrom_data[chrom_data$SequenceName %in% seqlevels(denyGR), ]
    chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$SequenceName), ]
    
    # Check if chromosome order is the same
    if (!all.equal(seqlevels(denyGR), chrom_data$SequenceName)) {
      print(paste("Chromosome order does not match for", genome_id, "genome."))
      break
    }
    seqlengths(denyGR) <- chrom_data$SequenceLength
    isCircular(denyGR) <- chrom_data$circular
    genome(denyGR)     <- genome_id
    
  } else if(genome_id=="T2T"){ 
    chrom_data <- GenomeInfoDb::getChromInfoFromNCBI("T2T-CHM13v2.0") # T2T not on UCSC
    main_chroms <- chrom_data[!grepl("_", chrom_data$SequenceName),]
    chrom_data <- chrom_data[chrom_data$SequenceName %in% seqlevels(denyGR), ]
    chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$SequenceName), ]
    
    # Check if chromosome order is the same
    if (!all.equal(seqlevels(denyGR), chrom_data$SequenceName)) {
      print(paste("Chromosome order does not match for", genome_id, "genome."))
      break
    }
    seqlengths(denyGR) <- chrom_data$SequenceLength
    isCircular(denyGR) <- chrom_data$circular
    genome(denyGR)     <- genome_id
  } else{
    
    chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
    main_chroms <- chrom_data[!grepl("_", chrom_data$chrom),]
    chrom_data <- chrom_data[chrom_data$chrom %in% seqlevels(denyGR), ]
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
  # 
  # Reformat output file name
  #fileNameOut <- sub("blacklist", "Excludable", file, ignore.case = TRUE)
  #fileNameOut <- sub("black", "Excludable", fileNameOut, ignore.case = TRUE)

  # Save as Rds object. file extension is changed to rds
  saveRDS(object = denyGR, file = file.path(dir_in, sub("bed$", "rds", file)))
  # excludeGR <- readRDS(file = file.path(dir_in, sub("bed$", "rds", fileNameOut)))

  # This small section is for creating a table on excludable set files
  # Get number of regions
  length <- length(denyGR)
  # Get minimum, median, and maximum region in each set
  mmm <-
    paste(
      format( min(width(denyGR)), big.mark = "," ),
      format( median(width(denyGR)), big.mark = "," ),
      format( max(width(denyGR)), big.mark = "," ),
      sep = " : "
    )
  # find missing chromosomes in excludable set
  missing <- gsub("chr", "",
                  paste(
                    as.character( main_chroms$chrom[ !(main_chroms$chrom %in% names( split( denyGR, seqnames(denyGR))))]),
                    collapse = ", "
                  )
  )
  if (missing != "") {missing = missing}
  else{missing="None"}

  # Append relevant information to dataframe
  d <- data.frame(file,
                   genome_id,
                   length,
                   mmm,
                   missing)

  colnames(d) <- c('Object',
                    'Assembly',
                    'Number of regions',
                    'Min : median : max of set',
                    'Chromosomes with no regions'
                   )
  # Add information to existing dataframe
  exclude_information <- rbind(exclude_information, d)
}

Source <- c(
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/ce10-blacklist.v2.bed.gz?raw=true', # ce10.Boyle_from_Excludable.v2.Excludable.bed
  'http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans', # ce10.Kundaje.ce10-Excludable.bed
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/ce11-blacklist.v2.bed.gz?raw=true', # ce11.Boyle_from_Excludable.v2.Excludable.bed
  'https://raw.githubusercontent.com/adomingues/redl_domingues_et_al_dev_2020/main/blacklisted.bed', # danRer10.Domingues.Excludable.bed
  'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8183574/bin/NIHMS1695157-supplement-Supplementary_Table_1-19.zip', # danRer10.Yang.Excludable.bed
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm3-blacklist.v2.bed.gz', # dm3.Boyle_from_Excludable.v2.Excludable.bed
  'http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/', # dm3.Kundaje.dm3-Excludable.bed
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz', # dm6.Boyle_from_Excludable.v2.Excludable.bed
  'https://www.encodeproject.org/files/ENCFF200UUD/', # hg19.Bernstein.Mint_Excludable_hg19.bed
  'https://www.encodeproject.org/files/ENCFF001TDO/', # hg19.Birney.wgEncodeDacMapabilityConsensusExcludable.bed
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz?raw=true', #hg19.Boyle_from_Excludable.v2.Excludable.bed
  'https://www.encodeproject.org/files/ENCFF001THR/', # hg19.Crawford.wgEncodeDukeMapabilityRegionsExcludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/mm10.full.blacklist.bed', # hg19.Lareau_full.Excludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/peaks/hg19_peaks.narrowPeak', # hg19.Lareau_MT.excludable.bed
  'https://www.encodeproject.org/files/ENCFF055QTV/', # hg19.Wold.hg19mitoExcludable.bed
  'https://www.encodeproject.org/files/ENCFF039QTN/', # hg19.Yeo.eCLIP_Excludableregions.hg19.bed
  'https://www.encodeproject.org/files/ENCFF023CZC/', # hg38.Bernstein.Mint_Excludable_GRCh38.bed
  'https://www.encodeproject.org/files/ENCFF356LFX/', # hg38.Kundaje.GRCh38_unified_Excludable.bed
  'https://www.encodeproject.org/files/ENCFF419RSJ/' , # hg38.Kundaje.GRCh38.Excludable.bed
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz?raw=true' , # hg38.Kundaje.with_Boyle.v2.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed' , # hg38.Lareau_full.Excludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/encodeBlacklist/hg38.encode.blacklist.bed' , # hg38.Lareau.MT_excludable_set.bed
  'https://www.encodeproject.org/files/ENCFF220FIN/' , # hg38.Reddy.wgEncodeDacMapabilityConsensusExcludable.hg38.bed
  'https://raw.githubusercontent.com/ewimberley/peakPass/main/excludedlists/hg38/peakPass60Perc_sorted.bed' , # hg38.Wimberley.peakpass_excludable_set.bed
  'https://www.encodeproject.org/files/ENCFF940NTE/' , # hg38.Wold.hg38mitoExcludable.bed
  'https://www.encodeproject.org/files/ENCFF269URO/' , # hg38.Yeo.eCLIP_Excludableregions.hg38liftover.bed.fixed.bed
  'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz?raw=true' , # mm10.Boyle_from_Excludable.v2.Excludable.bed
  'https://www.encodeproject.org/files/ENCFF790DJT/' , # mm10.Hardison.Excludable.full.bed
  'https://www.encodeproject.org/files/ENCFF226BDM/' , # mm10.Hardison.psuExcludable.mm10.bed
  'https://www.encodeproject.org/files/ENCFF999QPV/' , # mm10.Kundaje.anshul.Excludable.mm10.bed
  'https://www.encodeproject.org/files/ENCFF547MET/' , # mm10.Kundaje.mm10.Excludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/mm10.full.blacklist.bed' , # mm10.Lareau_full.Excludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/peaks/mm10_peaks.narrowPeak', # mm10.Lareau_MT.excludable.bed
  'https://www.encodeproject.org/files/ENCFF759PJK/' , # mm10.Wold.mm10mitoExcludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/mm9.full.blacklist.bed' , # mm9.Lareau_full.Excludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/peaks/mm9_peaks.narrowPeak', # mm9.Lareau_MT.excludable.bed
  'https://www.encodeproject.org/files/ENCFF299EZH/',  # mm9.Wold.mm9mitoExcludable.bed
  'New, ExcludeRanges', # T2T.Dozmorov-Ogata.excludable.bed
  'https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/peaks/chm13v2.0_peaks.narrowPeak', # T2T.Lareau_MT.excludable.bed
  'https://raw.githubusercontent.com/sklasfeld/GreenscreenProject/main/data/arabidopsis_blacklist_20inputs.bed', #TAIR10.Klasfeld_from_Excludable.Excludable.bed
  'https://github.com/sklasfeld/GreenscreenProject/blob/main/data/arabidopsis_greenscreen_20inputs.bed', #TAIR10.Klasfeld_from_Greenscreen.Excludable.bed
  'https://raw.githubusercontent.com/ewimberley/peakPass/main/excludedlists/tair10/predicted_excluded_list_sorted_0.6.bed' #TAIR10.Wimberley_peakpass.Excludable.bed
)

Year <- c(
  2018, # ce10.Boyle_from_Excludable.v2.Excludable.bed
  2012, # ce10.Kundaje.ce10-Excludable.bed
  2018, # ce11.Boyle_from_Excludable.v2.Excludable.bed
  2020, # danRer10.Domingues.Excludable.bed
  2021, # danRer10.Yang.Excludable.bed
  2018, # dm3.Boyle_from_Excludable.v2.Excludable.bed
  2012, # dm3.Kundaje.dm3-Excludable.bed
  2018, # dm6.Boyle_from_Excludable.v2.Excludable.bed
  2019, # hg19.Bernstein.Mint_Excludable_hg19.bed
  2011, # hg19.Birney.wgEncodeDacMapabilityConsensusExcludable.bed
  2018, # hg19.Boyle_from_Excludable.v2.Excludable.bed
  2011, # hg19.Crawford.wgEncodeDukeMapabilityRegionsExcludable.bed
  2017, # hg19.Lareau_full.Excludable.bed
  2017, # hg19.Lareau_MT.excludable.bed
  2016, # hg19.Wold.hg19mitoExcludable.bed
  2019, # hg19.Yeo.eCLIP_Excludableregions.hg19.bed
  2019, # hg38.Bernstein.Mint_Excludable_GRCh38.bed
  2020, # hg38.Kundaje.GRCh38_unified_Excludable.bed
  2016, # hg38.Kundaje.GRCh38.Excludable.bed
  2018, # hg38.Kundaje.with_Boyle.v2.bed
  2017, # hg38.Lareau_full.Excludable.bed
  2016, # hg38.Lareau.MT_excludable_set.bed
  2016, # hg38.Reddy.wgEncodeDacMapabilityConsensusExcludable.hg38.bed
  2021, # hg38.Wimberley.peakpass_excludable_set.bed
  2016, # hg38.Wold.hg38mitoExcludable.bed
  2019, # hg38.Yeo.eCLIP_Excludableregions.hg38liftover.bed.fixed.bed
  2018, # mm10.Boyle_from_Excludable.v2.Excludable.bed
  2016, # mm10.Hardison.Excludable.full.bed
  2016, # mm10.Hardison.psuExcludable.mm10.bed
  2016, # mm10.Kundaje.anshul.Excludable.mm10.bed
  2016, # mm10.Kundaje.mm10.Excludable.bed
  2017, # mm10.Lareau_full.Excludable.bed
  2017, # mm10.Lareau_MT.excludable.bed
  2016, # mm10.Wold.mm10mitoExcludable.bed
  2017, # mm9.Lareau_full.Excludable.bed
  2017, # mm9.Lareau_MT.excludable.bed
  2016, # mm9.Wold.mm9mitoExcludable.bed
  2022, # T2T.Dozmorov-Ogata.excludable.bed
  2022, # T2T.Lareau_MT.excludable.bed
  2021, #TAIR10.Klasfeld_from_Excludable.Excludable.bed
  2021, #TAIR10.Klasfeld_from_Greenscreen.Excludable.bed
  2021 #TAIR10.Wimberley_peakpass.Excludable.bed
)

exclude_information <- exclude_information[-1,]
exclude_information <- cbind(exclude_information, 'Year created or updated'=Year)
exclude_information <- cbind(exclude_information, Source)


# create empty dataframe, used for creating .csv file
gap_information <- matrix(NA, nrow = 1, ncol = 5)
colnames(gap_information) <- c('Object',
                               'Assembly',
                               'Number of regions',
                               'Min : median : max of set',
                               'Chromosomes with no regions'
                               )

# The following example demonstrates how the UCSC gaps data were processed
# All genomes
genomes <- c("hg19", "hg38", "mm9", "mm10")
# Process each genome
for (genome_id in genomes) {
  # Get chromosome info
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
  # Get genome-specific gaps table
  mySession <- browserSession()
  genome(mySession) <- genome_id
  # gaps <- getTable(ucscTableQuery(mySession, track = "gap"))
  query <- ucscTableQuery(mySession, table = "gap")
  gaps <- getTable(query)
  # Process each gap type
  for (gap_type in sort(unique(gaps$type))) {
    # gap_type <- "heterochromatin"
    gaps_selected <- gaps[gaps$type == gap_type, ]
    gapsGR <- makeGRangesFromDataFrame(gaps_selected, keep.extra.columns = TRUE)
    # Select seqinfo data for the gaps object
    chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome=genome_id)
    main_chroms <- chrom_data[!grepl("_", chrom_data$chrom),]
    chrom_data_subset <- chrom_data[chrom_data$chrom %in% seqlevels(gapsGR), ]
    chrom_data_subset <- chrom_data_subset[match(seqlevels(gapsGR), chrom_data_subset$chrom), ]
    if (!all(seqlevels(gapsGR) == chrom_data_subset$chrom)) {
      stop("Chromosome names do not match.")
    }
    # Assign seqinfo data
    seqlengths(gapsGR) <- chrom_data_subset$size
    isCircular(gapsGR) <- chrom_data_subset$circular
    genome(gapsGR)     <- genome_id
    # Construct file name, e.g., hg19.UCSC.gap_centromere.rds
    fileNameOut <- paste0(genome_id, ".UCSC.", gap_type, ".rds")
    # Save as Rds object
    saveRDS(object = gapsGR, file = file.path(dir_in, fileNameOut))

    # Get number of regions
    length <- length(gapsGR)
    # Get minimum, median, and maximum region in each set
    mmm <-
      paste(
        format( min(width(gapsGR)), big.mark = "," ),
        format( median(width(gapsGR)), big.mark = "," ),
        format( max(width(gapsGR)), big.mark = "," ),
        sep = " : "
      )
    # find missing chromosomes in excludable set
    missing <- gsub("chr", "",
                    paste(
                      as.character( main_chroms$chrom[ !(main_chroms$chrom %in% names( split( gapsGR, seqnames(gapsGR))))]),
                      collapse = ", "
                    )
    )
    if (missing != "") {missing = missing}
    else{missing="None"}

    df <- data.frame(fileNameOut, genome_id, length, mmm, missing)
    colnames(df) <- c('Object',
                      'Assembly',
                      'Number of regions',
                      'Min : median : max of set',
                      'Chromosomes with no regions'
                      )
    # Add information to existing dataframe
    gap_information <- rbind(gap_information, df)
  }
}

gapsGR <- genomation::readBed('~/Documents/decluderanges_data/package_data/gap_data/hg38.UCSC.centromere.bed')
fileNameOut <- 'hg38.UCSC.centromere.rds'
# Get number of regions
length <- length(gapsGR)
# Get minimum, median, and maximum region in each set
mmm <-
  paste(
    format( min(width(gapsGR)), big.mark = "," ),
    format( median(width(gapsGR)), big.mark = "," ),
    format( max(width(gapsGR)), big.mark = "," ),
    sep = " : "
  )
# find missing chromosomes in excludable set
missing <- gsub("chr", "",
                paste(
                  as.character( main_chroms$chrom[ !(main_chroms$chrom %in% names( split( gapsGR, seqnames(gapsGR))))]),
                  collapse = ", "
                )
)
if (missing != "") {missing = missing}else{
  missing="None"}
saveRDS(object = gapsGR, file = file.path(dir_in, fileNameOut))
df <- data.frame(fileNameOut, genome_id, length, mmm, missing)
colnames(df) <- c('Object',
                  'Assembly',
                  'Number of regions',
                  'Min : median : max of set',
                  'Chromosomes with no regions')
# Add information to existing dataframe
gap_information <- rbind(gap_information, df)

Source <- c( 'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm9&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm9&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm9&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema',
'http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1383614117_kpDvASW8YWQmcbxXyfWl9hv4fHnj&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_centromeres&hgta_ctDesc=table+browser+query+on+centromeres&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED'
)

Year <- c(
          rep(2020, 7),
          rep(2018, 5),
          rep(2021, 7),
          rep(2007, 3),
          2021,
          2022
)

gap_information <- gap_information[-1,]
gap_information <- cbind(gap_information, 'Year created or updated'=Year)
gap_information <- cbind(gap_information, Source)
gap_information <- gap_information[order(gap_information$Object), ]

information <- rbind(exclude_information, gap_information)

writexl::write_xlsx(
  information,
  file.path('inst', 'extdata', 'table_excluderanges_and_gaps.xlsx')
)

# writexl::write_xlsx(
#   information,
#   file.path(save_dir, 'Table_S1.xlsx')
# )

# # Number of samples per cell/tissue type?
# library(ggplot2)
# mtx_to_plot <- as.data.frame(table(gaps$type))
# colnames(mtx_to_plot) <- c("Type", "Number")
# mtx_to_plot <- mtx_to_plot[order(mtx_to_plot$Number), ]
# mtx_to_plot$Type <- factor(mtx_to_plot$Type, levels = mtx_to_plot$Type)
# ggplot(mtx_to_plot, aes(x = Number, y = Type, fill = Type)) +
#   geom_bar(stat="identity") +
#   theme_bw() + theme(legend.position = "none")
# ggsave("man/figures/excluderanges_hg19_gaps_number.png", width = 5, height = 2.5)



