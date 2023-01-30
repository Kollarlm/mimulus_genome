# Genespace
# Mimulus guttatus inland annual and coastal perennial genomes (L1 + S1)
# August 2022ß

#### Load in dependencies first!
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("rtracklayer")
# 
# # Install GENESPACE If it is not this version it will not run
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("jtlovell/GENESPACE@v0.9.3")

# GENESPACE
library(GENESPACE)

#Set Path
runwd <- file.path("/Users/lesliekollar/Desktop/mimulus_genome_R/mimulus_genome/Mimulus_guttatus_nanopore_genome_v1.2")
setwd(runwd)

gids <- c("S1","L1") #Whatever you list FIRST here will be the reference genome in the pan-genome function.

# Initialize GENESPACE
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  outgroup = NULL,
  ploidy = 1, # 8 and 14 generations inbred so consider it haploid. This can be changed for different ploidys for each species using the "rep()" function.
  diamondMode = "fast",
  orthofinderMethod = "fast",
  wd = runwd, 
  orthofinderInBlk = TRUE,
  overwrite = F,
  verbose = T,
  nCores = 4,
  minPepLen = 50,
  gffString = "gff", # How willl GENESPACE find the gff. I set it to look for gff but some people use other strings.
  pepString = "protein", #How will GENESPACE find the protein.fa files? I am using protein because it is in the file name. 
  #path2orthofinder = "/opt/anaconda3/envs/orthofinder/bin/orthofinder", # here is where you will change the path to your orthofinder executable if it works for you. 
  #/Users/lesliekollar/Applications/orthofinder/OrthoFinder/orthofinder/", 
  path2diamond = "diamond",
  path2mcscanx = "MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

# Parsing data
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "mRNA", 
  gffIdColumn = "ID",
  #gffStripText = "ID=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "ID=",
  troubleshoot = TRUE)


# Ran orthofinde on MSU's HPCC. 
# /opt/anaconda3/envs/orthofinder/bin/orthofinder -f /Users/lesliekollar/Desktop/Mimulus_guttatus_nanopore_genome/peptide -t 1 -a 1 -X -o /Users/lesliekollar/Desktop/Mimulus_guttatus_nanopore_genome/orthofinder
# gpar <- run_orthofinder(gsParam = gpar) -- will work if you call it in the terminal.


# Setting synteny parameters
gpar$params$nCores <- 1

#gpar <- set_syntenyParams(gpar, onlyOgAnchors = FALSE, blkSize = 2, nGaps = 20) # These are relaxed parameters as suggested by Issue 21 on the github page. 
# Normally, this would be... gpar <- synteny(gsParam = gpar)

gpar <- set_syntenyParams(gsParam = gpar)

# Running Synteny
#gpar <- synteny(gpar, overWrite = TRUE) # The overwrite function was suggested as a fix to the error. This didnt work even after clearing the results directory.

gpar <- synteny(gsParam = gpar)

tiff(filename = "S1_L1_plot.tiff",height = 20, width = 30, units="cm",
     compression = "lzw", res = 300)

plot_riparian(gpar,
              blackBg = FALSE,
              gapProp = 0.01)
dev.off()

## Plotting CHr8 only

regs <- data.frame(
  chr = c("chr8","chr8"),
  genome = c("S1","L1"),
  start= c(0, 0), #Walking out the start sites (i.e. getting smaller)
  end= c(7624865, 6576983),
  cols = c("gold", "green"))


tiff(filename = "S1_L1_chr8.tiff",height = 20, width = 30, units="cm",
     compression = "lzw", res = 300)

plot_riparian(gpar,
              onlyTheseRegions = regs)
dev.off()

# regs <- data.table(
#   genome= c("L1", "S1", "IM62"),
#   chr = c("chr8","chr8"),
#   start= c(824557, 820545), #Walking out the start sites (i.e. getting smaller)
#   end= c(7624865, 6576983), #Walking out the end sites (i.e. getting bigger)
#   cols = c("green", "yellow"))
# 
# plot_riparianHits(gpar, onlyTheseRegions = regs)

# ripSouceDat <- plot_riparianHits(
#   gpar, 
#   refGenome = "L1",
#   genomeIDs = c("L1", "S1"),
#   labelTheseGenomes = c("L1", "S1"),
#   gapProp = .01,
#   refChrCols = c("red", "orange", "yellow", "green", "blue", "purple", "pink", "red", "orange", "yellow", "green", "blue", "purple", "pink"),
#   blackBg = FALSE, 
#   returnSourceData = T, 
#   verbose = F)

pg <- pangenome(gpar)

pg_file <- fwrite(pg, "pg_file.csv")





