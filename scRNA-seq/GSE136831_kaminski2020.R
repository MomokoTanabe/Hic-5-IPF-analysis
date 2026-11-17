# _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
# IPF scRNA-seq SciAdv2020 Kaminski
# _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
# GSE136831
# Control n=28 (65916 cells), IPF n=32 (109825 cells)


# =========================================
# library
# =========================================
library(Seurat)
library(dplyr)
library(readr)
library(R.utils)
library(ggplot2)
library(patchwork)
library(Matrix)


# =========================================
# 1. Sample loading
# =========================================
# Extract the downloaded RDS.gz files
# Run the following command in the terminal:
# gunzip GSE136831_RAW/*.rds.gz

# List the extracted .rds files
samples <- list.files(
  "GSE136831_RAW",
  pattern = "\\.rds$",
  full.names = TRUE
)

# Load the files into a list
seurat_list <- lapply(samples, readRDS)

# Assign sample names to the list elements
names(seurat_list) <- gsub("\\.rds$", "", basename(samples))

# The RDS files contain UMI count and read count data.
# Use UMI counts for exonic regions.


# =========================================
# 2. Create and save Seurat objects one sample at a time
# =========================================
# Check the count matrix structure
example_sample <- seurat_list[[1]]
str(example_sample$umicount$exon, max.level = 2)
names(example_sample$umicount$exon)

# Use all UMI counts for exonic regions
# Process one sample at a time to reduce memory usage.

# Template for processing individual samples
# Replace "sample1" with the corresponding sample identifier.

rds_file <- "GSE136831_RAW/sample1.dgecounts.rds"  # Path to the input RDS file
sample_list <- readRDS(rds_file)  # Load the RDS file
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")

saveRDS(
  seu,
  file = "Seurat_samples/sample1_Seurat.rds"
)  # Replace "sample1" with the corresponding sample identifier.

rm(seu, counts)
gc()  # Release memory


# -----------------------------------------
# Sample: GSM4058900_133C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058900_133C-a.dgecounts.rds"
sample_list <- readRDS(rds_file) 
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058900_133C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058901_137C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058901_137C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058901_137C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058901_137C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058901_137C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058901_137C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058902_034C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058902_034C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058902_034C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058903_218C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058903_218C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058903_218C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058903_218C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058903_218C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058903_218C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058903_219C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058903_219C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058903_219C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058903_219C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058903_219C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058903_219C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058904_226C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058904_226C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058904_226C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058904_226C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058904_226C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058904_226C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058905_244C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058905_244C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058905_244C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058906_098C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058906_098C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058906_098C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058906_098C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058906_098C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058906_098C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058907_465C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058907_465C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058907_465C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058908_396C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058908_396C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058908_396C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058909_296C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058909_296C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058909_296C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058910_208C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058910_208C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058910_208C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058911_022C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058911_022C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058911_022C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058911_022C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058911_022C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058911_022C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058911_222C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058911_222C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058911_222C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058912_160C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058912_160C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058912_160C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058913_092C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058913_092C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058913_092C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058914_439C-b
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058914_439C-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058914_439C-b.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058914_439C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058914_439C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058914_439C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058915_065C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058915_065C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058915_065C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058916_388C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058916_388C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058916_388C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058917_92C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058917_92C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058917_92C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058917_192C-a
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058917_192C-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058917_192C-a.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058918_483C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058918_483C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058918_483C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058919_001C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058919_001C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058919_001C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058920_002C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058920_002C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058920_002C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058921_003C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058921_003C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058921_003C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058922_454C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058922_454C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058922_454C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058923_253C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058923_253C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058923_253C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058924_484C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058924_484C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058924_484C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058925_081C
# -----------------------------------------

rds_file <- "GSE136831_RAW/GSM4058925_081C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(
  seu,
  file = "rds_files/GSM4058925_081C.rds"
)
rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058926_137C
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058926_137C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058926_137C.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058927_084C
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058927_084C.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058927_084C.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058946_210I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058946_210I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058946_210I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058947_135I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058947_135I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058947_135I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058947_135I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058947_135I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058947_135I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058948_138I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058948_138I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058948_138I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058948_138I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058948_138I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058948_138I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058949_145I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058949_145I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058949_145I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058950_157I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058950_157I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058950_157I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058950_157I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058950_157I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058950_157I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058950_157I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058950_157I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058950_157I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058951_166I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058951_166I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058951_166I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058952_174I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058952_174I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058952_174I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058953_179I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058953_179I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058953_179I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058954_209I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058954_209I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058954_209I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058955_212I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058955_212I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058955_212I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058956_214I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058956_214I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058956_214I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058957_221I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058957_221I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058957_221I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058958_222I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058958_222I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058958_222I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058958_222I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058958_222I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058958_222I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058959_225I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058959_225I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058959_225I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058959_225I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058959_225I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058959_225I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058960_228I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058960_228I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058960_228I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058960_228I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058960_228I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058960_228I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058961_29I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058961_29I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058961_29I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058962_051I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058962_051I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058962_051I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058962_051I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058962_051I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058962_051I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058963_025I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058963_025I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058963_025I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058964_010I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058964_010I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058964_010I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058965_021I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058965_021I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058965_021I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058966_022I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058966_022I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058966_022I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058967_041I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058967_041I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058967_041I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058968_040I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058968_040I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058968_040I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058968_040I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058968_040I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058968_040I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058969_47I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058969_47I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058969_47I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058969_47I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058969_47I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058969_47I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058970_49I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058970_49I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058970_49I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058970_49I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058970_49I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058970_49I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058971_053I-d
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058971_053I-d.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058971_053I-d.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058971_053I-n
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058971_053I-n.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058971_053I-n.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058971_053I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058971_053I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058971_053I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058972_59I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058972_59I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058972_59I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058973_063I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058973_063I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058973_063I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058974_034I-a
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058974_034I-a.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058974_034I-a.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058975_123I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058975_123I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058975_123I.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058976_158I-b
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058976_158I-b.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058976_158I-b.rds"
)

rm(seu, counts)
gc()


# -----------------------------------------
# Sample: GSM4058977_177I
# -----------------------------------------
rds_file <- "GSE136831_RAW/GSM4058977_177I.dgecounts.rds"
sample_list <- readRDS(rds_file)
counts <- sample_list$umicount$exon$all
seu <- CreateSeuratObject(counts = counts)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(
  seu,
  file = "rds_files/GSM4058977_177I.rds"
)

rm(seu, counts)
gc()


# =========================================
# 3. Load Seurat objects and generate QC summary
# =========================================
library(purrr)
library(Seurat)
library(biomaRt)
library(tibble)

# Load individual Seurat objects
seurat_files <- list.files(
  "rds_files/", 
  pattern = "\\.rds$", 
  full.names = TRUE
)

seurat_objs <- setNames(
  lapply(seurat_files, readRDS),
  tools::file_path_sans_ext(basename(seurat_files))
)

# Recalculate mitochondrial gene percentages using Ensembl IDs.
# The original count matrices use Ensembl gene IDs rather than gene symbols.
# Therefore, mitochondrial genes must first be identified by mapping
# Ensembl gene IDs to HGNC symbols.


# Retrieve Ensembl gene ID-to-HGNC symbol annotations
mart <- useMart(
  "ensembl", 
  dataset = "hsapiens_gene_ensembl"
)

annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
)

# Identify Ensembl IDs corresponding to mitochondrial genes
mt_ensembl <- annot$ensembl_gene_id[
  grep("^MT-", annot$hgnc_symbol)
]

# Recalculate the percentage of mitochondrial gene counts
seurat_objs <- lapply(seurat_objs, function(seu) {
  mt_present <- intersect(rownames(seu), mt_ensembl)
  
  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu, 
    features = mt_present
  )
  
  return(seu)
})

## Generate a QC summary for each sample
qc_summary <- map_dfr(names(seurat_objs), function(sample_name) {
  seu <- seurat_objs[[sample_name]]
  meta <- seu@meta.data
  
  tibble(
    Sample = sample_name,
    Batch = stringr::str_extract(sample_name, "[A-Za-z]$"), 
    Cells = nrow(meta),
    MeanGenes = round(mean(meta$nFeature_RNA), 1),
    MedianGenes = round(median(meta$nFeature_RNA), 1),
    PercentMT = round(median(meta$percent.mt), 2)
  )
})

# Display and save the QC summary
qc_summary

write.csv(
  qc_summary, 
  "QC_summary_with_batch.csv", 
  row.names = FALSE
)

# Select samples for downstream analysis
# Control n=28, IPF n=32


# =========================================
# 4. QC: Visualization and filtering
# =========================================

# For each sample:
# 1. Retrieve the Seurat object.
# 2. Visualize QC metrics using violin plots.
# 3. Filter cells based on the sample-specific thresholds.
# 4. Normalize the data and identify highly variable features.
# 5. Update the Seurat object in the list.

# -----------------------------------------
# Sample: GSM4058900_133C-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058900_133C-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058900_133C-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058901_137C-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058901_137C-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058901_137C-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058902_034C
# -----------------------------------------

seu <- seurat_objs[["GSM4058902_034C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 20
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058902_034C"]] <- seu


# -----------------------------------------
# Sample: GSM4058903_219C-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058903_219C-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058903_219C-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058904_226C-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058904_226C-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058904_226C-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058905_244C
# -----------------------------------------

seu <- seurat_objs[["GSM4058905_244C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058905_244C"]] <- seu


# -----------------------------------------
# Sample: GSM4058906_098C-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058906_098C-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058906_098C-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058907_465C
# -----------------------------------------

seu <- seurat_objs[["GSM4058907_465C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058907_465C"]] <- seu


# -----------------------------------------
# Sample: GSM4058908_396C
# -----------------------------------------

seu <- seurat_objs[["GSM4058908_396C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2000 &
    percent.mt < 20
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058908_396C"]] <- seu


# -----------------------------------------
# Sample: GSM4058909_296C
# -----------------------------------------

seu <- seurat_objs[["GSM4058909_296C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058909_296C"]] <- seu


# -----------------------------------------
# Sample: GSM4058910_208C
# -----------------------------------------

seu <- seurat_objs[["GSM4058910_208C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058910_208C"]] <- seu


# -----------------------------------------
# Sample: GSM4058911_022C-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058911_022C-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058911_022C-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058912_160C
# -----------------------------------------

seu <- seurat_objs[["GSM4058912_160C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058912_160C"]] <- seu


# -----------------------------------------
# Sample: GSM4058913_092C
# -----------------------------------------

seu <- seurat_objs[["GSM4058913_092C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058913_092C"]] <- seu


# -----------------------------------------
# Sample: GSM4058914_439C-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058914_439C-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058914_439C-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058915_065C
# -----------------------------------------

seu <- seurat_objs[["GSM4058915_065C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058915_065C"]] <- seu


# -----------------------------------------
# Sample: GSM4058916_388C
# -----------------------------------------

seu <- seurat_objs[["GSM4058916_388C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058916_388C"]] <- seu


# -----------------------------------------
# Sample: GSM4058917_192C-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058917_192C-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058917_192C-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058918_483C
# -----------------------------------------

seu <- seurat_objs[["GSM4058918_483C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058918_483C"]] <- seu


# -----------------------------------------
# Sample: GSM4058919_001C
# -----------------------------------------

seu <- seurat_objs[["GSM4058919_001C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058919_001C"]] <- seu



# -----------------------------------------
# Sample: GSM4058920_002C
# -----------------------------------------

seu <- seurat_objs[["GSM4058920_002C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058920_002C"]] <- seu


# -----------------------------------------
# Sample: GSM4058921_003C
# -----------------------------------------

seu <- seurat_objs[["GSM4058921_003C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058921_003C"]] <- seu


# -----------------------------------------
# Sample: GSM4058922_454C
# -----------------------------------------

seu <- seurat_objs[["GSM4058922_454C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058922_454C"]] <- seu


# -----------------------------------------
# Sample: GSM4058923_253C
# -----------------------------------------

seu <- seurat_objs[["GSM4058923_253C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058923_253C"]] <- seu


# -----------------------------------------
# Sample: GSM4058924_484C
# -----------------------------------------

seu <- seurat_objs[["GSM4058924_484C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058924_484C"]] <- seu


# -----------------------------------------
# Sample: GSM4058925_081C
# -----------------------------------------

seu <- seurat_objs[["GSM4058925_081C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058925_081C"]] <- seu


# -----------------------------------------
# Sample: GSM4058926_137C
# -----------------------------------------

seu <- seurat_objs[["GSM4058926_137C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058926_137C"]] <- seu


# -----------------------------------------
# Sample: GSM4058927_084C
# -----------------------------------------

seu <- seurat_objs[["GSM4058927_084C"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058927_084C"]] <- seu


# -----------------------------------------
# Sample: GSM4058946_210I
# -----------------------------------------

seu <- seurat_objs[["GSM4058946_210I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058946_210I"]] <- seu


# -----------------------------------------
# Sample: GSM4058947_135I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058947_135I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058947_135I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058948_138I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058948_138I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058948_138I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058949_145I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058949_145I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058949_145I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058950_157I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058950_157I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058950_157I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058951_166I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058951_166I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058951_166I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058952_174I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058952_174I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058952_174I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058953_179I
# -----------------------------------------

seu <- seurat_objs[["GSM4058953_179I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058953_179I"]] <- seu


# -----------------------------------------
# Sample: GSM4058954_209I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058954_209I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058954_209I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058955_212I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058955_212I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058955_212I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058956_214I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058956_214I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3500 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058956_214I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058957_221I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058957_221I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058957_221I-a"]] <- seu



# -----------------------------------------
# Sample: GSM4058958_222I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058958_222I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058958_222I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058959_225I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058959_225I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058959_225I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058960_228I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058960_228I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058960_228I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058961_29I
# -----------------------------------------

seu <- seurat_objs[["GSM4058961_29I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058961_29I"]] <- seu


# -----------------------------------------
# Sample: GSM4058962_051I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058962_051I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058962_051I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058963_025I
# -----------------------------------------

seu <- seurat_objs[["GSM4058963_025I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058963_025I"]] <- seu


# -----------------------------------------
# Sample: GSM4058964_010I
# -----------------------------------------

seu <- seurat_objs[["GSM4058964_010I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058964_010I"]] <- seu


# -----------------------------------------
# Sample: GSM4058965_021I
# -----------------------------------------

seu <- seurat_objs[["GSM4058965_021I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058965_021I"]] <- seu


# -----------------------------------------
# Sample: GSM4058966_022I
# -----------------------------------------

seu <- seurat_objs[["GSM4058966_022I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058966_022I"]] <- seu


# -----------------------------------------
# Sample: GSM4058967_041I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058967_041I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058967_041I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058968_040I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058968_040I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058968_040I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058969_47I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058969_47I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058969_47I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058970_49I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058970_49I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058970_49I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058971_053I-n
# -----------------------------------------

seu <- seurat_objs[["GSM4058971_053I-n"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058971_053I-n"]] <- seu


# -----------------------------------------
# Sample: GSM4058972_59I
# -----------------------------------------

seu <- seurat_objs[["GSM4058972_59I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058972_59I"]] <- seu


# -----------------------------------------
# Sample: GSM4058973_063I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058973_063I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058973_063I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058974_034I-a
# -----------------------------------------

seu <- seurat_objs[["GSM4058974_034I-a"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058974_034I-a"]] <- seu


# -----------------------------------------
# Sample: GSM4058975_123I
# -----------------------------------------

seu <- seurat_objs[["GSM4058975_123I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058975_123I"]] <- seu


# -----------------------------------------
# Sample: GSM4058976_158I-b
# -----------------------------------------

seu <- seurat_objs[["GSM4058976_158I-b"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058976_158I-b"]] <- seu


# -----------------------------------------
# Sample: GSM4058977_177I
# -----------------------------------------

seu <- seurat_objs[["GSM4058977_177I"]]

VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 15
)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

seurat_objs[["GSM4058977_177I"]] <- seu


# Save the processed Seurat objects
saveRDS(
  seurat_objs,
  "rds_files/seurat_objs.rds"
)

# Reload the processed Seurat objects
seurat_objs <- readRDS("rds_files/seurat_objs.rds")


# =========================================
# 5. Integration of Seurat objects
# =========================================

# -----------------------------------------
# Control samples
# -----------------------------------------

# Define sample names for the Control group
control.samples <- c(
  "GSM4058900_133C-a",
  "GSM4058901_137C-a",
  "GSM4058902_034C",
  "GSM4058903_219C-a",
  "GSM4058904_226C-b",
  "GSM4058905_244C",
  "GSM4058906_098C-b",
  "GSM4058907_465C",
  "GSM4058908_396C",
  "GSM4058909_296C",
  "GSM4058910_208C",
  "GSM4058911_022C-a",
  "GSM4058912_160C",
  "GSM4058913_092C",
  "GSM4058914_439C-b",
  "GSM4058915_065C",
  "GSM4058916_388C",
  "GSM4058917_192C-a",
  "GSM4058918_483C",
  "GSM4058919_001C",
  "GSM4058920_002C",
  "GSM4058921_003C",
  "GSM4058922_454C",
  "GSM4058923_253C",
  "GSM4058924_484C",
  "GSM4058925_081C",
  "GSM4058926_137C",
  "GSM4058927_084C"
)

# Check whether all specified samples are available
stopifnot(all(control.samples %in% names(seurat_objs)))

# Extract Seurat objects for the Control group
control.list <- seurat_objs[control.samples]

# Add sample IDs to the metadata of each object
for (i in seq_along(control.list)) {
  sample_name <- control.samples[i]
  control.list[[i]]$sample_id <- sample_name
}

# Split the Control samples into groups of five to avoid memory issues
control.chunks <- split(
  control.samples,
  ceiling(seq_along(control.samples) / 5)
)

# Check the number of chunks
length(control.chunks)  # 6 (the last chunk contains 3 samples)

# Integrate and save each chunk separately
control.integrated.list <- list()

for (i in seq_along(control.chunks)) {
  message("Processing Control chunk ", i)
  
  sub.list <- control.list[control.chunks[[i]]]
  
  features <- SelectIntegrationFeatures(
    object.list = sub.list,
    nfeatures = 5000
  )
  
  anchors <- FindIntegrationAnchors(
    object.list = sub.list,
    anchor.features = features,
    dims = 1:30
  )
  
  integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:30
  )
  
  # Save each integrated chunk for subsequent integration
  saveRDS(
    integrated,
    paste0("rds_files/control_integrated_chunk", i, ".rds")
  )
  
  control.integrated.list[[i]] <- integrated
}

# Reload the individually integrated chunks
chunk.files <- list.files(
  "rds_files",
  pattern = "^control_integrated_chunk[0-9]+\\.rds$",
  full.names = TRUE
)

# Ensure that chunks are loaded in numerical order
chunk.files <- chunk.files[
  order(as.integer(sub(
    ".*chunk([0-9]+)\\.rds$",
    "\\1",
    chunk.files
  )))
]

control.integrated.list <- lapply(chunk.files, readRDS)

# Integrate the previously integrated chunks
features <- SelectIntegrationFeatures(
  object.list = control.integrated.list,
  nfeatures = 5000
)

anchors <- FindIntegrationAnchors(
  object.list = control.integrated.list,
  anchor.features = features,
  dims = 1:30
)

# Save and reload the integration anchors
saveRDS(anchors, "rds_files/anchors_control.rds")

anchors <- readRDS("rds_files/anchors_control.rds")

# Generate the final integrated Seurat object
control.integrated.final <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# Save the final integrated Seurat object for the Control group
saveRDS(
  control.integrated.final,
  "rds_files/control_integrated_final.rds"
)


# -----------------------------------------
# IPF samples
# -----------------------------------------

# Define sample names for the IPF group
IPF.samples <- c(
  "GSM4058946_210I",
  "GSM4058947_135I-a",
  "GSM4058948_138I-a",
  "GSM4058949_145I-a",
  "GSM4058950_157I-a",
  "GSM4058951_166I-a",
  "GSM4058952_174I-a",
  "GSM4058953_179I",
  "GSM4058954_209I-a",
  "GSM4058955_212I-a",
  "GSM4058956_214I-a",
  "GSM4058957_221I-a",
  "GSM4058958_222I-b",
  "GSM4058959_225I-b",
  "GSM4058960_228I-b",
  "GSM4058961_29I",
  "GSM4058962_051I-a",
  "GSM4058963_025I",
  "GSM4058964_010I",
  "GSM4058965_021I",
  "GSM4058966_022I",
  "GSM4058967_041I-b",
  "GSM4058968_040I-b",
  "GSM4058969_47I-b",
  "GSM4058970_49I-a",
  "GSM4058971_053I-n",
  "GSM4058972_59I",
  "GSM4058973_063I-b",
  "GSM4058974_034I-a",
  "GSM4058975_123I",
  "GSM4058976_158I-b",
  "GSM4058977_177I"
)

# Check whether all specified samples are available
stopifnot(all(IPF.samples %in% names(seurat_objs)))

# Extract Seurat objects for the IPF group
IPF.list <- seurat_objs[IPF.samples]

# Add sample IDs to the metadata of each object
for (i in seq_along(IPF.list)) {
  sample_name <- IPF.samples[i]
  IPF.list[[i]]$sample_id <- sample_name
}

# Split the IPF samples into groups of five
IPF.chunks <- split(
  IPF.samples,
  ceiling(seq_along(IPF.samples) / 5)
)

# Check the number of chunks
length(IPF.chunks)  # 7 (the last chunk contains 2 samples)

# Integrate and save each chunk separately
IPF.integrated.list <- list()

for (i in seq_along(IPF.chunks)) {
  message("Processing IPF chunk ", i)
  
  sub.list <- IPF.list[IPF.chunks[[i]]]
  
  features <- SelectIntegrationFeatures(
    object.list = sub.list,
    nfeatures = 5000
  )
  
  anchors <- FindIntegrationAnchors(
    object.list = sub.list,
    anchor.features = features,
    dims = 1:30
  )
  
  integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:30
  )
  
  # Save each integrated chunk for subsequent integration
  saveRDS(
    integrated,
    paste0("rds_files/IPF_integrated_chunk", i, ".rds")
  )
  
  IPF.integrated.list[[i]] <- integrated
}

# Reload the individually integrated chunks
chunk.files <- list.files(
  "rds_files",
  pattern = "^IPF_integrated_chunk[0-9]+\\.rds$",
  full.names = TRUE
)

# Ensure that chunks are loaded in numerical order
chunk.files <- chunk.files[
  order(as.integer(sub(
    ".*chunk([0-9]+)\\.rds$",
    "\\1",
    chunk.files
  )))
]

IPF.integrated.list <- lapply(chunk.files, readRDS)

# Integrate the previously integrated chunks
features <- SelectIntegrationFeatures(
  object.list = IPF.integrated.list,
  nfeatures = 5000
)

anchors <- FindIntegrationAnchors(
  object.list = IPF.integrated.list,
  anchor.features = features,
  dims = 1:30
)

# Save and reload the integration anchors
saveRDS(anchors, "rds_files/anchors_IPF.rds")

anchors <- readRDS("rds_files/anchors_IPF.rds")

# Generate the final integrated Seurat object
IPF.integrated.final <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# Save the final integrated Seurat object for the IPF group
saveRDS(
  IPF.integrated.final,
  "rds_files/IPF_integrated_final.rds"
)


# -----------------------------------------
# Integration of Control and IPF groups
# -----------------------------------------

# Load the final integrated Seurat objects for each group
control.integrated.final <- readRDS(
  "rds_files/control_integrated_final.rds"
)

IPF.integrated.final <- readRDS(
  "rds_files/IPF_integrated_final.rds"
)

# Add group labels to the metadata
control.integrated.final$group <- "Control"
IPF.integrated.final$group <- "IPF"

# Normalize, identify highly variable features, scale, and perform PCA
# separately for each Seurat object to reduce memory usage.

# Control group
control.integrated.final <- NormalizeData(
  control.integrated.final
)

control.integrated.final <- FindVariableFeatures(
  control.integrated.final,
  selection.method = "vst",
  nfeatures = 5000
)

control.integrated.final <- ScaleData(
  control.integrated.final
)

control.integrated.final <- RunPCA(
  control.integrated.final,
  features = VariableFeatures(control.integrated.final)
)

# IPF group
IPF.integrated.final <- NormalizeData(
  IPF.integrated.final
)

IPF.integrated.final <- FindVariableFeatures(
  IPF.integrated.final,
  selection.method = "vst",
  nfeatures = 5000
)

IPF.integrated.final <- ScaleData(
  IPF.integrated.final
)

IPF.integrated.final <- RunPCA(
  IPF.integrated.final,
  features = VariableFeatures(IPF.integrated.final)
)

# Select integration features shared between the two groups
combined.list <- list(
  control.integrated.final,
  IPF.integrated.final
)

features <- SelectIntegrationFeatures(
  object.list = combined.list,
  nfeatures = 5000
)

# Identify integration anchors using reciprocal PCA (RPCA)
anchors <- FindIntegrationAnchors(
  object.list = combined.list,
  anchor.features = features,
  dims = 1:30,
  reduction = "rpca"
)

# Integrate the Control and IPF groups
seurat.integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# Save the final integrated Seurat object
saveRDS(
  seurat.integrated,
  "rds_files/control_IPF_integrated.rds"
)

DefaultAssay(seurat.integrated)
# "integrated" ← OK


# =========================================
# 6. Scaling, clustering, and UMAP visualization
# =========================================

# Scale the integrated data
seurat.integrated <- ScaleData(seurat.integrated)

# Perform PCA
seurat.integrated <- RunPCA(seurat.integrated)

# Visualize the standard deviations of the principal components
ElbowPlot(seurat.integrated)

# The elbow point was identified at PC 10.

# Perform graph-based clustering
# The current analysis uses the first 10 principal components
# and a resolution of 0.5.
#
# In the original paper, the parameters for fibroblast and
# myofibroblast analysis were:
# nPCs = 25, k.param = 10, resolution = 0.7.

seurat.integrated <- FindNeighbors(
  seurat.integrated,
  dims = 1:10
)

seurat.integrated <- FindClusters(
  seurat.integrated,
  resolution = 0.5
)

# Perform UMAP dimensionality reduction
# The current analysis uses the first 10 principal components.
#
# The original paper used spread = 0.1 for UMAP visualization.

seurat.integrated <- RunUMAP(
  seurat.integrated,
  dims = 1:10
)

# Visualize the UMAP plot
DimPlot(
  seurat.integrated,
  reduction = "umap",
  label = TRUE
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )

# Check the dimensions of the integrated Seurat object
dim(seurat.integrated)  # Features (genes) x cells

# Check the number of cells in each group
table(seurat.integrated$group)  # Control and IPF

# Note:
# TGFB1I1 is present in the RNA assay but was not included
# among the 5,000 highly variable features.
# Therefore, TGFB1I1 will be manually added to the feature set.


# -----------------------------------------
# Add TGFB1I1 to VariableFeatures
# -----------------------------------------

# Retrieve the Ensembl gene ID for TGFB1I1
tgfb1i1_ensembl <- bitr(
  "TGFB1I1",
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = org.Hs.eg.db
)$ENSEMBL

tgfb1i1_ensembl

# Confirm that TGFB1I1 is present in the RNA assay
tgfb1i1_id <- "ENSG00000140682"

tgfb1i1_id %in% rownames(seurat.integrated[["RNA"]])

# Add TGFB1I1 to the existing VariableFeatures
var_genes <- VariableFeatures(seurat.integrated)

# Add TGFB1I1 only if it is not already included
if (!tgfb1i1_id %in% var_genes) {
  var_genes <- c(var_genes, tgfb1i1_id)
}

# Update VariableFeatures in the Seurat object
VariableFeatures(seurat.integrated) <- var_genes

# Recalculate scaling, PCA, and UMAP

# Scale the selected features
seurat.integrated <- ScaleData(
  seurat.integrated,
  features = var_genes
)

# Recalculate PCA
seurat.integrated <- RunPCA(
  seurat.integrated,
  features = var_genes
)

# Recalculate UMAP
seurat.integrated <- RunUMAP(
  seurat.integrated,
  dims = 1:25
)

# Save the updated integrated Seurat object
saveRDS(
  seurat.integrated,
  "rds_files/control_IPF_integrated.rds"
)


# -----------------------------------------
# Visualization of TGFB1I1 expression
# -----------------------------------------

# Visualize the overall UMAP
DimPlot(
  seurat.integrated,
  reduction = "umap",
  label = TRUE
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )

# Visualize TGFB1I1 expression across all cells
tgfb1i1_id <- "ENSG00000140682"

FeaturePlot(
  seurat.integrated,
  features = tgfb1i1_id
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )

# Set the default assay to RNA for gene expression visualization
DefaultAssay(seurat.integrated) <- "RNA"

# Visualize TGFB1I1 expression in Control cells
control_cells <- subset(
  seurat.integrated,
  subset = group == "Control"
)

FeaturePlot(
  control_cells,
  features = tgfb1i1_id
) +
  ggtitle("Control")

# Visualize TGFB1I1 expression in IPF cells
ipf_cells <- subset(
  seurat.integrated,
  subset = group == "IPF"
)

FeaturePlot(
  ipf_cells,
  features = tgfb1i1_id
) +
  ggtitle("IPF")

# Compare TGFB1I1 expression between Control and IPF
# using identical color scales and point sizes
library(patchwork)

# Control
p1 <- FeaturePlot(
  seurat.integrated[, seurat.integrated$group == "Control"],
  features = tgfb1i1_id,
  pt.size = 0.1,
  min.cutoff = 0,
  max.cutoff = 2,
  raster = FALSE
) +
  ggtitle("Control")

# IPF
p2 <- FeaturePlot(
  seurat.integrated[, seurat.integrated$group == "IPF"],
  features = tgfb1i1_id,
  pt.size = 0.1,
  min.cutoff = 0,
  max.cutoff = 2,
  raster = FALSE
) +
  ggtitle("IPF")

# Combine the plots
(p1 | p2) +
  plot_annotation(title = "TGFB1I1") +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )

# -----------------------------------------
# Downsampling to balance cell numbers
# -----------------------------------------

# Randomly downsample IPF cells to match the number of Control cells
set.seed(123)

cells.control <- colnames(seurat.integrated)[
  seurat.integrated$group == "Control"
]

cells.ipf <- colnames(seurat.integrated)[
  seurat.integrated$group == "IPF"
]

downsample.ipf <- sample(
  cells.ipf,
  size = length(cells.control)
)

# Subset the Seurat object to include all Control cells
# and the downsampled IPF cells
seurat.down <- subset(
  seurat.integrated,
  cells = c(cells.control, downsample.ipf)
)

# Check the number of cells in each group
table(seurat.down$group)

# Save the downsampled Seurat object
saveRDS(
  seurat.down,
  "rds_files/seurat.down.rds"
)


# =========================================
# 7. FeaturePlot
# =========================================

p <- FeaturePlot(
  seurat.down,
  features = tgfb1i1_id,
  split.by = "group",
  pt.size = 0.1,
  min.cutoff = 0,
  max.cutoff = 2
)

p + plot_annotation(
  title = "TGFB1I1") &
  theme(
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  ) &
  labs(x = "UMAP_1", y = "UMAP_2")



# =========================================
# 8. Cell type annotation
# =========================================

# -----------------------------------------
# Prepare gene expression data for CellMarker 2.0
# Convert gene identifiers to gene symbols
# -----------------------------------------

# Replace missing or empty gene symbols
missing_symbols <- is.na(rownames(expr_all)) | rownames(expr_all) == ""

rownames(expr_all)[missing_symbols] <- paste0(
  "Gene",
  seq_len(sum(missing_symbols))
)

# Make gene symbols unique
rownames(expr_all) <- make.unique(rownames(expr_all))

# Convert the expression matrix to a sparse matrix
expr_sparse <- as(expr_all, "dgCMatrix")

# Add the expression matrix as a new assay to the Seurat object
seurat.integrated[["expr_all"]] <- CreateAssayObject(
  counts = expr_sparse
)

# Set the default assay to the gene symbol-based expression assay
DefaultAssay(seurat.integrated) <- "expr_all"


# -----------------------------------------
# Identify cluster marker genes
# -----------------------------------------

library(Seurat)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Identify positive marker genes for each cluster
cluster_markers <- FindAllMarkers(
  seurat.integrated,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Extract the top 20 marker genes with the highest average log2 fold change
# for each cluster
top20_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

# Save the marker gene list as a CSV file
write.csv(
  top20_per_cluster,
  file = "top20_genes_per_cluster.csv",
  row.names = FALSE
)


# -----------------------------------------
# Load CellMarker 2.0 reference data
# -----------------------------------------

# Load CellMarker 2.0 human and mouse marker gene data
sc_marker <- read_csv(
  "Single cell markers_HumanMouse_CellMarker2.0.csv"
)

# Select human lung marker genes
marker_lung <- sc_marker %>%
  filter(
    tissue_type == "Lung" &
      species == "Human"
  )

# Create a list of marker genes for each cell type
marker_list <- split(
  marker_lung$marker,
  marker_lung$cell_name
)


# -----------------------------------------
# Annotate cell types based on marker genes
# -----------------------------------------

# Compare the top 20 marker genes of each cluster
# with the CellMarker 2.0 reference data
annotated_clusters <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(
    n = 20,
    wt = avg_log2FC
  ) %>%
  summarise(
    annot = {
      genes <- unique(gene)
      
      scores <- sapply(
        marker_list,
        function(x) length(intersect(genes, x))
      )
      
      if (all(scores == 0)) {
        as.character(unique(cluster))
      } else {
        names(which.max(scores))
      }
    }
  )

# Review the annotations using the CellMarker 2.0 reference data
annotated_clusters


# -----------------------------------------
# Add cell type annotations to the Seurat object
# -----------------------------------------

# Map cluster IDs to their assigned cell type labels
seurat.integrated$annotated_clusters <- plyr::mapvalues(
  seurat.integrated$seurat_clusters,
  from = annotated_clusters$cluster,
  to = annotated_clusters$annot
)

# Convert annotations to character strings
seurat.integrated$annotated_clusters <- as.character(
  seurat.integrated$annotated_clusters
)

# Refine selected cell type labels
seurat.integrated$annotated_clusters[
  seurat.integrated$annotated_clusters == "13"
] <- "Ciliated cell"

seurat.integrated$annotated_clusters[
  seurat.integrated$annotated_clusters == "Cancer stem cell"
] <- "Proliferative cell"

# Convert annotations to a factor
seurat.integrated$annotated_clusters <- as.factor(
  seurat.integrated$annotated_clusters
)


# -----------------------------------------
# Visualize cell type annotations
# -----------------------------------------

DimPlot(
  seurat.integrated,
  group.by = "annotated_clusters",
  cols = c(
    "Fibroblast" = "red",
    "Alveolar epithelial cell Type 2" = "green",
    "Endothelial cell" = "lightgreen",
    "Ciliated cell" = "gold",
    "Proliferative cell" = "orange",
    "Alveolar macrophage" = "lightpink",
    "Macrophage" = "pink",
    "Dendritic cell" = "lightblue",
    "Conventional dendritic cell 2(cDC2)" = "blue",
    "Natural killer cell" = "purple",
    "T cell" = "yellow",
    "B cell" = "hotpink",
    "Monocyte" = "cyan"
  ),
  label = FALSE,
  label.size = 3
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )


# -----------------------------------------
# Visualize TGFB1I1 expression by cell type
# -----------------------------------------

# Violin plot grouped by cell type and split by condition
VlnPlot(
  object = seurat.integrated,
  features = "TGFB1I1",
  group.by = "annotated_clusters",
  split.by = "group",
  cols = c(
    "Fibroblast" = "red",
    "Alveolar epithelial cell Type 2" = "green",
    "Endothelial cell" = "lightgreen",
    "Ciliated cell" = "gold",
    "Proliferative cell" = "orange",
    "Alveolar macrophage" = "lightpink",
    "Macrophage" = "pink",
    "Dendritic cell" = "lightblue",
    "Conventional dendritic cell 2(cDC2)" = "blue",
    "Natural killer cell" = "purple",
    "T cell" = "yellow",
    "B cell" = "hotpink",
    "Monocyte" = "cyan"
  ),
  pt.size = 0.1
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 5
    )
  )


# -----------------------------------------
# Compare TGFB1I1 expression between Control and IPF
# -----------------------------------------

# Generate separate violin plots for each condition
plots <- VlnPlot(
  seurat.integrated,
  features = "TGFB1I1",
  group.by = "annotated_clusters",
  split.by = "group",
  pt.size = 0,
  combine = FALSE
)

# Set colors for Control and IPF
plots <- lapply(
  plots,
  function(p) {
    p + scale_fill_manual(
      values = c(
        "Control" = "skyblue",
        "IPF" = "tomato"
      )
    )
  }
)

# Combine the plots
wrap_plots(plots)


# -----------------------------------------
# Stacked violin plots of fibrosis-related genes
# -----------------------------------------

VlnPlot(
  object = seurat.integrated,
  features = c(
    "TGFB1I1",
    "ACTA2",
    "COL1A1",
    "COL1A2",
    "CTHRC1"
  ),
  group.by = "annotated_clusters",
  stack = TRUE
)


# =========================================
# 9. Alternative extraction and gene symbol conversion
# =========================================

# Load required packages
library(Seurat)
library(Matrix)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)


# ------------------------------------------------------------
# Extract all data layers from the RNA assay
# ------------------------------------------------------------

DefaultAssay(seurat.integrated) <- "RNA"

data_layers <- grep(
  "^data",
  Layers(seurat.integrated[["RNA"]]),
  value = TRUE
)

cat("Number of detected data layers:", length(data_layers), "\n")


# ------------------------------------------------------------
# Extract expression data from each layer
# ------------------------------------------------------------

expr_list <- lapply(data_layers, function(lay) {
  GetAssayData(
    seurat.integrated,
    assay = "RNA",
    layer = lay
  )
})


# ------------------------------------------------------------
# Identify the union of genes across all layers
# ------------------------------------------------------------

all_genes <- Reduce(
  union,
  lapply(expr_list, rownames)
)


# ------------------------------------------------------------
# Align gene rows across all layers by zero-padding
# ------------------------------------------------------------

expr_list_padded <- lapply(expr_list, function(mat) {
  missing <- setdiff(all_genes, rownames(mat))
  
  if (length(missing) > 0) {
    missing_mat <- Matrix(
      0,
      nrow = length(missing),
      ncol = ncol(mat),
      dimnames = list(missing, colnames(mat))
    )
    
    mat <- rbind(mat, missing_mat)
  }
  
  mat[all_genes, , drop = FALSE]
})


# ------------------------------------------------------------
# Combine all layers by columns
# ------------------------------------------------------------

expr_all <- do.call(cbind, expr_list_padded)

cat("Number of genes after combining layers:", nrow(expr_all), "\n")


# ------------------------------------------------------------
# Convert Ensembl gene IDs to gene symbols
# ------------------------------------------------------------

ensembl_ids <- rownames(expr_all)

# Retrieve Ensembl ID-to-symbol annotations
symbol_table <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

# Retain the first entry for each Ensembl ID
symbol_table <- symbol_table[
  !duplicated(symbol_table$ENSEMBL),
]

symbol_map <- setNames(
  symbol_table$SYMBOL,
  symbol_table$ENSEMBL
)

# Map Ensembl IDs to gene symbols
symbol_ids <- symbol_map[ensembl_ids]

# Retain Ensembl IDs when no gene symbol is available
symbol_ids[is.na(symbol_ids)] <- ensembl_ids[is.na(symbol_ids)]


# ------------------------------------------------------------
# Replace row names with gene symbols
# ------------------------------------------------------------

rownames(expr_all) <- symbol_ids

cat("Number of genes after gene symbol conversion:", nrow(expr_all), "\n")


# ------------------------------------------------------------
# Add the converted expression matrix as a new assay
# ------------------------------------------------------------

seurat.integrated[["RNA_symbol"]] <- CreateAssayObject(
  counts = expr_all
)

DefaultAssay(seurat.integrated) <- "RNA_symbol"

seurat.integrated

# Save the Seurat object with the RNA_symbol assay
saveRDS(
  seurat.integrated,
  "rds_files/control_IPF_integrated_RNAsymbol.rds"
)

# Reload the saved object
seurat.integrated_RNAsymbol <- readRDS(
  "rds_files/control_IPF_integrated_RNAsymbol.rds"
)


# ------------------------------------------------------------
# Identify cluster marker genes and annotate cell types
# ------------------------------------------------------------

# Identify positive marker genes for each cluster
cluster_markers <- FindAllMarkers(
  seurat.integrated_RNAsymbol,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


# Select the top 20 marker genes per cluster
top20_per_cluster <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

# Save the top marker genes as a CSV file
write.csv(
  top20_per_cluster,
  file = "top20_genes_per_cluster_#2.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# Load CellMarker 2.0 and extract human lung markers
# ------------------------------------------------------------

sc_marker <- readr::read_csv(
  "Single cell markers_HumanMouse_CellMarker2.0.csv"
)

# Restrict the marker database to human lung tissue
marker_lung <- sc_marker %>%
  filter(
    tissue_type == "Lung" &
      species == "Human"
  )

# Generate a list of marker genes for each cell type
marker_list <- split(
  marker_lung$marker,
  marker_lung$cell_name
)


# ------------------------------------------------------------
# Annotate clusters based on CellMarker 2.0
# ------------------------------------------------------------

# Compare the top 20 marker genes of each cluster
# with the cell-type marker lists
annotated_clusters <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(
    n = 20,
    wt = avg_log2FC
  ) %>%
  summarise(
    annot = {
      genes <- unique(gene)
      
      scores <- sapply(
        marker_list,
        function(x) length(intersect(genes, x))
      )
      
      if (all(scores == 0)) {
        as.character(unique(cluster))
      } else {
        names(which.max(scores))
      }
    }
  )


# Cluster annotations obtained from CellMarker 2.0:
#
# Cluster 0  : Macrophage
# Cluster 1  : Macrophage
# Cluster 2  : Dendritic cell
# Cluster 3  : Alveolar macrophage
# Cluster 4  : Monocyte
# Cluster 5  : Macrophage
# Cluster 6  : Monocyte
# Cluster 7  : Natural killer cell
# Cluster 8  : T cell
# Cluster 9  : Alveolar epithelial cell Type 2
# Cluster 10 : Fibroblast
# Cluster 11 : Endothelial cell
# Cluster 12 : B cell
# Cluster 13 : Ciliated cell
# Cluster 14 : Conventional dendritic cell 2 (cDC2)
# Cluster 15 : Proliferative cell
# Cluster 16 : B cell
# Cluster 17 : T cell
# Cluster 18 : Proliferative cell


# ------------------------------------------------------------
# Add cell-type annotations to the Seurat object
# ------------------------------------------------------------

# Assign cell-type labels based on cluster IDs
seurat.integrated_RNAsymbol$annotated_clusters <- plyr::mapvalues(
  seurat.integrated_RNAsymbol$seurat_clusters,
  from = annotated_clusters$cluster,
  to = annotated_clusters$annot
)

# Convert annotations to character vectors before relabeling
seurat.integrated_RNAsymbol$annotated_clusters <- as.character(
  seurat.integrated_RNAsymbol$annotated_clusters
)

# Rename selected cell-type labels
seurat.integrated_RNAsymbol$annotated_clusters[
  seurat.integrated_RNAsymbol$annotated_clusters == "13"
] <- "Ciliated cell"

seurat.integrated_RNAsymbol$annotated_clusters[
  seurat.integrated_RNAsymbol$annotated_clusters == "Cancer stem cell"
] <- "Proliferative cell"

seurat.integrated_RNAsymbol$annotated_clusters[
  seurat.integrated_RNAsymbol$annotated_clusters == "Follicular B cell"
] <- "B cell"

seurat.integrated_RNAsymbol$annotated_clusters[
  seurat.integrated_RNAsymbol$annotated_clusters == "Conventional dendritic cell 2(cDC2)"
] <- "Conventional dendritic cell 2"

# Convert annotations to a factor
seurat.integrated_RNAsymbol$annotated_clusters <- as.factor(
  seurat.integrated_RNAsymbol$annotated_clusters
)


# Save the annotated Seurat object
saveRDS(
  seurat.integrated_RNAsymbol,
  "rds_files/control_IPF_integrated_RNAsymbol.rds"
)

# Reload the annotated object
seurat.integrated_RNAsymbol <- readRDS(
  "rds_files/control_IPF_integrated_RNAsymbol.rds"
)


# ============================================================
# 10. Visualization of cell-type annotations and gene expression
# ============================================================

# ------------------------------------------------------------
# UMAP visualization by cell type
# ------------------------------------------------------------

# Define colors for each cell type
cell_colors <- c(
  "Fibroblast" = "red",
  "Alveolar epithelial cell Type 2" = "green",
  "Endothelial cell" = "lightgreen",
  "Ciliated cell" = "gold",
  "Proliferative cell" = "orange",
  "Dendritic cell" = "lightblue",
  "Conventional dendritic cell 2" = "skyblue",
  "Alveolar macrophage" = "pink",
  "Macrophage" = "hotpink",
  "Monocyte" = "cyan",
  "Natural killer cell" = "purple",
  "T cell" = "yellow",
  "B cell" = "magenta"
)

# Specify the order of cell types in the legend
cell_order <- c(
  # Fibrosis-related cells
  "Fibroblast",
  "Proliferative cell",
  
  # Epithelial cells
  "Alveolar epithelial cell Type 2",
  "Ciliated cell",
  
  # Vascular cells
  "Endothelial cell",
  
  # Immune cells
  "Monocyte",
  "Macrophage",
  "Alveolar macrophage",
  "Dendritic cell",
  "Conventional dendritic cell 2",
  "Natural killer cell",
  "T cell",
  "B cell"
)

# Generate a UMAP plot colored by cell type
DimPlot(
  seurat.integrated_RNAsymbol,
  group.by = "annotated_clusters",
  cols = cell_colors,
  order = cell_order,
  label = FALSE,
  label.size = 10
) +
  scale_color_manual(
    values = cell_colors[cell_order],
    breaks = cell_order
  ) +
  theme_classic() +
  ggtitle("UMAP of Lung Cells by Cell Type") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20)
  ) +
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  )


# ------------------------------------------------------------
# Violin plot of TGFB1I1 expression by cell type
# ------------------------------------------------------------

# Specify the order of cell types on the x-axis
seurat.integrated_RNAsymbol$annotated_clusters <- factor(
  seurat.integrated_RNAsymbol$annotated_clusters,
  levels = c(
    "Fibroblast",
    "Proliferative cell",
    
    # Epithelial cells
    "Alveolar epithelial cell Type 2",
    "Ciliated cell",
    
    # Vascular cells
    "Endothelial cell",
    
    # Immune cells
    "Monocyte",
    "Macrophage",
    "Alveolar macrophage",
    "Dendritic cell",
    "Conventional dendritic cell 2",
    "Natural killer cell",
    "T cell",
    "B cell"
  )
)

# Generate violin plots split by experimental group
plots <- VlnPlot(
  seurat.integrated_RNAsymbol,
  features = "TGFB1I1",
  group.by = "annotated_clusters",
  split.by = "group",
  pt.size = 0,
  combine = FALSE
)

# Set colors for Control and IPF
plots <- lapply(plots, function(p) {
  p +
    scale_fill_manual(
      values = c(
        "Control" = "skyblue",
        "IPF" = "tomato"
      )
    ) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 15)
    )
})

# Display the violin plots
wrap_plots(plots)


# ------------------------------------------------------------
# Stacked violin plots of fibrosis-related gene expression
# ------------------------------------------------------------

# Create an abbreviated cell-type annotation column
seurat.integrated_RNAsymbol$annotatedclus_abbr <-
  seurat.integrated_RNAsymbol$annotated_clusters

# Replace cell-type names with abbreviations
seurat.integrated_RNAsymbol$annotatedclus_abbr <- plyr::revalue(
  seurat.integrated_RNAsymbol$annotatedclus_abbr,
  c(
    "Fibroblast" = "Fibro",
    "Proliferative cell" = "Prolif.",
    
    # Epithelial cells
    "Alveolar epithelial cell Type 2" = "AT2",
    "Ciliated cell" = "Ciliated",
    
    # Vascular cells
    "Endothelial cell" = "Endo",
    
    # Immune cells
    "Monocyte" = "Mono",
    "Macrophage" = "MΦ",
    "Alveolar macrophage" = "AM",
    "Dendritic cell" = "DC",
    "Conventional dendritic cell 2" = "cDC2",
    "Natural killer cell" = "NK",
    "T cell" = "T",
    "B cell" = "B"
  )
)

# Specify the order of abbreviated cell types
seurat.integrated_RNAsymbol$annotatedclus_abbr <- factor(
  seurat.integrated_RNAsymbol$annotatedclus_abbr,
  levels = c(
    "B",
    "T",
    "NK",
    "cDC2",
    "DC",
    "AM",
    "MΦ",
    "Mono",
    "Endo",
    "Ciliated",
    "AT2",
    "Prolif.",
    "Fibro"
  )
)

# Set the order of experimental groups
seurat.integrated_RNAsymbol$group <- factor(
  seurat.integrated_RNAsymbol$group,
  levels = c("IPF", "Control")
)

# Generate stacked violin plots
p <- VlnPlot(
  object = seurat.integrated_RNAsymbol,
  features = c(
    "TGFB1I1",
    "ACTA2",
    "COL1A1",
    "COL1A2",
    "CTHRC1"
  ),
  group.by = "annotatedclus_abbr",
  split.by = "group",
  pt.size = 0,
  cols = c("tomato", "skyblue"),
  stack = TRUE,
  combine = TRUE
)

# Adjust plot appearance
p +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 15),
    strip.text.x = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  )


# =========================================
# 11. Fibroblast cluster analysis
# =========================================

# -----------------------------------------
# Extract the fibroblast cluster
# -----------------------------------------

# Set cell identities to annotated cell types
Idents(seurat.integrated_RNAsymbol) <- "annotated_clusters"

# Extract fibroblasts
fibro <- subset(
  seurat.integrated_RNAsymbol,
  idents = "Fibroblast"
)

# Check the number of cells in each group
table(fibro$group)
# IPF: 4548; Control: 1063


# -----------------------------------------
# Data normalization and scaling
# -----------------------------------------

# Identify highly variable features
fibro <- FindVariableFeatures(
  fibro,
  selection.method = "vst",
  nfeatures = 2000
)

# Scale the data
fibro <- ScaleData(fibro)


# -----------------------------------------
# Dimensionality reduction by PCA
# -----------------------------------------

fibro <- RunPCA(
  fibro,
  features = VariableFeatures(fibro)
)

ElbowPlot(fibro)


# -----------------------------------------
# Neighbor detection and clustering
# -----------------------------------------

# Identify neighboring cells using the first 15 PCs
fibro <- FindNeighbors(
  fibro,
  dims = 1:15
)

# Perform graph-based clustering
fibro <- FindClusters(
  fibro,
  resolution = 0.2
)


# -----------------------------------------
# UMAP visualization
# -----------------------------------------

fibro <- RunUMAP(
  fibro,
  dims = 1:15
)

DimPlot(fibro, label = TRUE)

FeaturePlot(
  fibro,
  features = "TGFB1I1"
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )

# Set the group factor levels
fibro$group <- factor(
  fibro$group,
  levels = c("Control", "IPF")
)


# -----------------------------------------
# Identification of marker genes
# -----------------------------------------

fibro.markers <- FindAllMarkers(
  fibro,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(fibro.markers)

# Extract the top 20 marker genes for each cluster
top20_per_fibrocluster <- fibro.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

# Save the marker gene list
write.csv(
  top20_per_fibrocluster,
  file = "top20_genes_per_fibrocluster.csv",
  row.names = FALSE
)


# -----------------------------------------
# Cell type annotation using CellMarker 2.0
# -----------------------------------------

# Load CellMarker 2.0 marker gene data
sc_marker <- readr::read_csv(
  "Single cell markers_HumanMouse_CellMarker2.0.csv"
)

# Select human lung cell markers
marker_lung <- sc_marker %>%
  filter(
    tissue_type == "Lung" &
      species == "Human"
  )

# Generate marker gene lists for each cell type
marker_list <- split(
  marker_lung$marker,
  marker_lung$cell_name
)

# Annotate clusters based on the overlap between
# cluster marker genes and CellMarker 2.0 marker genes
fibro.annotated_clusters <- fibro.markers %>%
  group_by(cluster) %>%
  top_n(
    n = 20,
    wt = avg_log2FC
  ) %>%
  summarise(
    annot = {
      genes <- unique(gene)
      scores <- sapply(
        marker_list,
        function(x) length(intersect(genes, x))
      )
      
      if (all(scores == 0)) {
        as.character(unique(cluster))
      } else {
        names(which.max(scores))
      }
    }
  )

# Display the annotation results
fibro.annotated_clusters

# Add cell type annotations to the Seurat object
fibro$fibro.annotated_clusters <- plyr::mapvalues(
  fibro$seurat_clusters,
  from = fibro.annotated_clusters$cluster,
  to = fibro.annotated_clusters$annot
)

# Check the number of cells in each annotated cluster
table(fibro$fibro.annotated_clusters)

# Rename cell type labels
fibro$fibro.annotated_clusters <- as.character(
  fibro$fibro.annotated_clusters
)

fibro$fibro.annotated_clusters[
  fibro$fibro.annotated_clusters == "5"
] <- "Unknown"

fibro$fibro.annotated_clusters[
  fibro$fibro.annotated_clusters == "Alveolar cell Type 1"
] <- "Alveolar Type 1-like cell"

fibro$fibro.annotated_clusters[
  fibro$fibro.annotated_clusters == "Alveolar macrophage"
] <- "Alveolar macrophage-like cell"

fibro$fibro.annotated_clusters[
  fibro$fibro.annotated_clusters == "Epithelial cell"
] <- "Epithelial-like cell"

fibro$fibro.annotated_clusters[
  fibro$fibro.annotated_clusters == "Smooth muscle cell"
] <- "Smooth muscle-like cell"

fibro$fibro.annotated_clusters <- as.factor(
  fibro$fibro.annotated_clusters
)


# -----------------------------------------
# Save and reload the annotated fibroblast object
# -----------------------------------------

saveRDS(
  fibro,
  "rds_files/control_IPF_fibro.rds"
)

fibro <- readRDS(
  "rds_files/control_IPF_fibro.rds"
)


# -----------------------------------------
# UMAP visualization by annotated cell type
# -----------------------------------------

cell_colors <- c(
  "Myofibroblast" = "red",
  "Alveolar Type 1-like cell" = "green",
  "Epithelial-like cell" = "gold",
  "Alveolar macrophage-like cell" = "pink",
  "Smooth muscle-like cell" = "darkred",
  "Unknown" = "gray"
)

cell_order <- c(
  "Unknown",
  "Alveolar macrophage-like cell",
  "Smooth muscle-like cell",
  "Epithelial-like cell",
  "Alveolar Type 1-like cell",
  "Myofibroblast"
)

DimPlot(
  fibro,
  group.by = "fibro.annotated_clusters",
  cols = cell_colors,
  order = cell_order,
  label = FALSE,
  label.size = 10
) +
  theme_classic() +
  ggtitle("UMAP of Lung Fibroblasts by Cell Type") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20)
  ) +
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  )


# -----------------------------------------
# Violin plot of TGFB1I1 expression
# -----------------------------------------

# Create a separate column for abbreviated cell type labels
fibro$fibro.abbr <- fibro$fibro.annotated_clusters

# Abbreviate cell type labels
fibro$fibro.abbr <- plyr::revalue(
  fibro$fibro.abbr,
  c(
    "Myofibroblast" = "Myo",
    "Alveolar Type 1-like cell" = "AT1-like",
    "Epithelial-like cell" = "Epi-like",
    "Smooth muscle-like cell" = "SMC-like",
    "Alveolar macrophage-like cell" = "AM-like",
    "Unknown" = "Unknown"
  )
)

# Set the display order of cell types
fibro$fibro.abbr <- factor(
  fibro$fibro.abbr,
  levels = c(
    "Myo",
    "AT1-like",
    "Epi-like",
    "SMC-like",
    "AM-like",
    "Unknown"
  )
)

# Generate violin plots of TGFB1I1 expression
plots <- VlnPlot(
  fibro,
  features = "TGFB1I1",
  group.by = "fibro.abbr",
  split.by = "group",
  pt.size = 0,
  combine = FALSE
)

# Set colors for Control and IPF groups
plots <- lapply(
  plots,
  function(p) {
    p +
      scale_fill_manual(
        values = c(
          "Control" = "skyblue",
          "IPF" = "tomato"
        )
      ) +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)
      )
  }
)

# Display violin plots
patchwork::wrap_plots(plots)


# -----------------------------------------
# Changes in cluster composition
# (Stacked bar plot)
# -----------------------------------------

library(dplyr)
library(ggplot2)

# Calculate the number of cells in each cluster
# for each condition
cluster_counts <- fibro@meta.data %>%
  group_by(group, fibro.abbr) %>%
  summarise(
    CellCount = n()
  ) %>%
  ungroup()

# Calculate the proportion of each cluster
# within each condition
cluster_props <- cluster_counts %>%
  group_by(group) %>%
  mutate(
    Proportion = CellCount / sum(CellCount)
  )

# Define colors for each fibroblast cell type
cell_colors.abbr <- c(
  "Myo" = "red",
  "AT1-like" = "green",
  "Epi-like" = "gold",
  "AM-like" = "pink",
  "SMC-like" = "darkred",
  "Unknown" = "gray"
)

# Generate a stacked bar plot
ggplot(
  cluster_props,
  aes(
    x = group,
    y = Proportion,
    fill = fibro.abbr
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_fill_manual(
    values = cell_colors.abbr
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    y = "Cell proportion",
    x = "Condition"
  ) +
  theme_classic() +
  theme(
    # Axis title size
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    
    # Axis tick label size
    axis.text.x = element_text(
      size = 18,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 18),
    
    # Legend formatting
    legend.text = element_text(size = 18),
    legend.title = element_blank()
  )


# -----------------------------------------
# Proportion of TGFB1I1-positive cells
# -----------------------------------------

# Define TGFB1I1-positive cells based on an expression threshold
fibro$TGFB1I1_pos <- fibro@assays$RNA_symbol@data[
  "TGFB1I1",
] > 0.25

# Calculate the proportion of TGFB1I1-positive cells
# in each annotated cluster and condition
tgfb1i1_props <- fibro@meta.data %>%
  group_by(
    group,
    fibro.annotated_clusters
  ) %>%
  summarise(
    n_total = n(),
    n_pos = sum(TGFB1I1_pos)
  ) %>%
  mutate(
    Proportion = n_pos / n_total
  )


# -----------------------------------------
# FeaturePlot of TGFB1I1 expression
# -----------------------------------------

# Downsample IPF cells to match the number of Control cells
set.seed(123)

cells.control <- colnames(fibro)[
  fibro$group == "Control"
]

cells.ipf <- colnames(fibro)[
  fibro$group == "IPF"
]

downsample.ipf_fibro <- sample(
  cells.ipf,
  size = length(cells.control)
)

# Extract Control cells and downsampled IPF cells
seurat.down_fibroblast <- subset(
  fibro,
  cells = c(
    cells.control,
    downsample.ipf_fibro
  )
)

# Check the number of cells in each group
table(seurat.down_fibroblast$group)
# Control: 1063; IPF: 1063


# -----------------------------------------
# Save and reload the downsampled object
# -----------------------------------------

saveRDS(
  seurat.down_fibroblast,
  "rds_files/seurat.down_fibroblast.rds"
)

seurat.down_fibroblast <- readRDS(
  "rds_files/seurat.down_fibroblast.rds"
)


# -----------------------------------------
# TGFB1I1 FeaturePlot
# -----------------------------------------

# Set the default assay to RNA_symbol
DefaultAssay(seurat.down_fibroblast) <- "RNA_symbol"

# Set the group factor levels
seurat.down_fibroblast$group <- factor(
  seurat.down_fibroblast$group,
  levels = c("Control", "IPF")
)

# Generate split FeaturePlots
p <- FeaturePlot(
  seurat.down_fibroblast,
  features = "TGFB1I1",
  split.by = "group",
  pt.size = 0.5,
  min.cutoff = 0,
  max.cutoff = 2
)

# Display the plots
p +
  plot_annotation(
    title = "TGFB1I1",
    caption = "Control vs IPF"
  ) &
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  ) &
  theme(
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20)
  )
  

# =========================================
# 12. Comparison of TGFB1I1-positive cell proportions
# and expression levels in myofibroblasts
# =========================================

library(dplyr)
library(ggplot2)
library(tidyr)

# -----------------------------------------
# Extract myofibroblasts
# -----------------------------------------

myo <- subset(
  fibro,
  subset = fibro.annotated_clusters == "Myofibroblast"
)


# -----------------------------------------
# Extract TGFB1I1 expression and sample metadata
# -----------------------------------------

myo_df <- FetchData(
  myo,
  vars = c("TGFB1I1", "group", "sample_id")
) %>%
  rename(
    condition = group,
    sample = sample_id
  ) %>%
  mutate(
    TGFB1I1 = replace_na(TGFB1I1, 0)
  )


# -----------------------------------------
# Calculate TGFB1I1-positive cell proportions
# and mean expression levels for each sample
# -----------------------------------------

myo_summary <- myo_df %>%
  group_by(sample, condition) %>%
  summarise(
    prop_positive = mean(TGFB1I1 > 0),
    mean_expr = ifelse(
      all(TGFB1I1 == 0),
      0,
      mean(TGFB1I1[TGFB1I1 > 0])
    ),
    .groups = "drop"
  )


# -----------------------------------------
# Calculate group means and standard errors
# -----------------------------------------

myo_summary_stats <- myo_summary %>%
  group_by(condition) %>%
  summarise(
    prop_mean = mean(prop_positive),
    prop_sem = sd(prop_positive) / sqrt(n()),
    expr_mean = mean(mean_expr),
    expr_sem = sd(mean_expr) / sqrt(n()),
    .groups = "drop"
  )


# -----------------------------------------
# Statistical analysis
# -----------------------------------------

# Test for normality using the Shapiro-Wilk test

shapiro.test(
  myo_summary$prop_positive[
    myo_summary$condition == "Control"
  ]
)
# W = 0.62093, p-value = 2.25e-06

shapiro.test(
  myo_summary$prop_positive[
    myo_summary$condition == "IPF"
  ]
)
# W = 0.94721, p-value = 0.1423

shapiro.test(
  myo_summary$mean_expr[
    myo_summary$condition == "Control"
  ]
)
# W = 0.81367, p-value = 0.0008254

shapiro.test(
  myo_summary$mean_expr[
    myo_summary$condition == "IPF"
  ]
)
# W = 0.77472, p-value = 2.346e-05

# Non-normal distributions were observed;
# therefore, the Wilcoxon rank-sum test was used.

# Compare the proportions of TGFB1I1-positive cells
prop_test <- wilcox.test(
  prop_positive ~ condition,
  data = myo_summary
)

prop_test$p.value
# 0.00800582

# Compare mean TGFB1I1 expression levels in positive cells
expr_test <- wilcox.test(
  mean_expr ~ condition,
  data = myo_summary
)

expr_test$p.value
# 0.6384173


# -----------------------------------------
# Visualization
# -----------------------------------------

library(gridExtra)

# Recalculate Wilcoxon test results for plotting
prop_test <- wilcox.test(
  prop_positive ~ condition,
  data = myo_summary
)

expr_test <- wilcox.test(
  mean_expr ~ condition,
  data = myo_summary
)


# -----------------------------------------
# Proportion of TGFB1I1-positive cells
# -----------------------------------------

p1 <- ggplot(
  myo_summary,
  aes(
    x = condition,
    y = prop_positive * 100,
    fill = condition
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.4
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    height = 0,
    size = 1
  ) +
  labs(
    y = "TGFB1I1+ cells (%)",
    x = ""
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(myo_summary$prop_positive * 100) * 1.05,
    label = paste0(
      "p = ",
      signif(prop_test$p.value, 3)
    )
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "skyblue",
      "IPF" = "tomato"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"
  )


# -----------------------------------------
# Mean TGFB1I1 expression in positive cells
# -----------------------------------------

p2 <- ggplot(
  myo_summary,
  aes(
    x = condition,
    y = mean_expr,
    fill = condition
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.4
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    height = 0,
    size = 1
  ) +
  labs(
    y = "Mean TGFB1I1 expression (positive cells)",
    x = ""
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(myo_summary$mean_expr) * 1.05,
    label = paste0(
      "p = ",
      signif(expr_test$p.value, 3)
    )
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "skyblue",
      "IPF" = "tomato"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"
  )


# -----------------------------------------
# Display plots side by side
# -----------------------------------------

gridExtra::grid.arrange(
  p1,
  p2,
  ncol = 2
)


# =========================================
# 13. Comparison of TGFB1I1-positive cell proportions
# and expression levels in fibroblasts
# =========================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# Extract TGFB1I1 expression, group, and sample information
fibro_df <- FetchData(
  fibro,
  vars = c("TGFB1I1", "group", "sample_id")
) %>%
  rename(
    condition = group,
    sample = sample_id
  ) %>%
  mutate(
    TGFB1I1 = replace_na(TGFB1I1, 0)
  )

# Calculate the proportion of TGFB1I1-positive cells
# and mean expression among positive cells for each sample
fibro_summary <- fibro_df %>%
  group_by(sample, condition) %>%
  summarise(
    prop_positive = mean(TGFB1I1 > 0),
    mean_expr = ifelse(
      all(TGFB1I1 == 0),
      0,
      mean(TGFB1I1[TGFB1I1 > 0])
    ),
    .groups = "drop"
  )

# Calculate mean and SEM for each condition
fibro_summary_stats <- fibro_summary %>%
  group_by(condition) %>%
  summarise(
    prop_mean = mean(prop_positive),
    prop_sem = sd(prop_positive) / sqrt(n()),
    expr_mean = mean(mean_expr),
    expr_sem = sd(mean_expr) / sqrt(n()),
    .groups = "drop"
  )


# -----------------------------------------
# Statistical analysis
# -----------------------------------------

# Test for normality (Shapiro-Wilk test)
shapiro.test(
  fibro_summary$prop_positive[fibro_summary$condition == "Control"]
)
# W = 0.93434, p-value = 0.08837

shapiro.test(
  fibro_summary$prop_positive[fibro_summary$condition == "IPF"]
)
# W = 0.9707, p-value = 0.5191

shapiro.test(
  fibro_summary$mean_expr[fibro_summary$condition == "Control"]
)
# W = 0.88221, p-value = 0.005364

shapiro.test(
  fibro_summary$mean_expr[fibro_summary$condition == "IPF"]
)
# W = 0.89513, p-value = 0.004667

# Compare groups using the Wilcoxon rank-sum test

# Comparison of TGFB1I1-positive cell proportions
prop_test <- wilcox.test(
  prop_positive ~ condition,
  data = fibro_summary
)
prop_test$p.value
# 0.01003332

# Comparison of mean TGFB1I1 expression in positive cells
expr_test <- wilcox.test(
  mean_expr ~ condition,
  data = fibro_summary
)
expr_test$p.value
# 0.9696215


# -----------------------------------------
# Visualization
# -----------------------------------------

# TGFB1I1-positive cell proportions
p1 <- ggplot(
  fibro_summary,
  aes(x = condition, y = prop_positive * 100, fill = condition)
) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  labs(
    y = "TGFB1I1+ cells (%)",
    x = ""
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(fibro_summary$prop_positive * 100) * 1.05,
    label = paste0("p = ", signif(prop_test$p.value, 4))
  ) +
  scale_fill_manual(
    values = c("Control" = "skyblue", "IPF" = "tomato")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"
  )


# Mean TGFB1I1 expression in positive cells
p2 <- ggplot(
  fibro_summary,
  aes(x = condition, y = mean_expr, fill = condition)
) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  labs(
    y = "Mean TGFB1I1 expression (positive cells)",
    x = ""
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(fibro_summary$mean_expr) * 1.05,
    label = paste0("p = ", signif(expr_test$p.value, 4))
  ) +
  scale_fill_manual(
    values = c("Control" = "skyblue", "IPF" = "tomato")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"
  )


# Display plots side by side
gridExtra::grid.arrange(p1, p2, ncol = 2)




# ============================================================
# 14. Analysis of ACTA2, COL1A1, COL1A2, and CTHRC1 Expression
#     in Fibroblasts
# ============================================================

#--------------------------
# 14.1. Proportions of Positive Cells and Gene Expression
#       Levels in Fibroblasts
#--------------------------

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)


#--------------------------
# Retrieve Gene Expression Data
#--------------------------

# Retrieve expression values of TGFB1I1, ACTA2, COL1A1,
# COL1A2, and CTHRC1, together with group and sample information.

fibro_df <- FetchData(
  fibro,
  vars = c(
    "TGFB1I1", "ACTA2", "COL1A1", "COL1A2", "CTHRC1",
    "group", "sample_id"
  )
) %>%
  mutate(
    TGFB1I1 = replace_na(TGFB1I1, 0),
    ACTA2   = replace_na(ACTA2, 0),
    COL1A1  = replace_na(COL1A1, 0),
    COL1A2  = replace_na(COL1A2, 0),
    CTHRC1  = replace_na(CTHRC1, 0)
  ) %>%
  rename(
    condition = group,
    sample = sample_id
  )


#--------------------------
# Calculate the Proportion of Positive Cells
# and Mean Expression Levels for Each Sample
#--------------------------

fibro_summary <- fibro_df %>%
  group_by(sample, condition) %>%
  summarise(
    
    # Proportion of positive cells
    prop_TGFB1I1 = mean(TGFB1I1 > 0),
    prop_ACTA2   = mean(ACTA2 > 0),
    prop_COL1A1  = mean(COL1A1 > 0),
    prop_COL1A2  = mean(COL1A2 > 0),
    prop_CTHRC1  = mean(CTHRC1 > 0),
    
    # Mean expression levels among positive cells
    mean_TGFB1I1 = ifelse(
      all(TGFB1I1 == 0), 0,
      mean(TGFB1I1[TGFB1I1 > 0])
    ),
    
    mean_ACTA2 = ifelse(
      all(ACTA2 == 0), 0,
      mean(ACTA2[ACTA2 > 0])
    ),
    
    mean_COL1A1 = ifelse(
      all(COL1A1 == 0), 0,
      mean(COL1A1[COL1A1 > 0])
    ),
    
    mean_COL1A2 = ifelse(
      all(COL1A2 == 0), 0,
      mean(COL1A2[COL1A2 > 0])
    ),
    
    mean_CTHRC1 = ifelse(
      all(CTHRC1 == 0), 0,
      mean(CTHRC1[CTHRC1 > 0])
    ),
    
    .groups = "drop"
  )


#--------------------------
# Calculate Group Means and Standard Errors
# of the Mean (SEM)
#--------------------------

fibro_summary_stats <- fibro_summary %>%
  group_by(condition) %>%
  summarise(
    
    # Mean proportions of positive cells
    prop_TGFB1I1mean = mean(prop_TGFB1I1),
    prop_ACTA2mean   = mean(prop_ACTA2),
    prop_COL1A1mean  = mean(prop_COL1A1),
    prop_COL1A2mean  = mean(prop_COL1A2),
    prop_CTHRC1mean  = mean(prop_CTHRC1),
    
    # SEMs of the proportions of positive cells
    prop_TGFB1I1_sem = sd(prop_TGFB1I1) / sqrt(n()),
    prop_ACTA2_sem   = sd(prop_ACTA2) / sqrt(n()),
    prop_COL1A1_sem  = sd(prop_COL1A1) / sqrt(n()),
    prop_COL1A2_sem  = sd(prop_COL1A2) / sqrt(n()),
    prop_CTHRC1_sem  = sd(prop_CTHRC1) / sqrt(n()),
    
    # Mean expression levels among positive cells
    expr_TGFB1I1mean = mean(mean_TGFB1I1),
    expr_ACTA2mean   = mean(mean_ACTA2),
    expr_COL1A1mean  = mean(mean_COL1A1),
    expr_COL1A2mean  = mean(mean_COL1A2),
    expr_CTHRC1mean  = mean(mean_CTHRC1),
    
    # SEMs of mean expression levels among positive cells
    expr_TGFB1I1sem = sd(mean_TGFB1I1) / sqrt(n()),
    expr_ACTA2sem   = sd(mean_ACTA2) / sqrt(n()),
    expr_COL1A1sem  = sd(mean_COL1A1) / sqrt(n()),
    expr_COL1A2sem  = sd(mean_COL1A2) / sqrt(n()),
    expr_CTHRC1sem  = sd(mean_CTHRC1) / sqrt(n()),
    
    .groups = "drop"
  )


#--------------------------
# Statistical Analysis
#--------------------------

# Test for normality using the Shapiro-Wilk test.

# TGFB1I1
shapiro.test(fibro_summary$prop_TGFB1I1[fibro_summary$condition == "Control"])
# W = 0.93434, p-value = 0.08837

shapiro.test(fibro_summary$prop_TGFB1I1[fibro_summary$condition == "IPF"])
# W = 0.9707, p-value = 0.5191

shapiro.test(fibro_summary$mean_TGFB1I1[fibro_summary$condition == "Control"])
# W = 0.88221, p-value = 0.005364

shapiro.test(fibro_summary$mean_TGFB1I1[fibro_summary$condition == "IPF"])
# W = 0.89513, p-value = 0.004667


# ACTA2
shapiro.test(fibro_summary$prop_ACTA2[fibro_summary$condition == "Control"])
# W = 0.69227, p-value = 3.032e-06

shapiro.test(fibro_summary$prop_ACTA2[fibro_summary$condition == "IPF"])
# W = 0.85964, p-value = 0.0006758

shapiro.test(fibro_summary$mean_ACTA2[fibro_summary$condition == "Control"])
# W = 0.77116, p-value = 4.508e-05

shapiro.test(fibro_summary$mean_ACTA2[fibro_summary$condition == "IPF"])
# W = 0.80791, p-value = 5.759e-05


# COL1A1
shapiro.test(fibro_summary$prop_COL1A1[fibro_summary$condition == "Control"])
# W = 0.96275, p-value = 0.426

shapiro.test(fibro_summary$prop_COL1A1[fibro_summary$condition == "IPF"])
# W = 0.98239, p-value = 0.8654

shapiro.test(fibro_summary$mean_COL1A1[fibro_summary$condition == "Control"])
# W = 0.84966, p-value = 0.001134

shapiro.test(fibro_summary$mean_COL1A1[fibro_summary$condition == "IPF"])
# W = 0.97543, p-value = 0.6603


# COL1A2
shapiro.test(fibro_summary$prop_COL1A2[fibro_summary$condition == "Control"])
# W = 0.90965, p-value = 0.0224

shapiro.test(fibro_summary$prop_COL1A2[fibro_summary$condition == "IPF"])
# W = 0.83168, p-value = 0.0001705

shapiro.test(fibro_summary$mean_COL1A2[fibro_summary$condition == "Control"])
# W = 0.73562, p-value = 1.264e-05

shapiro.test(fibro_summary$mean_COL1A2[fibro_summary$condition == "IPF"])
# W = 0.69267, p-value = 6.766e-07


# CTHRC1
shapiro.test(fibro_summary$prop_CTHRC1[fibro_summary$condition == "Control"])
# W = 0.793, p-value = 0.0001038

shapiro.test(fibro_summary$prop_CTHRC1[fibro_summary$condition == "IPF"])
# W = 0.91285, p-value = 0.01336

shapiro.test(fibro_summary$mean_CTHRC1[fibro_summary$condition == "Control"])
# W = 0.8623, p-value = 0.002038

shapiro.test(fibro_summary$mean_CTHRC1[fibro_summary$condition == "IPF"])
# W = 0.84852, p-value = 0.0003855


# Based on the normality test results, the Wilcoxon
# rank-sum test was used consistently for group comparisons.


#--------------------------
# Comparison of the Proportions of Positive Cells
# Between Control and IPF Groups
#--------------------------

prop_tests <- list(
  
  TGFB1I1 = wilcox.test(
    prop_TGFB1I1 ~ condition,
    data = fibro_summary
  ),
  
  ACTA2 = wilcox.test(
    prop_ACTA2 ~ condition,
    data = fibro_summary
  ),
  
  COL1A1 = wilcox.test(
    prop_COL1A1 ~ condition,
    data = fibro_summary
  ),
  
  COL1A2 = wilcox.test(
    prop_COL1A2 ~ condition,
    data = fibro_summary
  ),
  
  CTHRC1 = wilcox.test(
    prop_CTHRC1 ~ condition,
    data = fibro_summary
  )
)

# Check p-values
lapply(prop_tests, function(x) x$p.value)

# TGFB1I1: 0.01003332
# ACTA2:   0.001263208
# COL1A1:  6.50028e-07
# COL1A2:  5.522939e-05
# CTHRC1:  4.804139e-05


#--------------------------
# Comparison of Mean Expression Levels Among
# Positive Cells Between Control and IPF Groups
#--------------------------

expr_tests <- list(
  
  TGFB1I1 = wilcox.test(
    mean_TGFB1I1 ~ condition,
    data = fibro_summary
  ),
  
  ACTA2 = wilcox.test(
    mean_ACTA2 ~ condition,
    data = fibro_summary
  ),
  
  COL1A1 = wilcox.test(
    mean_COL1A1 ~ condition,
    data = fibro_summary
  ),
  
  COL1A2 = wilcox.test(
    mean_COL1A2 ~ condition,
    data = fibro_summary
  ),
  
  CTHRC1 = wilcox.test(
    mean_CTHRC1 ~ condition,
    data = fibro_summary
  )
)

# Check p-values
lapply(expr_tests, function(x) x$p.value)


#--------------------------
# Plotting
#--------------------------

# Function to generate plots of the proportions
# of positive cells.
plot_proportion <- function(data, variable, gene, test) {
  
  ggplot(
    data,
    aes(
      x = condition,
      y = .data[[variable]] * 100,
      fill = condition
    )
  ) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(
      width = 0.1,
      fill = "white",
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.1,
      height = 0,
      size = 1
    ) +
    labs(
      y = paste0(gene, "+ cells (%)"),
      x = ""
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(data[[variable]] * 100) * 1.05,
      label = paste0(
        "p = ",
        signif(test$p.value, 4)
      )
    ) +
    scale_fill_manual(
      values = c(
        "Control" = "skyblue",
        "IPF" = "tomato"
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    )
}


# Function to generate plots of mean expression levels
# among positive cells.
plot_expression <- function(data, variable, gene, test) {
  
  ggplot(
    data,
    aes(
      x = condition,
      y = .data[[variable]],
      fill = condition
    )
  ) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(
      width = 0.1,
      fill = "white",
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.1,
      height = 0,
      size = 1
    ) +
    labs(
      y = paste0(
        "Mean ", gene, " expression (positive cells)"
      ),
      x = ""
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(data[[variable]]) * 1.05,
      label = paste0(
        "p = ",
        signif(test$p.value, 4)
      )
    ) +
    scale_fill_manual(
      values = c(
        "Control" = "skyblue",
        "IPF" = "tomato"
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    )
}


#--------------------------
# Proportions of Positive Cells
#--------------------------

p_hic5 <- plot_proportion(
  fibro_summary, "prop_TGFB1I1", "TGFB1I1",
  prop_tests$TGFB1I1
)

p_acta2 <- plot_proportion(
  fibro_summary, "prop_ACTA2", "ACTA2",
  prop_tests$ACTA2
)

p_col1a1 <- plot_proportion(
  fibro_summary, "prop_COL1A1", "COL1A1",
  prop_tests$COL1A1
)

p_col1a2 <- plot_proportion(
  fibro_summary, "prop_COL1A2", "COL1A2",
  prop_tests$COL1A2
)

p_cthrc1 <- plot_proportion(
  fibro_summary, "prop_CTHRC1", "CTHRC1",
  prop_tests$CTHRC1
)


#--------------------------
# Mean Expression Levels Among Positive Cells
#--------------------------

p_hic5_expr <- plot_expression(
  fibro_summary, "mean_TGFB1I1", "TGFB1I1",
  expr_tests$TGFB1I1
)

p_acta2_expr <- plot_expression(
  fibro_summary, "mean_ACTA2", "ACTA2",
  expr_tests$ACTA2
)

p_col1a1_expr <- plot_expression(
  fibro_summary, "mean_COL1A1", "COL1A1",
  expr_tests$COL1A1
)

p_col1a2_expr <- plot_expression(
  fibro_summary, "mean_COL1A2", "COL1A2",
  expr_tests$COL1A2
)

p_cthrc1_expr <- plot_expression(
  fibro_summary, "mean_CTHRC1", "CTHRC1",
  expr_tests$CTHRC1
)


#--------------------------
# Arrange Plots
#--------------------------

# Proportions of positive cells
prop_plots <- wrap_plots(
  p_hic5, p_acta2, p_col1a1, p_col1a2, p_cthrc1,
  ncol = 2
)

# Mean expression levels among positive cells
expr_plots <- wrap_plots(
  p_hic5_expr, p_acta2_expr, p_col1a1_expr,
  p_col1a2_expr, p_cthrc1_expr,
  ncol = 2
)

# Display the plots
prop_plots
expr_plots



#--------------------------
# 14.2. Proportions of ACTA2-, COL1A1-, COL1A2-, and
#       CTHRC1-Positive Cells Among TGFB1I1-Positive Fibroblasts
#--------------------------

#--------------------------
# Calculate the Proportion of Positive Cells
# Among TGFB1I1-Positive Fibroblasts
#--------------------------

# Function to calculate the proportion of cells positive
# for a target gene among TGFB1I1-positive fibroblasts.
calculate_positive_proportion <- function(data, gene) {
  
  hic5_pos <- data$TGFB1I1 > 0
  target_pos <- data[[gene]] > 0
  
  data %>%
    mutate(
      Hic5_pos = hic5_pos,
      target_pos = target_pos
    ) %>%
    group_by(sample, condition) %>%
    summarise(
      proportion = sum(Hic5_pos & target_pos) / sum(Hic5_pos),
      .groups = "drop"
    )
}


# Calculate proportions for each target gene.
prop_ACTA2_in_Hic5 <- calculate_positive_proportion(
  fibro_df, "ACTA2"
)

prop_COL1A1_in_Hic5 <- calculate_positive_proportion(
  fibro_df, "COL1A1"
)

prop_COL1A2_in_Hic5 <- calculate_positive_proportion(
  fibro_df, "COL1A2"
)

prop_CTHRC1_in_Hic5 <- calculate_positive_proportion(
  fibro_df, "CTHRC1"
)


#--------------------------
# Calculate Group Means and Standard Errors
# of the Mean (SEM)
#--------------------------

calculate_proportion_stats <- function(data) {
  
  data %>%
    group_by(condition) %>%
    summarise(
      mean_proportion = mean(proportion),
      sem_proportion = sd(proportion) / sqrt(n()),
      .groups = "drop"
    )
}


prop_ACTA2_in_Hic5_stats <- calculate_proportion_stats(
  prop_ACTA2_in_Hic5
)

prop_COL1A1_in_Hic5_stats <- calculate_proportion_stats(
  prop_COL1A1_in_Hic5
)

prop_COL1A2_in_Hic5_stats <- calculate_proportion_stats(
  prop_COL1A2_in_Hic5
)

prop_CTHRC1_in_Hic5_stats <- calculate_proportion_stats(
  prop_CTHRC1_in_Hic5
)


#--------------------------
# Statistical Analysis
#--------------------------

# Function to perform the Shapiro-Wilk normality test
# and Wilcoxon rank-sum test.
perform_proportion_tests <- function(data) {
  
  list(
    shapiro_control = shapiro.test(
      data$proportion[data$condition == "Control"]
    ),
    shapiro_ipf = shapiro.test(
      data$proportion[data$condition == "IPF"]
    ),
    wilcoxon = wilcox.test(
      proportion ~ condition,
      data = data
    )
  )
}


# Perform statistical tests for each target gene.
tests_ACTA2 <- perform_proportion_tests(prop_ACTA2_in_Hic5)
tests_COL1A1 <- perform_proportion_tests(prop_COL1A1_in_Hic5)
tests_COL1A2 <- perform_proportion_tests(prop_COL1A2_in_Hic5)
tests_CTHRC1 <- perform_proportion_tests(prop_CTHRC1)


# Check Shapiro-Wilk test results.
tests_ACTA2$shapiro_control
# W = 0.82557, p-value = 0.001658

tests_ACTA2$shapiro_ipf
# W = 0.86202, p-value = 0.001119

tests_COL1A1$shapiro_control
# W = 0.93066, p-value = 0.1419

tests_COL1A1$shapiro_ipf
# W = 0.93947, p-value = 0.08799

tests_COL1A2$shapiro_control
# W = 0.92715, p-value = 0.1207

tests_COL1A2$shapiro_ipf
# W = 0.68183, p-value = 8.579e-07

tests_CTHRC1$shapiro_control
# W = 0.7136, p-value = 4.087e-05

tests_CTHRC1$shapiro_ipf
# W = 0.9549, p-value = 0.2283


# Check Wilcoxon rank-sum test results.
tests_ACTA2$wilcoxon
# p-value = 0.6014697

tests_COL1A1$wilcoxon
# p-value = 0.003709152

tests_COL1A2$wilcoxon
# p-value = 0.01237655

tests_CTHRC1$wilcoxon
# p-value = 0.0002862374


#--------------------------
# Plotting
#--------------------------

# Function to generate violin plots showing the proportion
# of target gene-positive cells among TGFB1I1-positive fibroblasts.
plot_positive_proportion <- function(data, gene, test) {
  
  ggplot(
    data,
    aes(
      x = condition,
      y = proportion * 100,
      fill = condition
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.4
    ) +
    geom_boxplot(
      width = 0.1,
      fill = "white",
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.1,
      height = 0,
      size = 1
    ) +
    labs(
      y = paste0(
        gene, "+ among TGFB1I1+ fibroblasts (%)"
      ),
      x = ""
    ) +
    annotate(
      "text",
      x = 1.5,
      y = 120,
      label = paste0(
        "p = ",
        signif(test$p.value, 4)
      )
    ) +
    scale_fill_manual(
      values = c(
        "Control" = "skyblue",
        "IPF" = "tomato"
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    )
}


#--------------------------
# Generate Individual Plots
#--------------------------

p_acta2_in_hic5 <- plot_positive_proportion(
  prop_ACTA2_in_Hic5,
  "ACTA2",
  tests_ACTA2$wilcoxon
)

p_col1a1_in_hic5 <- plot_positive_proportion(
  prop_COL1A1_in_Hic5,
  "COL1A1",
  tests_COL1A1$wilcoxon
)

p_col1a2_in_hic5 <- plot_positive_proportion(
  prop_COL1A2_in_Hic5,
  "COL1A2",
  tests_COL1A2$wilcoxon
)

p_cthrc1_in_hic5 <- plot_positive_proportion(
  prop_CTHRC1_in_Hic5,
  "CTHRC1",
  tests_CTHRC1$wilcoxon
)


# Display individual plots.
p_acta2_in_hic5
p_col1a1_in_hic5
p_col1a2_in_hic5
p_cthrc1_in_hic5



# ============================================================
# 14.3. Comparison of ACTA2, COL1A1, COL1A2, and CTHRC1
#       expression between Hic-5-positive and
#       Hic-5-negative IPF fibroblasts
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#--------------------------
# 14.3.1. Prepare sample-level expression data
#--------------------------

# Genes of interest
genes <- c("ACTA2", "COL1A1", "COL1A2", "CTHRC1")

# Select IPF fibroblasts and classify Hic-5 status
fibro_IPF <- fibro_df %>%
  filter(condition == "IPF") %>%
  mutate(
    Hic5_status = ifelse(TGFB1I1 > 0, "Hic5+", "Hic5-")
  )

# Calculate mean expression for each gene
# in each sample and Hic-5 status
calculate_sample_expression <- function(data, gene) {
  
  expression_data <- data %>%
    group_by(sample, Hic5_status) %>%
    summarise(
      expression = mean(.data[[gene]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Reshape data to obtain paired values
  paired_data <- expression_data %>%
    pivot_wider(
      names_from = Hic5_status,
      values_from = expression
    ) %>%
    filter(
      !is.na(`Hic5-`),
      !is.na(`Hic5+`)
    )
  
  return(paired_data)
}

# Prepare sample-level data for each gene
expr_ACTA2 <- calculate_sample_expression(fibro_IPF, "ACTA2")
expr_COL1A1 <- calculate_sample_expression(fibro_IPF, "COL1A1")
expr_COL1A2 <- calculate_sample_expression(fibro_IPF, "COL1A2")
expr_CTHRC1 <- calculate_sample_expression(fibro_IPF, "CTHRC1")


#--------------------------
# 14.3.2. Paired Wilcoxon signed-rank tests
#--------------------------

# Perform paired Wilcoxon signed-rank test
perform_paired_wilcoxon <- function(data) {
  
  wilcox.test(
    data$`Hic5-`,
    data$`Hic5+`,
    paired = TRUE
  )
}

test_ACTA2 <- perform_paired_wilcoxon(expr_ACTA2)
test_COL1A1 <- perform_paired_wilcoxon(expr_COL1A1)
test_COL1A2 <- perform_paired_wilcoxon(expr_COL1A2)
test_CTHRC1 <- perform_paired_wilcoxon(expr_CTHRC1)

# Display statistical test results
test_ACTA2
test_COL1A1
test_COL1A2
test_CTHRC1

# Extract p-values
p_values <- tibble(
  Gene = genes,
  p_value = c(
    test_ACTA2$p.value,
    test_COL1A1$p.value,
    test_COL1A2$p.value,
    test_CTHRC1$p.value
  )
)

print(p_values)

# Previously reported p-values in Figure S3:
# ACTA2:  0.002551
# COL1A1: 3.05e-05
# COL1A2: 0.0001186
# CTHRC1: 6.326e-05


#--------------------------
# 14.3.3. Prepare data for visualization
#--------------------------

# Convert sample-level data into long format
prepare_plot_data <- function(data, gene) {
  
  data %>%
    select(sample, `Hic5-`, `Hic5+`) %>%
    pivot_longer(
      cols = c(`Hic5-`, `Hic5+`),
      names_to = "Hic5_status",
      values_to = "expression"
    ) %>%
    mutate(
      Hic5_status = factor(
        Hic5_status,
        levels = c("Hic5-", "Hic5+")
      ),
      Gene = gene
    )
}

plot_data_ACTA2 <- prepare_plot_data(expr_ACTA2, "ACTA2")
plot_data_COL1A1 <- prepare_plot_data(expr_COL1A1, "COL1A1")
plot_data_COL1A2 <- prepare_plot_data(expr_COL1A2, "COL1A2")
plot_data_CTHRC1 <- prepare_plot_data(expr_CTHRC1, "CTHRC1")


#--------------------------
# 14.3.4. Generate violin plots
#--------------------------

plot_gene_expression <- function(data, gene, p_value) {
  
  # Determine y-axis position for statistical annotation
  y_max <- max(data$expression, na.rm = TRUE)
  y_position <- y_max * 1.15
  
  # Convert p-value to significance label
  p_label <- case_when(
    p_value < 0.0001 ~ "****",
    p_value < 0.001  ~ "***",
    p_value < 0.01   ~ "**",
    p_value < 0.05   ~ "*",
    TRUE             ~ "ns"
  )
  
  ggplot(
    data,
    aes(x = Hic5_status, y = expression, fill = Hic5_status)
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      color = "black",
      alpha = 0.7
    ) +
    geom_point(
      position = position_jitter(width = 0.08),
      size = 1.5,
      alpha = 0.7
    ) +
    annotate(
      "segment",
      x = 1, xend = 2,
      y = y_position, yend = y_position
    ) +
    annotate(
      "segment",
      x = 1, xend = 1,
      y = y_position - y_max * 0.04,
      yend = y_position
    ) +
    annotate(
      "segment",
      x = 2, xend = 2,
      y = y_position - y_max * 0.04,
      yend = y_position
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_position + y_max * 0.04,
      label = p_label,
      size = 5
    ) +
    scale_fill_manual(
      values = c(
        "Hic5-" = "grey80",
        "Hic5+" = "tomato"
      )
    ) +
    scale_x_discrete(
      labels = c(
        "Hic5-" = "Hic-5−",
        "Hic5+" = "Hic-5+"
      )
    ) +
    labs(
      x = NULL,
      y = paste("Mean", gene, "expression")
    ) +
    coord_cartesian(
      ylim = c(0, y_position + y_max * 0.15)
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "none"
    )
}

# Generate plots for each gene
p_ACTA2 <- plot_gene_expression(
  plot_data_ACTA2, "ACTA2", test_ACTA2$p.value
)

p_COL1A1 <- plot_gene_expression(
  plot_data_COL1A1, "COL1A1", test_COL1A1$p.value
)

p_COL1A2 <- plot_gene_expression(
  plot_data_COL1A2, "COL1A2", test_COL1A2$p.value
)

p_CTHRC1 <- plot_gene_expression(
  plot_data_CTHRC1, "CTHRC1", test_CTHRC1$p.value
)

# Display all four plots
(p_ACTA2 | p_COL1A1) /
  (p_COL1A2 | p_CTHRC1)

  

# ============================================================
# 15. Reactome pathway enrichment analysis of the unknown
#     fibroblast cluster
# ============================================================

#--------------------------
# 15.1. Identification of upregulated DEGs
#--------------------------

# Set fibroblast cluster identities
Idents(fibro) <- "fibro.annotated_clusters"

# Identify upregulated DEGs in the Unknown cluster
# compared with all other fibroblast clusters
degs <- FindMarkers(
  fibro,
  ident.1 = "Unknown",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


#--------------------------
# 15.2. Conversion of gene symbols to Entrez IDs
#--------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

# Extract gene symbols from the DEG results
gene_symbols <- rownames(degs)

# Convert gene symbols to Entrez Gene IDs
entrez_ids <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


#--------------------------
# 15.3. Reactome pathway enrichment analysis
#--------------------------

library(ReactomePA)

# Perform Reactome pathway enrichment analysis
reactome_res <- enrichPathway(
  gene         = entrez_ids$ENTREZID,
  organism     = "human",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

# Display enrichment results
head(reactome_res)


#--------------------------
# 15.4. Visualization of Reactome enrichment results
#--------------------------

# Dot plot
dotplot(reactome_res, showCategory = 10) +
  ggtitle("Reactome Pathway Enrichment in Unknown Cluster") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

# Bar plot
barplot(reactome_res, showCategory = 10) +
  ggtitle("Reactome Pathway Enrichment in Unknown Cluster") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.margin = margin(5, 5, 5, 10)
  )




