#!/usr/bin/env Rscript
# read in the data
pb2_csvPath <- "/home/dgladish/projects/rrg-tperkins/dgladish/Projects/MSc/Brand_Perkins/T-ALL_Clonal/CITE-seq/data/PB2_rna_CITE.csv"
th1_csvPath <- "/home/dgladish/projects/rrg-tperkins/dgladish/Projects/MSc/Brand_Perkins/T-ALL_Clonal/CITE-seq/data/TH1_rna_CITE.csv"
th2_csvPath <- "/home/dgladish/projects/rrg-tperkins/dgladish/Projects/MSc/Brand_Perkins/T-ALL_Clonal/CITE-seq/data/TH2_rna_CITE.csv"

pb2_raw <- read.csv(pb2_csvPath)
rownames_pb2 <- pb2_raw$X
rownames(pb2_raw) <- rownames_pb2
pb2_raw$X <- NULL

th1_raw <- read.csv(th1_csvPath)
rownames_th1 <- th1_raw$X
rownames(th1_raw) <- rownames_th1
th1_raw$X <- NULL

th2_raw <- read.csv(th2_csvPath)
rownames_th2 <- th2_raw$X
rownames(th2_raw) <- rownames_th2
th2_raw$X <- NULL

# add sample-specific suffixes
colnames(pb2_raw) <- paste0(colnames(pb2_raw), "-PB2")
colnames(th1_raw) <- paste0(colnames(th1_raw), "-TH1")
colnames(th2_raw) <- paste0(colnames(th2_raw), "-TH2")

# merge the TH1 and TH2 matrices
merged_norm <- merge(th1_raw, th2_raw, by="row.names", all=TRUE)
rm(th1_csvPath)
rm(th2_csvPath)

rownames(merged_norm) <- merged_norm$Row.names
merged_norm$Row.names <- NULL

# subset the normal cell data frame
set.seed(123)

norm_names <- colnames(merged_norm)
print(head(norm_names))
print(length(norm_names))

sampled_norm_names <- sample(norm_names, 350)
print(length(sampled_norm_names))

merged_norm_subset <- merged_norm[, sampled_norm_names]
print(ncol(merged_norm_subset))

# merge the leukemia sample and TH1+TH2 subset
merged_all <- merge(th1_raw, th2_raw, by="row.names", all=TRUE)

# create the vector of normal cells from the TH1+TH2 subset
norm_names_subset <- colnames(merged_norm_subset)

# run SCEVAN, specifying pipelineCNA(norm_cell = norm_names_subset)
results <- SCEVAN::pipelineCNA(, norm_cell = norm_names_subset, par_cores = 1, SUBCLONES = TRUE, plotTree = TRUE, sample="")
saveRDS()
