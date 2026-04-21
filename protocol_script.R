# ================================
# Step 129
# ================================

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("edgeR")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("Biostrings")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(dplyr)
library(reshape2)
library(edgeR)
library(Biostrings)
library(ggplot2)

# ================================
# Step 130
# ================================

# load oligo count matrix
count.mat = read.csv("GSE315146_TestPool_count_updated.csv", header=T, row.names = 1) %>% as.matrix()

# generate a sample sheet that has sample information
sampleSheet = data.frame(sampleID = colnames(count.mat), repl=c("ePCR", "Maxi", rep(c("Rep1", "Rep2", "Rep3"), each=2)), pulldown = c(NA, NA, rep(c("Input", "IP"), 3)))

# calculate total number of counts for every sample
sampleSheet$countedNum = colSums(count.mat)

# ================================
# Step 131
# ================================

# flat count data as a data.frame
oligo.count = reshape2::melt(count.mat, varnames = c("oligoID", "sampleID"), value.name = "raw")

oligo.count = mutate(oligo.count, tileID = gsub("_[^_]+$", "", oligoID), barcodeID = gsub("^.*_", "", oligoID))

oligo.count = left_join(oligo.count, sampleSheet[, c("sampleID", "repl", "pulldown")], by="sampleID")

# calculate RPM (Reads Per Million reads) of oligo
oligo.count$rpm <- oligo.count$raw / (sampleSheet$countedNum[match(oligo.count$sampleID, sampleSheet$sampleID)]/10^6)

# set the parameters
TILE_NUM = length(unique(oligo.count$tileID))
BARCODE_NUM = 15
REPLICATE_NUM = 3

# ================================
# Step 132
# ================================

# missing oligo and tile numbers
oligo.count.zero = subset(oligo.count, raw==0)

temp = reshape2::dcast(oligo.count.zero, sampleID~., value.var = "oligoID")

sampleSheet$mOligoNum <- temp$"."[match(sampleSheet$sampleID, temp$sampleID)]

sampleSheet$mOligoRate <- sampleSheet$mOligoNum / (TILE_NUM*BARCODE_NUM)

temp = reshape2::dcast(oligo.count.zero, sampleID~., value.var = "tileID", fun.aggregate = function(x) length(unique(x)))

sampleSheet$mTileNum <- temp$"."[match(sampleSheet$sampleID, temp$sampleID)]

sampleSheet$mTileRate <- sampleSheet$mTileNum / TILE_NUM

print(sampleSheet)

# make a plot of distribution of oligos
p = ggplot(oligo.count, aes(log10(rpm+1))) + facet_wrap(~sampleID)
p = p + geom_density()
p = p + xlab("Log10 (RPM+1)") + ylab("Density")
p = p + ggtitle("RPM of oligos")
p

# ================================
# Step 133
# ================================

# if necessary, keep the relevant fRIP samples only
oligo.count = subset(oligo.count, sampleID %in% c("R1_IN", "R1_IP", "R2_IN", "R2_IP", "R3_IN", "R3_IP"))
sampleSheet = subset(sampleSheet, sampleID %in% c("R1_IN", "R1_IP", "R2_IN", "R2_IP", "R3_IN", "R3_IP"))

# Define functions for filtering and pooling

oligoFilter <- function(oligo.count, inputFilterBy="rpm", inputTh=0, ipFilterBy="rpm", ipTh=0) {
  # this function filters oligonucleotides based on either "raw" or "rpm" values.
  
  oligo.count$oligoID <- paste(oligo.count$tileID, oligo.count$barcodeID, sep="_")  
  
  # Input filter
  temp <- subset(oligo.count, pulldown=="Input" & get(inputFilterBy)>=inputTh)
  temp <- table(temp$oligoID)
  temp <- names(temp)[temp==REPLICATE_NUM] # oligoIDs should be expressed in all the replicates.
  oligo.count <- subset(oligo.count, oligoID %in% temp)
  
  # IP filter
  temp <- subset(oligo.count, pulldown=="IP" & get(ipFilterBy)>=ipTh)
  temp <- table(temp$oligoID)
  temp <- names(temp)[temp==REPLICATE_NUM] # oligoIDs should be expressed in all the replicates.
  oligo.count <- subset(oligo.count, oligoID %in% temp)
  
  cat("Proportion of oligos after filtering: ", length(unique(oligo.count$oligoID))/(TILE_NUM*BARCODE_NUM), "\n")
  oligo.count$oligoID <- NULL
  return(oligo.count)
}

# barcode pooling function
countPooler <- function(oligo.count) {
  
  oligo.pooled <- dcast(oligo.count, tileID~sampleID, value.var="raw", fun.aggregate = sum)
  rownames(oligo.pooled) <- oligo.pooled$tileID
  pooled.count <- data.matrix(oligo.pooled[,-1]) # Output as a matrix 
  
  cat("Proportion of tiles after filtering and pooling: ", nrow(pooled.count)/TILE_NUM, "\n")
  return(pooled.count)
}

# Perform filtering and pooling

pooled.count <- oligoFilter(oligo.count, inputFilterBy="rpm", inputTh=1, ipFilterBy="rpm", ipTh=1) %>% countPooler()

# ================================
# Step 134
# ================================

# Normalization
edgeR.obj <- DGEList(pooled.count)
edgeR.obj <- calcNormFactors(edgeR.obj, method="TMM") 

# Voom to use count data for linear modeling
pooled.model <- model.matrix(data= sampleSheet, ~ 0 + pulldown + repl)
y <- voom(edgeR.obj, pooled.model, plot = T)

# Linear model fit by limma
fit <- lmFit(y, pooled.model)

# Linear model contrast
contr <- makeContrasts(pulldownIP - pulldownInput, levels = colnames(coef(fit)))
out <- contrasts.fit(fit, contr)

# Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
out <- eBayes(out)

# Results
tile.res <- topTable(out, sort.by = "P", n = Inf) %>% tibble::rownames_to_column("tileID")
# hist(tile.res$P.Value) # Optional: check the distribution of p-values.

# Enriched tiles
tile.res$enriched <- "Insig"
tile.res$enriched[which(tile.res$adj.P.Val<0.01 & tile.res$logFC>0)] <- "Enriched"
tile.res$enriched[which(tile.res$adj.P.Val<0.01 & tile.res$logFC<0)] <- "Depleted"
table(tile.res$enriched)

# ================================
# Step 135
# ================================

# Read a fasta file of tile sequences.
tile.seq = readDNAStringSet("GSE315146_tiles.fasta")

# Write a fasta file for the ordered tiles.
temp <- tile.seq[tile.res$tileID[order(tile.res$t, decreasing = T)]]
writeXStringSet(temp, "PUM2_ranked.fasta")

