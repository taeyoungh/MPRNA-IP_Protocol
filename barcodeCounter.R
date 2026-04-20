#! /usr/bin/Rscript --vanilla

suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))

############################
# Description
# Last update: 1/04/2026

# Commnad Line Usage example for two input files
# barcodeCounter.R -i sample1.mapped,sample2.mapped -p OLIGO_POOL.fasta -t 2,2 -o output_prefix -n sample1,sample2
############################


############################
# 1. Command line arguments
############################

option_list <- list(make_option(c("-i", "--inputFiles"), type="character", default=NULL, help="Mapped file name. For multiple files, put comma-separated list"),
                    make_option(c("-p", "--oligoPool"), type="character", default=NULL, help="Fasta file of oligo pool"),
                    make_option(c("-t", "--mmThreshold"), type="character", default=NULL, help="Mismatch threshold for counting, For multiple files, put comma-separated list"),
                    make_option(c("-o", "--output"), type="character", default="barcodeCounter", help="Output file name"),
                    make_option(c("-n", "--sampleNames"), type="character", default=NULL, help="Sample names, For multiple files, put comma-separated list"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt) # print log
if (is.null(opt$sampleNames) | is.null(opt$inputFiles) | is.null(opt$oligoPool) | is.null(opt$mmThreshold)) {
  stop("One of required arguments is missing\n", call.=FALSE)
}

############################
# 2. Preprocessing
############################

sampleNames <- unlist(strsplit(opt$sampleNames, split=","))
inputFiles <- unlist(strsplit(opt$inputFiles, split=","))
mmThreshold <- setNames(as.numeric(unlist(strsplit(opt$mmThreshold, split=","))), sampleNames)

if (length(sampleNames) != length(inputFiles) | length(sampleNames) != length(mmThreshold) | length(inputFiles) != length(mmThreshold)) {
  stop("Arguments are different in size\n", call.=FALSE)
}

# read oligoPool fasta file
temp <- read.table(opt$oligoPool, sep="\n", header=F)[,1]
temp <- strsplit(gsub("^>|\\..*$", "", temp[seq(1, length(temp), 2)]), split="_")
oligoPool <- data.frame(TileID = paste(sapply(temp, "[[", 1), sapply(temp, "[[", 2), sep="_"),
                        BarcodeID = sapply(temp, "[[", 3))

############################
# 3. Process samples
############################

mismatch.summary <- list()

for (i in 1:length(sampleNames)) {
  cat(paste0("Sample ", i, ": ", sampleNames[i], "\n"))

  # read a mapped file.
  cat(paste0("\tReading ", inputFiles[i], "\n"))
  align <- read.table(inputFiles[i], header=F, sep="\t", colClasses = c("NULL", "NULL", "character", "character", "numeric")) # don't read qname and barcodeSeq.
  colnames(align) <- c("BarcodeID", "TileID", "MismatchNum")
  n <- nrow(align)

  # summary by mismatch number
  cat(paste0("\tSummarizing mismatches...\n"))
  mismatch.table <- as.data.frame.table(table(align$MismatchNum))
  colnames(mismatch.table) <- c("MismatchNum", "Number") # -1 means unmapped.
  mismatch.table <- mutate(mismatch.table, Prop = Number/n)
  mismatch.summary[[sampleNames[i]]] <- mismatch.table

  # count matrix
  cat(paste0("\tMaking count table... \n"))
  temp <- subset(align, MismatchNum<=mmThreshold[i] & MismatchNum>=0) %>% dplyr::count(TileID, BarcodeID, name=sampleNames[i])
  oligoPool <- left_join(oligoPool, temp, by=c("TileID", "BarcodeID"))
  oligoPool[[sampleNames[i]]][is.na(oligoPool[[sampleNames[i]]])] <- 0 # set missing oligo count as 0
}

############################
# 4. Outputs
############################

mismatch.summary = bind_rows(mismatch.summary, .id="Sample")
write.table(mismatch.summary, file=paste0(opt$output, "_qc.txt"), row.names = F, col.names = T, quote=F, sep="\t")

oligoPool = mutate(oligoPool, oligoID=paste(TileID, BarcodeID, sep="_")) %>% dplyr::select(-TileID, -BarcodeID)
rownames(oligoPool) = oligoPool$oligoID
write.table(dplyr::select(oligoPool, -oligoID), file=paste0(opt$output, "_count.csv"), row.names = T, col.names = NA, quote=F, sep=",")
