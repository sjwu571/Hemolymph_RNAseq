# Time point A as an example

# Download htseqcount files to a folder 

setwd("./RNA-seq/A")

filenames <- grep("htseqcount",list.files(),value=TRUE)
filenames

#?DESeqDataSetFromHTSeqCount

sample_names <- read.table("./biosample_result_A.txt", sep = "\t", header = T)
head(sample_names)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           file_name = filenames, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time))
sample_table

# create a DEseq2 object
library(DESeq2)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = getwd(),
                                  design = ~ sex)

dds


# filter/preprocess 


# filtering
library(edgeR)

x <- counts(dds)

minCPM = 0.5
nLibraries = 1
# save the CPM counts
#write.csv(cpm(DGEList(counts = x)), file = "CPM_matrix.csv")
# CMP counts for one gene
#plot(cpm(DGEList(counts = x))["G14641",], col = sample_table$sex)

# keep genes based on minCPM in nlibraries, idep rule
x <- x[ which( apply( cpm(DGEList(counts = x)),  1, function(y) sum(y>=minCPM)) >=  nLibraries ) , ]
dim(x)

# save 
saveRDS(x, file = "counts_A_minCPM0.5_nlib1.rds")

# load filtered count matrix

x <- readRDS("./counts_A_minCPM0.5_nlib1.rds")

dim(x)


# DEG

# Run DESeq2--------------------

library(DESeq2) # D.E.G.
FC <- 2 # Fold-change cutoff
FDR <- 0.05 # FDR cutoff
alpha <- 0.1 # independent filtering, default

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_names,
  design = ~sex)
dds <- DESeq2::DESeq(dds)


# Extract results--------------------

# Comparison 1 of 1:  F-M
res <- DESeq2::results(dds, 
                       contrast = c("sex", "F", "M"),
                       independentFiltering = TRUE,
                       alpha = alpha)
# Examine results 
summary(res)
plotMA(res)
#plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])

res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated
write.table(res, file = "DEG_table_fdr05.txt", quote = F)
