#### DE by time points ####
### F: T1 v T0, T2 v T1, T3 v T2
### M: T1 v T0, T2 v T1, T3 v T2

# 
filenames <- grep("htseqcount",list.files(),value=TRUE)
filenames

sample_names <- read.table("./biosamples_combined_2.txt", sep = "\t", header = T)
head(sample_names)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           file_name = filenames, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time),
                           ind = factor(sample_names$ind))
head(sample_table)

# only F
sample_table <- subset(sample_table, sex=="F")
dim(sample_table) # 18 bc 2 missing females at T0

# create a DEseq2 object, using time and replicates as covariate
# library(DESeq2)
# dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = getwd(),
#                                    design = ~ ind + time)
# 
# dds
# 
# # filtering
# library(edgeR)
# 
# x <- counts(dds)
# 
# minCPM = 0.5
# nLibraries = 1
# # save the CPM counts
# #write.csv(cpm(DGEList(counts = x)), file = "CPM_matrix.csv")
# # CMP counts for one gene
# #plot(cpm(DGEList(counts = x))["G14641",], col = sample_table$sex)
# 
# # keep genes based on minCPM in nlibraries, idep rule
# x <- x[ which( apply( cpm(DGEList(counts = x)),  1, function(y) sum(y>=minCPM)) >=  nLibraries ) , ]
# dim(x)
# # rownames(x)
# # dds_filtered <- dds[rownames(x),]
# # dds <- dds_filtered
# # x["G27464",]
# 
# # save 
# saveRDS(x, file = "counts_combined_F_only_minCPM0.5_nlib1.rds")

# load

x <- readRDS("~/counts_combined_F_only_minCPM0.5_nlib1.rds")

dim(x)

# Run DESeq2--------------------

library(DESeq2) # D.E.G.
FC <- 2 # Fold-change cutoff
FDR <- 0.05 # FDR cutoff
alpha <- 0.1 # independent filtering, default

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~ ind + time)

dds <- DESeq2::DESeq(dds)

# Extract results--------------------

res <- DESeq2::results(dds, 
                       contrast = c("time", "T1", "T0"), # T0 as baseline
                       independentFiltering = TRUE,
                       alpha = alpha)
# Examine results 
summary(res)
plotMA(res)

res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

write.table(res, file = "DEG_table_T1vT0.txt", quote = F)

### contrast T2 v T1 ###

res <- DESeq2::results(dds, 
                       contrast = c("time", "T2", "T1"),
                       independentFiltering = TRUE,
                       alpha = alpha)

# Examine results 
summary(res)
plotMA(res)
#plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
# if fgsea, do not subset
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated
write.table(res, file = "DEG_table_T2vT1.txt", quote = F)

### contrast T3 v T2 ###

res <- DESeq2::results(dds, 
                       contrast = c("time", "T3", "T2"),
                       independentFiltering = TRUE,
                       alpha = alpha)

# Examine results 
summary(res)
plotMA(res)
#plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
# if fgsea, do not subset
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated
write.table(res, file = "DEG_table_T3vT2.txt", quote = F)



#### Males #####

sample_names <- read.table("./biosamples_combined_2.txt", sep = "\t", header = T)
head(sample_names)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           file_name = filenames, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time),
                           ind = factor(sample_names$ind))
head(sample_table)

# only M
sample_table <- subset(sample_table, sex=="M")
dim(sample_table) 

# # create a DEseq2 object
# library(DESeq2)
# dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = getwd(),
#                                   design = ~ time + ind)
# 
# dds
# 
# 
# # filter/preprocess 
# 
# 
# # filtering
# library(edgeR)
# 
# x <- counts(dds)
# 
# minCPM = 0.5
# nLibraries = 1
# # save the CPM counts
# #write.csv(cpm(DGEList(counts = x)), file = "CPM_matrix.csv")
# # CMP counts for one gene
# #plot(cpm(DGEList(counts = x))["G14641",], col = sample_table$sex)
# 
# # keep genes based on minCPM in nlibraries, idep rule
# x <- x[ which( apply( cpm(DGEList(counts = x)),  1, function(y) sum(y>=minCPM)) >=  nLibraries ) , ]
# dim(x)
# # rownames(x)
# # dds_filtered <- dds[rownames(x),]
# # dds <- dds_filtered
# # x["G27464",]
# 
# # save 
# saveRDS(x, file = "counts_combined_M_only_minCPM0.5_nlib1.rds")
# 
# # load

x <- readRDS("./counts_combined_M_only_minCPM0.5_nlib1.rds")

dim(x)

# DEG

# Run DESeq2--------------------

library(DESeq2) # D.E.G.
FC <- 2 # Fold-change cutoff
FDR <- 0.05 # FDR cutoff
alpha <- 0.1 # independent filtering, default

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~ time + ind
)
dds <- DESeq2::DESeq(dds)


# Extract results--------------------

# T1 v T0
res <- DESeq2::results(dds, 
                       contrast = c("time", "T1", "T0"),
                       independentFiltering = TRUE,
                       alpha = alpha
)
# Examine results 
summary(res)
plotMA(res)

res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

write.table(res, file = "DEG_table_T1vT0.txt", quote = F)

### contrast T2 v T1 ###

res <- DESeq2::results(dds, 
                       contrast = c("time", "T2", "T1"),
                       independentFiltering = TRUE,
                       alpha = alpha)

# Examine results 
summary(res)
plotMA(res)
#plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
# if fgsea, do not subset
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated
write.table(res, file = "DEG_table_T2vT1.txt", quote = F)

### contrast T3 v T2 ###

res <- DESeq2::results(dds, 
                       contrast = c("time", "T3", "T2"),
                       independentFiltering = TRUE,
                       alpha = alpha)

# Examine results 
summary(res)
plotMA(res)
#plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
# if fgsea, do not subset
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated
write.table(res, file = "DEG_table_T3vT2.txt", quote = F)


