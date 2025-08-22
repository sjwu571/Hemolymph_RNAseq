#### Temporal analysis #####

# Bioconductor.pck <- c("SummarizedExperiment", "S4Vectors", "DESeq2",
#                       "Mfuzz", "ComplexHeatmap")
# BiocManager::install(version="3.18")

# Select.package.Bioc <- "MultiRNAflow"
# if(!require(package=Select.package.Bioc,
#             quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)){
#   BiocManager::install(pkgs=Select.package.Bioc)
# }## if(!require(package=Select.package.Bioc, quietly=TRUE, character.only=TRUE))

library(MultiRNAflow)

# Removed the two females missing t0, altogether.
# named individual replicate for each sample, r1, r2... r10.
# read in data, unfiltered count matrix with colnames changed
hemo_combined <- read.csv("./counts_matrix_combined_Mfuzz.csv", header = T)
colnames(hemo_combined)

SEresleuk500 <- DATAprepSE(RawCounts=hemo_combined,
                           Column.gene=1,
                           Group.position=1,
                           Time.position=2,
                           Individual.position=3)

SEresNORMleuk500 <- DATAnormalization(SEres=SEresleuk500,
                                      Normalization="vst",
                                      Blind.rlog.vst=FALSE,
                                      Plot.Boxplot=FALSE,
                                      Colored.By.Factors=TRUE,
                                      Color.Group=NULL,
                                      path.result=NULL)

SEresMfuzzLeuk500 <- MFUZZanalysis(SEresNorm=SEresNORMleuk500,
                                   DataNumberCluster=NULL,
                                   Method="hcpc", #default, or kmeans
                                   Membership=0.5, # prob of gene belong to each cluster
                                   Min.std=0.1,
                                   Plot.Mfuzz=TRUE,
                                   path.result=NULL,
                                   Name.folder.mfuzz="./")

# 0 genes excluded.
# 9819 genes excluded.
# 0 genes excluded.
# 11036 genes excluded.
# 
#res_Mfuzz <- SEresMfuzzLeuk500@metadata$MFUZZ
# I want to extract all gene ID for each cluster, both females and males
Iwant <- SEresMfuzzLeuk500@metadata$Results$UnsupervisedAnalysis$Mfuzz$Result.Mfuzz
head(Iwant)
str(Iwant)
# how many genes in each cluster for females
table(Iwant$Cluster_F)
table(Iwant$Cluster_M)
dim(Iwant)
summary(Iwant$Cluster_F)
str(Iwant)

# 
Cluster1_F <- Iwant[Iwant$Cluster_F==1 & !is.na(Iwant$Cluster_F),]
Cluster2_F <- Iwant[Iwant$Cluster_F==2 & !is.na(Iwant$Cluster_F),]
Cluster3_F <- Iwant[Iwant$Cluster_F==3 & !is.na(Iwant$Cluster_F),]
Cluster4_F <- Iwant[Iwant$Cluster_F==4 & !is.na(Iwant$Cluster_F),]

# ALl the gene ID in Cluster 1 for females
rownames(Cluster1_F)

# Males

Cluster1_M <- Iwant[Iwant$Cluster_M==1 & !is.na(Iwant$Cluster_M),]
Cluster2_M <- Iwant[Iwant$Cluster_M==2 & !is.na(Iwant$Cluster_M),]
Cluster3_M <- Iwant[Iwant$Cluster_M==3 & !is.na(Iwant$Cluster_M),]


head(Cluster2_M)

# make df and export

my_list <- list(
  Cluster1_F = Cluster1_F,
  Cluster1_M = Cluster1_M,
  Cluster2_F = Cluster2_F,
  Cluster2_M = Cluster2_M,
  Cluster3_F = Cluster3_F,
  Cluster3_M = Cluster3_M,
  Cluster4_F = Cluster4_F
)

head(my_list$Cluster1_F$Gene.Name)

getwd()
save(my_list, file = "List_geneID_Mfuzz_Prob_0.5.RData")
