# SPE score and SEL score
## Introduction
We developed an approach to systematically analyze the regulatory networks between enhancers and promoters across different cell types at the genome scale based on the E-P interactions predicted by the ABC model [1]. 
ABC E-P links for more celltypes can be downloaded [here](https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz).

## Steps
```r
setwd("Workdir")
rm(list=ls())
library(ggplot2)
library(tidyr)
library(factoextra)
source("./allfunction.R")
```
### 1. Creat E-P link binary matrix

![image](https://github.com/xmuhuanglab/eSelectivity/blob/main/image/step1.jpg)

```r
# Create a binarymatrix where rows represent all EP links and columns represent all cell types.
datapath <- "./data/"
suffix <- "_EP.txt"
allEP_binary <- binaryMatrix(datapath, # The filepath of EP link files for each cell
                             suffix # The suffix following the cell type in the filename
)

# save
saveRDS(allEP_binary, file = "allEP_binary.rds")
```
### 2. Separate E-P link binary matrix for each gene
```r
allEP_binary <- readRDS(file = "allEP_binary.rds")
pathout <- "./gene_EP/"
if (!dir.exists(pathout)){dir.create(pathout)}
sepBinaryMatrix(allEP_binary, # The binary matrix of all EP (the output of step1) 
                pathout # The output path of EP binary matrix of each gene
)
```
### 3. Caculate SPE score and SEL score for each gene

![image](https://github.com/xmuhuanglab/eSelectivity/blob/main/image/step2.jpg)

```r
allEP_binary <- readRDS(file = "allEP_binary.rds")
seppath <- "./gene_EP/"

score_matrix <- calculateScore(allEP_binary, # The binary matrix of all EP (the output of step1)
                               seppath # The filepath of EP binary matrix of each gene(the output of step2)
)
							   
write.table(score_matrix, "SPEscore_SELscore_matrix_raw.txt", sep="\t", quote = F, row.names = F, col.names = T)
```
### 4. K-means cluster and plot
```r
score_matrix <- read.delim("SPEscore_SELscore_matrix_raw.txt", skip=0, header=T, check.names = FALSE,sep = "\t",stringsAsFactors = F)

# cluster
set.seed(1234)
km_result <- kmeans(scale(score_matrix[,2:3]),3, nstart = 50)
# merge cluster result and raw matrix
score_matrix_cluster = cbind(score_matrix, cluster = km_result$cluster)

write.table(score_matrix_cluster,file = "spescore_selscore_matrix_kmeans_cluster.txt", sep="\t", quote = F, row.names = F, col.names = T)

# plot
score_matrix_cluster <- read.delim("spescore_selscore_matrix_kmeans_cluster.txt", skip=0, header=T, check.names = FALSE,sep = "\t",stringsAsFactors = F)
score_matrix_cluster$cluster <- factor(score_matrix_cluster$cluster, levels = c(1:3))

palette <- c("#B5C3C4","#BAD873","#C16EC6")

p <- ggplot(score_matrix_cluster, aes(x=SPE_score, y=SEL_score, color = cluster)) +
  geom_point(size = 2,shape = 20) +
  scale_color_manual(values = palette) +
  labs(x="SPE score",y="SEL score") +
  theme_bw() +
  theme(
    legend.position = c(0.2,0.7),
    axis.text.x =element_text(size=14,vjust = 0.5, hjust = 0.5, angle = 0), axis.text.y=element_text(size=14),
    axis.title.x =element_text(size=14,vjust = 0.5, hjust = 0.5), axis.title.y=element_text(size=14),
    legend.text=element_text(size=14),legend.title=element_text(size=14),
    panel.grid=element_blank(),panel.border = element_rect(size=0.5)
  )
p
ggsave(p,filename = "spescore_selscore_cluster_dotplot.jpg",height = 4,width = 4,dpi = 300)
```

![image](https://github.com/xmuhuanglab/eSelectivity/blob/main/image/step4.jpg)


## Reference
1. [Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR, Patwardhan TA, Jones TR, Nguyen TH, Ulirsch JC, Lekschas F et al: Genome-wide enhancer maps link risk variants to disease genes. Nature 2021, 593(7858):238-243.](https://doi.org/10.1038/s41586-021-03446-x
        
        
        
        )
