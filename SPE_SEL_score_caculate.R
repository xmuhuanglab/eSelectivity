library(ggplot2)

path_in <- "./SPE_SEL/data/";
path_out <- "./SPE_SEL/";

file_ins <- list.files(path = path_in, pattern = '.*_EP.txt');
file_ins

### 1. E-P link binary matrix

# Create a data frame where rows represent all EP links and columns represent all cell types.

allEP_name <- c()
celltypes <- c()

for(i in 1:length(file_ins)){
  EP <- read.delim(paste0(path_in,file_ins[i]), skip=0, header=T, check.names = FALSE, sep = "\t", stringsAsFactors = F)
  EP_name <- unique(paste(EP[,1],EP[,2],sep = '#'))
  allEP_name <- union(EP_name,allEP_name)
  celltypes[i] <- gsub("_EP.txt","",file_ins[i])
}

allEP_name <- sort(unique(allEP_name))

EPall <- data.frame(matrix(0, nrow = length(allEP_name), ncol = length(celltypes)))
rownames(EPall) <- allEP_name
colnames(EPall) <- celltypes

for(i in 1:length(file_ins)){
  EP <- read.delim(paste0(path_in,file_ins[i]), skip=0, header=T, check.names = FALSE, sep = "\t", stringsAsFactors = F)
  EP_name <- unique(paste(EP[,1],EP[,2],sep = '#'))
  EPall[,i] <- match(rownames(EPall),EP_name)
}

EPall[!is.na(EPall)] <- 1
EPall[is.na(EPall)] <- 0

# Sort the rows based on gene name and enhancer chromosome position

c1 <- ncol(EPall)

EPall <- separate(as.data.frame(rownames(EPall)), 1, sep = "#", into = c("gene","eh"))
EPall <- separate(as.data.frame(EPall$eh), 1, sep = "_", into = c("chr","start","end"))
EPall[,(c1+4):(c1+5)] <- apply(EPall[,(c1+4):(c1+5)], 2, as.numeric)
EP <- EPall[order(EPall$gene, EPall$chr, EPall$start), 1:(c1+2)]

# save
saveRDS(EP, file = paste0(path_out,"EPall.rds"))

### 2. E-P link binary matrix for each gene

path_in <- "./SPE_SEL/"
path_out <- "./SPE_SEL/gene_EP/"
if (!dir.exists(path_out)){dir.create(path_out)}

EPall <- readRDS(file = paste0(path_in, "EPall.rds"))
c1 <- ncol(EPall)

gene_list <- as.character(unique(EPall[,"gene"]))

for(i in 1:length(gene_list)){#i=1
  EP_matrix <- EPall[which(EPall[,"gene"] %in% gene_list[i]),1:(c1-2)]
  
  # keep multiple enhancer-gene links (enhancer > 1)
  if(nrow(EP_matrix) >= 2){
    write.table(EP_matrix, file=paste0(path_out,gene_list[i],".txt"), sep="\t", quote = F, row.names = T, col.names = T)
  }
}

### 3. SPE score and SEL score for each gene

path_in <- "./SPE_SEL/gene_EP/"
path_out <- "./SPE_SEL/"
file_ins <- list.files(path = path_in, pattern = '.*.txt');
file_ins

EPall <- readRDS(file = paste0(path_out,"EPall.rds"))
c1 <- ncol(EPall)
celltypes <- colnames(EPall)[1:c1]

score_matrix <- data.frame(gene = "N", SPE_score = 0, SEL_score = 0, stringsAsFactors = F)

for(i in 1:length(file_ins)){#i=1
  gene <- gsub(".txt","",file_ins[i])
  
  each_EP <- read.delim(paste0(path_in,file_ins[i]), skip=0, header=T, check.names = FALSE,sep = "\t",stringsAsFactors = F)
  each_EP_t <- t(each_EP)
  
  cell_num <- nrow(each_EP_t)
  eh_num <- ncol(each_EP_t)
  
  spe_score <- mean(cell_num - colSums(each_EP_t) + 1) / cell_num
  sel_score <- sd(rowSums(each_EP_t) / eh_num) / mean(rowSums(each_EP_t) / eh_num)
  
  score_matrix[i,] <- c(gene, spe_score, sel_score)
  
  progress_value <- (i/length(file_ins))*100
  if(progress_value %% 5 == 0){print(paste0(progress_value,"%"))}
  
}

write.table(score_matrix, file = paste0(path_out,"SPEscore_SELscore_matrix_raw.txt"), sep="\t", quote = F, row.names = F, col.names = T)

### 4. Kmeans clustering

# cluster
library(factoextra)
score_matrix <- read.delim(paste(path_out,"SPEscore_SELscore_matrix_raw.txt",sep=''), skip=0, header=T, check.names = FALSE,sep = "\t",stringsAsFactors = F)


set.seed(1234)
km_result <- kmeans(scale(score_matrix[,2:3]),3, nstart = 50)

fviz_cluster(object = km_result, data = score_matrix[,2:3],
             ellipse.type = "euclid", star.plot =T,repel = T,
             geom = ("print"), palette = 'jco', main = "",
             ggtheme = theme_minimal()) +
  theme(axis.title = element_blank())

score_matrix_cluster = cbind(score_matrix, cluster = km_result$cluster)

score_matrix_cluster$cluster = gsub("1","Var",score_matrix_cluster$cluster)
score_matrix_cluster$cluster = gsub("2","Spe",score_matrix_cluster$cluster)
score_matrix_cluster$cluster = gsub("3","Con",score_matrix_cluster$cluster)

write.table(score_matrix_cluster,file = paste(path_out,"spescore_selscore_matrix_kmeans_cluster.txt",sep=''), sep="\t", quote = F, row.names = F, col.names = T)

# plot
library(ggplot2)
score_matrix_cluster <- read.delim(paste(path_out,"spescore_selscore_matrix_kmeans_cluster.txt",sep=''), skip=0, header=T, check.names = FALSE,sep = "\t",stringsAsFactors = F)
score_matrix_cluster$cluster <- factor(score_matrix_cluster$cluster, levels = c("Spe","Var","Con"))

palette <- c("#B5C3C4","#BAD873","#C16EC6")

p <- ggplot(data=score_matrix_cluster,aes(x=SPE_score, y=SEL_score, color = cluster)) +
  geom_point(size = 3,shape = 20) +
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

ggsave(p,filename = paste(path_out,"spescore_selscore_cluster_dotplot.jpg",sep = ''),height = 4,width = 4,dpi = 300)
