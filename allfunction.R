binaryMatrix <- function(datapath, # The filepath of EP link files for each cell
                         suffix # The suffix following the cell type in the filename
                         ){
  file_ins <- list.files(path = datapath, pattern = paste0('.*',suffix);

  allEP_name <- c()
  celltypes <- c()

  for(i in 1:length(file_ins)){
    EP <- read.delim(paste0(datapath,file_ins[i]), skip=0, header=T, check.names = FALSE, sep = "\t", stringsAsFactors = F)
    EP_name <- unique(paste(EP[,1],EP[,2],sep = '#'))
    allEP_name <- union(EP_name,allEP_name)
    celltypes[i] <- gsub(suffix,"",file_ins[i])
    }

  allEP_name <- sort(unique(allEP_name))

  EPall <- data.frame(matrix(0, nrow = length(allEP_name), ncol = length(celltypes)))
  rownames(EPall) <- allEP_name
  colnames(EPall) <- celltypes

  for(i in 1:length(file_ins)){
    EP <- read.delim(paste0("datapath",file_ins[i]), skip=0, header=T, check.names = FALSE, sep = "\t", stringsAsFactors = F)
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

  return(EP)
}


sepBinaryMatrix <- function(allEP_binary, # The binary matrix of all EP (the output of step1)
                            pathout # The output path of EP binary matrix of each gene
                          ){
  c1 <- ncol(allEP_binary)

  gene_list <- as.character(unique(allEP_binary[,"gene"]))

  for(i in 1:length(gene_list)){
    EP_matrix <- allEP_binary[which(allEP_binary[,"gene"] %in% gene_list[i]),1:(c1-2)]

    # keep multiple enhancer-gene links (enhancer > 1)
    if(nrow(EP_matrix) >= 2){
      write.table(EP_matrix, file=paste0(sep_pathout,gene_list[i],".txt"), sep="\t", quote = F, row.names = T, col.names = T)
    }

	progress_value <- (i/length(file_ins))*100
	if(progress_value %% 5 == 0){print(paste0(progress_value,"%"))}
  }
}



calculateScore <- function(allEP_binary, # The binary matrix of all EP (the output of step1)
                           seppath, # The filepath of EP binary matrix of each gene(the output of step2)
                          ){
  c1 <- ncol(allEP_binary)
  celltypes <- colnames(allEP_binary)[1:c1]
  file_ins <- list.files(path = seppath, pattern = '.*.txt');

  score_matrix <- data.frame(gene = "N", SPE_score = 0, SEL_score = 0, stringsAsFactors = F)

  for(i in 1:length(file_ins)){#i=1
    gene <- gsub(".txt","",file_ins[i])

    each_EP <- read.delim(paste0(seppath,file_ins[i]), skip=0, header=T, check.names = FALSE,sep = "\t",stringsAsFactors = F)
    each_EP_t <- t(each_EP)

    cell_num <- nrow(each_EP_t)
    eh_num <- ncol(each_EP_t)

    spe_score <- mean(cell_num - colSums(each_EP_t) + 1) / cell_num
    sel_score <- sd(rowSums(each_EP_t) / eh_num) / mean(rowSums(each_EP_t) / eh_num)

    score_matrix[i,] <- c(gene, spe_score, sel_score)

    progress_value <- (i/length(file_ins))*100
    if(progress_value %% 5 == 0){print(paste0(progress_value,"%"))}
    }

  return(score_matrix)
}