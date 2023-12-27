options(repos = list(CRAN="http://cran.rstudio.com/"))
install.packages(c('dplyr','data.table','ggplot2'), repos = "http://cran.us.r-project.org")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", type = "source", checkBuilt = TRUE)

library(dplyr)
library(DESeq2)
library(ggplot2)
library(data.table)

# prepare the data for the model:
DE_prepare <- function(input_file_path){
  df <- fread(input_file_path, data.table=FALSE) 
  s <- colnames(df)[-1]
  # parse specifically to the file pattern:
  samples <- gsub("[[:digit:]]", "", s) # the name of the sample without the rep number
  replicates <- gsub('con', '',gsub('treat', '', s)) # the replicate number
  sample_df <- data.frame(sample=samples,replicate=replicates)
  row.names(sample_df) <- colnames(df)[-1]
  counts <- as.matrix(df[,2:ncol(df)])
  rownames(counts) <- df$gene_symbol
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sample_df,
                                design = ~ sample)
  
  return(dds)
}


get_DE_comparisons <- function(input_file_path, comparisons, outdir){
  input_file_path <- as.character(input_file_path)
  results_file <- fread(input_file_path, data.table = FALSE) 
  colnames(results_file)[1] <- "gene_symbol"
  dds <- DE_prepare(input_file_path)
  comp_results <- list()
  dds <- DESeq(dds, betaPrior=TRUE)
  for (comp in comparisons){
    comp_samples = strsplit(comp,"_vs_")[[1]]
    nominator <- comp_samples[1]
    denominator  <- comp_samples[2]
    # how to define specific comparison:
    test <- c("sample", nominator, denominator)
    comp_name <- comp 
    res <- results(dds, contrast=test)
    normalized_counts_df <- counts(dds, normalized = TRUE)
    res_df <- as.data.frame(res)
    res_df$gene_sym <- rownames(res)
    
    comp_results[[comp_name]] <- res_df
    }
  all_comparisons <- bind_rows(comp_results, .id = "comparison")
  rownames(all_comparisons) <- results_file$gene_symbol
  all_comparisons <- tibble::rownames_to_column(all_comparisons, "gene_symbol")
  file_name <- "results_de.csv" 
  out_file_de <- paste(outdir,file_name,sep= "/")
  write.csv(all_comparisons, out_file_de, quote = F, row.names = F)
  
  #plot a volcano plot of the data:
  ggplot(data=all_comparisons, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
  
  # filter genes with fold change above 2 and p-value below 0.05. 
  # Write a table of genes to a file. The table will include 3 columns:
  # gene name, fold change, p value.
  
  all_comparisons$FC <- with(all_comparisons, 2^all_comparisons$log2FoldChange)
  filtered_genes <- subset(all_comparisons, abs(FC)>2 & pvalue<0.05)
  filtered_genes <- tibble::rownames_to_column(filtered_genes, "gene_name")
  file_name <- "filtered_genes_de.csv" 
  out_file_filtered <- paste(outdir,file_name,sep= "/")
  write.csv(filtered_genes[ , c("gene_name", "FC", "pvalue")]  , out_file_filtered, quote = F, row.names = F)
}

# Rscript DEseq_analysis.R  input_file_path comparisons
args = commandArgs(trailingOnly=TRUE)
input_file_path <- as.character(args[1])
comparisons <- as.character(args[2]) # 'treat_vs_con'
outdir <- dirname(normalizePath(input_file_path))
get_DE_comparisons(input_file_path, comparisons, outdir)

