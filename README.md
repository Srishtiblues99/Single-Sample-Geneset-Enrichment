# Single-Sample-Geneset-Enrichment method is an extension of the GSEA method, working at the level of a single sample rather than a sample population as in the original GSEA application. The score derived from ssGSEA reflects the degree to which the input gene signature is coordinately up- or downregulated within a sample.
Gene Set Enrichment Analysis (GSEA) is a method for interpreting gene expression data. The method derives its power by focusing on gene sets, that is, groups of genes that share common biological function, chromosomal location, or regulation. We demonstrate how GSEA yields insights into several cancer-related data sets. Notably, where single-gene analysis finds little similarity between two independent studies of patient survival in a specific cancer, GSEA reveals many biological pathways in common.
# Load the packages to read, average and plot the data
library(matrixStats)
library(circlize)

library(ComplexHeatmap)
library(data.table)

# SSGSEA function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X) 
  gene_sets = lapply(gene_sets, function(genes)
  {which(row_names %in% genes)})
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    
    es_sample = sapply(gene_sets, function(gene_set_idx) {
     
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
   
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  if (norm) es = es / diff(range(es))
  
  
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

# Read the log2 normalized data into an object for further analysis
 data= readRDS("log_ttest.RDS") 
data[1:5,1:5] 

# Read the gene set matrix into an object
data = as.matrix(data)
gene_set= read.csv("C:/Users/shristi/Downloads/markers2Sep.csv")
head(gene_set)
gene_sets= as.list(as.data.frame(gene_set))
print("genes get ready")
res= ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)

# Transpose the res object for viewing purposes
res1 = t(res)
head(res1)
mat = as.matrix(res1)

# Zscore the ssgsea output for comparative analysis
for (i in 1:nrow(mat)){
  
  vec = as.numeric(mat[i,])
  mat[i,1:ncol(mat)] = (vec-mean(vec))/sd(vec)
  
}

# Heatmap function
Heatmap(t(mat),col=colorRamp2(c(-2,0,2),c("purple","white","light pink")))
