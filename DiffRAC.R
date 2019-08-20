# Infers the change in the response ratios using the customized model matrix and DESeq2, from the experimental design, the read count matrix, and the formula provided by the user
# Returns the customized model matrix and a DESeq dds object

DiffRAC <- function(design, counts, formula, num, denom)
{
  # Load the libraries
  library(DESeq2)
  library(plyr)
  
  # Inspect the design to require at least two replicates per condition
  if (sum(count(design, vars = colnames(design))$freq < 2) > 0)
  {
    stop("Error: At least two replicates per condition are required")
  }
  
  # Define the number of samples
  n <- nrow(design)
  
  # Create the design matrix
  ident <- diag(n) 
  design_mat <- rbind(ident, ident)
  readTypeB <- c(rep(0, nrow(ident)), rep(1, nrow(ident))) 
  design_mat <- cbind(design_mat, readTypeB)
  model_mat <- model.matrix(formula, data = model.frame(formula, data = design))
  model_mat <- model_mat[, -1, drop=F]
  model_mat <- rbind(matrix(0, nrow = nrow(model_mat), ncol = ncol(model_mat)), model_mat)
  design_mat <- cbind(design_mat, model_mat)
  colnames(design_mat)[1:n] <- paste("s", (1:n), sep="")
  design_mat[,1] <- 1
  colnames(design_mat)[1] <- "(Intercept)"
  row.names(design_mat) <- c(paste(row.names(design), ".", denom, sep=""), paste(row.names(design), ".", num, sep=""))
  
  # Reorder the count table
  count_input <- data.frame(counts[, row.names(design_mat)], row.names = row.names(counts))
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count_input,
                                colData = design_mat,
                                design = design_mat)
  dds <- DESeq(dds, full=design_mat, betaPrior = F)
  return(dds)
}
