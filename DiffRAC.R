# Infers the change in the response ratios using the customized model matrix and DESeq2, from the experimental design, the read count matrix, and the formula provided by the user
# The read types to be used as numerator and denominator of the count ratio must be specified using their suffixes
# Optionally, for large sample sizes, a condition-specific analysis can be performed, instead of a sample-specific investigation.
# Returns the customized model matrix and a DESeq dds object

DiffRAC <- function(design, counts, formula, num, denom, mode)
{
  # Load the libraries
  library(DESeq2)
  library(plyr)
  
  # Inspect the design to require at least two replicates per condition
  if (sum(plyr::count(design, vars = colnames(design))$freq < 2) > 0)
  {
    stop("Error: At least two replicates per condition are required")
  }
  
  # Define the number of samples
  n <- nrow(design)
  
  # Create the design matrix, for the sample-specific mode
  ident <- diag(n) 
  design_mat <- rbind(ident, ident)
  readTypeNum <- c(rep(0, nrow(ident)), rep(1, nrow(ident))) 
  design_mat <- cbind(design_mat, readTypeNum)
  model_mat <- model.matrix(formula, data = model.frame(formula, data = design))
  model_mat <- model_mat[, -1, drop=F]
  ext_model_mat <- rbind(matrix(0, nrow = nrow(model_mat), ncol = ncol(model_mat)), model_mat)
  design_mat <- cbind(design_mat, ext_model_mat)
  colnames(design_mat)[1:n] <- paste("s", (1:n), sep="")
  design_mat[,1] <- 1
  colnames(design_mat)[1] <- "(Intercept)"
  row.names(design_mat) <- c(paste(row.names(design), ".", denom, sep=""), paste(row.names(design), ".", num, sep=""))
  
  if (mode == "condition") # Condition-specific mode
  {
    colnames(model_mat) <- paste(colnames(model_mat), "_cond", sep="")
    design_mat_cond <- design_mat[, -(2:n)]
    design_mat_cond <- cbind(design_mat[,1, drop=F], rbind(model_mat, model_mat), design_mat_cond[, 2:ncol(design_mat_cond)])
    row.names(design_mat_cond) <- row.names(design_mat)
    design_mat <- design_mat_cond
  } else if (mode != "sample")
  {
    stop("Error: Mode not recognized")
  }
  
  if(sum(!(row.names(design_mat) %in% colnames(counts))) > 0)
  {
    stop("Error: The samples in the experimental design and count matrix are not the same")
  }
  
  # Reorder the count table
  count_input <- data.frame(counts[, row.names(design_mat)], row.names = row.names(counts))
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count_input,
                                colData = design_mat,
                                design = design_mat)
  dds <- DESeq(dds, full=design_mat, betaPrior = F)
  res <- list(design_mat, dds)
  names(res) <- c("model_mat", "dds")
  return(res)
}

