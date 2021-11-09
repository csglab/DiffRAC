# Infers the change in the response ratios using the customized model matrix and DESeq2, from the experimental design, the read count matrix, and the formula provided by the user
# The read types to be used as numerator and denominator of the count ratio must be specified using their suffixes
# Optionally, for large sample sizes, a condition-specific analysis can be performed, instead of a sample-specific investigation.
# Returns the customized model matrix and a DESeq dds object

# formula: The model formula for samples
# design: The design data frame. Each row is one sample, and each column is a variable
# count_num: The count matrix or count data frame for the numerator type (exonic reads)
# count_denom: The count matrix or count data frame for the denominator type (intronic reads)
# mode: Either "condition" or "sample"
# bias: The "bias" constant
DiffRAC.initialize <- function(formula, design, counts_num, counts_denom, mode, bias=1)
{
  # Load the libraries
  library(DESeq2)
  library(plyr)
  
  # check that the count tables have the same rows and columns
  if( ncol(counts_num) != ncol(counts_denom) |
      sum( colnames(counts_num) != colnames(counts_denom) ) > 0 )
  {
    stop("Incompatible columns in count tables.")
  }
  if( nrow(counts_num) != nrow(counts_denom) |
      sum( rownames(counts_num) != rownames(counts_denom) ) > 0 )
  {
    stop("Incompatible rows in count tables.")
  }
  
  ## ensure that the design data frame has the same rows as the columns of the count tables
  # choose only rows that correspond to a column in the count matrix
  design <- design[ rownames(design) %in% colnames(counts_num), , drop=F]
  if( nrow(design) != ncol(counts_num) )
  {
    stop("Incompatible or incomplete design matrix: some samples with count data were not found in design.")
  }
  # put the rows in the same order as columns of the count matrix
  design <- design[ match( colnames(counts_num), rownames(design) ) , , drop=F]
  
  # First, create the experiment design matrix (regardless of intron/exon status)
  model_mat <- model.matrix(
    formula,
    data = model.frame(formula, data = design)) # will drop rows with NAs
  # If the design matrix includes an intercept, drop it
  ## Note: I changed it from model_mat[,-1], because the user might specify a formula that leads to no intercept, or a model that leads to another variable that is entirely constant and therfore co-linear with intercept
  keep <- !( apply(model_mat,2,sd)==0 ) # any columns that have no variance (including intercept) will be dropped
  model_mat <- model_mat[, keep, drop=F]

  # Define the number of samples
  n <- nrow(model_mat)
  Intercept <- rep(1, 2*n)
  Numerator <- c( rep(0, n), rep(1, n) )

  # a zero matrix with the same dimensions as the model matrix
  Zero <- matrix( 0, nrow=nrow(model_mat), ncol=ncol(model_mat) )
  colnames(Zero) <- paste0(colnames(model_mat),":Ratio")
  colnames(model_mat) <- paste0(colnames(model_mat),":Denominator")
  
  if( mode == "condition" ) # Condition-specific mode
  {
    design_mat <- cbind(
      Intercept,
      rbind( model_mat, model_mat * bias ),
      Numerator,
      rbind( Zero, model_mat ) )
  } else if( mode == "sample" ) # Sample-specific mode
  {
    ident <- diag(n)
    colnames(ident) <- rownames(model_mat)
    
    design_mat <- cbind(
      Intercept,
      rbind( ident[,-1, drop=F], ident[,-1, drop=F] * bias ),
      Numerator,
      rbind( Zero, model_mat ) )
  } else { stop("Mode not recognized.") }
    
  # fix the column and row names
  colnames(design_mat)[1] <- "(Intercept)"
  rownames(design_mat) <- c(
    paste0(rownames(model_mat),":Denominator"),
    paste0(rownames(model_mat),":Numerator") )
    
  # Create the merged count table
  counts <- cbind( counts_denom, counts_num )
  colnames(counts) <- rownames(design_mat)
  
  return( list(design_mat=design_mat,counts=counts) )
}

# Modifies the bias term for a design matrix that is already constructed
# ratio: the ratio of the new bias constant to the previous bias constant
DiffRAC.modifyBias <- function( design_mat, mode, ratio )
{
  n <- nrow(design_mat)/2
  
  if( mode=="condition" )
  {
    design_mat[ (n+1):(n*2), 2:(ncol(design_mat)/2) ] <-
      design_mat[ (n+1):(n*2), 2:(ncol(design_mat)/2) ] * ratio
  } else if( mode=="sample" )
  {
    design_mat[ (n+1):(n*2), 2:n ] <-
      design_mat[ (n+1):(n*2), 2:n ] * ratio
  } else { stop("Mode not recognized.") }
  
  return(design_mat)
}

DiffRAC <- function( formula, design, counts_num, counts_denom, mode="condition", bias=1, optimizeBias=F )
{

  cat("\nInitializing DiffRAC framework...\n")
  # first, initialize the design matrix and counts, with an initial bias constant of 1
  drc <- DiffRAC.initialize( formula, design, counts_num, counts_denom, mode )

  # create a design matrix that does not consider the effect of experimental variables on stability
  n <- nrow(drc$design_mat)/2
  if( mode=="condition" )
  {
    trn_design_mat <- drc$design_mat[ , 1:(ncol(drc$design_mat)/2+1) ]
  } else if( mode=="sample" )
  {
    trn_design_mat <- drc$design_mat[ , 1:(n+1) ]
  } else { stop("Mode not recognized.") }
  
  
  cat("\nEstimating size factors and dispersions...\n")
  # initialize DESeq2, with a dummy design that simply takes into account the read types
  dds <- DESeqDataSetFromMatrix(countData = drc$counts,
                                colData = drc$design_mat,
                                design = ~ Numerator )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds,modelMatrix=trn_design_mat) # when estimating dispersions, consider the effect of experimental variables (or the sample effect) on transcription
  
  if( optimizeBias ) # requested to optimize the bias constant
  {
    ## Optimize bias
    # This works by fitting a model that does not have the "Numerator" model variables, with varying biase constants
    # In other words, we first try to find a bias constant that explains as much of the data as possible assuming only transcriptional regulation
    # If "sample" mode, the transcriptional regulation means sample-specific effects, and if "condition" mode, it means effect of experimental variables on transcription
    # Once the bias constant is optimized, then, the additional effect of post-transcriptional regulation is added to the model and is tested given that bias constant
    
    cat("\nOptimizing the bias constant...\n")
    
    optimLRT <- function(x)
    {
      full <- trn_design_mat # the model that contains only the transcriptional effects
      full[ (n+1):(n*2) , -c(1,ncol(trn_design_mat)) ] <-
        full[ (n+1):(n*2) , -c(1,ncol(trn_design_mat)) ] * x # modify the bias constant of this model
      reduced <- full[,c(1,ncol(trn_design_mat))] # only the intercept and the read-type variables
      
      lrt <- nbinomLRT( dds, full = full, reduced = reduced )
      
      stat <- sum( as.data.frame(results(lrt))$stat, na.rm = T )
      
      cat(paste0(x," : ",stat,"\n"))
      
      return(-stat)
    }
    
    optBias <- optimize(optimLRT,interval=c(0,1),tol = 0.001)
    
    bias <- optBias$minimum # the bias constant is now the value provided to the model
  }
  
  # modify the design matrix based on the bias term (either optimized or provided as argument)
  cat(paste0("The bias constant is ",bias,"\n"))
  if (bias != 1) { # no need to re-estimate dispersions if bias == 1
    drc$design_mat <- DiffRAC.modifyBias(drc$design_mat,mode,bias)
    cat("\nRe-estimating dispersion...\n")
    dds <- estimateDispersions(dds,modelMatrix=drc$design_mat)
  }
  
  cat("\nFitting model parameters...\n")

  # Run a Walt test with the bias-modified design matrix
  dds <- nbinomWaldTest(dds, modelMatrix=drc$design_mat, betaPrior = F)
  res <- list(
    model_mat=drc$design_mat,
    counts=drc$counts,
    dds=dds)
  
  return(res)
}

