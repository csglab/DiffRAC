# DiffRAC: A flexible R function for comparing response proportions in high-throughput sequencing count data

A R function to compare sequencing count ratios using a customized model matrix and DESeq2, from an experimental design, a formula and a read count table. In the example presented here, differential mRNA stability is being inferred from the ratios of exonic and intronic read counts.

# Input data

The input data example shown here represents an experimental design with two variables/treatments/conditions of interest. This function is designed to accomodate any number of experimental variables and the inclusion of interaction terms. 

## Experimental design

A simple `design` data.frame containing the sample names in the row.names, and the experimental variables in the columns with appropriate headers. The variables may be factors or numerical, and the user needs to make sure that the data has the intended class. A minimum of two replicates per condition (or per combination of conditions) is required. Sample names need to be given as row.names. For example:

| samples     | V1  | V2  |
| ----------- | --- | --- |
| cond1_rep1  |	1   | 0   |
| cond1_rep2  | 1   | 0   |
| cond2_rep1  |	1   | 1   |
| cond2_rep2  |	1   | 1   |
| cond3_rep1  |	0   | 0   |
| cond3_rep2  |	0   | 0   |
| cond4_rep1  |	0   | 1   |
| cond4_rep1  | 0   | 1   |

Please avoid the use of "-" and "." in the sample names. Names should also not start with a number.

## Count table

A `counts` data.frame including counts from the two types of reads (e.g. exonic and intronic read counts for the inference of differential mRNA stability) such as the one below. The count column titles must match the sample names in the design table, with the addition suffixes, denoting the two types of read counts (here ".e" and ".i" denote exonic and intronic read counts). The suffixes will define the reads that will be used as the numerator or denominator of the ratio. The suffixes must be separated from the sample name by a period ".". Again, please avoid the use of "-" and "." in the sample names and inside the suffixes. Names should not start with a number. The row names must contain gene IDs.

| Gene_ID | cond1_rep1.e  | cond1_rep1.i  | cond1_rep2.e  | cond1_rep2.i  | cond2_rep1.e  | cond2_rep1.i  | cond2_rep2.e  | cond2_rep2.i  | cond3_rep1.e  | cond3_rep1.i  | cond3_rep2.e  | cond3_rep2.i  | cond4_rep1.e  | cond4_rep1.i  | cond4_rep2.e  | cond4_rep2.i  |
| ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 22746   | 29  | 3   | 27  | 2   | 47  | 3   | 37  | 2   | 52  | 26  | 50  | 25  | 70  | 26  | 60  | 25  |
| 235623  | 122 | 92  | 90  | 18  | 299 | 45  | 454 | 177 | 145 | 115 | 113 | 41  | 322 | 68  | 477 | 200 |
| 238690  | 18  | 14  | 6   | 8   | 71  | 22  | 60  | 16  | 41  | 37  | 29  | 31  | 94  | 45  | 83  | 39  |
| 330369  | 5   | 35  | 4   | 17  | 149 | 55  | 276 | 149 | 28  | 58  | 27  | 40  | 172 | 78  | 299 | 172 |
| ...     | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |


## Formula

A `formula` indicating the relationship between the predictor and outcome variables. The variable names must be the same as in the design matrix. For example: 

\~ V1 + V2 + V1:V2

## Count ratio numerator and denominator

The type of read counts that represent the numerator and denominator of the ratio are given to DiffRAC as the `num` and  `denom` parameters. These strings must be the identical to the two suffixes in the count table. In the current example num = "e" and denom = "i".

## Sample-specific or condition-specific analyses

`mode="sample"` should be used, except in the case of large sample sizes. The `mode="condition"` parameter can then be used, to consider condition-specific changes instead of sample-specific changes for the read type that is the ratio denominator. This will significantly decrease the run time.

# The DiffRAC function

```r

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

```

# Usage

```r
DiffRAC_res <- DiffRAC(design, counts, formula, num, denom, mode)
```

# Output

DiffRAC returns a list with two elements:

The first element is the customized model matrix created by DiffRAC.

```r
DiffRAC_res$model_mat
```

The second element is a DESeq `dds` object (output from DESeqDataSetFromMatrix and DESeq functions). The user can then use their own contrasts and follow the DESeq2 documentation to get differential estimates (as a log2 fold-change) and identify differential events (such as differentially stabilized genes in this example).

```r
DiffRAC_res$dds
 ```
Please note that in case of factor variables, the variable names will be different from the inputed design and formula to the resulting dds object. For example, for a variable V1 with levels 0 and 1, a column V1 will be returned. For a variable V1 with levels "control" and "treatment", a column V1treatment will be returned. In the case of a variable V1 with levels "a", "b", and "c", then the design will contain V1b and V1c columns, each representing the change in stability relative to reference level "a". For numeric variables, the name will remain the same.

Differential events be identified by filtering for padj < 0.1, for example.

# Example

An example dataset is provided at `./examples`. Move into this folder, and load the read counts:

```r
counts <- as.matrix(read.table("test_counts.txt", header = T, row.names=1))
```

Then load the experimental design:

```r
design <- as.data.frame(read.table("test_design.txt", header = T, row.names=1))
```

Define the formula:

```r
formula <- as.formula(~ V1)
```

Define the numerator and denominator of the count ratio:

```r
num <- "e"
```

```r
denom <- "i"
```

Define the mode ("sample" should be used, except for large sample sizes):

```r
mode <- "sample"
```

And run DiffRAC:

```r
source("../DiffRAC.R")
DiffRAC_res <- DiffRAC(design, counts, formula, num, denom, mode)
```

To get differential estimates (differential stability estimates here), for example in the condition V1 level 1 vs 0:

```r
res <- as.data.frame(results(DiffRAC_res$dds, name="V1"))
head(res)
```
