# DiffRAC: A flexible R framework for comparing response proportions in high-throughput sequencing count data

A R framework to compare sequencing count ratios using a customized model matrix and DESeq2, from an experimental design, a formula and a read count tables. 

An example analysis can be found in the *examples* directory. 

# Inputs 

## design

The design data frame. Each row is one sample, and each column is an experimental variable. Sample names should be indicated as row.names, and the experimental variables as column names. Sample names have to match the column names in the count matrices (but not necessarily in the same order). The variables may be factors or numerical, and the user needs to make sure that the data has the intended class. There must be no NA values in the design. A minimum of two replicates per condition (or per combination of conditions) for the variable of interest should be used. For example:

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

Another valid example:

| samples | cell_type      |                                
| ------- | -------------- |                              
| s1      |	cellLineA      |                          
| s2      |	cellLineA      |                            
| s3      |	cellLineB      |                         
| s4      |	cellLineB      |      

Which is equivalent to:

| samples | cellLineA |                                
| ------- | --------- |                              
| s1      |	0         |                          
| s2      |	0         |                            
| s3      |	1         |                         
| s4      |	1         |

Please avoid the use of "-" and "." in the sample names. Names should also not start with a number.

## formula

The model formula for samples

A `formula` indicating the relationship between the predictor and outcome variables. The variable names must be the same as in the design matrix. For example: 

\~ V1 + V2 + V1:V2


## count_num

The count matrix or count data frame for the numerator type (exonic reads)

user would enter two different matrices for the exonic or intronic reads, but the rows and columns have to be in the same order, and the column names would have to match the sample names in the design matrix (although not necessarily with the same order)

A `counts` data.frame including counts from the two types of reads (e.g. exonic and intronic read counts for the inference of differential mRNA stability) such as the one below. The count column titles must match the sample names in the design table, with the addition suffixes, denoting the two types of read counts (here ".e" and ".i" denote exonic and intronic read counts). The suffixes will define the reads that will be used as the numerator or denominator of the ratio. The suffixes must be separated from the sample name by a period ".". Again, please avoid the use of "-" and "." in the sample names and inside the suffixes. Names should not start with a number. The row names must contain gene IDs.

| Gene_ID | cond1_rep1.e  | cond1_rep1.i  | cond1_rep2.e  | cond1_rep2.i  | cond2_rep1.e  | cond2_rep1.i  | cond2_rep2.e  | cond2_rep2.i  | cond3_rep1.e  | cond3_rep1.i  | cond3_rep2.e  | cond3_rep2.i  | cond4_rep1.e  | cond4_rep1.i  | cond4_rep2.e  | cond4_rep2.i  |
| ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 22746   | 29  | 3   | 27  | 2   | 47  | 3   | 37  | 2   | 52  | 26  | 50  | 25  | 70  | 26  | 60  | 25  |
| 235623  | 122 | 92  | 90  | 18  | 299 | 45  | 454 | 177 | 145 | 115 | 113 | 41  | 322 | 68  | 477 | 200 |
| 238690  | 18  | 14  | 6   | 8   | 71  | 22  | 60  | 16  | 41  | 37  | 29  | 31  | 94  | 45  | 83  | 39  |
| 330369  | 5   | 35  | 4   | 17  | 149 | 55  | 276 | 149 | 28  | 58  | 27  | 40  | 172 | 78  | 299 | 172 |
| ...     | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

## count_denom

The count matrix or count data frame for the denominator type (intronic reads)

# The read types to be used as numerator and denominator of the count ratio must be specified using their suffixes
## Count ratio numerator and denominator

The type of read counts that represent the numerator and denominator of the ratio are given to DiffRAC as the `num` and  `denom` parameters. These strings must be the identical to the two suffixes in the count table. In the current example num = "e" and denom = "i".

## mode

Optionally, for large sample sizes, a condition-specific analysis can be performed, instead of a sample-specific investigation.

Either "condition" or "sample"

`mode="sample"` should be used, except in the case of large sample sizes. The `mode="condition"` parameter can then be used, to consider condition-specific changes instead of sample-specific changes for the read type that is the ratio denominator. This will significantly decrease the run time.

# bias: The "bias" constant

optionally perform an optimization to obtain the bias term
alternatively, the bias term can be supplied by the user (or be left as 1, which is the default).

# The DiffRAC function



# Usage

```r
DiffRAC_res <- DiffRAC(design, counts, formula, num, denom, mode)
```

# Output

Returns the customized model matrix and a DESeq dds object

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
