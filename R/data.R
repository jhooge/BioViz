#' iPSC derived motor neuron cultures from C9ORF72 carriers
#'
#' This dataset includes -log10(pvalues) and -log2(foldchanges)
#' of genes associated with Amytrophic Lateral Sclerosis. It has been
#' derived by Stephen Turner from the GEO GSE52202 dataset by aligning it to hg19
#' and analysing it with DESeq2.
#'
#' @format A data frame with 16406 rows and 4 variables:
#' \describe{
#'   \item{id}{HGNC Gene Symbol}
#'   \item{fc}{Log2 Foldchange}
#'   \item{pvalue}{Adjusted pvalues}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse52202}
#' @source \url{https://gist.github.com/stephenturner/806e31fce55a8b7175af}
"iPSC"

#' Classification Data with Five Classes
#' 
#' This dataset includes benchmark data of five classes for the benchmarking of 
#' classification models.
#' It includes 250 samples and 15 variables and an additional label
#' vector in the first column for supervised classification. 10 of these variables (V2, 
#' V2, ..., V10) are informative with respect to the class labels and 5
#' are noise variables (U1, U2, ..., U5) which do not correlate with the label vector.
#' Each of the variables in this dataset are drawn normal distributed 
#' with different means and a fixed standard deviation of 1.
#'
#' @format A data frame with 250 rows and 16 columns:
#' \describe{
#'   \item{Class}{Label or Response Vector}
#'   \item{V2, ..., V10}{Informative Variables}
#'   \item{U1, ... , U5}{Uninformative Variables}
#' }
"fiveClass"