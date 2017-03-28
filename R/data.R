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
