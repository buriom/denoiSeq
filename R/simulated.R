#' simulated data
#'
#'  A dataset containing simulated data based on parameter values N = 1,2,...,50
#'  , p = 0.0017  and f =0.1,0.2,...,0.5. The values of N were repeated 15 times
#'   to generate  750 genes.  This dataset  contains 750 observational genes
#'   with 5 experimental samples for each condition,  summarised as a 750 by 10
#'   integer matrix. The first 428 genes are not differentially expressed
#'   between the two conditions whereas the last 322 genes are. The gene counts
#'    were generated in accordance to the probability distribution derived in
#'    Ndifon et al.
#' @format A matrix with 750 rows and 10 columns:
#' \describe{
#'   \item{Conditions}{There are 5 columns for each of the conditions A and B.}
#'   \item{Transcripts}{There are 750 distinct genes without  names.}
#' }

"simdat"
