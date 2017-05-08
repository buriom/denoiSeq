#' simulated data
#'
#'  A dataset containing simulated data based on parameter values N = 1,2,...,50, p = 0.0017
#'  and f = 0.5. The values of N were repeated 15 times to generate  750 genomic events. The
#'  first 375 genomic events are strictly not different while the last 375 are a mixture of both
#'   non different and  different (obtained by reversing
#'   the alignment for the second condition to 50,49, . . .,1, hence values like 25 and 26 are indistinguishable).
#'
#' @format A matrix with 750 rows and 10 columns:
#' \describe{
#'   \item{Conditions}{There are 5 columns for each of the conditions A and B.}
#'   \item{Transcripts}{There are 750 distinct transcripts without  names.}
#' }

"simdat"
