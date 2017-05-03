# _____________________________________________________________________________

# Determining differential expression
# _____________________________________________________________________________

# To obtain parameter values sampled at each
# step of Gibbs sampling.
getN_A <- function(step,RD){
  rezult <- RD@output
  return(rezult[[1]][[step]]$N_A)
}
getN_B <- function(step,RD){
  rezult <- RD@output
  return(rezult[[1]][[step]]$N_B)
}
getp <- function(step,RD){
  rezult <- RD@output
  return(rezult[[1]][[step]]$p)
}
getf <- function(step,RD){
  rezult <- RD@output
  return(rezult[[1]][[step]]$f)
}


#' Calculate the test statistic.
#'
#' Given two samples of a parameter N_i from each of the two conditions,
#'  it calculates the area of the distribution of the log2 transformed differences that falls within the region of practical
#'  equivalance (ROPE) defaulted to [-0.5,0.5].
#'
#' @param N_A A vector of samples of N_i for condition A.
#' @param N_B A vector of samples of N_i for condition B.
#' @param rope-limit A float that delimits the range of the region of practical
#'  equivalance, ROPE. A default value of 0.5 is used.
#' @return A value between 0 and 1 that can be interpreted as the probability that there is no difference.

#' @export
DEstat <- function(N_A, N_B, rope_limit=0.5) {
    # log2 difference of samples of the same
    # parameter across the 2 conditions
    dif <- log2(N_A) - log2(N_B)
    # region of practical equivalence
    rope <- sum(dif > -rope_limit & dif < rope_limit)
    return(rope/length(N_A))
}

#' Differential Expression Analysis based on a bottom up model
#'
#' To extract the results of Bayesian inference returned by denoiSeq function.
#'
#' @param RD A readsData object with the output slot filled.
#' @param steps  An integer representing the number of iterations.
#' @param burnin An integer for the number of iterations to be considered as burn in values
#' @param rope-limit A float that delimits the range of the region of practical
#'  equivalance, ROPE. A default value of 0.5 is used.
#'
#' @return A dataframe with 3 columns namely; the log2 fold change (log2FC), the standard error of the log2 fold change (lgfcSE) and the test static (stat).

#' @export
results <- function(RD, steps, burnin = floor(steps/3),  rope_limit = 0.5) {
  N_Asamples <- t(sapply(1:steps,getN_A,RD=RD))
  N_Bsamples <- t(sapply(1:steps,getN_B,RD=RD))
  N_Asamples <- tail(N_Asamples, steps - burnin)
  N_Bsamples <- tail(N_Bsamples, steps - burnin)
  m <- ncol(N_Asamples)
  values <- mapply(DEstat, split(t(N_Asamples), 1:m), split(t(N_Bsamples), 1:m), rope_limit = rope_limit )
  lfc_mat <-   N_Bsamples/N_Asamples
  lfc_mean <- apply(lfc_mat,2,mean)
  lfc_SE <- apply(lfc_mat,2,sd)/sqrt(nrow(lfc_mat))
  df <- data.frame("log2FC(B vs A)" = lfc_mean,lfcSE = lfc_SE, stat = values)
  rownames(df) <- RD@geneNames
  return(df)
}

