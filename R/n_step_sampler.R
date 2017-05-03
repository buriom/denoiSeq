# _____________________________________________________________________________

# Calling the Gibbs sampling algorithm for
# multiple steps.
# _____________________________________________________________________________


# To calculate the size factors used in
# normalizing the counts
size_factors <- function(counts) {
  # To eliminate rows with zeros.
    counts <- counts[apply(counts,1,prod)>0,]
    geomean <- (apply(counts, 1, prod))^(1/ncol(counts))
    s_j <- apply(counts/geomean, 2, median)
    return(s_j/max(s_j))
}

#' To estimate model parameters for PCR sequencing using Bayesian Inference.
#'
#' Uses the Gibb's sampling algorithm to sample from the joint distribution of the model parameters.
#'
#' The counts in each column are used to estimate the size factors which are
#' in turn used to normalise the counts. The function then uses the
#' rows of the matrix to estimate parameter N_i for each transcript and uses the
#' entire dataset (combined from both conditions)  to estimate \code{p} and
#' \code{f}. The result is an estimate for each of the \code{m} \code{N_i}
#' parameters for each condition and estimates for \code{p} and \code{f}.
#'
#' @param RD A readsData object with atleast the counts slot filled.
#' @param steps  An integer representing the number of iterations.
#' @param tuningSteps An integer representing the number of iterations to be
#'   used for tuning the step sizes.
# ' @param init.vals A list containing initial values for all N's (in both
# '   conditions), p and f.
# ' @param step.size A list containing step sizes values for all the N's, p and
# '   f.

#' @return The same readsData object but with the output slot filled by a  list of 2 lists; a list named samples which contains parameter samples
#' for each of  \code{N_i},\code{p} and \code{f}, at each iteration, and a second list called
#'  stepsize which  contains the tuned step sizes.

#' @export
denoiseq <- function(RD,steps, tuningSteps = floor(steps/3)) {
  #unpacking the counts
  counts <- RD@counts
  m = nrow(counts)
  #subsetting to obtain matrices for each condition
  counts_A = counts[, RD@replicates$A]
  n_A <- ncol(counts_A)
  counts_B = counts[,RD@replicates$B]
  n_B <- ncol(counts_B)

  # Calculating the size factors used in
  # normalizing the counts
  propotns <- size_factors(counts)
  #initialising a list for the initial values
  parm_samples <- list(RD@initValues)
  #initialising a list for the step sizes
  stepsize_vectr <- list(RD@stepSizes)
  #tuning
  for (i in 2:tuningSteps) {
    results <- gibbsampling2(counts, counts_A,
                             counts_B, parm_samples[[i - 1]], propotns,
                             stepsize_vectr[[i - 1]], m, n_A, n_B)
    parm_samples[[i]] <- results$parms
    stepsize_vectr[[i]] <- results$stepsize
  }
  #setting tuned stepsizes
  step.size <- results$stepsize
  for (i in (tuningSteps + 1):steps) {
    parm_samples[[i]] <- gibbsampling(counts,
                                      counts_A, counts_B, parm_samples[[i -
                                                                          1]], propotns, step.size, m, n_A,
                                      n_B)
  }
  RD <- setOutput(RD,list(samples = parm_samples, stepsize = stepsize_vectr))
  return(RD)
}

