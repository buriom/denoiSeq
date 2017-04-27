# _____________________________________________________________________________

# Calling the Gibbs sampling algorithm for
# multiple steps.
# _____________________________________________________________________________


# To calculate the size factors used in
# normalizing the counts
size_factors <- function(counts) {
  # To eliminate rows with zeros.
    counts <- counts[-(which(apply(counts,1,prod)==0)),]
    geomean <- (apply(counts, 1, prod))^(1/ncol(counts))
    s_j <- apply(counts/geomean, 2, median)
    return(s_j/max(s_j))
}

#' To estimate model parameters for PCR sequencing using Bayesian Inference.
#'
#' Uses the Gibb's sampling algorithm to sample from the joint distribution of the model parameters.
#'
#' The counts in each column are used to estimate the size factors which are
#' in turn used to normalise the counts in each column. The function then uses the
#' rows of the matrix to estimate parameter N_i for each transcript and uses the
#' entire dataset (combined from both conditions)  to estimate \code{p} and
#' \code{f}. The result is an estimate for each of the \code{m} \code{N_i}
#' parameters for each condition and estimates for \code{p} and \code{f}.
#'
#' @param counts An m by n integer matrix of counts from two experimental conditions
#' @param steps  An integer representing the number of iterations.
#' @param n_A  An integer the number of columns of counts in condition A.
#' @param tuning.steps An integer representing the number of iterations to be
#'   used for tuning the step sizes.
#' @param init.vals A list containing initial values for all N's (in both
#'   conditions), p and f.
#' @param step.size A list containing step sizes values for all the N's, p and
#'   f.

#' @return A list of 2 lists,i.e samples which contains parameter samples
#' for each of  \code{N_i},\code{p} and \code{f} for each iteration, and
#'  stepsize which  is a list of the tuned step sizes.

#' @export
nstep_gibbsampling <- function(counts,steps, n_A = ncol(counts)/2,
      tuning.steps = floor(steps/3),
    init.vals = list(N_A = rep(1, m), N_B = rep(1,
        m), p = 0.0001, f = 0.01), step.size = list(stepsizeN_A = rep(1,
        m), stepsizeN_B = rep(1, m), stepsize_p = 5e+05, stepsize_f = (1000))) {
    m = nrow(counts)
    n_B = ncol(counts) - n_A
    counts_A = counts[, 1:n_A]
    counts_B = counts[,(n_A + 1):ncol(counts)]
    # Calculating the size factors used in
    # normalizing the counts
    propotns <- size_factors(counts)
    #if(!(is.null(rownames(counts)))){
    names(init.vals$N_A) <- names(init.vals$N_B) <- rownames(counts)
      #}
    parms <- init.vals
    parm_samples <- list(parms)
    stepsize_vectr <- list(step.size)
    #tuning
    for (i in 2:tuning.steps) {
        results <- gibbsampling2(counts, counts_A,
            counts_B, parm_samples[[i - 1]], propotns,
            stepsize_vectr[[i - 1]], m, n_A, n_B)
        parm_samples[[i]] <- results$parms
        stepsize_vectr[[i]] <- results$stepsize
    }
    #setting tuned stepsizes
    step.size <- results$stepsize
    for (i in (tuning.steps + 1):steps) {
        parm_samples[[i]] <- gibbsampling(counts,
            counts_A, counts_B, parm_samples[[i -
                1]], propotns, step.size, m, n_A,
            n_B)
    }
    return(list(samples = parm_samples, stepsize = stepsize_vectr))
}



theta0 <- list(N_A=rep(1,71),N_A=rep(1,71), p=0.0001, f=0.01)

#step sizes for each of the parameters
Ap=rep(5e5,1); Af=(1e3)

A <- list(stepsizeN_A=rep(1,71),stepsizeN_B=rep(1,71),stepsize_p=Ap,stepsize_f=Af)

ptm <- proc.time()
steps <- 50;tune <- 15
rezult <- nstep_gibbsampling (ERCC,steps,tuning.steps = tune)
proc.time() - ptm
