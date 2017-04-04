# _____________________________________________________________________________

# Calling the Gibbs sampling algorithm for
# multiple steps.
# _____________________________________________________________________________


# To calculate the size factors used in
# normalizing the counts
ratios <- function(counts) {
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
#' @param counts An m by (n+k) matrix of counts from two experimental conditions
#' @param n_A  The number of columns of counts in condition A.
#' @param n_B The number of columns of counts in condition B.
#' @param counts_A An m by n matrix of counts from condition A.
#' @param counts_B An m by k matrix of counts from condition B.
#' @param steps  An integer representing the number of iterations.
#' @param tuning.steps An integer representing the number of iterations to be
#'   used for tuning the step sizes.
#' @param init.vals A list containing initial values for all the N's (for both
#'   conditions), p and f.
#' @param step.size A list containing step sizes values for all the N's, p and
#'   f.

#' @return A list of samples  for each of  \code{N_i} ,\code{p} and \code{f}.

nstep_gibbsampling <- function(counts, n_A = ncol(counts)/2,
    n_B = ncol(counts) - n_A, m = nrow(counts),
    counts_A = counts[, 1:n_A], counts_B = counts[,
        (n_A + 1):ncol(counts)], steps, tuning.steps = floor(steps/3),
    init.vals = list(N_A = rep(1, m), N_B = rep(1,
        m), p = 0.005, f = 0.1), step.size = list(ss_A = rep(1,
        m), ss_B = rep(1, m), ss_p = 5e+05, ss_f = (1000))) {
    # Calculating the size factors used in
    # normalizing the counts
    propotns <- ratios(counts)
    names(init.vals$N_A) <- names(init.vals$N_B) <- rownames(counts)
    parms <- init.vals
    parm_samples <- list(parms)
    stepsize_vectr <- list(step.size)
    for (i in 2:tuning.steps) {
        results <- gibbsampling2(counts, counts_A,
            counts_B, parm_samples[[i - 1]], propotns,
            stepsize_vectr[[i - 1]], m, n_A, n_B)
        parm_samples[[i]] <- results$parms
        stepsize_vectr[[i]] <- results$stepsize
    }
    step.size <- results$stepsize
    for (i in (tuning.steps + 1):steps) {
        parm_samples[[i]] <- gibbsampling(counts,
            counts_A, counts_B, parm_samples[[i -
                1]], propotns, step.size, m, n_A,
            n_B)
    }
    return(list(samples = parm_samples, stepsize = stepsize_vectr))
}
