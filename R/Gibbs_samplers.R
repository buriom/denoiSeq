# _____________________________________________________________________________

# Implementation of the Gibbs sampling algorithm to sample from the
# joint posterior distribution of N, p and f.
# _____________________________________________________________________________


# To implement a one step Gibbs sampler
gibbsampling <- function(counts, counts_A, counts_B, parms, propotns, stepsize,
                         m, n_A, n_B) {
  # To update p's
  parms$p <- MH_p(c(rep(parms$N_A, n_A), rep(parms$N_B, n_B)), parms$f,
                  parms$p, rep(propotns, each = m), as.vector(counts), stepsize$stepsize_p)
  # To update f's
  parms$f <- MH_f(c(rep(parms$N_A, n_A), rep(parms$N_B, n_B)), parms$f,
                  parms$p, rep(propotns, each = m), as.vector(counts), stepsize$stepsize_f)
  # To update N's
  parms$N_A <- as.vector(mapply(MH_N, split(parms$N_A, 1:m), parms$f,
                                parms$p, split(rep(propotns[1:5], each = m), 1:m), split(counts_A,
                                1:m), split(stepsize$stepsizeN_A[1:m], 1:m)))
  parms$N_B <- as.vector(mapply(MH_N, split(parms$N_B, 1:m), parms$f,
                                parms$p, split(rep(propotns[6:10], each = m), 1:m), split(counts_B,
                                1:m), split(stepsize$stepsizeN_B[1:m], 1:m)))
  return(parms)
}

# To implement a one step Gibbs sampler, while tuning the step size
gibbsampling2 <- function(counts, counts_A, counts_B, parms, propotns,
                          stepsize, m, n_A, n_B) {
  # To update p's
  results_p <- MH_p2(c(rep(parms$N_A, n_A), rep(parms$N_B, n_B)), parms$f,
                     parms$p, rep(propotns, each = m), as.vector(counts), stepsize$stepsize_p)
  parms$p <- results_p[1]
  stepsize$stepsize_p <- results_p[2]
  # To update f's
  results_f <- MH_f2(c(rep(parms$N_A, n_A), rep(parms$N_B, n_B)), parms$f,
                     parms$p, rep(propotns, each = m), as.vector(counts), stepsize$stepsize_f)
  parms$f <- results_f[1]
  stepsize$stepsize_f <- results_f[2]
  # To update N's
  results_NA <- mapply(MH_N2, split(parms$N_A, 1:m), parms$f, parms$p,
                       split(rep(propotns[1:5], each = m), 1:m), split(counts_A, 1:m),
                       split(stepsize$stepsizeN_A[1:m], 1:m))
  parms$N_A <- results_NA[1, ]
  stepsize$stepsizeN_A <- results_NA[2, ]
  results_NB <- mapply(MH_N2, split(parms$N_B, 1:m), parms$f, parms$p,
                       split(rep(propotns[6:10], each = m), 1:m), split(counts_B, 1:m),
                       split(stepsize$stepsizeN_B[1:m], 1:m))
  parms$N_B <- results_NB[1, ]
  stepsize$stepsizeN_B <- results_NB[2, ]
  return(list(parms = parms, stepsize = stepsize))
}
