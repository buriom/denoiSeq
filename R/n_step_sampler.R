# _____________________________________________________________________________

# Calling the Gibbs sampling algorithm for multiple steps.
# _____________________________________________________________________________


# To calculate the size factors used in normalizing the counts
size_factors <- function(counts) {
    # To eliminate rows with zeros.
    counts <- counts[apply(counts, 1, prod) > 0, ]
    geomean <- (apply(counts, 1, prod)) ^ (1 / ncol(counts))
    s_j <- apply(counts / geomean, 2, median)
    return(s_j / max(s_j))
}

#' Differential expression analysis based on a bottom-up model (a superposition of a binomial and negative binomial distribution)
#'
#' The denoiseq function perfoms default analysis by first normalising the counts and then estimating the model
#' parameters using Bayesian inference. Size factors are estimated from count matrix and used for the normalisation.
#' The  Gibb's sampling algorithm is then used to sample from the joint distribution of the model parameters.
#'
#' The denoiSeq package is based on a bottom up model for PCR amplified sequencing developed by Ndifon et al (2012)
#'  which culminates into a  superposition of the binomial and negative binomial as the modelled distribution for the gene counts.
#' The derived distribution has three main parameters, i.e \code{N,p} and \code{f} which represent the initial gene amount  before
#' amplification, the amplification efficiency and the dilution rate, respectively.
#'
#' This package  uses Bayesian inference to
#'  estimate the model parameters. The counts in each column are used to estimate the size factors (Anders and Huber,2010) which are
#' in turn used to normalise the counts. For an \code{m} by \code{n} count matrix, we have
#' \code{2m} \code{N_i} parameters to estimate (\code{m} for the same genes in each of the two conditions) in addition to parameters \code{p}
#' and \code{f}. denoiseq then uses the
#' rows in each condition to estimate parameter N_i for each gene in each condition, and the
#' entire dataset combined from both conditions  to estimate \code{p} and
#' \code{f}. The result is an estimate for each of the \code{m} \code{N_i}
#' parameters for each condition and estimates for \code{p} and \code{f}.
#'
#' For differential expression analysis, the parameters of interest are the \code{N_i}s for each gene, \code{i}.
#' @param RDobject A readsData object with atleast the counts slot filled.
#' @param steps  An integer representing the number of iterations.
#' @param tuningSteps An integer representing the number of iterations to be
#'   used for tuning the step sizes. Defaulted to a third of steps.
#'
#' @return The same readsData object but with the output slot filled by a  list of 2 lists; a list named samples which contains parameter samples
#' for each of  \code{N_i},\code{p} and \code{f}, at each iteration, and a second list called
#'  stepsize which  contains the tuned step sizes.
#' @examples
#' #pre -filtering to remove lowly expressed genes
#' ERCC <- ERCC[rowSums(ERCC)>0,]
#' RD <- new('readsData',counts = ERCC)
#' steps <- 30
#' #30 steps are just for illustration here. Atleast 5000 steps are adequate.
#' BI <- denoiseq(RD,steps)
#'
#' @export
denoiseq <- function(RDobject, steps, tuningSteps = floor(steps / 3)) {
    # unpacking the counts
    counts <- RDobject@counts
    m <-  nrow(counts)
    # subsetting to obtain matrices for each condition
    counts_A <-  counts[, RDobject@replicates$A]
    n_A <- ncol(counts_A)
    counts_B  <-  counts[, RDobject@replicates$B]
    n_B <- ncol(counts_B)
    # Calculating the size factors used in normalizing the counts
    propotns <- size_factors(counts)
    # initialising a list for the initial values
    parm_samples <- list(RDobject@initValues)
    # initialising a list for the step sizes
    stepsize_vectr <- list(RDobject@stepSizes)
    # tuning
    for (i in 2:tuningSteps) {
        results <- gibbsampling2(counts, counts_A, counts_B, parm_samples[[i - 1]], propotns,
                                 stepsize_vectr[[i - 1]], m, n_A, n_B)
        parm_samples[[i]] <- results$parms
        stepsize_vectr[[i]] <- results$stepsize
    }
    # setting tuned stepsizes
    step.size <- results$stepsize
    for (i in (tuningSteps + 1):steps) {
        parm_samples[[i]] <- gibbsampling(counts, counts_A, counts_B, parm_samples[[i - 1]],
                                          propotns, step.size, m, n_A, n_B)
    }
    RDobject <- setOutput(RDobject, list(samples = parm_samples, stepsize = stepsize_vectr))
    return(RDobject)
}
