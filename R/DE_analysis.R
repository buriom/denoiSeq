# _____________________________________________________________________________

# Determining differential expression
# _____________________________________________________________________________

# To obtain parameter values sampled at each
# step of Gibbs sampling.
getN_A <- function(rezult, step) return(rezult[[1]][[step]]$N_A)
getN_B <- function(rezult, step) return(rezult[[1]][[step]]$N_B)
getp <- function(rezult, step) return(rezult[[1]][[step]]$p)
getf <- function(rezult, step) return(rezult[[1]][[step]]$f)


#' Get samples of a parameter
#'
#' To extract samples of individual parameters from the list output of nstepgibbs
#' produced by  multiple runs of the gibbs sampler.
#'
#' @param parm A string of the parameter name i.e  p, f or transcript name.
#' @param rezult A list parameter values at each step
#' @param steps An integer representing number of iterations of the gibbs sampler.
#' @param condition An integer representing the two experimental conditions; typically, 1 or 2.
#'
#' @return A vector of parameter samples

get_samples_of <- function(parm, rezult, steps = length(rezult),
    condition = 1) {
    if (parm == "p") {
        p_samples <- t(lapply(1:steps, getp, rezult = rezult))
        return(p_samples)
    } else if (parm == "f") {
        f_samples <- t(lapply(1:steps, getf, rezult = rezult))
        return(f_samples)
    } else if (parm %in% names(rezult[[1]][[1]]$N_A) &
        condition == 1) {
        N_Asamples <- t(mapply(getN_A, 1:steps,
            MoreArgs = list(rezult = rezult)))
        N <- N_Asamples[, parm]
        return(N)
    } else if (parm %in% names(rezult[[1]][[1]]$N_B) &
        condition == 2) {
        N_Bsamples <- t(mapply(getN_B, 1:steps,
            MoreArgs = list(rezult = rezult)))
        N <- N_Bsamples[, parm]
        return(N)
    } else {
        print("Unknown parameter")
    }
}

#' Differential Expression test
#'
#' Given two samples of a parameter N_i representing each of the two conditions,
#'  it determines significance of the differences between the two conditionss based on
#'  a cutoff value.
#'
#' @param N_A A vector of samples of N_i for condition A.
#' @param N_B A vector of samples of N_i for condition B.
#' @param cutoff A threshold value that determines the significance of the difference.
#'
#' @return A logical value about the significance of difference.

DEtest <- function(N_A, N_B, cutoff) {
    # log2 difference of samples of the same
    # parameter across the 2 conditions
    dif <- log2(N_A) - log2(N_B)
    # region of practical equivalence
    rope <- sum(dif > -0.5 & dif < 0.5)
    return((rope/length(N_A)) < cutoff)
}

#' Differential Expression Analysis based on a bottom up model
#'
#' To determine all the differentially expressed genes across the two conditions
#'  in a given dataset.
#'
#' @param burnin An integer for the number of iterations to be considered as burn in values
#' @inheritParams get_samples_of
#' @inheritParams DEtest
#'
#' @return A logival vector highlighting significantly expressed transcripts.

DE_analys <- function(steps, burnin, rezult, cutoff) {
    N_Asamples <- t(mapply(getN_A, 1:steps, MoreArgs = list(rezult = rezult)))
    N_Bsamples <- t(mapply(getN_B, 1:steps, MoreArgs = list(rezult = rezult)))
    N_Asamples <- tail(N_Asamples, steps - burnin)
    N_Bsamples <- tail(N_Bsamples, steps - burnin)
    m <- ncol(N_Asamples)
    values <- mapply(DEtest, split(t(N_Asamples),
        1:m), split(t(N_Bsamples), 1:m), cutoff = cutoff)
    sgf_values <- (colnames(N_Bsamples))[values ==
        TRUE]
    return(sgf_values)
}
