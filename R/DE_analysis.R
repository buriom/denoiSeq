# _____________________________________________________________________________

# Determining differential expression
# _____________________________________________________________________________

# To obtain parameter values sampled at each step of Gibbs sampling.
getN_A <- function(step, RDobject) {
    rezult <- RDobject@output
    return(rezult[[1]][[step]]$N_A)
}
getN_B <- function(step, RDobject) {
    rezult <- RDobject@output
    return(rezult[[1]][[step]]$N_B)
}
getp <- function(step, RDobject) {
    rezult <- RDobject@output
    return(rezult[[1]][[step]]$p)
}
getf <- function(step, RDobject) {
    rezult <- RDobject@output
    return(rezult[[1]][[step]]$f)
}


# Calculate the test statistic.

DEstat <- function(N_A, N_B, rope_limit = 0.5) {
    # log2 difference of samples of the same parameter across the 2 conditions
    dif <- log2(N_A) - log2(N_B)
    # region of practical equivalence
    rope <- sum(dif > -rope_limit & dif < rope_limit)
    return(rope/length(N_A))
}

#' Extract results from denoiseq analysis
#'
#' Extracts the results of Bayesian inference returned by denoiSeq function and
#'  computes the summary and test statistics.
#'
#' @param RDobject A readsData object with the output slot filled.
#' @param steps  An integer representing the number of iterations.
#' @param burnin An integer for the number of iterations to be considered as burn in values
#' @param rope-limit A float that delimits the range of the region of practical
#'  equivalance, ROPE. A default value of 0.5 is used.
#'
#' @return A dataframe with 3 columns namely; the log2 fold change (log2FC), the standard error of the log2 fold change (lgfcSE) and the test static (stat).
#'
#' @examples
#' #pre -filtering to remove lowly expressed genes
#' ERCC <- ERCC[rowSums(ERCC)>10,]
#' RD <- new('readsData',counts = ERCC)
#' steps <- 100
#' #100 steps are not adequate. Just for demonstration here.
#' BI <- denoiseq(RD,steps)
#' rez <- results(BI,steps)
#' head(rez)
#' @export
results <- function(RDobject, steps, burnin = floor(steps/3), rope_limit = 0.5) {
    N_Asamples <- t(sapply(1:steps, getN_A, RDobject = RDobject))
    N_Bsamples <- t(sapply(1:steps, getN_B, RDobject = RDobject))
    N_Asamples <- tail(N_Asamples, steps - burnin)
    N_Bsamples <- tail(N_Bsamples, steps - burnin)
    m <- ncol(N_Asamples)
    values <- mapply(DEstat, split(t(N_Asamples), 1:m), split(t(N_Bsamples), 1:m),
        rope_limit = rope_limit)
    lfc_mat <- N_Bsamples/N_Asamples
    lfc_mean <- apply(lfc_mat, 2, mean)
    lfc_SE <- apply(lfc_mat, 2, sd)/sqrt(nrow(lfc_mat))
    df <- data.frame(`log2FC(B vs A)` = lfc_mean, lfcSE = lfc_SE, stat = values)
    rownames(df) <- RDobject@geneNames
    return(df)
}
