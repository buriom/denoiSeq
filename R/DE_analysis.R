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
#' The test statistic is calculated  by first log2 transforming the samples of the relevant parameters
#' (i.e N_iA and N_iB), subtracting one sample from the other, and then determining the propotion of this
#' distribution of differences that lies in the region of practical equivalence (ROPE). We then select
#'  genes whose test statistic value  is below a given threshold as differentially expressed.
#'  From simulated data, a typical threshold value is around 0.3. However, one
#'  can also arrange the stat column in ascending order and select the most differentially expressed genes.
#'
#' @param RDobject A readsData object with a filled output slot.
#' @param steps  An integer representing the number of iterations.
#' @param burnin An integer for the number of iterations to be considered as burn in values.
#' @param rope_limit A float that delimits the range of the region of practical
#'  equivalance, ROPE. A default value of 0.5 is used.
#'
#' @return A dataframe with 3 columns namely; the log2 fold change (log2FC), the standard error of the log2 fold change (lgfcSE) and the test static (stat).
#'
#' @examples
#' #pre -filtering to remove lowly expressed genes
#' ERCC <- ERCC[rowSums(ERCC)>15,]
#' RD <- new('readsData',counts = ERCC)
#' steps <- 100
#' #100 steps are not adequate. Just for illustration here.
#' BI <- denoiseq(RD,steps)
#' rez <- results(BI,steps)
#' head(rez)
#' dim(rez)
#' #Determine significant genes
#' sgf <- rez[rez$stat<0.3,]
#' head(sgf)
#' dim(sgf)
#' #Re-ordering according to most differentially expressed
#' rez <- rez[with(rez,order(stat)),]
#' head(rez,10)
#' @importFrom utils tail
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
    df <- data.frame(`log2FC(B vs A)` = lfc_mean, lfcSE = lfc_SE, ROPE_propn = values)
    rownames(df) <- RDobject@geneNames
    return(df)
}
