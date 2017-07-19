# _____________________________________________________________________________

# Determining differential expression
#  _____________________________________________________________________________

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
    return(rope / length(N_A))
}

#' Compute the test statistic
#'
#' Extracts  posterior samples of the  parameters which are returned by denoiSeq function and
#'  computes the summary and test statistics.
#'
#' To calculate the  test statistic, this function  first log2 transforms the  posterior samples of the two relevant parameters
#' i.e N_iA and N_iB. It then subtracts posterior samples of one of the parameters  from the other and
#' determines the proportion of this
#' distribution of differences that lies in the region of practical equivalence (ROPE) (Kruschke, 2011).
#' The genes can then be arranged in an ascending order of  the ROPE_propn column and we can select the most
#' differentially expressed genes as those whose ROPE_propn is less than a particular threshold value.
#'
#' Using both real and  simulated data, optimal values between 0.0007 and 0.4 were obtained for the threshold.
#'
#' @param RDobject A readsData object with a filled output slot.
#' @param steps  An integer representing the number of iterations.
#' @param burnin An integer for the number of iterations to be considered as burn in values. A default value
#' equivalent to a third of steps is used.
#' @param rope_limit A float that delimits the range of the region of practical
#'  equivalence, ROPE. A default value of 0.5 is used.
#'
#' @return A dataframe with 3 columns namely; the log2 fold change (log2FC), the standard error of the log2
#' fold change (lgfcSE) and the test static ( ROPE_propn).
#'
#' @examples
#' #pre -filtering to remove lowly expressed genes
#' ERCC <- ERCC[rowSums(ERCC)>0,]
#' RD <- new('readsData',counts = ERCC)
#' steps <- 30
#' #30 steps are just for illustration here. At least 5000 steps are adequate.
#' BI <- denoiseq(RD,steps)
#' rez <- results(BI,steps)
#' head(rez)
#'
#' #Re-ordering according to most differentially expressed
#' rez <- rez[with(rez,order( ROPE_propn)),]
#' head(rez,10)
#'
#' #Determine significant genes using a threshold of 0.38.
#' sgf <- rez[rez$ ROPE_propn<0.38,]
#' head(sgf)
#' dim(sgf)
#' @importFrom utils tail
#' @export

results <- function(RDobject, steps, burnin = floor(steps / 3), rope_limit = 0.5) {
    N_Asamples <- t(sapply(1:steps, getN_A, RDobject = RDobject))
    N_Bsamples <- t(sapply(1:steps, getN_B, RDobject = RDobject))
    N_Asamples <- tail(N_Asamples, steps - burnin)
    N_Bsamples <- tail(N_Bsamples, steps - burnin)
    m <- ncol(N_Asamples)
    values <- mapply(DEstat, split(t(N_Asamples), 1:m), split(t(N_Bsamples), 1:m), rope_limit = rope_limit)
    lfc_mat <- N_Bsamples / N_Asamples
    lfc_mean <- apply(lfc_mat, 2, mean)
    lfc_SE <- apply(lfc_mat, 2, sd) / sqrt(nrow(lfc_mat))
    df <- data.frame(lfc_mean, lfc_SE, values)
    colnames(df) <- c(" log2FC(B / A)", "lfcSE  ", "ROPE_propn")
    rownames(df) <- RDobject@geneNames
    return(df)
}
