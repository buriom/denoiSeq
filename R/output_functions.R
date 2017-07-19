# _____________________________________________________________________________

# Samples' extraction _____________________________________________________________________________


#' Get posterior samples  of a parameter
#'
#' Extracts posterior samples of individual parameters contained in the  output slot of
#' the readsData object returned by denoiSeq.
#'
#' @param parm A parameter name string i.e  p, f or gene name.
#' @param RDobject A readsData object with a filled output slot.
#' @param steps An integer representing number of iterations  used while calling denoiseq.
#' @param condition A character (either  A or B) representing the two experimental conditions.
#'
#' @return A vector  of parameter samples, of length equal to steps.
#'
#' @examples
#' #pre-filtering to remove lowly expressed genes
#' ERCC <- ERCC[rowSums(ERCC)>0,]
#' RD <- new('readsData',counts = ERCC)
#' steps <- 30
#' #30 steps are just for illustration here. Atleast 5000 steps are adequate.
#' BI <- denoiseq(RD,steps)
#' samples <- getSamplesOf(BI,'ERCC-00051',steps)
#' plot(samples,type='l', main = 'History plot of ERCC-00051')
#'
#' @export
getSamplesOf <- function(RDobject, parm, steps, condition = "A") {
    N_Asamples <- t(sapply(1:steps, getN_A, RDobject = RDobject))
    N_Bsamples <- t(sapply(1:steps, getN_B, RDobject = RDobject))
    colnames(N_Asamples) <- colnames(N_Bsamples) <- RDobject@geneNames
    if (is.character(parm)) {
        if (parm == "p") {
            p_samples <- t(lapply(1:steps, getp, RDobject = RDobject))
            return(unlist(p_samples))
        } else if (parm == "f") {
            f_samples <- t(lapply(1:steps, getf, RDobject = RDobject))
            return(unlist(f_samples))
        } else if (parm %in% RDobject@geneNames & condition == "A") {
            N <- N_Asamples[, parm]
            return(unlist(N))
        } else if (parm %in% RDobject@geneNames & condition == "B") {
            N <- N_Bsamples[, parm]
            return(unlist(N))
        } else {
            print("Unknown parameter")
        }
    } else if (!is.character(parm)) {
        if (condition == "A") {
            N_Asamples <- t(sapply(1:steps, getN_A, RDobject = RDobject))
            N <- N_Asamples[, parm]
            return(unlist(N))
        } else if (condition == "B") {
            N_Bsamples <- t(sapply(1:steps, getN_B, RDobject = RDobject))
            N <- N_Bsamples[, parm]
            return(unlist(N))
        } else {
            print("Unknown parameter")
        }
    }
}

#' Get  values of the tuned step sizes.
#'
#' Extracts the tuned step sizes for sampling each parameter from  the return value of  denoiseq.
#'
#' @param RDobject A readsData object with a filled output slot.
#' @return A list of the tuned step sizes for sampling each of the parameters.
#' @examples
#' #pre -filtering to remove lowly expressed genes
#' ERCC <- ERCC[rowSums(ERCC)>0,]
#' RD <- new('readsData',counts = ERCC)
#' steps <- 30
#' #30 steps are just for illustration here. Atleast 5000 steps are adequate.
#' BI <- denoiseq(RD,steps)
#' tunedStepSize(BI)
#'
#' @export
#'
tunedStepSize <- function(RDobject) {
    rezult <- RDobject@output
    tune <- length(rezult[[2]])
    return(rezult[[2]][[tune]])
}
