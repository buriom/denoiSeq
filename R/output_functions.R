# _____________________________________________________________________________

# Samples' extraction
# _____________________________________________________________________________


#' Get sampled values  of a parameter
#'
#' To extract samples of individual parameters contianed in the  output slot
#' produced by  denoiSeq.
#'
#' @param parm A parameter name string i.e  p, f or transcript name.
#' @param RD A readsData object with a filled output slot.
#' @param steps An integer representing number of iterations of the one step gibbs sampler.
#' @param condition An character (either  A or B) representing the two experimental conditions.
#'
#' @return A vector of parameter samples

#' @export
getSamplesOf <- function(RD, parm, steps,
    condition = "A") {
    N_Asamples <- t(sapply(1:steps,getN_A,RD=RD))
    N_Bsamples <- t(sapply(1:steps,getN_B,RD=RD))
    colnames(N_Asamples) <- colnames(N_Bsamples) <- RD@geneNames
    if(is.character(parm)){
      if (parm == "p") {
      p_samples <- t(lapply(1:steps, getp, RD = RD))
      return(unlist(p_samples))
    } else if (parm == "f") {
      f_samples <- t(lapply(1:steps, getf, RD = RD))
      return(unlist(f_samples))
    } else if (parm %in% RD@geneNames &
               condition == "A") {
      N <- N_Asamples[, parm]
      return(unlist(N))
    }
     else if (parm %in% RD@geneNames &
               condition == "B") {
      N <- N_Bsamples[, parm]
      return(unlist(N))
     }
      else {
        print("Unknown parameter")
      }
    }
    else if(!is.character(parm)){
      if (condition == 1) {
        N_Asamples <- t(sapply(1:steps,getN_A,RD=RD))
        N <- N_Asamples[, parm]
        return(unlist(N))

      } else if(condition ==2){
        N_Bsamples <- t(sapply(1:steps,getN_B,RD=RD))
        N <- N_Bsamples[, parm]
        return(unlist(N))
      }
      else {
        print("Unknown parameter")
      }
    }

}

#' Get tuned values of the step size.
#'
#' To extract samples of individual parameters contianed in the  output slot
#' produced by  denoiSeq.
#'
#' @param RD A readsData object with a filled output slot.
#' @return A list of the tuned step sizes of each of the parameter.

#' @export
#'
tunedStepSize <- function(RD ){
  rezult <- RD@output
  tune <- length(rezult[[2]])
  return (rezult[[2]][[tune]])
}
