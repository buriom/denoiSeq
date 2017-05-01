#_______________________________________________________________________________

# Create the base readsData class
#_______________________________________________________________________________

#' An S4 class to represent summarised counts and the output of Bayesian Inference.
#'
#' @slot counts  A positive integer matrix containing summarised counts for each
#'  genomic event (genes, exons, transcripts, etc)  in the two conditions.
#' @slot geneNames A character vector containing the names of the genomic event.
#' @slot initValues A list containing initial values for each parameter. Defaulted to \emph{N_A} = rep(1, nrow(counts)),\emph{N_B} = rep(1, nrow(counts)),\emph{p}= 0.0001, \emph{f} = 0.01
#' @slot stepSizes A list containing step sizes for sampling each parameter. Defaulted to stepsizeN_A = rep(1, nrow(counts)),stepsizeN_B = rep(1, nrow(counts)),stepsize_p= 1e3, stepsizeN_A = 5e7
#'
#' @name readsData
readsData <- setClass(
  # Set the name for the class
  "readsData",

  # Define the slots
  slots = c(
    counts = "matrix",
    geneNames = "character",
    replicates   = "list",
    initValues   = "list",
    stepSizes   = "list",
    output = "list"
  )
)

#Initializing
setMethod ("initialize", signature  = "readsData",
           definition = function (.Object,

                                  counts,
                                  geneNames = rownames(counts),
                                  replicates = list(A = 1:(ncol(counts)/2),B= (ncol(counts)/2+1):ncol(counts)),
                                  stepSizes = list(stepsizeN_A = rep(1, nrow(counts)),
                                                    stepsizeN_B = rep(1,nrow(counts)),
                                                    stepsize_p = 5e+07, stepsize_f = 1e3),
                                  initValues = list(N_A = rep(1, nrow(counts)),
            N_B = rep(1,nrow(counts)), p = 0.0001, f = 0.01)){
             .Object@counts <- counts
             .Object@replicates <- replicates
             .Object@stepSizes <- stepSizes
             .Object@geneNames <- geneNames
             .Object@initValues <- initValues
             validObject(.Object)
             return (.Object)
           })

validity=function(object)
{
  rval <- NULL
  if( nrow(object@counts) != length(object@geneNames) ){
    return("reads and geneNamess don't match")
  }
  else if(sum(object@counts < 0)>0){
    return("counts cannot be negative")
  }
  else return(TRUE)
}

setValidity( "readsData",validity)
#__________________________________________________________

#initValues
#__________________________________________________________

#' Mutator method for the initValues slot of the readsData object.
#'
#' Alters the value of the initValues slot.
#'
#' @param theObject a readsData object
#' @param initval A list of initial values for each parameter.
#' @return The same readsData object with the initValues slot updated.
#'
setGeneric(name="setInitValues",
           def=function(theObject,initval)
           {
             standardGeneric("setInitValues")
           }
)

#method for setting the ouput slot
setMethod(f="setInitValues",
          signature="readsData",
          definition=function(theObject,initval)
          {
            theObject@initValues <- initval
            return(theObject)
          }
)
#______________________________________________________________

#stepSizes
#_______________________________________________________________
#' Mutator method for the stepSizes slot of the readsData object.
#'
#' Alters the value of the stepSizes slot.
#'
#' @param theObject a readsData object.
#' @param stepSizesval A list of step sizes for each parameter.
#' @return The same readsData object with the stepSizes slot updated.
#'
setGeneric(name="setStepSizes",
           def=function(theObject,stepSizesval)
           {
             standardGeneric("setStepSizes")
           }
)

#method for setting the ouput slot
setMethod(f="setStepSizes",
          signature="readsData",
          definition=function(theObject,stepSizesval)
          {
            theObject@stepSizes <- stepSizesval
            return(theObject)
          }
)

# replicates <- list(A=1:5,B=6:10)
# load("/home/buri/denoiSeq/data/ERCC.RData")
# geneNames <- row.names(ERCC)
# initValues = list(N_A = rep(1, 71),
# N_B = rep(1,71), p = 0.0001, f = 0.01)
# CD <- new("readsData", counts = ERCC,geneNames = geneNames,initValues = initValues)
