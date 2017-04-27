######################################################################
# Create the base readsData class
#
# This is used to represent the most basic agent in a simulation.
readsData <- setClass(
  # Set the name for the class
  "readsData",

  # Define the slots
  slots = c(
    counts = "matrix",
    annotation = "character",
    replicates   = "list",
    init.values   = "list",
    step.sizes   = "list"
  ),

  # Set the default values for the slots. (optional)
  # prototype=list(
  #   init.values = list(N_A = rep(1, nrow(counts)),
  #                      N_B = rep(1,nrow(counts)), p = 0.0001, f = 0.01),
  #   step.sizes = list(stepsizeN_A = rep(1, nrow(counts)),
  #                     stepsizeN_B = rep(1,nrow(counts)),
  #                     stepsize_p = 5e+07, stepsize_f = 1e3)
  # ),

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  {
    if(sum(object@counts)<0) {
      return("counts cannot be negative")
    }
    return(TRUE)
  }
)

setMethod ("initialize", signature  = "readsData",
           definition = function (.Object,
                                  replicates,
                                  counts,
                                  annotation,
                                  step.sizes = list(stepsizeN_A = rep(1, nrow(counts)),
                                                    stepsizeN_B = rep(1,nrow(counts)),
                                                    stepsize_p = 5e+07, stepsize_f = 1e3),
                                  init.values = list(N_A = rep(1, nrow(counts)),
            N_B = rep(1,nrow(counts)), p = 0.0001, f = 0.01)){
             .Object@counts <- counts
             .Object@replicates <- replicates
             .Object@step.sizes <- step.sizes
             .Object@annotation <- annotation
             .Object@init.values <- init.values
             return (.Object)
           })


replicates <- list(A=1:5,B=6:10)
annotations <- row.names(ERCC)
init.values = list(N_A = rep(1, 71),
N_B = rep(1,71), p = 0.0001, f = 0.01)
CD <- new("readsData", counts = ERCC,replicates=replicates,annotation = annotations)

denoiSeq <- function(CD,steps, tuning.steps = floor(steps/3)) {
  #unpacking the counts
  counts <- CD@counts
  m = nrow(counts)
  #subsetting to obtain matrices for each condition
  counts_A = counts[, CD@replicates$A]
  n_A <- ncol(counts_A)
  counts_B = counts[,CD@replicates$B]
  n_B <- ncol(counts_B)

  # Calculating the size factors used in
  # normalizing the counts
  propotns <- size_factors(counts)
  if(is.null(rownames(counts))){
    names(init.valuess$N_A) <- names(init.valuess$N_B) <- CD@annotation
  }
  #initialising a list for the initial values
  parm_samples <- list(CD@init.values)
  #initialising a list for the step sizes
  stepsize_vectr <- list(CD@step.sizes)
  #tuning
  for (i in 2:tuning.steps) {
    results <- gibbsampling2(counts, counts_A,
                             counts_B, parm_samples[[i - 1]], propotns,
                             stepsize_vectr[[i - 1]], m, n_A, n_B)
    parm_samples[[i]] <- results$parms
    stepsize_vectr[[i]] <- results$stepsize
  }
  #setting tuned stepsizes
  step.size <- results$stepsize
  for (i in (tuning.steps + 1):steps) {
    parm_samples[[i]] <- gibbsampling(counts,
                                      counts_A, counts_B, parm_samples[[i -
                                                                          1]], propotns, step.size, m, n_A,
                                      n_B)
  }
  return(list(samples = parm_samples, stepsize = stepsize_vectr))
}

