#_______________________________________________________________________________

# Create the base readsData class
#_______________________________________________________________________________

#To define the new class
readsData <- setClass(
  # Set the name for the class
  "readsData",

  # Define the slots
  slots = c(
    counts = "matrix",
    annotation = "character",
    replicates   = "list",
    init.values   = "list",
    step.sizes   = "list",
    output = "list"
  ),
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

#Initializing
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

#method for setting the ouput slot
setGeneric(name="setoutput",
           def=function(theObject,outputval)
           {
             standardGeneric("setoutput")
           }
)

setMethod(f="setoutput",
          signature="readsData",
          definition=function(theObject,outputval)
          {
            theObject@output <- outputval
            return(theObject)
          }
)

#method for getting the contents of the output slot
setGeneric(name="getoutput",
           def=function(theObject)
           {
             standardGeneric("getoutput")
           }
)

setMethod(f="getoutput",
          signature="readsData",
          definition=function(theObject)
          {
            return(theObject@output)
          }
)
#__________________________________________________________

#annotation
#__________________________________________________________
#method for setting the annotation slot
setGeneric(name="setannotation",
           def=function(theObject,annotval)
           {
             standardGeneric("setannotation")
           }
)

setMethod(f="setannotation",
          signature="readsData",
          definition=function(theObject,annotval)
          {
            theObject@annotation <- annotval
            return(theObject)
          }
)

#method for getting the contents of the output slot
setGeneric(name="getannotation",
           def=function(theObject)
           {
             standardGeneric("getannotation")
           }
)

setMethod(f="getannotation",
          signature="readsData",
          definition=function(theObject)
          {
            return(theObject@annotation)
          }
)
#______________________________________________________________

#counts
#_______________________________________________________________
#method for setting the ouput slot
setGeneric(name="setcounts",
           def=function(theObject,countsval)
           {
             standardGeneric("setcounts")
           }
)

setMethod(f="setcounts",
          signature="readsData",
          definition=function(theObject,countsval)
          {
            theObject@counts <- countsval
            return(theObject)
          }
)

#method for getting the contents of the output slot
setGeneric(name="getcounts",
           def=function(theObject)
           {
             standardGeneric("getcounts")
           }
)

setMethod(f="getcounts",
          signature="readsData",
          definition=function(theObject)
          {
            return(theObject@counts)
          }
)


replicates <- list(A=1:5,B=6:10)
annotation <- row.names(ERCC)
init.values = list(N_A = rep(1, 71),
N_B = rep(1,71), p = 0.0001, f = 0.01)
CD <- new("readsData", counts = ERCC,replicates=replicates,annotation = annotation)

