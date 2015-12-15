#' Creates strings of sample pair comparisons
#' 
#' Takes in a vector of strings of sample names and returns strings of 
#' all the comparisons, either within a testing group or between testing groups.
#'  This is not intended to be called directly by the user.
#' 
#' @param samples1 A vector of sample names from one group
#' @param samples2 A vector of sample names from the other group, if we want to
#'  specify between-group comparisons
#' @param type 'within' or 'between'. 'within' returns all the sample pairs 
#' within samples1, 'between' returns all the sample pairs
#' between samples1 and samples2
#' @return A vector of sample pair comparisons of the form 'sample1 vs sample2'
#'  for use with the M3D functions
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}

determineGroupComps <- function(samples1,samples2=NULL,type){
  len <- length(samples1)
  comps <- vector()
  if (type=='within'){
    if ((len)==1){
      stop('There is only one sample in this group, you need replicates')
    }
    comps <- unlist(lapply(1:(len-1), function(i){
      remainder <- samples1[(i+1):len]
      paste(samples1[i], ' vs ', remainder)
    }))
    return(comps)
  }
  else if (type=='between'){  
    comps <- unlist(lapply(1:len, function(i){
      paste(samples1[i], ' vs ', samples2)
    }))
  }
  
  else {
    stop('Type is not valid, choose within or between')
  }
}
