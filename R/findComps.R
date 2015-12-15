#' Finds columns in the M3D test-statistic matrix
#'
#' Returns the columns of the test-statistic matrix that refer to specific
#'  samples. This is not intended to be called directly by the user.

#' @param MMD A matrix containing the M3D test-statistic, the difference the
#'  full and methylation blind metrics,
#' for each region in the CpGs object. Each column is a comparison between two
#'  samples, which are described in the column names.
#' @param samples A vector of sample pairs of the form 'sample1 vs sample2' as
#'  returned from determineGroupComps
#' @return Returns the indices of the M3D test-statistic components that contain
#'  the sample pair comparisons in 'samples'
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}

findComps <- function(MMD,samples){
  
  id <- lapply(samples, function(s){
    which(colnames(MMD)==s)
  })
  return(unlist(id))
}

#' Toy data for the package - 1000 CpG regions to be tested in a GRanges object
#'
#' @name CpGsDemo
#' @docType data
#' @author Tom Mayo
#' @keywords data
NULL

#' Toy data for the package - the output of the M3D_Wrapper function.
#'
#' @name MMDlistDemo
#' @docType data
#' @author Tom Mayo
#' @keywords data
NULL

#' Toy data for the package - methylation data for cytosines sites within the testing regions only, in an rrbs object.
#'
#' @name rrbsDemo
#' @docType data
#' @author Tom Mayo
#' @keywords data
NULL

#' Toy data for the package - the output of the pvals function.
#'
#' @name PDemo
#' @docType data
#' @author Tom Mayo
#' @keywords data
NULL
