#' Finds the median
#'
#' Returns the median of a list of values with corresponding frequencies. This
#'  is not intended to be called directly by the user.
#' 
#' @param values A vector of the unique values that occur
#' @param freqs A vector of the number of occurrences of each value
#' @return Returns the median value of the data comprising each entry in values
#'  repeated the corresponding entry in freqs number of times, as a numeric.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
medianFreq <- function(values,freqs){
  if (all(freqs==0)){
    med<-0
    return(med)
  }
  nonzeros<-which(freqs!=0)   #take out the zero frequencies
  values<-values[nonzeros]
  freqs<-freqs[nonzeros]
  ord<-order(values)  #sort
  values<-values[ord]
  freqs<-freqs[ord]
  # find where the median must lie, different for odd and even counts
  if (sum(freqs) %% 2 == 1) {
    medInd <- ceiling(sum(freqs)/2) # find halfway index
    cumfreqs<-cumsum(freqs) # find where this hits
    med<-unique(values[which(cumfreqs==min(cumfreqs[cumfreqs>=medInd]))])
  } else {
    medInd1 <- sum(freqs)/2 # find first halfway index
    medInd2 <- medInd1 + 1 # find second
    cumfreqs<-cumsum(freqs) # find where this hits
    med1 <-unique(values[which(cumfreqs==min(cumfreqs[cumfreqs>=medInd1]))])
    med2 <-unique(values[which(cumfreqs==min(cumfreqs[cumfreqs>=medInd2]))])
    med <- (med1 + med2) / 2
  }
  return(med)
}
