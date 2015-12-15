#' Reads in ENCODE RRBS data
#'
#' Reads in RRBS data in bed file format from the ENCODE consortium and outputs
#'  an rrbs data structure. Adapted from readBismark in the BiSeq package.
#' 
#' @param files A character pointing the the rrbs files downloads from the
#'  ENCODE database.
#' @param colData Samples' names plus additional sample information as
#'  character, data.frame or DataFrame.
#' @param eData Experiment data to describe the work. This is used to create 
#' the BSraw object as in the BiSeq package.
#' @return Returns a BSraw object storing methylation and coverage data -
#'  the underlying structure for this package.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @export
#' @examples
#' \donttest{
#' # download the files and change the working directory
#' # to that location
#' files <- c('wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bed.gz',
#' 'wgEncodeHaibMethylRrbsH1hescHaibSitesRep2.bed.gz',
#' 'wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed.gz',
#' 'wgEncodeHaibMethylRrbsK562HaibSitesRep2.bed.gz')  
#' group <- factor(c('H1-hESC','H1-hESC','K562','K562'))
#' samples <- c('H1-hESC1','H1-hESC2','K562-1','K562-2')
#' colData <- DataFrame(group,row.names= samples)
#' rrbs <- readENCODEdata(files,colData)}

readENCODEdata <- function(files, colData,eData=NaN) {
  if (nrow(colData) != length(files)) {
    stop("Row number of colData must equal length of files.")
  }
  
  methData = list()
  
  for (i in 1:length(files)) {
    cat(paste("Processing sample ", rownames(colData)[i], " ... \n", sep=""))
    
    bismark <- scan(files[i], skip=1, sep="\t",
                    what=list("character", integer(), NULL, NULL, 
                              integer(), 'character',
                              NULL,NULL,NULL,NULL,integer()))  
    message('Begin sorting data')
    bismarkEntries <- c(1,2,5,6,11)
    Order <- with(bismark,order(bismark[[6]],bismark[[1]],bismark[[2]]))
    for (j in 1:length(bismarkEntries)){   
      entry <- bismarkEntries[j]
      bismark[[entry]] <- bismark[[entry]][Order]
    }
    rm(Order)
    message('Data sort complete')
    methData[[i]] = GRanges(
      seqnames=bismark[[1]],
      strand=bismark[[6]],
      ranges=IRanges(start=bismark[[2]], width=1),
      methylated=as.integer(round(0.01*bismark[[5]]*bismark[[11]])),
      reads=bismark[[5]])
    
    rm(bismark)
  }
  
  cat("Building BSraw object.\n")
  
  fData <- methData[[1]]
  
  if(length(methData) > 1){
    for(i in seq(along=methData)[-1]){
      fData <- unique(c(fData, methData[[i]]))
    }
  }
  elementMetadata(fData) <- NULL
  names(fData) <- as.character(1:length(fData))
  
  tReads <- matrix(integer(length = length(fData) * 
                             length(methData)), nrow=length(fData))
  mReads <- matrix(integer(length = length(fData) * 
                             length(methData)), nrow=length(fData))
  
  for(i in seq(along=methData)){
    mtch <- findOverlaps(fData, methData[[i]])
    mtch.m <-  as.matrix(mtch)
    ind1 <- mtch.m[, 1]
    ind2 <- mtch.m[, 2]
    tReads[ind1, i] <- elementMetadata(methData[[i]])$reads[ind2]
    mReads[ind1, i] <- elementMetadata(methData[[i]])$methylated[ind2]
  }
  
  colnames(tReads) <- rownames(colData)
  colnames(mReads) <- rownames(colData)
  rownames(tReads) <- names(fData)
  rownames(mReads) <- names(fData)
  if(is.na(eData)){
    eData <- SimpleList(exp='Demo')
  }
  
  rrbs = BSraw(
    exptData=eData,
    colData = colData,
    rowRanges = fData,
    totalReads = tReads,
    methReads = mReads)
  
  return(rrbs)
}
