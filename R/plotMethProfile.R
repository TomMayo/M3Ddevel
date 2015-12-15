#' Plots methylation profiles over a specific region
#' 
#' Plots a smoothed methylation profile for each of the two testing groups. 
#' Within each group, the mean of 
#' methylation level is taken, smoothed and plotted, along with the individual 
#' values. 
#' 
#' @param rrbs An rrbs object containing methylation and coverage data as 
#' created using the BiSeq pacakge
#' @param CpGs A GRanges object with each row being a testing region
#' @param group1 The name of the first testing group
#' @param group2 The name of the second testing group
#' @param CpGindex The index within the CpGs object of the region we are 
#' plotting
#' @return NULL, the function plots the profiles
#' @export
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' # plot the 9th region in the Toy Data Set
#' data(rrbsDemo)
#' data(CpGsDemo)
#' plotMethProfile(rrbsDemo, CpGsDemo, 'H1-hESC', 'K562', 9)

plotMethProfile<- function(rrbs,CpGs,group1,group2,CpGindex){
  Ov <- findOverlaps(CpGs[CpGindex],rowRanges(rrbs))
  CpGisland <- rowRanges(rrbs)[subjectHits(Ov)]
  if (length(CpGisland)==0){   # in case we pass an empty island
    warning('This island is empty, returning without plotting')
    return()
  }
  # STEP 1. Get the Cystosine sites
  CytoSites <- start(ranges(CpGisland))
  Ovlaps <- findOverlaps(CpGisland,rowRanges(rrbs))
  sampleNames1 <- rownames(colData(rrbs))[colData(rrbs)$group==group1]   
  # STEP 2. Get the methylation data. 
  #Average over the groups. Find the proportions
  sampleNames2 <- rownames(colData(rrbs))[colData(rrbs)$group==group2]
  methylData1 <- methReads(rrbs)[subjectHits(Ovlaps),sampleNames1]
  methylData2 <- methReads(rrbs)[subjectHits(Ovlaps),sampleNames2]
  coverage1 <- totalReads(rrbs)[subjectHits(Ovlaps),sampleNames1]
  coverage2 <- totalReads(rrbs)[subjectHits(Ovlaps),sampleNames2]
  if (all(coverage1==0) | all(coverage2==0)){
    message('One sample has no coverage')
    return()
  }
  methylGroup1 <- rowSums(methylData1)   # group the data
  methylGroup2 <- rowSums(methylData2)
  covGroup1 <- rowSums(coverage1)
  covGroup2 <- rowSums(coverage2)
  methLvl1 <- methylGroup1/covGroup1 # proportions
  methLvl2 <- methylGroup2/covGroup2
  # STEP 3. Smooth the data
  h = 80 # bandwidth  
  numCs <- length(methylData1[,1])
  numSamples <- length(c(sampleNames1,sampleNames2))
  # Weights: each columns is the weights for that C
  Weights <- unlist(lapply(1:numCs^2, function(i){
    row <- ceiling(i/numCs)
    col <- ((i-1) %% (numCs)) + 1
    if (abs(CytoSites[row]-CytoSites[col]) <= h){
      1 - (abs(CytoSites[row]-CytoSites[col]))/h
    } else {
      0}}))
  Weights <- matrix(Weights, nrow=numCs, ncol=numCs, byrow=TRUE)
  smoothMethylData1 <- t(t(methylData1)%*%Weights)/t(t(coverage1)%*%Weights)
  smoothMethylData2 <- t(t(methylData2)%*%Weights)/t(t(coverage2)%*%Weights)
  smoothGroup1 <- rowMeans(smoothMethylData1,na.rm=TRUE)
  smoothGroup2 <- rowMeans(smoothMethylData2, na.rm=TRUE)
  # STEP 4 - plot, order the cytosines first
  ord <- order(CytoSites) 
  CytoSites <- CytoSites[ord]
  methLvl1 <-methLvl1[ord]
  methLvl2 <- methLvl2[ord]
  smoothGroup1 <- smoothGroup1[ord]
  smoothGroup2 <- smoothGroup2[ord]
  
  plot(CytoSites,methLvl1,col='blue',
       ylim=range(c(methLvl1,methLvl2),na.rm=TRUE),ann=FALSE)
  lines(CytoSites,smoothGroup1,col='blue',lty=1, lwd=2)
  points(CytoSites,methLvl2,col='red',pch=22)
  lines(CytoSites,smoothGroup2,col='red',lty=1,lwd=2)
  title(ylab='Methylation Level')
  xlabel = paste(unique(seqnames(CpGisland)),unique(strand(CpGisland)),sep=', ')
  title(xlab=xlabel)
  title(main=paste("Island",CpGindex,sep=' '), font.main=4)
}
