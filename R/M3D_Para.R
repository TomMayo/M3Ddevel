#' Computes the components of the M3D test-statistic over all regions for all
#'  sample-pairs.
#' 
#' Parallel implementation of M3D_Wrapper function. 
#' Returns the two components of the M3D test-statistic - the MMD
#'  (Gretton et al. 2006) for the full data and the coverge
#' only data, respectively - for all regions and all samples pairs, as a matrix.
#' 
#' @param rrbs An rrbs object containing methylation and coverage data as
#'  created using the BiSeq pacakge
#' @param CpGs  A GRanges object detailing the testing regions.
#' @param overlaps The overlaps between the list of testing regions and the
#'  methylation data. This is obtained using the
#' function findOverlaps(CpGs,rrbs) for a GRanges object CpGs detailing the
#'  testing regions.
#' @param num.cores Integer giving the number or cores to use. Defaults to the 
#'  maximum available
#' @return This returns the two components of the M3D test-statistic for each
#'  region over all sample pairs as a matrix. 
#' Subtracting them gives the M3D test-statistic. This is processed with the
#'  function pvals.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @references Gretton, A., Borgwardt, K. M., Rasch, M., Scholkopf, B., Smola,
#'  A. J. (2006). A kernel method for the two-sample-problem. In Advances in
#'   neural information processing systems (pp. 513-520).
#' @examples
#' \donttest{data(rrbsDemo)
#' data(CpGsDemo)
#' M3D_list <- M3D_Para(rrbsDemo,CpGsDemo)
#' head(M3d_list$Full-M3D_list$Coverage)}
#' @export


M3D_Para <- function(rrbs, CpGs, overlaps,num.cores=NaN){
    # use max cores if num cores not specified
    if (is.na(num.cores)){
        num.cores <- detectCores()
    }
    
    numIslands <- length(CpGs)
    islsPerCore <- ceiling(numIslands/num.cores) # number of islands per core
    
    rrbsList <- list()
    ovList <- list()
    indList <- list()
    for (i in 1:num.cores){
        if (i==num.cores){
            inds <- (((i-1)*islsPerCore)+1):numIslands
        } else {
            inds <- (((i-1)*islsPerCore)+1):(i*islsPerCore)
        }
        
        ovs <- overlaps[queryHits(overlaps) %in% inds]
        rrbsList[[i]]<-rrbs[subjectHits(ovs)]
        len <- length(subjectHits(ovs))
        q.map <- as.integer(queryHits(ovs)-(queryHits(ovs)[1]-1))
        s.map <- as.integer(1:len)
        ovs <- new('Hits',queryHits=q.map,subjectHits=s.map,
                   queryLength=len,subjectLength=len)
        ovList[[i]] <- ovs
        indList[[i]] <- inds
    }
    paraList <- rep(TRUE,num.cores)
    
    # run M3D in parallel
    MMDret <- mcmapply(function(x, y, z) M3D_Wrapper(x,y,z), rrbsList, ovList, 
                       paraList, mc.cores=num.cores, mc.preschedule=FALSE)
    
    # reconstruct the MMD matrices to return - col names first
    nSamples = length(colnames(methReads(rrbs)))
    if (nSamples==2){
        samplesIdx <- c(1,2)
        numPairs <- 1
    } else {
        a <- unlist(lapply(1:(nSamples-1), function(i) rep(i,(nSamples-i))))
        b <- unlist(lapply(1:(nSamples-1), function(i) (i+1):nSamples))   
        samplesIdx <- cbind(a,b)
        numPairs <- length(samplesIdx[,1])
    }
    ColumnNames <- vector()
    for (pairInd in 1:numPairs){
        if (numPairs==1){
            pair <- c(1,2)
        } else {
            pair <- samplesIdx[pairInd,]
        }
        sample1 <- colnames(methReads(rrbs))[pair[1]]
        sample2 <- colnames(methReads(rrbs))[pair[2]]
        ColumnNames <- c(ColumnNames,paste(sample1, ' vs ', sample2))
    }    
    
    #now matrix
    fullmx <- matrix(NaN,ncol=length(ColumnNames)*2, nrow=numIslands)
    if (is.matrix(MMDret)){
        for (i in 1:num.cores){
            temp <- matrix(MMDret[,i], ncol=length(ColumnNames)*2, nrow=length(indList[[i]]))
            fullmx[indList[[i]],] <- temp
        }
    } else if (is.list(MMDret)) {
        startRow <- 1
        for (i in 1:num.cores){
            numCols <- length(MMDret[[i]][,1])
            stopRow <- startRow+numCols-1
            fullmx[startRow:stopRow,] <- MMDret[[i]]
            startRow <- stopRow+1
        }
    } else {
        stop('The parallel call returns neither a list nor a matrix')
    }
    end <- length(fullmx[1,])/2
    MMD <- fullmx[,1:end]
    MMDCoverage <- fullmx[,(end+1):(2*end)]
    colnames(MMD) <- ColumnNames
    colnames(MMDCoverage) <- ColumnNames
    ret <- list(MMD,MMDCoverage)
    names(ret) <- c('Full','Coverage')    
    return(ret)
}
