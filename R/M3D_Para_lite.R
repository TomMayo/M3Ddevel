#' Computes the components of the M3D test-statistic over all regions for all
#'  sample-pairs.
#'
#' Parallel implementation of M3D_Wrapper_lite function.
#' Returns the M3D test-statistic, without the compenents, for all regions
#' with the cross group sample pairs averaged to save memory (in column one),
#' as a matrix.
#'
#' @param rrbs An rrbs object containing methylation and coverage data as
#'  created using the BiSeq pacakge
#' @param CpGs  A GRanges object detailing the testing regions.
#' @param overlaps The overlaps between the list of testing regions and the
#'  methylation data. This is obtained using the
#' function findOverlaps(CpGs,rrbs) for a GRanges object CpGs detailing the
#'  testing regions.
#' @param group1 The name of the first group for the comparison. This is stored
#'  in colData(rrbs). Default finds first unique group in colData(rrbs).
#' @param group2 The name of the second group for the comparison. This is stored
#'  in colData(rrbs). Default finds second unique group in colData(rrbs).
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


M3D_Para_lite <- function(rrbs, CpGs, overlaps, group1=NaN,group2=NaN, num.cores=NaN){
    # use max cores if num cores not specified
    if (is.na(group1)){
        group1 <- as.character(unique(colData(rrbs)$group)[1])
    }
    if (is.na(group2)){
        group2 <- as.character(unique(colData(rrbs)$group)[2])
    }
    if (is.na(num.cores)){
        num.cores <- parallel::detectCores()
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
    group1List <- rep(group1,num.cores)
    group2List <- rep(group2,num.cores)

    # run M3D in parallel
    MMDret <- parallel::mcmapply(function(x, y, z, w) M3D_Wrapper_lite(x,y,z,w),
                                 rrbsList, ovList, group1List, group2List,
                                 mc.cores = num.cores, mc.preschedule = FALSE)


    # reconstruct the MMD matrices to return - col names first
    nSamples = sum(colData(rrbs)$group==group1,colData(rrbs)$group==group2)
    if (nSamples==2){
        samplesIdx <- c(1,2)
        numPairs <- 1
    } else {
        a <- unlist(lapply(1:(nSamples-1), function(i) rep(i,(nSamples-i))))
        b <- unlist(lapply(1:(nSamples-1), function(i) (i+1):nSamples))
        samplesIdx <- cbind(a,b)
        numPairs <- length(samplesIdx[,1])
    }
    ### make colnames - vectorize
    ColumnNames <- unlist(lapply(1:numPairs, function(pairInd){
        if (numPairs==1){
            pair <- c(1,2)
        } else {
            pair <- samplesIdx[pairInd,]
        }
        sample1 <- colnames(methReads(rrbs))[pair[1]]
        sample2 <- colnames(methReads(rrbs))[pair[2]]
        return(paste(sample1, ' vs ', sample2))
    }))

    l <- length(ColumnNames)
    temp <- matrix(rep(0,2*l),nrow=2,ncol=l) # to work with findcomps (hack)
    colnames(temp) <- ColumnNames
    ### code from pvals
    samples1 <- rownames(colData(rrbs))[colData(rrbs)[,]==group1]
    samples2 <- rownames(colData(rrbs))[colData(rrbs)[,]==group2]
    within1 <- determineGroupComps(samples1,type='within')
    within2 <- determineGroupComps(samples2,type='within')
    within <- c(within1,within2)
    idsWithin <- c(findComps(temp,within1),findComps(temp,within2))
    ColumnNames <- c('MeanBetween',ColumnNames[idsWithin])


    #now matrix
    M3D_stat <- matrix(NaN,ncol=length(ColumnNames), nrow=numIslands)
    if (is.matrix(MMDret)){
        for (i in 1:num.cores){
            temp <- matrix(MMDret[,i], ncol=length(ColumnNames), nrow=length(indList[[i]]))
            M3D_stat[indList[[i]],] <- temp
        }
    } else if (is.list(MMDret)) {
        startRow <- 1
        for (i in 1:num.cores){
            numCols <- length(MMDret[[i]][,1])
            stopRow <- startRow+numCols-1
            startRow <- stopRow+1
        }
    } else {
        stop('The parallel call returns neither a list nor a matrix')
    }
    colnames(M3D_stat) <- ColumnNames
    return(M3D_stat)
}
