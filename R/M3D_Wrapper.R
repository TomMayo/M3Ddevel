#' Computes the components of the M3D test-statistic over all regions for all
#'  sample-pairs.
#'
#' Returns the two components of the M3D test-statistic - the MMD
#'  (Gretton et al. 2006) for the full data and the coverge
#' only data, respectively - for all regions and all samples pairs, as a matrix.
#'
#' @param rrbs An rrbs object containing methylation and coverage data as
#'  created using the BiSeq pacakge
#' @param overlaps The overlaps between the list of testing regions and the
#'  methylation data. This is obtained using the
#' function findOverlaps(CpGs,rrbs) for a GRanges object CpGs detailing the
#'  testing regions.
#' @param para Set to true if called via M3D_Para
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
#' overlaps <- findOverlaps(CpGsDemo,rrbsDemo)
#' M3D_list <- M3D_Wrapper(rrbsDemo,overlaps)
#' head(M3d_list$Full-M3D_list$Coverage)}
#' @export

M3D_Wrapper <- function(rrbs, overlaps, para=FALSE){
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


    # store objects here to avoid slow lookups, potential memory hazard
    m_reads <- methReads(rrbs)
    t_reads <- totalReads(rrbs)
    locs_all <- start(ranges(rowRanges(rrbs)))
    q_hits <- queryHits(overlaps)
    sub_hits <- subjectHits((overlaps))

    islands <- unique(queryHits(overlaps))
    CSites <- rowRanges(rrbs)
    MMD <- matrix(NA,nrow=length(islands),ncol=numPairs)
    MMDCoverage <- matrix(NA,nrow=length(islands),ncol=numPairs)
    col_sample_names <- colnames(m_reads)

    ### make colnames - vectorize
    ColumnNames <- unlist(lapply(1:numPairs, function(pairInd){
        if (numPairs==1){
            pair <- c(1,2)
        } else {
            pair <- samplesIdx[pairInd,]
        }
        sample1 <- colnames(m_reads)[pair[1]]
        sample2 <- colnames(m_reads)[pair[2]]
        return(paste(sample1, ' vs ', sample2))
    }))

    # loop over islands, then over samples
    pb <- txtProgressBar(min=1,max=length(islands),style=3)
    for (i in 1:length(islands)) {
        island <- islands[i]
        methIndices <- sub_hits[q_hits==island]
        meth_isl <- m_reads[methIndices,]
        total_isl <- t_reads[methIndices,]
        unmeth_isl <- total_isl - meth_isl

        # compute the location matrix
        ## define easy labels
        locs <- locs_all[methIndices]
        locs <- locs - min(locs)
        # make location matrix
        G <- locs%*%t(locs)
        L <- ncol(G)
        nor <- rep(G[seq(1,L^2,L+1)],L)
        # locMx is the squared distance between sites
        locMx <- -2*G + matrix(nor, nrow=L, byrow=TRUE) + matrix(nor, nrow=L)
        locInds <- which(locMx!=0)

        for (pairInd in 1:numPairs){
            if (numPairs==1){
                pair <- c(1,2)
            } else {
                pair <- samplesIdx[pairInd,]
            }
            sample1 <- col_sample_names[pair[1]]
            sample2 <- col_sample_names[pair[2]]

            # this is the meth data
            methData <- meth_isl[,c(sample1,sample2)]
            # this is the total data
            unmethData <- unmeth_isl[,c(sample1,sample2)]
            testData <- matrix(0, nrow = dim(methData)[1], ncol = 4)
            testData[,1:2] <- methData
            testData[,3:4] <- unmethData

            # compute the MMDs
            res <- M3D_Single(testData,locMx,locInds,method='MinusCovMMD')
            MMD[i,pairInd] <- res[[1]]
            MMDCoverage[i,pairInd] <- res[[2]]
        }
        setTxtProgressBar(pb,i)
    }
    close(pb)
    colnames(MMD) <- ColumnNames
    colnames(MMDCoverage) <- ColumnNames
    ret <- list(MMD,MMDCoverage)
    names(ret) <- c('Full','Coverage')
    if (para==FALSE){
        return(ret)
    } else {
        return(cbind(MMD,MMDCoverage))
    }
}
