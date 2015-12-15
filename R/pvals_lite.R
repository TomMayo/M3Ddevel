#' Computes p-values
#'
#' Returns p-values for each region reflecting the probability of observing the
#'  mean test-statistic
#' of the between group comparisons among the inter-replicate comparisons. The
#' 'lite' version works with the outcome of M3D_Wrapper_lite and the 'lite'
#' parallel version.
#'
#' @param rrbs An rrbs object containing methylation and coverage data as
#'  created using the BiSeq pacakge
#' @param CpGs A GRanges object with each row being a testing region
#' @param M3D_stat_lite A matrix containing the M3D test-statistic, as returned
#' from the M3D_Wrapper_lite function. This has the mean of the MMDs for the
#' comparisons between groups in the first column, and the within groups MMDs
#' in the others, for each region in the CpGs object.
#' @param group1 The name of the first group for the comparison. This is stored
#'  in colData(rrbs)
#' @param group2 The name of the second group for the comparison. This is stored
#'  in colData(rrbs)
#' @param smaller Determines whether the p-value is computed whether the
#'  test-statistic is greater or lesser than inter-replicate
#' values. For our purposes, it should be set to FALSE.
#' @param comparison Details which groups we are using to define our empirical
#'  testing distribution. The default is to
#' use all of them, however, should the user find one group contains unusually
#'  high variability, then that group can be selected.
#' Values are either 'allReps', 'Group1' or 'Group2'.
#' @param method Determines which method is used to calculate p-values.
#' 'empirical' uses the empirical distribution directly, without modelling.
#'  This is the default. 'model', fits an exponential distribution to the tail
#'  of our null distribution.
#' @param closePara Sets a threshold for how close the exponential curve should
#' fit the empirical distribution in the 'model' method. If the method produces
#' errors, consider raising this parameter.
#' @return Returns a list P, with 2 entries. 'FDRmean' is the Benjamini-Hochberg
#'  adjusted p-values. The unadjusted p-values are stored in 'Pmean'.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @export
#' @import BiSeq
#' @import GenomicRanges
#' @import IRanges
#' @examples
#' \donttest{data(rrbsDemo)
#' data(CpGsDemo)
#' data(MMDlistDemo)
#' M3D <- MMDlistDemo$Full-MMDlistDemo$Coverage
#' M3D <- cbind(rowMeans(M3D[,2:5],na.rm=TRUE),M3D[,1],M3D[,6])
#' colnames(M3D) <- c('MeanBetween','H1-hESC1  vs  H1-hESC2','K562-1  vs  K562-2')
#' group1 <- unique(colData(rrbsDemo)$group)[1]
#' group2 <-unique(colData(rrbsDemo)$group)[2]
#' PDemo <- pvals(rrbsDemo, CpGsDemo, M3D,
#'             group1, group2, smaller=FALSE,comparison='allReps')}


pvals_lite <- function(rrbs, CpGs, M3D_stat_lite, group1, group2,
                  smaller=FALSE,comparison='allReps',method='empirical',closePara=0.005){
    # get the names of the columns to compare from M3D_stat_lite object
    samples1 <- rownames(colData(rrbs))[colData(rrbs)[,]==group1]
    samples2 <- rownames(colData(rrbs))[colData(rrbs)[,]==group2]
    within1 <- determineGroupComps(samples1,type='within')
    within2 <- determineGroupComps(samples2,type='within')
#     within <- c(within1,within2)
#     between <- determineGroupComps(samples1,samples2,type='between')
#
#     ids1 <- findComps(M3D_stat_lite,within1)
#     ids2 <- findComps(M3D_stat_lite,within2)
#     ids3 <- findComps(MMD,between)
#
#     # get the MMD in order of between groups, then within groups
#     Within <- cbind(MMD[,within1],MMD[,within2])
#     colnames(Within) <- c(within1,within2)
#     Between <- MMD[,between]
#     AvsB <- cbind(Between,Within)
    Within <- M3D_stat_lite[,2:length(M3D_stat_lite[1,])]

    # find the total counts over each island, per sample
    nCpGs <- length(CpGs)
    ovlaps <- findOverlaps(CpGs,rowRanges(rrbs))
    islandList <- unique(queryHits(ovlaps))
    nIslands <- length(islandList)

    if (method=='model' | method=='model-con'){
        # fit an exponential to the tail of the curve, using the 95th percentile
        wndSze <- 0.01
        cutoff <- quantile(Within,0.95,na.rm=TRUE)
        cutoff <- floor(cutoff/wndSze)*wndSze
        # create a sliding window and count the occurences in the bins
        top <- ceiling(max(Within,na.rm=TRUE)/wndSze)*wndSze
        base <- seq(cutoff, top, wndSze)
        tot <- length(Within[!is.na(Within)])
        cdf <- unlist(lapply(base,function(i) {
            sum(Within<=i,na.rm=TRUE)/tot
        }))
        inds <- cdf!=1 # remove the log(0) terms
        cdf <- cdf[inds]
        base <- base[inds]

        # lambda method - search over lambdas, choose the
        # closest without underestimating the probabilities
        # (assume lambda between 1 and 150)

        lambdaTests <- seq(1,150,0.1)
        fit <- rep(NaN,length(lambdaTests))
        pass <- rep(TRUE,length(lambdaTests))
        countClose <- rep(TRUE,length(lambdaTests))
        for (i in 1:length(lambdaTests)){
            # sum squared fit
            l<- lambdaTests[i]
            fvect <- exp(-l*base)-1+cdf
            if (any(fvect<0)){
                pass[i] <- FALSE
            }
            countClose[i] <- sum(abs(fvect)>closePara)
            #       fit[i] <- sum(fvect*fvect)
            fit[i] <- abs(prod(fvect))
        }
        close <- which(countClose-min(countClose)<=2)
        passes <- which(pass==TRUE)
        fit <- fit[close]
        passes <- passes[close]
        lambdaTests <- lambdaTests[close]
        ind <- which(fit==min(fit,na.rm=TRUE))
        lambda <- lambdaTests[ind]
        cat('lambda = ',lambda)

        # get the pvalues (under the exponential distribution with lambda)
        # use the mean of the inter-group M3D stats
        Pmean <- lapply(1:length(CpGs), function(i) {
            testStat <- M3D_stat_lite[i,1]
            if (testStat <= 0) {
                testStat <- 0
            }
            exp(-lambda* testStat)
        })
    } else if (method=='empirical'){
        # use the empirical distribution of the test-statistic as null distribution
        Pmean <- lapply(1:nIslands, function(j) {
            if (smaller==TRUE){
                if (comparison=='Group1'){
                    length(which(M3D_stat_lite[,within1] <=M3D_stat_lite[j,1]))/
                        length(within1)/nIslands
                } else if (comparison=='Group2')  {
                    length(which(M3D_stat_lite[,within2]<=M3D_stat_lite[j,1]))/
                        length(within2)/nIslands
                } else {
                    length(which(Within <= M3D_stat_lite[j,1]))/length(Within)
                }

            } else {

                if (comparison=='Group1'){
                    length(which(M3D_stat_lite[,within1]>=M3D_stat_lite[j,1]))/
                        length(within1)/nIslands
                } else if (comparison=='Group2')  {
                    length(which(M3D_stat_lite[,within2]>=M3D_stat_lite[j,1]))/
                        length(within2)/nIslands
                } else {
                    length(which(Within >= M3D_stat_lite[j,1]))/length(Within)
                }
            }
        })
    }
    Pmean <- unlist(Pmean)
    FDRmean <- p.adjust(as.vector(Pmean),method='BH')
    P <- list(Pmean,FDRmean)
    names(P) <- c('Pmean','FDRmean')
    return(P)
}
