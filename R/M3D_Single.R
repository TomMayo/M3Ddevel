#' Computes the components of the M3D test-statistic over one region for 2
#'  samples
#' 
#' Returns the two components of the M3D test-statistic - the MMD
#'  (Gretton et al. 2006) for the full data and the coverge
#' only data, respectively. This is not intended to be called directly
#'  by the user.
#' 
#' @param testData Contains the methylation data over two samples for a given
#'  region.
#' @param locMx The matrix of distances between the CpG sites.
#' @param locInds The indices of the non-zero entries of locMx.
#' @param method This specifies whether to return the full MMD and the
#'  methylation blind MMD (focusing only on the coverage)
#' or just the former. if method ='MinusCovMMD' it is both, all other values
#'  return just the full MMD.
#' @return This returns the value of the MMD for the region between the two
#'  samples as a numeric. If method
#' is set to 'MinusCovMMD', a list is returned of the full MMD and the coverage
#'  only MMD. Subtracting them gives the M3D test-statistic.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @references Gretton, A., Borgwardt, K. M., Rasch, M., Scholkopf, B., Smola,
#'  A. J. (2006). A kernel method for the two-sample-problem. In Advances in
#'   neural information processing systems (pp. 513-520).

M3D_Single <- function(testData,locMx,locInds,method='MinusCovMMD'){
    
    # get the methylation data, to calculate the MMD 
    # (this is how many times each one appears)
    meth11 <- testData[,1]%*%t(testData[,1])
    meth12 <- testData[,1]%*%t(testData[,2])
    meth22 <- testData[,2]%*%t(testData[,2])
    unMeth11 <- testData[,3]%*%t(testData[,3])
    unMeth12 <- testData[,3]%*%t(testData[,4])
    unMeth22 <- testData[,4]%*%t(testData[,4])
    # number of observations of 'x', first variable
    m <- sum(testData[,1]+testData[,3])
    # number of observations of 'y', second variable
    n <- sum(testData[,2]+testData[,4])
    
    if(m==0 | n==0){
        biased <- NaN
        biasedCov <- NaN
        res <- list(biased,biasedCov)
        return(res)
    }
    
    # get the summing weights
    weights11 <- meth11 + unMeth11
    weights12 <- meth12 + unMeth12
    weights22 <- meth22 + unMeth22
    
    # estimate sigma from the data
    mdist <- medianFreq(locMx[locInds],weights12[locInds]) 
    sigma <- sqrt(mdist/2) # we have sigma = sqrt(mdist/2) 
    if (sigma==0 | is.na(sigma)){
        sigma <- 1
    }
    
    # apply the rbf to the locMx
#     rbfMx <- explocMx^(-1/2/sigma^2)  
    rbfMx <- exp(-1/2/sigma^2 * locMx)  

    
    # now sum, biased one first
    sumKxx = sum(weights11 * rbfMx)
    sumKyy = sum(weights22 * rbfMx)
    sumKxy = sum(weights12 * rbfMx)
    
    # biased estimate
    biased <- sqrt(sumKxx/(m*m) + sumKyy/(n*n) -2*sumKxy/(m*n) )
    
    if (method!='MinusCovMMD'){
        biased <- sqrt(sumKxx/(m*m) + sumKyy/(n*n) -2*sumKxy/(m*n))
        return(biased)
    } else {
        total1 <- testData[,1]+testData[,3]
        total2 <- testData[,2]+testData[,4]
        
        weights11Cov <- total1%*%t(total1)
        weights12Cov <- total1%*%t(total2)
        weights22Cov <- total2%*%t(total2)
        
        sumKxxCov = sum(weights11Cov * rbfMx)
        sumKyyCov = sum(weights22Cov * rbfMx)
        sumKxyCov = sum(weights12Cov * rbfMx)
        
        biasedCov <- sqrt(sumKxxCov/(m*m) + sumKyyCov/(n*n)
                          - 2*sumKxyCov/(m*n))
        res <- list(biased,biasedCov)
        return(res)
    }
}
