#' Computes the components of the M3D test-statistic over one region for 2
#'  samples
#'
#' Returns the two components of the M3D test-statistic - the MMD
#'  (Gretton et al. 2006) for the full data and the coverge
#' only data, respectively. This is not intended to be called directly
#'  by the user.
#'
#' @param test_data Contains the methylation data over two samples for a given
#'  region.
#' @param loc_mx The matrix of distances between the CpG sites.
#' @param locInds The indices of the non-zero entries of loc_mx.
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

M3D_Single <- function(test_data,loc_mx,locInds,method='MinusCovMMD'){

    # get the methylation data, to calculate the MMD
    # (this is how many times each one appears)
    t1 <- test_data[,1]
    t2 <- test_data[,2]
    t3 <- test_data[,3]
    t4 <- test_data[,4]
    meth11 <- tcrossprod(t1, t1)
    meth12 <- tcrossprod(t1, t2)
    meth22 <- tcrossprod(t2, t2)
    unMeth11 <- tcrossprod(t3, t3)
    unMeth12 <- tcrossprod(t3, t4)
    unMeth22 <- tcrossprod(t4, t4)
    # number of observations of 'x', first variable

    m <- sum(t1+t3)
    # number of observations of 'y', second variable
    n <- sum(t2+t4)

    if(m==0 | n==0){
        biased <- NaN
        biased_cov <- NaN
        res <- list(biased,biased_cov)
        return(res)
    }

    # get the summing weights
    weights11 <- meth11 + unMeth11
    weights12 <- meth12 + unMeth12
    weights22 <- meth22 + unMeth22

    # estimate sigma from the data
    mdist <- median_freq(loc_mx[locInds],weights12[locInds])
    sigma <- sqrt(mdist/2) # we have sigma = sqrt(mdist/2)
    if (sigma==0 | is.na(sigma)){
        sigma <- 1
    }

    # apply the rbf to the loc_mx
    #     rbfMx <- exploc_mx^(-1/2/sigma^2)
    rbfMx <- exp(-1/2/sigma^2 * loc_mx)


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
        total1 <- t1+t3
        total2 <- t2+t4

        weights11Cov <- tcrossprod(total1,total1)
        weights12Cov <- tcrossprod(total1, total2)
        weights22Cov <- tcrossprod(total2, total2)

        sumKxxCov = sum(weights11Cov * rbfMx)
        sumKyyCov = sum(weights22Cov * rbfMx)
        sumKxyCov = sum(weights12Cov * rbfMx)

        biased_cov <- sqrt(sumKxxCov/(m*m) + sumKyyCov/(n*n)
                           - 2*sumKxyCov/(m*n))
        res <- list(biased,biased_cov)
        return(res)
    }
}
