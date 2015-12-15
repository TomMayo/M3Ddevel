\name{M3D-package}
\alias{M3D-package}
\alias{M3D}
\docType{package}
\title{
Non-parametric statistical testing
}
\description{
This package identifies statistically significantly differentially
    methylated regions of CpGs. It uses kernel methods, sepcifically the Maxmimum Mean
    Discrepancy (Gretton et al. 2006), to measure differences in methylation profiles, and relates
    these to inter-replicate changes, whilst accounting for variation in
    coverage profiles.
}
\details{
\tabular{ll}{
Package: \tab M3D\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-07-17\cr  
License: \tab Artistic License 2.0 \cr
}
This package works on RRBS data as processed by the BiSeq package. The starting point is an rrbs object, a class defined by the BiSeq package (Hebestreit et al. 2013), and a GRanges object outlining the regions to test. The maximum mean discrepancy (MMD) (Gretton et al. 2006) is calculated over each region for each pair of samples, once with respect to methylation levels and once respecting only coverage. These two values are subtracted to form a test-statistic and between-group values are compared to inter-replicate values to provide p-values. These reflect the empirical probability of observing the between-group methylation differences among the replicates.

Function list:

determineGroupComps: returns a vector of the sample comparisons
findComps: returns the indices of the M3D test-statistic that corresponding to particular samples
M3D_Single: Computes the two components of the M3D test-statistic over 1 island for 1 sample pair.
M3D_Wrapper: Computes the two components of the M3D test-statistic over all sample pairs over all islands.
medianFreq: Returns the median of data summarised by unique values and the frequency with which they occur. 
pvals: Returns empirical p-values for the regions based on the M3D test-statistic.
}
\author{
Tom Mayo

Maintainer: Tom Mayo <t.mayo@ed.ac.uk>
}
\references{
Gretton, A., Borgwardt, K. M., Rasch, M., Scholkopf, B., Smola, A. J. (2006). A kernel method for the two-sample-problem. In Advances in neural information processing systems (pp. 513-520).

Hebestreit, K., Dugas, M., Klein, H. U. (2013). Detection of significantly differentially methylated regions in targeted bisulfite sequencing data. Bioinformatics, 29(13), 1647-1653.
}

\keyword{ package }
%\seealso{
%\code{\link[BiSeq:BiSeq-package]{BiSeq}}
%}

