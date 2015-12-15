library('M3Ddevel')
context('Hidden functions')

test_that('determineGroupComps works properly', {
    samples1 <- c('a','b','c','d')
    ans<-c("a  vs  b", "a  vs  c", "a  vs  d", "b  vs  c", "b  vs  d", "c  vs  d")
    expect_equal(M3Ddevel:::determineGroupComps(samples1,type='within'),ans)
    samples2 <- c('y','z')
    ans<-c("a  vs  y", "a  vs  z", "b  vs  y", "b  vs  z", "c  vs  y", "c  vs  z",
           "d  vs  y", "d  vs  z")
    expect_equal(M3Ddevel:::determineGroupComps(samples1,samples2,'between'),ans)
    samples1 <- 'a'
    samples2 <- 'b'
    ans <- 'a  vs  b'
    expect_equal(M3Ddevel:::determineGroupComps(samples1,samples2,'between'), ans)
})

test_that('findComps works properly', {
    data('MMDlistDemo')
    expect_equal(M3Ddevel:::findComps(MMDlistDemo[[1]],'K562-1  vs  K562-2'),6)
    expect_equal(M3Ddevel:::findComps(MMDlistDemo[[1]],'H1-hESC1  vs  H1-hESC2'),1)
    samples <- c("H1-hESC1  vs  H1-hESC2", "H1-hESC1  vs  K562-1",
                 "H1-hESC1  vs  K562-2", "H1-hESC2  vs  K562-1",
                 "H1-hESC2  vs  K562-2", "K562-1  vs  K562-2")
    expect_equal(M3Ddevel:::findComps(MMDlistDemo[[1]],samples),1:6)
})

test_that('M3D_Single works properly', {
    library('BiSeq')
    data('rrbsDemo')
    data('CpGsDemo')
    rrbs <- rrbsDemo
    overlaps <- GenomicRanges::findOverlaps(CpGsDemo,rrbs)
    numPairs <- 6
    nSamples <- 4
    a <- unlist(lapply(1:(nSamples-1), function(i) rep(i,(nSamples-i))))
    b <- unlist(lapply(1:(nSamples-1), function(i) (i+1):nSamples))
    samplesIdx <- cbind(a,b)

    methIndices <- GenomicRanges::subjectHits(overlaps[
        GenomicRanges::queryHits(overlaps)==1])

    sample1 <- "H1-hESC1"
    sample2 <- "H1-hESC2"
    methData <- BiSeq::methReads(rrbs)[methIndices,c(sample1,sample2)]
    totalData <- BiSeq::totalReads(rrbs)[methIndices,c(sample1,sample2)]
    unmethData <- totalData-methData
    testData <- data.frame(meth=methData,unmeth=unmethData)
    CSites <- rowRanges(rrbs)
    locs <- start(ranges(CSites[methIndices]))
    locs <- locs - min(locs)
    # make location matrix
    G <- locs%*%t(locs)
    L <- ncol(G)
    nor <- rep(G[seq(1,L^2,L+1)],L)
    # locMx is the squared distance between sites
    locMx <- -2*G + matrix(nor, nrow=L, byrow=TRUE) + matrix(nor, nrow=L)
    locInds <- which(locMx!=0)
    expect_equal(M3Ddevel:::M3D_Single(testData,locMx,locInds,method='MinusCovMMD'),
                list(0.133288,0.1324517), tolerance=1.0e-3)

    sample1 <- "K562-1"
    sample2 <- "H1-hESC2"
    methData <- methReads(rrbs)[methIndices,c(sample1,sample2)]
    totalData <- totalReads(rrbs)[methIndices,c(sample1,sample2)]
    unmethData <- totalData-methData
    testData <- data.frame(meth=methData,unmeth=unmethData)
    expect_equal(M3Ddevel:::M3D_Single(testData,locMx,locInds,method='MinusCovMMD'),
                list(0.2358175,0.2618952), tolerance=1.0e-3)

    methIndices <- subjectHits(overlaps[queryHits(overlaps)==15])
    locs <- start(ranges(CSites[methIndices]))
    locs <- locs - min(locs)
    # make location matrix
    G <- locs%*%t(locs)
    L <- ncol(G)
    nor <- rep(G[seq(1,L^2,L+1)],L)
    # locMx is the squared distance between sites
    locMx <- -2*G + matrix(nor, nrow=L, byrow=TRUE) + matrix(nor, nrow=L)
    locInds <- which(locMx!=0)
    sample1 <- "H1-hESC1"
    sample2 <- "K562-2"
    methData <- methReads(rrbs)[methIndices,c(sample1,sample2)]
    totalData <- totalReads(rrbs)[methIndices,c(sample1,sample2)]
    unmethData <- totalData-methData
    testData <- data.frame(meth=methData,unmeth=unmethData)
    expect_equal(M3Ddevel:::M3D_Single(testData,locMx,locInds,method='MinusCovMMD'),
                list(0.06645972,0.06732423), tolerance=1.0e-3)

    methIndices <- subjectHits(overlaps[queryHits(overlaps)==121])
    locs <- start(ranges(CSites[methIndices]))
    locs <- locs - min(locs)
    # make location matrix
    G <- locs%*%t(locs)
    L <- ncol(G)
    nor <- rep(G[seq(1,L^2,L+1)],L)
    # locMx is the squared distance between sites
    locMx <- -2*G + matrix(nor, nrow=L, byrow=TRUE) + matrix(nor, nrow=L)
    locInds <- which(locMx!=0)
    sample1 <- "K562-1"
    sample2 <- "K562-2"
    methData <- methReads(rrbs)[methIndices,c(sample1,sample2)]
    totalData <- totalReads(rrbs)[methIndices,c(sample1,sample2)]
    unmethData <- totalData-methData
    testData <- data.frame(meth=methData,unmeth=unmethData)
    expect_equal(M3Ddevel:::M3D_Single(testData,locMx,locInds,method='MinusCovMMD'),
                list(0.06144124,0.06245846), tolerance=1.0e-3)
})

test_that('medianFreq works properly', {
    expect_equal(M3Ddevel:::medianFreq(1:10,1:10),7)
    expect_equal(M3Ddevel:::medianFreq(1:10,rep(1,10)),5.5)
    expect_equal(M3Ddevel:::medianFreq(2:6,c(1,2,1,45,2)),5)
})

test_that('median_freq works properly', {
    expect_equal(M3Ddevel:::median_freq(1:10,1:10),7)
    expect_equal(M3Ddevel:::median_freq(1:10,rep(1,10)),5.5)
    expect_equal(M3Ddevel:::median_freq(2:6,c(1,2,1,45,2)),5)
})


