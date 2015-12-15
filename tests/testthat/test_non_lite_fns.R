library('M3Ddevel')
context('Non-lite functions')

test_that('pvals works properly', {
    data('rrbsDemo')
    data('CpGsDemo')
    data('MMDlistDemo')
    data('PDemo')
    group1 <- 'H1-hESC'
    group2 <- 'K562'
    M3D <- MMDlistDemo$Full-MMDlistDemo$Coverage
    expect_equal(pvals(rrbsDemo, CpGsDemo, M3D, group1, group2,
                      smaller=FALSE,comparison='allReps',method='empirical',
                      closePara=0.005),PDemo, tolerance=1.0e-4)
})
