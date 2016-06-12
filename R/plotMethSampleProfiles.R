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
#' @param plot_title Optional title for the plot, overrides automatic title
#' @return NULL, the function plots the profiles
#' @export
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' # plot the 9th region in the Toy Data Set
#' data(rrbsDemo)
#' data(CpGsDemo)
#' plotMethProfile(rrbsDemo, CpGsDemo, 'H1-hESC', 'K562', 9)

plotMethSampleProfiles <- function(rrbs, CpGs, group1, group2, CpGindex, 
                                   samples=NaN, plot_title=NaN){
    Ov <- findOverlaps(CpGs[CpGindex],rowRanges(rrbs))
    CpGisland <- rowRanges(rrbs)[subjectHits(Ov)]
    if (length(CpGisland)==0){   # in case we pass an empty island
        warning('This island is empty, returning without plotting')
        return()
    }
    # STEP 1. Get the Cystosine sites
    CytoSites <- start(ranges(CpGisland))
    Ovlaps <- findOverlaps(CpGisland,rowRanges(rrbs))
    # order them to get better plots
    ord <- order(CytoSites)
    CytoSites <- CytoSites[ord]

    # default is all samples from the two groups
    if (is.na(samples)){
        samples <- which(colData(rrbs)$group %in% c(group1,group2))
    }
    
    # STEP 2. Get the methylation data. Average over the groups. Find the proportions
    # initilise lists
    sample_names <- list()
    methyl_group <- list()
    cov_group <- list()

    # get data to speed up
    colour_groups <- c('dodgerblue3', 'red')
    all_samps <- colnames(methReads(rrbs))
    meth_data <- methReads(rrbs)[subjectHits(Ovlaps)[ord],]
    cov_data <- totalReads(rrbs)[subjectHits(Ovlaps)[ord],]
    cov_test <- list()
    cov_test_sum <- FALSE
    meth_lvl <- list()
    smooth_group <- list()
    lower_meth <- 1 # initialise the lower bound for plot
    upper_meth <- 0 # and for upper bound
    # weight matrix for smoothing
    # STEP 3. Smooth it. Using the method in BiSeq paper.
    h = 80 # bandwidth
    numCs <- length(CytoSites)
    numSamples <- length(samples)
    Weights <- matrix(0,numCs,numCs) # each columns is the weights for that C
    for (col in 1:numCs){
        for (row in 1:numCs){
            if (abs(CytoSites[row]-CytoSites[col]) <= h){
                Weights[row,col] <- 1 - (abs(CytoSites[row]-CytoSites[col]))/h
            }
        }
    }
    
    for (i in 1:length(samples)){
        samp_num <- samples[i]
        sample_temp <- all_samps[samp_num]
        group_temp <- as.character(colData(rrbs)$group[samp_num])
        if(group_temp==group1){
            colour_groups[[i]] <- 'dodgerblue3'
        } else { colour_groups[[i]] <- 'red'}
        sample_names[[i]] <- sample_temp
        methyl_group[[i]] <- meth_data[, sample_temp]
        cov_group[[i]] <- cov_data[, sample_temp]
        cov_test[[i]] <- all(cov_group[[i]]==0)
        cov_test_sum <- cov_test_sum | cov_test[[i]]
        meth_lvl[[i]] <- methyl_group[[i]]/cov_group[[i]]
        smooth_group[[i]] <- t(t(methyl_group[[i]]) %*% Weights) / 
            t(t(cov_group[[i]]) %*% Weights)
        temp_range <- range(meth_lvl[[i]], na.rm = TRUE, finite = TRUE)
        lower_meth <- min(lower_meth, temp_range[1])
        upper_meth <- max(upper_meth, temp_range[2])
    }
    
    
    if (cov_test_sum==TRUE){
        message('At least one sample has no coverage: ', 
                paste(sample_names[[which(cov_test==TRUE)]], sep=', '))
        return()
    }

    
    # STEP 4 Plot.
    meth_range <- c(lower_meth, upper_meth)
    plot(CytoSites, meth_lvl[[1]], ylim = meth_range, 
         col = colour_groups[[1]], ann = FALSE)
    lines(CytoSites, smooth_group[[1]], col = colour_groups[[1]], lty=1, lwd=2)
    for (i in 2:length(samples)){
        points(CytoSites, meth_lvl[[i]], col=colour_groups[[i]], pch=22)
        lines(CytoSites, smooth_group[[i]],col=colour_groups[[i]], lty=1,lwd=2)
    }
    title(ylab='Methylation Level')
    xlabel = paste(unique(seqnames(CpGisland)),sep=', ')
    title(xlab=xlabel)
    if (is.na(plot_title)){
        title(main=paste("Island",CpGindex, group1 , 'vs', group2 ,sep=' '), font.main=3)
    } else {
        title(main=plot_title, font.main=3)
    }    #legend(min(CytoSites),1,c(as.character(group1),as.character(group2)), cex=0.8, col=c("blue","red"), pch=21:22)
}