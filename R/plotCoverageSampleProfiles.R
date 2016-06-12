#' Plots coverage profiles over a specific region
#' 
#' Plots the profile for each of the two testing groups. 
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

plotCoverageSampleProfiles <- function(rrbs, CpGs, group1, group2, CpGindex, 
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
    
    # default is all samples from the two groups
    if (is.na(samples)){
        samples <- which(colData(rrbs)$group %in% c(group1,group2))
    }
    
    # STEP 2. Get the methylation data. Average over the groups. Find the proportions
    # initilise lists
    sample_names <- list()
    cov_group <- list()
    
    # get data to speed up
    colour_groups <- c('dodgerblue3', 'red')
    all_samps <- colnames(methReads(rrbs))
    cov_data <- totalReads(rrbs)[subjectHits(Ovlaps),]

    upper_cov <- 0 # and for upper bound
   
    
    for (i in 1:length(samples)){
        samp_num <- samples[i]
        sample_temp <- all_samps[samp_num]
        group_temp <- as.character(colData(rrbs)$group[samp_num])
        if(group_temp==group1){
            colour_groups[[i]] <- 'dodgerblue3'
        } else { colour_groups[[i]] <- 'red'}
        sample_names[[i]] <- sample_temp
        cov_group[[i]] <- cov_data[, sample_temp]
        temp_range <- range(cov_group[[i]], na.rm = TRUE, finite = TRUE)
        upper_cov <- max(upper_cov, temp_range[2])
    }

    
    # STEP 4 Plot.
    
    plot(CytoSites, cov_group[[1]], ylim = c(0, upper_cov), 
         col = colour_groups[[1]], ann = FALSE)
    lines(CytoSites, cov_group[[1]], col = colour_groups[[1]], lty=1, lwd=2)
    for (i in 2:length(samples)){
        points(CytoSites, cov_group[[i]], col=colour_groups[[i]], pch=22)
        lines(CytoSites, cov_group[[i]],col=colour_groups[[i]], lty=1,lwd=2)
    }
    title(ylab='Coverage')
    xlabel = paste(unique(seqnames(CpGisland)),sep=', ')
    title(xlab=xlabel)
    if (is.na(plot_title)){
        title(main=paste("Island",CpGindex, group1 , 'vs', group2 ,sep=' '), font.main=3)
    } else {
        title(main=plot_title, font.main=3)
    }    #legend(min(CytoSites),1,c(as.character(group1),as.character(group2)), cex=0.8, col=c("blue","red"), pch=21:22)
}