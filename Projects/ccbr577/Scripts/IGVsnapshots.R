# IGVsnapshots.R
# Generate an IGV batch file to look at all regions with peaks (just do one snapshot if they overlap)
# Randy Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc


vars <- commandArgs(TRUE)

tmp <- strsplit(vars, '=')

vars <- as.list(sapply(tmp, `[`, 2))
names(vars) <- sapply(tmp, `[`, 1)

##### testing #####
if(FALSE)
{
    vars$wtIns <- '../Data/Hughes_TAF3_WT_q20ct2.bed'
    vars$wtPeaks <- '../Results/Hughes_TAF3_WT_q20ct2_peaks.bed'
    vars$mtIns <- '../Data/Hughes_TAF3mutM880A_q20ct2.bed'
    vars$mtPeaks <- '../Results/Hughes_TAF3mutM880A_q20ct2_peaks.bed'
}


##### Defaults #####


##### Read in data #####

## wtIns <- read.table(vars$wtIns, skip = 1, stringsAsFactors = FALSE, header = FALSE)
wtPeaks <- read.table(vars$wtPeaks, skip = 1, stringsAsFactors = FALSE, header = FALSE)
names(wtPeaks) <- c('chr', 'start', 'end', 'genes', 'V5', 'strand', 'thickStart', 'thickEnd')

## mtIns <- read.table(vars$mtIns, skip = 1, stringsAsFactors = FALSE, header = FALSE)
mtPeaks <- read.table(vars$mtPeaks, skip = 1, stringsAsFactors = FALSE, header = FALSE)
names(mtPeaks) <- c('chr', 'start', 'end', 'genes', 'V5', 'strand', 'thickStart', 'thickEnd')

Peaks <- list(wt = wtPeaks,
              mt = mtPeaks)


##### Get snapshot coordinates #####

snapshots <- matrix(nrow = 0, ncol = 3, dimnames = list(character(), c("chr", "start", "end")))

allchrs <- unique(c(wtPeaks$chr, mtPeaks$chr))
chr.curr <- allchrs[1]

while(dim(Peaks$wt)[1] > 0 & dim(Peaks$mt)[1] > 0)
{
    # reset pick
    pick <- ''
    nopick <- ''

    # move to next chromosome if we need to
    if(Peaks$wt$chr[1] != chr.curr & Peaks$mt$chr[1] != chr.curr)
    {
        allchrs <- allchrs[-1]

        chr.curr <- allchrs[1]

        # this should be the next one ... could be buggy!
        if(Peaks$wt$chr[1] != allchrs[1] | Peaks$mt$chr[1] != allchrs[1])
            stop('Possible problem in sorting of bed files detected')
    }

    # find next peak start
    if(dim(Peaks$wt)[1] == 0)
    {
        pick <- 'mt'
        nopick <- 'wt'
    }

    if(dim(Peaks$mt)[1] == 0 & pick == '')
    {
        pick <- 'wt'
        nopick <- 'mt'
    }

    if(Peaks$wt$chr[1] != chr.curr & pick == '')
    {
        pick <- 'mt'
        nopick <- 'wt'
    }

    if(Peaks$mt$chr[1] != chr.curr & pick == '')
    {
        pick <- 'wt'
        nopick <- 'mt'
    }

    if(Peaks$wt$start[1] < Peaks$mt$start[1] & pick == '')
    {
        pick <- 'wt'
        nopick <- 'mt'
    }

    if(Peaks$wt$start[1] > Peaks$mt$start[1] & pick == '')
    {
        pick <- 'mt'
        nopick <- 'wt'
    }

    # identify overlaps between the two
    if(dim(Peaks[[nopick]])[1] > 0)
    {
        if(Peaks[[pick]]$end[1] > Peaks[[nopick]]$start[1])
        {
            both <- TRUE
            end <- max(c(Peaks[[nopick]]$end[1], Peaks[[pick]]$end[1]))
        }else{
            both <- FALSE
            end <- Peaks[[pick]]$end[1]
        }

        both <- FALSE
    }

    # construct location for snapshot
    snapshots <- rbind(snapshots,
                       c(chr.curr,
                         as.character(Peaks[[pick]]$start[1]),
                         as.character(end)))

    # remove included peaks
    Peaks[[pick]] <- Peaks[[pick]][-1,]

    if(both)
        Peaks[[nopick]] <- Peaks[[nopick]][-1,]
}


##### Write file #####

cat("snapshotDirectory ~/Documents/Work/CCBR-Commons/Projects/ccbr577/Results/peaks/", '\n', file = 'IGVbatch.txt')

for(i in 1:dim(snapshots)[1])
{
    cat("goto ", snapshots[i,1], ":", snapshots[i,2], "-", snapshots[i,3], "\n",
        file = 'IGVbatch.txt', sep = '', append = TRUE)
    cat("snapshot\n", file = 'IGVbatch.txt', append = TRUE)
}
