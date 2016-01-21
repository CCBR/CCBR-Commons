# presentation.R
# code for July 14 presentation to CCBR
# Randall Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc

# run the test code (including portions inside if(FALSE) statements) in densityScan and compareDensities first
geneInfo <- data.frame(MGI = unique(translate$MGI),
                       len = 0,
                       stringsAsFactors = FALSE)

# calculate gene lengths for all genes -- this bit is copied and modified from densityScan
for(g in geneInfo$MGI)
{
    # get all IDs
    tmp <- subset(translate, MGI == g)
    ensmu <- unique(tmp$Ensembl)

    # find chromosomal coordinates of g
    genes.sub <- subset(genes, info %in% ensmu)


    # find length of region covered by g (kb)
    # be careful if we have any overlapping regions
    overlap <- sapply(genes.sub$start, function(x) x > genes.sub$start & x < genes.sub$end)

    while(any(overlap))
    {
        # will will drop this row, but first merge it with the first region that it overlaps with
        fix <- (1:dim(genes.sub)[1])[apply(overlap, 2, any)][1]
        mergewith <- which(overlap[,fix])[1]

        genes.sub$end[mergewith] <- with(genes.sub, max(end[mergewith], end[fix]))
        genes.sub <- genes.sub[-fix,]

        # check if we have any more overlaps
        overlap <- sapply(genes.sub$start, function(x) x > genes.sub$start & x < genes.sub$end)
    }

    geneInfo$len[geneInfo$MGI == g] <- with(genes.sub, sum(end - start)) / 1000
}

geneInfo <- merge(geneInfo, data.frame(MGI = rownames(Z1),
                                       wt = as.vector(Z1), # wild type
                                       mt = as.vector(Z2), # mutant
                                       stringsAsFactors = FALSE), all = TRUE)

geneInfo$wtPeak <- with(geneInfo, !is.na(wt) & wt > 0)
geneInfo$mtPeak <- with(geneInfo, !is.na(mt) & mt > 0)

geneInfo$wt[is.na(geneInfo$wt)] <- 0
geneInfo$mt[is.na(geneInfo$mt)] <- 0

if(FALSE)
{
    save(geneInfo, file = '~/Documents/Papers & Abstracts/Randy/Presentations/15.07.13 CCBR Mtg/geneInfo.RData')
    load('~/Documents/Papers & Abstracts/Randy/Presentations/15.07.13 CCBR Mtg/geneInfo.RData')
}

# take a look at observation of peaks

wtDens <- with(geneInfo, ksmooth(len, wtPeak, bandwidth = 200))
with(geneInfo, plot(len, wtPeak))
lines(wtDens)
plot(wtDens, type = 'l', ylim = c(0, 0.06), xlim = c(0,1500))

mtDens <- with(geneInfo, ksmooth(len, mtPeak, bandwidth = 200))
with(geneInfo, plot(len, mtPeak))
lines(mtDens)
plot(mtDens, type = 'l', ylim = c(0, 0.06), xlim = c(0,1500))

# is there an association between len and peak density?
# I tried looking at log(density) as well as dropping influential points, but there are always more influential
# points that pop out...this isn't a very good system to model with linear regression
wt.dens <- glm(wt ~ len, data = subset(geneInfo, wtPeak)) # significant association, but beta is extremely small
mt.dens <- glm(mt ~ len, data = subset(geneInfo, mtPeak)) # significant association, but beta is extremely small

# is there an association between len and presence of peak?
wt.pres <- glm(wtPeak ~ len, family = binomial, data = geneInfo) # no association
mt.pres <- glm(mtPeak ~ len, family = binomial, data = geneInfo) # some association - very small beta: OR = 1.0006
