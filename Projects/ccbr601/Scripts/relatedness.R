# relatedness.R
# Calculation of relatedness in NPHS2 cohort
# Randy Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc


library(ALDsuite)
library(ALDdata)


######### read in raw data #########

dat <- read.table('../Data/Copy of Wenkosi_chr 1FinalReport.txt', skip = 9, header = TRUE,
                  stringsAsFactors = FALSE, na.strings = '-')

dat2 <- read.table('../Data/kareshma 2_chromosome 1 FinalReport.txt', skip = 9, header = TRUE,
                   stringsAsFactors = FALSE, na.strings = '-')


# These two files partially overlap, and George wanted to add a couple individuals from the second file
# ...these individuals have some overlapping data in the first file, though.
# For now, I will do a full merge...for any overlapping mismatches I will keep the information in dat

dat <- merge(dat, dat2, all = TRUE)

# Since GC_Score is spelled GCScore in dat2, we can identify where they came from
# Case1: Neither GC_Score or GCScore are missing
#        These observations overlap and agree -- keep them!
# Case2: Either GC_Score or GCScore are missing
#        Either there is no overlap, or there was a discrepancy between the two -- drop these!
# Case3: Both GC_Score and GCScore are missing
#        We don't have any of these

# many of these might be labeled as missing in one file -- drop missings
dat <- subset(dat, !is.na(Allele1) & !is.na(Allele2))

# start by finding duplicate SampleID/SNPName combinations
dat$rowID <- with(dat, paste(SNPName, SampleID, sep = ':'))

tmp <- table(dat$rowID)
tmp <- tmp[tmp == 2] # --> there are only 9 discrepancies among 3 markers (affecting 18 rows)

# drop dups
dat <- subset(dat, !rowID %in% names(tmp))

rm(dat2)

# make sure we drop NPHS2 and the surrounding region
nphs2 <- 179535000
dat <- subset(dat, !(Chr == 1 & Position > nphs2 - 1.5e6 & Position < nphs2 + 1.5e6))


######### MAP file #########

# get unique markers
map <- unique(subset(dat, select = c(Chr, SNPName, Position)))

# calculate genetic position and reorder columns
map$genpos <- with(map, gen.calc(Chr, Position)$gen.pos) / 100 # want this in Morgans
map <- map[order(map$Position),c(1,2,4,3)]

write.table(map, file = '../Data/nphs2.map', sep = '\t', quote = FALSE, na = '0', row.names = FALSE,
            col.names = FALSE)


######### PED file #########

ped <- data.frame(family = NA,
                  indiv = unique(dat$SampleID),
                  father = NA,
                  mother = NA,
                  sex = 'other',
                  pheno = 1,
                  stringsAsFactors = FALSE)

ped$family <- paste('F', ped$indiv, sep = '')
ped$father <- paste('D', ped$indiv, sep = '')
ped$mother <- paste('M', ped$indiv, sep = '')

# create a matrix of geontypes to be cbound to ped
pedmat <- matrix(NA, nrow = dim(ped)[1], ncol = dim(map)[1],
                 dimnames = list(ped$indiv, map$SNPName))

for(j in map$SNPName)
{
    tmp <- subset(dat, SNPName == j, select = c(SampleID, Allele1, Allele2))
    pedmat[tmp$SampleID,j] <- with(tmp, paste(Allele1, Allele2, sep = '\t'))
}

pedmat[is.na(pedmat)] <- '0\t0'
tmp <- apply(pedmat, 1, paste, collapse = '\t')

write.table(cbind(ped, tmp), file = '../Data/nphs2.ped', sep = '\t', quote = FALSE, na = '0',
            row.names = FALSE, col.names = FALSE)


######### Run PLINK #########

system('./checkIBD_F')


######### People Cherie is interested in #########

## ids <- c('NC046', 'SA015', 'SA036', 'SA048', 'SA054', 'SA088', 'SA089')
## ids <- c('WNK03426', 'WNK03427', 'WNK03385', 'WNK03379', 'WNK03367', 'WNK03347', 'WNK08583')
## ids <- c('NR043226', 'NR043225', 'NR043206', 'NR043200', 'NR043188', 'NR043168', 'NR051111')
ids <- c('CL300052', 'CL300045', 'CL300049', 'CL300043', 'CL300062', 'CL300041', 'WNK03347',
         'WNK03367', 'WNK03374', 'WNK03379', 'WNK03385', 'WNK03392', 'WNK03426', 'WNK03427')


######### Look at results before sending #########

tmp <- read.table('../Results/plink.het', header = TRUE, stringsAsFactors = FALSE)

# see F for the ids above
v260eHom <- lapply(ids, function(i) tmp[grep(i, tmp$IID),'F'])
names(v260eHom) <- ids

hist(sapply(tmp$F, function(x) max(x,0)), xlab = 'F', main = 'Histogram of F')

tmp <- read.table('../Results/plink.genome', header = TRUE)
with(tmp, plot(Z0, Z1, xlim = 0:1, ylim = 0:1))
abline(1, -1)

# using the IDs below...
tmp$flag1 <- FALSE
tmp$flag1[unlist(sapply(ids, grep, x = tmp$IID1))] <- TRUE
tmp$flag2 <- FALSE
tmp$flag2[unlist(sapply(ids, grep, x = tmp$IID2))] <- TRUE

tmp$flag <- tmp$flag1 & tmp$flag2
with(subset(tmp, flag), plot(Z0, Z1, xlim = 0:1, ylim = 0:1))
abline(1, -1)
