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


######### MAP file #########

# get unique markers
map <- unique(subset(dat, select = c(Chr, SNPName, Position)))

# calculate genetic position and reorder columns
map$genpos <- with(map, gen.calc(Chr, Position)$gen.pos) / 100 # want this in Morgans
map <- map[order(map$Position),c(1,2,4,3)]

write.table(map, file = '../Data/nphs2.map', sep = '\t', quote = FALSE, na = '0', row.names = FALSE,
            col.names = FALSE)


######### PED file #########

ped <- data.frame(family = unique(dat$SampleID),
                  indiv = NA,
                  father = NA,
                  mother = NA,
                  sex = 'other',
                  pheno = 1,
                  stringsAsFactors = FALSE)

ped$indiv <- paste('I', ped$family, sep = '')
ped$father <- paste('D', ped$indiv, sep = '')
ped$mother <- paste('M', ped$indiv, sep = '')

# this will be painfully slow, but it will save a lot of coding/thinking time.... :P
for(j in map$SNPName)
{
    tmp <- subset(dat, SNPName == j, select = c(SampleID, Allele1, Allele2))
    names(tmp) <- c('SampleID', paste(j, 1:2, sep = '.'))

    # merge by familyID so that the columns remain ordered correctly...minor point since they are all unique
    ped <- merge(ped, tmp, all = TRUE, by.x = 'family', by.y = 'SampleID')
}

write.table(ped, file = '../Data/nphs2.ped', sep = '\t', quote = FALSE, na = '0', row.names = FALSE,
            col.names = FALSE)


######### Run PLINK #########

system('./checkIBD_F')


######### Look at results before sending #########

tmp <- read.table('../Results/plink.het', header = TRUE)

hist(sapply(tmp$F, function(x) max(x,0)), xlab = 'F', main = 'Histogram of F')
