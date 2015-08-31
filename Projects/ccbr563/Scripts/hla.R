# hla.R
# Preliminary analysis of HLA data for Denise Whitby
# Randall Johnson
# CCR Collaborative Bioinformatics Resource
# Advanced Biomedical Computing Center at Frederick National Laboratory
# Leidos Biomedical Research, Inc

library(XLConnect)
library(lpSolve)
library(WriteXLS)

#############
# Data Prep #
#############

##### HLA data from Carrington lab #####
dat <- readWorksheet(loadWorkbook('../Data/Cameroon_KS_HLA.xlsx'), sheet = 'Sheet1')
dat2 <- rbind(readWorksheet(loadWorkbook('../Data/HLAtypingResults_plate5_150824.xlsx'), sheet = 'Sheet1'),
              readWorksheet(loadWorkbook('../Data/HLAtypingResults_plate6_150824.xlsx'), sheet = 'Sheet1'))

dat <- merge(dat, dat2, by.x = names(dat),
             by.y = c('HGAL', 'PID', 'A_1', 'A_2', 'B_1', 'B_2', 'C_1', 'C_2', 'DQA1_1', 'DQA1_2', 'DQB1_1',
                      'DQB1_2', 'DRB1_1', 'DRB1_2', 'DRB3_1', 'DRB3_2', 'DRB4_1', 'DRB4_2', 'DRB5_1', 'DRB5_2',
                      'DPB1_1', 'DPB1_2'), all = TRUE)

dat$PID <- gsub('_VOS', '', dat$PID)

clin <- read.csv('../Data/ksbase.csv', stringsAsFactors = FALSE,
                 na.strings = c('Legitimate Skip', '-1', '88', '99'))
clin$case <- with(clin, substr(ID, 1, 1) == 'K')

dat <- merge(dat, subset(clin, select = c('ID', 'case')), by.x = 'PID', by.y = 'ID')

names(dat) <- gsub('Export.', '', names(dat), fixed = TRUE)
rownames(dat) <- dat$PID

##### Already have matches for these #####
matched <- readLines('hla_matched.txt')
matched <- unlist(strsplit(matched, ':', fixed = TRUE))

dat <- subset(dat, !PID %in% matched)

##### Load HLA data from my R package #####
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaA_broad.RData')
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaB_broad.RData')
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaBw.RData')
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaC_broad.RData')
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaCgrp.RData')
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaDQA_broad.RData')
load('~/Documents/Work/R-packages/HLA/HLA/data/hlaDRB1_broad.RData')


##### Figure out genotypes at different levels of accuracy #####

# allele group level (2 digit typing)
agroup <- function(var)
    sapply(strsplit(var, ":"), `[`, 1)

dat$C1A1.agroup <- agroup(dat$C1A1)
dat$C1A2.agroup <- agroup(dat$C1A2)

dat$C1B1.agroup <- agroup(dat$C1B1)
dat$C1B2.agroup <- agroup(dat$C1B2)

dat$C1C1.agroup <- agroup(dat$C1C1)
dat$C1C2.agroup <- agroup(dat$C1C2)

dat$DQA.A1.agroup <- agroup(dat$DQA.A1)
dat$DQA.A2.agroup <- agroup(dat$DQA.A2)

dat$DRB1.A1.agroup <- agroup(dat$DRB1.A1)
dat$DRB1.A2.agroup <- agroup(dat$DRB1.A2)

# broad level (old typing)
dat$C1A1.broad <- hlaA.broad[dat$C1A1.agroup]
dat$C1A2.broad <- hlaA.broad[dat$C1A2.agroup]

dat$C1B1.broad <- hlaB.broad[dat$C1B1.agroup]
dat$C1B2.broad <- hlaB.broad[dat$C1B2.agroup]

dat$C1C1.broad <- hlaC.broad[dat$C1C1.agroup]
dat$C1C2.broad <- hlaC.broad[dat$C1C2.agroup]

dat$DQA.A1.broad <- hlaDQA.broad[dat$DQA.A1.agroup]
dat$DQA.A2.broad <- hlaDQA.broad[dat$DQA.A2.agroup]

dat$DRB1.A1.broad <- hlaDRB1.broad[dat$DRB1.A1.agroup]
dat$DRB1.A2.broad <- hlaDRB1.broad[dat$DRB1.A2.agroup]

# update specific protein level (4 digit typing)
# not all individuals have 4 digit typing...for those that don't, make NA
dat$C1A1 <- with(dat, ifelse(nchar(C1A1) == 5, C1A1, NA))
dat$C1A2 <- with(dat, ifelse(nchar(C1A2) == 5, C1A2, NA))

dat$C1B1 <- with(dat, ifelse(nchar(C1B1) == 5, C1B1, NA))
dat$C1B2 <- with(dat, ifelse(nchar(C1B2) == 5, C1B2, NA))

dat$C1C1 <- with(dat, ifelse(nchar(C1C1) == 5, C1C1, NA))
dat$C1C2 <- with(dat, ifelse(nchar(C1C2) == 5, C1C2, NA))

dat$DQA.A1 <- with(dat, ifelse(nchar(DQA.A1) == 5, DQA.A1, NA))
dat$DQA.A2 <- with(dat, ifelse(nchar(DQA.A2) == 5, DQA.A2, NA))

dat$DRB1.A1 <- with(dat, ifelse(nchar(DRB1.A1) == 5, DRB1.A1, NA))
dat$DRB1.A2 <- with(dat, ifelse(nchar(DRB1.A2) == 5, DRB1.A2, NA))


##### Figure out Bw/C group status #####

# HLA-Bw4
dat$bw4.1 <- dat$C1A1.agroup %in% subset(hlaBw, hla == 'A' & bwGroup == 'Bw4')$allele
dat$bw4.2 <- dat$C1A2.agroup %in% subset(hlaBw, hla == 'A' & bwGroup == 'Bw4')$allele
dat$bw4.3 <- (dat$C1B1.agroup %in% subset(hlaBw, hla == 'B' & (bwGroup == 'Bw4' | is.na(bwGroup)))$allele |
              dat$C1B1 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw4')$allele) &
             !dat$C1B1 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw6' | is.na(bwGroup))$allele
dat$bw4.4 <- (dat$C1B2.agroup %in% subset(hlaBw, hla == 'B' & (bwGroup == 'Bw4' | is.na(bwGroup)))$allele |
              dat$C1B2 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw4')$allele) &
             !dat$C1B2 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw6' | is.na(bwGroup))$allele

dat$bw4 <- with(dat, bw4.1 | bw4.2 | bw4.3 | bw4.4)

# HLA-Bw6
dat$bw6.1 <- (dat$C1B1.agroup %in% subset(hlaBw, hla == 'B' & (bwGroup == 'Bw6' | is.na(bwGroup)))$allele |
              dat$C1B1 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw6')$allele) &
             !dat$C1B1 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw4' | is.na(bwGroup))$allele
dat$bw6.2 <- (dat$C1B2.agroup %in% subset(hlaBw, hla == 'B' & (bwGroup == 'Bw6' | is.na(bwGroup)))$allele |
              dat$C1B2 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw6')$allele) &
             !dat$C1B2 %in% subset(hlaBw, hla == 'B' & bwGroup == 'Bw4' | is.na(bwGroup))$allele

dat$bw6 <- with(dat, bw6.1 | bw6.2)

# HLA-C group 1
dat$cg1.1 <- (dat$C1C1.agroup %in% subset(hlaCgrp, group == 'C1')$allele |
              dat$C1C1 %in% subset(hlaCgrp, group == 'C1')$allele) &
             !dat$C1C1 %in% subset(hlaCgrp, group == 'C2')$allele
dat$cg1.2 <- (dat$C1C2.agroup %in% subset(hlaCgrp, group == 'C1')$allele |
              dat$C1C2 %in% subset(hlaCgrp, group == 'C1')$allele) &
             !dat$C1C2 %in% subset(hlaCgrp, group == 'C2')$allele

dat$cg1 <- with(dat, cg1.1 | cg1.2)

# HLA-C group 2
dat$cg2.1 <- (dat$C1C1.agroup %in% subset(hlaCgrp, group == 'C2')$allele |
              dat$C1C1 %in% subset(hlaCgrp, group == 'C2')$allele) &
             !dat$C1C1 %in% subset(hlaCgrp, group == 'C1')$allele
dat$cg2.2 <- (dat$C1C2.agroup %in% subset(hlaCgrp, group == 'C2')$allele |
              dat$C1C2 %in% subset(hlaCgrp, group == 'C2')$allele) &
             !dat$C1C2 %in% subset(hlaCgrp, group == 'C1')$allele

dat$cg2 <- with(dat, cg2.1 | cg2.2)


#####################
# Goodness of Match #
#####################

ncase <- sum(dat$case)
ncont <- sum(!dat$case)
nvars <- ncase * ncont

# position of cases/controls to help identify rows/columns associated with the individual
dat$pos[dat$case] <- 1:ncase
dat$pos[!dat$case] <- 1:ncont


##### Set up integer linear program variables #####

# vector to maximize (i.e. goodness of match)
c <- rep(0, nvars)

# contraints (i.e. make sure we match each case to one unique control)
A <- matrix(0, nrow = ncase + ncont, ncol = nvars,
            dimnames = list(1:(ncase + ncont), 1:nvars))

# right hand side of contraints (i.e. one control per case, one or fewer matches per control)
b <- rep(1, ncase + ncont)


##### Assign values to c and A #####
compare.two.genes <- function(a1, a2, b1, b2)
{
    # individual in {a, b}
    # allele in {1, 2}

    comp <- c(sum((a1 == b1) + (a2 == b2), na.rm = TRUE),
              sum((a1 == b2) + (a2 == b1), na.rm = TRUE))

    return(max(comp))
}

for(i in 1:ncase)
{
    # label row for the ith case
    cse <- which(dat$pos == i & dat$case)
    rownames(A)[i] <- dat$PID[cse]

    # we want one control for this case
    A[i,((i-1)*ncont + 1):(i*ncont)] <- 1

    for(j in 1:ncont)
    {
        # we only want to match this control with at most one case
        A[ncase + j,(i-1)*ncont + j] <- 1

        # this is the column we are working with
        column <- (i - 1)*ncont + j

        # label column for the i:jth pairing
        cnt <- which(dat$pos == j & !dat$case)
        colnames(A)[column] <- with(dat, paste(PID[cse], PID[cnt], sep = ":"))
        names(c)[column] <- with(dat, paste(PID[cse], PID[cnt], sep = ":"))

        # Broad match score
        c[column] <- c[column] + with(dat, sum(compare.two.genes(C1A1.broad[cse], C1A2.broad[cse],
                                                                 C1A1.broad[cnt], C1A2.broad[cnt]),
                                               compare.two.genes(C1B1.broad[cse], C1B2.broad[cse],
                                                                 C1B1.broad[cnt], C1B2.broad[cnt]),
                                               compare.two.genes(C1C1.broad[cse], C1C2.broad[cse],
                                                                 C1C1.broad[cnt], C1C2.broad[cnt]),
                                               compare.two.genes(DRB1.A1.broad[cse], DRB1.A2.broad[cse],
                                                                 DRB1.A1.broad[cnt], DRB1.A2.broad[cnt]),
                                               compare.two.genes(DQA.A1.broad[cse], DQA.A2.broad[cse],
                                                                 DQA.A1.broad[cnt], DQA.A2.broad[cnt]),
                                               na.rm = TRUE))

        # allele group match bonus
        c[column] <- c[column] + 0.5 * with(dat, sum(compare.two.genes(C1A1.agroup[cse], C1A2.agroup[cse],
                                                                       C1A1.agroup[cnt], C1A2.agroup[cnt]),
                                                     compare.two.genes(C1B1.agroup[cse], C1B2.agroup[cse],
                                                                       C1B1.agroup[cnt], C1B2.agroup[cnt]),
                                                     compare.two.genes(C1C1.agroup[cse], C1C2.agroup[cse],
                                                                       C1C1.agroup[cnt], C1C2.agroup[cnt]),
                                                     compare.two.genes(DRB1.A1.agroup[cse], DRB1.A2.agroup[cse],
                                                                       DRB1.A1.agroup[cnt], DRB1.A2.agroup[cnt]),
                                                     compare.two.genes(DQA.A1.agroup[cse], DQA.A2.agroup[cse],
                                                                       DQA.A1.agroup[cnt], DQA.A2.agroup[cnt]),
                                                     na.rm = TRUE))

        # specific protein match bonus
        c[column] <- c[column] + 0.25 * with(dat, sum(compare.two.genes(C1A1[cse], C1A2[cse],
                                                                        C1A1[cnt], C1A2[cnt]),
                                                      compare.two.genes(C1B1[cse], C1B2[cse],
                                                                        C1B1[cnt], C1B2[cnt]),
                                                      compare.two.genes(C1C1[cse], C1C2[cse],
                                                                        C1C1[cnt], C1C2[cnt]),
                                                      compare.two.genes(DRB1.A1[cse], DRB1.A2[cse],
                                                                        DRB1.A1[cnt], DRB1.A2[cnt]),
                                                      compare.two.genes(DQA.A1[cse], DQA.A2[cse],
                                                                        DQA.A1[cnt], DQA.A2[cnt]),
                                                      na.rm = TRUE))

        # Bw4/6 or C group mismatch penalty
        c[column] <- c[column] - 2 * with(dat, sum(bw4[cse] != bw4[cnt],
                                                   bw6[cse] != bw6[cnt],
                                                   cg1[cse] != cg1[cnt],
                                                   cg2[cse] != cg2[cse],
                                                   na.rm = TRUE))
    }
}


#########
# Match #
#########

set.seed(928374)
matches <- lp(direction = 'max', objective.in = c, const.mat = A,
              ## const.dir = c(rep("==", ncase), rep("<=", ncont)),
              const.dir = rep("<=", ncase + ncont), # no more than one match per person
              const.rhs = b, all.bin = TRUE)

hist(matches$objective[matches$solution == 1])
summary(matches$objective[matches$solution == 1])

make.one.sheet <- function(mats)
{
    cases <- sapply(strsplit(mats, ":"), `[`, 1)
    conts <- sapply(strsplit(mats, ":"), `[`, 2)

    retval <- subset(dat[as.vector(rbind(cases, conts)),],
                     select = c('PID', 'C1A1', 'C1A2', 'C1B1', 'C1B2', 'C1C1', 'C1C2',
                         'DQA.A1', 'DQA.A2', 'DRB1.A1', 'DRB1.A2', 'C1A1.broad', 'C1A2.broad',
                         'C1B1.broad', 'C1B2.broad', 'C1C1.broad', 'C1C2.broad', 'DQA.A1.broad',
                         'DQA.A2.broad', 'DRB1.A1.broad', 'DRB1.A2.broad', 'bw4', 'bw6', 'cg1',
                         'cg2'))
    retval$score <- NA
    retval$score[1:length(mats) * 2 - 1] <- c[mats]

    return(retval)
}

# matching of actual cases and controls ~ 1:1 matching
actual <- with(matches, objective[solution == 1])
actual.df <- make.one.sheet(names(actual))
picked.df <- make.one.sheet(names(actual[order(-actual)][1:16]))

write.table(names(actual), 'hla_matched.txt', row.names = FALSE, col.names = FALSE,
            quote = FALSE, append = TRUE)
