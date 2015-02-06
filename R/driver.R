# driver.R
# power calculations for Whitby lab
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# SAIC-Frederick, Inc
# Created August 7, 2013
# Last Modified September 26, 2014

library(Epi)
library(MASS)
library(Rtools)
library(gdata)

# Transmission
power.cont.tab(n1 = 74, n2 = 228, p2 = .2, sig.level = 0.05 / 10000, power = 0.8)

# KSHV
power.cont.tab(n1 = 304, n2 = 802, p2 = 0.2, sig.level = 0.05 / 10000, power = 0.8)


########
# Data #
########

# Read in data
full <- read.csv('data/ksbase.csv', na.strings = c('', 'Legitimate Skip', '-1', '88', '99'))
dat <- full

# removing "Other, specify" variables
dat$A4a <- NULL
dat$A5a <- NULL
dat$A6i <- NULL
dat$A23a <- NULL
dat$A24a <- NULL
dat$A25a <- NULL
dat$A26a <- NULL
dat$A29A <- NULL
dat$A32A <- NULL
dat$B1e <- NULL
dat$B2g <- NULL
dat$B3a1 <- NULL
dat$B3d1 <- NULL
dat$B4f <- NULL
dat$B4i <- NULL
dat$B4f1 <- NULL
dat$B5b10 <- NULL
dat$B5b12 <- NULL
dat$B5e <- NULL
dat$B7a <- NULL
dat$B15a <- NULL
dat$B16a <- NULL
dat$B20i <- NULL
dat$B26g <- NULL
dat$B26h <- NULL
dat$B26j1 <- NULL
dat$B26k1 <- NULL
dat$B26l1 <- NULL
dat$B26l1a <- NULL
dat$B33g <- NULL
dat$C1i <- NULL
dat$C6i <- NULL
dat$C8i <- NULL
dat$C9a1 <- NULL
dat$D1a <- NULL
dat$E1i <- NULL
dat$G4i <- NULL

# excluded from the analysis
dat$A2 <- NULL # inverviewer
dat$A1 <- NULL # interview date
## dat$B3d <- NULL # Type of measure - size of the drink from the previous question (see email on 1/2/14)
dat$B9a9 <- NULL # specify
dat$B26c1a <- NULL # Pills: NA (??)
dat$B26f1a <- NULL # Injectinos: NA (??)
dat$B26i1a <- NULL # hormones IUD: NA (??)
dat$C3a <- NULL # Which other skin diseases (not sure what to do about this one??)
dat$C7a <- NULL # What was the reason for your longest stay at the hospital? (These are quite variable and will have very little overlap between individuals. Not sure they are going to give us a lot of information unless we collapse some of the categories.)
dat$C8b1 <- NULL # same as previous 2
dat$D1b <- NULL # same as previous 3
dat$G5 <- NULL # appointment for followup
dat$G6 <- NULL # validated by
dat$G7 <- NULL # interviewer code

# "didn't answer" or "don't know" questions
dat$A6h <- NULL
dat$A17h <- NULL
dat$A27a <- NULL
dat$B2f <- NULL
dat$B4h <- NULL
dat$B5b11 <- NULL
dat$B6e <- NULL
dat$B26f <- NULL
dat$B33f <- NULL
dat$C1h <- NULL
dat$C3a2 <- NULL
dat$C5h <- NULL
dat$C6g <- NULL
dat$C6h <- NULL
dat$C6b2 <- NULL
dat$C6b3 <- NULL
dat$C6c4 <- NULL
dat$C6c5 <- NULL
dat$C7a1 <- NULL
dat$C8h <- NULL
dat$C8d2 <- NULL
dat$E1g <- NULL
dat$E1h <- NULL
dat$E1b5 <- NULL
dat$E1b6 <- NULL
dat$G4h <- NULL

# dates (see 'VariableNames.xls' for names...too many to bother changing right now)
mos <- c(paste(0, 1:9, sep = ''), 10:12)
names(mos) <- c('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP',
                'OCT', 'NOV', 'DEC')
make.date <- function(x, drop = '31DEC1959:23:59:59')
{
    x <- as.character(x)
    retval <- paste(substr(x, 6, 9), mos[substr(x, 3, 5)], substr(x, 1, 2), sep = '-')
    return(as.numeric(cal.yr(ifelse(x == drop, NA, retval))))
}

#dat$A1 <- make.date(dat$A1) # interview date
dat$A4 <- make.date(dat$A4) # birth date
dat$C8d1 <- make.date(dat$C8d1)
dat$F1 <- make.date(dat$F1)
dat$G3 <- make.date(dat$G3)

# Female specific variables -> make NA for males
dat$A19 <- ifelse(dat$A3 == 1, NA, dat$A19)
dat$A19a <- ifelse(dat$A3 == 1, NA, dat$A19a)
dat$A19b <- ifelse(dat$A3 == 1, NA, dat$A19b)
dat$A19c <- ifelse(dat$A3 == 1, NA, dat$A19c)
dat$B23 <- ifelse(dat$A3 == 1, NA, dat$B23)
dat$B26a <- ifelse(dat$A3 == 1, NA, dat$B26a)
dat$B26b <- ifelse(dat$A3 == 1, NA, dat$B26b)
dat$B26c <- ifelse(dat$A3 == 1, NA, dat$B26c)
dat$B26d <- ifelse(dat$A3 == 1, NA, dat$B26d)
dat$B26e <- ifelse(dat$A3 == 1, NA, dat$B26e)

# Male specific variables -> make NA for females
dat$B22 <- ifelse(dat$A3 == 2, NA, dat$B22)

# convert categorical variables into dummy variables
dummify.me <- quote({
    for(l in unique(dat[[var]])) # cover observed levels
    {
        if(!is.na(l))
            dat[[paste(var, l, sep = '_')]] <- dat[[var]] == l
    }

    dat[[var]] <- NULL
})

var <- 'A5'; eval(dummify.me)
var <- 'A7'; eval(dummify.me)
var <- 'A23'; eval(dummify.me)
var <- 'A24'; eval(dummify.me)
var <- 'A25'; eval(dummify.me)
var <- 'A26'; eval(dummify.me)
var <- 'A29'; eval(dummify.me)
var <- 'A32'; eval(dummify.me)
var <- 'B7'; eval(dummify.me)
var <- 'B15'; eval(dummify.me)
var <- 'B16'; eval(dummify.me)
var <- 'B18a'; eval(dummify.me)
var <- 'B19'; eval(dummify.me)
var <- 'C9a'; eval(dummify.me)
var <- 'F6'; eval(dummify.me)
var <- 'G2'; eval(dummify.me)

# KS variables
# Ovall - 1 = KSHV infected
#       - 2 = KSHV uninfected
# K8_1_OD = some score for K8 region of KSHV??
# K8 = Negative or positive for K8 region of KSHV?
# ORF73_OD = some score for ORF73 region of KSHV?
# ORF73_Score = Negative or positive for ORF73 region of KSHV?
# Orf = exactly the same as ORF73_Score
dat$Orf <- NULL

# convert to numeric
dat$K8_1_Score <- as.numeric(dat$K8_1_Score)
dat$ORF73_Score <- as.numeric(dat$ORF73_Score)

# don't know what these are
dat$UniqueKey <- NULL
dat$RECSTATUS <- NULL
dat$ExperimentID <- NULL
dat$Treat <- NULL
dat$Plate_ID <- NULL
dat$Test_Date <- NULL
dat$Batch <- NULL
dat$Sample_ID <- NULL
dat$Comments <- NULL
dat$Sort <- NULL

# found a few incorrect KS diagnosis dates...adding dummy data for now
dat$F1[is.na(dat$F1) & substr(dat$ID, 1, 1) == 'K'] <- 0
dat$F1[!is.na(dat$F1) & substr(dat$ID, 1, 1) == 'T'] <- NA

####################
# Cluster Analysis #
####################

# comparing indivdiuals to individuals over all variables
# -> eigenvectors measure how similar individuals are to eachother in different dimensions
cormat <- cor(t(dat[,!names(dat) %in% c('ID', 'F1', 'G1')]), use = 'pairwise.complete.obs')

eig <- eigen(cormat)
rownames(eig$vector) <- dat$ID

set.seed(345789)
randord <- sample(1:dim(dat)[1])

# color by HIV status (no clear correlation here)
parcoord(eig$vectors[randord,1:10], col = dat$G1[randord])
plot(eig$vectors[,1], eig$vectors[,4], col = dat$G1)
plot(eig$vectors[,1], jitter(rep(0, dim(eig$vectors)[1])), col = dat$G1)

# color by KS status
parcoord(eig$vectors[randord,1:10], col = is.na(dat$F1[randord]) + 1)
plot(eig$vectors[,1], jitter(rep(0, dim(eig$vectors)[1])), col = is.na(dat$F1[randord]) + 1)

# first PC accounts for a lot of the variation
(cumsum(eig$values) / sum(abs(eig$values)))[1:10]

dat$match <- eig$vector[,1]

write.csv(eig$vector[,1], file = 'cluster.csv')


##########################
# Select pool of samples #
##########################

######### First pass #########

# hiv negative / KS positive individuals .. among Males
cases <- as.character(subset(dat, G1 == 2 & !is.na(F1) & A3 == 1)$ID)

# hiv positive / KS positive individuals .. among Males
cases <- c(cases, as.character(subset(dat, G1 == 1 & !is.na(F1) & A3 == 1)$ID))

dat$cases <- dat$ID %in% cases

# take a look at where the cases cluster (black = KS-, red = KS+, circle=Female, diamond = Male)
with(subset(dat, !cases), plot(match, jitter(rep(0, length(match))),
                                       col = 'black', pch = ifelse(A3 == 1, 18, 20)))
with(subset(dat, cases), points(match, jitter(rep(0, length(match))),
                                        col = 'red', pch = ifelse(A3 == 1, 18, 20)))
abline(v = -0.02915) # less than this is group 1


# group 3 is between these two
abline(v = -0.0289)
abline(v = -0.0287)

# care needs to be taken with the remaining two cases!
# match on HIV status explicitly on this last group because all KS+ in gorup3 are HIV+
dat$pool <- ifelse(dat$match < -0.02915 & dat$A3 == 1, 1,
            ifelse(dat$match < -0.02890 & dat$A3 == 1, 2,
            ifelse(dat$match < -0.02870 & dat$A3 == 1 & dat$G1 == 1, 3, NA)))

dat$HIV <- dat$G1 == 1
dat$KS <- !is.na(dat$F1)

# see what we get
table(dat$KS, dat$pool, dat$HIV)

# sort and write file
dat <- dat[order(dat$pool, dat$HIV, dat$KS),]

write.csv(subset(dat, select = c('ID', 'HIV', 'KS', 'pool')), file = 'pools.csv', row.names = FALSE)

######### Second pass #########

# drop KSHV- indivdiuals
dat$pool <- ifelse(dat$Ovall == 2, NA, dat$pool)

# drop individuals who are obvious outliers in the previous PCA
# -- take everything less than -0.0287, male, KSHV+
dat.sub <- subset(dat, match < -0.0287 & A3 == 1 & Ovall == 1)

cormat <- cor(t(dat.sub[,!names(dat.sub) %in% c('ID', 'F1', 'G1', 'match', 'cases', 'pool', 'HIV', 'KS')]),
              use = 'pairwise.complete.obs')

eig <- eigen(cormat)
rownames(eig$vector) <- dat.sub$ID
dat.sub$match2 <- eig$vectors[,1]

set.seed(345789)
randord <- sample(1:dim(dat.sub)[1])

# color by HIV status (no clear correlation here)
parcoord(eig$vectors[randord,1:10], col = dat$G1[randord])
plot(eig$vectors[,1], eig$vectors[,2], col = dat$G1)
plot(eig$vectors[,1], jitter(rep(0, dim(eig$vectors)[1])), col = dat$G1)

# color by KS status
parcoord(eig$vectors[randord,1:10], col = is.na(dat$F1[randord]) + 1)

plot(eig$vectors[,1], jitter(rep(0, dim(eig$vectors)[1])), col = is.na(dat$F1[randord]) + 1)

with(dat.sub, plot(match, match2, col = ifelse(cases, 'red', 'black')))

# color by pool from previous PCA
with(dat.sub, plot(match, match2, col = ifelse(is.na(pool), 4, pool), pch = ifelse(cases, 18, 21)))

# first PC accounts for a lot of the variation
(cumsum(eig$values) / sum(abs(eig$values)))[1:10]

# subpools
dat$subpool <- character(dim(dat)[1])

dat$subpool <- ifelse(dat$ID %in% subset(dat.sub, pool == 1)$ID, 'a', NA)
dat$subpool <- ifelse(dat$ID %in% subset(dat.sub, pool == 2)$ID,
                      ifelse(dat$match < -0.029095, 'b', 'c'), dat$subpool)
dat$subpool <- ifelse(dat$ID %in% subset(dat.sub, pool == 3)$ID, 'd', dat$subpool)

# look at this again with the subpool information added in
with(dat, plot(match, jitter(rep(0, dim(dat)[1])), col = ifelse(subpool == 'c', 4, pool),
               pch = ifelse(cases, 18, 21), ylim = c(-.1, .1),
               xlim = c(min(match, na.rm = TRUE), max(match[pool == 3], na.rm = TRUE))))

#################################################
# Evaluate which variables are most influential #
#################################################

vars <- names(dat)[!names(dat) %in% c('ID', 'F1', 'G1', 'match', 'cases', 'pool', 'HIV', 'KS', 'subpool')]

pvals <- rep(NA, length(vars))
names(pvals) <- vars

for(v in vars)
{
    if(all(is.na(dat[[v]][!is.na(dat$subpool)])))
    {
        pvals[v] <- 1
        next
    }

    null <- eval(parse(text = paste("lm(", v, " ~ 1, data = subset(dat, !is.na(subpool)))", sep = '')))
    alt <- try(eval(parse(text = paste("lm(", v, " ~ match, data = subset(dat, !is.na(subpool)))", sep = ''))))

    if(class(alt) == 'try-error') # too many missing data! -- analyze if missingness is a factor
    {
        null <- eval(parse(text = paste("glm(is.na(", v, ") ~ 1, data = subset(dat, !is.na(subpool)),",
                           "family = binomial)")))
        alt <- eval(parse(text = paste("glm(is.na(", v, ") ~ match, data = subset(dat, !is.na(subpool)),",
                          "family = binomial)", sep = '')))

        model <- anova(null, alt, test = 'Chisq')

        pvals[v] <- model[["Pr(>Chi)"]][2]
    }else{
        model <- anova(null, alt)

        pvals[v] <- model[["Pr(>F)"]][2]
    }
}

high <- unique(sapply(strsplit(names(pvals)[which(pvals < 0.01)], '_', fixed = TRUE), function(x) x[1]))
moderate <- unique(sapply(strsplit(names(pvals)[which(pvals > 0.01 & pvals < 0.05)], '_', fixed = TRUE),
                          function(x) x[1]))
trend <- unique(sapply(strsplit(names(pvals)[which(pvals > 0.05 & pvals < 0.1)], '_', fixed = TRUE),
                       function(x) x[1]))

# looks OK
translation <- read.xls('data/VariableNames.xlsx', stringsAsFactors = FALSE)
translation[translation$NewName %in% high, 2:3]
translation[translation$NewName %in% moderate, 2:3]
translation[translation$NewName %in% trend, 2:3]

#######################
# Pick final controls #
#######################

dat$pick <- dat$cases & !is.na(dat$subpool)

tmp <- which(dat$pick)

for(i in tmp)
{
    # all where pool is defined are male
    hiv <- dat$HIV[i]
    mtch <- dat$match[i]
    pl <- dat$pool[i]
    for(j in 1:4) # this is the ratio of cases:controls we want
    {
        pick.new <- as.character(with(subset(dat, !pick & HIV == hiv & !KS & pool == pl),
                                      ID[which.min(abs(mtch - match))]))

        dat$pick[dat$ID == pick.new] <- TRUE
    }
}

# we may not quite maximize the number of controls to pick from if we run out of controls in a particular pool
# if so, sample some more in other pools

total.sample <- 525

how.many <- total.sample - sum(dat$pick)

if(how.many > 0)
{
    set.seed(298347)

    pick.new <- with(dat, which(!pick & !KS & !is.na(pool)))

    dat$pick[sample(pick.new, how.many)] <- TRUE
}

dat.sub <- subset(dat, pick, select = c('ID', 'pool', 'HIV', 'KS'))

write.table(dat.sub, file = 'pools.txt', row.names = FALSE)

write.csv(subset(dat.sub, select = c('ID', 'HIV', 'KS')), file = 'matched.csv', row.names = FALSE)
