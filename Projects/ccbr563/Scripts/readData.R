# readData.R
# Read and format the clinical data for downstream analysis
# Randall Johnson
# Leidos Biomedical Research, Inc


library(Epi)
library(magrittr)

# Read in data
full <- read.csv('../Data/ksbase.csv', na.strings = c('', 'Legitimate Skip', '-1', '88', '99'), stringsAsFactors = FALSE)
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

# a few extra variables I will want
dat$case <- with(dat, substr(ID, 1, 1) == 'K')

dat$female <- dat$A3 == 2

months <- 1:12
names(months) <- c('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC')

dat$dob <-  dat$A4
            ## with(dat, paste(substr(A4, 6, 9),  # date of birth
            ##                 months[substr(A4, 3, 5)],
            ##                 substr(A4, 1, 2), sep = '-')) %>%
            ## cal.yr()

dat$hiv <- dat$G1 == 1

rownames(dat) <- dat$ID
