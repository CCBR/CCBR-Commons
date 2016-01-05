# reconstructAgeMatching.R
# Randall Johnson
# Leidos Biomedical Research, Inc


#library(lpSolve)
library(lpSolveAPI)

##### Read in data #####

source('readData.R')
rownames(dat) <- dat$ID

##### Set up Integer programming #####

# matching groups
# K=KS+, k=KS-
# F=female, M=male
# H=HIV+, h=hiv-
KF  <- with(dat, which( case &  female))
KM  <- with(dat, which( case & !female))

KFH <- with(dat, which( case &  female &  hiv))
KMH <- with(dat, which( case & !female &  hiv))
KFh <- with(dat, which( case &  female & !hiv))
KMh <- with(dat, which( case & !female & !hiv))

kFH <- with(dat, which(!case &  female &  hiv))
kMH <- with(dat, which(!case & !female &  hiv))
kFh <- with(dat, which(!case &  female & !hiv))
kMh <- with(dat, which(!case & !female & !hiv))

# possible matches
permute <- function(case, cont)
    return(cbind(rep(case, each = length(cont)),
                 rep(cont, length(case))))

matches <- list()

# matching as described in the paper ("The study was designed so that four HIV-positive controls and two HIV-negative controls were matched to each KS case based on sex and 5-year age brackets.")
matches$paper <- rbind(permute(dat$ID[KF], dat$ID[kFH]),
                       permute(dat$ID[KF], dat$ID[kFh]),
                       permute(dat$ID[KM], dat$ID[kMH]),
                       permute(dat$ID[KM], dat$ID[kMh]))

# matching including HIV status as a constraint
matches$onHIV <- rbind(permute(dat$ID[KFH], dat$ID[kFH]),
                       permute(dat$ID[KMH], dat$ID[kMH]),
                       permute(dat$ID[KFh], dat$ID[kFh]),
                       permute(dat$ID[KMh], dat$ID[kMh]))


# weighting vector to maximize
w <- list()

for(i in 1:2)
{
    w[[i]] <- rep(0, dim(matches[[i]])[1])
    names(w[[i]]) <- paste(matches[[i]][,1], matches[[i]][,2], sep = ':')

    # calculate difference in date of birth (in years)
    for(j in 1:dim(matches[[i]])[1])
    {
        w[[i]][j] <- with(dat, abs(dob[ID == matches[[i]][j,1]] - dob[ID == matches[[i]][j,2]]))^-1 # invert and maximize

        # if they are within 5 years of eachother, that is good enough (this takes care of shared birthdays that score Inf)
        if(w[[i]][j] > 1)
            w[[i]][j] <- 1
    }
}
names(w) <- c('paper', 'onHIV')


# constraints for matching as described by paper
rnames <- list()


paper <- make.lp(nrow = length(unique(matches$paper[,1])) + # each case matched to up to 4 HIV+ controls
                        length(unique(matches$paper[,1])) + # each case matched to up to 2 HIV- controls
                        length(unique(matches$paper[,2])),  # each control matched to only 1 case or fewer
                 ncol = dim(matches$paper)[1])

rnames$paper <- c(rep(unique(matches$paper[,1]), 2), unique(matches$paper[,2]))

for(j in 1:dim(paper)[2])
{
    case <- which(rnames$paper == matches$paper[j,1])
    cont <- which(rnames$paper == matches$paper[j,2])

    contID <- rnames$paper[cont]
    contHasHIV <- with(dat, hiv[ID == contID])

    # first set of rows is for HIV+ matches, second set is for HIV- matches (see definition of paper)
    set.mat(paper, case[ifelse(contHasHIV, 1, 2)], j, 1)

    # each control matched to only 1 case or fewer
    set.mat(paper, cont, j, 1)
}

tmp <- length(get.constr.value(paper))
set.constr.type(paper, types = c(rep('<=', 2*length(unique(matches$paper[,1]))),
                                 rep('=', length(unique(matches$paper[,2])))),
                1:tmp)
set.constr.value(paper, rhs = c(rep(4, length(unique(matches$paper[,1]))),  # each case matched by up to 4 HIV+ controls
                                rep(2, length(unique(matches$paper[,1]))),  # each case matched by up to 2 HIV- controls
                                rep(1, length(unique(matches$paper[,2])))), # each control matched to exactly 1 case
                 constraints = 1:tmp)

set.objfn(paper, -w$paper, 1:length(w$paper))

# constraints for additionally matching by HIV status
onHIV <- make.lp(nrow = length(unique(matches$onHIV[,1])) + # each case matched to up to X controls
                        length(unique(matches$onHIV[,2])),  # each control matched to exactly 1 case
                 ncol = dim(matches$onHIV)[1])

rnames$onHIV <- c(unique(matches$onHIV[,1]), unique(matches$onHIV[,2]))

for(j in 1:dim(onHIV)[2])
{
    case <- which(rnames$onHIV == matches$onHIV[j,1])
    cont <- which(rnames$onHIV == matches$onHIV[j,2])

    # each case matched to up to 5 controls
    set.mat(onHIV, case, j, 1)

    # each control matched to only 1 case or fewer
    set.mat(onHIV, cont, j, 1)
}

tmp <- length(get.constr.value(onHIV))
set.constr.type(onHIV, types = c(rep('<=', length(unique(matches$onHIV[,1]))), # each case matched by up to X controls
                                 rep('=', length(unique(matches$onHIV[,2])))), # each control matched to exactly 1 case
                1:tmp)
set.constr.value(onHIV, rhs = ifelse(rnames$onHIV %in% dat$ID[c(kFH, kMH, kMh)], 1,  # controls
                              ifelse(rnames$onHIV %in% dat$ID[c(KFH, KMH)], 5, 14)), # KMh cases are disporportionately matched
                 constraints = 1:tmp)

set.objfn(onHIV, -w$onHIV, 1:length(w$onHIV))

##### Run Matching #####

set.seed(298347)

solve(paper)
tmp <- get.variables(paper) * w$paper
hist(tmp[tmp != 0 & tmp < 1])

picks <- data.frame(case = matches$paper[as.logical(get.variables(paper)),1],
                    cont = matches$paper[as.logical(get.variables(paper)),2],
                    stringsAsFactors = FALSE)
picks$ageDiff <- abs(dat[picks$case,"dob"] - dat[picks$cont, "dob"])
picks$caseHIV <- dat[picks$case, 'hiv']
picks$contHIV <- dat[picks$cont, 'hiv']
picks$female <- dat[picks$case, 'female']
write.csv(picks, file = 'matched_on_gender_age.csv')

solve(onHIV)

tmp <- get.variables(onHIV) * w$onHIV
hist(tmp[tmp != 0])

picks <- data.frame(case = matches$onHIV[as.logical(get.variables(onHIV)),1],
                    cont = matches$onHIV[as.logical(get.variables(onHIV)),2],
                    stringsAsFactors = FALSE)
picks$ageDiff <- abs(dat[picks$case,"dob"] - dat[picks$cont, "dob"])
picks$hiv <- dat[picks$case, 'hiv']
picks$female <- dat[picks$case, 'female']
write.csv(picks, file = 'matched_on_gender_hiv_age.csv')
