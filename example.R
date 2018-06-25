## this is an example how online updating cumulative statistic with CUEE
## can be obtained

## simulated_survdata is generated using the genData function provided

source("./bczph.cuee.R")
survDat <- read.csv("./simulated_survdata.csv")

dims <- 3 
cuma <- previnvvar <- matrix(0, dims, dims)
prevbeta <- cumu <- prevcumscore <- prevcumprod <- matrix(0, dims, 1)

## block size 5000
survDat$block <- ceiling((1:nrow(survDat)) / 5000)

updateStats <- c()
pb <- txtProgressBar(style = 3, min = 0, max = max(survDat$block))

for (i in 1:max(survDat$block)){
    blockDat <- survDat[survDat$block == i, ]
    
    bczphFit <- bczph.cuee(Surv(survtime, status) ~ age + Sex + Black,
                        cumA = cuma, cumU = cumu, prevBeta = prevbeta,
                        prevInvVar = previnvvar, prevCumScore = prevcumscore,
                        prevCumProd = prevcumprod, newdata = blockDat,
                        transform = "km")
    updateStats[i] <- bczphFit$Stat
    
    cuma <- bczphFit$cumA
    cumu <- bczphFit$cumU
    prevbeta <- bczphFit$beta
    previnvvar <- bczphFit$invVar
    prevcumscore <- bczphFit$cumScore
    prevcumprod <- bczphFit$cumProd
    
    Sys.sleep(0.1); setTxtProgressBar(pb, i)
}

plot(updateStats, type = "l", ylim = c(0, 10))
abline(a = qchisq(0.95, 3), b = 0, col = "red", lty = 2)
