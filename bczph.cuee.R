bczph.cuee <-
    function(formula,
             cumA = NA,
             cumU = NA,
             prevBeta = NA,
             prevInvVar = NA,
             prevCumScore = NA,
             prevCumProd = NA,
             newdata,
             transform = "km")
    {
        ## local fit to get betahat
        newFit <-
            do.call(coxph, list(formula = formula, data = newdata))
        ########################################################################
        ## scaled schoenfeld residual of the local fit, not very useful
        ## just need event times, betahat and its inverse variance
        ## schoenfeld residual
        varnames <- names(newFit$coefficients)
        nvar <- length(varnames)
        ssresid <- resid(newFit, "schoenfeld")
        sresid <- ssresid %*% newFit$var
        ndead <- length(ssresid) / nvar
        times <- as.numeric(rownames(ssresid))
        if (is.character(transform)) {
            tname <- transform
            ttimes <- switch(
                transform,
                'identity' = times,
                'rank'    = rank(times),
                'log'     = log(times),
                'km' = {
                    temp <- survfitKM(factor(rep(1, nrow(
                        newFit$y
                    ))), newFit$y,
                    se.fit = FALSE)
                    # A nuisance to do left cont KM
                    t1 <- temp$surv[temp$n.event > 0]
                    t2 <- temp$n.event[temp$n.event > 0]
                    km <- rep(c(1, t1), c(t2, 0))
                    if (is.null(attr(sresid, 'strata')))
                        1 - km
                    else
                        (1 - km[sort.list(sort.list(times))])
                },
                stop("Unrecognized transform")
            )
        } else {
            tname <- deparse(substitute(transform))
            if (length(tname) > 1)
                tname <- 'user'
            ttimes <- transform(times)
        }
        xx <- ttimes - mean(ttimes)
        ########################################################################
        betaHat <- newFit$coefficients
        invVarHat <- solve(newFit$var)
        ## if this is the first block, initialize all components to zero
        if (sum(is.na(prevBeta))) {
            cumA <- matrix(0, nvar, nvar)
            cumU <- matrix(0, nvar, 1)
            prevBeta <- matrix(0, nvar, 1)
            prevInvVar <- matrix(0, nvar, nvar)
            prevCumScore <- matrix(0, nvar, 1)
            prevCumProd <- matrix(0, nvar, 1)
        }
        ## get betaCheck, equation 21 in technometrics paper
        betaCheck <- solve(prevInvVar + invVarHat,
                           prevCumProd + invVarHat %*% betaHat)

        ## get the score vector evaluated at betaCheck
        ctrl <- coxph.control(iter.max = 0)
        checkFit <- do.call(coxph, list(formula = formula, data = newdata,
                                init = betaCheck, control = ctrl))
        score <- prevCumScore + colSums(coxph.detail(checkFit)$score)
        invVarCheck <- solve(checkFit$var)

        ## betaTilde
        betaTilde <- solve(prevInvVar + invVarCheck,
                           prevCumProd + invVarCheck %*% betaCheck + score)

        ## tildeFit to get the A and U parts of the cumulative test statistic
        tildeFit <- do.call(coxph, list(formula = formula, data = newdata,
                                        init = betaTilde, control = ctrl))
        V <- tildeFit$var
        cA <- sum(xx ^ 2) * V / ndead
        nssresid <- resid(tildeFit, "schoenfeld")
        ## schoenfeld residual
        nsresid <- nssresid %*% V
        cU <- c(xx %*% nsresid)
        cumA <- cumA + cA
        cumU <- cumU + cU
        cumZ <- c(t(cumU) %*% solve(cumA, cumU))
        output <- list(
            cumA = cumA,
            cumU = cumU,
            beta = betaTilde,
            invVar = prevInvVar + invVarCheck,
            cumScore = score,
            cumProd = prevCumProd + invVarCheck %*% betaCheck,
            Stat = cumZ,
            p = 1 - pchisq(cumZ, ncol(sresid))
        )
        class(output) <- "bczph.cuee"
        return(output)
    }
