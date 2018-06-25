bczph.cee <-
    function(formula,
             cumA = 0,
             cumU = 0,
             prevBeta = 0,
             prevInvVar = 0,
             newdata,
             transform = "km")
    {
        newfit <- do.call(coxph, list(formula = formula, data = newdata))
        ## scaled schoenfeld residual
        ssresid <- resid(newfit, "schoenfeld")
        ## schoenfeld residual
        sresid <- ssresid %*% newfit$var
        varnames <- names(newfit$coefficients)
        nvar <- length(varnames)
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
                        newfit$y
                    ))), newfit$y,
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
        cumInvVar <- prevInvVar + solve(newfit$var)
        cumVar <- solve(cumInvVar)
        cumBeta <- cumVar %*% (prevInvVar %*% prevBeta +
                                   solve(newfit$var, %*% newfit$coefficients))
        cumFit <-
            do.call(
                coxph,
                list(
                    formula = formula,
                    data = newdata,
                    init = cumBeta,
                    control = coxph.control(iter.max = 0)
                )
            )
        V <- cumFit$var
        cA <- sum(xx ^ 2) * V / ndead
        nssresid <- resid(cumFit, "schoenfeld")
        ## schoenfeld residual
        nsresid <- nssresid %*% V
        cU <- c(xx %*% nsresid)
        cumA <- cumA + cA
        cumU <- cumU + cU
        cumZ <- c(t(cumU) %*% solve(cumA,%*% cumU))
        output <- list(
            cumA = cumA,
            cumU = cumU,
            beta = cumBeta,
            invVar = cumInvVar,
            Stat = cumZ,
            p = 1 - pchisq(cumZ, ncol(sresid))
        )
        class(output) <- "bczph.cee"
        return(output)
    }
