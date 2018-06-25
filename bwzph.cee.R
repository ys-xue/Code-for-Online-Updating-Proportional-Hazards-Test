bwzph.cee <-
    function(formula,
             wA,
             wU,
             newdata,
             transform = "km",
             width = 5,
             wbeta,
             wvar) {
        newfit <- do.call(coxph, list(formula = formula, data = newdata))
        wbeta[[length(wbeta) + 1]] <- newfit$coefficients
        wvar[[length(wvar) + 1]] <- newfit$var
        if (length(wbeta) > width) {
            wbeta <- tail(wbeta, width)
            wvar <- tail(wvar, width)
        }
        tosum <-
            lapply(1:length(wbeta), function(x)
                solve(wvar[[x]], wbeta[[x]]))
        invvar <- lapply(1:length(wvar), function(x)
            solve(wvar[[x]]))
        temp <- Reduce("+", tosum)
        cuminvVar <- Reduce("+", invvar)
        cumBeta <- solve(cuminvVar, temp)
        cumfit <- do.call(
            coxph,
            list(
                formula = formula,
                data = newdata,
                init = cumBeta,
                control = coxph.control(iter.max = 0)
            )
        )
        # scaled Schoenfeld residuals
        sresid <- resid(cumfit, "schoenfeld")
        nsresid <- sresid %*% cumfit$var
        varnames <- names(cumfit$coefficients)
        nvar <- length(varnames)
        ndead <- length(sresid) / nvar
        if (nvar == 1) {
            times <- as.numeric(names(sresid))
        } else {
            times <- as.numeric(rownames(sresid))
        }
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
        }
        else {
            tname <- deparse(substitute(transform))
            if (length(tname) > 1)
                tname <- 'user'
            ttimes <- transform(times)
        }
        ttimes <- ttimes - mean(ttimes)

        V <- cumfit$var
        cA <- sum(ttimes ^ 2) * V / ndead
        cU <- c(ttimes %*% nsresid)
        wA[[length(wA) + 1]] <- cA
        wU[[length(wU) + 1]] <- cU
        if (length(wA) > width) {
            wA <- tail(wA, width)
            wU <- tail(wU, width)
        }
        cumA <- Reduce("+", wA)
        #print(cumA)
        cumU <- Reduce("+", wU)
        cumZ <- c(t(cumU) %*% solve(cumA, cumU))
        output <- list(
            WindowA = wA,
            WindowU = wU,
            WindowBeta = wbeta,
            WindowVar = wvar,
            cumBeta = cumBeta,
            Stat = cumZ,
            p = 1 - pchisq(cumZ, ncol(sresid))
        )
        class(output) <- "bwzph.cee"
        return(output)
    }
