genData <- function(lambda, beta, n, tmax, epsilon) {
    beta <- do.call(matrix, args = list(beta, length(beta)))
    rtmax <- runif(n, 0, tmax)
    vepsilon <- rbinom(n, 1, epsilon)
    rtmax <- tmax * vepsilon + (1 - vepsilon) * rtmax
    rtmax <- pmin(rtmax, tmax)
    p <- length(beta)
    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, 0.5)
    x3 <- rbinom(n, 1, 0.1)
    x <- cbind(x1, x2, x3)
    hazard <- lambda * exp(x %*% beta)
    tfail <- rexp(n, rate = 1)
    tfail <- tfail / hazard
    status <- (tfail <= rtmax) + 0
    time <- tfail * status + rtmax * (1 - status)
    return(data.frame(
        survtime = time,
        status = status,
        age = x1,
        Sex = x2,
        Black = x3
    ))
}



genDataFrailty <- function(lambda, beta, n,  tmax, epsilon, sd) {
        # Use epsilon<tmax> + (1-epsilon)U(0, tmax) as random censoring time
        rtmax <- runif(n, 0, tmax)
        vepsilon <- rbinom(n, 1, epsilon)
        rtmax <- tmax * vepsilon + (1 - vepsilon) * rtmax
        # rtmax <- runif(n, ltuning * tmax, utuning * tmax)
        rtmax <- pmin(rtmax, tmax)
        p <- length(beta)
        x1 <- rnorm(n)
        x2 <- rbinom(n, 1, 0.5)
        x3 <- rbinom(n, 1, 0.1)
        x <- cbind(x1, x2, x3)
        # frailty
        w <- rnorm(n, 0, sd)
        hazard <- lambda * exp(x %*% beta + w)
        # attention
        tfail <- rexp(n, rate = 1)
        tfail <- tfail / hazard
        # why normalize using hazard?
        status <- (tfail <= rtmax) + 0
        time <- tfail * status + rtmax * (1 - status)
        return(data.frame(
            survtime = time,
            status = status,
            age = x1,
            Sex = x2,
            Black = x3
        ))
    }
