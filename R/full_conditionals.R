# _____________________________________________________________________________

# Implementation of the full conditionals for parameters N, p and f.
# _____________________________________________________________________________
# To calculate the summation in the probability distribution
summation_term <- function(N, f, p, propotns, x) {
    j <- 0:min(c(x, N))
    ind <- lgamma(N + x - j + 1) - lgamma(N - j + 1) - lgamma(j + 1) - lgamma(x - j + 1) + j * (log(p + (1 - p) * f * propotns) -
        log(1 - p) - log(1 - f * propotns)) - log(N + 0.005 + x - j)
    u <- max(ind)
    z <- sum(exp(ind - u))
    return(u + log(z))
}

# To calculate the full conditional for N.
fullcond_N <- function(N, f, p, propotns, x) {
    # summation
    b <- mapply(summation_term, x = x, propotns = propotns, MoreArgs = list(N = N, p = p, f = f))
    # multiplying by the prior
    fc <- sum(log(N) + N * (log(p) + log(1 - propotns * f) - log(p + (1 - p) * propotns * f))) + sum(b)
    return(fc)
}

# To calculate the full conditional for p.
fullcond_p <- function(N, f, p, propotns, x) {
    # summation
    b <- mapply(summation_term, x = x, propotns = propotns, N = N, MoreArgs = list(p = p, f = f))
    # multiplying by the prior alpha0 <- 1 beta0 <- 1 (alpha0 - 1) * log(p) + (beta0 - 1) *log(1 - p) +
    fc <- sum(x * log(1 - p) + N * log(p) - (N + x) * (log(p + (1 - p) * propotns * f))) + sum(b)
    return(fc)
}

# To calculate the full conditional for f.
fullcond_f <- function(N, f, p, propotns, x) {
    # summation
    b <- mapply(summation_term, x = x, propotns = propotns, N = N, MoreArgs = list(p = p, f = f))
    # multiplying by the prior theta0 <- 1 gamma0 <- 1 (theta0 - 1) * log(f) + (gamma0 - 1) *log(1 - f) +
    fc <- sum(x * log(f) + N * log(1 - propotns * f) - (N + x) * (log(p + (1 - p) * propotns * f))) + sum(b)
    return(fc)
}
