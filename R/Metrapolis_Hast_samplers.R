# _____________________________________________________________________________

# Implementation of the Metrapolis-Hastings algorithm for sampling from each of
# the full conditionals
# _____________________________________________________________________________

# Normal Metrapolis-Hastings samplers


#' @importFrom stats runif median sd dbeta rbeta
# To implement Metrapolis-Hastings for N
MH_N <- function(N_b, f, p, propotns, x, step_size) {
    k <- runif(2, 0, 1)
    u1 <- k[1]
    u2 <- k[2]
    N_c <- exp(log(N_b) + step_size * 2 * (u1 - 1 / 2))
    r <- exp(fullcond_N(N_c, f, p, propotns, x) - fullcond_N(N_b, f, p, propotns, x))
    if (is.nan(r) || u2 > r) {
        N_a <- N_b
    } else {
        N_a <- N_c
    }
    return(N_a)
}

# To implement Metrapolis Hastings for p
MH_p <- function(N, f, p_b, propotns, x, step_size) {
    u <- runif(1, 0, 1)
    p_c <- rbeta(1, step_size * p_b, step_size * (1 - p_b))
    r <- exp(fullcond_p(N, f, p_c, propotns, x) - fullcond_p(N, f, p_b, propotns, x) + dbeta(p_b, step_size * p_c, step_size *
        (1 - p_c), log = TRUE) - dbeta(p_c, step_size * p_b, step_size * (1 - p_b), log = TRUE))
    if (u < r) {
        p_a <- p_c
    } else {
        p_a <- p_b
    }
    return(p_a)
}

# To implement Metrapolis-Hastings for f
MH_f <- function(N, f_b, p, propotns, x, step_size) {
    u <- runif(1, 0, 1)
    f_c <- rbeta(1, step_size * f_b, step_size * (1 - f_b))
    r <- exp(fullcond_f(N, f_c, p, propotns, x) - fullcond_f(N, f_b, p, propotns, x) + dbeta(f_b, step_size * f_c, step_size *
        (1 - f_c), log = TRUE) - dbeta(f_c, step_size * f_b, step_size * (1 - f_b), log = TRUE))
    if (u < r) {
        f_a <- f_c
    } else {
        f_a <- f_b
    }
    return(f_a)
}

# ____________________________________________________________________________

# Metrapolis-Hastings modified to tune the step size
# _____________________________________________________________________________
#
#  To implement Metrapolis-Hastings for N
MH_N2 <- function(N_b, f, p, propotns, x, step_size) {
    k <- runif(2, 0, 1)
    u1 <- k[1]
    u2 <- k[2]
    N_c <- exp(log(N_b) + step_size * 2 * (u1 - 1 / 2))
    r <- exp(fullcond_N(N_c, f, p, propotns, x) - fullcond_N(N_b, f, p, propotns, x))
    # tuning to obtain around 42% acceptance rate
    if (is.nan(r) || u2 > r) {
        N_a <- N_b
        step_size <- step_size / 1.007
    } else {
        N_a <- N_c
        step_size <- step_size * 1.01
    }
    return(c(N_a, step_size))
}


# To implement Metrapolis-Hastings for p
MH_p2 <- function(N, f, p_b, propotns, x, step_size) {
    u <- runif(1, 0, 1)
    p_c <- rbeta(1, step_size * p_b, step_size * (1 - p_b))
    r <- exp(fullcond_p(N, f, p_c, propotns, x) - fullcond_p(N, f, p_b, propotns, x) + dbeta(p_b, step_size * p_c, step_size *
        (1 - p_c), log = TRUE) - dbeta(p_c, step_size * p_b, step_size * (1 - p_b), log = TRUE))
    if (u < r) {
        p_a <- p_c
        step_size <- step_size / 1.01
    } else {
        p_a <- p_b
        step_size <- step_size * 1.007
    }
    return(c(p_a, step_size))
}

# To implement Metrapolis-Hastings for f
MH_f2 <- function(N, f_b, p, propotns, x, step_size) {
    u <- runif(1, 0, 1)
    f_c <- rbeta(1, step_size * f_b, step_size * (1 - f_b))
    r <- exp(fullcond_f(N, f_c, p, propotns, x) - fullcond_f(N, f_b, p, propotns, x) + dbeta(f_b, step_size * f_c, step_size *
        (1 - f_c), log = TRUE) - dbeta(f_c, step_size * f_b, step_size * (1 - f_b), log = TRUE))
    if (u < r) {
        f_a <- f_c
        step_size <- step_size / 1.01
    } else {
        f_a <- f_b
        step_size <- step_size * 1.007
    }
    return(c(f_a, step_size))
}
