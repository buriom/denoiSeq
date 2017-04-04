
context("test the summation term in the likelihood")

N <- 2 # initial copy number
counts <- c(3,5) # observed final copy number
p <- 0.2 # amplicationefficiency
f <- 0.5 # downsampling rate

propotns <- c(1,1)

test_that("calculation of the summation term in the likelihood is correct", {
  expect_equal(sum(mapply(summation_term, counts = counts, propotns = propotns,
                          MoreArgs = list(N = N, p = p, f = f))), log(2100/16))
})
