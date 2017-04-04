# colnames(N_samples2) <- rownames(taqman)

chainplot <- function(parm, steps, cond = 1) {
  if (parm == "p") {
    p_samples <- t(lapply(1:steps, getp, rezult = rezult))
    p <- head(p_samples[1, ], steps)
    plot(plot(p, type = "l", main = "History plot of p."))
  } else if (parm == "f") {
    p_samples <- t(lapply(1:steps, getf, rezult = rezult))
    p <- head(p_samples[1, ], steps)
    plot(plot(p, type = "l", main = "History plot of f."))
  } else if (parm %in% rownames(counts_A) & cond == 1) {
    N_Asamples <- t(lapply(1:steps, getN_A, rezult = rezult))
    N <- N_Asamples[, "parm"]
    plot(head(N), type = "l", main = paste("History plot of", "parm", "in condition 1"))
  } else if (parm %in% rownames(counts_A) & cond == 2) {
    N_Bsamples <- t(lapply(1:steps, getN_B, rezult = rezult))
    N <- N_Bsamples[, "parm"]
    plot(head(N), type = "l", main = paste("History plot of", "parm", "in condition 2"))
  } else {
    print("Unknown parameter")
  }
}

aucplot <- function(parm, steps, cond = 1) {
  if (parm == "p") {
    p_samples <- t(lapply(1:steps, getp, rezult = rezult))
    p <- head(p_samples[1, ], steps)
    acf(p, lag.max = 100, type = c("correlation"), plot = T)
  } else if (parm == "f") {
    f_samples <- t(lapply(1:steps, getf, rezult = rezult))
    f <- head(f_samples[1, ], steps)
    acf(f, lag.max = 100, type = c("correlation"), plot = T)
  } else if (parm %in% rownames(counts_A) & cond == 1) {
    N_Asamples <- t(lapply(1:steps, getN_A, rezult = rezult))
    N <- N_Asamples[, "parm"]
    acf(N, lag.max = 100, type = c("correlation"), plot = T)
  } else if (parm %in% rownames(counts_A) & cond == 2) {
    N_Bsamples <- t(lapply(1:steps, getN_B, rezult = rezult))
    N <- N_Bsamples[, "parm"]
    acf(N, lag.max = 100, type = c("correlation"), plot = T)
  } else {
    print("Unknown parameter")
  }
}
