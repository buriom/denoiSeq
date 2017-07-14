
<!-- README.md is generated from README.Rmd. Please edit that file -->
An R package for Differential Expression Analysis. Given count data from two experimental conditions, denoiSeq helps one determine which transcripts are differentially expressed across the two conditions using Bayesian inference of the parameters of a bottom-up model for PCR amplification developed in "Chromatin conformation governs T cell receptor J beta gene segment usage", by Ndifon et al.

To use the package, one needs to create a `readsData` object and invoke the `denoiseq` function on it. The results are obtained from the return value of denoiseq using the `results` function which then computes the test statistic used in differential analysis.

    RD <- new("readsData", counts = ERCC)  #creating the readsData object

    steps <- 3000  #steps for MCMC

    BI <- denoiseq(RD, steps)  #invoking denoiseq on the readsData object

    rez <- results(BI,steps)  #computing the test statistic

This package can be istalled from CRAN using install.packages("denoiSeq") or from github using devtools::install\_github("buriom/denoiSeq").
