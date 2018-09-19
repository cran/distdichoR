#'@export
print.distdicho <- function(x, ...) {
    # t.test
    if (!is.null(x$ttest)) {
        cat("======================================================\n")
        cat("===              t-Test                            ===\n")
        cat("======================================================\n")
        print(x$ttest)
    }


    # proportion test
    cat("======================================================\n")
    cat("===              Distributional method             ===\n")
    cat("======================================================\n")

    # methods 
    cat("Distributional estimates for the comparison of proportions", ifelse(x$arguments$tail ==
        "upper", "above", "below"), "the cut-point", x$arguments$cutpoint, "\n")

    # variances
    if (!is.null(x$arguments$ratio)) {
        if (x$arguments$ratio == 0) {
              if (x$arguments$correction) {
                   cat("Standard error computed with correction for large effects\n\n")
            } else {
                cat("Standard error computed with correction for unknown variance ratio\n\n")
            }} else {
            cat("Standard error computed under the hypothesis that the ratio of variances is equal to",
                x$arguments$ratio, "\n\n")
        }
    } else {
        cat("Standard error computed under the hypothesis that the ratio of variances is equal to 1 \n\n")
    }
    if (!is.null(x$arguments$alpha))
        cat("Alpha:", x$arguments$alpha, "\n\n")

    # table
    tab1 <- data.frame(Group = x$data.name, Obs = c(x$parameter["n1"], x$parameter["n2"]), Mean = c(x$parameter["m1"],
        x$parameter["m2"]), Std.Dev = c(x$parameter["s1"], x$parameter["s2"]), Dist.prop. = c(x$prop["proportion1"],
        x$prop["proportion2"]))

    tab2 <- data.frame(Stat = c("Diff. prop", "Risk ratio", "Odds ratio"), Estimate = c(x$dist.estimates["difference in proportions"],
        x$dist.estimates["risk ratio (distributional estimate)"], x$dist.estimates["odds ratio (distributional estimate)"]),
        Std.Err = c(x$se["standard error 'propdiff'"], x$se["standard error 'rr'"], x$se["standard error 'or'"]),
        CI.lower = c(x$ci["lower limit 'propdiff'"], x$ci["lower limit 'rr'"], x$ci["lower limit 'or'"]),
        CI.upper = c(x$ci["upper limit 'propdiff'"], x$ci["upper limit 'rr'"], x$ci["upper limit 'or'"]))
   
     # print table
    print(tab1, row.names = F)
    cat("\n------------------------------------------------------\n")
    print(tab2, row.names = F)
    cat("\n------------------------------------------------------\n")

    # confidence intervals
    cat("*", x$arguments$conf.level * 100, "percent confidence interval\n")

    # bootstrap
    if (!is.null(x$arguments$bootci)) {
        cat(ifelse(x$arguments$bootci == TRUE, "* confidence interval calculated using bootstrap standard error",
            "* confidence interval calculated using distributional standard error"), "\n\n")
    }
    cat("------------------------------------------------------")
}

