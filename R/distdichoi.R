#' nomal data (immdediate form, allowing unequal variances)
#'
#' Immediate form of the distributional method for dichotomising normal data
#' allowing for assumptions of unequal variances
#' (based on Sauzet et al. 2014 and Peacock et al. 2012).
#'
#'
#' distdichoi takes no data, but the number of observations as well as the mean and standard deviations of both groups.
#' It first returns the results of a two-group unpaired t-test (allowing for unequal variances in the unequal variance cases).
#' Followed by the distributional estimates and their standard errors (see Sauzet et al. 2014 and Peacock et al. 2012)
#' for a difference in proportions, risk ratio and odds ratio. It also provides the distributional confidence intervals for the statistics estimated
#' (this assumes an asymptotic normal distribution of estimates and might not be valid for small sample sizes (see Sauzet et al. 2014 for details)).
#' Estimates are calculated using either assumption of equal variances in both groups (default R = 1) or assumption of
#' unequal variance ratio (R != 1 & R !=0 for known variance ratio and R=0 for correction for unknown variance ratio).
#'
#'
#'
#'
#'@param n1 A number specifying the number of observations in the exposed group.
#'@param m1 A number specifying the mean of the exposed group.
#'@param s1 A number specifying the standard deviation of the exposed group.
#'@param n2 A number specifying the number of observations in the unexposed (reference) group
#'@param m2 A number specifying the mean of the unexposed (reference) group.
#'@param s2 A number specifying the standard deviation of the unexposed (reference) group.
#'@param cp A numeric value specifying the cut point under or over which the distributional proportions are computed.
#'@param tail A character string specifying the tail of the distribution in which the proportions are computed.
#'            Must be either 'lower' (default) or 'upper'.
#'@param R A numeric value indicating the true ratio of variances (R = Var(group1)/Var(group2).
#'         A value of 0 specifies that the true ratio of variances is unknown.
#'@param conf.level Confidence level of the interval.
#'
#' @return A list with class 'distdicho' containing the following components:
#' \item{data.name}{The names of the data.}
#' \item{arguments}{A list with the specified arguments.}
#' \item{parameter}{The mean, standard error and number of observations for both groups.}
#' \item{prop}{The estimated proportions below / above the cut point for both groups.}
#' \item{dist.estimates}{The difference in proportions, risk ratio and odds ratio of the groups.}
#' \item{se}{The estimated standard error of the difference in proportions, the risk ratio and the odds ratio.}
#' \item{ci}{The confidence intervals of the difference in proportions, the risk ratio and the odds ratio.}
#' \item{method}{A character string indicating the used method.}
#' \item{ttest}{A list containing the results of a t-test.}
#'
#'
#'@seealso \code{\link[distdichoR]{distdicho}}, \code{\link[distdichoR]{distdichogen}}, \code{\link[distdichoR]{distdichoigen}}, \code{\link[distdichoR]{regdistdicho}}
#'
#'@references
#' Peacock J.L., Sauzet O., Ewings S.M., Kerry S.M. Dichotomising continuous data while retaining statistical power using a distributional approach.  Statist. Med; 2012;26:3089-3103.
#' Sauzet, O., Peacock, J. L. Estimating dichotomised outcomes in two groups with unequal variances: a distributional approach.  Statist. Med; 2014 33 4547-4559 ;DOI: 10.1002/sim.6255.
#' Peacock, J.L., Bland, J.M., Anderson, H.R.: Preterm delivery: effects of socioeconomic factors, psychological stress, smoking, alcohol, and caffeine. BMJ 311(7004), 531-535 (1995).
#'
#'@examples
#'# Immediate form of distdicho
#' distdichoi(n1 = 494, m1 = 3267.4, s1 = 441.3,
#'            n2 = 983, m2 = 3452, s2 = 435.9,
#'            cp = 2500, tail = 'upper')
#'
#'## Proportions of low birth weight babies among smoking and non-smoking mothers
#'## (data from Peacock et al. 1995). Returns distributional estimates, standard 
#'## errors and distributional confidence intervals for differences in proportions,
#'## RR and OR of babies having a birth weight under 2500g (low birth weight LBW)
#'## for group smoker (mother smokes) over the odds of LBW in group non-smoker 
#'## (mother doesn't smoke)
#'# distdicho and distdichoi are returning the same results
#' bw_smoker <- bwsmoke$birthwt[bwsmoke$smoke == 'smoker']
#' bw_nonsmoker <- bwsmoke$birthwt[bwsmoke$smoke == 'non-smoker']
#' distdicho(x = bw_smoker, y = bw_nonsmoker, cp = 2500)
#' distdichoi(n1 = length(bw_smoker[!is.na(bw_smoker)]), 
#'            m1 = mean(bw_smoker, na.rm = TRUE), 
#'            s1 = sd(bw_smoker, na.rm = TRUE),
#'            n2 = length(bw_nonsmoker[!is.na(bw_smoker)]), 
#'            m2 = mean(bw_nonsmoker, na.rm = TRUE), 
#'            s2 = sd(bw_nonsmoker, na.rm = TRUE), 
#'            cp = 2500)
#'
#'@export
distdichoi <- function(n1, m1, s1, n2, m2, s2, cp = 0, tail = c("lower", "upper"), R = 1, conf.level = 0.95) {
    
    # 1 verify arguments
    tail <- match.arg(tail)
    
    if (length(n1) != 1 || !is.numeric(n1) || n1%%1 != 0 || n1 <= 0 || is.infinite(n1)) 
        stop("'n1' must be a single positive number")
    
    if (length(n2) != 1 || !is.numeric(n2) || n2%%1 != 0 || n2 <= 0 || is.infinite(n2)) 
        stop("'n2' must be a single positive number")
    
    if (length(m1) != 1 || !is.numeric(m1) || is.infinite(m1)) 
        stop("'m1' must be a single number")
    
    if (length(m2) != 1 || !is.numeric(m2) || is.infinite(m2)) 
        stop("'m2' must be a single number")
    
    if (length(s1) != 1 || !is.numeric(s1) || s1 <= 0 || is.infinite(s1)) 
        stop("'s1' must be a single positive number")
    
    if (length(s2) != 1 || !is.numeric(s2) || s2 <= 0 || is.infinite(s2)) 
        stop("'s2' must be a single positive number")
    
    if (length(cp) != 1 || !is.numeric(cp) || is.infinite(cp)) 
        stop("'cp' must be a single number")
    
    if (length(R) != 1 || !is.numeric(R) || R < 0 || is.infinite(R)) 
        stop("'R' must be a single non-negative number")
    
    if (length(conf.level) != 1 || !is.numeric(conf.level) || conf.level < 0 || conf.level > 1 || is.infinite(conf.level)) 
        stop("'conf.level' must be a single number between 0 and 1")
    
    dname <- c("x", "y")
    
    # 3 calculate sd
    if (R == 0) {
        sd1 <- s1
        sd2 <- s2
    } else {
        sd1 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * R * s2^2)/(n1 + n2 - 2))
        sd2 <- sqrt(((n1 - 1) * R^(-1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
    }
    
    # 4 calculate proportions
    if (tail == "upper") {
        prop1 <- 1 - stats::pnorm((cp - m1)/sd1)
        prop2 <- 1 - stats::pnorm((cp - m2)/sd2)
        
    } else {
        prop1 <- stats::pnorm((cp - m1)/sd1)
        prop2 <- stats::pnorm((cp - m2)/sd2)
    }
    
    # 5 calculate propdiff, distrr, distor
    propdiff <- prop1 - prop2
    distrr <- prop1/prop2
    distor <- prop1 * (1 - prop2)/((1 - prop1) * prop2)
    
    # 6 calculate se of propdiff, distrr, distor
    if (R != 0) {
        sediff <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * n2))
        selogrr <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * prop1^2 * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * prop2^2 * n2))
        selogor <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * prop1^2 * (1 - prop1)^2 * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * prop2^2 * (1 - prop2)^2 * n2))
    } else {
        if (tail == "upper") {
            if ((cp - m1) > 0) {
                corr1 <- (0.05/sqrt(n1)) * ((cp - m1)/sd1) * ((cp - m1)/sd1 + 2.6)
                corr2 <- (0.45/sqrt(n1)) * abs((cp - m1)/sd1)^(5/2)
                corr3 <- corr2
            } else {
                corr1 <- (0.05/sqrt(n1)) * ((cp - m1)/sd1) * ((cp - m1)/sd1 - 2.6)
                corr2 <- (0.055/sqrt(n1)) * ((cp - m1)/sd1) * ((cp - m1)/sd1 - 2.6)
                corr3 <- (0.45/sqrt(n1)) * abs((cp - m1)/sd1)^(5/2)
            }
        } else {
            if ((cp - m1) > 0) {
                corr1 <- (0.05/sqrt(n1)) * ((cp - m1)/sd1) * ((cp - m1)/sd1 - 2.6)
                corr2 <- (0.055/sqrt(n1)) * ((cp - m1)/sd1) * ((cp - m1)/sd1 - 2.6)
                corr3 <- (0.45/sqrt(n1)) * abs((cp - m1)/sd1)^(5/2)
            } else {
                corr1 <- -(0.05/sqrt(n1)) * ((cp - m1)/sd1) * ((cp - m1)/sd1 + 2.6)
                corr2 <- (0.45/sqrt(n1)) * abs((cp - m1)/sd1)^(5/2)
                corr3 <- corr2
            }
        }
        sediff <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * n2)) + corr1
        selogrr <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * prop1^2 * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * prop2^2 * n2)) + corr2
        selogor <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * prop1^2 * (1 - prop1)^2 * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * prop2^2 * (1 - prop2)^2 * n2)) + corr2
    }
    
    # 7 recalculate se of rr and or
    alselogrr <- sqrt((exp(selogrr^2) - 1) * exp(2 * log(distrr) + selogrr^2))
    alselogor <- sqrt((exp(selogor^2) - 1) * exp(2 * log(distor) + selogor^2))
    
    
    # 8 confidence intervals
    lev <- 1 - ((1 - conf.level)/2)
    
    ciinf <- propdiff - stats::qnorm(lev) * sediff
    cisup <- propdiff + stats::qnorm(lev) * sediff
    
    ciinfrr <- exp(log(distrr) - stats::qnorm(lev) * selogrr)
    cisuprr <- exp(log(distrr) + stats::qnorm(lev) * selogrr)
    
    ciinfor <- exp(log(distor) - stats::qnorm(lev) * selogor)
    cisupor <- exp(log(distor) + stats::qnorm(lev) * selogor)
    
    # 9 calculate a t-test
    var.equal <- if (R == 1) 
        TRUE else FALSE
    t.testi <- function(n1, m1, s1, n2, m2, s2, alternative = c("two.sided", "less", "greater"), var.equal = FALSE, mu = 0, conf.level = 0.95) {
        alternative <- match.arg(alternative)
        
        dname <- "immediate form (no data)"
        method <- paste(if (!var.equal) 
            "Welch", "Two Sample t-test")
        estimate <- c(m1, m2)
        names(estimate) <- c("mean of x", "mean of y")
        if (var.equal) {
            df <- n1 + n2 - 2
            v <- 0
            if (n1 > 1) 
                v <- v + (n1 - 1) * s1^2
            if (n2 > 1) 
                v <- v + (n2 - 1) * s2^2
            v <- v/df
            stderr <- sqrt(v * (1/n1 + 1/n2))
        } else {
            stderrx <- sqrt(s1^2/n1)
            stderry <- sqrt(s2^2/n2)
            stderr <- sqrt(stderrx^2 + stderry^2)
            df <- stderr^4/(stderrx^4/(n1 - 1) + stderry^4/(n2 - 1))
        }
        if (stderr < 10 * .Machine$double.eps * max(abs(m1), abs(m2))) 
            stop("data are essentially constant")
        tstat <- (m1 - m2 - mu)/stderr
        
        pval <- 2 * stats::pt(-abs(tstat), df)
        alpha <- 1 - conf.level
        cint <- stats::qt(1 - alpha/2, df)
        cint <- tstat + c(-cint, cint)
        
        cint <- mu + cint * stderr
        names(tstat) <- "t"
        names(df) <- "df"
        names(mu) <- "difference in means"
        attr(cint, "conf.level") <- conf.level
        rval <- list(statistic = tstat, parameter = df, p.value = pval, conf.int = cint, estimate = estimate, null.value = mu, alternative = alternative, method = method, data.name = dname)
        class(rval) <- "htest"
        return(rval)
    }
    ttest <- t.testi(n1, m1, s1, n2, m2, s2, var.equal = var.equal)
    
    
    # 11 put all together in a list
    cutpoint <- cp
    names(cutpoint) <- "cutpoint"
    tail <- tail
    names(tail) <- "tail"
    ratio <- R
    names(ratio) <- "R"
    conf.level <- conf.level
    names(conf.level) <- "conf.level"
    correction <- FALSE
    arguments <- list(cutpoint = cutpoint, tail = tail, ratio = ratio, conf.level = conf.level, correction = correction)
    parameter <- c(n1, m1, s1, n2, m2, s2)
    names(parameter) <- c("n1", "m1", "s1", "n2", "m2", "s2")
    method <- paste("Distributional approach with", ifelse(R == "1", "equal", ifelse(R == 0, "unknown", "unequal")), "variances")
    prop <- c(prop1, prop2)
    names(prop) <- c("proportion1", "proportion2")
    estimates <- c(propdiff, distrr, distor)
    names(estimates) <- c("difference in proportions", "risk ratio (distributional estimate)", "odds ratio (distributional estimate)")
    se <- c(sediff, alselogrr, alselogor)
    names(se) <- c("standard error 'propdiff'", "standard error 'rr'", "standard error 'or'")
    ci.diff <- c(ciinf, cisup)
    names(ci.diff) <- c("lower limit 'propdiff'", "upper limit 'propdiff'")
    ci.rr <- c(ciinfrr, cisuprr)
    names(ci.rr) <- c("lower limit 'rr'", "upper limit 'rr'")
    ci.or <- c(ciinfor, cisupor)
    names(ci.or) <- c("lower limit 'or'", "upper limit 'or'")
    ci <- c(ci.diff, ci.rr, ci.or)
    
    res <- list(data.name = dname, arguments = arguments, parameter = parameter, prop = prop, dist.estimates = estimates, se = se, ci = ci, method = method, ttest = ttest)
    
    
    # 12 printing the list
    class(res) <- "distdicho"
    return(res)
}






