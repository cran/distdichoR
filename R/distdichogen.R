#' normal, skew-normal or gamma distributed data
#'
#' distdichogen first returns the results of a two-group unpaired t-test.
#' Followed by the distributional estimates and their standard errors (see Sauzet et al. 2014 and Peacock et al. 2012)
#' for a difference in proportions, risk ratio and odds ratio. It also provides the distributional confidence intervals for the statistics estimated.
#' distdicho_gen takes normal (dist = 'normal'), skew normal (dist = 'sk_normal') and gamma (dist = 'gamma') distributed data.
#' The data can either be given as two variables, which provide the outcome in each group or specified as a formula
#' of the form lhs ~ rhs where lhs is a numeric variable giving the data values and
#' rhs a factor with two levels giving the corresponding exposed and unexposed groups.
#' In all cases, it is assumed that there are only two groups.
#'
#'
#'
#'
#'
#'
#'@param x A numeric vector of data values.
#'@param y A numeric vector of data values.
#'@param cp A numeric value specifying the cut point under which the distributional proportions are computed.
#'@param tail A character string specifying the tail of the distribution in which the proportions are computed.
#'            Must be either 'lower' (default) or 'upper'.
#'@param conf.level Confidence level of the interval.
#'@param dist A character string specifying the distribution of the data. Must be either 'normal' (default), 'sk_normal or 'gamma'.
#'@param bootci A logical variable indicating whether bootstrap bias-corrected confidence intervals are calculated
#'        instead of distributional ones.
#'@param nrep A numeric value, specifies the number of bootstrap replications (nrep must be higher than the number of observations).
#'@param formula A formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and
#'         rhs a factor with two levels giving the corresponding exposed and unexposed groups.
#'@param data An optional matrix or data frame containing the variables
#'            in the formula. By default, the variables are taken from \code{environment(formula)}.
#'@param exposed A character string specifying the grouping value of the exposed group.
#'@param ... Further arguments to be passed to or from methods.            
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
#' @seealso \code{\link[distdichoR]{distdicho}}, \code{\link[distdichoR]{distdichoi}}, \code{\link[distdichoR]{distdichoigen}}, \code{\link[distdichoR]{regdistdicho}}
#'
#'@references
#' Peacock J.L., Sauzet O., Ewings S.M., Kerry S.M. Dichotomising continuous data while retaining statistical power using a distributional approach.  Statist. Med; 2012; 26:3089-3103.
#' Sauzet, O., Peacock, J. L. Estimating dichotomised outcomes in two groups with unequal variances: a distributional approach.  Statist. Med; 2014 33 4547-4559 ;DOI: 10.1002/sim.6255.
#' Sauzet, O., Ofuya, M., Peacock, J. L. Dichotomisation using a distributional approach when the outcome is skewed BMC Medical Research Methodology 2015, 15:40; doi:10.1186/s12874-015-0028-8.
#' Peacock, J.L., Bland, J.M., Anderson, H.R.: Preterm delivery: effects of socioeconomic factors, psychological stress, smoking, alcohol, and caffeine. BMJ 311(7004), 531-535 (1995).
#'
#'
#'@examples
#'## Proportions of low birth weight babies among smoking and non-smoking mothers
#'## (data from Peacock et al. 1995). Returns distributional estimates, standard 
#'## errors and distributional confidence intervals for differences in proportions,
#'## RR and OR of babies having a birth weight under 2500g (low birth weight)
#'## for group smoker (mother smokes) over the odds of LBW in group non-smoker 
#'## (mother doesn't smoke)
#'# Formula interface
#'distdichogen(birthwt ~ smoke, cp = 2500, data = bwsmoke, exposed = 'smoker',
#'              dist = 'sk_normal')
#'# Data stored in two vectors
#' bw_smoker <- bwsmoke$birthwt[bwsmoke$smoke == 'smoker']
#' bw_nonsmoker <- bwsmoke$birthwt[bwsmoke$smoke == 'non-smoker']
#' distdichogen(x = bw_smoker, y = bw_nonsmoker, 
#'               cp = 2500, tail = 'lower', dist = 'sk_normal')
#'
#'
#'## Body Mass Index (BMI) and parity. Returns distributional estimates, standard
#'## errors and distributional confidence intervals for difference in proportions,
#'## RR and OR of obese mothers (BMI of >30kg/m^2) for group_par=1 (multiparity) 
#'## over the odds of obesity in group_par=0 (primiparity)
#'distdichogen(bmi ~ group_par, cp = 30, data = bmi, exposed = '1',
#'              tail = 'upper', dist = 'sk_normal')
#'
#'
#'
#'
#'@export
distdichogen <- function(x, ...) {
    UseMethod("distdichogen")
}

#'@rdname distdichogen
#'@export
distdichogen.default <- function(x, y, cp = 0, tail = c("lower", "upper"), conf.level = 0.95, dist = c("normal", "sk_normal", "gamma"), bootci = FALSE, nrep = 2000, ...) {
    
    ### Checking if sn is installed
    if (!is.element("sn", utils::installed.packages()[,1])) {
      utils::install.packages("sn")
    }
    base::require("sn")  
    ###
    
    dist <- match.arg(dist)
    # normal distribution ###############################################################################################
    if (dist == "normal") {
        
        # 1 verify arguments
        tail <- match.arg(tail)
        
        if (length(x) <= 1 || !is.numeric(x) || any(is.infinite(x))) 
            stop("'x' must be a numeric vector")
        
        if (length(y) <= 1 || !is.numeric(y) || any(is.infinite(y))) 
            stop("'y' must be a numeric vector")
        
        if (length(cp) != 1 || !is.numeric(cp) || is.infinite(cp)) 
            stop("'cp' must be a single number")
        
        if (length(conf.level) != 1 || !is.numeric(conf.level) || conf.level < 0 || conf.level > 1 || is.infinite(conf.level)) 
            stop("'conf.level' must be a single number between 0 and 1")
        
        if (length(nrep) != 1 || !is.numeric(nrep) || nrep%%1 != 0 || nrep <= 0 || is.infinite(nrep)) 
            stop("'nrep' must be a single positive number")
        
        
        dname <- c(deparse(substitute(x)), deparse(substitute(y)))
        
        
        # 2 compute n, m, s
        xok <- !is.na(x)
        yok <- !is.na(y)
        
        x <- x[xok]
        y <- y[yok]
        
        n1 <- length(x)
        n2 <- length(y)
        s1 <- stats::sd(x)
        s2 <- stats::sd(y)
        m1 <- mean(x)
        m2 <- mean(y)
        
        
        # 3 calculate sd
        sd1 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
        sd2 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
        
        
        
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
        
        
        # 6 calculate se of Propdiff, distrr, distor
        sediff <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * n2))
        selogrr <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * prop1^2 * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * prop2^2 * n2))
        selogor <- sqrt(exp(-(cp - m1)^2/sd1^2)/(2 * pi * prop1^2 * (1 - prop1)^2 * n1) + exp(-(cp - m2)^2/sd2^2)/(2 * pi * prop2^2 * (1 - prop2)^2 * n2))
        
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
        ttest <- stats::t.test(x, y, var.equal = TRUE, conf.level = conf.level)
        
        
        # 10 use bootstrap ci
        if (bootci == TRUE) {
            data <- data.frame(var1 = c(x, y), var2 = c(rep(1, length(x)), rep(2, length(y))))
            boot.fkt <- function(data, indices, varr = 1, cp, tail = "lower") {
                d <- data[indices, ]
                x <- d[d$var2 == 1, 1]
                y <- d[d$var2 == 2, 1]
                
                n1 <- length(x)
                n2 <- length(y)
                s1 <- stats::sd(x)
                s2 <- stats::sd(y)
                m1 <- mean(x)
                m2 <- mean(y)
                
                if (varr == 0) {
                  sd1 <- s1
                  sd2 <- s2
                } else {
                  sd1 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * varr * s2^2)/(n1 + n2 - 2))
                  sd2 <- sqrt(((n1 - 1) * varr^(-1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
                }
                
                if (tail == "upper") {
                  prop1 <- 1 - stats::pnorm((cp - m1)/sd1)
                  prop2 <- 1 - stats::pnorm((cp - m2)/sd2)
                } else {
                  prop1 <- stats::pnorm((cp - m1)/sd1)
                  prop2 <- stats::pnorm((cp - m2)/sd2)
                }
                
                r <- numeric(3)
                r[1] <- prop1 - prop2
                r[2] <- prop1/prop2
                r[3] <- prop1 * (1 - prop2)/((1 - prop1) * prop2)
                return(r)
            }
            boot.out <- boot::boot(data = data, statistic = boot.fkt, R = nrep, strata = data$var2, cp = cp, tail = tail, varr = 1)
            
            boot.ci.diff <- boot::boot.ci(boot.out = boot.out, index = 1, conf = conf.level, type = "bca")
            boot.ci.rr <- boot::boot.ci(boot.out = boot.out, index = 2, conf = conf.level, type = "bca")
            boot.ci.or <- boot::boot.ci(boot.out = boot.out, index = 3, conf = conf.level, type = "bca")
            # modify ci
            ciinf <- boot.ci.diff$bca[4]
            ciinfrr <- boot.ci.rr$bca[4]
            ciinfor <- boot.ci.or$bca[4]
            
            cisup <- boot.ci.diff$bca[5]
            cisuprr <- boot.ci.rr$bca[5]
            cisupor <- boot.ci.or$bca[5]
        }
        
        # 11 put all together in a list
        cutpoint <- cp
        names(cutpoint) <- "cutpoint"
        tail <- tail
        names(tail) <- "tail"
        conf.level <- conf.level
        names(conf.level) <- "conf.level"
        bootci <- bootci
        names(bootci) <- "bootci"
        arguments <- list(cutpoint = cutpoint, tail = tail, conf.level = conf.level, bootci = bootci)
        parameter <- c(n1, m1, s1, n2, m2, s2)
        names(parameter) <- c("n1", "m1", "s1", "n2", "m2", "s2")
        method <- paste("Distributional approach with equal variances")
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
        
        # skew normal distribution ###############################################################################################
    } else if (dist == "sk_normal") {
        
        # 1 verify arguments
        tail <- match.arg(tail)
        
        if (length(x) <= 1 || !is.numeric(x) || any(is.infinite(x))) 
            stop("'x' must be a numeric vector")
        
        if (length(y) <= 1 || !is.numeric(y) || any(is.infinite(y))) 
            stop("'y' must be a numeric vector")
        
        if (length(cp) != 1 || !is.numeric(cp) || is.infinite(cp)) 
            stop("'cp' must be a single number")
        
        if (length(conf.level) != 1 || !is.numeric(conf.level) || conf.level < 0 || conf.level > 1 || is.infinite(conf.level)) 
            stop("'conf.level' must be a single number between 0 and 1")
        
        if (length(nrep) != 1 || !is.numeric(nrep) || nrep%%1 != 0 || nrep <= 0 || is.infinite(nrep)) 
            stop("'nrep' must be a single positive number")
        
        dname <- c(deparse(substitute(x)), deparse(substitute(y)))
        
        # 2 compute n, m, s
        xok <- !is.na(x)
        yok <- !is.na(y)
        
        x <- x[xok]
        y <- y[yok]
        
        n1 <- length(x)
        n2 <- length(y)
        s1 <- stats::sd(x)
        s2 <- stats::sd(y)
        m1 <- mean(x)
        m2 <- mean(y)
        
        
        # 3 calculate sd, muz, w, alpha
        alpha <- stats::coef(sn::selm(c(x, y) ~ 1), "dp")["alpha"]
        muz <- sqrt(2/pi) * alpha/sqrt(1 + alpha^2)
        sd1 <- sd2 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
        w1 <- w2 <- sd1/sqrt(1 - muz^2)
        alpha21 <- m1 - w1 * muz
        alpha22 <- m2 - w2 * muz
        
        
        
        
        # 4 calculate proportions
        if (tail == "upper") {
            prop1 <- 1 - sn::psn(cp, alpha21, w1, alpha)
            prop2 <- 1 - sn::psn(cp, alpha22, w2, alpha)
            
        } else {
            prop1 <- sn::psn(cp, alpha21, w1, alpha)
            prop2 <- sn::psn(cp, alpha22, w2, alpha)
        }
        
        
        
        # 5 calculate propdiff, distrr, distor
        propdiff <- prop1 - prop2
        distrr <- prop1/prop2
        distor <- prop1 * (1 - prop2)/((1 - prop1) * prop2)
        
        
        # 6 calculate se of Propdiff, distrr, distor
        sediff <- sqrt((w1^2/n1) * (1 - muz^2) * ((2 * (-exp(-(cp - alpha21)^2/(2 * w1^2))/sqrt(2 * pi * w1^2)) * stats::pnorm(alpha * (cp - alpha21)/w1))^2 + (2 * (-exp(-(cp - alpha22)^2/(2 * 
            w1^2))/sqrt(2 * pi * w1^2)) * stats::pnorm(alpha * (cp - alpha22)/w1))^2))
        selogrr <- sqrt((w1^2/n1) * (1 - muz^2) * ((2 * (-exp(-(cp - alpha21)^2/(2 * w1^2))/sqrt(2 * pi * w1^2)) * stats::pnorm(alpha * (cp - alpha21)/w1))^2/prop1^2 + (2 * (-exp(-(cp - alpha22)^2/(2 * 
            w1^2))/sqrt(2 * pi * w1^2)) * stats::pnorm(alpha * (cp - alpha22)/w1)/prop2)^2))
        selogor <- sqrt((w1^2/n1) * (1 - muz^2) * ((2 * (-exp(-(cp - alpha21)^2/(2 * w1^2))/sqrt(2 * pi * w1^2)) * stats::pnorm(alpha * (cp - alpha21)/w1))^2/(prop1 * (1 - prop1))^2 + (2 * (-exp(-(cp - 
            alpha22)^2/(2 * w1^2))/sqrt(2 * pi * w1^2)) * stats::pnorm(alpha * (cp - alpha22)/w1)/(prop2 * (1 - prop2)))^2))
        
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
        ttest <- stats::t.test(x, y, var.equal = TRUE, conf.level = conf.level)
        
        
        # 10 use bootstrap ci
        if (bootci == TRUE) {
            data <- data.frame(var1 = c(x, y), var2 = c(rep(1, length(x)), rep(2, length(y))))
            boot.fkt2 <- function(data, indices, cp, tail = "lower") {
                d <- data[indices, ]
                x <- d[d$var2 == 1, 1]
                y <- d[d$var2 == 2, 1]
                n1 <- length(x)
                n2 <- length(y)
                s1 <- stats::sd(x)
                s2 <- stats::sd(y)
                m1 <- mean(x)
                m2 <- mean(y)
                
                alpha <- stats::coef(sn::selm(c(x, y) ~ 1), "dp")["alpha"]
                muz <- sqrt(2/pi) * alpha/sqrt(1 + alpha^2)
                sd1 <- sd2 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
                w1 <- w2 <- sd1/sqrt(1 - muz^2)
                alpha21 <- m1 - w1 * muz
                alpha22 <- m2 - w2 * muz
                
                if (tail == "upper") {
                  prop1 <- 1 - sn::psn(cp, alpha21, w1, alpha)
                  prop2 <- 1 - sn::psn(cp, alpha22, w2, alpha)
                } else {
                  prop1 <- sn::psn(cp, alpha21, w1, alpha)
                  prop2 <- sn::psn(cp, alpha22, w2, alpha)
                }
                
                r <- numeric(3)
                r[1] <- prop1 - prop2
                r[2] <- prop1/prop2
                r[3] <- prop1 * (1 - prop2)/((1 - prop1) * prop2)
                return(r)
            }
            boot.out <- boot::boot(data = data, statistic = boot.fkt2, R = nrep, strata = data$var2, cp = cp, tail = tail)
            boot.ci.diff <- boot::boot.ci(boot.out = boot.out, index = 1, conf = 0.95, type = "bca")
            boot.ci.rr <- boot::boot.ci(boot.out = boot.out, index = 2, conf = 0.95, type = "bca")
            boot.ci.or <- boot::boot.ci(boot.out = boot.out, index = 3, conf = 0.95, type = "bca")
            # modify ci
            ciinf <- boot.ci.diff$bca[4]
            ciinfrr <- boot.ci.rr$bca[4]
            ciinfor <- boot.ci.or$bca[4]
            
            cisup <- boot.ci.diff$bca[5]
            cisuprr <- boot.ci.rr$bca[5]
            cisupor <- boot.ci.or$bca[5]
        }
        
        
        # 11 put all together in a list
        cutpoint <- cp
        names(cutpoint) <- "cutpoint"
        tail <- tail
        names(tail) <- "tail"
        alpha <- alpha
        names(alpha) <- alpha
        conf.level <- conf.level
        names(conf.level) <- "conf.level"
        bootci <- bootci
        names(bootci) <- "bootci"
        arguments <- list(cutpoint = cutpoint, tail = tail, alpha = alpha, conf.level = conf.level, bootci = bootci)
        parameter <- c(n1, m1, s1, n2, m2, s2)
        names(parameter) <- c("n1", "m1", "s1", "n2", "m2", "s2")
        method <- paste("Distributional approach when the outcome is skewed")
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
        
        # gamma distribution ###############################################################################################
    } else {
        # 1 verify arguments
        tail <- match.arg(tail)
        
        if (length(x) <= 1 || !is.numeric(x) || any(is.infinite(x))) 
            stop("'x' must be a numeric vector")
        
        if (length(y) <= 1 || !is.numeric(y) || any(is.infinite(y))) 
            stop("'y' must be a numeric vector")
        
        if (length(cp) != 1 || !is.numeric(cp) || is.infinite(cp)) 
            stop("'cp' must be a single number")
        
        if (length(conf.level) != 1 || !is.numeric(conf.level) || conf.level < 0 || conf.level > 1 || is.infinite(conf.level)) 
            stop("'conf.level' must be a single number between 0 and 1")
        
        if (length(nrep) != 1 || !is.numeric(nrep) || nrep%%1 != 0 || nrep <= 0 || is.infinite(nrep)) 
            stop("'nrep' must be a single positive number")
        
        dname <- c(deparse(substitute(x)), deparse(substitute(y)))
        
        # 2 compute n, m, s
        xok <- !is.na(x)
        yok <- !is.na(y)
        
        x <- x[xok]
        y <- y[yok]
        
        n1 <- length(x)
        n2 <- length(y)
        s1 <- stats::sd(x)
        s2 <- stats::sd(y)
        m1 <- mean(x)
        m2 <- mean(y)
        
        
        # 3 calculate sd, alpha
        sd12 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
        v <- (n1 + n2 - 1) * (s1^2 + s2^2)/(n1 + n2 -2)
        alpha <- 0.5 * (m1 + m2)^2/v
        
        
        
        # 4 calculate proportions
        if (tail == "upper") {
            prop1 <- 1 - stats::pgamma(cp, rate = alpha/m1, shape = alpha)
            prop2 <- 1 - stats::pgamma(cp, rate = alpha/m2, shape = alpha)
        } else {
            prop1 <- stats::pgamma(cp, rate = alpha/m1, shape = alpha)
            prop2 <- stats::pgamma(cp, rate = alpha/m2, shape = alpha)
        }
        
        
        
        # 5 calculate propdiff, distrr, distor
        propdiff <- prop1 - prop2
        distrr <- prop1/prop2
        distor <- prop1 * (1 - prop2)/((1 - prop1) * prop2)
        
        
        # 6 calculate se of Propdiff, distrr, distor
        
        sediff <- sqrt(1/(alpha)) * ((alpha * cp)^alpha/gamma(alpha)) * sqrt(1/n1 * (1/m1)^(2 * alpha) * exp(-2 * alpha * cp/m1) + 1/n2 * (1/m2)^(2 * alpha) * exp(-2 * alpha * cp/m2))
        selogrr <- sqrt(1/(alpha)) * ((alpha * cp)^alpha/gamma(alpha)) * sqrt(1/(n1 * prop1) * (1/m1)^(2 * alpha) * exp(-2 * alpha * cp/m1) + 1/(n2 * prop2) * (1/m2)^(2 * alpha) * exp(-2 * 
            alpha * cp/m2))
        selogor <- sqrt(1/(alpha)) * ((alpha * cp)^alpha/gamma(alpha)) * sqrt(1/(n1 * prop1 * (1 - prop1)) * (1/m1)^(2 * alpha) * exp(-2 * alpha * cp/m1) + 1/(n2 * prop2 * (1 - prop2)) * 
            (1/m2)^(2 * alpha) * exp(-2 * alpha * cp/m2))
        
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
        ttest <- stats::t.test(x, y, var.equal = TRUE, conf.level = conf.level)
        
        
        # 10 use bootstrap ci
        if (bootci == TRUE) {
            data <- data.frame(var1 = c(x, y), var2 = c(rep(1, length(x)), rep(2, length(y))))
            boot.fkt3 <- function(data, indices, cp, tail = "lower") {
                d <- data[indices, ]
                x <- d[d$var2 == 1, 1]
                y <- d[d$var2 == 2, 1]
                n1 <- length(x)
                n2 <- length(y)
                s1 <- stats::sd(x)
                s2 <- stats::sd(y)
                m1 <- mean(x)
                m2 <- mean(y)
                
                sd12 <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + n2 - 2))
                
                v <- (n1 + n2 - 1) * (s1^2 + s2^2)/(n1 + n2 -2)
                alpha <- 0.5 * (m1 + m2)^2/v
                
                if (tail == "upper") {
                  prop1 <- 1 - stats::pgamma(cp, rate = alpha/m1, shape = alpha)
                  prop2 <- 1 - stats::pgamma(cp, rate = alpha/m2, shape = alpha)
                } else {
                  prop1 <- stats::pgamma(cp, rate = alpha/m1, shape = alpha)
                  prop2 <- stats::pgamma(cp, rate = alpha/m2, shape = alpha)
                }
                
                r <- numeric(3)
                r[1] <- prop1 - prop2
                r[2] <- prop1/prop2
                r[3] <- prop1 * (1 - prop2)/((1 - prop1) * prop2)
                return(r)
            }
            boot.out <- boot::boot(data = data, statistic = boot.fkt3, R = nrep, strata = data$var2, cp = cp, tail = tail)
            boot.ci.diff <- boot::boot.ci(boot.out = boot.out, index = 1, conf = 0.95, type = "bca")
            boot.ci.rr <- boot::boot.ci(boot.out = boot.out, index = 2, conf = 0.95, type = "bca")
            boot.ci.or <- boot::boot.ci(boot.out = boot.out, index = 3, conf = 0.95, type = "bca")
            # modify ci
            ciinf <- boot.ci.diff$bca[4]
            ciinfrr <- boot.ci.rr$bca[4]
            ciinfor <- boot.ci.or$bca[4]
            
            cisup <- boot.ci.diff$bca[5]
            cisuprr <- boot.ci.rr$bca[5]
            cisupor <- boot.ci.or$bca[5]
        }
        
        
        
        # 11 put all together in a list
        cutpoint <- cp
        names(cutpoint) <- "cutpoint"
        tail <- tail
        names(tail) <- "tail"
        alpha <- alpha
        names(alpha) <- alpha
        conf.level <- conf.level
        names(conf.level) <- "conf.level"
        bootci <- bootci
        names(bootci) <- "bootci"
        arguments <- list(cutpoint = cutpoint, tail = tail, alpha = alpha, conf.level = conf.level, bootci = bootci)
        parameter <- c(n1, m1, s1, n2, m2, s2)
        names(parameter) <- c("n1", "m1", "s1", "n2", "m2", "s2")
        method <- paste("Distributional approach when the outcome is gamma distributed")
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
}

#'@rdname distdichogen
#'@export
distdichogen.formula <- function(formula, data, exposed, ...) {
    if (missing(formula) || (length(formula) != 3L) || (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$exposed <- NULL
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (levels(g)[1] != exposed && levels(g)[2] != exposed)
      stop("argument exposed is invalid")
    if (levels(g)[1] != exposed) {
        g <- factor(g, levels(g)[c(2, 1)])
    }
    DNAME <- levels(g)
    if (nlevels(g) != 2L) 
        stop("grouping factor must have exactly 2 levels")
    DATA <- stats::setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("distdichogen", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}




