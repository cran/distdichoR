#' normal, skew-normal or gamma distributed data (via linear regression)
#'
#' Provides adjusted distributional estimates for the comparison of proportions for a dichotomised dependent continuous variable
#' derived from a linear regression of the continuous outcome on the grouping variable and other covariates as described in Sauzet et al. 2015.
#'
#'
#'
#' regdistdicho returns the distributional estimates and their standard errors (see Sauzet et al. 2014 and Peacock et al. 2012)
#' for a difference in proportions, risk ratio and odds ratio. It also provides the distributional confidence intervals for the statistics estimated. 
#' The estimation is based on the marginal means of a linear regression of the outcome on the grouping variable and other covariates.
#'
#'
#'@param mod A linear model of the form lm(lhs ~ rhs) where lhs is a numeric variable giving the data values and
#'           rhs is the grouping variable and other covariates.
#'@param group_var A character string specifying the name of the grouping variable.
#'@param cp A numeric value specifying the cut point under which the distributional proportions are computed.
#'@param tail A character string specifying the tail of the distribution in which the proportions are computed,
#'            must be either 'lower' (default) or 'upper'.
#'@param conf.level Confidence level of the interval.    
#'@param dist A character string specifying the distribution of the error variable in the linear regression, must be either 
#'            'normal' (default), 'sk_normal or 'gamma'. 
#'@param alpha A numeric value specifying further parameter of the skew normal / gamma distribution.       
#'
#'
#'
#'@return A list with class 'distdicho' containing the following components:
#' \item{data.name}{The names of the data.}
#' \item{arguments}{A list with the specified arguments.}
#' \item{parameter}{The marginal mean, standard error and number of observations for both groups.}
#' \item{prop}{The estimated proportions below / above the cut point for both groups.}
#' \item{dist.estimates}{The difference in proportions, risk ratio and odds ratio of the groups.}
#' \item{se}{The estimated standard error of the difference in proportions, the risk ratio and the odds ratio.}
#' \item{ci}{The confidence intervals of the difference in proportions, the risk ratio and the odds ratio.}
#' \item{method}{A character string indicating the used method.}
#'
#'@seealso \code{\link[distdichoR]{distdicho}}, \code{\link[distdichoR]{distdichoi}}, \code{\link[distdichoR]{distdichogen}}, \code{\link[distdichoR]{distdichoigen}}
#'
#'@references
#' Peacock J.L., Sauzet O., Ewings S.M., Kerry S.M. Dichotomising continuous data while retaining statistical power using a distributional approach.  2012 Statist. Med; 26:3089-3103.
#' Sauzet, O., Peacock, J. L. Estimating dichotomised outcomes in two groups with unequal variances: a distributional approach. 2014 Statist. Med; 33 4547-4559 ;DOI: 10.1002/sim.6255.
#' Sauzet, O., Brekenkamp, J., Brenne, S. , Borde, T., David, M., Razum, O., Peacock, J.L. 2015. A distributional approach to obtain adjusted differences in population at risk with a comparison with
#' other regressions methods using perinatal data.  In preparation.
#' Peacock, J.L., Bland, J.M., Anderson, H.R.: Preterm delivery: effects of socioeconomic factors, psychological stress, smoking, alcohol, and caffeine.  BMJ 311(7004), 531-535 (1995).
#'
#'
#'@examples
#' ## Proportions of low birth weight babies among smoking and non-smoking mothers
#' ## (data from Peacock et al. 1995)
#' mod_smoke <- lm(birthwt ~ smoke + gest, data = bwsmoke)
#' regdistdicho(mod = mod_smoke, group_var = 'smoke', cp = 2500, tail = 'lower')
## Ensure that package lsmeans (and nlme for multilevel regression) is installed


#'@export
regdistdicho <- function(mod, group_var, cp = 0, tail = c("lower", "upper"), conf.level = 0.95, dist = c("normal", "sk_normal", "gamma"), alpha=1) {
  
    ### Checking if lsmeans is installed
    if (!is.element("lsmeans", utils::installed.packages()[,1])) {
      utils::install.packages("lsmeans")
    }
    base::require("lsmeans")
    ###

    
    ## Normal Regression
    if (substr(mod$call[1], 1, 3)=="lm") {
      
      # Type of variable
      # Checking if grouping variable is a factor with less than 2 categories
      if (!(is.factor(mod$model[, group_var])) || length(levels(mod$model[, group_var])) < 2) {
        stop("Grouping variable has to be a factor with two or more categories!")
      }
    
      # Finding the reference categorie
      v <- c()
      for(i in 1:length(levels(mod$model[,group_var]))) {
        v[i] <- paste(group_var,as.character(levels(mod$model[,group_var]))[i], sep = "")
      }
      c <- c()
      count <- 1
      for (j in 1:length(v)) {
        for(i in 1:length(names(mod$coefficients))) {
          if(v[j] == as.character(names(mod$coefficients))[i]) {
            c[count] <- j
            count <- count+1
          }
        }
      }
      refnumber <- (1:length(levels(mod$model[, group_var])))[-c]
    
      # If grouping varibale is a factor of more than 2 categories 
      if (length(levels(mod$model[, group_var])) > 2) {
        
        # Setting up counting variables and a list for the outputs
        wert <- 1
        res_list <- rep( list(list()), (length(levels(mod$model[, group_var])) - 1 ) )
        
        # Going through every possible pairing of the categories 
        for(i in c) {
          
          # Residual Std Dev
          rse <- summary(mod)$sigma
        
          # Marginal Effects
          me <- lsmeans::lsmeans(mod, group_var)
          me.1 <- summary(me)$lsmean[i]
          me.2 <- summary(me)$lsmean[refnumber]
        
          # No. of observations
          n.1 <- table(mod$model[, group_var])[i]
          n.2 <- table(mod$model[, group_var])[refnumber]
        
          # Distdichoi and saving it to the list 
          res <- do.call(distdichoigen, list(n1 = n.1, m1 = me.1, s1 = rse, n2 = n.2, m2 = me.2, s2 = rse, cp = cp, tail = tail, conf.level = conf.level, dist = dist, alpha = alpha))
          res$data.name <- c(levels(mod$model[, group_var])[i], levels(mod$model[, group_var])[refnumber])
          #names(res$ttest$estimate) <- c(paste("mean", "of", levels(mod$model[, group_var])[i], sep=" "), paste("mean", "of", levels(mod$model[, group_var])[refnumber], sep=" "))
          res$ttest <- c()
          res_list[[wert]] <- res
          wert <- wert+1
            
        }
        return(res_list)
      }
      if (length(levels(mod$model[, group_var])) == 2) {
  
        # Residual Std Dev
        rse <- summary(mod)$sigma
    
    
        # Marginal Effects
        me <- lsmeans::lsmeans(mod, group_var)
        me.1 <- summary(me)$lsmean[c]
        me.2 <- summary(me)$lsmean[refnumber]
    
        # No. of observations
        n.1 <- table(mod$model[, group_var])[c]
        n.2 <- table(mod$model[, group_var])[refnumber]
    
        res <- do.call(distdichoigen, list(n1 = n.1, m1 = me.1, s1 = rse, n2 = n.2, m2 = me.2, s2 = rse, cp = cp, tail = tail, conf.level = conf.level, dist = dist, alpha = alpha))
        res$data.name <- c(levels(mod$model[, group_var])[c], levels(mod$model[, group_var])[refnumber])
        #names(res$ttest$estimate) <- c(paste("mean", "of", levels(mod$model[, group_var])[c], sep=" "), paste("mean", "of", levels(mod$model[, group_var])[refnumber], sep=" "))  
        res$ttest <- c()
        return(res)
      }
    }
    
    
    
    ## Multilevel
    
    ### Checking if nlme is installed
    if (!is.element("nlme", utils::installed.packages()[,1])) {
      utils::install.packages("nlme")
    }
    base::require("nlme")
    ###
    
    if (substr(mod$call[1], 1, 3)=="lme") {
      # Type of variable
      # Checking if grouping variable is a factor with less than 2 categories
      if (!(is.factor( mod$data[[which(names(mod$data) == group_var)]])) || length(table(mod$data[, group_var])) < 2) {
        stop("Grouping variable has to be a factor with two or more categories!")
      }
      
      # Finding the reference categorie
      v <- c()
      for(i in 1:length(table(mod$data[, group_var]))) {
        v[i] <- paste(group_var,as.character(names(table(mod$data[, group_var])))[i], sep = "")
      }
      c <- c()
      count <- 1
      for (j in 1:length(v)) {
        for(i in 1:length(colnames(mod$varFix))) {
          if(v[j] == as.character(colnames(mod$varFix))[i]) {
            c[count] <- j
            count <- count+1
          }
        }
      }
      refnumber <- (1:length(table(mod$data[, group_var])))[-c]
      
      # If grouping varibale is a factor of more than 2 categories 
      if (!(is.factor( mod$data[[which(names(mod$data) == group_var)]])) || length(table(mod$data[, group_var])) > 2) {
        
        # Setting up counting variables and a list for the outputs
        wert <- 1
        res_list <- rep( list(list()), (length(table(mod$data[, group_var])) - 1 ) )
        
        # Going through every possible pairing of the categories 
        for(i in c) {
          
          # Residual Std Dev
          rse <- sqrt(as.numeric(nlme::VarCorr(mod)[1,1])+as.numeric(nlme::VarCorr(mod)[2,1]))
          
          # Marginal Effects
          me <- lsmeans::lsmeans(mod, group_var)
          me.1 <- summary(me)$lsmean[i]
          me.2 <- summary(me)$lsmean[refnumber]
          
          # No. of observations
          n.1 <- table(mod$data[, group_var])[i]
          n.2 <- table(mod$data[, group_var])[refnumber]
          
          # Distdichoi and saving it to the list 
          res <- do.call(distdichoigen, list(n1 = n.1, m1 = me.1, s1 = rse, n2 = n.2, m2 = me.2, s2 = rse, cp = cp, tail = tail, conf.level = conf.level, dist = dist, alpha = alpha))
          res$data.name <- c(names(table(mod$data[, group_var]))[i], names(table(mod$data[, group_var]))[refnumber])
          #names(res$ttest$estimate) <- c(paste("mean", "of", names(table(mod$data[, group_var]))[i], sep=" "), paste("mean", "of", names(table(mod$data[, group_var]))[refnumber], sep=" "))
          res$ttest <- c()
          res_list[[wert]] <- res
          wert <- wert+1
          
        }
        return(res_list)
      }
      if (!(is.factor( mod$data[[which(names(mod$data) == group_var)]])) || length(table(mod$data[, group_var])) == 2) {
        
        # Residual Std Dev
        rse <- sqrt(as.numeric(nlme::VarCorr(mod)[1,1])+as.numeric(nlme::VarCorr(mod)[2,1]))
        
        
        # Marginal Effects
        me <- lsmeans::lsmeans(mod, group_var)
        me.1 <- summary(me)$lsmean[c]
        me.2 <- summary(me)$lsmean[refnumber]
        
        # No. of observations
        n.1 <- table(mod$data[, group_var])[c]
        n.2 <- table(mod$data[, group_var])[refnumber]
        
        res <- do.call(distdichoigen, list(n1 = n.1, m1 = me.1, s1 = rse, n2 = n.2, m2 = me.2, s2 = rse, cp = cp, tail = tail, conf.level = conf.level, dist = dist, alpha = alpha))
        res$data.name <- c(names(table(mod$data[, group_var]))[c], names(table(mod$data[, group_var]))[refnumber])
        #names(res$ttest$estimate) <- c(paste("mean", "of", names(table(mod$data[, group_var]))[c], sep=" "), paste("mean", "of", names(table(mod$data[, group_var]))[refnumber], sep=" "))  
        res$ttest <- c()
        return(res)
      }
    } 
}