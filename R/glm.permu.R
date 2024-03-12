#' @title Permutation Function
#' @description  In practical applications, the assumptions underlying generalized linear models frequently face violations, including incorrect specifications of the outcome variable's distribution or omitted predictors. These deviations can render the results of standard generalized linear models unreliable. As the sample size increases, what might initially appear as minor issues can escalate to critical concerns. To address these challenges, we adopt a permutation-based inference method tailored for generalized linear models. This approach offers robust estimations that effectively counteract the mentioned problems, and its effectiveness remains consistent regardless of the sample size.
#' @param variables a data frame of all the variables.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param npermu the number of permutation times. The default value is 1000.
#' @param CI.percent the confidence level. The default value is 0.95.
#' @importFrom stats glm
#' @return a vector.

### "permu" is for permutation only
permu = function(variables, family, npermu = 1000, CI.percent = .95)
{

  if(npermu < 100){stop("The number of permutation times (npermu) needs to be at least 100.")}

  names(variables)[1] = "y"
  orig_fit = glm(y~., family = family, data = variables)
  orig_beta = summary(orig_fit)$coef[2,1]

  sample_size = length(variables[,1])
  #---------------permutation--------------
  count = 0
  beta.per = rep(NA, npermu) # To store the beta1 after each permutation

  x1_orig = variables$x1

  for (j in 1:npermu){ # for each permutation

    variables$x1 = x1_orig
    idx = sample(1:sample_size) # permutate the index of x1
    variables$x1 = variables$x1[idx] # x1 after the permutation
    fit.per = glm(y~., family = family, data = variables)
    beta.per[j] = summary(fit.per)$coef[2,1] # get the beta1 after permutation
    if (abs(beta.per[j]) >= abs(orig_beta)) {
      count = count + 1
    }

  }

  ### get the decision
  beta.order = beta.per[order(beta.per)] # order all the beta1
  lower = beta.order[npermu*((1 - CI.percent)/2)]
  upper = beta.order[npermu*(1 - (1-CI.percent)/2)]

  ### get the new p-value
  p = count/npermu

  return(c(orig_beta, p, lower, upper))

}


#' @title Permutation-Based Inference for Generalized Linear Models
#' @description  In practical applications, the assumptions underlying generalized linear models frequently face violations, including incorrect specifications of the outcome variable's distribution or omitted predictors. These deviations can render the results of standard generalized linear models unreliable. As the sample size increases, what might initially appear as minor issues can escalate to critical concerns. To address these challenges, we adopt a permutation-based inference method tailored for generalized linear models. This approach offers robust estimations that effectively counteract the mentioned problems, and its effectiveness remains consistent regardless of the sample size.
#' @details In the big data era, the need to reevaluate traditional statistical methods is paramount due to the challenges posed by vast datasets. While larger samples theoretically enhance accuracy and hypothesis testing power without increasing false positives, practical concerns about inflated Type-I errors persist. The prevalent belief is that larger samples can uncover subtle effects, necessitating dual consideration of p-value and effect size. Yet, the reliability of p-values from large samples remains debated.
#'
#' The fact is that larger samples can exacerbate minor issues into significant errors, leading to false conclusions when assumption violations exist. In response, a permutation-based test is introduced to counterbalance the effects of sample size and assumption discrepancies by neutralizing them between actual and permuted data. This approach effectively stabilizes nominal Type I error rates across various sample sizes, thereby ensuring robust statistical inferences even amidst breached conventional assumptions in big data.
#'
#' There are many situations can lead to the assumption violations in generalized linear models such as a scenario of distribution misspecification and a scenario involving unobserved predictors.
#'
#' For example, consider the problem of fitting a Poisson regression to analyze a dataset comprising one outcome variable \eqn{y} and one predictor \eqn{x_1}. The objective is to determine the statistical significance of the predictor’s association with the outcome variable, primarily through the p-value of the regression coefficient for the predictor. In the first scenario, the actual distribution of the outcomes diverges from the Poisson distribution that the model presumes. In the second scenario, outcomes are influenced by an unobserved predictor \eqn{x_2}. Under both situations, the Type I error rates cannot be accurately estimated, and their biases increase as the sample size grows.
#'
#' To utilize an interaction term, a more complex model is required, which cannot be directly applied using this function.
#'
#' @param outcome a vector of the response variable.
#' @param predictors a data frame of all the predictors.
#' @param family a description of the error distribution and link function to be used in the model. We can handle all families supported by glm function.
#' @param npermu the number of permutation times. The default value is 1000.
#' @param CI.percent the confidence level. The default value is 0.95.
#' @return a data frame of estimates of regression coefficients with their permutation p-values and permutation confidence intervals.
#' @export
#' @examples
#' set.seed(0)
#' x1 = rnorm(10, 0, 1)
#' x2 = rnorm(10, 0, 2)
#' x3 = rnorm(10, 0, 3)
#' lambda = exp(x3)
#' y = rpois(10, lambda)
#' X = as.data.frame(cbind(x1, x2))
#' glm.fit = glm(y~., "poisson", data = cbind(y, X))
#' summary(glm.fit)$coef
#' glm.permu(y, X, "poisson")

### "glm.permu" is for shifting the column of predictors based on the "permu" function
glm.permu = function(outcome, predictors, family, npermu = 1000, CI.percent = 0.95)
{
  names_orig = names(predictors)

  nvariable = dim(predictors)[2]

  for (i in 1:nvariable){
    names(predictors)[i] = paste0("x", i)
  }

  result = matrix(NA, nvariable, 4, dimnames = list(names_orig, c("coef", "p-value", "CI_lower", "CI_upper")))
  data.frame(result)

  for (i in 1:nvariable){
    x.shift = predictors[,1]
    predictors[,1] = predictors[, i]
    predictors[, i] = x.shift

    variables = as.data.frame(cbind(outcome, predictors))

    permutation = permu(variables, family, npermu, CI.percent)
    result[i,] = permutation
  }

  return(result)
}
