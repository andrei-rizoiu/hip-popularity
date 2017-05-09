require(pracma)
require(nloptr)

#' Calculates, based on a structure of lambda history ( pairs (Time, Count) ) 
#' the average value of the number of viewcounts (lambda) at the given moment t.
#' The lambda_history can either by a past history of event counts (observed or 
#' simulated) or NULL (default) if NULL, an artificial series is generated with 
#' parameters until time t-1, and then value at time t is generated. The 
#' external influence should be a vector of the length higher or equal than the
#' desired period + 1 (for moment 0 of time). NULL is no influence (a value of
#' zero).
#' 
#' @param params (gamma, eta, K, beta, c, theta, mu1, mu2, ...) and alpha: simulation 
#'   parameters
#' @param ext_infl: the series containing S(t) - the external influence in our 
#'   model
predict_theoretical_lambda <- function(params = c(gamma = 100, eta = 100, K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, mu1 = 1), t, 
                                       lambda_history = NULL, ext_infl = NULL, alpha = 2.016) {
  res <- .check_fix_mus_ext_infl(params = params, ext_infl = ext_infl, t = t)
  params <- res$params
  ext_infl <- res$ext_infl
  
  if (t == 0) {
    # remember that the first element is influence at time 0, 
    # ext_infl[i] is the influence at time i-1
    infl_val <- 0
    for (i in 1:length(ext_infl)) {
      muname <- sprintf("mu%d", i)
      infl_val <- infl_val + params[[muname]] * ext_infl[[i]][1]
    }
    return(params$gamma + infl_val) 
  }
  
  if (! is.null(lambda_history)) {
    # the event has 2 components: (magnitude_i, time_i)
    mat_event = matrix(unlist(lambda_history), ncol = 2, byrow = F)
    timei = mat_event[,1]
    lambdai = mat_event[,2]
    
    # eliminate values which do not interest us (scope of time is above our t)
    lambdai <- lambdai[timei < t]
    timei <- timei[timei < t]
    
    # calculate the time difference between each lambdai and the current moment
    diffi <- t - timei
    
    # calculate the rest of the formula
    # 1. time decay part
    result <- (diffi + params$c) ^ (-1-params$theta)
    # 2. multiply the values of lambda
    result <- result * lambdai
    # 3. add the external influence
    infl_val <- 0
    for (i in 1:length(ext_infl)) {
      muname <- sprintf("mu%d", i)
      infl_val <- infl_val + params[[muname]] * ext_infl[[i]][t+1]
    }
    # 4. sum everything up and finish calculations
    result <- params$eta  + infl_val + params$K * ((alpha - 1) / (alpha - params$beta - 1)) * sum(result)
  } else {
    series <- generate_simulated_data(params = params, time = t, ext_infl = ext_infl)
    result <- series$Count[series$Time == t]
  }
  
  (result)
}

############################ Functions for testing the fitting on simulated data
#' Function used to generate viewcount data (with or without noise), starting 
#' from an initial value. it calls predic_theoretical_lambda to get next value 
#' based on past values. The length of the generated list is time
#' 
#' @param params (gamma, eta, mu, K, beta, c, theta) and alpha: simulation 
#'   parameters
#' @param time: how many timeslices should be simulated
#' @param ext_infl: the series containing S(t) - the external influence in our 
#'   model
#' @param noise: should uniform noise be added to the generated series
#' @param multiplicative: if noise, should the noise be a proportion of the 
#'   generated values, or the generated noise should simply be added to the 
#'   values E.g.: multiplicative = true: counts = counts + counts*factor, else:
#'   counts = counts + factor
#' @param factor: the factor of the noise
generate_simulated_data <- function(params = c(gamma = 100, eta = 100, K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, mu1 = 1), time,
                                    ext_infl = NULL, alpha = 2.016, noise = FALSE, factor = 0.3, multiplicative = TRUE, prefix = NULL) {
  # add the initial event
  event_list <- data.frame(matrix(data = NA, nrow = time+1, ncol = 2), row.names = NULL)
  colnames(event_list) <- c("Time", "Count")
  event_list$Time <- 0:time
  if (is.null(prefix))
    prefix <- predict_theoretical_lambda(params = params, t = 0, ext_infl = ext_infl, alpha = alpha)
  
  # initialize event count list
  upper <- length(prefix)
  if (upper > time + 1) upper <- time + 1
  event_list$Count[1:upper] <- prefix[1:upper]
  next_t <- upper
  
  # go through each moment of time
  if ( next_t <= time) {
    for (t in next_t:time) {
      event_list$Count[t+1] <- predict_theoretical_lambda(params = params, t = t, lambda_history = event_list, ext_infl = ext_infl, alpha = alpha)
    }
  }
  
  # add noise to the generated series
  if (noise) {
    if (multiplicative) {
      # generate gaussian proportions between -factor to factor
      factor = abs(factor)
      noise <- runif(n = time, min = -1 * factor, max = factor) * event_list$Count[-1]
    } else {
      factor = abs(factor)
      noise <- runif(n = time, min = -1 * factor, max = factor)
    }
    
    # add the noise to the generated event list
    event_list$Count[-1] <- event_list$Count[-1] + noise
  }
  
  (event_list)
}

#' this function will need to be minimized in order to fit the theoretical 
#' parameters to the real observed values
#' 
#' @param params a list or vector with the models parameters.
#' @param real_lambda_history the observed event counts to which to compare when
#'   calculating error.
#' @param ext_infl: the series containing S(t) - the external influence in our 
#'   model
#' @param alpha (default 2.016) the value of alpha
#' @param return_error_vector (default F) if True, return the error vector 
#'   instead of just the sum
#' @param lowerBound (default NULL) - a vector of values which are the lower 
#'   bounds for the parameters. If a given parameter is lower than its lower 
#'   bound, than the error function is set to a very large value.
#' @param upperBound (default NULL) - same as lowerBound, just for the upper
#'   bound.
#' @param alpha_regularizer (default 0) - metaparameter for regularizing the
#'   values of parameters. If zero, no regularization.
#' @param param_scaling (default NULL) - vector of scaling values for each
#'   parameter. If NULL, a dummy vector of ones will be used.
error_function <- function(params, real_lambda_history, ext_infl = NULL, alpha = 2.016, 
                           return_error_vector = FALSE, disable_gradient_params = NULL,
                           lowerBound = NULL, upperBound = NULL,
                           alpha_regularizer = 0, param_scaling = NULL) {
  
  # obtain the parameters from the list of of parameters. This will be coming 
  # from the optimization procedure process parameters
  ## also correct the parameters to be a names list.
  params <- .correct_names(params)
  
  # check if we got bounds
  if (is.null(lowerBound)) lowerBound <- rep(-Inf, times = length(params))
  if (is.null(upperBound)) upperBound <- rep(Inf, times = length(params))
  
  # check if we got params scaling
  if (is.null(param_scaling))
    param_scaling <- rep(x = 1, times = length(params))
  param_scaling <- unlist(param_scaling)
  
  # sanity check - check if the normalization params are the same length as the params
  # at this stage, they should be compatible
  if (length(params) != length(param_scaling))
    stop(sprintf("[error_function]: You have %d parameters and %s normalization parameters. Cannot work in these conditions!", length(params), length(param_scaling)))
  
  # calculate the theoretical values of lambda, based on our parameters and the
  # known external influence. this boils down to simulating the series with the
  # parameters
  theoretical_lambda <- generate_simulated_data(params = params, 
                                                time = max(real_lambda_history$Time), 
                                                ext_infl = ext_infl, 
                                                alpha = alpha)
  
  # calculate the error
  result <- (real_lambda_history$Count - theoretical_lambda$Count) ^2 / 2
  if (return_error_vector) {
    return(result)  
  }
  result <- sum(result)
  
  # calculate the regularization
  #   disable_reg <- c("theta", "c")
  disable_reg <- c(5, 6)
  regularized_sum <- (as.numeric(params) / param_scaling) ^ 2
  regularized_sum[disable_reg] <- 0
  reg <- 0.5 * alpha_regularizer * sum(regularized_sum, na.rm = T)
  result <- result + reg
  
  # L-BFGS-B doesn't like Inf values for this function, therefore we will replace Inf with a ridiculously large value
  maxValue <- .Machine$double.xmax - 1
  
  if ( !is.finite(result) ) result <- maxValue
  if ( is.nan(result) ) result <- maxValue
  if ( result > maxValue ) result <- maxValue
  
  # check bounds
  for (i in 1:length(params)) {
    if ( params[[i]] < lowerBound[i]) result <- maxValue
    if ( params[[i]] > upperBound[i]) result <- maxValue
  }
  
  return(result)
}

###################### The next section is for calculating the gradient in a closed form ####################################
# the closed form derivatives are (in latex format):
# \fp{J}{K} = \sum_{t=1}^T e(t) \fp{\lambda(t)}{K} ; 
#     \fp{\lambda(t)}{K} = \frac{1 - \alpha}{\alpha - \beta - 1} \sum_{\tau=1}^{t} \hat\lambda(t-\tau) (\tau + c)^{-(1+\theta)}
# \fp{J}{\beta} = \sum_{t=1}^T e(t) \fp{\lambda(t)}{\beta} ;
#     \fp{\lambda(t)}{\beta} = \frac{1 - \alpha}{\left( \alpha - \beta - 1 \right)^2} \sum_{\tau=1}^{t} \hat\lambda(t-\tau) (\tau + c)^{-(1+\theta)}
# \fp{J}{c} = \frac{1 - \alpha}{\alpha - \beta - 1} (1+\theta) \sum_{t=1}^T e(t) \fp{\lambda(t)}{c} ;
#     \fp{\lambda(t)}{c} = \sum_{\tau=1}^{t} \hat\lambda(t-\tau) (\tau + c)^{-(2+\theta)}
# \fp{J}{\theta} = \frac{1 - \alpha}{\alpha - \beta - 1} \sum_{t=1}^T e(t) \sum_{\tau=1}^{t} \hat\lambda(t-\tau) \ln (\tau+c) (\tau + c)^{-(1+\theta)}
##
##
#' Calculates the gradient of the \lambda(t), dependent on the values of the 
#' \lambda(t) and of the gradient of \lambda(0:t-1) at previous times. This 
#' gradient is necessary in the calculation of the gradient of the error 
#' function for fitting. The formula of calculation is given in version 3 of the
#' work documents.
#' 
#' @param params a list or vector with the models parameters.
#' @param alpha parameter of the model, non variable for now.
#' @param t the time at which to calculate the gradient.
#' @param lambda_values the series of values for lambda (needed in the
#'   calculation). If NULL, then the series is generated using the parameters.
#' @param previous_grad_lambda the matrix containing the gradient values at
#'   previous moments of time. If NULL, then this function will reconstruct the
#'   matrix by recursivelly calling itself for previous moments of time.
#' @param ext_infl: the series containing S(t) - the external influence in our 
#'   model
#' @return the matrix, with t+1 rows and nvar columns, giving the values of the
#'   gradient for each variable, at each moment of time <= the indicated time.
grad_lambda <- function(params = c(gamma = 100, eta = 100, K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, mu1 = 1), t, 
                        lambda_values = NULL, previous_grad_lambda = NULL, ext_infl = NULL, alpha = 2.016) {
  
  # process parameters and external influence
  res <- .check_fix_mus_ext_infl(params = params, ext_infl = ext_infl, t = t)
  params <- res$params
  ext_infl <- res$ext_infl
  
  # there is a pathological case, in which no big calculations are necesary: t = 0
  # simply test and return what we were asked
  if (t == 0) {
    calculation <- data.frame(t = 0, gamma = 1, eta = 0, K = 0, beta = 0, c = 0, theta = 0, q = 0, row.names = NULL )
    for (i in 1:length(ext_infl)) {
      muname <- sprintf("mu%d", i)
      calculation[[muname]] <- ext_infl[[i]][1]
    }
    return(calculation)
  }
  
  # check if we were given the lambda_values
  if ( is.null(lambda_values) ) {
    lambda_values <- generate_simulated_data(params = params, time = t, ext_infl = ext_infl, alpha = alpha)
  }
  
  # check if we have the previous values of the gradient of lambda
  if ( is.null(previous_grad_lambda)) {
    # call ourselves with one less time step, to get the previous values
    previous_grad_lambda <- NULL
    for (time in 0:(t-1)) {
      previous_grad_lambda <- grad_lambda(params = params, t = time, lambda_values = lambda_values, 
                                          previous_grad_lambda = previous_grad_lambda, ext_infl = ext_infl, alpha = alpha)
    }
  }
  
  # keep only what we want
  prev_grad <- previous_grad_lambda[previous_grad_lambda$t < t,]
  
  # eliminate values which do not interest us (scope of time is above our t)
  lambdai <- lambda_values$Count[lambda_values$Time < t]
  timei <- lambda_values$Time[lambda_values$Time < t]
  
  # calculate the time difference between each lambdai and the current moment
  taui <- t - timei
  
  # calculate the rest of the formula
  # 1. time decay part
  decay <- (taui + params$c) ^ (-1-params$theta)
  decay_c <- (taui + params$c) ^ (-2-params$theta) # the decay has a slighter different form for parameter c
  
  # 2. get the previous grad lambda values and sum them up
  prev_grad <- prev_grad * decay
  prev_grad <- colSums(prev_grad, na.rm = T)[-1] #first column is Time
  
  # 3. multiply the values of real lambda
  result <-  decay * lambdai #valid for q
  result_c <- -(1 + params$theta) * sum(decay_c * lambdai) # the form for parameter c
  # for parameter theta, we have an additional term ln(taui+c)
  result_theta <- -1 * sum(result * log(taui+params$c))
  result_q <- sum(result)
  
  # calculate final values
  val_q = params$K * ((alpha - 1) / (alpha - params$beta - 1))
  grad_theta <- val_q * (result_theta + prev_grad["theta"])
  grad_c <-  val_q * (result_c + prev_grad["c"])
  grad_eta = 1 + val_q * prev_grad["eta"]
  grad_q = result_q + val_q * prev_grad["q"]
  grad_K = ((alpha - 1) / (alpha - params$beta - 1)) * grad_q
  grad_beta = (params$K * (alpha - 1)) / ((alpha - params$beta - 1)^2) * grad_q
  grad_gamma = val_q * prev_grad["gamma"]
  
  calculation <- data.frame(t = t, gamma = grad_gamma, eta = grad_eta, K = grad_K, 
                            beta = grad_beta, c = grad_c, theta = grad_theta, q = grad_q, 
                            row.names = NULL )
  # calculate and add the mus
  for (i in 1:length(ext_infl)) {
    muname <- sprintf("mu%d", i)
    calculation[[muname]] <- ext_infl[[i]][t+1] + val_q * prev_grad[muname]
  }
  
  prev_grad <- previous_grad_lambda[previous_grad_lambda$t < t,]
  calculation <- rbind(prev_grad, calculation)
  
  return(calculation)
}

# using the results of the previous function, this next function will calculate
# the closed form gradient of the error function, over the parameters
error_function_gradient <- function(params = c(gamma = 100, eta = 100, K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, mu1 = 1), real_lambda_history, 
                                    ext_infl = NULL, alpha = 2.016, return_error_vector = FALSE, disable_gradient_params = "beta",
                                    lowerBound = NULL, upperBound = NULL,
                                    alpha_regularizer = 0, param_scaling = NULL) {
  # obtain the parameters from the list of parameters. This will be coming from
  # the optimization procedure process parameters
  params <- .correct_names(params)
  
  # check if we got bounds
  if (is.null(lowerBound)) lowerBound <- rep(-Inf, times = length(params))
  if (is.null(upperBound)) upperBound <- rep(Inf, times = length(params))
  
  # check if we got params scaling
  if (is.null(param_scaling))
    param_scaling <- rep(x = 1, times = length(params))
  param_scaling <- unlist(param_scaling)
  
  # what is the maximum extent of time?
  max_time <- max(real_lambda_history$Time)
  
  # generate the articial series corresponding to the parameters
  lambda_values <- generate_simulated_data(params = params, time = max_time, ext_infl = ext_infl, alpha = alpha)
  
  # calculate the grad of the lambda, for the given parameters, using the function defined above.
  lambda_grad <- grad_lambda(params = params, t = max_time, lambda_values = lambda_values, ext_infl = ext_infl, alpha = alpha)
  
  # calculating e(t)
  error_t <-  lambda_values$Count - real_lambda_history$Count
  
  return_val <- colSums(error_t * lambda_grad)
  # we don't want the t and q parameters 
  return_val <- return_val[!(names(return_val) %in% c("t", "q"))]
  
  # calculate the regularization
  #   disable_reg <- c("theta", "c")
  disable_reg <- c(5, 6)
  #   reg <- alpha_regularizer * (as.numeric(params) / param_scaling) ## wrong, but fitting results are with this
  reg <- alpha_regularizer / param_scaling ## correct
  reg[disable_reg] <- 0
  return_val <- return_val + reg
  
  # L-BFGS-B doesn't like Inf values for this function, therefore we will replace Inf with a ridiculously large value
  maxValue <- .Machine$double.xmax - 1
  
  # check bounds
  for (i in 1:length(params)) {
    if ( params[[i]] < lowerBound[i]) return_val[i] <- -1 * maxValue
    if ( params[[i]] > upperBound[i]) return_val[i] <- maxValue
  }
  
  return_val[is.na(return_val)] <- maxValue
  return_val[is.nan(return_val)] <- maxValue
  return_val[return_val > maxValue] <- maxValue
  return_val[!is.finite(return_val)] <- maxValue
  
  if (!is.null(disable_gradient_params)) {
    #   return_val["eta"] <- 0 # this is how you exclude a variable from fitting
    for ( i in disable_gradient_params) return_val[i] <- 0
  }
  
  return(return_val)
}

# calculate the gradient of the function, using the finite differences. Do the
# calculation in parallel
parallel_gradient <- function(.params, ...) { # Now use the cluster 
  dp = cbind(rep(0,length(.params)),diag(.params * 1e-8));   
  #   Fout = parCapply(.cl, dp, function(x) fn(.params + x, ...)); # Parallel 
  Fout = apply(dp, MARGIN = 2, function(x) error_function(.params + x, ...)); # Parallel 
  return((Fout[-1]-Fout[1])/diag(dp[,-1]));                  #
}

# artificial data simulation
fit_artificial_data <- function(params = c(gamma = 100, eta = 100, K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, mu1 = 1), ext_infl = NULL,
                                alpha = 2.016, time = 100, maxit = 100, noise = FALSE, factor = 0.3, multiplicative = TRUE) {
  # process parameters
  params <- .check_fix_mus_ext_infl(params = params, ext_infl = ext_infl)$params
  
  # first set the parameters and compute various measures
  n <- get_n(K = params$K, beta = params$beta, c = params$c, theta = params$theta) # no longer defined
  q <- params$K * ((alpha - 1) / (alpha - params$beta - 1))
  pars <- data.frame(params, q = q, n = n, row.names = "initial")
  #   pars <- data.frame(gamma = params$gamma, eta = params$eta, K = params$K, beta = params$beta,  
  #                      q = q, c = params$c, theta = params$theta, n = n)
  #   rownames(pars) <- "initial"
  # generate the simulate series and fit the curve using our fitting method
  gen_series <- generate_simulated_data(params = params, 
                                        time = time, 
                                        ext_infl = ext_infl,
                                        alpha = alpha,
                                        noise = noise, 
                                        factor = factor, 
                                        multiplicative = multiplicative)$Count
  
  results_gen_series <- fit_series(data_series = gen_series, 
                                   initial_params = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 
                                   ext_infl = ext_infl, alpha = alpha, method = "nloptr" )$model
  #   #   print(gen_series)
  #   results_gen_series <- optim( par = unlist(.correct_names(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1))), ## easy to fit starting from perfect parameters, no?: unlist(params), 
  #                                fn = error_function,
  #                                gr = error_function_gradient,
  #                                real_lambda_history = gen_series, 
  #                                ext_infl = ext_infl,
  #                                control = list(maxit = maxit),
  #                                method = "BFGS")
  #   # check convergence
  #   if (results_gen_series$convergence != 0) {
  #     warning(sprintf("The fitting did not converge. The message is: '%s'. Try increasing the number of iterations using the maxit parameter!", results_gen_series$message))
  #   }
  
  # compute the measures using the fitted values and put them into the same data frame
  q_fit <- results_gen_series$par["K"] * ((alpha - 1) / (alpha - results_gen_series$par["beta"] - 1))
  n_fit <- get_n(K = results_gen_series$par["K"], beta = results_gen_series$par["beta"], c = results_gen_series$par["c"], theta = results_gen_series$par["theta"])
  pars <- rbind(pars, data.frame( as.list(results_gen_series$par), q = q_fit, n = n_fit, row.names = "fitted"))
  
  # print the results
  print(pars, digits=4)
  
  return(results_gen_series)
}

# this function calculates the value of q, based on K and beta.
get_q <- function(K = 0.024, beta = 0.5, alpha = 2.016) {
  result <- K * ((alpha - 1) / (alpha - beta - 1))
  
  return(result)
}

#' The next functions are to fit arbitrary series of real data. Internally, it 
#' constructs the real data series + timestamp, and calls the optimization 
#' procedure (with the method in the parameter).
#' 
#' @param data_series - an array of numbers which are to be fitted
#' @param init_params ( default c(gamma = 100, eta = 100, K = 0.0024, beta = 1, 
#'   c = -0.73, theta = 3.35, mu1 = 1) )- the initial parameters in the search 
#'   space
#' @param lowerBound (default NULL) - a vector of values which are the lower 
#'   bounds for the parameters. If a given parameter is lower than its lower 
#'   bound, than the error function is set to a very large value (for optim, 
#'   nloptr does it itself)
#' @param upperBound (default NULL) - same as lowerBound, just for the upper 
#'   bound.
#' @param alpha_regularizer (default 0) - metaparameter for regularizing the 
#'   values of parameters. If zero, no regularization.
#' @param param_scaling (default NULL) - vector of scaling values for each 
#'   parameter. If NULL, a dummy vector of ones will be used.
#' @params disable_optim_params (default NULL) - parameters given in this array 
#'   will not be optimized and will remain at the initial values
#' @params method (options: c("nloptr", "optim") - default "nloptr") the method
#'   used to implement the gradient descent. "nloptr" is based on the "nloptr"
#'   library, it is faster and has included bound checks.
#'   
#' @return a vector of number (the generated counts from the fitted parameters) 
#'   and the fitted model as returned by the optimization procedure.
fit_series <- function(data_series, initial_params = c(gamma = 100, eta = 10, K = 0.0024, beta = 0.1, c = 0.1, theta = 3.35, mu1 = 1), 
                       ext_infl = NULL, alpha = 2.016, disable_relaxation_kernel = F, lowerBound = NULL, upperBound = NULL,
                       alpha_regularizer = 0, param_scaling = NULL, disable_optim_params = "beta", method = "nloptr") {
  # process parameters
  initial_params <- .check_fix_mus_ext_infl(params = initial_params, ext_infl = ext_infl)$params
  
  # first compose the data structure that we need for fitting
  counts <- data.frame(Time = 1:length(data_series)-1, Counts = data_series)
  
  # if we are to disable the relaxation kernel, we set K to zero and make sure it is not updated in the gradient
  disable_gradient_params = disable_optim_params
  if (disable_relaxation_kernel) {
    initial_params["K"] <- 0
    disable_gradient_params <- c(disable_gradient_params, "K", "beta", "theta", "c")
  }
  
  ## check method
  if (! method %in% c("nloptr", "optim")) {
    warning(sprintf("Method %s not recognised. Defaulting to \"nloptr\"", method))
    method <- "nloptr"
  }
  
  ##################### DONE checking. Start fitting.
  if (method == "nloptr") {
    ## alternativelly, use the library "nloptr" lbfgs implementation (implements
    ## L-BFGS). Note that the lbfgs function in "nloptr" does the bound checking on
    ## its own, so no need to pass it to the error function hack
    fitted_model <- NULL
    require(nloptr)
    tryCatch(
      fitted_model <- lbfgs(x0 = unlist(initial_params),
                            fn = error_function, 
                            gr = error_function_gradient, 
                            lower = lowerBound, upper = upperBound,
                            control = list(xtol_rel = "1e-8", maxeval = 1000, check_derivatives = F), ## TODO: maybe lower maxeval for faster (less accurate fitting)
                            real_lambda_history = counts,  ## from here, our function parameters
                            alpha = alpha,
                            ext_infl = ext_infl,
                            disable_gradient_params = disable_gradient_params,
                            alpha_regularizer = alpha_regularizer, param_scaling = param_scaling) ,
      error = function(e) warning(paste("There was an error: ", e)))
    
    # check convergence
    if (!is.null(fitted_model) && fitted_model$convergence < 0) {
      warning(sprintf("The fitting did not converge. The message is: '%s'. Try increasing the number of iterations using the maxit parameter!", fitted_model$message))
      fitted_model <- NULL ## fallback to optim
    }
    
    ## if something went wrong, fallback to optim
    if (is.null(fitted_model)) {
      warning("\"nloptr\" did not succeed! Trying fallback optim!")
      method <- "optim"
    }
  }
  
  if (method == "optim") {
    ## next, fit the parameters of our model - using R supplied BFGS
    fitted_model <- NULL
    tryCatch(
      fitted_model <- optim( par = initial_params,
                             fn = error_function,
                             gr = error_function_gradient,
                             real_lambda_history = counts, 
                             alpha = alpha,
                             ext_infl = ext_infl,
                             disable_gradient_params = disable_gradient_params,
                             lowerBound = lowerBound, upperBound = upperBound,
                             alpha_regularizer = alpha_regularizer, param_scaling = param_scaling,
                             control = list(maxit = 5000),
                             method = "BFGS") ,
      error = function(e) warning(paste("There was an error: ", e)))
    
    # check convergence
    if (!is.null(fitted_model) & fitted_model$convergence != 0) {
      warning(sprintf("The fitting did not converge. The message is: '%s'. Try increasing the number of iterations using the maxit parameter!", fitted_model$message))
    }
  }
  
  if (is.null(fitted_model)) return(list(fitted_counts = rep(NA, times = nrow(counts), model = NA)))
  
  ## if here, all went alright.
  ## in case of nloptr lbfgs, we need to reinstate the model's parameter names.
  names(fitted_model$par) <- names(initial_params)
  
  # generate a artificial sequence with the fitted parameters
  fitted_counts <- generate_simulated_data(params = fitted_model$par, time = nrow(counts)-1, ext_infl = ext_infl, alpha = alpha)$Count
  
  # and return the fitted counts and model
  result <- list(fitted_counts = fitted_counts, model = fitted_model)
  return(result)
}

#' Internal usage method only. It checks that the parameters "mu_i" are in sync 
#' with the external influences received by the algorihms. The purpose is to 
#' have as many parameters as internal influence series. If t is provided 
#' (default NULL), than the external influences will be tailored to the length 
#' required, if needed by adding zeros.
#' 
#' @param params - list or numeric vector, containing parameters. 6 default 
#'   parameters (gamma, eta, K, beta, c, theta) and a variable number of mu 
#'   parameters (mu1, mu2 ...)
#' @param ext_infl (default NULL)- a list of series which are the external
#'   influence
#' @param t (default NULL) - the time length for which to tailor the external
#'   influence
#'   
#' @return a list containing the new parameters and the new list of external
#'   influences.
.check_fix_mus_ext_infl <- function(params, ext_infl = NULL, t = NULL, lBound = 0, uBound = 505.90) {
  # process parameters
  params <- .correct_names(params)
  
  # deal with the case of no external influence (we got a NULL parameter)
  if (is.null(ext_infl)) {
    ext_infl <- list(rep(x = 0, times = 2))
  }
  
  # by default we have only one series of external influence. Make sure that if
  # we have more, than the number of mus is identical the parameters are:
  # (gamma, eta, K, beta, c, theta, mu1, mu2, ...). 6 fixed ones + as many mus
  # as necessary
  
  # if less parameters than external influences => add parameters
  if (length(ext_infl) + 6 > length(params)) {
    # we need more mus. add with warning
    #     warning(sprintf("Insufficiant mu parameters. I got %d external series and only %d mu params. Adding the rest.", length(ext_infl), length(params)-6))
    existing_mus <- length(params) - 6
    for (i in 1:(length(ext_infl) + 6 - length(params))) {
      existing_mus <- existing_mus + 1
      varname <- sprintf("mu%d", existing_mus)
      params[[varname]] <- runif(n = 1, min = lBound, max = uBound) #initially 1
    }
  }
  
  # if more parameters than external influences => add dummy external influence
  if (length(ext_infl) + 6 < length(params)) {
    needed_series <- length(params) - 6 - length(ext_infl)
    for (i in 1:needed_series) {
      ext_infl <- c(ext_infl, list(NA))
    }
  }
  
  # if t is provided, make sure that each external influence is tailored to the desired length t+1 (time starts at 0)
  if (!is.null(t)) {
    # the external influence should be a vector of the same length as the desired period
    ext_infl <- lapply(X = ext_infl, 
                       FUN = function(x) {
                         x <- unlist(x)
                         if (length(x) > t+1) x <- x[1:(t+1)]
                         if (length(x) < t+1) x <- c(x, rep(x = NA, times = t + 1 - length(x)))
                         x[!is.finite(x)] = 0
                         
                         return(x)
                       })
  }
  
  return(list(params = params, ext_infl = ext_infl))
}

.get_n <- function(params, alpha = 2.016, mmin = 1) {
  return(get_n(K = params$K, beta = params$beta, c = params$c, theta = params$theta, alpha = alpha, mmin = mmin)) 
}

# the get_n function calculates the average number of events generated by an event
# useful to detect super-critical regimes
get_n <- function(K = 0.024, alpha = 2.016, beta = 0.5, mmin = 1, c = 0.001, theta = 0.2) {
  
  if (!is.finite(K * alpha * beta * mmin * c * theta))
    return(NA)
  
  if (beta >= (alpha - 1) ) {
    warning("The closed expression calculated by this function does NOT hold for beta >= alpha - 1")
    
    return (NA)
  }
  
  if (theta <= 0 ) {
    warning(sprintf("The closed expression calculated by this function does NOT hold for theta <= 0 (K=%.4f, beta=%.2f, theta=%.2f)", K, beta, theta))
    
    return (NA)
  }
  
  n0 = (K * (alpha - 1)) / (alpha - 1 - beta)
  int = 1 / ( theta * c^theta)
  #     if (theta > 0.001)
  #         int1 = integrate(f=function(t) 1/((t+c)^(1+theta)),lower=0,upper=Inf)$value
  #     else
  #         int1 = int
  #     
  #     if ( abs(int - int1) >= 0.00001) {
  #         warning(sprintf("Integral (val=%.4f) and closed form (val=%.4f) solution in calculating n do not give the same result! Params: theta: (%.4f", int1, int, theta))
  #     }
  
  n = n0 * int
  
  return (n) 
}

#' Calculates the time-decaying kernel corresponding to a given event/set of 
#' events, at the givent time. The mernel function is the probability that an 
#' event is spawned at the current time, considering the previous events given 
#' as parameters.
#' 
#' @param event - the event is a 2 elements list (mi, ti), magnitude and arrival
#'   time of the event.
#' @param t - time at which to calculate the the probability.
#' @param params - list of parameters (see above). Note that gamma, eta and mui
#'   are not necessary here.
kernel_fct <- function(event, t, params = list( K = 0.024, beta = 0.5, c = 0.001, theta = 0.2), mmin = 1, alpha = 2.016 ) {
  params <- .correct_names(params)
  
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]
  
  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- params$K
  
  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ params$beta
  
  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + params$c)^(1 + params$theta)
  
  val = fun_f * fun_ro * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  
  (val)
}

#' Calculate the total number of events generated by a single event, of 
#' magnitude 1 at time 0. This is done by numerical integration of the generated
#' curve and it allows, for super-critical videos, to determine the 
#' characteristic time and for sub-critical the maturity time.
#' 
#' @return a list with 4 components: endo - the endogenous reponse (the number 
#'   of events as a response to a single event of mag 1 at time 0), n - 
#'   branching factor, the characteristic time and the maturity time at the
#'   given percentage (controlled by the "omega_perc_maturity_time" param)
get_endogenous_response <- function(params = c(gamma = 1, eta = 0, K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, mu1 = 1), 
                                    alpha = 2.016, mmin = 1, max_simulation_time = 10000, omega_perc_maturity_time = 0.95) {
  params <- .correct_names(params)
  n <- .get_n(params = params, alpha = alpha, mmin = mmin)
  # remember: one event of mag 1 at time 0
  params$gamma <- 1
  params$eta <- 0
  series <- generate_simulated_data(params = params, time = max_simulation_time, alpha = alpha)$Count
  
  endo <- NULL
  tryCatch( 
    endo <- quadl(f = function(prs) { sapply(X = prs, FUN = function(x) return(series[(floor(x)+1)])) }, 
                  xa = 0, xb =  max_simulation_time),
    error = function(e) warning(paste("There was an error: ", e)))
  if (!is.numeric(endo)) endo <- NA
  
  ## if integration failed, try the default method
  if (is.na(endo)) {
    tryCatch( 
      endo <- integrate(f = function(prs) { sapply(X = prs, FUN = function(x) return(series[(floor(x)+1)])) }, 
                        lower = 0, upper = max_simulation_time, subdivisions=max(50000, max_simulation_time) )$value,
      error = function(e) warning(paste("There was an error: ", e)))
    if (!is.numeric(endo)) endo <- NA
  }
  
  ## if still bad, just use a poor approximation, not to return NA
  if (is.na(endo)) {
    endo <- sum(series, na.rm = T)
  }
  
  ## determine the characteristic time by substracting consequtive values from series
  ## in other words I want the first non-negative difference
  characteristic_time <- which(series[2:length(series)] - series[1:(length(series)-1)] > 0)[1]
  maturity_time <- NA
  
  ## maturity time only exists in sub-critical regimes (in other words, when the
  ## characteristic time does not)
  if (is.na(characteristic_time)) {
    maturity_time <- which(cumsum(series) >= omega_perc_maturity_time * endo)[1]
  }
  
  return( list(endo = endo, n = n, characteristic_time = characteristic_time, maturity_time = maturity_time, impulse_response = list(series) ))
}

#' This is an internal function required to correct the names of the parameters.
#' It is needed because the "nloptr lbfgs" loses the parameter names. Even when
#' passed an initial parameter named array, it loses the names when recomputing
#' them (for example when exploring the parameter space.) The purpose is to
#' check if the parameter array has the required names. If not, they will be
#' attributed the default names needed for this library. The parameter names is
#' assumed to be in order: 
#'      c("gamma", "eta", "K", "beta", "c", "theta", "mu1", "mu2", ...)
.correct_names <- function(params) {
  ## first make sure the params are a list
  params <- as.list(unlist(params))
  ## second make sure they have names attached
  if (is.null(names(params))) {
    par_names <- c("gamma", "eta", "K", "beta", "c", "theta")
    lng <- length(par_names)
    for (i in (lng+1):length(params)) {
      par_names <- c(par_names, sprintf("mu%d", i-lng))
    }
    names(params) <- par_names
  }
  
  return(params)
}
