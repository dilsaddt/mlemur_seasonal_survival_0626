library(nimble)

# vectorized version of dbinom
dbinom_vector <- nimbleFunction(
  run = function( x = double(1)
                  , size = double(1)
                  , prob = double(1) 
                  , log = integer(0, default = 0)
  ) {
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = prob, size = size, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dbinom_vector", dbinom_vector, envir = .GlobalEnv)

# vectorized version of rbinom
rbinom_vector <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , size = double(1)
                  , prob = double(1)
  ) {
    returnType(double(1))
    return(rbinom(length(size), prob = prob, size = size))
  })
assign("rbinom_vector", rbinom_vector, envir = .GlobalEnv)

# vectorized version of dbinom with vectorized input parameter size and scalar input parameter prob
dbinom_vector_inMixed <- nimbleFunction(
  run = function( x = double(1)
                  , size = double(1)
                  , prob = double(0) 
                  , log = integer(0, default = 0)
  ) {
    Prob <- nimNumeric(length = length(size), value = prob)
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = Prob, size = size, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dbinom_vector_inMixed", dbinom_vector_inMixed, envir = .GlobalEnv)

# vectorized version of rbinom with vectorized input parameter size and scalar input parameter prob
rbinom_vector_inMixed <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , size = double(1)
                  , prob = double(0)
  ) {
    Prob <- nimNumeric(length = length(size), value = prob)
    returnType(double(1))
    return(rbinom(length(size), prob = Prob, size = size))
  })
assign("rbinom_vector_inMixed", rbinom_vector_inMixed, envir = .GlobalEnv)


# vectorized version of dbern
dbern_vector <- nimbleFunction(
  run = function( x = double(1)
                  , prob = double(1) 
                  , log = integer(0, default = 0)
  ) {
    size <- nimNumeric(length = length(prob), value = 1)
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = prob, size = size, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dbern_vector", dbern_vector, envir = .GlobalEnv)

# vectorized version of rbern
rbern_vector <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , prob = double(1)
  ) {
    size <- nimNumeric(length = length(prob), value = 1)
    returnType(double(1))
    return(rbinom(length(prob), prob = prob, size = size))
  })
assign("rbern_vector", rbern_vector, envir = .GlobalEnv)

# vectorized version of dbern with scalar input parameter prob that will be repeated 'length'-times
dbern_vector_inScalar <- nimbleFunction(
  run = function( x = double(1)
                  , prob = double(0)
                  , length = double(0)
                  , log = integer(0, default = 0)
  ) {
    
    size <- nimNumeric(length = length, value = 1)
    Prob <- nimNumeric(length = length, value = prob)
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = prob, size = size, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dbern_vector_inScalar", dbern_vector_inScalar, envir = .GlobalEnv)

# vectorized version of rbern with scalar input parameter prob that will be repeated 'length'-times
rbern_vector_inScalar <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , prob = double(0)
                  , length = double(0)
  ) {
    size <- nimNumeric(length = length, value = 1)
    Prob <- nimNumeric(length = length, value = prob)
    returnType(double(1))
    return(rbinom(length(Prob), prob = Prob, size = size))
  })
assign("rbern_vector_inScalar", rbern_vector_inScalar, envir = .GlobalEnv)


# vectorized version of dunif with scalar input parameters min and max that will be repeated 'length'-times
dunif_vector_inScalar <- nimbleFunction(
  run = function( x = double(1)
                  , min = double(0)
                  , max = double(0)
                  , length = double(0)
                  , log = integer(0, default = 0)
                  
  ) {
    
    Min <- nimNumeric(length = length, value = min)
    Max <- nimNumeric(length = length, value = max)
    returnType(double(0))
    logProb <- sum(dunif(x, min = Min, max = Max, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dunif_vector_inScalar", dunif_vector_inScalar, envir = .GlobalEnv)

# vectorized version of runif with scalar input parameters min and max that will be repeated 'length'-times
runif_vector_inScalar <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , min = double(0)
                  , max = double(0)
                  , length = double(0)
  ) {
    Min <- nimNumeric(length = length, value = min)
    Max <- nimNumeric(length = length, value = max)
    returnType(double(1))
    return(runif(length(Min), min = Min, max = Max))
  })
assign("runif_vector_inScalar", runif_vector_inScalar, envir = .GlobalEnv)


# vectorized version of dunif with min and max being vectors
dunif_vector <- nimbleFunction(
  run = function( x = double(1)
                  , min = double(1)
                  , max = double(1)
                  , log = integer(0, default = 0)
                  
  ) {
    returnType(double(0))
    logProb <- sum(dunif(x, min = min, max = max, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dunif_vector", dunif_vector, envir = .GlobalEnv)

# vectorized version of runif with min and max being vectors
runif_vector <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , min = double(1)
                  , max = double(1)
  ) {
    returnType(double(1))
    return(runif(length(min), min = min, max = max))
  })
assign("runif_vector", runif_vector, envir = .GlobalEnv)


# vectorized version of dpois with scalar input parameter lambda that will be repeated 'length'-times
dpois_vector_inScalar <- nimbleFunction(
  run = function( x = double(1)
                  , lambda = double(0)
                  , length = double(0)
                  , log = integer(0, default = 0)
  ) {
    Lambda <- nimNumeric(length = length, value = lambda)
    returnType(double(0))
    logProb <- sum(dpois(x, lambda = Lambda, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dpois_vector_inScalar", dpois_vector_inScalar, envir = .GlobalEnv)

# vectorized version of rpois with scalar input parameter lambda that will be repeated 'length'-times
rpois_vector_inScalar <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , lambda = double(0)
                  , length = double(0)
  ) {
    Lambda <- nimNumeric(length = length, value = lambda)
    returnType(double(1))
    return(rpois(length(Lambda), lambda = Lambda))
  })
assign("rpois_vector_inScalar", rpois_vector_inScalar, envir = .GlobalEnv)



# vectorized version of dgamma with scalar input parameters shape and rate that will be repeated 'length'-times
dgamma_vector_inScalar <- nimbleFunction(
  run = function( x = double(1)
                  , shape = double(0)
                  , rate = double(0) 
                  , length = double(0)
                  , log = integer(0, default = 0)
                  
  ) {
    
    Shape <- nimNumeric(length = length, value = shape)
    Rate <- nimNumeric(length = length, value = rate)
    returnType(double(0))
    logProb <- sum(dgamma(x, shape = Shape, rate = Rate, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dgamma_vector_inScalar", dgamma_vector_inScalar, envir = .GlobalEnv)

# vectorized version of rgamma with scalar input parameters shape and rate that will be repeated 'length'-times
rgamma_vector_inScalar <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , shape = double(0)
                  , rate = double(0)
                  , length = double(0)
  ) {
    Shape <- nimNumeric(length = length, value = shape)
    Rate <- nimNumeric(length = length, value = rate)
    returnType(double(1))
    return(rgamma(length(Shape), shape = Shape, rate = Rate))
  })
assign("rgamma_vector_inScalar", rgamma_vector_inScalar, envir = .GlobalEnv)


# vectorized version of dnorm
dnorm_vector <- nimbleFunction(
  run = function( x = double(1)
                  , mean = double(1)
                  , sd = double(1) 
                  , log = integer(0, default = 0)
  ) {
    returnType(double(0))
    logProb <- sum(dnorm(x, mean = mean, sd = sd, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dnorm_vector", dnorm_vector, envir = .GlobalEnv)

# vectorized version of rnorm
rnorm_vector <- nimbleFunction(
  run = function(n = integer(0, default = 1)
                 , mean = double(1)
                 , sd = double(1)
  ) {
    returnType(double(1))
    return(rnorm(length(mean), mean = mean, sd = sd))
  })
assign("rnorm_vector", rnorm_vector, envir = .GlobalEnv)


# vectorized version of dnorm with scalar input parameters mean and sd that will be repeated 'length'-times
dnorm_vector_inScalar <- nimbleFunction(
  run = function( x = double(1)
                  , mean = double(0)
                  , sd = double(0)
                  , length = double(0)
                  , log = integer(0, default = 0)
                  
  ) {
    
    Mean <- nimNumeric(length = length, value = mean)
    SD <- nimNumeric(length = length, value = sd)
    returnType(double(0))
    logProb <- sum(dnorm(x, mean = Mean, sd = SD, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })
assign("dnorm_vector_inScalar", dnorm_vector_inScalar, envir = .GlobalEnv)

# vectorized version of rnorm with scalar input parameters mean and sd that will be repeated 'length'-times
rnorm_vector_inScalar <- nimbleFunction(
  run = function( n = integer(0, default = 1)
                  , mean = double(0)
                  , sd = double(0)
                  , length = double(0)
  ) {
    Mean <- nimNumeric(length = length, value = mean)
    SD <- nimNumeric(length = length, value = sd)
    returnType(double(1))
    return(rnorm(length(Mean), mean = Mean, sd = SD))
  })
assign("rnorm_vector_inScalar", rnorm_vector_inScalar, envir = .GlobalEnv)

# function to calculate conditional entry probabilities:
nimF.nu <- nimbleFunction(
  run = function(beta.b = double(1), n.occ = integer(0)){
    # define dimensions of vectors
    b <- nimNumeric(length = n.occ, value = NA) # vector of entry probabilities
    nu <- nimNumeric(length = n.occ, value = NA) # vector of conditional entry probabilities
    
    # first occasion
    nu[1] <- beta.b[1]/sum(beta.b[1:n.occ])
    b[1] <- nu[1]
    
    # all subsequent occasions
    for(t in 2:n.occ){
      b[t] <- beta.b[t]/sum(beta.b[1:n.occ])
      nu[t] <- b[t]/(1-sum(b[1:(t-1)]))
    }
    
    returnType(double(1)) # we return a numeric vector
    return(nu)
  }
)
assign("nimF.nu", nimF.nu, envir = .GlobalEnv)


# function to estimate occasion t at which the focal individual entered?
nimF.entered <- nimbleFunction(
  run = function(n.occ = integer(0)
                 , u = double(1) # vector containing a 1 at the occasions during which the focal individual was alive
  ){
    entered <- nimNumeric(value = 0, length = n.occ)
    if(sum(u)>0){
      entered[min(which(u[1:n.occ]==1))] <- 1
    }
    
    returnType(double(1)) # we return an numeric vector
    return(entered)
  }
)
assign("nimF.entered", nimF.entered, envir = .GlobalEnv)