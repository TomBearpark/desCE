simulModel <- function(N = 1000, g = 0.5, p = c(0.5, 0.5), tau = 1,
                       s2n = 1) {
  
  ######################################################
  # N: sample size
  # g: fraction in group 1
  # p: pscore in both groups c(p_0, p_1)
  # tau: treatment shift
  # s2n: we normalize the variances of the groups to be 
  #      1, as such s2n parameterizes the mean.
  ######################################################
  
  # grouping structure: 100*g% gets in one group
  Ng1 <- ceiling(g*N)
  Ng0 <- N - Ng1
  G <- c(rep(0, Ng0), rep(1, Ng1))
  
  # assign treatment depending on group
  D <- c(stats::rbinom(Ng0, 1, p[1]), 
         stats::rbinom(Ng1, 1, p[2]))
  
  Y0 <- c(rnorm(Ng0, 0, sigma),
          rnorm(Ng1, s2n, sigma))
  Y1 <- Y0 + tau
  Y1het <- Y0 + tau * (1 + G)
  
  Y <- D * Y1 + (1 - D) * Y0
  Yhet <- D * Y1het + (1 - D) * Y0
  
  df <- data.frame(Y=Y, Yhet=Yhet, D=D, G=G, Y1=Y1, Y1het=Y1het, Y0=Y0)
  return(df)
}

regEst <- function(df, tEff = "homosk") {
  
  if (tEff == "homosk") {
    b <- c(fixest::feols(Y ~ 1 + D, data = df, vcov="hetero")$coefficients[[2]],
           fixest::feols(Y ~ 1 + D + G, data = df, vcov="hetero")$coefficients[[2]],
           fixest::feols(Y ~ 1 + D + G*D, data = df, vcov="hetero")$coefficients[[2]],
           fixest::feols(Y ~ 1 + D + G + G*D, data = df, vcov="hetero")$coefficients[[2]])
  }

  if (tEff == "hetero") {
    res1 <- fixest::feols(Yhet ~ 1 + D, data = df, vcov="hetero")
    res2 <- fixest::feols(Yhet ~ 1 + D + G, data = df, vcov="hetero")
    res3 <- fixest::feols(Yhet ~ 1 + D + G*D, data = df, vcov="hetero")
    res4 <- fixest::feols(Yhet ~ 1 + D + G + G*D, data = df, vcov="hetero")
    
    b <- c(res1$coefficients[[2]], res1$coefficients[[2]],
           res2$coefficients[[2]], res2$coefficients[[2]] + res2$coefficients[[3]],
           res3$coefficients[[2]], res3$coefficients[[2]] + res3$coefficients[[3]],
           res4$coefficients[[2]], res4$coefficients[[2]] + res4$coefficients[[4]])
    
  }
   
   return(b)
}

gfe <- function(df, nGroups, out.var, covs.var, treat.var=NULL,
                tol = 1e-08, nInitGuesses = 10, typeGuess = "normal", 
                nCores = 1, seed = 8894) {
  
  # consider adding error checking here
  
  # prepare data (outcome, covariates, treatment)
  y <- df[[out.var]]
  X <- as.matrix(df[covs.var])
  #D <- as.matrix(df[treat.var])
  ncovs <- ncol(X)
  
  # start from different initial guesses
  if (nCores == 1) {
    set.seed(seed)
    res <- lapply(c(1:nInitGuesses), 
                  function(i) algorithmRun(y, X, ncovs, nGroups, typeGuess, tol, df))
  } else {
    res <- parallel::mclapply(c(1:nInitGuesses),
                              function(i) algorithmRun(y, X, ncovs, nGroups, typeGuess, tol, df),
                              mc.cores = nCores, mc.set.seed = seed)
  }

  # pick result with lowest loss function
  bestGuess <- which.min(unlist(lapply(res, "[[", "lossFinal")))
  out <- res[[bestGuess]]

  return(list(theta=out$theta, alpha=out$alpha, g=out$g, nIters = out$nIters,
              tolFinal = out$tolFinal, algPaths = out$algPaths))
}

algorithmRun <- function(y, X, ncovs, nGroups, typeGuess, tol, df) {
  
  guess <- drawInitialGuess(y, ncovs, nGroups, typeGuess)
  alpha0 <- guess$alpha
  theta0 <- guess$theta
  
  eps <- 1e08           # initial tolerance
  epsList <- c()        # store tolerance
  thetaList <- list()   # store thetas
  alphaList <- list()   # store alphas
  groupList <- list()   # store groups
  s <- 0
  
  while (eps >= tol) {
    
    # Step 1: group observations 
    g <- groupingGet(y, X, theta0, alpha0)
    
    # Step 2: optimize theta and alpha with RSS loss    
    df$g <- g
    modelEst <- fixest::feols(Y ~ X | g, data=df)
    theta1 <- modelEst$coefficients
    alpha1 <- as.matrix(fixest::fixef(modelEst)$g)

    # compute loss function and store outcomes
    eps <- max(abs(c(theta0 - theta1, alpha0 - alpha1)))
    theta0 <- as.matrix(theta1)
    alpha0 <- alpha1
    epsList <- c(epsList, eps)
    
    s <- s + 1
    thetaList[[s]] <- theta1
    alphaList[[s]] <- alpha1
    groupList[[s]] <- g
  }
  
  loss <- sum(modelEst$residuals^2)

  algPaths <- list(eps = epsList, theta = thetaList, alpha = alphaList, g = groupList)
  
  return(list(theta=theta1, alpha=alpha1, g=g, nIters = s,
              tolFinal = tol, lossFinal = loss, algPaths = algPaths))
}

drawInitialGuess <- function(y, ncovs, nGroups, type = "norm") {
  
  # initial guesses
  sdY <- sd(y)
  if (type == "unif") {
    theta0 <- as.matrix(stats::runif(ncovs, min=-2*sdY, max=2*sdY))
    alpha0 <- as.matrix(stats::runif(nGroups, min=-2*sdY, max=2*sdY))
  } else {
    theta0 <- as.matrix(stats::rnorm(ncovs, 0, sdY))
    alpha0 <- as.matrix(stats::rnorm(nGroups, 0, 1))
  }
  return(list(theta=theta0, alpha=alpha0))
  
}

groupingGet <- function(y, X, theta, alpha) {
  losses <- apply(alpha, 1, function(a) lossStep1get(a, y, X, theta))
  return(apply(losses, 1, which.min))
}

lossStep1get <- function(a, y, X, theta) {
  return((y - X %*% theta - a)^2)
}

