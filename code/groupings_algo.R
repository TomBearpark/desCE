library(tidyverse)
library(fixest)
set.seed(1)

N     <- 10000
beta  <- 1
delta <- 5

df <- tibble(i = 1:N, 
             
             # Treatment effect: can make heterogeneous later 
             beta = beta, 
             
             # Randomly assign covariate
             G = sample(c(0, 1), N, replace = TRUE), 
             
             # Assign treatment, random conditional on. G
             D = stats::rbinom(N, 1, 1 / (G + 2)), 
             
             # Calculate the potential outcomes 
             y0 = rnorm(N, mean = delta * G), 
             y1 = y0 + beta, 
             
             # Calculate the observed outcomes. 
             y = D * y1 + (1-D) * y0
             )

# Check what we would recover if we had full information 
feols(y ~ D, df, vcov = "hetero")
feols(y ~ D | G, df, vcov = "hetero")
feols(y ~ D + G, df, vcov = "hetero")
feols(y ~ -1 + D + G + D*G, df, vcov = "hetero")

# What we would observe in the real world... 
df.obs <- df %>% select(i, y, D)

# B&M ---------------------------------------------------------------------

run_gfe <- function(beta.guess, 
                    alpha.guess, 
                    df.obs, 
                    tol = 0.000000001, max.iter = 1000, verbose = FALSE){
  
  # Matrices to store results in 
  beta.mat  <- matrix(nrow = max.iter, ncol = 1)
  alpha.mat <- matrix(nrow = max.iter, ncol = 2)
  
  beta.mat[1, ]  <- beta.guess
  alpha.mat[1, ] <- alpha.guess
  
  div <- tol + 1
  s   <-  0
  
  while(abs(div) > tol){
    
    if(verbose) print(div)
    
    s <- s + 1
    if(s > max.iter) stop("failed")
    
    beta <- beta.mat[s, ]
    alpha <- alpha.mat[s, ]
    
    # 2. Assignment step
    df.obs$e1 <- (df.obs$y - (beta[1]*df.obs$D + alpha[1]))^2
    df.obs$e2 <- (df.obs$y - (beta[1]*df.obs$D + alpha[2]))^2
    
    df.obs$g <- ifelse(df.obs$e1 < df.obs$e2, 0, 1)
    
    # 3. Update step 
    m <- feols(y ~ D | g, data = df.obs)
    
    # Store new coefficients
    beta.mat[s+1, ]  <- coef(m)[1]
    alpha.mat[s+1, ] <- drop(fixest::fixef(m)$g)
    
    # Check how close we are to convergence 
    div <- max(c(beta.mat[s,]-beta.mat[s+1,], alpha.mat[s,]-alpha.mat[s+1,]))
  }
  
  beta.mat <- beta.mat[1:s+1, ]
  
  return(list(div = div, df = df.obs, beta.mat = beta.mat))
}

guess_params <- function(beta.params, alpha.params){
  
  return(
    list(
      beta = rnorm(length(beta.params$mean), 
                         beta.params$mean, beta.params$sd), 
      alpha = rnorm(length(alpha.params$mean), 
                          alpha.params$mean, alpha.params$sd)
    )
  )
}

run_guesses_gfe <- function(beta.params, alpha.params, 
                            num.guesses = 100){
  
  attempts <- map(
    1:num.guesses, 
    function(ii){
      message(paste0('guess number ', ii))
      guess <- 
        guess_params(beta.params = beta.params, alpha.params = alpha.params)
      out <- run_gfe(guess$beta, guess$alpha, df.obs)
      out
    }
  )
  
  attempts
}

select_best_attempt <- function(attempts){
  # Select either the minimum div attempt, or the first div that attains the min
  divs <- map_dbl(seq_along(attempts), \(x) attempts[[x]]$div)
  attempts[[which(divs==min(divs))[1]]]
}


# Code testing
beta.params  <- list(mean = 2, sd = 3)
alpha.params <- list(mean = c(0, 0), sd = c(3, 3))

guess <- guess_params(beta.params, alpha.params)
out   <- run_gfe(guess$beta, guess$alpha, df.obs)
out$div

# Overall run
attempts <- run_guesses_gfe(beta.params, alpha.params)

map_dbl(1:100, 
        \(x) attempts[[x]]$beta.mat[length(attempts[[x]]$beta.mat)]) %>% 
  density %>% 
  plot

out <- select_best_attempt(attempts)

feols(y ~ D , data = out$df)
feols(y ~ D | g, data = out$df)


