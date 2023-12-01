pacman::p_load(tidyverse, fixest)

# cv funcs ----------------------------------------------------------------

# Split a data set into training and test data. 
split <- function(df, id.vars, test.prop){
  
  # Get unique values of each id variable
  id.vals  <- lapply(id.vars, \(x) unique(df[[x]]))
  # Randomly sample across these unique values
  test.ids <- lapply(
    id.vals, \(x) sample(x, size = (test.prop^(1/length(id.vars))) * length(x), 
                                          replace = FALSE))
  df$test <- 1 
  for(ii in seq_along(id.vars)) {
    
    df$test <- ifelse(df[[ id.vars[[ii]] ]] %in% test.ids[[ii]], 
                          1 * df$test, 0)
  }
  df
}  

# Calculate the RMSE between two vectors
rmse <- function(predictions, targets){
  sqrt(mean((predictions - targets)^2))
}

# Loop over iteratiosn (1:K) and a list of models, evaluating their RMSE
# out of sample performance 
crossvalidate <- function(df, K, models, dependent, id.vars, test.prop){
  
  map_dfr(
    1:K, 
    function(kk){
      data <- split(df, id.vars, test.prop)
      training_set <- data[data$test != 1,]
      testing_set  <- data[data$test == 1,]
      
      map_dfr(
        models, 
        function(model){
          m         <- feols(as.formula(model), data = training_set)
          predicted <- predict(m, testing_set)
          tibble(model = model, 
                 rmse = rmse(predicted, testing_set[[dependent]]))
        }
      )
    }
  ) 
}

# Helpful wrapper
run_cv <- function(df, dep.var, FEs, id.vars, test.prop, K){
  
  models <- paste0(dep.var, FEs)
  cv.r <- crossvalidate(df, K = K, model = models, dependent = dep.var, 
                        id.vars = id.vars, test.prop = test.prop)
  
  cv.r %>% 
    group_by(model) %>% 
    summarize(rmse = mean(rmse)) %>% 
    arrange(rmse) 
}

# Choose the model with the lowest RMSE 
best.model <- function(cv.out){
  cv.out <- arrange(cv.out, rmse)
  as.formula(cv.out$model[[1]])
}

# lasso -------------------------------------------------------------------
run_double_lasso <- function(pred.vars, X.pen, X.nopen, standardize = TRUE, 
                             maxit = 1000L){
  
  X.mat    <- cbind(X.pen, X.nopen)
  dim.nopen <- ifelse(is.null(X.nopen), 0, ncol(X.nopen))
  penalties <- c(rep(1, ncol(X.pen)), rep(0, dim.nopen))
  
  map(
    names(pred.vars), 
    function(var){
      message(paste0("running for ", var))
      
      lasso.cv <- glmnet::cv.glmnet(x=X.mat, y=pred.vars[[var]], 
                                    alpha=1, standardize = standardize,
                                    penalty.factor = penalties, 
                                    maxit = maxit)
      
      glmnet::glmnet(x=X.mat, y=pred.vars[[var]], alpha=1, 
                     lambda = c(lasso.cv$lambda.min),
                     penalty.factor = penalties, standardize = standardize)      
    }
  ) %>% 
    set_names(names(pred.vars))
}


extract_coefs <- function(l.out, tol, sel.nopen){
  map_dfr(names(l.out), 
          function(xx){
            l <- l.out[[xx]]
            tibble(term = rownames(l$beta), beta = l$beta[,1], var = xx)
          }) %>% 
    mutate(
      forced = ifelse(term %in% sel.nopen, 1, 0), 
      chosen = ifelse(abs(beta) > tol, 1, 0), 
      selected = ifelse(chosen == 1 | forced == 1, 1, 0))
}

get_selected_controls <- function(l.out, tol, sel.nopen){
  extract_coefs(l.out=l.out, tol=tol, sel.nopen=sel.nopen) %>% 
    dplyr::filter(selected == 1) %>% 
    dplyr::pull(term) %>% 
    unique()
}

get_post_selection_regdf <- function(df, outcome, treatments, 
                                     X, 
                                     selected_controls, 
                                     treatment.names=NULL){
  selected_controls <- selected_controls[!selected_controls %in% treatments]
  
  reg.df <- bind_cols(y = df[[outcome]], 
                          dplyr::select(df, all_of(treatments)), 
                          set_names(as_tibble (X[, selected_controls]), 
                                    selected_controls)
                      )  %>% 
    janitor::clean_names()
  
  if(!is.null(treatment.names)){
    names(reg.df)[names(reg.df) %in% treatments] <- treatment.names
  }
  reg.df
}

.f <- function(...) as.formula(paste0(...))
