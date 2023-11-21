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
