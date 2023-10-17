## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(Hmisc, fixest, tidyverse, janitor, broom, marginaleffects, 
               useful)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out <- paste0(db, "/out/presentation/")

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1,
         time2   = time^2, time3 = time^3, time4 = time^4,
         year2 = year^2) %>% 
  rename(y = growthWDI, country = countryname) 

df.sim <- df %>% 
  select(time1 = time, iso, y, x = temp1, time2 = time2, year, year2, 
         continent, 
         gdpCAP_wdi) %>% 
  mutate(y = 100 * y) %>% 
  # Add first differences
  group_by(iso) %>% 
  mutate(fd.y = y - dplyr::lag(y), 
         fd.x = x - dplyr::lag(x), 
         l.x   = dplyr::lag(x), 
         x2 = x^2, x1 = x) %>% 
  ungroup()  %>% 
  mutate(iso = as.factor(iso), 
         timeCont = paste0(time1,continent))

N  <- nrow(df.sim)
Ni <- length(unique(df.sim$iso))

reg0 <- feols(data = df.sim , y ~  x1 + x2 | iso + time1 + iso[time1], 
              panel.id = c('time1', 'iso'), cluster = ~iso)
rf0 <- predict_poly(reg0, "x", 0, 30, 14) 
plot_rf_poly(rf0)

# funcs --------------------------------------------------------------------

# Add folds, ensuring we can estimate all group covariates 
split <- function(df.sim, id.vars, test.prop){
  
  id.vals  <- lapply(id.vars, \(x) unique(df.sim[[x]]))
  test.ids <- lapply(id.vals, 
                     \(x) sample(x, 
                                 size = sqrt(test.prop) * length(x), 
                                 replace = FALSE))
  df.sim$test <- 1 
  for(ii in seq_along(id.vars)) {
    
    df.sim$test <- ifelse(df.sim[[ id.vars[[ii]] ]] %in% test.ids[[ii]], 
                          1 * df.sim$test, 0)
  }
  df.sim
}  

rmse <- function(predictions, targets){
  sqrt(mean((predictions - targets)^2))
}

crossvalidate <- function(df.sim, K, models, dependent, id.vars, test.prop){
  
  map_dfr(
    1:K, 
    function(kk){
      data <- split(df.sim, id.vars, test.prop)
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

run_cv <- function(dep.var, FEs, id.vars, test.prop, K){

  models <- paste0(dep.var, FEs)
  cv.r <- crossvalidate(df.sim, K = K, model = models, dependent = dep.var, 
                        id.vars = id.vars, test.prop = test.prop)
  
  cv.r %>% 
    group_by(model) %>% 
    summarize(rmse = mean(rmse)) %>% 
    arrange(rmse) 
}

best.model <- function(cv.out){
  cv.out <- arrange(cv.out, rmse)
  as.formula(cv.out$model[[1]])
}

# run ---------------------------------------------------------------------


# Step 1: estimate \hat f(x), using CV

# Step 1.1: choose best FEs

FEs <- list(' ~ 1 | iso',
            ' ~ 1 | time1 + iso',
            ' ~ 1 | time1 + iso + iso[time1]',
            ' ~ 1 | time1 + iso + iso[time1] + iso[time2]')

id.vars <- c("iso", "time1")
test.prop <- .2

results.x <- run_cv("x", FEs, id.vars, test.prop, K = 100)
results.y <- run_cv("y", FEs, id.vars, test.prop, K = 100)

m.x <- feols(best.model(results.x), df.sim, combine.quick = FALSE)
m.y <- feols(best.model(results.y), df.sim, combine.quick = FALSE)
etable(m.x, m.y, fitstat = c("aic", "bic", "r2"))

# Step 1.2, estimate \hat f(x)
df.sim$hat.x1   <- predict(m.x, df.sim)
df.sim$tilde.x1 <- residuals(m.x)

df.sim$hat.y1   <- predict(m.y, df.sim)
df.sim$tilde.y1 <- residuals(m.y)

# illustrate what happened
df.sim %>% 
  filter(iso %in% c("USA", "CHN")) %>% ggplot(aes(x = year)) + 
  geom_line(aes(y = hat.x1)) + geom_line(aes(y = x), color = 'red') + 
  facet_wrap(~iso)

# Step 2: estimate the relationship between x and y
m2 <- feols(tilde.y1 ~ tilde.x1, df.sim, vcov = 'hetero')
m2


# heterogeneity -----------------------------------------------------------

## non-linear x ----------------------------------------------------------

gg <- 10

df.sim <- df.sim %>% mutate(xbin = Hmisc::cut2(x, g = gg, levels.mean = T))
cuts   <- cut2(df.sim$x, g = gg, onlycuts = T)

m.nlx <- feols(tilde.y1 ~ i(xbin, tilde.x1) + xbin, df.sim, vcov = 'hetero')
m.nlx %>% coefplot(keep = 'tilde')

# Add stuff so we can plot
pdf.nlx <- m.nlx %>% tidy() %>% 
  filter(str_detect(term, "tilde")) %>% 
  mutate(center = as.numeric(str_remove_all(term, "xbin::|:tilde.x1| "))) %>% 
  mutate(max = cuts[-1], min = cuts[-length(cuts)])

rf0$r <- NA
rf0$r[1] <- 0
for(tt in 2:29){
  estimate <- pdf.nlx %>% 
    mutate(tt = tt) %>% filter(between(tt, min, max)) %>% 
    pull(estimate) %>% as.numeric()
  
  rf0$r[tt] <- rf0$r[tt-1] + 1 * estimate
}
rf0$r <- rf0$r - rf0$r[rf0$temp == 14]

rects <- df.sim %>% 
  group_by(xbin) %>% summarize(minX = min(x), 
                               maxX = max(x), 
                               meanGDPPC = log(mean(gdpCAP_wdi)))

ggplot(rf0) + 
  geom_line(aes(x = temp, y = r), color = 'red') + 
  geom_line(aes(x = temp, y = response), color = 'blue') + 
  geom_hline(yintercept = 0) + 
  geom_rect(data = rects, aes(ymin=-Inf, ymax=Inf, xmin=minX,
                            xmax=maxX, fill=meanGDPPC), 
            alpha =0.5) +
  scale_fill_viridis_c(direction = -1)
  

# clustering on y ---------------------------------------------------------
df.sim$ky <- kmeans(df.sim$tilde.x1, centers = 5)$cluster

feols(tilde.y1 ~ i(ky, tilde.x1), df.sim, vcov = 'hetero') %>% 
  coefplot(keep = 'tilde')

df.sim %>% group_by(ky) %>% dplyr::summarize(temp = mean(x))

