## Implement double LASSO model selection on K=0,1,2 in BHM reg
rm(list = ls(all = TRUE))

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")
pacman::p_load(fixest, tidyverse, janitor, broom, useful, hdm)
theme_set(theme_bw())

if (Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
  code <- '~/Documents/GitHub/desCE/'
} else if (Sys.info()['user'] == "ux310uq-gl443t") {
  root <- "D:/Dropbox (Princeton)/projects/"  
} else if (Sys.info()['user'] == "fpalomba") {
  root <- "/Users/fpalomba/Dropbox (Princeton)/projects/"  
}
set.seed(123)
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out  <- paste0(db, "/out/draft/")

if(Sys.info()['user'] == "tombearpark") source(file.path(code, '/utils/cvFuncs.R'))

# load and clean data -----------------------------------------------------

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1 = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1, time1 = time, 
         time2   = time^2, time3 = time^3, time4 = time^4) %>% 
  rename(y = growthWDI, 
         x1 = temp1, x2 = temp2, 
         w1 = precip1, w2 = precip2)

# baseline reg ------------------------------------------------------------

reg0 <- feols(data = df, y ~ x1 + x2 + w1 + w2 |
                  iso + time1 + iso[time1] + iso[time2], 
                panel.id = c('time1', 'iso'), 
                cluster = ~iso)

etable(reg0, keep = c("x1", "x2"))

# double selection --------------------------------------------------------
# hdm does not allow user to choose which coefficients penalize and which not
# so we will use glmnet. first we check that we get the same answers when 
# lambda is the same
y <- df$y
d1 <- df$x1
d2 <- df$x2
covs.cont <- c("w1", "w2")
X <- model.matrix(y ~ w1 + w2 + as.factor(iso) + as.factor(time1) + 
                    (time1 + time2 + time3 + time4) * as.factor(iso), df)[, -1]

# follow Gelman's advice to standardize dummies https://statmodeling.stat.columbia.edu/2009/07/11/when_to_standar/
Xstd <- X
Xstd[,covs.cont] <- apply(X[,covs.cont], 2, function(x) x / (2*sd(x)))


sel.penalized <- (unlist(lapply(purrr::map(stringr::str_split(colnames(Xstd), ":"), 2), function(l) length(l))) == 1)
time.fe <- unlist(lapply(stringr::str_split(colnames(Xstd), "as.factor\\(time1\\)"), function(l) length(l) == 2)) |
           unlist(lapply(stringr::str_split(colnames(Xstd), "as.factor\\(time2\\)"), function(l) length(l) == 2)) |
           unlist(lapply(stringr::str_split(colnames(Xstd), "as.factor\\(time3\\)"), function(l) length(l) == 2)) |
           unlist(lapply(stringr::str_split(colnames(Xstd), "as.factor\\(time4\\)"), function(l) length(l) == 2))
sel.penalized[colnames(Xstd) %in% c("time1", "time2", "time3", "time4")] <- TRUE
sel.penalized[time.fe] <- TRUE


X.nopen <- Xstd[,!sel.penalized]
X.pen <- Xstd[,sel.penalized]

# run lasso for each treatment equation and for outcome equation
lasso.d1.cv <- glmnet::cv.glmnet(x=cbind(X.pen, X.nopen), y=d1, alpha=1, standardize = FALSE,
                                 penalty.factor = c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen))))

lasso.d1 <- glmnet::glmnet(x=cbind(X.pen, X.nopen), y=d1, alpha=1, lambda = c(lasso.d1.cv$lambda.1se, lasso.d1.cv$lambda.min),
                           penalty.factor = c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen))), standardize = FALSE)


lasso.d2.cv <- glmnet::cv.glmnet(x=cbind(X.pen, X.nopen), y=d2, alpha=1, standardize = FALSE,
                                 penalty.factor = c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen))))

lasso.d2 <- glmnet::glmnet(x=cbind(X.pen, X.nopen), y=d2, alpha=1, lambda = c(lasso.d2.cv$lambda.1se, lasso.d2.cv$lambda.min),
                           penalty.factor = c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen))), standardize = FALSE)


lasso.y.cv <- glmnet::cv.glmnet(x=cbind(X.pen, X.nopen), y=y, alpha=1, maxit = 1000, standardize = FALSE,
                                 penalty.factor = c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen))))

lasso.y <- glmnet::glmnet(x=cbind(X.pen, X.nopen), y=y, alpha=1, lambda = c(lasso.y.cv$lambda.1se, lasso.y.cv$lambda.min),
                           penalty.factor = c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen))), standardize = FALSE)


# store coefficients different from 0
# according to 1stdv rule 
lam.std.sel <- (lasso.d1$beta[ , 1] > 1e-04) | (lasso.d2$beta[ , 1] > 1e-04) | (lasso.y$beta[ , 1] > 1e-04) | sel.penalized == FALSE

# according to min mse 
lam.min.sel <- (lasso.d1$beta[ , 2] > 1e-04) | (lasso.d2$beta[ , 2] > 1e-04) | (lasso.y$beta[ , 2] > 1e-04) | sel.penalized == FALSE

sum(lam.std.sel)
sum(lam.min.sel)

post.lasso <- list()
post.lasso[["lam.std"]] <- lm(y ~ -1 + X[, lam.std.sel])
post.lasso[["lam.min"]] <- lm(y ~ -1 + X[, lam.min.sel])


rf0 <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
rf0 %>%  
  plot_rf_poly(facet.var = 'id')

