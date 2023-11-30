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
y  <- df$y
d1 <- df$x1
d2 <- df$x2
covs.cont <- c("w1", "w2")
df <- df %>% 
  mutate(iso = as.factor(iso), yearFE = as.factor(time1))

X <- model.matrix(y ~ w1 + w2 + iso + yearFE + 
                    (time1 + time2 + time3 + time4) * iso, df)[, -1]

# follow Gelman's advice to standardize dummies https://statmodeling.stat.columbia.edu/2009/07/11/when_to_standar/
Xstd <- X
Xstd[,covs.cont] <- apply(X[,covs.cont], 2, function(x) x / (2*sd(x)))

# We only want to penalise time trend related stuff
sel.penalized <- str_subset(colnames(Xstd), "time")
sel.nopen     <- str_subset(colnames(Xstd), "time", negate = TRUE)

stopifnot(length(sel.penalized) + length(sel.nopen) == dim(X)[2])

X.pen   <- Xstd[, sel.penalized]
X.nopen <- Xstd[, sel.nopen]

X.mat <- cbind(X.pen, X.nopen)
penalties <- c(rep(1, ncol(X.pen)), rep(0, ncol(X.nopen)))

# run lasso for each treatment equation and for outcome equation
# Temperature
lasso.d1.cv <- glmnet::cv.glmnet(x=X.mat, y=d1, alpha=1, standardize = FALSE,
                                 penalty.factor = penalties)

lasso.d1 <- glmnet::glmnet(x=X.mat, y=d1, alpha=1, 
                           lambda = c(lasso.d1.cv$lambda.1se, lasso.d1.cv$lambda.min),
                           penalty.factor = penalties, standardize = FALSE)

# Temperature square
lasso.d2.cv <- glmnet::cv.glmnet(x=X.mat, y=d2, alpha=1, standardize = FALSE,
                                 penalty.factor = penalties)
                                 
lasso.d2 <- glmnet::glmnet(x=X.mat, y=d2, alpha=1, 
                           lambda = c(lasso.d2.cv$lambda.1se, lasso.d2.cv$lambda.min),
                           penalty.factor = penalties, standardize = FALSE)

# Outcome
lasso.y.cv <- glmnet::cv.glmnet(x=X.mat, y=y, alpha=1, maxit = 5000, standardize = FALSE,
                                penalty.factor = penalties)

lasso.y <- glmnet::glmnet(x=X.mat, y=y, alpha=1, 
                          lambda = c(lasso.y.cv$lambda.1se, lasso.y.cv$lambda.min),
                          penalty.factor = penalties, standardize = FALSE)

# store coefficients different from 0
# according to 1stdv rule 
lam.std.sel <- (lasso.d1$beta[ , 1] > 1e-04) | (lasso.d2$beta[ , 1] > 1e-04) | (lasso.y$beta[ , 1] > 1e-04) | sel.penalized == FALSE

# according to min mse 
lam.min.sel <- (lasso.d1$beta[ , 2] > 1e-04) | (lasso.d2$beta[ , 2] > 1e-04) | (lasso.y$beta[ , 2] > 1e-04) | sel.penalized == FALSE

sum(lam.std.sel)
sum(lam.min.sel)


# create response function ------------------------------------------------
reg.df <- bind_cols(y = df$y, 
                    xtreat1 = df$x1, xtreat2 = df$x2, 
                    as_tibble(X[, lam.std.sel])
)  %>% 
  janitor::clean_names()

# Check i get whats going on. 
var.df <- tibble(var = names(lam.std.sel), include = lam.std.sel)
var.df %>% pull(var)


# Total number of potential parameters
# Ni -1 (ref) * 3 (unit FE, 2& time controls), Nt-1, 2 precip
max.K <- 4
n.unitFE   <- (length(unique(df$iso)))-1
n.timeFEs  <- (length(unique(df$time1)))-1
n.trends   <- max.k * length(unique(df$iso))
n.controls <- length(c("w1", "w2"))
n.unitFE + n.timeFEs + n.trends + n.controls

# Which do we keep?
var.df <- var.df %>% 
  mutate(fe= str_replace_all(var, "as.factor\\(time", "timeFE"), 
         fe = str_replace_all(fe, "1\\)", "")) %>% 
  mutate(fe = str_remove_all(fe, "as.factor\\(iso\\)")) %>% 
  mutate(fe.type = case_when(
    str_detect(fe, "FE") ~ "time", 
    str_detect(fe, "time1") ~ "trend 1", 
    str_detect(fe, "time2") ~ "trend 2",
    str_detect(fe, "time3") ~ "trend 3",
    str_detect(fe, "time4") ~ "trend 4",
    str_detect(fe, "w") ~ 'precip', 
    !str_detect(fe, "time") ~ "unit"
  ))

var.df %>% 
  group_by(fe.type) %>% 
  tally()
var.df %>% 
  filter(include) %>% 
  group_by(fe.type) %>% 
  tally()

view(var.df)



names(reg.df) <- str_remove_all(names(reg.df), "as_factor_")

vars <- names(reg.df)[-1] %>% 
  paste0(collapse = "+")

m.dml <- feols(as.formula(paste0("y ~ ", vars)), bind_cols(reg.df, iso= df$iso), 
               cluster = 'iso')

rf.dml <- predict_poly(m.dml, "xtreat", 0, 30, 14, ci_level = 95, id.col = "DML") 

rf.dml %>%  
  plot_rf_poly(facet.var = 'id')



rf0 <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
rf0 %>%  
  plot_rf_poly(facet.var = 'id')

