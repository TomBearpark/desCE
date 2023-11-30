## Implement double LASSO model selection on K=0,1,2 in BHM reg
rm(list = ls(all = TRUE))

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")
pacman::p_load(fixest, tidyverse, janitor, broom, useful, glmnet)
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
                  iso + time1, 
                panel.id = c('time1', 'iso'), 
                cluster = ~iso)

reg1 <- feols(data = df, y ~ x1 + x2 + w1 + w2 |
                iso + time1 + iso[time1], 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

reg2 <- feols(data = df, y ~ x1 + x2 + w1 + w2 |
                iso + time1 + iso[time1] + iso[time2], 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

etable(reg0,reg1,reg2,keep = c("x1", "x2"))

# double selection --------------------------------------------------------
# hdm does not allow user to choose which coefficients penalize and which not
# so we will use glmnet. first we check that we get the same answers when 
# lambda is the same
y  <- df$y
d1 <- df$x1
d2 <- df$x2

df <- df %>% 
  mutate(iso = as.factor(iso), yearFE = as.factor(time1))

X <- model.matrix(y ~ w1 + w2 + iso + yearFE + 
                    (time1 + time2 + time3 + time4) * iso, df)[, -1]

covs.to.std <- c(c("w1", "w2"), str_subset(colnames(X), "time"))

# follow Gelman's advice to standardize dummies https://statmodeling.stat.columbia.edu/2009/07/11/when_to_standar/
Xstd <- X
# ACtually no we dont
# Xstd[,covs.to.std] <- apply(X[,covs.to.std], 2, function(x) x / (2*sd(x)))

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
standardize <- TRUE
lasso.d1.cv <- glmnet::cv.glmnet(x=X.mat, y=d1, alpha=1, standardize = standardize,
                                 penalty.factor = penalties)

lasso.d1 <- glmnet::glmnet(x=X.mat, y=d1, alpha=1, 
                           lambda = c(lasso.d1.cv$lambda.1se, lasso.d1.cv$lambda.min),
                           penalty.factor = penalties, standardize = standardize)

# Temperature square
lasso.d2.cv <- glmnet::cv.glmnet(x=X.mat, y=d2, alpha=1, standardize = standardize,
                                 penalty.factor = penalties)
                                 
lasso.d2 <- glmnet::glmnet(x=X.mat, y=d2, alpha=1, 
                           lambda = c(lasso.d2.cv$lambda.1se, lasso.d2.cv$lambda.min),
                           penalty.factor = penalties, standardize = standardize)

# Outcome
lasso.y.cv <- glmnet::cv.glmnet(x=X.mat, y=y, alpha=1, standardize = standardize,
                                penalty.factor = penalties, maxit = 1000L)

plot(lasso.y.cv)

lasso.y <- glmnet::glmnet(x=X.mat, y=y, alpha=1, 
                          lambda = c(lasso.y.cv$lambda.1se, lasso.y.cv$lambda.min),
                          penalty.factor = penalties, standardize = standardize)

# store coefficients different from 0
# according to 1stdv rule 
tol <- 1e-8

lam.std.sel <- 
  bind_rows(
    tibble(term = rownames(lasso.d1$beta), beta = lasso.d1$beta[,1], 
           var = "d1"), 
    tibble(term = rownames(lasso.d2$beta), beta = lasso.d2$beta[,1], 
           var = "d2"), 
    tibble(term = rownames(lasso.y$beta), beta = lasso.y$beta[,1], 
           var = "y")) %>% 
  mutate(chosen = ifelse(beta > tol, 1, 0), 
         forced = ifelse(term %in% sel.nopen, 1, 0), 
         selected = ifelse(chosen == 1 | forced == 1, 1, 0))
  
# according to min mse 
lam.min.sel <-   bind_rows(
  tibble(term = rownames(lasso.d1$beta), beta = lasso.d1$beta[,2], 
         var = "d1"), 
  tibble(term = rownames(lasso.d2$beta), beta = lasso.d2$beta[,2], 
         var = "d2"), 
  tibble(term = rownames(lasso.y$beta), beta = lasso.y$beta[,2], 
         var = "y")) %>% 
  mutate(chosen = ifelse(beta > tol, 1, 0), 
         forced = ifelse(term %in% sel.nopen, 1, 0), 
         selected = ifelse(chosen == 1 | forced == 1, 1, 0))

# some diagnositics -------------------------------------------------------

lam.std.sel %>% group_by(var) %>% 
  summarize(selected = sum(selected), 
            forced = sum(forced), 
            chosen = sum(chosen))

lam.min.sel %>% group_by(var) %>% 
  summarize(selected = sum(selected), 
            forced = sum(forced), 
            chosen = sum(chosen))

countries <- lam.min.sel %>% 
  filter(selected == 1, str_detect(term, "time")) %>% 
  mutate(iso = substr(term, 4, 6)) %>% 
  group_by(iso, term) %>% summarize(selected = sum(selected)) %>% 
  ungroup() %>% 
  arrange(term) %>% 
  mutate(term = str_split(term, ":", simplify = TRUE)[,2]) %>% 
  group_by(iso) %>% tally() %>% 
  arrange(-n) 
  
unique(df$iso)[!unique(df$iso) %in% countries$iso]


# create response function ------------------------------------------------

include.std <- lam.std.sel %>% filter(selected == 1) %>% pull(term) %>% unique()
include.min <- lam.min.sel %>% filter(selected == 1) %>% pull(term) %>% unique()

reg.df.std <- bind_cols(y = df$y, xtreat1 = df$x1, xtreat2 = df$x2, 
                    as_tibble(X[, include.std]))  %>% 
  janitor::clean_names()

reg.df.min <- bind_cols(y = df$y, xtreat1 = df$x1, xtreat2 = df$x2, 
                        as_tibble(X[, include.min]))  %>% 
  janitor::clean_names()

# Check i get whats going on. 
var.df <- bind_rows(tibble(var = include.std, type = "std"), 
                    tibble(var = include.min, type = "min"))
                    
# Total number of potential parameters
# Ni -1 (ref) * 3 (unit FE, 2& time controls), Nt-1, 2 precip
max.K <- 4
n.unitFE   <- (length(unique(df$iso)))-1
n.timeFEs  <- (length(unique(df$time1)))-1
n.trends   <- max.K * length(unique(df$iso))
n.controls <- length(c("w1", "w2"))
n.unitFE + n.timeFEs + n.trends + n.controls

# Which do we keep?
var.df <- var.df %>% 
  mutate(fe = var) %>% 
  mutate(fe.type = case_when(
    str_detect(fe, "yearFE") ~ "time", 
    str_detect(fe, "time1") ~ "trend 1", 
    str_detect(fe, "time2") ~ "trend 2",
    str_detect(fe, "time3") ~ "trend 3",
    str_detect(fe, "time4") ~ "trend 4",
    str_detect(fe, "w1|w2") ~ 'precip', 
    !str_detect(fe, "time|year") ~ "unit"
  ))

var.df %>% 
  group_by(fe.type, type) %>% 
  tally() %>% arrange(type)

vars.std <- names(reg.df.std)[-1] %>% 
  paste0(collapse = "+")

vars.min <- names(reg.df.min)[-1] %>% 
  paste0(collapse = "+")

m.dml.std <- feols(as.formula(paste0("y ~ ", vars.std)), 
                   bind_cols(reg.df.std, iso= df$iso), 
               cluster = 'iso')

m.dml.min <- feols(as.formula(paste0("y ~ ", vars.min)), 
                   bind_cols(reg.df.min, iso= df$iso), 
                   cluster = 'iso')

# plot --------------------------------------------------------------------

rf.dml.std <- predict_poly(m.dml.std, "xtreat", 0, 30, 14, ci_level = 95, id.col = "DML STD")
rf.dml.min <- predict_poly(m.dml.min, "xtreat", 0, 30, 14, ci_level = 95, id.col = "Double lasso")
rf.bhm     <- predict_poly(reg2, "x", 0, 30, 14, ci_level = 95, id.col = "BHM")
rf.nps     <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "NPS")
rf.cv      <- predict_poly(reg1, "x", 0, 30, 14, ci_level = 95, id.col = "Cross Validation") 

bind_rows(rf.dml.min, rf.cv) %>% 
  left_join(select(rf.bhm, temp, bhm = response)) %>% 
  left_join(select(rf.nps, temp, nps = response))  %>%
  plot_rf_poly(facet.var = 'id') + 
  geom_line(aes(x = temp, y = bhm), color = "red", linetype= "dashed") + 
  geom_line(aes(x = temp, y = nps), color = "green", linetype= "dashed") + 
  xlab("Temperature (C)") 

ggsave(paste0(dir.out, "double_selection_rf.png"), height = 3, width = 7)


# predicted values for USA and CHN under each model -----------------------

df$pred.dml  <- predict(m.dml.min)
df$pred.cv   <- predict(reg1)
df$pred.bhm  <- predict(reg2)
df$pred.nps  <- predict(reg0)

vars.nox <- names(reg.df.min)[-c(1:3)] %>% 
  paste0(collapse = "+")


mX.dml.min <- feols(as.formula(paste0("xtreat1 ~ ", vars.nox)), 
                   bind_cols(reg.df.min, iso= df$iso), 
                   cluster = 'iso')

regX0 <- feols(data = df, x1 ~ w1 + w2 |
                iso + time1, 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

regX1 <- feols(data = df, x1 ~  w1 + w2 |
                iso + time1 + iso[time1], 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

regX2 <- feols(data = df, x1 ~ + w1 + w2 |
                iso + time1 + iso[time1] + iso[time2], 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

df$pred.dmlX  <- predict(mX.dml.min)
df$pred.cvX   <- predict(regX1)
df$pred.bhmX  <- predict(regX2)
df$pred.npsX  <- predict(regX0)


df %>% 
  filter(iso == "CHN") %>% 
  select(iso, year, y, starts_with("pred")) %>% 
  pivot_longer(starts_with("pred")) %>% 
  mutate(var = ifelse(str_detect(name, "X"), "X", "Y")) %>% 
  mutate(name = str_remove_all(name, "X")) %>% 
  bind_rows(df %>% filter(iso == "CHN") %>% 
              select(year, Y="y", X="x1") %>% 
              pivot_longer(c(Y, X)) %>% mutate(var = name) %>% 
              mutate(name = "Value") 
  ) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = value, color = name)) + 
  facet_wrap(~var, scales = 'free')
