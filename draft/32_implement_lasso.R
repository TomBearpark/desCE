## Implement double LASSO model selection on K=0,1,2 in BHM reg

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

if(Sys.info()['user'] == "tombearpark") source(file.path(code, '/utils/cvFuncs.R'))

set.seed(123)
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out  <- paste0(db, "/out/draft/")

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

# baseline regs for comparisons -------------------------------------------

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

# prepare data --------------------------------------------------------

df <- df %>% mutate(iso = as.factor(iso), yearFE = as.factor(time1))

X <- model.matrix(y ~ w1 + w2 + iso + yearFE + 
                    (time1 + time2 + time3 + time4) * iso, df)[, -1]

# We only want to penalise time trend related stuff
sel.penalized <- str_subset(colnames(X), "time")
sel.nopen     <- str_subset(colnames(X), "time", negate = TRUE)

stopifnot(length(sel.penalized) + length(sel.nopen) == dim(X)[2])

X.pen   <- X[, sel.penalized]
X.nopen <- X[, sel.nopen]

# run double lasso --------------------------------------------------------
pred.vars <- df[, c("y", "x1", "x2")]
# run selection for each treatment equation and for outcome equation
l.out <- run_double_lasso(pred.vars = pred.vars, 
                          X.pen = X.pen, X.nopen = X.nopen,
                          standardize = TRUE)

# extract relevant stuff --------------------------------------------------
tol <- 1e-8
lam.min.sel <- extract_coefs(l.out, tol, sel.nopen = sel.nopen)
include.min <- get_selected_controls(l.out, tol, sel.nopen)

# some diagnositics -------------------------------------------------------

lam.min.sel %>% 
  group_by(var) %>% 
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
  
# Which countries didn't get time trends?
unique(df$iso)[!unique(df$iso) %in% countries$iso]

tibble(fe = include.min, type = "min") %>% 
  mutate(fe.type = case_when(
    str_detect(fe, "yearFE") ~ "time", 
    str_detect(fe, "time1") ~ "trend 1", 
    str_detect(fe, "time2") ~ "trend 2",
    str_detect(fe, "time3") ~ "trend 3",
    str_detect(fe, "time4") ~ "trend 4",
    str_detect(fe, "w1|w2") ~ 'precip', 
    !str_detect(fe, "time|year|w1|w2") ~ "unit"
  )) %>% group_by(fe.type, type) %>% tally() %>% arrange(type)

# create response function ------------------------------------------------

reg.df.min <- bind_cols(y = df$y, xtreat1 = df$x1, xtreat2 = df$x2, 
                        as_tibble(X[, include.min]))  %>% 
  janitor::clean_names()

vars.min <- paste0(names(reg.df.min)[-1], collapse = "+")
  
m.dml.min <- feols(as.formula(paste0("y ~ ", vars.min)), 
                   bind_cols(reg.df.min, iso = df$iso), 
                   cluster = 'iso')

# plot --------------------------------------------------------------------

rf.dml.min <- predict_poly(m.dml.min, "xtreat", 0, 30, 14, ci_level = 95, 
                           id.col = "Double lasso")
rf.nps     <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "NPS")
rf.cv      <- predict_poly(reg1, "x", 0, 30, 14, ci_level = 95, 
                           id.col = "Cross Validation") 
rf.bhm     <- predict_poly(reg2, "x", 0, 30, 14, ci_level = 95, id.col = "BHM")

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
df$pred.nps  <- predict(reg0)
df$pred.cv   <- predict(reg1)
df$pred.bhm  <- predict(reg2)

vars.nox <- paste0(names(reg.df.min)[-c(1:3)], collapse = "+") 
  
mX.dml.min <- feols(as.formula(paste0("xtreat1 ~ ", vars.nox)), 
                   bind_cols(reg.df.min, iso= df$iso), 
                   cluster = 'iso')

regX0 <- feols(data = df, x1 ~ w1 + w2 | iso + time1, 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

regX1 <- feols(data = df, x1 ~  w1 + w2 | iso + time1 + iso[time1], 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

regX2 <- feols(data = df, x1 ~ + w1 + w2 | iso + time1 + iso[time1] + iso[time2], 
              panel.id = c('time1', 'iso'), 
              cluster = ~iso)

df$pred.dmlX  <- predict(mX.dml.min)
df$pred.cvX   <- predict(regX1)
df$pred.bhmX  <- predict(regX2)
df$pred.npsX  <- predict(regX0)


country.toplot <- "GBR"
df %>% 
  filter(iso == !!country.toplot) %>% 
  select(iso, year, y, starts_with("pred")) %>% 
  pivot_longer(starts_with("pred")) %>% 
  mutate(var = ifelse(str_detect(name, "X"), "X", "Y")) %>% 
  mutate(name = str_remove_all(name, "X")) %>% 
  bind_rows(df %>% filter(iso == !!country.toplot) %>% 
              select(year, Y="y", X="x1") %>% 
              pivot_longer(c(Y, X)) %>% mutate(var = name) %>% 
              mutate(name = "Value") 
  ) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = value, color = name)) + 
  facet_wrap(~var, scales = 'free') + 
  ggtitle(country.toplot)
