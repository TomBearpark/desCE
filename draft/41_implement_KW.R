## Implement CV model selection on K=0,1,2 in BHM reg

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, haven, janitor, broom, useful)
theme_set(theme_classic())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
  code <- '~/Documents/GitHub/desCE/'
}else{
  stop("error")
}
set.seed(123)
root      <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(root, "/data//")
dir.out  <- paste0(root, "/out/draft/")

source(file.path(code, '/utils/cvFuncs.R'))

# load data ---------------------------------------------------------------

df.kw.raw <- read_dta(paste0(dir.data, 
                             "/Replication_Data_KW/panel_short.dta")
                      ) %>% mutate(t = year - 1900, t2 = t^2) %>% 
  arrange(ID, year) %>% 
  group_by(ID) %>% 
  mutate(
    temp_1 = T_CRU, 
    prec_1 = P_CRU, 
    temp_2 = T_CRU*T_CRU, 
    prec_2 = P_CRU*P_CRU, 
    dT = T_CRU - dplyr::lag(T_CRU, 1),
    dP = P_CRU - dplyr::lag(P_CRU, 1)
  ) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(l_temp            = dplyr::lag(temp_1), 
         d_temp            = temp_1 - l_temp, 
         temp_int_d_temp   = temp_1 * d_temp, 
         l_d_temp          = dplyr::lag(d_temp), 
         temp_int_l_d_temp = temp_1 * l_d_temp,
         
         l_prec            = dplyr::lag(prec_1), 
         d_prec            = prec_1 - l_prec, 
         prec_int_d_prec   = prec_1 * d_prec, 
         l_d_prec          = dplyr::lag(d_prec), 
         prec_int_l_d_prec = prec_1 * l_d_prec
  ) %>% 
  ungroup()

df <- df.kw.raw %>% 
  select(y = dlgdp_pc_usd, x1 = temp_1, x2 = temp_2, 
         w1 = prec_1, w2 = prec_2, FES=StructChange, 
         ID, time1 = t, iso)  %>% 
  mutate(time1, time2 = time1^2, time3 = time1^3, time4 = time1^4, 
         const = 1) %>% 
  na.omit()


# baseline regression for comparison --------------------------------------

rf <- map_dfr(
  0:4, 
  function(K){
    fe <- 'ID+time1'
    if(K>0) fe <- paste0(fe, "+", 
                         paste0("ID[time", 1:K, "]", collapse = "+"))
    reg <- feols(data = df, 
                 as.formula(paste0('y ~ x1 + x2 + w1 + w2 |', fe)), 
                    panel.id = c('time1', 'ID'), cluster = "iso")
    predict_poly(reg, "x", 0, 30, 14, ci_level = 95, id.col = paste0("K=", K)) 
  }
)

bind_rows(rf) %>% 
  group_by(id) %>%
  mutate(opt = which(response == max(response))) %>% 
  
  plot_rf_poly(facet.var = 'id', nrow=1)  +
  geom_vline(aes(xintercept = opt), linetype = "dashed", alpha = .3)  +
  xlab("Temperature (C)") + ylab("")
  

ggsave(file.path(dir.out, "KW_K_poly.png"), 
       width = 10, height = 3)


# CV ----------------------------------------------------------------------

FEs <- list(
  ' ~ w1 + w2 | FES + ID',
  
  ' ~ w1 + w2 | FES + time1 + ID',
  
  ' ~ w1 + w2 | FES + time1 + ID + ID[time1]',
  ' ~ w1 + w2 | FES + time1 + ID + ID[time1] + ID[time2]', 
  ' ~ w1 + w2 | FES + time1 + ID + ID[time1] + ID[time2] + ID[time3]',
  ' ~ w1 + w2 | FES + time1 + ID + ID[time1] + ID[time2] + ID[time3] + ID[time4]'
  
)

# Remove structural change variable which makes things complicated
FEs <- str_remove_all(FEs, "FES")

# Parameters for CV

id.vars   <- c("ID", "time1")
test.prop <- .2
K         <- 100

split(df, id.vars, .2)

results.x1 <- run_cv(df, "x1", FEs, id.vars, test.prop, K = K, na.rm = TRUE)
results.x2 <- run_cv(df, "x2", FEs, id.vars, test.prop, K = K, na.rm = TRUE)
results.y  <- run_cv(df, "y",  FEs, id.vars, test.prop, K = K, na.rm = TRUE)

# Plot result
bind_rows(
  results.x1 %>% mutate(order = row_number(), var = 'Temperature'), 
  results.x2 %>% mutate(order = row_number(), var = 'Temperature Squared'), 
  results.y  %>% mutate(order = row_number(), var = 'GDP-PC Growth'), 
) %>% 
  mutate(model = str_split(model, "\\|", simplify = TRUE)[,2], 
         Model = model) %>% 
  mutate(Model = ifelse(!str_detect(model, "time1\\]"), "K=0", Model),
         Model = ifelse(str_detect(model, "time1\\]"),  "K=1", model), 
         Model = ifelse(str_detect(model, "time2\\]"),  "K=2", Model), 
         Model = ifelse(str_detect(model, "time3\\]"),  "K=3", Model), 
         Model = ifelse(str_detect(model, "time4\\]"),  "K=4", Model), 
  ) %>%
  filter(order == 1)


# LASSO -------------------------------------------------------------------
df <- df %>% mutate(ID = as.factor(ID), yearFE = as.factor(time1)) 

X <- model.matrix(y ~ w1 + w2 + ID + yearFE + 
                    (time1 + time2 + time3 + time4) * ID, df)[, -1]

# We only want to penalise time trend related stuff
sel.penalized <- str_subset(colnames(X), "time")
sel.nopen     <- str_subset(colnames(X), "time", negate = TRUE)

stopifnot(length(sel.penalized) + length(sel.nopen) == dim(X)[2])

X.pen   <- X[, sel.penalized]
X.nopen <- X[, sel.nopen]

pred.vars <- df[, c("y", "x1", "x2")]
# run selection for each treatment equation and for outcome equation
l.out <- run_double_lasso(pred.vars = pred.vars, 
                          X.pen = X.pen, X.nopen = X.nopen,
                          standardize = TRUE)

tol <- 1e-8
lam.min.sel <- extract_coefs(l.out, tol, sel.nopen = sel.nopen)
include.min <- get_selected_controls(l.out, tol, sel.nopen)

## create response function ------------------------------------------------

reg.df.min <- bind_cols(y = df$y, xtreat1 = df$x1, xtreat2 = df$x2, 
                        as_tibble(X[, include.min]))  %>% 
  janitor::clean_names()

vars.min <- paste0(names(reg.df.min)[-1], collapse = "+")

m.dml.min <- feols(as.formula(paste0("y ~ ", vars.min)), 
                   bind_cols(reg.df.min, ID = df$ID, iso = df$iso), 
                   cluster = 'iso')

rf.dml.min <- predict_poly(m.dml.min, "xtreat", 0, 30, 14, ci_level = 95, 
                           id.col = "Double lasso") %>% 
  mutate(facet = id)

# Show CV selected and LASSO selected values
rfs <- rf.dml.min %>% 
  bind_rows(
    rf %>% filter(id=="K=0") %>% mutate(id = "CV: K=0")
  )

rfs %>% 
  rename(Model = id) %>% 
  ggplot() + 
  geom_hline(yintercept = 0, color = "grey", alpha = 0.9)+
  geom_line(aes(x = temp, y = response)) +
  geom_ribbon(aes(x = temp, ymin = lower, ymax = upper), alpha=.4) +
  facet_wrap(~Model) +
  xlab("Temperature (C)")  + 
  theme(legend.position = 'none')

ggsave(paste0(dir.out, "KW_double_selection_rf.png"), height = 3, width = 7)
