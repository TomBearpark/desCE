## Implement CV model selection on K=0,1,2 in BHM reg

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, haven, janitor, broom, useful)
theme_set(theme_bw())

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

reg0.1 <- feols(data = df, y ~ x1 + x2 + w1 + w2 |
                   ID + time1 + ID[time1] + ID[time2], 
                panel.id = c('time1', 'ID'), cluster = "iso")

rf0.1 <- predict_poly(reg0.1, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
bind_rows(rf0.1) %>% plot_rf_poly(facet.var = 'id')


# CV ----------------------------------------------------------------------

FEs <- list(
  ' ~ w1 + w2 | FES + iso',
  ' ~ w1 + w2 | FES + time1 + iso',
  ' ~ w1 + w2 | FES + time1 + iso + iso[time1]',
  ' ~ w1 + w2 | FES + time1 + iso + iso[time1] + iso[time2]', 
  ' ~ w1 + w2 | FES + time1 + iso + iso[time1] + iso[time2] + iso[time3]',
  ' ~ w1 + w2 | FES+ time1 + iso + iso[time1] + iso[time2] + iso[time3] + iso[time4]'
  # ,  ' ~ w1 + w2 | continent^time1 + iso'
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
  ggplot() + 
  geom_col(aes(x = order, y = rmse, fill = model)) + 
  facet_wrap(~var, scales = 'free') + 
  xlab("Ranking")
