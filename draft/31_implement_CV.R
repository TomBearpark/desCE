## Implement CV model selection on K=0,1,2 in BHM reg

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, janitor, broom, useful)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
  code <- '~/Documents/GitHub/desCE/'
}else{
  stop("error")
}
set.seed(123)
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out  <- paste0(db, "/out/draft/")

source(file.path(code, '/utils/cvFuncs.R'))

# load and clean data -----------------------------------------------------

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1 = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1,
         time2   = time^2, time3 = time^3, time4 = time^4,
         year2   = year^2) %>% 
  rename(y = growthWDI, country = countryname) 

df.sim <- df %>% 
  select(time1 = time, iso, y, x = temp1, time2, time3, time4, 
         year, year2, 
         continent, 
         gdpCAP_wdi, w1 = precip1, w2 = precip2) %>% 
  mutate(y = 100 * y) %>% 
  # Add first differences
  group_by(iso) %>% 
  mutate(fd.y = y - dplyr::lag(y), 
         fd.x = x - dplyr::lag(x), 
         l.x  = dplyr::lag(x), 
         x2   = x^2, x1 = x) %>% 
  ungroup()  %>% 
  mutate(iso = as.factor(iso), 
         timeCont = paste0(time1,continent))

N  <- nrow(df.sim)
Ni <- length(unique(df.sim$iso))

# baseline regression for comparison --------------------------------------

reg0.1 <- feols(data = df.sim, y ~ x1 + x2 + w1 + w2 |
                iso + time1 + iso[time1] + iso[time2], 
              panel.id = c('time1', 'iso'), cluster = "iso")


rf0.1 <- predict_poly(reg0.1, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
bind_rows(rf0.1) %>% plot_rf_poly(facet.var = 'id')

# CV ---------------------------------------------------------------------

FEs <- list(
  # ' ~ w1 + w2 | iso',
            ' ~ w1 + w2 | time1 + iso',
            ' ~ w1 + w2 | time1 + iso + iso[time1]',
            ' ~ w1 + w2 | time1 + iso + iso[time1] + iso[time2]', 
            ' ~ w1 + w2 | time1 + iso + iso[time1] + iso[time2] + iso[time3]',
            ' ~ w1 + w2 | time1 + iso + iso[time1] + iso[time2] + iso[time3] + iso[time4]'
            # ,  ' ~ w1 + w2 | continent^time1 + iso'
            )

# Parameters for CV
df.sim$time1continent <- paste0(df.sim$time1, df.sim$continent)
id.vars   <- c("iso", "time1")
test.prop <- .2
K         <- 100

split(df.sim, id.vars, .2)

results.x1 <- run_cv(df.sim, "x1", FEs, id.vars, test.prop, K = K)
results.x2 <- run_cv(df.sim, "x2", FEs, id.vars, test.prop, K = K)
results.y  <- run_cv(df.sim, "y", FEs, id.vars, test.prop, K = K)

# Plot result
bind_rows(
  results.x1 %>% mutate(order = row_number(), var = 'Temperature'), 
  results.x2 %>% mutate(order = row_number(), var = 'Temperature Squared'), 
  results.y  %>% mutate(order = row_number(), var = 'GDP-PC Growth'), 
) %>% 
  mutate(model = str_split(model, "\\|", simplify = TRUE)[,2], 
         Model = model) %>% 
  mutate(Model = ifelse(str_detect(model, "time1\\]"),  "K=1", model), 
         Model = ifelse(!str_detect(model, "time1\\]"), "K=0", Model), 
         Model = ifelse(str_detect(model, "time2\\]"),  "K=2", Model), 
         Model = ifelse(str_detect(model, "time3\\]"),  "K=3", Model), 
         Model = ifelse(str_detect(model, "time4\\]"),  "K=4", Model), 
         ) %>%
  ggplot() + 
  geom_col(aes(x = order, y = rmse, fill = Model)) + 
  facet_wrap(~var, scales = 'free') + 
  xlab("Ranking")

ggsave(paste0(dir.out, "cv_result_BHM.png"), height = 3, width = 8)

# bonus plots -------------------------------------------------------------

m.x1 <- feols(best.model(results.x1), df.sim, combine.quick = FALSE)
m.x2 <- feols(best.model(results.x2), df.sim, combine.quick = FALSE)
m.y  <- feols(best.model(results.y), df.sim, combine.quick = FALSE)

etable(m.x1, m.x2, m.y, fitstat = c("aic", "bic", "r2"))

# Step 1.2, estimate \hat f(x)
for(var in c('x1', 'x2', 'y')) {
  df.sim[paste0('hat.', var)]   <- predict(get(paste0("m.", var)), df.sim)
  df.sim[paste0('tilde.', var)] <- residuals(get(paste0("m.", var)))
}

# illustrate what happened
df.sim %>% 
  filter(iso %in% c("USA", "CHN")) %>% 
  ggplot(aes(x = year)) + 
  geom_line(aes(y = hat.x1)) + geom_line(aes(y = x), color = 'red') + 
  facet_wrap(~iso)

# Step 2: estimate the relationship between x and y
m2   <- feols(tilde.y ~ tilde.x1, df.sim, cluster = 'iso')
m2.2 <- feols(tilde.y ~ tilde.x1 + tilde.x2, df.sim, cluster = 'iso')
etable(m2, m2.2)

# Compare to BHM version
bind_rows(
  rf0.1, 
  predict_poly(m2.2, "tilde", 0, 30, 14, id.col = "Tilde")  
) %>% plot_rf_poly(fill.var = 'id')

