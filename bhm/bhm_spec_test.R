## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")
pacman::p_load(fixest, tidyverse, useful)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out <- paste0(db, "/out/presentation/")

# load data  --------------------------------------------------------------

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1) %>% 
  rename(y = growthWDI, country = countryname) 
for(kk in 1:8){
  df[paste0("time", kk)] <- df$time^kk
}

# compare models ----------------------------------------------------------

reg0 <- feols(data = df , 
             y ~  temp1 + temp2 + precip1 + precip2
          | country + time, 
        panel.id = c('time', 'country'), 
        cluster = ~country)

reg1 <- feols(data = df , 
              y ~  temp1 + temp2 + precip1 + precip2
              | country + time + country[time], 
              panel.id = c('time', 'country'), 
              cluster = ~country)

reg2 <- feols(data = df , 
              y ~  temp1 + temp2 + precip1 + precip2
              | country + time + country[time] + country[time2], 
              panel.id = c('time', 'country'), 
              cluster = ~country)

reg3 <- feols(data = df , 
              y ~  temp1 + temp2 + precip1 + precip2
              | country + time + country[time] + country[time2] + 
                country[time3], 
              panel.id = c('time', 'country'), 
              cluster = ~country)

reg8 <- feols(data = df , 
              y ~  temp1 + temp2 + precip1 + precip2
              | country + time + 
                country[time] + country[time2] + 
                country[time3] + country[time4] + 
                country[time5] + country[time6] + 
                country[time7] + country[time8], 
              panel.id = c('time', 'country'), 
              # vcov = 'hetero'
              cluster = ~country
              )

etable(reg0, reg1, reg2, reg3, reg8,
       fitstat = c("n", "r2", "wr2", "wald", "f.stat", "aic", "bic", "rmse"))

bind_rows(
  predict_poly(reg0, "temp", 0, 35, 14, ci_level = 95, id.col = "0"), 
  predict_poly(reg1, "temp", 0, 35, 14, ci_level = 95, id.col = "1"), 
  predict_poly(reg2, "temp", 0, 35, 14, ci_level = 95, id.col = "2"), 
  predict_poly(reg3, "temp", 0, 35, 14, ci_level = 95, id.col = "3"), 
  predict_poly(reg8, "temp", 0, 35, 14, ci_level = 95, id.col = "8")
) %>% 
  plot_rf_poly(facet.var = 'id')


bind_rows(
  predict_poly(reg0, "temp", 0, 35, 14, ci_level = 95, id.col = "0"), 
  predict_poly(reg1, "temp", 0, 35, 14, ci_level = 95, id.col = "1"), 
  predict_poly(reg2, "temp", 0, 35, 14, ci_level = 95, id.col = "2")
) %>% 
  plot_rf_poly(facet.var = 'id', fill.color = "#fdb927", 
               line.color = "#552583", 
               add_theme = TRUE) + 
  xlab("Temperature (C)")

ggsave(paste0(dir.out, "rf_comparison.png"), height = 3, width = 8)

# diagnostic on trends ----------------------------------------------------

m1 <- feols(data = df , 
      y ~ temp1 + temp2 + precip1 + precip2 + i(country, time)
      |   country + time, 
      panel.id = c('time', 'country'), 
      cluster = ~country)
etable(m1, reg1, keep = "temp")
coefplot(m1, keep = "country")

broom::tidy(m1) %>% 
  filter(str_detect(term, "country")) %>% 
  ggplot() + geom_density(aes(x = statistic))

# rejection of a model? ---------------------------------------------------

# Compare regression residuals with and without the trends
r <- df %>% 
  select(year, country, 
         y, temp1, temp2) %>% 
  mutate(r0 = resid(reg0), 
         r1 = resid(reg1), 
         r2 = resid(reg2))

var(r$y)
var(r$r0)
var(r$r1)
var(r$r2)

r %>% 
  pivot_longer(c("y", starts_with("r")), values_to = "residual", 
               names_to = "model") %>% 
  ggplot() + 
  geom_point(aes(x = residual, y = temp1)) + 
  facet_wrap(~model)
  
cluster_resids <- function(resids, r){
  
  r %>% 
    mutate(k1 = kmeans(resids, 1)$cluster, 
           k2 = kmeans(resids, 2)$cluster, 
           k3 = kmeans(resids, 3)$cluster, 
           k4 = kmeans(resids, 4)$cluster, 
           ) %>% 
    pivot_longer(cols = starts_with("k"), 
                 names_to = "type", values_to = "cluster") %>% 
    mutate(cluster = as.factor(cluster)) %>% 
    ggplot() + 
    geom_point(aes(x = temp1, y = y, color = cluster),
               alpha = .5) + 
    facet_wrap(~type, nrow = 1)
  
}

cluster_resids(r$y, r)
cluster_resids(r$r0, r)
cluster_resids(r$r1, r)
cluster_resids(r$r2, r)

