## 1. Replicate BHM with various time trends
## 2. Look at trends in temperature

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")
pacman::p_load(fixest, broom, tidyverse, useful)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out  <- paste0(db, "/out/draft/")

# load data  --------------------------------------------------------------

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(!is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1) %>% 
  rename(y = growthWDI, country = countryname) 

for(kk in 1:8){
  df[paste0("time", kk)] <- df$time^kk
}

df %>% 
  filter(iso %in% c("USA", "CHN")) %>% 
  select(iso, year, y, Temp = temp1) %>% 
  pivot_longer(cols = c(y, Temp)) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = value)) + 
  facet_wrap(~iso + name, scales = 'free')
ggsave(paste0(dir.out, "non_stationary.png"), height = 6, width = 10)

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
               add_theme = FALSE) + 
  xlab("Temperature (C)") + ylab("")

ggsave(paste0(dir.out, "rf_comparison.png"), height = 3, width = 8)


# trends in temperature  --------------------------------------------------

m1 <- feols(data = df , 
            temp1 ~ i(country, time, ref = "AFG")
            |  country , 
            panel.id = c('time', 'country'), 
            vcov = "hetero")
etable(m1, keep = "time")

broom::tidy(m1) %>% 
  filter(str_detect(term, "country")) %>% 
  mutate(signif = ifelse(p.value<0.5, 1, 0)) %>% 
  summarize(sum(signif), n())

                          