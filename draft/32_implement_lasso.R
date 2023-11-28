## Implement double LASSO model selection on K=0,1,2 in BHM reg

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

rf0 <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
rf0 %>%  
  plot_rf_poly(facet.var = 'id')

# double selection --------------------------------------------------------


