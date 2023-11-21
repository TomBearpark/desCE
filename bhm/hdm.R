pacman::p_load("hdm")
if(!require(pacman)) install.packages('pacman')
pacman::p_load(Hmisc, fixest, tidyverse, janitor, broom, marginaleffects, 
               useful, glmnet)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
  code <- '~/Documents/GitHub/desCE/'
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out <- paste0(db, "/out/presentation/")

source(file.path(code, '/utils/cvFuncs.R'))


# load and clean data -----------------------------------------------------

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1,
         
         time2   = time^2, time3 = time^3, time4 = time^4,
         year2 = year^2) %>% 
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
         l.x   = dplyr::lag(x), 
         x2 = x^2, x1 = x) %>% 
  ungroup()  %>% 
  mutate(iso = as.factor(iso), 
         timeCont = paste0(time1,continent))

N  <- nrow(df.sim)

Ni <- length(unique(df.sim$iso))

# baseline regression for comparison --------------------------------------

reg0 <- feols(data = df.sim , y ~ x1 + w1 + w2 + 
                i(iso, time1) + i(iso, time2) + i(iso) + i(time1), 
              panel.id = c('time1', 'iso'), cluster = ~iso)

rf0 <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
plot_rf_poly(rf0)

# lasso -------------------------------------------------------------------

x <- model.matrix(reg0)
d <- x[, c("x1")]
w <- x[, colnames(x)[!colnames(x) %in% c("x1")]]
y <- df.sim$y

one.step   <- rlasso(y ~ x)
structural <- rlasso(y ~ w)
propensity <- rlasso(d ~ w)

summary(one.step,   all = FALSE)
summary(structural, all = FALSE)
summary(propensity, all = FALSE)
