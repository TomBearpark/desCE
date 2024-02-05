if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")
pacman::p_load(fixest, broom, tidyverse, useful, marginaleffects, sandwich)
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
  mutate(temp1 = UDel_temp_popweight,temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1) %>% 
  rename(y = growthWDI, country = countryname) %>% 
  mutate(country = as.factor(country))

for(kk in 1:8) df[paste0("time", kk)] <- df$time^kk

# se ----------------------------------------------------------------------

tol <- 1*10^(-12)

# Paper replication
fm1 <- feols(data = df , y ~ temp1 + temp2 + precip1 + precip2
              | country + time + country[time] + country[time2], 
              panel.id = c('time', 'country'), 
              cluster = 'country')

# LM version
lm1 <- lm(as.formula(y ~ temp1 + temp2 + precip1 + precip2 
                     + as.factor(country) + as.factor(time1)
                     + country:time1 + country:time2),
          df)

# Check coefficients are the same
stopifnot(coef(fm1)-coef(lm1)[2:5] < tol)

coeftable(fm1)
coeftable(lm1)[2:5, ]


# Replicate fixest SEs from the LM ones. This doesn't replicate, but seems 
# closest we can get. Its what vignette says to do to replicate: 
# https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html

se(fm1, ssc = ssc(fixef.K = "full"))
sandwich::vcovCL(lm1, cluster = ~country, type = "HC1")[2:5, 2:5] %>% se()


# hc3 version -------------------------------------------------------------

# Doesn't run
sandwich::vcovCL(lm1, cluster = ~country, type = "HC3")

# Try jackknife version, which apparently is the same thing
lm.kn <- sandwich::vcovBS(lm1, cluster = ~country, type = "jackknife")
se(lm.kn[2:5, 2:5])
