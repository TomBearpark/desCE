if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")
pacman::p_load(fixest, broom, tidyverse, useful, marginaleffects, sandwich)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}

db       <- file.path(root, "BP_2023_fesearch")
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

df <- df %>% 
  select(country, starts_with('time'), y, temp1, temp2, precip1, precip2)

# se ----------------------------------------------------------------------

tol <- 1*10^(-12)

# Paper replication
fm1 <- feols(data = df , y ~ temp1 + temp2 + precip1 + precip2
              | country + time + country[time] + country[time2], 
              panel.id = c('time', 'country'), cluster = 'country')

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
# closest we can get. Its what vignette says should be the same: 
# https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html

se(fm1, ssc = ssc(fixef.K = "full"))
sandwich::vcovCL(lm1, cluster = ~country, type = "HC1")[2:5, 2:5] %>% se()

# hc3 version -------------------------------------------------------------

# Doesn't run
# sandwich::vcovCL(lm1, cluster = ~country, type = "HC3")

# Try jackknife version, which apparently is the same thing as HC3? 
lm.kn <-   sandwich::vcovBS(lm1, cluster = ~country, type = "jackknife")
lm.kn.c <- sandwich::vcovJK(lm1, cluster = ~country, center = "estimate")

se(lm.kn[2:5, 2:5])
se(lm.kn.c[2:5, 2:5])
se(fm1)
# Only about 1% difference, and the CR3 are apparently smaller
100*(se(lm.kn[2:5, 2:5]) - se(fm1)) / se(fm1)
100*(se(lm.kn.c[2:5, 2:5]) - se(fm1)) / se(fm1)


# summcluster  ------------------------------------------------------------

pacman::p_load('summclust')

m0 <- fm1

m1 <- feols(data = df , y ~ temp1 + temp2 + precip1 + precip2 + 
              i(country, time1, ref = "ZIM") + i(country, time2, ref = "ZIM")
             + i(country, ref = 'ZIM') + i(time1)
         , panel.id = c('time', 'country'), cluster = 'country')

m2 <- feols(data = df , y ~ temp1 + temp2 + precip1 + precip2 + 
              i(country, time1, ref = "ZIM") + i(country, time2, ref = "ZIM")
            | country  + time1, 
            panel.id = c('time', 'country'), 
            cluster = "country"
            )

etable(m0, m1, m2, keep = "temp")  

# Check we get same coefficient estimates 
# coef(fm1)[1:2]
# se(fm1)[1:2]
# coef(fm)[2:3]
# se(fm)[2:3]

# etable(fm, fm1, keep = "temp")

# Sumcluster
vars <- c("temp1", "temp2")
s0 <- summclust(m0, cluster = ~country, params = vars) 
s1 <- summclust(m1, cluster = ~country, params = vars)
s2 <- summclust(m2, cluster = ~country, params = vars)

summary(s1)
summary(s2)
# Check everything gives same coefficient estimates 
coef(fm1)[vars]
s1$coef_estimates[vars]
s2$coef_estimates[vars]

vcov(fm1)[vars, vars]
s1$vcov[vars, vars]
s2$vcov[vars, vars]


