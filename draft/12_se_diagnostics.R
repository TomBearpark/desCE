if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")

pacman::p_load(fixest, broom, tidyverse, useful, marginaleffects, sandwich, 
               summclust)
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
  mutate(temp1 = UDel_temp_popweight,
         temp2 = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1) %>% 
  rename(y = growthWDI, country = countryname) %>% 
  mutate(country = as.factor(country)) 

for(kk in 1:8) df[paste0("time", kk)] <- df$time^kk

df <- df %>% 
  select(country, starts_with('time'), y, temp1, temp2, precip1, precip2)


# regs --------------------------------------------------------------------
m.k0 <- feols(y ~ temp1 + temp2 + precip1 + precip2 |
                country + time, 
              df, cluster = "country")

m.k0.nofe <- feols(y ~ temp1 + temp2 + precip1 + precip2+
                i(country, ref = "ZIM") + i(time), 
              df, cluster = "country")

etable(m.k0, m.k0.nofe, keep = "temp")
tempvars <- c("temp1", "temp2")

# Compare default SEs 
se.k0      <- se(vcov(m.k0))[tempvars]
se.k0.nofe <- se(vcov(m.k0.nofe))[tempvars]

(se.k0 - se.k0.nofe) / se.k0

# Compare CR3 version from summclust
s0      <- summclust(m.k0, cluster = ~country, params = tempvars) 
s0.nofe <- summclust(m.k0.nofe, cluster = ~country, params = tempvars) 

se.summclust.k0      <- se(s0$vcov)[tempvars]
se.summclust.k0.nofe <- se(s0.nofe$vcov)[tempvars]

(se.summclust.k0 - se.summclust.k0.nofe) / se.summclust.k0
