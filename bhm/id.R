## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, janitor, broom, marginaleffects, useful)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out <- paste0(db, "/out/presentation/")

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
  select(time1 = time, iso, y, x = temp1, time2 = time2, year, year2) %>% 
  mutate(y = 100 * y)

N  <- nrow(df.sim)
Ni <- length(unique(df.sim$iso))


# regs --------------------------------------------------------------------


