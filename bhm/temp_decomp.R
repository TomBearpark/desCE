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

for(kk in 1:12){
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


mm <- feols(temp1 ~ 1 | country + year + country[time] + country[time2], df)
df %>% 
  mutate(temp.res = residuals(mm)) %>% 
  filter(iso %in% c("USA", "CHN", "ITA", "GBR")) %>% 
  select(iso, year, y, temp.res = temp.res) %>% 
  pivot_longer(cols = c(temp.res)) %>% 
  ggplot(aes(x = year, y = value)) + 
  geom_hline(yintercept = 0, color = 'grey') + 
  geom_line() + 
  geom_smooth() +
  facet_wrap(~iso + name, scales = 'free') 

fixef(mm)$year %>% plot
