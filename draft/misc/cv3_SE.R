rm(list = ls())

if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")

pacman::p_load(fixest, tidyverse)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  root <- "/Users/fpalomba/Dropbox (Princeton)/projects/"
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

# Default version
m.k0 <- feols(y ~ temp1 + temp2 + precip1 + precip2 + i(time) |
                country  + country[time] + country[time2], 
              df, 
              cluster = "country")


vcov3 = function(object, cluster){
  
  data    = fixest:::fetch_data(object)
  data    = as.data.frame(data) 
  cluster = data[,cluster]
  
  unique_cluster = as.vector(unique(cluster))
  G       = length(unique(cluster))
  fml_all = object$fml_all
  
  obj_call = object$call
  obj_call$data = quote(data[cluster != g,]) 
  
  beta_centered = lapply(c(unique_cluster), 
                         function(g) coef(eval(obj_call)) - coef(object)) 
  vcov = lapply(seq_along(unique_cluster), 
                function(g) tcrossprod(as.matrix(beta_centered[[g]])))
  
  vcov = Reduce("+", vcov) * (G-1)/G 
  vcov
}


v.cr3 <- vcov3(m.k0, "country")

100 * (sqrt(diag(vcov(m.k0))) - sqrt(diag(v.cr3))) / sqrt(diag(v.cr3))
