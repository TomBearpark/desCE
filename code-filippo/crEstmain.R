rm(list = ls())

if(!require(pacman)) install.packages('pacman')
if(!require(useful)) devtools::install_github("TomBearpark/useful")

pacman::p_load(fixest, broom, tidyverse, useful, marginaleffects, sandwich, 
               summclust)
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

# EHW HC3 
m.k0 <- feols(y ~ temp1 + temp2 + precip1 + precip2 |
                    country + time + country[time], 
              df, cluster = "country")

# CR0 
m.k0.nossc <- feols(y ~ temp1 + temp2 + precip1 + precip2 |
                    country + time + country[time], 
              df, cluster = "country", ssc=ssc(adj=FALSE, fixef.K = "full", cluster.adj=FALSE))

# CR1 (equation 8 in MacKinnon. Nielsen, and Webb (2023, JoE))
m.k0.cr1 <- feols(y ~ temp1 + temp2 + precip1 + precip2 |
                      country + time + country[time], 
                    df, cluster = "country", ssc=ssc(adj=FALSE, fixef.K = "nested", cluster.adj=TRUE)) # defaul cluster


tempvars <- c("temp1", "temp2")

# Compare default SEs 
se.k0  <- se(vcov(m.k0))[tempvars]
se.k0.nossc <- se(vcov(m.k0.nossc))[tempvars]
se.k0.cr1  <- se(vcov(m.k0.cr1))[tempvars]


# manual -----
df$time.f <- as.factor(df$time)

X <- model.matrix(y ~ -1 + temp1 + temp2 + precip1 + precip2 + country + time.f + country:time, df)
X <- X[,-ncol(X)]
eps <- m.k0$residuals
XX <- t(X) %*% X
XXinv <- solve(XX)

#df.adj <- nrow(df)/(nrow(df)-ncol(X))
#vcov.Het <- XXinv %*% t(X) %*% diag(eps^2) %*% X %*% XXinv*df.adj
#cbind(sqrt(diag(vcov.Het[tempvars,tempvars])), se.k0.het)

clusters <- unique(df$country)
G <- length(clusters)

eps <- as.matrix(m.k0$residuals)

D <- fastDummies::dummy_columns(df, select_columns = c("country"), remove_selected_columns = TRUE)
D <- as.matrix(D[,str_detect(colnames(D), 'country')])

S <- (eps %*% t(eps)) * (D %*% t(D))

# standard CR0 estimator (no adjustments)
vcov.CL <- XXinv %*% t(X) %*% S %*% X %*% XXinv 

# CR1 estimator (with cluster size and dof adjustment akin to HC1)
N <- nrow(X) # panel size
K <- ncol(X) # coeffs + fixed effects 

# Difference with fixest:
# we are clustering at the county level, thus county FE and country:time are nested
# however, fixest seems to count country:time as non-nested and also adds a + 1
# to get similar (yet not exact) results substitute Keff <- K + 1 - G 
Keff <- K - 2*G  

vcov.CL.adj <- vcov.CL * (G/(G-1)) * ((N-1)/(N-Keff)) 
cbind(sqrt(diag(vcov.CL[tempvars,tempvars])), se.k0.nossc, 
      sqrt(diag(vcov.CL.adj[tempvars,tempvars])), se.k0.cr1)


vcov3 = function(object, cluster){
  
  data = fixest:::fetch_data(object)
  data = as.data.frame(data) # need this guy to process tibbles
  cluster = data[,cluster]
  unique_cluster = as.vector(unique(cluster))
  G = length(unique(cluster))
  fml_all = object$fml_all
  
  obj_call = object$call
  obj_call$data = quote(data[cluster != g,]) # need to remove a cluster not keep only one
  
  beta_centered = lapply(c(unique_cluster), function(g) coef(eval(obj_call)) - coef(object)) 
  vcov = lapply(seq_along(unique_cluster), function(g) tcrossprod(as.matrix(beta_centered[[g]])))
  vcov = Reduce("+", vcov) * (G-1)/G # correction factor has changed
  vcov
}

rbind(se(vcov3(m.k0.cr1, "country"))[tempvars], 
      se.k0,
      se.k0.nossc,
      se.k0.cr1)



