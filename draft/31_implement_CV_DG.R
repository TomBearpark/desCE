## Implement CV model on D&G

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, haven, janitor, broom, useful, marginaleffects)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
  code <- '~/Documents/GitHub/desCE/'
}else{
  stop("error")
}
set.seed(123)
root      <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(root, "/data/112583-V1/DG_Reply_2012_Replication_File/")
dir.out  <- paste0(root, "/out/draft/")

source(file.path(code, '/utils/cvFuncs.R'))


# load data, clean up a bit -----------------------------------------------

df <- paste0(dir.data, 
             "DG_CORRECTED_SAMPLE_2012.dta") %>% 
  read_dta() %>% 
  mutate(
    dry_dd89 = dd89 * dry, 
    dry_dd89_sq = dd89 * dd89 * dry,
    dry_prcp =  prcp * dry, 
    dry_prcp_sq = prcp *  prcp * dry, 
    
    irr_dd89 = dd89 * (1-dry), 
    irr_dd89_sq = dd89 * dd89 * (1-dry),
    irr_prcp =  prcp * (1-dry), 
    irr_prcp_sq = prcp *  prcp * (1-dry)
  ) %>% 
  mutate(fland = as.numeric(fland), 
         usdaRegion = str_sub(uuyy, 1, 1), 
         censusRegion = str_sub(ddyy, 1, 1))

unique(df$usdaRegion) %>% sort
unique(df$censusRegion) %>% sort


# replicate column 1a -----------------------------------------------------
s.stat <- df %>% 
  group_by(dry) %>% 
  summarise(
    wy = weighted.mean(y, fland), 
    y  = mean(y), 
    x  = mean(dd89), 
    wx = weighted.mean(dd89, fland))


# xi: areg y dry dry_dd89 dry_dd89_sq dry_prcp dry_prcp_sq 
# irr_dd89 irr_dd89_sq irr_prcp irr_prcp_sq i.year x1-x9 [aweight=fland], a(fips) cluster(fips);

m0 <- feols(
  y ~ dry + 
    dry_dd89 + dry_dd89_sq + dry_prcp + dry_prcp_sq +
    irr_dd89 + irr_dd89_sq + irr_prcp + irr_prcp_sq +
    x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9  | 
    year + fips, 
  data = df, 
  cluster = "fips", 
  weights = ~fland, 
  panel.id = c('fips', 'year'), 
  
  ssc = ssc(fixef.K = "full")
)

# To replicate SEs from their paper, need to specify fixef.K = full
se(m0)['dry_dd89']
se(m0, ssc = ssc(fixef.K = "full"))['dry_dd89']

# Table values: 
# lincom 100*(dry_dd89 + 2*`m_dry_dd89'*dry_dd89_sq);
# lincom 100*(irr_dd89 + 2*`m_irr_dd89'*irr_dd89_sq);
m.x.dry <- weighted.mean(df$dry_dd89, df$fland)
m.x.irr <- weighted.mean(df$irr_dd89, df$fland)
hypotheses(m0, '100 * (dry_dd89 + 2 * m.x.dry * dry_dd89_sq )=0')
hypotheses(m0, '100 * (irr_dd89 + 2 * m.x.irr * irr_dd89_sq )=0')

# model without the heterogeneity, which they prefer ----------------------

run_feols <- function(FE, df){
  ff <- as.formula(
    paste0("y ~ dd89 + dd89_sq + prcp + prcp_sq +", 
           "x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9  | ", FE))
  feols(
    ff, 
    data = df, 
    cluster = "fips", 
    weights = ~fland, 
    panel.id = c('fips', 'year'), 
    ssc = ssc(fixef.K = "full")
  )
}

# TWFE
fe1 <- "year+fips"
m1.1 <- run_feols(fe1, df)

# USDA region
fe2 <- "year^usdaRegion+fips"
m1.2 <- run_feols(fe2, df)

# Census region
fe3 <- "year^censusRegion+fips"
m1.3 <- run_feols(fe3, df)

# State-by-year
fe4 <- "year^statefips+fips"
m1.4 <- run_feols(fe4, df)

# State by year usda region by year
fe5 <- "year^statefips+year^usdaRegion+fips"
m1.5 <- run_feols(fe4, df)

etable(m1.1, m1.2, m1.3, m1.4)
range <- quantile(df$dd89, c(0, .5, 1))

key <- tibble(order = 1:5, fe = c(fe1, fe2, fe3, fe4, fe5))

bind_rows(
  predict_poly(m1.1, "dd89", range[1], range[3], range[2], 
               step.length = 100, id.col = "year+fips") , 
  predict_poly(m1.2, "dd89", range[1], range[3], range[2], 
               step.length = 100, id.col = "year^usdaRegion+fips") , 
  predict_poly(m1.3, "dd89", range[1], range[3], range[2], 
               step.length = 100, id.col = "year^censusRegion+fips") , 
  predict_poly(m1.4, "dd89", range[1], range[3], range[2], 
               step.length = 100, id.col = "year^statefips+fips"), 
  predict_poly(m1.5, "dd89", range[1], range[3], range[2], 
               step.length = 100, id.col = fe5) 
) %>% 
  mutate(fe = id) %>% 
  left_join(key) %>% 
  mutate(fe = fct_reorder(id, order)) %>% 
  plot_rf_poly(facet.var = 'fe', nrow = 1)

# cv ----------------------------------------------------------------------
controls <- " ~ prcp + prcp_sq + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
FEs <- paste0(controls, "|", list(fe1, fe2, fe3, fe4))

# Parameters for CV

id.vars   <- c("fips", "year")
test.prop <- .25
K         <- 50

split(df, id.vars, .2)

results.x1 <- run_cv(df, "dd89", FEs, id.vars, test.prop, K = K, na.rm = TRUE)
results.x2 <- run_cv(df, "dd89_sq", FEs, id.vars, test.prop, K = K, na.rm = TRUE)
results.y  <- run_cv(df, "y",  FEs, id.vars, test.prop, K = K, na.rm = TRUE)

bind_rows(
  results.x1 %>% mutate(order = row_number(), var = 'dd'), 
  results.x2 %>% mutate(order = row_number(), var = 'dd Squared'), 
  results.y  %>% mutate(order = row_number(), var = 'profits'), 
) %>% 
  mutate(model = str_split(model, "\\|", simplify = TRUE)[,2], 
         Model = model) %>% 
  ggplot() + 
  geom_col(aes(x = order, y = rmse, fill = model)) + 
  facet_wrap(~var, scales = 'free') + 
  xlab("Ranking")

