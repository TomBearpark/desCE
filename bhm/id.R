## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(Hmisc, fixest, tidyverse, janitor, broom, marginaleffects, 
               useful)
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
         precip1 = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1,
         time2   = time^2, time3 = time^3, time4 = time^4,
         year2   = year^2) %>% 
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
         l.x  = dplyr::lag(x), 
         x2   = x^2, x1 = x) %>% 
  ungroup()  %>% 
  mutate(iso = as.factor(iso), 
         timeCont = paste0(time1,continent))

N  <- nrow(df.sim)
Ni <- length(unique(df.sim$iso))

# baseline regression for comparison --------------------------------------

reg0 <- feols(data = df.sim , y ~  x1 + x2 + w1 + w2 | 
                iso + time1 + iso[time1] + iso[time2], 
              panel.id = c('time1', 'iso'), cluster = ~iso)
rf0 <- predict_poly(reg0, "x", 0, 30, 14, ci_level = 95, id.col = "BHM") 
plot_rf_poly(rf0)

# Rescale time trends to prevent issues with inversions
df.simF <- df.sim %>% mutate(time1 = time1 / 10, 
                             time2 = time1^2, 
                             time3 = time1^3, 
                             time4 = time1^4, 
                             time5 = time1^5, 
                             time6 = time1^6)

# models ------------------------------------------------------------------

ws   <- paste0(" w", 1:2, collapse = " +")
twfe <- " + i(time1) + i(iso)"

FEs.ftest <- map(
  0:6, 
  function(k){
    if(k == 0) {
      trends <- ''
    }
    else{
      twfe <- paste0(twfe, " + ")
      trends <- paste0(" i(iso, time", 1:k, ")", collapse = " + ")
    }
    paste0(ws, twfe, trends, collapse = " +")
  })




# -------------------------------------------------------------------------
# Show ordering according to an F-test on the outcome model

models <- 
  map(FEs.ftest, 
    function(x){
      feols(as.formula(paste0("y ~ x1 + x2 + ", x)), 
            data = df.simF, cluster = "iso")
    })

coefplot(models, keep = "x1")
coefplot(models, keep = "x2")

# Compare AIC/BIC of models
map(models, \(x) BIC(x)) %>% unlist() %>%  plot()
map(models, \(x) AIC(x)) %>% unlist() %>%  plot()

# I dont understand how a model with pol6 6 time trend can be estimated 
# when some units only have 5 observations??
df.simF$pred.m6 <- predict(models[[7]], df.simF)

df.simF %>% 
  filter(iso == "AFG") %>% 
  select(year, starts_with("pred"), y)

tidy(models[[7]]) %>% 
  filter(str_detect(term, "AFG"))

tidy(models[[7]]) %>% filter(str_detect(term, "time6")) %>% 
  ggplot() + geom_histogram(aes(x = p.value))

# -------------------------------------------------------------------------
# Show ordering according to an F-test on the treatment model

models.x <- 
  map(FEs.ftest, 
      function(x){
        feols(as.formula(paste0("x1 ~", x)), 
              data = df.simF, se = "hetero")
      })

map(models.x, \(x) BIC(x)) %>% unlist() %>% plot()
map(models.x, \(x) AIC(x)) %>% unlist() %>% plot()

map(models.x, \(x) fitstat(x, "f")$f$stat) %>% unlist() %>% plot()


hypotheses(
  models.x[[1]], joint = "w"
)
anova(models.x[[1]], models.x[[2]])

fitstat(models.x[[1]], "f")


wald(models.x[[1]])

# CV ---------------------------------------------------------------------

FEs <- list(' ~ w1 + w2 | iso',
            ' ~ w1 + w2 | time1 + iso',
            ' ~ w1 + w2 | time1 + iso + iso[time1]',
            ' ~ w1 + w2 | time1 + iso + iso[time1] + iso[time2]', 
            ' ~ w1 + w2 | time1 + iso + iso[time1] + iso[time2] + iso[time3]',
            ' ~ w1 + w2 | time1 + iso + iso[time1] + iso[time2] + iso[time3] + iso[time4]', 
            ' ~ w1 + w2 | continent^time1 + iso'
)

# Parameters for CV
df.sim$time1continent <- paste0(df.sim$time1, df.sim$continent)
id.vars   <- c("iso")
test.prop <- .2
K         <- 100

split(df.sim, id.vars, .2)

results.x1 <- run_cv(df.sim, "x1", FEs, id.vars, test.prop, K = K)
results.x2 <- run_cv(df.sim, "x2", FEs, id.vars, test.prop, K = K)
results.y  <- run_cv(df.sim, "y", FEs, id.vars, test.prop, K = K)

bind_rows(
  results.x1 %>% mutate(order = row_number(), var = 'x1'), 
  results.x2 %>% mutate(order = row_number(), var = 'x2'), 
  results.y  %>% mutate(order = row_number(), var = 'y'), 
) %>% 
  mutate(model = str_split(model, "\\|", simplify = TRUE)[,2]) %>% 
  filter(!str_detect(model, "4|3")) %>% 
  ggplot() + 
  geom_col(aes(x = order, y = rmse, fill = model)) + 
  facet_wrap(~var, scales = 'free')

m.x1 <- feols(best.model(results.x1), df.sim, combine.quick = FALSE)
m.x2 <- feols(best.model(results.x2), df.sim, combine.quick = FALSE)
m.y  <- feols(best.model(results.y), df.sim, combine.quick = FALSE)

etable(m.x1, m.x2, m.y, fitstat = c("aic", "bic", "r2"))

# Step 1.2, estimate \hat f(x)
for(var in c('x1', 'x2', 'y')) {
  df.sim[paste0('hat.', var)]   <- predict(get(paste0("m.", var)), df.sim)
  df.sim[paste0('tilde.', var)] <- residuals(get(paste0("m.", var)))
}

# illustrate what happened
df.sim %>% 
  filter(iso %in% c("USA", "CHN")) %>% 
  ggplot(aes(x = year)) + 
  geom_line(aes(y = hat.x1)) + geom_line(aes(y = x), color = 'red') + 
  facet_wrap(~iso)

# Step 2: estimate the relationship between x and y
m2   <- feols(tilde.y ~ tilde.x1, df.sim, cluster = 'iso')
m2.2 <- feols(tilde.y ~ tilde.x1 + tilde.x2, df.sim, cluster = 'iso')
etable(m2, m2.2)

# Compare to BHM version
bind_rows(
  rf0, 
  predict_poly(m2.2, "tilde", 0, 30, 14, id.col = "Tilde")  
) %>% plot_rf_poly(fill.var = 'id')

# heterogeneity -----------------------------------------------------------

## non-linear x ----------------------------------------------------------

gg <- 3

df.sim <- df.sim %>% mutate(xbin = Hmisc::cut2(x, g = gg, levels.mean = T))
cuts   <- cut2(df.sim$x, g = gg, onlycuts = T)

m.nlx <- feols(tilde.y ~ i(xbin, tilde.x1) + xbin, df.sim, cluster = 'iso')
m.nlx %>% coefplot(keep = 'tilde')

# Add stuff so we can plot
pdf.nlx <- m.nlx %>% tidy() %>% 
  filter(str_detect(term, "tilde")) %>% 
  mutate(center = as.numeric(str_remove_all(term, "xbin::|:tilde.x1| "))) %>% 
  mutate(max = cuts[-1], min = cuts[-length(cuts)])

rf0$r <- rf0$cat <- NA
rf0$r[1] <- 0

for(tt in seq(2, 31, 1)){
  temp <- rf0$temp[tt]
  if(temp > max(pdf.nlx$max)){
    estimate <- pdf.nlx$estimate[pdf.nlx$max == max(pdf.nlx$max)]
  }else if (temp < min(pdf.nlx$min)) {
    estimate <- pdf.nlx$estimate[pdf.nlx$min == min(pdf.nlx$min)]
  }else{
    estimate <- pdf.nlx %>% 
      mutate(temp = temp) %>% filter(dplyr::between(temp, min, max)) %>% 
      pull(estimate) %>% as.numeric()
  }
  rf0$cat[tt] <- estimate
  rf0$r[tt] <- rf0$r[tt-1] + 1 * estimate
}
rf0$r <- rf0$r - rf0$r[rf0$temp == 14]
rf0
rects <- df.sim %>% 
  group_by(xbin) %>% summarize(minX = min(x), 
                               maxX = max(x), 
                               meanGDPPC = log(mean(gdpCAP_wdi)), 
                               year = mean(year))


ggplot(rf0) + 
  geom_line(aes(x = temp, y = r), color = 'red') + 
  geom_line(aes(x = temp, y = response), color = 'blue') + 
  geom_hline(yintercept = 0) + 
  geom_rect(data = rects, aes(ymin=-Inf, ymax=Inf, xmin=minX,
                            xmax=maxX, fill=year),
            alpha =0.5) +
  scale_fill_viridis_c(direction = -1)

# Why is there some monotonicity? Its due to time period heterogeneity i think
unique(df.sim$xbin)
ii <- df.sim %>% 
  select(xbin, iso, year) %>% 
  distinct() %>% 
  filter(xbin == "20.7691") %>% 
  pull(iso)

shp <- rnaturalearth::ne_download(returnclass = "sf") %>% 
  select(iso = ISO_A3)

shp %>% 
  mutate(check = ifelse(iso %in% ii, "poor", "other")) %>% 
  ggplot() + 
  geom_sf(aes(fill = check))
