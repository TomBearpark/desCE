## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, janitor, broom, marginaleffects, useful)
pacman::p_load(ggpubr, ggrepel)

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


# 1. are x and y trended -----------------------------------------------------
x.trends <- feols(x ~ i(iso) + i(iso, time1), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')

y.trends <- feols(y ~ i(iso) + i(iso, time1), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')
pdf <- bind_rows(
  mutate(tidy(x.trends), var = "Temperature"),  
  mutate(tidy(y.trends), var = "GDP-Growth")
) %>% 
  filter(str_detect(term, "time1"))

pdf %>% 
  ggplot() + 
  geom_rect(aes(xmin=-1.96, xmax=1.96, ymin=-Inf, ymax=Inf), 
            fill = "pink", 
            alpha = .01) + 
  geom_histogram(aes(x = statistic), alpha = .4) + 
  geom_vline(xintercept = 0, color = 'pink') + 
  facet_wrap(~var) + xlab("Time trend t-stat")

ggsave(paste0(dir.out, "time_trends.png"), height = 3, width = 7)

pdf %>% mutate(sig = 1*(abs(statistic)>1.96)) %>% group_by(var) %>% 
  summarize(mean(sig))

# Do the trends correlate
left_join(select(tidy(x.trends), term, xt = estimate), 
          select(tidy(y.trends), term, yt = estimate)) %>% 
  filter(!str_detect(term, "Intercept")) %>% 
  mutate(type = case_when(
    str_detect(term, "time1") & str_detect(term, "iso") ~ "trend", 
    !str_detect(term, "year") & str_detect(term, "iso") ~ "fe iso",
    .default = NA)) %>% 
  ggplot(aes(x = xt, y = yt)) + 
  geom_point() + 
  geom_smooth( method = 'lm', se = FALSE) + 
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")
    )) + 
  facet_wrap(vars(type), scales = 'free')
ggsave(paste0(dir.out, "time_trends_corr.png"), height = 3, width = 7)

#  2. are they quadratic trended  --------------------------------------------
x.trends2 <- feols(x ~ i(iso) + i(iso, time1) + i(iso, time2), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')

y.trends2 <- feols(y ~ i(iso) + i(iso, time1) + i(iso, time2), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')
pdf2 <- bind_rows(
  mutate(tidy(x.trends2), var = "Temperature"),  
  mutate(tidy(y.trends2), var = "GDP-Growth")
) %>% 
  filter(str_detect(term, "time")) %>% 
  mutate(term = case_when(str_detect(term, "time1") ~ "linear", 
                          str_detect(term, "time2") ~ "quadratic" )) 
pdf2 %>% 
  ggplot() + 
  geom_rect(aes(xmin=-1.96, xmax=1.96, ymin=-Inf, ymax=Inf), 
            fill = "pink", 
            alpha = .01) + 
  geom_histogram(aes(x = statistic), alpha = .4) + 
  geom_vline(xintercept = 0, color = 'pink') + 
  facet_wrap(vars(term, var)) + xlab("Time trend t-stat")

ggsave(paste0(dir.out, "time2_trends.png"), height = 3, width = 7)

pdf2 %>% mutate(sig = 1*(abs(statistic)>1.96)) %>% group_by(var, term) %>% 
  summarize(mean(sig))


# 3. is the trend there when we include time FEs --------------------------

x.trends3 <- 
  feols(x ~ i(time1, ref = 1) + i(iso) + i(iso, time1) + i(iso, time2), 
        df.sim, panel.id = c('iso', 'time1'), vcov = 'hetero')

y.trends3 <- 
  feols(y ~i(time1, ref = 1) + i(iso) + i(iso, time1) + i(iso, time2), 
        df.sim, panel.id = c('iso', 'time1'), vcov = 'hetero')

# Check we understand FWL
lm(resid(y.trends3) ~ resid(x.trends3))
feols(y ~ x + i(time1, ref = 1) + i(iso) + i(iso, time1) + i(iso, time2), 
      df.sim, panel.id = c('iso', 'time1'), 
      vcov = 'hetero')


pdf3 <- bind_rows(
  mutate(tidy(x.trends3), var = "Temperature"),  
  mutate(tidy(y.trends3), var = "GDP-Growth")) %>% 
  mutate(term = case_when(
    str_detect(term, "time1")& str_detect(term, "iso")  ~ "linear", 
    str_detect(term, "time2") & str_detect(term, "iso")  ~ "quadratic",
    str_detect(term, "time1") & !str_detect(term, "iso") ~ "fe time", 
    !str_detect(term, "time1") & str_detect(term, "iso") ~ "fe iso",)) %>% 
  filter(!str_detect(term, "Intercept"))

pdf3 %>% 
  filter(var == "Temperature") %>% 
  filter(str_detect(term, "linear")) %>% 
  ggplot() + geom_histogram(aes(x =statistic))

pdf3 %>% 
  ggplot() + 
  geom_rect(aes(xmin=-1.96, xmax=1.96, ymin=-Inf, ymax=Inf), 
            fill = "pink", 
            alpha = .01) + 
  geom_histogram(aes(x = statistic), alpha = .4) + 
  geom_vline(xintercept = 0, color = 'pink') + 
  facet_wrap(vars(var, term), nrow = 2) + 
  xlab("t-stat")
ggsave(paste0(dir.out, "time2_trends_wFE.png"), height = 4, width = 9)

# ggsave(paste0(dir.out, "time2_trends.png"), height = 3, width = 7)
left_join(select(tidy(x.trends3), term, xt = estimate),
          select(tidy(y.trends3), term, yt = estimate)) %>% 
  mutate(term = case_when(
    str_detect(term, "time1")& str_detect(term, "iso")  ~ "linear", 
    str_detect(term, "time2") & str_detect(term, "iso")  ~ "quadratic",
    str_detect(term, "time1") & !str_detect(term, "iso") ~ "fe time", 
    !str_detect(term, "time1") & str_detect(term, "iso") ~ "fe iso",)) %>% 
  filter(!str_detect(term, "Intercept")) %>% 
  ggplot(aes(x = xt, y = yt)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")
    )) + 
  facet_wrap(~term, scales = 'free')
ggsave(paste0(dir.out, "time2_trends_wFE_corr.png"), height = 4, width = 9)

pdf3 %>% mutate(sig = 1*(abs(statistic)>1.96)) %>% group_by(var, term) %>% 
  summarize(mean(sig))




# compare estimates -------------------------------------------------------
x.twfe <- feols(x ~ i(iso) + i(time1), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')

y.twfe <- feols(y ~ i(iso) + i(time1), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')


tibble(df.sim, 
       twfe = predict(x.twfe), 
       quad_timefe = predict(x.trends3)
       ) %>% 
  filter(iso %in% c("GBR", "USA", "CHN")) %>% 
  pivot_longer(cols = c(x, twfe, quad_timefe)) %>% 
  ggplot() + 
  geom_line(aes(x = time1, y = value, color = name)) + 
  facet_wrap(~iso, scales = 'free', ncol = 1) + 
  ggtitle("Temp fitted values")

tibble(df.sim, 
       twfe = predict(y.twfe), 
       quad_timefe = predict(y.trends3)
) %>% 
  filter(iso %in% c("GBR", "USA", "CHN")) %>% 
  pivot_longer(cols = c(y, twfe, quad_timefe)) %>% 
  ggplot() + 
  geom_line(aes(x = time1, y = value, color = name)) + 
  facet_wrap(~iso, scales = 'free', ncol = 1) + 
  ggtitle("GDP-PC gr fitted values")






