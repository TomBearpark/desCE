pacman::p_load(fixest, parallel, ggplot2, reshape2, wesanderson, latex2exp,
               fastDummies, tidyverse)

if (Sys.getenv("USER") == "fpalomba") {
  path <- "/Users/fpalomba/Dropbox (Princeton)/projects/BP_2023_fesearch/"
} else {
  path <- ""
}

source(paste0(path, "code/funs.R"))

Nsimul <- 50       # number of samples drawn
Nobs <- 1000       # number of obs per sample
g <- 0.5           # fraction in group 1
p <- c(1/3,2/3)    # pscore for each group (p_0, p_1)
tau <- 1           # treatment shift
s2n <- 1           # signal-to-noise ratio in potential outcomes
Ncores <- parallel::detectCores() - 1

# simulate bunch of datasets
set.seed(8894)

toplot <- NULL
for (Nobs in c(1000, 5000, 10000)) {
  
  for (s2n in c(1, 3, 5)) {
    dfSims <- lapply(c(1:Nsimul), function(i) simulModel(N=Nobs, g=g, p=p, tau=tau, s2n=s2n))
    res <- parallel::mclapply(dfSims, function(df) gfe(df, out.var="Y", covs.var="D", nGroups = 2),
                              mc.cores=Ncores, mc.set.seed=8894)
    
    toplot  <- rbind(toplot,
                    cbind(c(1:Nsimul), rep(Nobs, Nsimul), rep(s2n, Nsimul),
                          unlist(lapply(c(1:length(res)), function(i) res[[i]][["theta"]] ))))
  }
}

toplot <- as.data.frame(toplot)
names(toplot) <- c("simul", "Nobs", "s2n", "theta")
toplot$signal <- as.factor(toplot$s2n)
toplot$Nobs <- as.factor(toplot$Nobs)
theme_set(theme_bw())
p <- ggplot(toplot) +
  geom_vline(xintercept = tau, color = "#FF0000") +
  geom_histogram(aes(x=theta)) + 
  facet_grid(signal~Nobs, labeller = label_both) + 
  xlab(TeX("$\\tau$")) + ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot2::ggsave(p, file=paste0(path, "simuls/gfe_NxS.png"), width=10, height = 8, dpi="retina")

# run base regressions just to fumble
# source(paste0(path, "code/baseRegRun.R"))







