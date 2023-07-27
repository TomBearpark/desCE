pacman::p_load(fixest, parallel, ggplot2, reshape2, wesanderson, latex2exp,
               fastDummies, tidyverse)

if ((Sys.getenv("USER") == "fpalomba") & (Sys.info()["sysname"] == "Darwin")) {
  path <- "/Users/fpalomba/Dropbox (Princeton)/projects/BP_2023_fesearch/"
  Ncores <- parallel::detectCores() - 1
} else if ((Sys.getenv("USER") == "fpalomba") & (Sys.info()["sysname"] == "Linux")) {
  path <- "/scratch/network/fpalomba/BP_2023_fesearch/"
  Ncores <- 32
} else {
  path <- ""  
}

source(paste0(path, "code/funs.R"))

Nsimul <- 100      # number of samples drawn
Nobs <- 1000       # number of obs per sample
TT <- 1            # number of periods
gNum <- 5          # number of groups
p <- c(1/3,1/3,1/3,2/3,2/3)   # pscore for each group
tau <- 1           # treatment shift
s2n <- 1           # signal-to-noise ratio in potential outcomes (parametrizes distance between mean of groups)

# simulate bunch of datasets
set.seed(8894)
dfSims <- lapply(c(1:Nsimul), function(i) simulModel(N=Nobs, TT=10, gNum=gNum, p=p, tau=tau, s2n=s2n))

save(dfSims, file=paste0(path, "data/dfSims.rdata"))

res <- lapply(dfSims, function(df) gfe(df, out.var="Y", covs.var="D", nGroups = 5)) #, 
                           #mc.cores=Ncores, mc.set.seed=8894)

set.seed(8894)
toplot <- NULL
for (Tobs in c(1, 2, 5)) {
  
  for (s2n in c(1, 3, 5)) {
    dfSims <- lapply(c(1:Nsimul), function(i) simulModel(N=Nobs, TT=Tobs, gNum=g, p=p, tau=tau, s2n=s2n))
    res <- parallel::mclapply(dfSims, function(df) gfe(df, out.var="Y", covs.var="D", nGroups = 2),
                              mc.cores=Ncores, mc.set.seed=8894)
    
    toplot  <- rbind(toplot,
                    cbind(c(1:Nsimul), rep(Tobs, Nsimul), rep(s2n, Nsimul),
                          unlist(lapply(c(1:length(res)), function(i) res[[i]][["theta"]] ))))
  }
}

toplot <- as.data.frame(toplot)
names(toplot) <- c("simul", "Tobs", "s2n", "theta")
toplot$signal <- as.factor(toplot$s2n)
toplot$Tobs <- as.factor(toplot$Tobs)
theme_set(theme_bw())
p <- ggplot(toplot) +
  geom_vline(xintercept = tau, color = "#FF0000") +
  geom_histogram(aes(x=theta)) + 
  facet_grid(signal~Tobs, labeller = label_both) + 
  xlab(TeX("$\\tau$")) + ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot2::ggsave(p, file=paste0(path, "simuls/gfe_TxS.png"), width=10, height = 8, dpi="retina")

save.image(file = paste0(path, "simuls/workspace_gfe.RData"))


# run base regressions just to fumble
# source(paste0(path, "code/baseRegRun.R"))







