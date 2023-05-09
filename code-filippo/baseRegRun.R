# run basic regs
regResults <- parallel::mclapply(dfSims, function(df) regEst(df, tEff="homosk"), mc.cores=Ncores)
regResultsHet <- parallel::mclapply(dfSims, function(df) regEst(df, tEff="hetero"), mc.cores=Ncores)

# plot estimated coefs (heterogeneous TE)
df <- data.frame(simul = rep(c(1:Nsimul), 2), g = c(rep("G0", Nsimul), rep("G1", Nsimul)),
                 "b1"=c(unlist(lapply(regResultsHet, "[[", 1)), unlist(lapply(regResultsHet, "[[", 2))),
                 "b2"=c(unlist(lapply(regResultsHet, "[[", 3)), unlist(lapply(regResultsHet, "[[", 4))),
                 "b3"=c(unlist(lapply(regResultsHet, "[[", 5)), unlist(lapply(regResultsHet, "[[", 6))),
                 "b4"= c(unlist(lapply(regResultsHet, "[[", 7)), unlist(lapply(regResultsHet, "[[", 8))))

toplot <- reshape2::melt(df, id=c("simul", "g"))

toplot$variable <- factor(toplot$variable,
                          levels=c("b1", "b2", "b3", "b4"),
                          labels=c("Only D", "D and G", "D and G*D", "D, G, and G*D"))


theme_set(theme_bw())
p <- ggplot(toplot) + 
  geom_vline(xintercept = tau, color = "#FF0000") +
  geom_vline(xintercept = tau + 1, color = "#00A08A") +
  geom_histogram(aes(x=value, color=g, group=g), fill="white") + 
  facet_wrap(~variable, nrow=2) + 
  scale_color_manual(values=wes_palette("Darjeeling1"), labels=c("Group 0", "Group 1"), name="") + 
  xlab(TeX("$\\beta$")) + ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.justification = c(0, 1), legend.position = c(0, 1),
        legend.background = element_rect(fill='transparent', colour = NA),
        legend.box.background = element_rect(fill='transparent', colour = NA),
        legend.key = element_rect(colour = NA, fill = NA))
ggplot2::ggsave(p, file=paste0(path, "simuls/teHeterog.png"), width=8, height = 8, dpi="retina")

# plot estimated coefs (homogeneous TE)
df <- data.frame(simul = c(1:Nsimul),
                 b1 = unlist(lapply(regResults, "[[", 1)),
                 b2 = unlist(lapply(regResults, "[[", 2)),
                 b3 = unlist(lapply(regResults, "[[", 3)),
                 b4 = unlist(lapply(regResults, "[[", 4)))

toplot <- reshape2::melt(df, id=c("simul"))

toplot$variable <- factor(toplot$variable,
                          levels=c("b1", "b2", "b3", "b4"),
                          labels=c("Only D", "D and G", "D and G*D", "D, G, and G*D"))

theme_set(theme_bw())
p<- ggplot(toplot) + 
  geom_histogram(aes(x=value, color=variable), fill="white") + 
  facet_wrap(~variable, nrow=2) + 
  geom_vline(xintercept = tau) +
  scale_color_manual(values=wes_palette("Darjeeling1")) + 
  xlab(TeX("$\\beta$")) + ylab("") +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot2::ggsave(p, file=paste0(path, "simuls/tehomog.png"), width=8, height = 8, dpi="retina")
