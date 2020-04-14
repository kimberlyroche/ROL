library(compositions)
library(driver)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

simulation <- 2
condition <- 2

# ----------------------------------------------------------------------------------------------------
#   SIMULATION 1: CORRELATION BETWEEN FEATURES FROM FULLY RANDOM SAMPLES
# ----------------------------------------------------------------------------------------------------

if(simulation == 1) {
  # simulate these conditions
  #   1: approx. uniform compositions
  #   2: high vs. log relative abundance (more like what we'd observe in practice)

  # simulate a bunch of compositions WITH DIFFERENT TOTAL COUNTS
  # turn these into counts by simulating: (1) a random total number of reads from a Poisson
  #                                       (2) drawing counts from a multinomial
  # convert these to logratios
  N <- 200
  D <- 50
  pseudocount <- 0.5
  alr_mat <- matrix(NA, N, D-1) # columns are samples
  clr_mat <- matrix(NA, N, D)
  ilr_mat <- matrix(NA, N, D-1)
  for(i in 1:N) {
    if(condition == 1) {
      comps <- compositions::rDirichlet.rcomp(1, rep(1, D))
    } else {
      comps <- rexp(D, 1)
      comps <- comps/sum(comps)
    }
    total_reads <- rnbinom(1, mu=10000, size=10)
    counts <- rmultinom(1, total_reads, comps)
    alr_mat[i,] <- driver::alr(t(counts) + pseudocount)
    clr_mat[i,] <- driver::clr(t(counts) + pseudocount)
    ilr_mat[i,] <- driver::ilr(t(counts) + pseudocount)
  }
  
  # get the correlations between features
  correlations <- data.frame(value=c(), logratio=c())
  
  alr_correlations <- cor(alr_mat)
  correlations <- rbind(correlations,
                        data.frame(value=alr_correlations[upper.tri(alr_correlations, diag=F)], logratio="alr"))
  clr_correlations <- cor(clr_mat)
  correlations <- rbind(correlations,
                        data.frame(value=clr_correlations[upper.tri(clr_correlations, diag=F)], logratio="clr"))
  ilr_correlations <- cor(ilr_mat)
  correlations <- rbind(correlations,
                        data.frame(value=ilr_correlations[upper.tri(ilr_correlations, diag=F)], logratio="ilr"))
  
  correlations$logratio <- as.factor(correlations$logratio)
  
  p <- ggplot(correlations) +
    geom_density(aes(x=value, color=logratio)) +
    xlab("correlation")
  ggsave(paste0("corrdensity_sim",simulation,"_cond",condition,".png"), p, units="in", dpi=150, height=4, width=7)
}

# ----------------------------------------------------------------------------------------------------
#   SIMULATION 2: CORRELATION BETWEEN FEATURES FROM A TIME SERIES
#                 Does autocorrelation between samples produce (apparent) asymmetric correlation
#                 between features?
# ----------------------------------------------------------------------------------------------------

if(simulation == 2) {
  # simulate these conditions
  #   1: low autocorrelation
  #   2: high autocorrelation
  T <- 200
  D <- 50
  pseudocount <- 0.5
  alr_mat <- matrix(NA, T, D-1)
  clr_mat <- matrix(NA, T, D)
  ilr_mat <- matrix(NA, T, D-1)
  alr_props <- matrix(NA, T, D)
  clr_props <- matrix(NA, T, D)
  ilr_props <- matrix(NA, T, D)
  for(i in 1:T) {
    if(condition == 1) {
      comps <- compositions::rDirichlet.rcomp(1, rep(1, D))
    } else {
      comps <- rexp(D, 1)
      comps <- comps/sum(comps)
    }
    total_reads <- rnbinom(1, mu=10000, size=10)
    counts <- rmultinom(1, total_reads, comps)
    alr_draw <- driver::alr(t(counts) + pseudocount)
    clr_draw <- driver::clr(t(counts) + pseudocount)
    ilr_draw <- driver::ilr(t(counts) + pseudocount)
    if(i == 1) {
      alr_mat[i,] <- alr_draw
      clr_mat[i,] <- clr_draw
      ilr_mat[i,] <- ilr_draw
    } else {
      if(condition == 1) {
        ac <- 0.1
      } else {
        ac <- 0.9
      }
      alr_mat[i,] <- ac*alr_mat[i-1,] + (1-ac)*alr_draw
      clr_mat[i,] <- ac*clr_mat[i-1,] + (1-ac)*clr_draw
      ilr_mat[i,] <- ac*ilr_mat[i-1,] + (1-ac)*ilr_draw
    }
    alr_props[i,] <- driver::alrInv(alr_mat[i,])
    clr_props[i,] <- driver::clrInv(clr_mat[i,])
    ilr_props[i,] <- driver::ilrInv(ilr_mat[i,])
  }
  
  if(D <= 10) {
    alr_prop_df <- driver::gather_array(alr_props, "proportion", "sample", "taxon")
    alr_prop_df$taxon <- as.factor(alr_prop_df$taxon)
    clr_prop_df <- driver::gather_array(clr_props, "proportion", "sample", "taxon")
    clr_prop_df$taxon <- as.factor(clr_prop_df$taxon)
    ilr_prop_df <- driver::gather_array(ilr_props, "proportion", "sample", "taxon")
    ilr_prop_df$taxon <- as.factor(ilr_prop_df$taxon)
    
    coul = brewer.pal(4, "Spectral")
    coul = colorRampPalette(coul)(length(unique(prop_df$taxon)))
    
    #p.alr <- ggplot(alr_prop_df, aes(x=sample, y=proportion, fill=taxon)) +
    #  geom_bar(position="fill", stat="identity", width=1) +
    #  scale_fill_manual(values=coul) +
    #  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #  theme(legend.text=element_text(size=8)) +
    #  theme(legend.position="none")
    p.clr <- ggplot(clr_prop_df, aes(x=sample, y=proportion, fill=taxon)) + 
      geom_bar(position="fill", stat="identity", width=1) +
      scale_fill_manual(values=coul) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.text=element_text(size=8)) +
      theme(legend.position="none")
    #p.ilr <- ggplot(ilr_prop_df, aes(x=sample, y=proportion, fill=taxon)) + 
    #  geom_bar(position="fill", stat="identity", width=1) +
    #  scale_fill_manual(values=coul) +
    #  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #  theme(legend.text=element_text(size=8)) +
    #  theme(legend.position="none")
    
    #p <- grid.arrange(p.alr, p.clr, p.ilr, nrow=3)
    ggsave(paste0("timecourse_sim",simulation,"_cond",condition,".png"), p.clr, units="in", dpi=150, height=2, width=8)
  } else {
    # get the correlations between features
    correlations <- data.frame(value=c(), logratio=c())
    
    alr_correlations <- cor(alr_mat)
    correlations <- rbind(correlations,
                          data.frame(value=alr_correlations[upper.tri(alr_correlations, diag=F)], logratio="alr"))
    clr_correlations <- cor(clr_mat)
    correlations <- rbind(correlations,
                          data.frame(value=clr_correlations[upper.tri(clr_correlations, diag=F)], logratio="clr"))
    ilr_correlations <- cor(ilr_mat)
    correlations <- rbind(correlations,
                          data.frame(value=ilr_correlations[upper.tri(ilr_correlations, diag=F)], logratio="ilr"))
    
    correlations$logratio <- as.factor(correlations$logratio)
    
    p <- ggplot(correlations) +
      geom_density(aes(x=value, color=logratio)) +
      xlab("correlation")
    ggsave(paste0("corrdensity_sim",simulation,"_cond",condition,".png"), p, units="in", dpi=150, height=4, width=7)
    
  }
}

