#############################################
# Function to compute estimates for CRM/MBM #
#############################################

# Utility Functions
library(dplyr)

##############################################################################
######### Scenario 1: Incomplete Benchmark + dependent sample 2 ##############
##############################################################################
s2.2.single_p <- function(simsize = nsim, pplsize = targetN, 
                          pB = scenarios2$pB_vec[1], p1_sub = scenarios2$p1_sub_vec[1], p2 = scenarios2$p2_vec[1],
                          dep.factor = scenarios2$dep.factor_vec[1]){
  # Create datasets
  sample1.dt <- matrix(0, nrow = simsize, ncol = pplsize) 
  sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  
  # Label Benchmark members
  benchmark.vec <- c()
  benchmark.vec <- rbinom(pplsize, 1, pB)
  benchmark.dt <- t(replicate(simsize, benchmark.vec))
  # Simulate for each element
  for (i in 1:simsize){
    set.seed(1234+i)
    # Compute the n in Benchmark
    n_B <- sum(benchmark.dt[i,])
    # Incompletely select Benchmark people into samples
    sample1.dt[i,which(benchmark.dt[i,] == 1)] <- rbinom(n_B, 1, p1_sub)
    # Draw individuals for source 2 for CRM or MBM to generate proportion
    # Subgroup with benchmark labels
    sample2.dt[i,which(benchmark.dt[i,] == 1)] <- rbinom(n_B, 1, p2*dep.factor)
    # Subgroup without benchmark labels
    n_not_B <- pplsize - n_B
    sample2.dt[i,which(benchmark.dt[i,] == 0)] <- rbinom(n_not_B, 1, p2)
    
  }
  # Got two data sample, sample1.dt and sample2.dt
  ####################
  # --- Estimate --- #
  ####################
  # /// CRM /// #
  # Overlapping individuals, equivalent to labeling the overlapping
  sample1.2.dt <- sample1.dt*sample2.dt
  # Marginal counts for sources
  n1 <- apply(sample1.dt,1,sum)
  n2 <- apply(sample2.dt,1,sum)
  # Counts of overlapping
  m <- apply(sample1.2.dt,1,sum)
  # LP estimator 
  Nhat_LP <- n1*n2/m
  # Chapman estimator
  Nhat_Chap <- (n1+1)*(n2+1)/(m+1) - 1
  
  # Summary
  True.m.crm <- pplsize * pB * p1_sub * p2 * dep.factor 
  Mean.m.crm <- mean(m)
  Med.m.crm <- median(m)
  SD.m.crm <- sd(m)
  Min.m.crm <- min(m)
  Max.m.crm <- max(m)
  
  MSE.Nhat_LP <- mean((Nhat_LP-pplsize)^2)
  Mean.Nhat_LP <- mean(Nhat_LP)
  Med.Nhat_LP <- median(Nhat_LP)
  SD.Nhat_LP <- sd(Nhat_LP)
  LB95.Nhat_LP <- quantile(Nhat_LP, 0.025) %>% as.numeric()
  UB95.Nhat_LP <- quantile(Nhat_LP, 0.975) %>% as.numeric()
  Min.Nhat_LP <- min(Nhat_LP)
  Max.Nhat_LP <- max(Nhat_LP)
  Mean.Nhat_LP_Min.m <- mean(Nhat_LP[which(m == Min.m.crm)])
  Rel.Bias.Nhat_LP <- mean(Nhat_LP/pplsize)
  LB95.Rel.Bias.Nhat_LP <- quantile(Nhat_LP/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_LP <- quantile(Nhat_LP/pplsize, 0.975) %>% as.numeric()
  
  MSE.Nhat_Chap <- mean((Nhat_Chap-pplsize)^2)
  Mean.Nhat_Chap <- mean(Nhat_Chap)
  Med.Nhat_Chap <- median(Nhat_Chap)
  SD.Nhat_Chap <- sd(Nhat_Chap)
  LB95.Nhat_Chap <- quantile(Nhat_Chap, 0.025) %>% as.numeric()
  UB95.Nhat_Chap <- quantile(Nhat_Chap, 0.975) %>% as.numeric()
  Min.Nhat_Chap <- min(Nhat_Chap)
  Max.Nhat_Chap <- max(Nhat_Chap)
  Mean.Nhat_Chap_Min.m <- mean(Nhat_Chap[which(m == Min.m.crm)])
  Rel.Bias.Nhat_Chap <- mean(Nhat_Chap/pplsize)
  LB95.Rel.Bias.Nhat_Chap <- quantile(Nhat_Chap/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_Chap <- quantile(Nhat_Chap/pplsize, 0.975) %>% as.numeric()
  
  
  # /// MBM (Linked, wrong benchmark, wrong multiplier) /// #
  # Get multiplier
  p_m_link <- m/n2
  # 1. MLE
  Nhat_MBM_lk_mle <- n1/p_m_link
  # 2. MVUE
  Nhat_MBM_lk_mvue <- ((n1 -1)/p_m_link) + 1
  
  # Summary
  MSE.Nhat_MBM_lk_mle <- mean((Nhat_MBM_lk_mle-pplsize)^2)
  Mean.Nhat_MBM_lk_mle <- mean(Nhat_MBM_lk_mle)
  Med.Nhat_MBM_lk_mle <- median(Nhat_MBM_lk_mle)
  SD.Nhat_MBM_lk_mle <- sd(Nhat_MBM_lk_mle)
  LB95.Nhat_MBM_lk_mle <- quantile(Nhat_MBM_lk_mle, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_lk_mle <- quantile(Nhat_MBM_lk_mle, 0.975) %>% as.numeric()
  Min.Nhat_MBM_lk_mle <- min(Nhat_MBM_lk_mle)
  Max.Nhat_MBM_lk_mle <- max(Nhat_MBM_lk_mle)
  Mean.Nhat_MBM_lk_mle_Min.m <- mean(Nhat_MBM_lk_mle[which(m == Min.m.crm)])
  Rel.Bias.Nhat_MBM_lk_mle <- mean(Nhat_MBM_lk_mle/pplsize)
  LB95.Rel.Bias.Nhat_MBM_lk_mle <- quantile(Nhat_MBM_lk_mle/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_lk_mle <- quantile(Nhat_MBM_lk_mle/pplsize, 0.975) %>% as.numeric()
  
  MSE.Nhat_MBM_lk_mvue <- mean((Nhat_MBM_lk_mvue-pplsize)^2)
  Mean.Nhat_MBM_lk_mvue <- mean(Nhat_MBM_lk_mvue)
  Med.Nhat_MBM_lk_mvue <- median(Nhat_MBM_lk_mvue)
  SD.Nhat_MBM_lk_mvue <- sd(Nhat_MBM_lk_mvue)
  LB95.Nhat_MBM_lk_mvue <- quantile(Nhat_MBM_lk_mvue, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_lk_mvue <- quantile(Nhat_MBM_lk_mvue, 0.975) %>% as.numeric()
  Min.Nhat_MBM_lk_mvue <- min(Nhat_MBM_lk_mvue)
  Max.Nhat_MBM_lk_mvue <- max(Nhat_MBM_lk_mvue)
  Mean.Nhat_MBM_lk_mvue_Min.m <- mean(Nhat_MBM_lk_mvue[which(m == Min.m.crm)])
  Rel.Bias.Nhat_MBM_lk_mvue <- mean(Nhat_MBM_lk_mvue/pplsize)
  LB95.Rel.Bias.Nhat_MBM_lk_mvue <- quantile(Nhat_MBM_lk_mvue/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_lk_mvue <- quantile(Nhat_MBM_lk_mvue/pplsize, 0.975) %>% as.numeric()
  
  
  # /// MBM - (Unlinked, wrong benchmark, true multiplier) /// #
  # Find correct number of overlapping
  sample.B.2.dt <- benchmark.dt*sample2.dt
  # Marginal counts for sources
  n1 <- apply(sample1.dt,1,sum)
  n2 <- apply(sample2.dt,1,sum)
  # Counts of overlapping to get correct multiplier
  m.B2 <- apply(sample.B.2.dt,1,sum)
  # Get correct multiplier estimate
  p_m_unlk <- m.B2/n2
  
  # 1. MBM MLE: Undercounted Benchmark divided by true multiplier
  Nhat_MBM_unlk_mle <- n1/p_m_unlk
  # 2. MBM MVUE
  Nhat_MBM_unlk_mvue <- ((n1 -1)/p_m_unlk) + 1
  
  # Summary
  Mean.m.mbm.unlk <- mean(m.B2)
  Med.m.mbm.unlk <- median(m.B2)
  SD.m.mbm.unlk <- sd(m.B2)
  Min.m.mbm.unlk <- min(m.B2)
  Max.m.mbm.unlk <- max(m.B2)
  
  MSE.Nhat_MBM_unlk_mle <- mean((Nhat_MBM_unlk_mle-pplsize)^2)
  Mean.Nhat_MBM_unlk_mle <- mean(Nhat_MBM_unlk_mle)
  Med.Nhat_MBM_unlk_mle <- median(Nhat_MBM_unlk_mle)
  SD.Nhat_MBM_unlk_mle <- sd(Nhat_MBM_unlk_mle)
  LB95.Nhat_MBM_unlk_mle <- quantile(Nhat_MBM_unlk_mle, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_unlk_mle <- quantile(Nhat_MBM_unlk_mle, 0.975) %>% as.numeric()
  Min.Nhat_MBM_unlk_mle <- min(Nhat_MBM_unlk_mle)
  Max.Nhat_MBM_unlk_mle <- max(Nhat_MBM_unlk_mle)
  Rel.Bias.Nhat_MBM_unlk_mle <- mean(Nhat_MBM_unlk_mle/pplsize)
  LB95.Rel.Bias.Nhat_MBM_unlk_mle <- quantile(Nhat_MBM_unlk_mle/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_unlk_mle <- quantile(Nhat_MBM_unlk_mle/pplsize, 0.975) %>% as.numeric()
  
  
  MSE.Nhat_MBM_unlk_mvue <- mean((Nhat_MBM_unlk_mvue-pplsize)^2)
  Mean.Nhat_MBM_unlk_mvue <- mean(Nhat_MBM_unlk_mvue)
  Med.Nhat_MBM_unlk_mvue <- median(Nhat_MBM_unlk_mvue)
  SD.Nhat_MBM_unlk_mvue <- sd(Nhat_MBM_unlk_mvue)
  LB95.Nhat_MBM_unlk_mvue <- quantile(Nhat_MBM_unlk_mvue, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_unlk_mvue <- quantile(Nhat_MBM_unlk_mvue, 0.975) %>% as.numeric()
  Min.Nhat_MBM_unlk_mvue <- min(Nhat_MBM_unlk_mvue)
  Max.Nhat_MBM_unlk_mvue <- max(Nhat_MBM_unlk_mvue)
  Rel.Bias.Nhat_MBM_unlk_mvue <- mean(Nhat_MBM_unlk_mvue/pplsize)
  LB95.Rel.Bias.Nhat_MBM_unlk_mvue <- quantile(Nhat_MBM_unlk_mvue/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_unlk_mvue <- quantile(Nhat_MBM_unlk_mvue/pplsize, 0.975) %>% as.numeric()

  
  # Combine output
  df <- data.frame(# LP estimator
    MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP, Mean.Nhat_LP_Min.m, 
    Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
    # Chapman estimator
    MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap, Mean.Nhat_Chap_Min.m,
    Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
    
    # MBM_lk_mle estimator (linked)
    MSE.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle, Med.Nhat_MBM_lk_mle, SD.Nhat_MBM_lk_mle, LB95.Nhat_MBM_lk_mle, UB95.Nhat_MBM_lk_mle, Min.Nhat_MBM_lk_mle, Max.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle_Min.m, 
    Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle,
    # MBM_direct estimator (linked)
    MSE.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue, Med.Nhat_MBM_lk_mvue, SD.Nhat_MBM_lk_mvue, LB95.Nhat_MBM_lk_mvue, UB95.Nhat_MBM_lk_mvue, Min.Nhat_MBM_lk_mvue, Max.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue_Min.m,
    Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue,
    
    # MBM_unlk_mle estimator (unlinked)
    MSE.Nhat_MBM_unlk_mle, Mean.Nhat_MBM_unlk_mle, Med.Nhat_MBM_unlk_mle, SD.Nhat_MBM_unlk_mle, LB95.Nhat_MBM_unlk_mle, UB95.Nhat_MBM_unlk_mle, Min.Nhat_MBM_unlk_mle, Max.Nhat_MBM_unlk_mle,
    Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle,
    
    # MBM_unlk_mvue estimator (unlinked)
    MSE.Nhat_MBM_unlk_mvue, Mean.Nhat_MBM_unlk_mvue, Med.Nhat_MBM_unlk_mvue, SD.Nhat_MBM_unlk_mvue, LB95.Nhat_MBM_unlk_mvue, UB95.Nhat_MBM_unlk_mvue, Min.Nhat_MBM_unlk_mvue, Max.Nhat_MBM_unlk_mvue,
    Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue,
    
    # counts of overlap
    True.m.crm, Mean.m.crm, Med.m.crm, SD.m.crm, Min.m.crm, Max.m.crm,
    Mean.m.mbm.unlk, Med.m.mbm.unlk, SD.m.mbm.unlk, Min.m.mbm.unlk, Max.m.mbm.unlk
  )
  
  return(df)
}


s2.2.vary_p <- function(simsize = nsim, pplsize = targetN, n.scenario = nrow(scenarios2), 
                        pB = scenarios2$pB_vec, p1_sub = scenarios2$p1_sub_vec, p2 = scenarios2$p2_vec,
                        dep.factor = scenarios2$dep.factor_vec){
  # Start changing scenarios2 across the vector of pre-specified p
  # Create list, each element is one sub-scenario 
  s2.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- n.scenario
  for(l in 1:length.ls){
    s2.ls[[l]] <-  s2.2.single_p(simsize = simsize, pplsize = pplsize, pB = pB[l],  p1_sub = p1_sub[l], 
                                 p2 = p2[l],
                                 dep.factor = dep.factor[l])
  }
  # Collapse the list
  s2.comb <- do.call(rbind,s2.ls)
  # Add scenarios2
  scenario <- data.frame(pB, p1_sub, p2, dep.factor)
  # Final output
  final.s2.comb <- as.data.frame(cbind(scenario,s2.comb))
  return(final.s2.comb)
}



#######################################################
# Function to plot CRM vs MBM on different estimators #
#######################################################

simu2.2.Dt.CRM.MBM.diff.Est <- function(dt.combo = Simu2.s2.Nhat, p1_pB.i = pB.vec[1], p1_sub.i = p1_sub.vec[1], dep.factor.i = dep.factor.vec[1]){
  
  CRM.LP <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i & dep.factor == dep.factor.i, 
                   select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
  CRM.LP$Methods <- rep("CRM_LP", nrow(CRM.LP))
  colnames(CRM.LP) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  CRM.Chap <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i & dep.factor == dep.factor.i, 
                     select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
  CRM.Chap$Methods <- rep("CRM_Chap", nrow(CRM.Chap))
  colnames(CRM.Chap) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.lk_mle <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i & dep.factor == dep.factor.i, 
                       select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle))
  MBM.lk_mle$Methods <- rep("MBM_link_mle", nrow(MBM.lk_mle))
  colnames(MBM.lk_mle) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
 
  MBM.lk_mvue <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i & dep.factor == dep.factor.i, 
                        select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue))
  MBM.lk_mvue$Methods <- rep("MBM_link_mvue", nrow(MBM.lk_mvue))
  colnames(MBM.lk_mvue) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mle <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i & dep.factor == dep.factor.i, 
                         select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle))
  MBM.unlk_mle$Methods <- rep("MBM_unlk_mle", nrow(MBM.unlk_mle))
  colnames(MBM.unlk_mle) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mvue <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i & dep.factor == dep.factor.i, 
                          select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue))
  MBM.unlk_mvue$Methods <- rep("MBM_unlk_mvue", nrow(MBM.unlk_mvue))
  colnames(MBM.unlk_mvue) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, MBM.lk_mle, MBM.lk_mvue, MBM.unlk_mle, MBM.unlk_mvue)) # 
  
  return(CRM.MBM.given.p1)
}
  
simu2.2.plot.CRM.MBM.diff.Est <- function(dt.combo = subdt, p1_pB.i = pB.vec[1], p1_sub.i = p1_sub.vec[1], dep.factor.i = dep.factor.vec[1], ErrorBar = "Yes"){
  if (ErrorBar == "Yes"){
    p <- ggplot(dt.combo, aes(x = p2, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
      geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
      geom_line(aes(y = Rel.Bias.Nhat), size = 0.5) +
      geom_point(aes(y = Rel.Bias.Nhat)) +
      geom_hline(aes(yintercept=1), color = "black",
                 linetype="dashed", size=0.5) +
      scale_y_continuous(breaks = seq(0.4, 1.7, by = 0.05)) +
      # Force all plots have the same limits and breaks
      # scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.4, 1.4, by = 0.05)) + 
      ylab("Nhat/N") + 
      xlab("Sampling Probability for second sample (p2)") +
      ggtitle(paste("(Simu2) CRM vs MBM--Plot for", "pB =", p1_pB.i, "pB_sub =", p1_sub.i, "dep.level =", dep.factor.i, sep = " ")) +
      theme_classic() +
      theme(plot.title = element_text(size = 14),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.position = "bottom") 
  }
  if (ErrorBar == "No"){
    p <- ggplot(dt.combo, aes(x = p2, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
      # geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
      geom_line(aes(y = Rel.Bias.Nhat), size = 0.5) +
      geom_point(aes(y = Rel.Bias.Nhat)) +
      geom_hline(aes(yintercept=1), color = "black",
                 linetype="dashed", size=0.5) +
      scale_y_continuous(breaks = seq(0.4, 1.7, by = 0.05)) +
      # If force all plots have the same limits and breaks
      # scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.4, 1.4, by = 0.05)) +
      ylab("Nhat/N") + 
      xlab("Sampling Probability for second sample (p2)") +
      ggtitle(paste("(Simu2) CRM vs MBM--Plot for", "pB =", p1_pB.i, "pB_sub =", p1_sub.i, "dep.level =", dep.factor.i, sep = " ")) +
      theme_classic() +
      theme(plot.title = element_text(size = 14),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.position = "bottom") 
  }
  return(p)     
}


####################################################################################
# Function to plot MBM combined different undercounting level and dependence level #
####################################################################################

simu2.2.Dt.MBM.diff.Est <- function(dt.combo = Simu2.s2.Nhat, p2.i = p2.vec, 
                                    p1_pB.i = pB.vec, p1_sub = p1_sub.vec, 
                                    dep.factor = dep.factor.vec){

  MBM.lk_mle <- subset(dt.combo, p2 == p2.i & pB == p1_pB.i, 
                       select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle))
  MBM.lk_mle$Methods <- rep("MBM_link_MLE", nrow(MBM.lk_mle))
  colnames(MBM.lk_mle) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mle <- subset(dt.combo, p2 == p2.i & pB == p1_pB.i, 
                         select = c(pB, p1_sub, p2, dep.factor, Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle))
  MBM.unlk_mle$Methods <- rep("MBM_unlk_MLE", nrow(MBM.unlk_mle))
  colnames(MBM.unlk_mle) <- c("pB", "p1_sub", "p2", "dep.factor", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.given.p1 <- as.data.frame(rbind(MBM.lk_mle, MBM.unlk_mle)) # 
  
  return(MBM.given.p1)
}

simu2.2.plot.MBM.diff.Est <- function(dt.combo = subdt[which(subdt$Methods == "MBM_unlk_MLE"),], 
                                      p2.i = p2.vec,
                                      p1_pB.i = pB.vec, 
                                      p1_sub = p1_sub.vec, 
                                      dep.factor = dep.factor.vec,
                                      title = "(A)",
                                      linkage = "link"){
   
    dt.combo$p1_sub <- as.character(dt.combo$p1_sub)
    
    if(linkage == "link"){
      limit_vec <- c(0,4)
      break_vec <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
    }
    if(linkage == "unlink"){
      limit_vec <- c(0,4)
      break_vec <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
    }

    p <- ggplot(dt.combo, aes(x = dep.factor, y = Rel.Bias.Nhat, group = p1_sub, color = p1_sub)) + 
      geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
      geom_line(size = 0.5) +
      geom_point() +
      scale_colour_discrete(name  ="Completeness of the Benchmark") +
      geom_hline(aes(yintercept=1), color = "black",
                 linetype="dashed", size=0.5) +
      
      # If linked MBM (keep all others the same as this)
      # scale_y_continuous(limits = limit_vec, breaks = break_vec) +
      # If unlinked MBM (keep all others the same as this)
      scale_y_continuous(limits = limit_vec, breaks = break_vec) +
      ylab(expression(hat(N)/N)) + 
      xlab("Degree of Sources Dependence") +
      ggtitle(title) +
      theme_classic() +
      theme(plot.title = element_text(size = 14),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16),
            legend.position = "none") 
               
  return(p)     
}
