# Simulation 2 - Discuss the difference on proposed framework of CRM vs MBM
library(dplyr)
library(ggplot2)

# For MBM - Negative binomial
# Key points: Benchmark population is not a random sample, second sample is a random sample

# For CRM - Two binomial process
# Key points: Both samples are random

#########################################
# Function to compute estimates for CRM #
#########################################

# A single sub scenario
crm.1 <- function(simsize = nsim, pplsize = targetN, p1 =scenarios1$p1[1], p2 = scenarios1$p2[1]){
  # Create a list with length = number of scenarios = length(p1.vec)
  # Element:
  sample1.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  # Simulate for each element
  for (i in 1:simsize){
    set.seed(1234+i)
    # Draw individuals for source 1, equivalent to labeling the individual
    sample1.dt[i,] <- rbinom(pplsize, 1, p1)
    # Draw individuals for source 2, equivalent to labeling the individual
    sample2.dt[i,] <- rbinom(pplsize, 1, p2)
  }
  # Marginal counts for sources
  n1 <- apply(sample1.dt,1,sum)
  n2 <- apply(sample2.dt,1,sum)
  # Overlapping individuals, equivalent to labeling the overlapping
  sample1.2.dt<-sample1.dt*sample2.dt
  # Counts of overlapping
  m <- apply(sample1.2.dt,1,sum)
  ####################
  # --- Estimate --- #
  ###################
  # LP Estimator for CRM
  # Est Total N for CRM
  Nhat_LP <- n1*n2/m
  # Chapman Estimator for CRM
  # Est Total N for CRM
  Nhat_Chap <- (n1+1)*(n2+1)/(m+1) - 1
  # Summary
  True.p_m <- p1
  Mean.p_m_crm <- mean(m/n2)
  SD.p_m_crm <- sd(m/n2)
  LB95.p_m_crm <- quantile(m/n2, 0.025) %>% as.numeric()
  UB95.p_m_crm <- quantile(m/n2, 0.975) %>% as.numeric()
  
  True.m <- pplsize * p1 * p2
  Mean.m <- mean(m)
  Med.m <- median(m)
  SD.m <- sd(m)
  Min.m <- min(m)
  Max.m <- max(m)
  
  MSE.Nhat_LP <- mean((Nhat_LP-pplsize)^2)
  Mean.Nhat_LP <- mean(Nhat_LP)
  Med.Nhat_LP <- median(Nhat_LP)
  SD.Nhat_LP <- sd(Nhat_LP)
  LB95.Nhat_LP <- quantile(Nhat_LP, 0.025) %>% as.numeric()
  UB95.Nhat_LP <- quantile(Nhat_LP, 0.975) %>% as.numeric()
  Min.Nhat_LP <- min(Nhat_LP)
  Max.Nhat_LP <- max(Nhat_LP)
  Mean.Nhat_LP_Min.m <- mean(Nhat_LP[which(m == Min.m)])
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
  Mean.Nhat_Chap_Min.m <- mean(Nhat_Chap[which(m == Min.m)])
  Rel.Bias.Nhat_Chap <- mean(Nhat_Chap/pplsize)
  LB95.Rel.Bias.Nhat_Chap <- quantile(Nhat_Chap/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_Chap <- quantile(Nhat_Chap/pplsize, 0.975) %>% as.numeric()
  
  # To analyze the simulation results
  # One-sample t test 
  t.test.2side_LP <- t.test(Nhat_LP, mu = pplsize)
  pvalue.t.2side_LP <- round(t.test.2side_LP$p.value, digits = 5)
  t.test.gtr_LP <- t.test(Nhat_LP, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_LP <- round(t.test.gtr_LP$p.value, digits = 5)
  t.test.less_LP <- t.test(Nhat_LP, mu = pplsize, alternative = "less")
  pvalue.t.less_LP <- round(t.test.less_LP$p.value, digits = 5)
  
  t.test.2side_Chap <- t.test(Nhat_Chap, mu = pplsize)
  pvalue.t.2side_Chap <- round(t.test.2side_Chap$p.value, digits = 5)
  t.test.gtr_Chap <- t.test(Nhat_Chap, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_Chap <- round(t.test.gtr_Chap$p.value, digits = 5)
  t.test.less_Chap <- t.test(Nhat_Chap, mu = pplsize, alternative = "less")
  pvalue.t.less_Chap <- round(t.test.less_Chap$p.value, digits = 5)
  
  # Combine output
  df <- data.frame(# LP estimator
                   MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP,
                   pvalue.t.2side_LP, pvalue.t.gtr_LP, pvalue.t.less_LP, Mean.Nhat_LP_Min.m, 
                   Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
                   # Chapman estimator
                   MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap,
                   pvalue.t.2side_Chap, pvalue.t.gtr_Chap, pvalue.t.less_Chap, Mean.Nhat_Chap_Min.m,
                   Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
                   # plug-in estimate of p1
                   True.p_m, Mean.p_m_crm, SD.p_m_crm, LB95.p_m_crm, UB95.p_m_crm,
                   # counts of overlap
                   True.m, Mean.m, Med.m, SD.m, Min.m, Max.m)
  return(df)
}
# Changing p vector from both low to both high
crm.1.ls <- function(simsize = nsim, pplsize = targetN, p1 = scenarios1$p1 , p2 = scenarios1$p2){
  # Start changing scenarios across the vector of pre-specified p
  # Create list, each element is one sub-scenario of p1 and p2
  s1.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  if(length(p1) == length(p2)) {length.ls <- length(p1)}
  for(l in 1:length.ls){
    s1.ls[[l]] <-  crm.1(simsize = simsize, pplsize = pplsize, p1 = p1[l], p2 = p2[l])
  }
  # Collapse the list
  s1.comb <- do.call(rbind,s1.ls)
  # Add scenarios
  scenario <- data.frame(p1, p2)
  # Final output
  final.s1.comb <- as.data.frame(cbind(scenario,s1.comb))
  return(final.s1.comb)
}



#########################################
# Function to compute estimates for MBM #
#########################################

# A single sub scenario
mbm.1 <- function(simsize = nsim, pplsize = targetN, pB =scenarios1$pB[12], p2 = scenarios1$p2[12]){
  # Create a list with length = number of scenarios = length(p1.vec)
  # Element:
  # all.benchmart.dt <- matrix(0, nrow = sim.B, ncol = pplsize)
  benchmark.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  # For each fixed benchmark, run simsize times sample 2, draw fixed benchmark n.fix.B times
  # Nhat_MBM_nb_mvue <- c()
  # Nhat_MBM_nb_mle <- c()
  # for(j in 1:n.fix.B){
  # Label Benchmark
  # set.seed(2345+j)
  # all.benchmart.dt[j,] <- rbinom(pplsize, 1, pB)
  # benchmark.vec <- all.benchmart.dt[j,]
  # set.seed(1234)
  benchmark.vec <- rbinom(pplsize, 1, pB)
  benchmark.dt <- t(replicate(simsize, benchmark.vec))
  # Simulate for each element
  for (i in 1:simsize){
    set.seed(1234+i)
    # Draw individuals for source 2, equivalent to labeling the individual
    sample2.dt[i,] <- rbinom(pplsize, 1, p2)
  }
  # Marginal counts for sources
  nB <- apply(benchmark.dt,1,sum)
  n2 <- apply(sample2.dt,1,sum)
  # Overlapping individuals, equivalent to labeling the overlapping
  sampleB.2.dt <- benchmark.dt*sample2.dt
  # Counts of overlapping
  m <- apply(sampleB.2.dt,1,sum)
  
  ####################
  # --- Estimate --- #
  ###################
  # Estimator for MBM 
  # 1. Use NB random value generator
  p_m <- m/n2
  # MVUE
  Nhat_MBM_nb_mvue <- ((nB -1)/p_m) + 1
  # Nhat_MBM_nb_mvue_j <- ((nB -1)/p_m) + 1
  
  # 2. Use direct formula
  p_m <- m/n2
  # MLE
  Nhat_MBM_nb_mle <- nB/p_m
  # Nhat_MBM_nb_mle_j <- nB/p_m

  # # Mean_Nhat_MBM_nb_mvue_j <- mean(Nhat_MBM_nb_mvue_j)
  # Nhat_MBM_nb_mvue <- c(Nhat_MBM_nb_mvue, Nhat_MBM_nb_mvue_j)
  # # Mean_Nhat_MBM_nb_mle_j <- mean(Nhat_MBM_nb_mle_j)
  # Nhat_MBM_nb_mle <- c(Nhat_MBM_nb_mle, Nhat_MBM_nb_mle_j)
  # }
  
  # Summary
  True.p_m <- pB
  Mean.p_m_mbm <- mean(m/n2)
  SD.p_m_mbm <- sd(m/n2)
  LB95.p_m_mbm <- quantile(m/n2, 0.025) %>% as.numeric()
  UB95.p_m_mbm <- quantile(m/n2, 0.975) %>% as.numeric()
  
  True.m <- pplsize * pB * p2
  Mean.m <- mean(m)
  Med.m <- median(m)
  SD.m <- sd(m)
  Min.m <- min(m)
  Max.m <- max(m)
  
  MSE.Nhat_MBM_nb_mvue <- mean((Nhat_MBM_nb_mvue-pplsize)^2)
  Mean.Nhat_MBM_nb_mvue <- mean(Nhat_MBM_nb_mvue)
  Med.Nhat_MBM_nb_mvue <- median(Nhat_MBM_nb_mvue)
  SD.Nhat_MBM_nb_mvue <- sd(Nhat_MBM_nb_mvue)
  LB95.Nhat_MBM_nb_mvue <- quantile(Nhat_MBM_nb_mvue, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_nb_mvue <- quantile(Nhat_MBM_nb_mvue, 0.975) %>% as.numeric()
  Min.Nhat_MBM_nb_mvue <- min(Nhat_MBM_nb_mvue)
  Max.Nhat_MBM_nb_mvue <- max(Nhat_MBM_nb_mvue)
  Mean.Nhat_MBM_nb_mvue_Min.m <- mean(Nhat_MBM_nb_mvue[which(m == Min.m)])
  Rel.Bias.Nhat_MBM_nb_mvue <- mean(Nhat_MBM_nb_mvue/pplsize)
  LB95.Rel.Bias.Nhat_MBM_nb_mvue <- quantile(Nhat_MBM_nb_mvue/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_nb_mvue <- quantile(Nhat_MBM_nb_mvue/pplsize, 0.975) %>% as.numeric()

  MSE.Nhat_MBM_nb_mle <- mean((Nhat_MBM_nb_mle-pplsize)^2)
  Mean.Nhat_MBM_nb_mle <- mean(Nhat_MBM_nb_mle)
  Med.Nhat_MBM_nb_mle <- median(Nhat_MBM_nb_mle)
  SD.Nhat_MBM_nb_mle <- sd(Nhat_MBM_nb_mle)
  LB95.Nhat_MBM_nb_mle <- quantile(Nhat_MBM_nb_mle, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_nb_mle <- quantile(Nhat_MBM_nb_mle, 0.975) %>% as.numeric()
  Min.Nhat_MBM_nb_mle <- min(Nhat_MBM_nb_mle)
  Max.Nhat_MBM_nb_mle <- max(Nhat_MBM_nb_mle)
  Mean.Nhat_MBM_nb_mle_Min.m <- mean(Nhat_MBM_nb_mle[which(m == Min.m)])
  Rel.Bias.Nhat_MBM_nb_mle <- mean(Nhat_MBM_nb_mle/pplsize)
  LB95.Rel.Bias.Nhat_MBM_nb_mle <- quantile(Nhat_MBM_nb_mle/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_nb_mle <- quantile(Nhat_MBM_nb_mle/pplsize, 0.975) %>% as.numeric()
  
  
  # To analyze the simulation results
  # One-sample t test 
  t.test.2side_MBM_nb_mvue <- t.test(Nhat_MBM_nb_mvue, mu = pplsize)
  pvalue.t.2side_MBM_nb_mvue <- round(t.test.2side_MBM_nb_mvue$p.value, digits = 5)
  t.test.gtr_MBM_nb_mvue <- t.test(Nhat_MBM_nb_mvue, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_nb_mvue <- round(t.test.gtr_MBM_nb_mvue$p.value, digits = 5)
  t.test.less_MBM_nb_mvue <- t.test(Nhat_MBM_nb_mvue, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_nb_mvue <- round(t.test.less_MBM_nb_mvue$p.value, digits = 5)
  
  t.test.2side_MBM_nb_mle <- t.test(Nhat_MBM_nb_mle, mu = pplsize)
  pvalue.t.2side_MBM_nb_mle <- round(t.test.2side_MBM_nb_mle$p.value, digits = 5)
  t.test.gtr_MBM_nb_mle <- t.test(Nhat_MBM_nb_mle, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_nb_mle <- round(t.test.gtr_MBM_nb_mle$p.value, digits = 5)
  t.test.less_MBM_nb_mle <- t.test(Nhat_MBM_nb_mle, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_nb_mle <- round(t.test.less_MBM_nb_mle$p.value, digits = 5)
  
  # Combine output
  df <- data.frame(# MBM_nb_mvue estimator
                  MSE.Nhat_MBM_nb_mvue, Mean.Nhat_MBM_nb_mvue, Med.Nhat_MBM_nb_mvue, SD.Nhat_MBM_nb_mvue, LB95.Nhat_MBM_nb_mvue, UB95.Nhat_MBM_nb_mvue, Min.Nhat_MBM_nb_mvue, Max.Nhat_MBM_nb_mvue,
                  pvalue.t.2side_MBM_nb_mvue, pvalue.t.gtr_MBM_nb_mvue, pvalue.t.less_MBM_nb_mvue, Mean.Nhat_MBM_nb_mvue_Min.m,
                  Rel.Bias.Nhat_MBM_nb_mvue, LB95.Rel.Bias.Nhat_MBM_nb_mvue, UB95.Rel.Bias.Nhat_MBM_nb_mvue,
                  # MBM_nb_mle estimator
                  MSE.Nhat_MBM_nb_mle, Mean.Nhat_MBM_nb_mle, Med.Nhat_MBM_nb_mle, SD.Nhat_MBM_nb_mle, LB95.Nhat_MBM_nb_mle, UB95.Nhat_MBM_nb_mle, Min.Nhat_MBM_nb_mle, Max.Nhat_MBM_nb_mle,
                  pvalue.t.2side_MBM_nb_mle, pvalue.t.gtr_MBM_nb_mle, pvalue.t.less_MBM_nb_mle, Mean.Nhat_MBM_nb_mle_Min.m,
                  Rel.Bias.Nhat_MBM_nb_mle, LB95.Rel.Bias.Nhat_MBM_nb_mle, UB95.Rel.Bias.Nhat_MBM_nb_mle,
                  # plug-in estimate of pB
                  True.p_m, Mean.p_m_mbm, SD.p_m_mbm, LB95.p_m_mbm, UB95.p_m_mbm,
                  # counts of overlap
                  True.m, Mean.m, Med.m, SD.m, Min.m, Max.m)
  return(df)
}
# Changing p vector from both low to both high
mbm.1.ls <- function( simsize = nsim, pplsize = targetN, pB = scenarios1$pB , p2 = scenarios1$p2){
  # Start changing scenarios across the vector of pre-specified p
  # Create list, each element is one sub-scenario of pB and p2
  s1.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  if(length(pB) == length(p2)) {length.ls <- length(pB)}
  for(l in 1:length.ls){
    s1.ls[[l]] <-  mbm.1(simsize = simsize, pplsize = pplsize, pB = pB[l], p2 = p2[l])
  }
  # Collapse the list
  s1.comb <- do.call(rbind,s1.ls)
  # Add scenarios
  scenario <- data.frame(pB, p2)
  # Final output
  final.s1.comb <- as.data.frame(cbind(scenario,s1.comb))
  return(final.s1.comb)
}



#######################################################
# Function to plot CRM vs MBM on different estimators #
#######################################################

simu2.plot.CRM.MBM.diff.Est <- function(dt.crm = crm.1.Nhat, dt.mbm = mbm.1.Nhat, p1_pB = p1){
  
  CRM.LP <- subset(dt.crm, p1 == p1_pB, select = c(p1, p2, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
  CRM.LP$Methods <- rep("CRM_LP", nrow(CRM.LP))
  colnames(CRM.LP) <- c("p1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  CRM.Chap <- subset(dt.crm, p1 == p1_pB, select = c(p1, p2, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
  CRM.Chap$Methods <- rep("CRM_Chap", nrow(CRM.Chap))
  colnames(CRM.Chap) <- c("p1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  MBM.nb.mvue <- subset(dt.mbm, pB == p1_pB, select = c(pB, p2, Rel.Bias.Nhat_MBM_nb_mvue, LB95.Rel.Bias.Nhat_MBM_nb_mvue, UB95.Rel.Bias.Nhat_MBM_nb_mvue))
  MBM.nb.mvue$Methods <- rep("MBM_nb_mvue", nrow(MBM.nb.mvue))
  colnames(MBM.nb.mvue) <- c("p1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  MBM.nb.mle <- subset(dt.mbm, pB == p1_pB, select = c(pB, p2, Rel.Bias.Nhat_MBM_nb_mle, LB95.Rel.Bias.Nhat_MBM_nb_mle, UB95.Rel.Bias.Nhat_MBM_nb_mle))
  MBM.nb.mle$Methods <- rep("MBM_nb_mle", nrow(MBM.nb.mle))
  colnames(MBM.nb.mle) <- c("p1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, MBM.nb.mvue, MBM.nb.mle))
  
  p <- ggplot(CRM.MBM.given.p1, aes(x = p2, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
        geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
        geom_line(aes(y = Rel.Bias.Nhat), size = 0.8) +
        geom_point(aes(y = Rel.Bias.Nhat)) +
        geom_hline(aes(yintercept=1), color = "black",
                   linetype="dashed", size=0.5) +
        ylab("Nhat/N") + 
        xlab("Sampling Probability for second sample (p2)") +
        ggtitle(paste("CRM vs MBM (N = 5000) -- Plot for", "p1 (pB) =", p1_pB, sep = " ")) +
        theme_classic() +
        theme(plot.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16)) 
 return(p)     
}

#######################################################
# Function to plot CRM vs MBM on plug-in estimate p_m #
#######################################################

simu2.plot.CRM.MBM.diff.p_m <- function(dt.crm = crm.1.Nhat, dt.mbm = mbm.1.Nhat, p1_pB = p1){
  CRM.p_m <- subset(dt.crm, p1 == p1_pB, select = c(p1, p2, Mean.p_m_crm, SD.p_m_crm, LB95.p_m_crm, UB95.p_m_crm))
  CRM.p_m$Methods <- rep("CRM-Est(p1)", nrow(CRM.p_m))
  colnames(CRM.p_m) <- c("p1", "p2", "Mean.p_m","SD.p_m", "LB95.p_m", "UB95.p_m","Methods")
  MBM.p_m <- subset(dt.mbm, pB == p1_pB, select = c(pB, p2, Mean.p_m_mbm, SD.p_m_mbm, LB95.p_m_mbm, UB95.p_m_mbm))
  MBM.p_m$Methods <- rep("MBM-Est(pB)", nrow(MBM.p_m))
  colnames(MBM.p_m) <- c("p1", "p2", "Mean.p_m", "SD.p_m", "LB95.p_m", "UB95.p_m", "Methods")
  
  CRM.MBM.p_m <- as.data.frame(rbind(CRM.p_m, MBM.p_m))
  
  p <- ggplot(CRM.MBM.p_m, aes(x = p2, y = Mean.p_m, group = Methods, color = Methods)) + 
    geom_ribbon(aes(ymin = LB95.p_m, ymax = UB95.p_m), linetype=2, alpha=0.1) +
    geom_line(aes(y = Mean.p_m), size = 0.8) +
    geom_point(aes(y = Mean.p_m)) +
    geom_hline(aes(yintercept=p1_pB), color = "black",
               linetype="dashed", size=0.5) +
    ylab("Mean Plug-in Estimate of pB") + 
    xlab("Sampling Probability for second sample (p2)") +
    ggtitle(paste("CRM vs MBM -- Plot for", "Plug-in Est p1(pB) =", p1_pB, sep = " ")) +
    theme_classic() +
    theme(plot.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) 
  return(p)     
  }
  

#######################################################
# Function to plot CRM vs MBM on relative p1/p2 ratio #
#######################################################


simu2.plot.CRM.MBM.p1.to.p2 <- function(dt.crm = crm.1.Nhat, dt.mbm = mbm.1.Nhat, p1_pB = p1){
  CRM.LP <- subset(dt.crm, select = c(rel.p1.to.p2, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
  CRM.LP$Methods <- rep("CRM_LP", nrow(CRM.LP))
  colnames(CRM.LP) <- c("Ratio.p1.p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  CRM.Chap <- subset(dt.crm, select = c(rel.p1.to.p2, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
  CRM.Chap$Methods <- rep("CRM_Chap", nrow(CRM.Chap))
  colnames(CRM.Chap) <- c("Ratio.p1.p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  MBM.nb.mvue <- subset(dt.mbm,  select = c(rel.pB.to.p2, Rel.Bias.Nhat_MBM_nb_mvue, LB95.Rel.Bias.Nhat_MBM_nb_mvue, UB95.Rel.Bias.Nhat_MBM_nb_mvue))
  MBM.nb.mvue$Methods <- rep("MBM_nb_mvue", nrow(MBM.nb.mvue))
  colnames(MBM.nb.mvue) <- c("Ratio.p1.p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  MBM.nb.mle <- subset(dt.mbm, select = c(rel.pB.to.p2, Rel.Bias.Nhat_MBM_nb_mle, LB95.Rel.Bias.Nhat_MBM_nb_mle, UB95.Rel.Bias.Nhat_MBM_nb_mle))
  MBM.nb.mle$Methods <- rep("MBM_nb_mle", nrow(MBM.nb.mle))
  colnames(MBM.nb.mle) <- c("Ratio.p1.p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, MBM.nb.mvue, MBM.nb.mle))
  
  p <- ggplot(CRM.MBM.given.p1, aes(x = Ratio.p1.p2, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
    # geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
    # geom_line(aes(y = Rel.Bias.Nhat), size = 0.8) +
    geom_point(aes(y = Rel.Bias.Nhat)) +
    geom_hline(aes(yintercept=1), color = "black",
               linetype="dashed", size=0.5) +
    ylab("Nhat/N") + 
    xlab("Ratio of p1 to p2") +
    ggtitle("CRM vs MBM (N = 5000) -- Plot for ratio of p1 to p2") +
    theme_classic() +
    theme(plot.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) 
  return(p)     
}
