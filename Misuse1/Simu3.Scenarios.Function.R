# Utility Functions
library(dplyr)

#############################################
# Function to compute estimates for CRM/MBM #
#############################################

###########################################################################
######### Scenario 1: Incomplete Benchmark + random sample 2 ##############
###########################################################################
s3.1.single_p <- function(simsize = nsim, pplsize = targetN, 
                        pB = scenarios1$pB_vec[1], p1_sub = scenarios1$p1_sub_vec[1], p2 = scenarios1$p2_vec[1]){
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
    sample2.dt[i,] <- rbinom(pplsize, 1, p2)
    # Sample 2 drawing process is EQUIVALENT to follow (checked already)
    # # Draw individuals for source 2 for CRM or MBM to generate proportion
    # # sample1 counts
    # n_s1 <- sum(sample1.dt[i,])
    # n_not_s1 <- pplsize-n_s1
    # # Unbiased sample 2
    # sample2.dt[i,which(sample1.dt[i,] == 1)] <- rbinom(n_s1, 1, p2)
    # sample2.dt[i,which(sample1.dt[i,] == 0)] <- rbinom(n_not_s1, 1, p2)
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
  True.m.crm <- pplsize * pB * p1_sub * p2
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
  
  # # /// MBM (Unlinked) Undercounted Benchmark with post-sampling adjustment of multiplier /// #
  # # Find adjustment of the multiplier
  # # 1. Make sure sample2 is a random sample - Checked
  # # 2. Get the information of % of undercounting in Benchmark (deflation factor)
  # est_p1_sub <- c()
  # for(i in 1:simsize){
  #   ith.B.sample1.overlap <- sum(benchmark[i,] == 1 & sample1.dt[i,] == 0)
  #   ith.B <- sum(benchmark[i,])
  #   ith.est_p1_sub <- 1-ith.B.sample1.overlap/ith.B
  #   est_p1_sub <- c(est_p1_sub, ith.est_p1_sub)
  # }
  # # Find correct number of overlapping
  # sample.B.2.dt <- benchmark*sample2.dt
  # # Marginal counts for sources
  # n1 <- apply(sample1.dt,1,sum)
  # n2 <- apply(sample2.dt,1,sum)
  # # Counts of overlapping to get correct multiplier
  # nB2 <- apply(sample.B.2.dt,1,sum)
  # # Get correct multiplier estimate
  # p_m_unlink <- nB2/n2
  # # Post-sampling adjustment of the multiplier
  # adj_p_m_unlink <- est_p1_sub * p_m_unlink
  # # MBM estimate: Undercounted Benchmark divided by adjusted multiplier
  # Nhat_AdjMBM <- n1/adj_p_m_unlink
  # # Summary
  # MSE.Nhat_AdjMBM <- mean((Nhat_AdjMBM-pplsize)^2)
  # Mean.Nhat_AdjMBM <- mean(Nhat_AdjMBM)
  # SD.Nhat_AdjMBM <- sd(Nhat_AdjMBM)
  # LB95.Nhat_AdjMBM <- quantile(Nhat_AdjMBM, 0.025) %>% as.numeric()
  # UB95.Nhat_AdjMBM <- quantile(Nhat_AdjMBM, 0.975) %>% as.numeric()
  
  
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
  
  t.test.2side_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize)
  pvalue.t.2side_MBM_lk_mle <- round(t.test.2side_MBM_lk_mle$p.value, digits = 5)
  t.test.gtr_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_lk_mle <- round(t.test.gtr_MBM_lk_mle$p.value, digits = 5)
  t.test.less_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_lk_mle <- round(t.test.less_MBM_lk_mle$p.value, digits = 5)
  
  t.test.2side_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize)
  pvalue.t.2side_MBM_lk_mvue <- round(t.test.2side_MBM_lk_mvue$p.value, digits = 5)
  t.test.gtr_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_lk_mvue <- round(t.test.gtr_MBM_lk_mvue$p.value, digits = 5)
  t.test.less_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_lk_mvue <- round(t.test.less_MBM_lk_mvue$p.value, digits = 5)
  
  t.test.2side_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize)
  pvalue.t.2side_MBM_unlk_mle <- round(t.test.2side_MBM_unlk_mle$p.value, digits = 5)
  t.test.gtr_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_unlk_mle <- round(t.test.gtr_MBM_unlk_mle$p.value, digits = 5)
  t.test.less_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_unlk_mle <- round(t.test.less_MBM_unlk_mle$p.value, digits = 5)
  
  t.test.2side_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize)
  pvalue.t.2side_MBM_unlk_mvue <- round(t.test.2side_MBM_unlk_mvue$p.value, digits = 5)
  t.test.gtr_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_unlk_mvue <- round(t.test.gtr_MBM_unlk_mvue$p.value, digits = 5)
  t.test.less_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_unlk_mvue <- round(t.test.less_MBM_unlk_mvue$p.value, digits = 5)
  
  
  # Combine output
  df <- data.frame(# LP estimator
                    MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP,
                    pvalue.t.2side_LP, pvalue.t.gtr_LP, pvalue.t.less_LP, Mean.Nhat_LP_Min.m, 
                    Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
                    # Chapman estimator
                    MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap,
                    pvalue.t.2side_Chap, pvalue.t.gtr_Chap, pvalue.t.less_Chap, Mean.Nhat_Chap_Min.m,
                    Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
                    
                    # MBM_lk_mle estimator (linked)
                    MSE.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle, Med.Nhat_MBM_lk_mle, SD.Nhat_MBM_lk_mle, LB95.Nhat_MBM_lk_mle, UB95.Nhat_MBM_lk_mle, Min.Nhat_MBM_lk_mle, Max.Nhat_MBM_lk_mle,
                    pvalue.t.2side_MBM_lk_mle, pvalue.t.gtr_MBM_lk_mle, pvalue.t.less_MBM_lk_mle, Mean.Nhat_MBM_lk_mle_Min.m, 
                    Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle,
                    # MBM_direct estimator (linked)
                    MSE.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue, Med.Nhat_MBM_lk_mvue, SD.Nhat_MBM_lk_mvue, LB95.Nhat_MBM_lk_mvue, UB95.Nhat_MBM_lk_mvue, Min.Nhat_MBM_lk_mvue, Max.Nhat_MBM_lk_mvue,
                    pvalue.t.2side_MBM_lk_mvue, pvalue.t.gtr_MBM_lk_mvue, pvalue.t.less_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue_Min.m,
                    Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue,
                    
                    # MBM_unlk_mle estimator (unlinked)
                    MSE.Nhat_MBM_unlk_mle, Mean.Nhat_MBM_unlk_mle, Med.Nhat_MBM_unlk_mle, SD.Nhat_MBM_unlk_mle, LB95.Nhat_MBM_unlk_mle, UB95.Nhat_MBM_unlk_mle, Min.Nhat_MBM_unlk_mle, Max.Nhat_MBM_unlk_mle,
                    pvalue.t.2side_MBM_unlk_mle, pvalue.t.gtr_MBM_unlk_mle, pvalue.t.less_MBM_unlk_mle, 
                    Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle,
                    
                    # MBM_unlk_mvue estimator (unlinked)
                    MSE.Nhat_MBM_unlk_mvue, Mean.Nhat_MBM_unlk_mvue, Med.Nhat_MBM_unlk_mvue, SD.Nhat_MBM_unlk_mvue, LB95.Nhat_MBM_unlk_mvue, UB95.Nhat_MBM_unlk_mvue, Min.Nhat_MBM_unlk_mvue, Max.Nhat_MBM_unlk_mvue,
                    pvalue.t.2side_MBM_unlk_mvue, pvalue.t.gtr_MBM_unlk_mvue, pvalue.t.less_MBM_unlk_mvue, 
                    Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue,
                    
                    # counts of overlap
                    True.m.crm, Mean.m.crm, Med.m.crm, SD.m.crm, Min.m.crm, Max.m.crm,
                    Mean.m.mbm.unlk, Med.m.mbm.unlk, SD.m.mbm.unlk, Min.m.mbm.unlk, Max.m.mbm.unlk
                    )
  
  return(df)
}


s3.1.vary_p <- function(simsize = nsim, pplsize = targetN, n.scenario = nrow(scenarios1), 
                      pB = scenarios1$pB_vec, p1_sub = scenarios1$p1_sub_vec, p2 = scenarios1$p2_vec){
  # Start changing scenarios1 across the vector of pre-specified p
  # Create list, each element is one sub-scenario 
  s2.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- n.scenario
  for(l in 1:length.ls){
    s2.ls[[l]] <-  s3.1.single_p(simsize = simsize, pplsize = pplsize, pB = pB[l],  p1_sub = p1_sub[l], p2 = p2[l])
  }
  # Collapse the list
  s2.comb <- do.call(rbind,s2.ls)
  # Add scenarios2
  scenario <- data.frame(pB, p1_sub, p2)
  # Final output
  final.s2.comb <- as.data.frame(cbind(scenario,s2.comb))
  return(final.s2.comb)
}



#######################################################
# Function to plot CRM vs MBM on different estimators #
#######################################################

simu3.1.plot.CRM.MBM.diff.Est <- function(dt.combo = Simu3.s1.Nhat, p1_pB.i = pB.vec[1], p1_sub.i = p1_sub.vec[1]){
  
  CRM.LP <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
  CRM.LP$Methods <- rep("CRM_LP", nrow(CRM.LP))
  colnames(CRM.LP) <- c("pB", "p1_sub", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  CRM.Chap <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
  CRM.Chap$Methods <- rep("CRM_Chap", nrow(CRM.Chap))
  colnames(CRM.Chap) <- c("pB", "p1_sub", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.lk_mle <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2, Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle))
  MBM.lk_mle$Methods <- rep("MBM_link_mle", nrow(MBM.lk_mle))
  colnames(MBM.lk_mle) <- c("pB", "p1_sub", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  MBM.lk_mvue <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2, Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue))
  MBM.lk_mvue$Methods <- rep("MBM_link_mvue", nrow(MBM.lk_mvue))
  colnames(MBM.lk_mvue) <- c("pB", "p1_sub", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mle <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2, Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle))
  MBM.unlk_mle$Methods <- rep("MBM_unlk_mle", nrow(MBM.unlk_mle))
  colnames(MBM.unlk_mle) <- c("pB", "p1_sub", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mvue <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2, Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue))
  MBM.unlk_mvue$Methods <- rep("MBM_unlk_mvue", nrow(MBM.unlk_mvue))
  colnames(MBM.unlk_mvue) <- c("pB", "p1_sub", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, MBM.lk_mle, MBM.lk_mvue, MBM.unlk_mle, MBM.unlk_mvue)) 
  
  p <- ggplot(CRM.MBM.given.p1, aes(x = p2, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
    geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
    geom_line(aes(y = Rel.Bias.Nhat), size = 0.5) +
    geom_point(aes(y = Rel.Bias.Nhat)) +
    geom_hline(aes(yintercept=1), color = "black",
               linetype="dashed", size=0.5) +
    # scale_y_continuous(limits = c(0.4, 2.5), breaks = seq(0, 100, 10)) +
    ylab("Nhat/N") + 
    xlab("Sampling Probability for second sample (p2)") +
    ggtitle(paste("(Sim3) CRM vs MBM--Plot for", "pB =", p1_pB.i, "pB_sub =", p1_sub.i, sep = " ")) +
    theme_classic() +
    theme(plot.title = element_text(size = 14),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = "bottom") 
  return(p)     
}


######################################################################################
# Scenario 2 variant - Incomplete Benchmark + biased sample 2 (not a random sample 2)
######################################################################################
s3.2.single_p <- function(simsize = nsim, pplsize = targetN, 
                        pB = scenarios2$pB_vec[1], p1_sub = scenarios2$p1_sub_vec[1],
                        p2.X1 = scenarios2$p2.X1_vec[1], p2.X0 = scenarios2$p2.X0_vec[1], p.X1 = scenarios2$p.X1_vec[1]){
  # Create datasets
  sample1.dt <- matrix(0, nrow = simsize, ncol = pplsize) 
  sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  
  set.seed(1234)
  # Label Covariate X
  covariate.vec <- rbinom(pplsize, 1, p.X1)
  covariate.dt <- t(replicate(simsize, covariate.vec))
  n.X0 <- sum(covariate.vec == 0)
  n.X1 <- pplsize - n.X0
  
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
    # Draw individuals for source 2 for CRM or MBM to generate proportion given X status
    sample2.dt[i, which(covariate.vec == 0)] <- rbinom(n.X0, 1, p2.X0)
    sample2.dt[i, which(covariate.vec == 1)] <- rbinom(n.X1, 1, p2.X1)
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
  True.m.crm <- pplsize * pB * p1_sub * p2.X1
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
  
  
  # /// MBM (Linked, wrong benchmark, wrong multiplier) /// 
  # Get multiplier
  p_m_link <- m/n2
  # 1. MBM MLE
  Nhat_MBM_lk_mle <- n1/p_m_link
  # 2. MBM MVUE
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
  # Estimate
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
  
  t.test.2side_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize)
  pvalue.t.2side_MBM_lk_mle <- round(t.test.2side_MBM_lk_mle$p.value, digits = 5)
  t.test.gtr_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_lk_mle <- round(t.test.gtr_MBM_lk_mle$p.value, digits = 5)
  t.test.less_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_lk_mle <- round(t.test.less_MBM_lk_mle$p.value, digits = 5)
  
  t.test.2side_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize)
  pvalue.t.2side_MBM_lk_mvue <- round(t.test.2side_MBM_lk_mvue$p.value, digits = 5)
  t.test.gtr_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_lk_mvue <- round(t.test.gtr_MBM_lk_mvue$p.value, digits = 5)
  t.test.less_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_lk_mvue <- round(t.test.less_MBM_lk_mvue$p.value, digits = 5)
  
  t.test.2side_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize)
  pvalue.t.2side_MBM_unlk_mle <- round(t.test.2side_MBM_unlk_mle$p.value, digits = 5)
  t.test.gtr_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_unlk_mle <- round(t.test.gtr_MBM_unlk_mle$p.value, digits = 5)
  t.test.less_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_unlk_mle <- round(t.test.less_MBM_unlk_mle$p.value, digits = 5)
  
  t.test.2side_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize)
  pvalue.t.2side_MBM_unlk_mvue <- round(t.test.2side_MBM_unlk_mvue$p.value, digits = 5)
  t.test.gtr_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize, alternative = "greater")
  pvalue.t.gtr_MBM_unlk_mvue <- round(t.test.gtr_MBM_unlk_mvue$p.value, digits = 5)
  t.test.less_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize, alternative = "less")
  pvalue.t.less_MBM_unlk_mvue <- round(t.test.less_MBM_unlk_mvue$p.value, digits = 5)
  
  # Combine output
  df <- data.frame(# LP estimator
    MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP,
    pvalue.t.2side_LP, pvalue.t.gtr_LP, pvalue.t.less_LP, Mean.Nhat_LP_Min.m, 
    Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
    # Chapman estimator
    MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap,
    pvalue.t.2side_Chap, pvalue.t.gtr_Chap, pvalue.t.less_Chap, Mean.Nhat_Chap_Min.m,
    Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
    
    # MBM_lk_mle estimator (linked)
    MSE.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle, Med.Nhat_MBM_lk_mle, SD.Nhat_MBM_lk_mle, LB95.Nhat_MBM_lk_mle, UB95.Nhat_MBM_lk_mle, Min.Nhat_MBM_lk_mle, Max.Nhat_MBM_lk_mle,
    pvalue.t.2side_MBM_lk_mle, pvalue.t.gtr_MBM_lk_mle, pvalue.t.less_MBM_lk_mle, Mean.Nhat_MBM_lk_mle_Min.m, 
    Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle,
    # MBM_lk_mvue estimator (linked)
    MSE.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue, Med.Nhat_MBM_lk_mvue, SD.Nhat_MBM_lk_mvue, LB95.Nhat_MBM_lk_mvue, UB95.Nhat_MBM_lk_mvue, Min.Nhat_MBM_lk_mvue, Max.Nhat_MBM_lk_mvue,
    pvalue.t.2side_MBM_lk_mvue, pvalue.t.gtr_MBM_lk_mvue, pvalue.t.less_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue_Min.m,
    Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue,
    
    # MBM_unlk_mle estimator (unlinked)
    MSE.Nhat_MBM_unlk_mle, Mean.Nhat_MBM_unlk_mle, Med.Nhat_MBM_unlk_mle, SD.Nhat_MBM_unlk_mle, LB95.Nhat_MBM_unlk_mle, UB95.Nhat_MBM_unlk_mle, Min.Nhat_MBM_unlk_mle, Max.Nhat_MBM_unlk_mle,
    pvalue.t.2side_MBM_unlk_mle, pvalue.t.gtr_MBM_unlk_mle, pvalue.t.less_MBM_unlk_mle, 
    Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle,
    
    # MBM_unlk_mvue estimator (unlinked)
    MSE.Nhat_MBM_unlk_mvue, Mean.Nhat_MBM_unlk_mvue, Med.Nhat_MBM_unlk_mvue, SD.Nhat_MBM_unlk_mvue, LB95.Nhat_MBM_unlk_mvue, UB95.Nhat_MBM_unlk_mvue, Min.Nhat_MBM_unlk_mvue, Max.Nhat_MBM_unlk_mvue,
    pvalue.t.2side_MBM_unlk_mvue, pvalue.t.gtr_MBM_unlk_mvue, pvalue.t.less_MBM_unlk_mvue, 
    Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue,
    
    # counts of overlap
    True.m.crm, Mean.m.crm, Med.m.crm, SD.m.crm, Min.m.crm, Max.m.crm,
    Mean.m.mbm.unlk, Med.m.mbm.unlk, SD.m.mbm.unlk, Min.m.mbm.unlk, Max.m.mbm.unlk
  )
  
  
  return(df)
}

# Scenario 2 variant - sample 2 is not random
s3.2.vary_p <- function(simsize = nsim, pplsize = targetN, n.scenario = nrow(scenarios2), 
                      pB = scenarios2$pB_vec, p1_sub = scenarios2$p1_sub_vec, 
                      p2.X1 = scenarios2$p2.X1_vec, p2.X0 = scenarios2$p2.X0_vec, p.X1 = scenarios2$p.X1_vec){
  # Start changing scenarios2 across the vector of pre-specified p
  # Create list, each element is one sub-scenario 
  s3.2.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- n.scenario
  for(l in 1:length.ls){
    s3.2.ls[[l]] <-  s3.2.single_p(simsize = simsize, pplsize = pplsize, pB = pB[l],  p1_sub = p1_sub[l], 
                                   p2.X1 = p2.X1[l], p2.X0 = p2.X0[l], p.X1 = p.X1[l])
  }
  # Collapse the list
  s3.2.comb <- do.call(rbind,s3.2.ls)
  # Add scenarios2
  scenario <- data.frame(p.X1, pB, p1_sub, p2.X1, p2.X0)
  # Final output
  final.s3.2.comb <- as.data.frame(cbind(scenario,s3.2.comb))
  return(final.s3.2.comb)
}



#######################################################
# Function to plot CRM vs MBM on different estimators #
#######################################################

simu3.2.plot.CRM.MBM.diff.Est <- function(dt.combo = Simu3.s2.Nhat, p1_pB.i = scenarios2$pB_vec[1], p1_sub.i = scenarios2$p1_sub_vec[1]){
  
  CRM.LP <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2.X1, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
  CRM.LP$Methods <- rep("CRM_LP", nrow(CRM.LP))
  colnames(CRM.LP) <- c("pB", "p1_sub", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  CRM.Chap <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2.X1, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
  CRM.Chap$Methods <- rep("CRM_Chap", nrow(CRM.Chap))
  colnames(CRM.Chap) <- c("pB", "p1_sub", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.lk_mle <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2.X1, Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle))
  MBM.lk_mle$Methods <- rep("MBM_link_mle", nrow(MBM.lk_mle))
  colnames(MBM.lk_mle) <- c("pB", "p1_sub", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  MBM.lk_mvue <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2.X1, Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue))
  MBM.lk_mvue$Methods <- rep("MBM_link_mvue", nrow(MBM.lk_mvue))
  colnames(MBM.lk_mvue) <- c("pB", "p1_sub", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mle <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2.X1, Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle))
  MBM.unlk_mle$Methods <- rep("MBM_unlk_mle", nrow(MBM.unlk_mle))
  colnames(MBM.unlk_mle) <- c("pB", "p1_sub", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  MBM.unlk_mvue <- subset(dt.combo, pB == p1_pB.i & p1_sub == p1_sub.i, select = c(pB, p1_sub, p2.X1, Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue))
  MBM.unlk_mvue$Methods <- rep("MBM_unlk_mvue", nrow(MBM.unlk_mvue))
  colnames(MBM.unlk_mvue) <- c("pB", "p1_sub", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
  
  CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, MBM.lk_mle, MBM.lk_mvue, MBM.unlk_mle, MBM.unlk_mvue)) # 
  
  p <- ggplot(CRM.MBM.given.p1, aes(x = p2.X1, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
    geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
    geom_line(aes(y = Rel.Bias.Nhat), size = 0.5) +
    geom_point(aes(y = Rel.Bias.Nhat)) +
    geom_hline(aes(yintercept=1), color = "black",
               linetype="dashed", size=0.5) +
    ylab("Nhat/N") + 
    xlab("Sampling Probability for second sample (p2 given X=1)") +
    ggtitle(paste("(Sim3.2) CRM vs MBM--Plot for", "pB =", p1_pB.i, "pB_sub =", p1_sub.i, sep = " ")) +
    theme_classic() +
    theme(plot.title = element_text(size = 14),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) 
  return(p)     
}





##############################################################################
######### Scenario 3: Incomplete Benchmark + dependent sample 2 ##############
##############################################################################
s3.3.single_p <- function(simsize = nsim, pplsize = targetN, 
                          pB = scenarios3$pB_vec[1], p1_sub = scenarios3$p1_sub_vec[1], p2 = scenarios3$p2_vec[1],
                          dep.factor = scenarios3$dep.factor_vec[1]){
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
  
  # # To analyze the simulation results
  # # One-sample t test 
  # t.test.2side_LP <- t.test(Nhat_LP, mu = pplsize)
  # pvalue.t.2side_LP <- round(t.test.2side_LP$p.value, digits = 5)
  # t.test.gtr_LP <- t.test(Nhat_LP, mu = pplsize, alternative = "greater")
  # pvalue.t.gtr_LP <- round(t.test.gtr_LP$p.value, digits = 5)
  # t.test.less_LP <- t.test(Nhat_LP, mu = pplsize, alternative = "less")
  # pvalue.t.less_LP <- round(t.test.less_LP$p.value, digits = 5)
  # 
  # t.test.2side_Chap <- t.test(Nhat_Chap, mu = pplsize)
  # pvalue.t.2side_Chap <- round(t.test.2side_Chap$p.value, digits = 5)
  # t.test.gtr_Chap <- t.test(Nhat_Chap, mu = pplsize, alternative = "greater")
  # pvalue.t.gtr_Chap <- round(t.test.gtr_Chap$p.value, digits = 5)
  # t.test.less_Chap <- t.test(Nhat_Chap, mu = pplsize, alternative = "less")
  # pvalue.t.less_Chap <- round(t.test.less_Chap$p.value, digits = 5)
  # 
  # t.test.2side_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize)
  # pvalue.t.2side_MBM_lk_mle <- round(t.test.2side_MBM_lk_mle$p.value, digits = 5)
  # t.test.gtr_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize, alternative = "greater")
  # pvalue.t.gtr_MBM_lk_mle <- round(t.test.gtr_MBM_lk_mle$p.value, digits = 5)
  # t.test.less_MBM_lk_mle <- t.test(Nhat_MBM_lk_mle, mu = pplsize, alternative = "less")
  # pvalue.t.less_MBM_lk_mle <- round(t.test.less_MBM_lk_mle$p.value, digits = 5)
  # 
  # t.test.2side_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize)
  # pvalue.t.2side_MBM_lk_mvue <- round(t.test.2side_MBM_lk_mvue$p.value, digits = 5)
  # t.test.gtr_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize, alternative = "greater")
  # pvalue.t.gtr_MBM_lk_mvue <- round(t.test.gtr_MBM_lk_mvue$p.value, digits = 5)
  # t.test.less_MBM_lk_mvue <- t.test(Nhat_MBM_lk_mvue, mu = pplsize, alternative = "less")
  # pvalue.t.less_MBM_lk_mvue <- round(t.test.less_MBM_lk_mvue$p.value, digits = 5)
  # 
  # t.test.2side_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize)
  # pvalue.t.2side_MBM_unlk_mle <- round(t.test.2side_MBM_unlk_mle$p.value, digits = 5)
  # t.test.gtr_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize, alternative = "greater")
  # pvalue.t.gtr_MBM_unlk_mle <- round(t.test.gtr_MBM_unlk_mle$p.value, digits = 5)
  # t.test.less_MBM_unlk_mle <- t.test(Nhat_MBM_unlk_mle, mu = pplsize, alternative = "less")
  # pvalue.t.less_MBM_unlk_mle <- round(t.test.less_MBM_unlk_mle$p.value, digits = 5)
  # 
  # t.test.2side_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize)
  # pvalue.t.2side_MBM_unlk_mvue <- round(t.test.2side_MBM_unlk_mvue$p.value, digits = 5)
  # t.test.gtr_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize, alternative = "greater")
  # pvalue.t.gtr_MBM_unlk_mvue <- round(t.test.gtr_MBM_unlk_mvue$p.value, digits = 5)
  # t.test.less_MBM_unlk_mvue <- t.test(Nhat_MBM_unlk_mvue, mu = pplsize, alternative = "less")
  # pvalue.t.less_MBM_unlk_mvue <- round(t.test.less_MBM_unlk_mvue$p.value, digits = 5)
  # 
  
  # Combine output
  df <- data.frame(# LP estimator
    MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP, Mean.Nhat_LP_Min.m, 
    # pvalue.t.2side_LP, pvalue.t.gtr_LP, pvalue.t.less_LP, 
    Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
    # Chapman estimator
    MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap, Mean.Nhat_Chap_Min.m,
    # pvalue.t.2side_Chap, pvalue.t.gtr_Chap, pvalue.t.less_Chap, 
    Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
    
    # MBM_lk_mle estimator (linked)
    MSE.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle, Med.Nhat_MBM_lk_mle, SD.Nhat_MBM_lk_mle, LB95.Nhat_MBM_lk_mle, UB95.Nhat_MBM_lk_mle, Min.Nhat_MBM_lk_mle, Max.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle_Min.m, 
    # pvalue.t.2side_MBM_lk_mle, pvalue.t.gtr_MBM_lk_mle, pvalue.t.less_MBM_lk_mle, 
    Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle,
    # MBM_direct estimator (linked)
    MSE.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue, Med.Nhat_MBM_lk_mvue, SD.Nhat_MBM_lk_mvue, LB95.Nhat_MBM_lk_mvue, UB95.Nhat_MBM_lk_mvue, Min.Nhat_MBM_lk_mvue, Max.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue_Min.m,
    # pvalue.t.2side_MBM_lk_mvue, pvalue.t.gtr_MBM_lk_mvue, pvalue.t.less_MBM_lk_mvue, 
    Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue,
    
    # MBM_unlk_mle estimator (unlinked)
    MSE.Nhat_MBM_unlk_mle, Mean.Nhat_MBM_unlk_mle, Med.Nhat_MBM_unlk_mle, SD.Nhat_MBM_unlk_mle, LB95.Nhat_MBM_unlk_mle, UB95.Nhat_MBM_unlk_mle, Min.Nhat_MBM_unlk_mle, Max.Nhat_MBM_unlk_mle,
    # pvalue.t.2side_MBM_unlk_mle, pvalue.t.gtr_MBM_unlk_mle, pvalue.t.less_MBM_unlk_mle, 
    Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle,
    
    # MBM_unlk_mvue estimator (unlinked)
    MSE.Nhat_MBM_unlk_mvue, Mean.Nhat_MBM_unlk_mvue, Med.Nhat_MBM_unlk_mvue, SD.Nhat_MBM_unlk_mvue, LB95.Nhat_MBM_unlk_mvue, UB95.Nhat_MBM_unlk_mvue, Min.Nhat_MBM_unlk_mvue, Max.Nhat_MBM_unlk_mvue,
    # pvalue.t.2side_MBM_unlk_mvue, pvalue.t.gtr_MBM_unlk_mvue, pvalue.t.less_MBM_unlk_mvue, 
    Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue,
    
    # counts of overlap
    True.m.crm, Mean.m.crm, Med.m.crm, SD.m.crm, Min.m.crm, Max.m.crm,
    Mean.m.mbm.unlk, Med.m.mbm.unlk, SD.m.mbm.unlk, Min.m.mbm.unlk, Max.m.mbm.unlk
  )
  
  return(df)
}


s3.3.vary_p <- function(simsize = nsim, pplsize = targetN, n.scenario = nrow(scenarios3), 
                        pB = scenarios3$pB_vec, p1_sub = scenarios3$p1_sub_vec, p2 = scenarios3$p2_vec,
                        dep.factor = scenarios3$dep.factor_vec){
  # Start changing scenarios3 across the vector of pre-specified p
  # Create list, each element is one sub-scenario 
  s3.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- n.scenario
  for(l in 1:length.ls){
    s3.ls[[l]] <-  s3.3.single_p(simsize = simsize, pplsize = pplsize, pB = pB[l],  p1_sub = p1_sub[l], 
                                 p2 = p2[l],
                                 dep.factor = dep.factor[l])
  }
  # Collapse the list
  s3.comb <- do.call(rbind,s3.ls)
  # Add scenarios3
  scenario <- data.frame(pB, p1_sub, p2, dep.factor)
  # Final output
  final.s3.comb <- as.data.frame(cbind(scenario,s3.comb))
  return(final.s3.comb)
}



#######################################################
# Function to plot CRM vs MBM on different estimators #
#######################################################

simu3.3.Dt.CRM.MBM.diff.Est <- function(dt.combo = Simu3.s3.Nhat, p1_pB.i = pB.vec[1], p1_sub.i = p1_sub.vec[1], dep.factor.i = dep.factor.vec[1]){
  
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
  
simu3.3.plot.CRM.MBM.diff.Est <- function(dt.combo = subdt, p1_pB.i = pB.vec[1], p1_sub.i = p1_sub.vec[1], dep.factor.i = dep.factor.vec[1], ErrorBar = "Yes", overORunderEst = "Over"){
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
      ggtitle(paste("(Sim3.3) CRM vs MBM--Plot for", "pB =", p1_pB.i, "pB_sub =", p1_sub.i, "dep.level =", dep.factor.i, sep = " ")) +
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
      # Force all plots have the same limits and breaks
      # scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0.4, 1.4, by = 0.05)) +
      ylab("Nhat/N") + 
      xlab("Sampling Probability for second sample (p2)") +
      ggtitle(paste("(Sim3.3) CRM vs MBM--Plot for", "pB =", p1_pB.i, "pB_sub =", p1_sub.i, "dep.level =", dep.factor.i, sep = " ")) +
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









######################## xxxxxxxxxxxxxxxxxxxxx NOT RUN xxxxxxxxxxxxxxxxxxxxx ########################
# Link/Not link samples in MBM and Incomplete Benchmark and biased multiplier (to correct incompleteness issue in Benchmark)
# Scernario 3
# s3.single_p <- function(simsize = nsim, pplsize = targetN, 
#                         pB = scenarios3$pB_vec[1], p1_sub = scenarios3$p1_sub_vec[1], p2 = scenarios3$p2_vec[3]){
#   # Create datasets
#   benchmark <- matrix(0, nrow = simsize, ncol = pplsize)
#   sample1.dt <- matrix(0, nrow = simsize, ncol = pplsize) 
#   sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
#   # Simulate for each element
#   for (i in 1:simsize){
#     set.seed(1234+i)
#     # Draw individuals for Benchmark
#     benchmark[i,] <- rbinom(pplsize, 1, pB)
#     n_B <- sum(benchmark[i,])
#     n_not_B <- pplsize-n_B
#     # Incompletely select Benchmark people into samples
#     sample1.dt[i,which(benchmark[i,] == 1)] <- rbinom(n_B, 1, p1_sub)
#     sample1.dt[i,which(benchmark[i,] == 0)] <- rbinom(n_not_B, 1, 0)
#     # Draw individuals for source 2 for CRM or MBM to generate proportion
#     # sample1 counts
#     n_B_s1 <- sum(sample1.dt[i,] == 1 & benchmark[i,] == 1)
#     n_miss_B_s1 <- n_B-n_B_s1
#     # Unbiased sample 2 (??????????)
#     sample2.dt[i,which(sample1.dt[i,] == 1 & benchmark[i,] == 1)] <- rbinom(n_B_s1, 1, p2*p1_sub)
#     sample2.dt[i,which(sample1.dt[i,] == 0 & benchmark[i,] == 1)] <- rbinom(n_miss_B_s1, 1, p2*(1-p1_sub))
#     sample2.dt[i,which(benchmark[i,] == 0)] <- rbinom(n_not_B, 1, p2)
#   }
#   # Got two data sample, sample1.dt and sample2.dt
#   ####################
#   # --- Estimate --- #
#   ####################
#   # /// CRM /// #
#   # Overlapping individuals, equivalent to labeling the overlapping
#   sample1.2.dt <- sample1.dt*sample2.dt
#   # Marginal counts for sources
#   n1 <- apply(sample1.dt,1,sum)
#   n2 <- apply(sample2.dt,1,sum)
#   # Counts of overlapping
#   n11 <- apply(sample1.2.dt,1,sum)
#   # LP estimator 
#   Nhat_CRM <- n1*n2/n11
#   # Summary
#   avg_n11_CRM <- mean(n11)
#   avg_n1_CRM <- mean(n1)
#   avg_n2_CRM <- mean(n2)
#   # Statistics
#   MSE.Nhat_CRM <- mean((Nhat_CRM-pplsize)^2)
#   Mean.Nhat_CRM <- mean(Nhat_CRM)
#   SD.Nhat_CRM <- sd(Nhat_CRM)
#   LB95.Nhat_CRM <- quantile(Nhat_CRM, 0.975)
#   UB95.Nhat_CRM <- quantile(Nhat_CRM, 0.025)
#   
#   # /// MBM (Linked, wrong benchmark, wrong multiplier) /// #
#   # Get multiplier
#   p_m_link <- n11/n2
#   # MBM estimator
#   Nhat_MBM_link <- n1/p_m_link
#   # Summary
#   avg_multiplier_link <- mean(p_m_link)
#   avg_n11_MBM_link <- mean(n11)
#   avg_n1_MBM_link <- mean(n1)
#   avg_n2_MBM_link <- mean(n2)
#   # Statistics
#   MSE.Nhat_MBM_link <- mean((Nhat_MBM_link-pplsize)^2)
#   Mean.Nhat_MBM_link <- mean(Nhat_MBM_link)
#   SD.Nhat_MBM_link <- sd(Nhat_MBM_link)
#   LB95.Nhat_MBM_link <- quantile(Nhat_MBM_link, 0.975)
#   UB95.Nhat_MBM_link <- quantile(Nhat_MBM_link, 0.025)
#   
#   # /// MBM - (Unlinked, wrong benchmark, true multiplier) /// #
#   # Find correct number of overlapping
#   sample.B.2.dt <- benchmark*sample2.dt
#   # Marginal counts for sources
#   n1 <- apply(sample1.dt,1,sum)
#   n2 <- apply(sample2.dt,1,sum)
#   # Counts of overlapping to get correct multiplier
#   nB2 <- apply(sample.B.2.dt,1,sum)
#   # Get correct multiplier estimate
#   p_m_unlink <- nB2/n2
#   # MBM estimate: Undercounted Benchmark divided by true multiplier
#   Nhat_MBM <- n1/p_m_unlink
#   # Summary
#   avg_multiplier_MBM <- mean(p_m_unlink)
#   avg_nB2_MBM <- mean(nB2)
#   avg_n1_MBM<- mean(n1)
#   avg_n2_MBM <- mean(n2)
#   # Statistics
#   MSE.Nhat_MBM <- mean((Nhat_MBM-pplsize)^2)
#   Mean.Nhat_MBM <- mean(Nhat_MBM)
#   SD.Nhat_MBM <- sd(Nhat_MBM)
#   LB95.Nhat_MBM <- quantile(Nhat_MBM, 0.975)
#   UB95.Nhat_MBM <- quantile(Nhat_MBM, 0.025)
#   
#   # /// MBM (Unlinked) Undercounted Benchmark with post-sampling adjustment of multiplier /// #
#   # Find adjustment of the multiplier
#   # 1. Make sure sample2 is a random sample - Checked
#   # 2. Get the information of % of undercounting in Benchmark (deflation factor)
#   est_p1_sub <- c()
#   for(i in 1:simsize){
#     ith.B.sample1.overlap <- sum(benchmark[i,] == 1 & sample1.dt[i,] == 0)
#     ith.B <- sum(benchmark[i,])
#     ith.est_p1_sub <- 1-ith.B.sample1.overlap/ith.B
#     est_p1_sub <- c(est_p1_sub, ith.est_p1_sub)
#   }
#   # Find correct number of overlapping
#   sample.B.2.dt <- benchmark*sample2.dt
#   # Marginal counts for sources
#   n1 <- apply(sample1.dt,1,sum)
#   n2 <- apply(sample2.dt,1,sum)
#   # Counts of overlapping to get correct multiplier
#   nB2 <- apply(sample.B.2.dt,1,sum)
#   # Get correct multiplier estimate
#   p_m_unlink <- nB2/n2
#   # Post-sampling adjustment of the multiplier
#   adj_p_m_unlink <- est_p1_sub * p_m_unlink
#   # MBM estimate: Undercounted Benchmark divided by adjusted multiplier
#   Nhat_AdjMBM <- n1/adj_p_m_unlink
#   # Summary
#   avg_multiplier_AdjMBM <- mean(adj_p_m_unlink)
#   avg_nB2_AdjMBM <- mean(nB2)
#   avg_n1_AdjMBM<- mean(n1)
#   avg_n2_AdjMBM <- mean(n2)
#   # Statistics
#   MSE.Nhat_AdjMBM <- mean((Nhat_AdjMBM-pplsize)^2)
#   Mean.Nhat_AdjMBM <- mean(Nhat_AdjMBM)
#   SD.Nhat_AdjMBM <- sd(Nhat_AdjMBM)
#   LB95.Nhat_AdjMBM <- quantile(Nhat_AdjMBM, 0.975)
#   UB95.Nhat_AdjMBM <- quantile(Nhat_AdjMBM, 0.025)
#   
#   # Combine both estimates
#   df <- data.frame(# CRM
#     avg_n11_CRM, avg_n1_CRM, avg_n2_CRM, 
#     MSE.Nhat_CRM, Mean.Nhat_CRM, SD.Nhat_CRM, LB95.Nhat_CRM, UB95.Nhat_CRM, 
#     # MBM (linked)
#     avg_multiplier_link, avg_n11_MBM_link, avg_n1_MBM_link, avg_n2_MBM_link, 
#     MSE.Nhat_MBM_link, Mean.Nhat_MBM_link, SD.Nhat_MBM_link, LB95.Nhat_MBM_link, UB95.Nhat_MBM_link, 
#     # MBM (unlinked)
#     avg_multiplier_MBM, avg_nB2_MBM, avg_n1_MBM, avg_n2_MBM, 
#     MSE.Nhat_MBM, Mean.Nhat_MBM, SD.Nhat_MBM, LB95.Nhat_MBM, UB95.Nhat_MBM,
#     # MBM (unlinked)
#     avg_multiplier_AdjMBM, avg_nB2_AdjMBM, avg_n1_AdjMBM, avg_n2_AdjMBM, 
#     MSE.Nhat_AdjMBM, Mean.Nhat_AdjMBM, SD.Nhat_AdjMBM, LB95.Nhat_AdjMBM, UB95.Nhat_AdjMBM
#   )
#   df <- round(df, digits = 3)
#   return(df)
# }
# 
# 
# s3.vary_p <- function(simsize = nsim, pplsize = targetN, n.scenario = nrow(scenarios3), 
#                       pB = scenarios3$pB_vec, p1_sub = scenarios3$p1_sub_vec, p2 = scenarios3$p2_vec){
#   # Start changing scenarios3 across the vector of pre-specified p
#   # Create list, each element is one sub-scenario 
#   s3.ls <- list()
#   # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
#   length.ls <- n.scenario
#   for(l in 1:length.ls){
#     s3.ls[[l]] <-  s3.single_p(simsize = simsize, pplsize = pplsize, pB = pB[l],  p1_sub = p1_sub[l], p2 = p2[l])
#   }
#   # Collapse the list
#   s3.comb <- do.call(rbind,s3.ls)
#   # Add scenarios3
#   scenario <- data.frame(pB, p1_sub, p2)
#   # Final output
#   final.s3.comb <- as.data.frame(cbind(scenario,s3.comb))
#   return(final.s3.comb)
# }
# 
