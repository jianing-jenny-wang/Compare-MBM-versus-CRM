# Simulation 4 - Take a benchmark(census) as a nominal sample
library(dplyr)

############################################################
######## Scenario case 1 (Sample 2 is homogeneous) #########
############################################################

#########################################
# Function to compute estimates for CRM #
#########################################

crm.1 <- function(simsize = nsim, pplsize = targetN, p.X1 = scenarios1$p.X1[1], pB.X0 = scenarios1$pB.X0[1], pB.X1 = scenarios1$pB.X1[1], p2 = scenarios1$p2[1]){
# Create datasets
benchmark.dt <- matrix(0, nrow = simsize, ncol = pplsize)
sample1.dt <- matrix(0, nrow = simsize, ncol = pplsize)
sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)

set.seed(1234)
# Label Covariate X
covariate.vec <- rbinom(pplsize, 1, p.X1)
covariate.dt <- t(replicate(simsize, covariate.vec))
n.X0 <- sum(covariate.vec == 0)
n.X1 <- pplsize - n.X0

# Simulate for each element
for (i in 1:simsize){
  set.seed(1234+i)
  # Benchmark is a random sample
  benchmark.dt[i, which(covariate.vec == 0)] <- rbinom(n.X0, 1, pB.X0)
  benchmark.dt[i, which(covariate.vec == 1)] <- rbinom(n.X1, 1, pB.X1)
  # All benchmark people got sampled into the sample 1 for CRM and MBM
  sample1.dt[i,which(benchmark.dt[i,] == 1)] <- 1
  # Draw individuals for source 2 for CRM or MBM to generate proportion
  sample2.dt[i,] <- rbinom(pplsize, 1, p2)
}
# Got two data sample, sample1.dt and sample2.dt
# Marginal counts for sources
n1 <- apply(sample1.dt,1,sum)
n2 <- apply(sample2.dt,1,sum)
# Overlapping individuals, equivalent to labeling the overlapping
sample1.2.dt<-sample1.dt*sample2.dt
# # Create contingency table
# ppl.dt <- data.frame(covariate.vec,benchmark.vec, sample1.dt[i,], sample2.dt[i,])
# colnames(ppl.dt) <- c("X","Benchmark","Sample1","Sample2")
# ppl.dt$X <- ifelse(ppl.dt$X == 1, "A", "B")
# cont.tb <- table(ppl.dt$Sample1,ppl.dt$Sample2,ppl.dt$X) # For given X, sample 1 (row) by sample 2 (col)
# cont.tb.l <- as.data.frame(cont.tb)
# colnames(cont.tb.l) <- c("Sample1","Sample2","X","Counts")
####################
# --- Estimate --- #
####################
# /// CRM /// #
# Counts of overlapping
m <- apply(sample1.2.dt,1,sum)
# LP estimator 
# Est Total N for CRM
Nhat_LP <- n1*n2/m
# Chapman Estimator for CRM
# Est Total N for CRM
Nhat_Chap <- (n1+1)*(n2+1)/(m+1) - 1
# Stratified Estimate using LP and Chapman
Nhat_LP_Strata <- c()
Nhat_Chap_Strata <- c()
for(i in 1:simsize){
  sample.1 <- sample1.dt[i,]
  sample.2 <- sample2.dt[i,]
  sample.m <- sample1.2.dt[i,]
  n1.X1 <- sum(sample.1[which(covariate.vec == 1)])
  n1.X0 <- sum(sample.1[which(covariate.vec == 0)])
  n2.X1 <- sum(sample.2[which(covariate.vec == 1)])
  n2.X0 <-  sum(sample.2[which(covariate.vec == 0)])
  m.X1 <- sum(sample.m[which(covariate.vec == 1)])
  m.X0 <- sum(sample.m[which(covariate.vec == 0)])
  # Nhat.X1_LP <-  n1.X1*n2.X1/m.X1
  # Nhat.X0_LP <-  n1.X0*n2.X0/m.X0
  # Nhat_LP.X <- Nhat.X1_LP + Nhat.X0_LP
  # Nhat_LP_Strata <- c(Nhat_LP_Strata, Nhat_LP.X )
  Nhat.X1_Chap <-  (n1.X1+1)*(n2.X1+1)/(m.X1+1) - 1
  Nhat.X0_Chap <-  (n1.X0+1)*(n2.X0+1)/(m.X0+1) - 1
  Nhat_Chap.X <- Nhat.X1_Chap + Nhat.X0_Chap
  Nhat_Chap_Strata <- c(Nhat_Chap_Strata, Nhat_Chap.X)
}
# Summary
True.m <- pplsize * (1-p.X1) * pB.X0 * p2 + pplsize * p.X1 * pB.X1 * p2
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

# Stratified
MSE.Nhat_Chap_Strata <- mean((Nhat_Chap_Strata-pplsize)^2)
Mean.Nhat_Chap_Strata <- mean(Nhat_Chap_Strata)
Med.Nhat_Chap_Strata <- median(Nhat_Chap_Strata)
SD.Nhat_Chap_Strata <- sd(Nhat_Chap_Strata)
LB95.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata, 0.025) %>% as.numeric()
UB95.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata, 0.975) %>% as.numeric()
Min.Nhat_Chap_Strata <- min(Nhat_Chap_Strata)
Max.Nhat_Chap_Strata <- max(Nhat_Chap_Strata)
Mean.Nhat_Chap_Strata_Min.m <- mean(Nhat_Chap_Strata[which(m == Min.m)])
Rel.Bias.Nhat_Chap_Strata <- mean(Nhat_Chap_Strata/pplsize)
LB95.Rel.Bias.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata/pplsize, 0.025) %>% as.numeric()
UB95.Rel.Bias.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata/pplsize, 0.975) %>% as.numeric()

# Combine output
df <- data.frame(# LP estimator
                MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP,
                Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
                # Chapman estimator
                MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap,
                Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
                # Chapman estimator with stratification
                MSE.Nhat_Chap_Strata, Mean.Nhat_Chap_Strata, Med.Nhat_Chap_Strata, SD.Nhat_Chap_Strata, LB95.Nhat_Chap_Strata, UB95.Nhat_Chap_Strata, Min.Nhat_Chap_Strata, Max.Nhat_Chap_Strata,
                Rel.Bias.Nhat_Chap_Strata, LB95.Rel.Bias.Nhat_Chap_Strata, UB95.Rel.Bias.Nhat_Chap_Strata,
                # counts of overlap
                True.m, Mean.m, Med.m, SD.m, Min.m, Max.m)
}

# Changing p vector from both low to both high
crm.1.ls <- function(simsize = nsim, pplsize = targetN, p.X1 = scenarios1$p.X1, pB.X0 = scenarios1$pB.X0, pB.X1 = scenarios1$pB.X1, p2 = scenarios1$p2){
  # Start changing scenarios across the vector of pre-specified p
  # Create list, each element is one sub-scenario of pB and p2
  s1.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- length(p2)
  for(l in 1:length.ls){
    s1.ls[[l]] <-  crm.1(simsize = simsize, pplsize = pplsize, p.X1 = p.X1[l], pB.X0 = pB.X0[l], pB.X1 = pB.X1[l], p2 = p2[l])
  
  }
  # Collapse the list
  s1.comb <- do.call(rbind,s1.ls)
  # Add scenarios
  scenario <- data.frame(p.X1, pB.X0, pB.X1, p2)
  # Final output
  final.s1.comb <- as.data.frame(cbind(scenario,s1.comb))
  return(final.s1.comb)
}


####################################################
# Function to compute estimates for MBM (Unlinked) #
####################################################

mbm.1 <- function(simsize = nsim, pplsize = targetN, p.X1 = scenarios1$p.X1[2], pB.X0 = scenarios1$pB.X0[2], pB.X1 = scenarios1$pB.X1[2], p2 = scenarios1$p2[2]){
  # Create datasets
  sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  
  set.seed(1234)
  # Label Covariate X
  covariate.vec <- rbinom(pplsize, 1, p.X1)
  covariate.dt <- t(replicate(simsize, covariate.vec))
  n.X0 <- sum(covariate.vec == 0)
  n.X1 <- pplsize - n.X0
  # Label Benchmark members
  benchmark.vec <- c()
  benchmark.vec[which(covariate.vec == 0)] <- rbinom(n.X0, 1, pB.X0)
  benchmark.vec[which(covariate.vec == 1)] <- rbinom(n.X1, 1, pB.X1)
  benchmark.dt <- t(replicate(simsize, benchmark.vec))
  # Simulate for each element
  for (i in 1:simsize){
    set.seed(1234+i)
    # Draw individuals for source 2 for CRM or MBM to generate proportion
    sample2.dt[i,] <- rbinom(pplsize, 1, p2)
  }
  # Got two data sample, sample1.dt and sample2.dt
  # Marginal counts for sources
  nB <- apply(benchmark.dt,1,sum)
  n2 <- apply(sample2.dt,1,sum)
  # Overlapping individuals, equivalent to labeling the overlapping
  sample.B.2.dt<-benchmark.dt*sample2.dt
  # Counts of overlapping
  m <- apply(sample.B.2.dt,1,sum)
  ####################
  # --- Estimate --- #
  ####################
  # /// MBM /// #
  p_m <- m/n2
  # 1. Use MLE
  Nhat_MBM_lk_mle <- nB/p_m
  # 2. Use MVUE
  Nhat_MBM_lk_mvue <- ((nB -1)/p_m) + 1
  
  # 3. Unlinked sample 2, calculate the theoretical m rather than counting the overlap
  # Calculate the theoretical number of sample 2
  sample2.X1 <- covariate.dt*sample2.dt
  n2_unlink.X1 <- apply(sample2.X1, 1, sum) # Mean = pB.X1
  # Compute the theoretical number of benchmark in sample 2
  pB.X1_nB <- unique(nB)/n.X1
  set.seed(1234)
  m_unlink <- rbinom(simsize, n2_unlink.X1, pB.X1_nB)
  # Compute the multiplier
  p_m_unlink <- m_unlink/n2
  # MLE 
  Nhat_MBM_unlk_mle <- nB/p_m_unlink
  # MVUE
  Nhat_MBM_unlk_mvue <- ((nB-1)/p_m_unlink) + 1
  
  # Summary
  True.m <- pplsize * (1-p.X1) * pB.X0 * p2 + pplsize * p.X1 * pB.X1 * p2
  Mean.m <- mean(m)
  Med.m <- median(m)
  SD.m <- sd(m)
  Min.m <- min(m)
  Max.m <- max(m)
  
  MSE.Nhat_MBM_lk_mle <- mean((Nhat_MBM_lk_mle-pplsize)^2)
  Mean.Nhat_MBM_lk_mle <- mean(Nhat_MBM_lk_mle)
  Med.Nhat_MBM_lk_mle <- median(Nhat_MBM_lk_mle)
  SD.Nhat_MBM_lk_mle <- sd(Nhat_MBM_lk_mle)
  LB95.Nhat_MBM_lk_mle <- quantile(Nhat_MBM_lk_mle, 0.025) %>% as.numeric()
  UB95.Nhat_MBM_lk_mle <- quantile(Nhat_MBM_lk_mle, 0.975) %>% as.numeric()
  Min.Nhat_MBM_lk_mle <- min(Nhat_MBM_lk_mle)
  Max.Nhat_MBM_lk_mle <- max(Nhat_MBM_lk_mle)
  Mean.Nhat_MBM_lk_mle_Min.m <- mean(Nhat_MBM_lk_mle[which(m == Min.m)])
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
  Mean.Nhat_MBM_lk_mvue_Min.m <- mean(Nhat_MBM_lk_mvue[which(m == Min.m)])
  Rel.Bias.Nhat_MBM_lk_mvue <- mean(Nhat_MBM_lk_mvue/pplsize)
  LB95.Rel.Bias.Nhat_MBM_lk_mvue <- quantile(Nhat_MBM_lk_mvue/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_MBM_lk_mvue <- quantile(Nhat_MBM_lk_mvue/pplsize, 0.975) %>% as.numeric()
  

  # Combine output
  df <- data.frame(# MBM_lk_mle estimator
                    MSE.Nhat_MBM_lk_mle, Mean.Nhat_MBM_lk_mle, Med.Nhat_MBM_lk_mle, SD.Nhat_MBM_lk_mle, LB95.Nhat_MBM_lk_mle, UB95.Nhat_MBM_lk_mle, Min.Nhat_MBM_lk_mle, Max.Nhat_MBM_lk_mle,
                    Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle,
                    # MBM_lk_mvue estimator
                    MSE.Nhat_MBM_lk_mvue, Mean.Nhat_MBM_lk_mvue, Med.Nhat_MBM_lk_mvue, SD.Nhat_MBM_lk_mvue, LB95.Nhat_MBM_lk_mvue, UB95.Nhat_MBM_lk_mvue, Min.Nhat_MBM_lk_mvue, Max.Nhat_MBM_lk_mvue,
                    Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue,
                    
                    # # MBM_unlk_mle estimator
                    # MSE.Nhat_MBM_unlk_mle, Mean.Nhat_MBM_unlk_mle, Med.Nhat_MBM_unlk_mle, SD.Nhat_MBM_unlk_mle, LB95.Nhat_MBM_unlk_mle, UB95.Nhat_MBM_unlk_mle, Min.Nhat_MBM_unlk_mle, Max.Nhat_MBM_unlk_mle,
                    # # pvalue.t.2side_MBM_unlk_mle, pvalue.t.gtr_MBM_unlk_mle, pvalue.t.less_MBM_unlk_mle, Mean.Nhat_MBM_unlk_mle_Min.m,
                    # Rel.Bias.Nhat_MBM_unlk_mle, LB95.Rel.Bias.Nhat_MBM_unlk_mle, UB95.Rel.Bias.Nhat_MBM_unlk_mle,
                    # # MBM_unlk_mvue estimator
                    # MSE.Nhat_MBM_unlk_mvue, Mean.Nhat_MBM_unlk_mvue, Med.Nhat_MBM_unlk_mvue, SD.Nhat_MBM_unlk_mvue, LB95.Nhat_MBM_unlk_mvue, UB95.Nhat_MBM_unlk_mvue, Min.Nhat_MBM_unlk_mvue, Max.Nhat_MBM_unlk_mvue,
                    # # pvalue.t.2side_MBM_unlk_mvue, pvalue.t.gtr_MBM_unlk_mvue, pvalue.t.less_MBM_unlk_mvue, Mean.Nhat_MBM_unlk_mvue_Min.m,
                    # Rel.Bias.Nhat_MBM_unlk_mvue, LB95.Rel.Bias.Nhat_MBM_unlk_mvue, UB95.Rel.Bias.Nhat_MBM_unlk_mvue,
                    
                    # counts of overlap
                    True.m, Mean.m, Med.m, SD.m, Min.m, Max.m)
}

# Changing p vector from both low to both high
mbm.1.ls <- function(simsize = nsim, pplsize = targetN, p.X1 = scenarios1$p.X1, pB.X0 = scenarios1$pB.X0, pB.X1 = scenarios1$pB.X1, p2 = scenarios1$p2){
  # Start changing scenarios across the vector of pre-specified p
  # Create list, each element is one sub-scenario of pB and p2
  s1.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- length(p2)
  for(l in 1:length.ls){
    s1.ls[[l]] <-  mbm.1(simsize = simsize, pplsize = pplsize, p.X1 = p.X1[l], pB.X0 = pB.X0[l], pB.X1 = pB.X1[l], p2 = p2[l])
    
  }
  # Collapse the list
  s1.comb <- do.call(rbind,s1.ls)
  # Add scenarios
  scenario <- data.frame(p.X1, pB.X0, pB.X1, p2)
  # Final output
  final.s1.comb <- as.data.frame(cbind(scenario,s1.comb))
  return(final.s1.comb)
}


################################################################
######## Scenario case 2 (Sample 2 is conditional on X) #########
################################################################

crm.1.2 <- function(simsize = nsim, pplsize = targetN, p.X1 = scenarios2$p.X1[1], 
                    pB.X0 = scenarios2$pB.X0[1], pB.X1 = scenarios2$pB.X1[1], 
                    p2.X0 = scenarios2$p2.X0[1], p2.X1 = scenarios2$p2.X1[1]){
  # Create datasets
  benchmark.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  sample1.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  sample2.dt <- matrix(0, nrow = simsize, ncol = pplsize)
  
  set.seed(1234)
  # Label Covariate X
  covariate.vec <- rbinom(pplsize, 1, p.X1)
  covariate.dt <- t(replicate(simsize, covariate.vec))
  n.X0 <- sum(covariate.vec == 0)
  n.X1 <- pplsize - n.X0
  
  # Simulate for each element
  for (i in 1:simsize){
    set.seed(1234+i)
    # Benchmark is a random sample
    benchmark.dt[i, which(covariate.vec == 0)] <- rbinom(n.X0, 1, pB.X0)
    benchmark.dt[i, which(covariate.vec == 1)] <- rbinom(n.X1, 1, pB.X1)
    # All benchmark people got sampled into the sample 1 for CRM and MBM
    sample1.dt[i,which(benchmark.dt[i,] == 1)] <- 1
    # Draw individuals for source 2 for CRM to generate proportion
    sample2.dt[i, which(covariate.vec == 0)] <- rbinom(n.X0, 1, p2.X0)
    sample2.dt[i, which(covariate.vec == 1)] <- rbinom(n.X1, 1, p2.X1)
  }
  # Got two data sample, sample1.dt and sample2.dt
  # Marginal counts for sources
  n1 <- apply(sample1.dt,1,sum)
  n2 <- apply(sample2.dt,1,sum)
  # Overlapping individuals, equivalent to labeling the overlapping
  sample1.2.dt<-sample1.dt*sample2.dt
  ####################
  # --- Estimate --- #
  ####################
  # /// CRM /// #
  # Counts of overlapping
  m <- apply(sample1.2.dt,1,sum)
  # LP estimator 
  # Est Total N for CRM
  Nhat_LP <- n1*n2/m
  # Chapman Estimator for CRM
  # Est Total N for CRM
  Nhat_Chap <- (n1+1)*(n2+1)/(m+1) - 1
  # Stratified Estimate using LP and Chapman
  Nhat_LP_Strata <- c()
  Nhat_Chap_Strata <- c()
  for(i in 1:simsize){
    sample.1 <- sample1.dt[i,]
    sample.2 <- sample2.dt[i,]
    sample.m <- sample1.2.dt[i,]
    n1.X1 <- sum(sample.1[which(covariate.vec == 1)])
    n1.X0 <- sum(sample.1[which(covariate.vec == 0)])
    n2.X1 <- sum(sample.2[which(covariate.vec == 1)])
    n2.X0 <-  sum(sample.2[which(covariate.vec == 0)])
    m.X1 <- sum(sample.m[which(covariate.vec == 1)])
    m.X0 <- sum(sample.m[which(covariate.vec == 0)])
    # Nhat.X1_LP <-  n1.X1*n2.X1/m.X1
    # Nhat.X0_LP <-  n1.X0*n2.X0/m.X0
    # Nhat_LP.X <- Nhat.X1_LP + Nhat.X0_LP
    # Nhat_LP_Strata <- c(Nhat_LP_Strata, Nhat_LP.X )
    Nhat.X1_Chap <-  (n1.X1+1)*(n2.X1+1)/(m.X1+1) - 1
    Nhat.X0_Chap <-  (n1.X0+1)*(n2.X0+1)/(m.X0+1) - 1
    Nhat_Chap.X <- Nhat.X1_Chap + Nhat.X0_Chap
    Nhat_Chap_Strata <- c(Nhat_Chap_Strata, Nhat_Chap.X)
  }
  # Summary
  True.m <- pplsize * (1-p.X1) * pB.X1 * p2.X0 + pplsize * p.X1 * pB.X1 * p2.X1 
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
  
  # Stratified
  MSE.Nhat_Chap_Strata <- mean((Nhat_Chap_Strata-pplsize)^2)
  Mean.Nhat_Chap_Strata <- mean(Nhat_Chap_Strata)
  Med.Nhat_Chap_Strata <- median(Nhat_Chap_Strata)
  SD.Nhat_Chap_Strata <- sd(Nhat_Chap_Strata)
  LB95.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata, 0.025) %>% as.numeric()
  UB95.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata, 0.975) %>% as.numeric()
  Min.Nhat_Chap_Strata <- min(Nhat_Chap_Strata)
  Max.Nhat_Chap_Strata <- max(Nhat_Chap_Strata)
  Mean.Nhat_Chap_Strata_Min.m <- mean(Nhat_Chap_Strata[which(m == Min.m)])
  Rel.Bias.Nhat_Chap_Strata <- mean(Nhat_Chap_Strata/pplsize)
  LB95.Rel.Bias.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata/pplsize, 0.025) %>% as.numeric()
  UB95.Rel.Bias.Nhat_Chap_Strata <- quantile(Nhat_Chap_Strata/pplsize, 0.975) %>% as.numeric()
  
  # Combine output
  df <- data.frame(# LP estimator
    MSE.Nhat_LP, Mean.Nhat_LP, Med.Nhat_LP, SD.Nhat_LP, LB95.Nhat_LP, UB95.Nhat_LP, Min.Nhat_LP, Max.Nhat_LP,
    Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP,
    # Chapman estimator
    MSE.Nhat_Chap, Mean.Nhat_Chap, Med.Nhat_Chap, SD.Nhat_Chap, LB95.Nhat_Chap, UB95.Nhat_Chap, Min.Nhat_Chap, Max.Nhat_Chap,
    Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap,
    # Chapman estimator with stratification
    MSE.Nhat_Chap_Strata, Mean.Nhat_Chap_Strata, Med.Nhat_Chap_Strata, SD.Nhat_Chap_Strata, LB95.Nhat_Chap_Strata, UB95.Nhat_Chap_Strata, Min.Nhat_Chap_Strata, Max.Nhat_Chap_Strata,
    Rel.Bias.Nhat_Chap_Strata, LB95.Rel.Bias.Nhat_Chap_Strata, UB95.Rel.Bias.Nhat_Chap_Strata,
    # counts of overlap
    True.m, Mean.m, Med.m, SD.m, Min.m, Max.m)
}

# Changing p vector from both low to both high
crm.1.2.ls <- function(simsize = nsim, pplsize = targetN, p.X1 = scenarios2$p.X1, 
                       pB.X0 = scenarios2$pB.X0, pB.X1 = scenarios2$pB.X1, 
                       p2.X0 = scenarios2$p2.X0, p2.X1 = scenarios2$p2.X1){
  # Start changing scenarios across the vector of pre-specified p
  # Create list, each element is one sub-scenario of pB and p2
  s1.ls <- list()
  # Put in list (Dim = length(p.vex) with each element has dim = simsize * pplsize)
  length.ls <- length(pB.X0)
  for(l in 1:length.ls){
    s1.ls[[l]] <-  crm.1.2(simsize = simsize, pplsize = pplsize, p.X1 = p.X1[l], 
                         pB.X0 = pB.X0[l], pB.X1 = pB.X1[l], 
                         p2.X0 = p2.X0[l], p2.X1 = p2.X1[l])
    
  }
  # Collapse the list
  s1.comb <- do.call(rbind,s1.ls)
  # Add scenarios
  scenario <- data.frame(p.X1, pB.X0, pB.X1, p2.X0, p2.X1)
  # Final output
  final.s1.comb <- as.data.frame(cbind(scenario,s1.comb))
  return(final.s1.comb)
}




#######################################################
# Function to plot CRM vs MBM on different estimators #
#######################################################

simu3.plot.CRM.MBM.diff.Est <- function(dt.crm = simu3.crm.1.Nhat, dt.mbm = NA, p1_pB = pB.X1.vec[4], p.X1 = 0.3, title = "A"){
  if(!is.na(dt.mbm)) { # If with MBM
    CRM.LP <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
    CRM.LP$Methods <- rep("LP", nrow(CRM.LP))
    colnames(CRM.LP) <- c("p1.X1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    CRM.Chap <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
    CRM.Chap$Methods <- rep("Chapman", nrow(CRM.Chap))
    colnames(CRM.Chap) <- c("p1.X1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    CRM.Chap_Strata <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2, Rel.Bias.Nhat_Chap_Strata, LB95.Rel.Bias.Nhat_Chap_Strata, UB95.Rel.Bias.Nhat_Chap_Strata))
    CRM.Chap_Strata$Methods <- rep("Chapman (Stratified)", nrow(CRM.Chap_Strata))
    colnames(CRM.Chap_Strata) <- c("p1.X1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    
    MBM.lk_mle <- subset(dt.mbm, pB.X1 == p1_pB, select = c(pB.X1, p2, Rel.Bias.Nhat_MBM_lk_mle, LB95.Rel.Bias.Nhat_MBM_lk_mle, UB95.Rel.Bias.Nhat_MBM_lk_mle))
    MBM.lk_mle$Methods <- rep("MBM_lk_mle", nrow(MBM.lk_mle))
    colnames(MBM.lk_mle) <- c("p1.X1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    MBM.lk_mvue <- subset(dt.mbm, pB.X1 == p1_pB, select = c(pB.X1, p2, Rel.Bias.Nhat_MBM_lk_mvue, LB95.Rel.Bias.Nhat_MBM_lk_mvue, UB95.Rel.Bias.Nhat_MBM_lk_mvue))
    MBM.lk_mvue$Methods <- rep("MBM_lk_mvue", nrow(MBM.lk_mvue))
    colnames(MBM.lk_mvue) <- c("p1.X1", "p2", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    
    CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, CRM.Chap_Strata, MBM.lk_mle, MBM.lk_mvue))  
    
    p <- ggplot(CRM.MBM.given.p1, aes(x = p1_pB, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
      geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
      geom_line(aes(y = Rel.Bias.Nhat), size = 0.8) +
      geom_point(aes(y = Rel.Bias.Nhat)) +
      geom_hline(aes(yintercept=1), color = "black",
                 linetype="dashed", size=0.5) +
      scale_y_continuous(limits = c(0,8.5), breaks = c(0,2,4,6,8)) +
      ylab(expression(hat(N)/N)) + 
      xlab("Sampling Probability for sample 2") +
      ggtitle(title) +
      theme_classic() +
      theme(plot.title = element_text(size = 14),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1))
  }
  
  if(is.na(dt.mbm) & length(p1_pB)==1){ # If without MBM
    CRM.LP <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2.X1, p2.X0, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
    CRM.LP$Methods <- rep("LP", nrow(CRM.LP))
    colnames(CRM.LP) <- c("p1.X1", "p2.X1", "p2.X0", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    CRM.Chap <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2.X1, p2.X0, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
    CRM.Chap$Methods <- rep("Chapman", nrow(CRM.Chap))
    colnames(CRM.Chap) <- c("p1.X1", "p2.X1", "p2.X0", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    CRM.Chap_Strata <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2.X1, p2.X0, Rel.Bias.Nhat_Chap_Strata, LB95.Rel.Bias.Nhat_Chap_Strata, UB95.Rel.Bias.Nhat_Chap_Strata))
    CRM.Chap_Strata$Methods <- rep("Chapman (Stratified)", nrow(CRM.Chap_Strata))
    colnames(CRM.Chap_Strata) <- c("p1.X1", "p2.X1", "p2.X0", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    
    CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, CRM.Chap_Strata))
    # Compute an overall (weighted average) sampling probability of sample 2. Show this on x-axis
    CRM.MBM.given.p1$ave.p2 <- apply(cbind(CRM.MBM.given.p1$p2.X1, CRM.MBM.given.p1$p2.X0), 1, 
                                     FUN = function(x){return(x[1]*p.X1 + x[2]*(1-p.X1))})

    p <- ggplot(CRM.MBM.given.p1, aes(x = ave.p2, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
      geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
      geom_line(aes(y = Rel.Bias.Nhat), size = 0.8) +
      geom_point(aes(y = Rel.Bias.Nhat)) +
      geom_hline(aes(yintercept=1), color = "black",
                 linetype="dashed", size=0.5) +
      scale_y_continuous(limits = c(0,8.5), breaks = c(0,2,4,6,8)) +
      ylab(expression(hat(N)/N)) + 
      xlab("Sampling Probability of Sample 2") +
      ggtitle(title) +
      theme_classic() +
      theme(plot.title = element_text(size = 14),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.position = "bottom") +
      guides(color = guide_legend(nrow = 1))
    
    return(p)
  }
  
  if(is.na(dt.mbm) & length(p1_pB)>1){
    CRM.LP <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2.X1, Rel.Bias.Nhat_LP, LB95.Rel.Bias.Nhat_LP, UB95.Rel.Bias.Nhat_LP))
    CRM.LP$Methods <- rep("LP", nrow(CRM.LP))
    colnames(CRM.LP) <- c("p1.X1", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    CRM.Chap <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2.X1, Rel.Bias.Nhat_Chap, LB95.Rel.Bias.Nhat_Chap, UB95.Rel.Bias.Nhat_Chap))
    CRM.Chap$Methods <- rep("Chapman", nrow(CRM.Chap))
    colnames(CRM.Chap) <- c("p1.X1", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    CRM.Chap_Strata <- subset(dt.crm, pB.X1 == p1_pB, select = c(pB.X1, p2.X1, Rel.Bias.Nhat_Chap_Strata, LB95.Rel.Bias.Nhat_Chap_Strata, UB95.Rel.Bias.Nhat_Chap_Strata))
    CRM.Chap_Strata$Methods <- rep("Chapman (Stratified)", nrow(CRM.Chap_Strata))
    colnames(CRM.Chap_Strata) <- c("p1.X1", "p2.X1", "Rel.Bias.Nhat", "LB95.Rel.Bias.Nhat", "UB95.Rel.Bias.Nhat","Methods")
    
    CRM.MBM.given.p1 <- as.data.frame(rbind(CRM.LP, CRM.Chap, CRM.Chap_Strata))
    
    p <- ggplot(CRM.MBM.given.p1, aes(x = p1.X1, y = Rel.Bias.Nhat, group = Methods, color = Methods)) + 
      geom_ribbon(aes(ymin = LB95.Rel.Bias.Nhat, ymax = UB95.Rel.Bias.Nhat), linetype=2, alpha=0.1) +
      geom_line(aes(y = Rel.Bias.Nhat), size = 0.8) +
      geom_point(aes(y = Rel.Bias.Nhat)) +
      geom_hline(aes(yintercept=1), color = "black",
                 linetype="dashed", size=0.5) +
      ylab("Nhat/N") + 
      xlab("Sampling Probability for sample 1(X=1)") +
      ggtitle(paste("Heterogeneous sample 2: p2(X=1)=",p2.X1, sep = " ")) +
      theme_classic() +
      theme(plot.title = element_text(size = 14),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16),
            legend.position = "bottom") +
      guides(color = guide_legend(nrow = 2))
      
  }
  return(p)   
}


