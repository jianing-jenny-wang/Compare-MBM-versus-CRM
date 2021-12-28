# Investigation on problems from simulation
library(ggplot2)
### For CRM LP estimator and Chapman estimator ###
### Q1: Why small sample create bias on Mean.Nhat? ###
# Dig into comparison between LP and Chapman estimator for CRM (dual-system), on different sampling probabilities
# Ref: 
# 1. T1.5S1N30 Mauricio Sadinle 2008
# 2. T1.5S1N46, T1.5S1N46.1 both are MLE math derivation
# 3. Temp S1N85 Candace M. Lynn 2009 Master Thesis: Statement " The Lincoln-Petersen method, although simple, is subject to bias if m2 is small or zero, 
# and it can overestimate the population."
summary(Nhat_LP)
summary(Nhat_Chap)
bad.index <- which(Nhat_LP > quantile(Nhat_LP, 0.75))
m[bad.index]
comb.bad.index <- data.frame(m[bad.index], Nhat_LP[bad.index], n1[bad.index], n2[bad.index])

summary(m)
small.m <- which(m < quantile(m, 0.25))
comb.small.m <- data.frame(m[small.m], Nhat_LP[small.m], n1[small.m], n2[small.m])

# For each pair of p1 and p2, plot the density of Nhat
# Basic density
Nhat_LP.dt <- data.frame(Nhat_LP)
Min.Nhat_LP <- min(Nhat_LP)
Max.Nhat_LP <- max(Nhat_LP)
ggplot(Nhat_LP.dt, aes(x=Nhat_LP)) +  
  geom_density()  +  
  geom_vline(aes(xintercept=mean(Nhat_LP), color="Mean"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=median(Nhat_LP), color="Median"),
             linetype="dashed", size=1, show.legend = T) +
 geom_vline(aes(xintercept=pplsize, color="True"),
             linetype="dashed", size=1, show.legend = T) +
  scale_color_manual("Lines",
                     values=c("Mean"="blue","Median"="red","True"="black")) +
  ggtitle(paste("Plot for", "p1 =", p1, ",", "p2 =", p2, sep = " ")) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) 

Nhat_Chap.dt <- data.frame(Nhat_Chap)
Min.Nhat_Chap <- min(Nhat_Chap)
Max.Nhat_Chap <- max(Nhat_Chap)
ggplot(Nhat_Chap.dt, aes(x=Nhat_Chap)) +  
  geom_density()  +  
  geom_vline(aes(xintercept=mean(Nhat_Chap), color="Mean"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=median(Nhat_Chap), color="Median"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=pplsize, color="True"),
             linetype="dashed", size=1, show.legend = T) +
  scale_color_manual("Lines",
                     values=c("Mean"="blue","Median"="red","True"="black")) +
  ggtitle(paste("Plot for", "p1 =", p1, ",", "p2 =", p2, sep = " ")) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) 

# Finding: When m is small, such as m = c(4,5,6) in current run with seed (1234), the N are much over-estimated, over 10000 (twice than the true N 5000)
# Then: 
# Why this tends to overestimate when m is small?
# Answer: m is very small and p is very small
# What happen to LP estimator?
# Answer: check overleaf step 3 file for MLE of hypergeometric and proof for bias via jensen's inequality, ref are T1.5S1N46, T1.5S1N46.1 both are MLE math derivation

# Try to do a one-sample t-test (Two-sided if not sure Nhat is larger/smaller than trueN, one-side if sure Nhat is larger)
# Two sided
a <- t.test(Nhat_LP, mu = pplsize)
# One sided
t.test(Nhat_LP, mu = pplsize, alternative = "greater")
# Two sided
t.test(Nhat_Chap, mu = pplsize)
# One sided
t.test(Nhat_Chap, mu = pplsize, alternative = "greater")






### For MBM NB estimator or direct calculation estimator ###
# For each pair of p1 and p2, plot the density of Nhat
# Basic density
Nhat_MBM_nb.dt <- data.frame(Nhat_MBM_nb)
Min.Nhat_MBM_nb <- min(Nhat_MBM_nb)
Max.Nhat_MBM_nb <- max(Nhat_MBM_nb)
ggplot(Nhat_MBM_nb.dt, aes(x=Nhat_MBM_nb)) +  
  geom_density()  +  
  geom_vline(aes(xintercept=mean(Nhat_MBM_nb), color="Mean"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=median(Nhat_MBM_nb), color="Median"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=pplsize, color="True"),
             linetype="dashed", size=1, show.legend = T) +
  scale_color_manual("Lines",
                     values=c("Mean"="blue","Median"="red","True"="black")) +
  ggtitle(paste("Plot for", "pB =", pB, "(fixed Benchmark across simulation),", "p2 =", p2, sep = " ")) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) 

Nhat_MBM_direct.dt <- data.frame(Nhat_MBM_direct)
Min.Nhat_MBM_direct <- min(Nhat_MBM_direct)
Max.Nhat_MBM_direct <- max(Nhat_MBM_direct)
ggplot(Nhat_MBM_direct.dt, aes(x=Nhat_MBM_direct)) +  
  geom_density()  +  
  geom_vline(aes(xintercept=mean(Nhat_MBM_direct), color="Mean"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=median(Nhat_MBM_direct), color="Median"),
             linetype="dashed", size=1, show.legend = T) +
  geom_vline(aes(xintercept=pplsize, color="True"),
             linetype="dashed", size=1, show.legend = T) +
  scale_color_manual("Lines",
                     values=c("Mean"="blue","Median"="red","True"="black")) +
  ggtitle(paste("Plot for", "pB =", pB, "(fixed Benchmark across simulation),", "p2 =", p2, sep = " ")) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) 


# Note: if we want to compare two simulation results, use two-sample t test


