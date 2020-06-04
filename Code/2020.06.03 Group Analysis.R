##Load packages

library(readr)
library(dplyr)
library(ggplot2) 
library(gridExtra)
library(dgof)
library(tidyr)
library(stargazer)
library(xtable)
library(stats)


##Get data on avergae of treatment effects by group
group_data = read.csv("group_ates.csv")

# Clean data
group_data$X = NULL


#Plot histogram of the true ATE distribution and the distribution of average of treatment effects across groups for different methods
#The goal is to see we get the same distribution of treatment effects across groups

#S-Learner
d1 = gather(group_data,key="Method", value = "Estimate", true_ate, ate_s)
p1 = ggplot(d1, aes(x = Estimate, color = Method)) +
  geom_histogram(fill="white", position="dodge", binwidth = 0.1)
  
#T-Learner
d2 = gather(group_data,key="Method", value = "Estimate", true_ate, ate_t)
p2 = ggplot(d2, aes(x = Estimate, color = Method)) +
  geom_histogram(fill="white", position="dodge", binwidth = 0.1)


#X-Learner
d3 = gather(group_data,key="Method", value = "Estimate", true_ate, ate_x)
p3 = ggplot(d3, aes(x = Estimate, color = Method)) +
  geom_histogram(fill="white", position="dodge", binwidth = 0.1)

# Causal Forest (R-Learner)
d4 = gather(group_data,key="Method", value = "Estimate", true_ate, ate_cf)
p4 = ggplot(d4, aes(x = Estimate, color = Method)) +
  geom_histogram(fill="white", position="dodge", binwidth = 0.1)


#Plot the Histograms
grid.arrange(p1,p2,p3,p4)

#Analytically Test whether 
#the true ATE distribution and the distribution of average of treatment effects across groups for different methods
t1 = ks.test(group_data$true_ate, group_data$ate_s)
t2  = ks.test(group_data$true_ate, group_data$ate_t)
t3 = ks.test(group_data$true_ate, group_data$ate_x)
t4 = ks.test(group_data$true_ate, group_data$ate_cf)

#Summarize Results
km_df <- tibble(Method = c("S-Learner" , "T-Learner", "X-Learner", "Causal Forest"), Conclusion = rep("Reject", 4),pvalue = c(t1$p.value,t2$p.value,t3$p.value,t4$p.value))

####################################################################################################################################
#Clustering Analysis
#Examine whether the treatment effects cluster in a that resembeles the distribution across providers

#Get clusters for true ATE, X-Learner, and T-Learner
#We use 8(=number of providers) clusters
k1 = kmeans(group_data$true_ate,8)
k2 = kmeans(group_data$ate_t,8)
k3 = kmeans(group_data$ate_x,8)

#Plot the centers of the clusters for various methods
k_clusters = data.frame(True_ATE = sort(k1$centers), T_Learner = sort(k2$centers),
                          X_Learner = sort(k3$centers))
