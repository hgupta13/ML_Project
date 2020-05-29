################################################
# Date: 2020.05.29
# Course: ECON 293 - ML
# Project: Final Project
################################################

  rm(list = ls())
  options(scipen=999)
  
############################
# LOAD PACKAGES
############################
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(rmarkdown)
  library(kableExtra)
  library(utils)
  library(readxl)
  library(randomForest)
  library(glmnet)
  library(clubSandwich)
  library(sandwich)
  
############################
# LOAD DATA
############################
  df <- read_csv('sample.csv') %>%
    select(-X1) %>%
    mutate(gender = gender - 1) %>%
    filter(!is.na(NGO_donations),
           !is.na(student_wealth))
  itt_out <- list()
  
############################
# FUNCTIONS
############################

  
############################
# REPLICATE TABLE 3 OF ROMERO ET AL.
############################
  reg_fun <- paste0('IRT_score ~ ', 
                    # Group fixed effect
                    'factor(groupid) + ',
                    # Treatment
                    'W + ',
                    # Student fixed characteristics
                    paste('age', 'gender', 'student_wealth', 'factor(grade)', sep = ' + '),
                    ' + ',
                    # School fixed characteristics
                    paste('school_facilities', 'time_to_bank', 'enrollment_2015', 'rural', 'NGO_donations', sep = ' + '))

  
############################
# Calculate ITT OLS, cluster standard errors at school level
############################
  lm_out <- lm(reg_fun, data = df) 
  itt_unweighted <- coef_test(lm_out, vcov = "CR1", 
                              cluster = df$schoolid, test = "naive-t")['W', ] 
  itt_unweighted <- tibble(type = 'unweighted',
                    estimate = itt_unweighted$beta,
                    se = itt_unweighted$SE,
                    ci_low = itt_unweighted$beta - 1.96 * itt_unweighted$SE,
                    ci_high = itt_unweighted$beta + 1.96 * itt_unweighted$SE)
  
############################
# 2 - IPW Estimator
############################
  W <- as.matrix(df %>% select(W))
  
    # Manually create indicator variables for X$grade
    grs <- matrix(nrow = nrow(W), ncol = length(unique(df$grade)), 0L)
    grades <- unique(df$grade)
    for(i in 1:length(grades)) {
      grs[ , i] <- as.numeric(grades[[i]] == df$grade)
    }
    X <- cbind(grs,
               as.matrix(df %>% 
                           select(school_facilities, time_to_bank, enrollment_2015, 
                                  rural, NGO_donations, gender, age, student_wealth)))
  # Estimate propensity weights
  p <- glm(W ~ X, family = "binomial") %>%
      predict(type = 'response')
  
  # Append propensity weights to df
  df_ipw <- df %>%
    # Calculate propensity weights
    mutate(weight = (W / p) + ((1 - W) / (1 - p)))
  
  # Estimate IPW regression
  lm_out <- lm(reg_fun, data = df_ipw, weights = weight)
  itt_ipw <- coef_test(lm_out, vcov = "CR1", 
                              cluster = df_ipw$schoolid, test = "naive-t")['W', ]
  itt_ipw <- tibble(type = 'IPW',
                    estimate = itt_ipw$beta,
                    se = itt_ipw$SE,
                    ci_low = itt_ipw$beta - 1.96 * itt_ipw$SE,
                    ci_high = itt_ipw$beta + 1.96 * itt_ipw$SE)
  
############################
# 3 - Double robust estimator
############################
  aipw_fun <- paste0('IRT_score ~ ', 
                    # Group fixed effect
                    'factor(groupid) + ',
                    # School fixed characteristics
                    paste('school_facilities', 'time_to_bank', 'enrollment_2015', 'rural', 'NGO_donations', sep = ' + '),
                    # Interaction
                    ' + W * (',
                    # Student fixed characteristics
                    paste('age', 'gender', 'student_wealth', 'factor(grade)', sep = ' + '),
                    ')')
  lm_out <- lm(aipw_fun, data = df)   
    
    # Predict - all treated
    df_treatall <- df %>% 
      mutate(W = 1)
    y_treatall <- predict(lm_out, df_treatall)
    
    # Predict - all control
    df_treatnone <- df %>% 
      mutate(W = 0)
    y_treatnone <- predict(lm_out, df_treatnone)
    
    # Predict actual
    actual_pred = predict(lm_out, df)
    
    # Calculate AIPW estimate
    G <- y_treatall - y_treatnone +
      ((df$W - p) * (df$IRT_score - actual_pred)) / (p * (1 - p))
    tau.hat <- mean(G)
    se.hat <- sqrt(var(G) / (length(G) - 1))
    
    # Format output
    itt_aipw <- tibble(type = 'aipw',
                         estimate = tau.hat,
                         se = se.hat,
                         ci_low = tau.hat - 1.96 * se.hat,
                         ci_high = tau.hat + 1.96 * se.hat)

    
############################
# FORMAT OUTPUT - ITT
############################
  # 1 - Propensity score overlap
  plot_tib <- tibble(treatment = factor(W), score = p)
  ggplot(plot_tib) +
      geom_histogram(aes(x = score, y = stat(density), fill = treatment),
                     alpha = 0.3, position = 'identity') +
    labs(x = 'Propensity Score',
         y = 'Density',
         title = 'Logit Propensity Scores')
  
  # 2 - Table of ITT Estimates
  itt_table <- bind_rows(itt_unweighted,
                         itt_ipw,
                         itt_aipw)
  kable(itt_table, digits = 3)


  