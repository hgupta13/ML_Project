################################################
# Date: 6/1/2020
# Course: ECON 293
# Project: Final Project
# Description: HTE Estimation
################################################

  rm(list = ls())
  options(scipen=999)
  
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(rmarkdown)
  library(kableExtra)
  library(utils)
  library(readxl)
  library(grf)
  library(devtools)
  library(causalTree)
  library(gridExtra)
  library(causalTree)
  library(estimatr)
  
############################
# LOAD DATA
############################
  df <- read_csv('sample.csv') %>%
    select(-X1) %>%
    mutate(gender = gender - 1) %>%
    filter(!is.na(NGO_donations),
           !is.na(student_wealth))
  
############################
# SPLIT, PREPARE DATA
############################
  # Split into test, train samples. Stratify split by group and treatment
  # to ensure that the test and train datasets can calculate CATEs within each
  # randomization group
  df_train <- df %>%
    group_by(W, groupid) %>%
    sample_frac(0.6) %>%
    ungroup()
  df_test <- df[(df$studentid %in% df_train$studentid) == FALSE, ]
  
  # Remove student id from data frame, recast groupid as factor variable
  df_train <- df_train %>%
    select(-c(studentid))
  df_test <- df_test %>%
    select(-c(studentid))
  
  # Split data frame into outcome, features, treatment
  W <- df_train$W
  Y <- df_train$IRT_score
  X <- df_train %>% 
    select(-c(W, IRT_score))
  
  # Repeat for test set
  W_test <- df_test$W
  Y_test <- df_test$IRT_score
  X_test <- df_test %>% 
    select(-c(W, IRT_score))

############################
# S LEARNER
############################      
  # 1 - Calculate S-Learner (see slide 22 of Lecture 4)
  s_learn <- regression_forest(cbind(W, X), Y)
  pred_s_0 <- predict(s_learn, cbind(0, X))$predictions
  pred_s_1 <- predict(s_learn, cbind(1, X))$predictions
  pred_s_oob <- predict(s_learn)$predictions
  pred_s_0[W == 0] <- pred_s_oob[W == 0] 
  pred_s_1[W == 1] <- pred_s_oob[W == 1]
  pred_s <- pred_s_1 - pred_s_0 
  
  # Predict on test data
  tauhat_s <- predict(s_learn, cbind(1, X_test))$predictions - 
    predict(s_learn, cbind(0, X_test))$predictions
  
############################
# T LEARNER
############################   
  # 2 - Calculate T-Learner (see slide 21 of Lecutre 4)
  tf0 <- regression_forest(X[W==0,], Y[W==0])
  tf1 <- regression_forest(X[W==1,], Y[W==1])
  tf.preds.0 <- predict(tf0, X)$predictions
  tf.preds.1 <- predict(tf1, X)$predictions
  tf.preds.0[W==0] <- predict(tf0)$predictions #OOB
  tf.preds.1[W==1] <- predict(tf1)$predictions #OOB
  pred_t <- tf.preds.1 - tf.preds.0
  
  # Predict on test data
  tauhat_t <- predict(tf1, X_test)$predictions - 
    predict(tf0, X_test)$predictions
  
############################
# X LEARNER
############################ 
  # A - Predict Y_i from X_i when W_i == 0 (use T-learner forest 0),
  # Learn \tau_1 by predicting delta from X_i when W_i == 1
  yhat0 = predict(tf0, X[W==1,])$predictions
  xf1 = regression_forest(X[W==1,], Y[W==1]-yhat0)
  xf.preds.1 = predict(xf1, X)$predictions
  xf.preds.1[W==1] = predict(xf1)$predictions
  
  # B - Swap: Predict Y_i from X_i when W_i == 1 (use T-learner forest 1),
  # Learn \tau_0
  yhat1 = predict(tf1, X[W==0,])$predictions
  xf0 = regression_forest(X[W==0,], yhat1-Y[W==0])
  xf.preds.0 = predict(xf0, X)$predictions
  xf.preds.0[W==0] = predict(xf0)$predictions
  
  # C - Estimate the propensity score - regression forest
  # of W on X
  propf = regression_forest(X, W)  # , tune.parameters = TRUE)
  ehat = predict(propf)$predictions
  
  # D - Estimate \hat{\tau}(x)
  pred_x = (1 - ehat) * xf.preds.1 + ehat * xf.preds.0
  
  # E - Predict
  ehat.test <- predict(propf, X_test)$predictions
  xf.preds.1.test <- predict(xf1, X_test)$predictions
  xf.preds.0.test <- predict(xf0, X_test)$predictions
  tauhat_xl_test <- (1 - ehat.test) * xf.preds.1.test + ehat.test * xf.preds.0.test
  
############################
# CAUSAL FOREST
############################ 
  cf <- causal_forest(X, Y, W, num.trees = dim(X)[1])
  pred_cf <- predict(cf)$predictions
  tauhat_cf <- predict(cf, newdata = X_test)$predictions
  
############################
# WITHIN-GROUP ATE
############################
  # Estimate ATE within each matched pair of schools
  group_ate <- df %>%
    # Collapse to pair-treatment group level
    group_by(groupid, W) %>%
    mutate(avg_score = mean(IRT_score),
           sq_error = (avg_score - IRT_score) ^ 2) %>%
    summarize(avg_score = mean(avg_score),
              count = n(),
              sum_error = sum(sq_error)) %>%
    ungroup() %>%
    # Collapse to pair level to estimate pair treatment effect, error
    mutate(avg_score = if_else(W == 0, -1 * avg_score, avg_score),
           st_err = sum_error / (count - 1)) %>%
    group_by(groupid) %>%
    summarize(treatment_effect = sum(avg_score),
              st_err = sum(st_err),
              count = sum(count)) %>%
    ungroup() %>%
    mutate(st_err = sqrt(st_err))
  
  # Map into test dataset
  tauhat_group_ate <- df_test %>% 
    inner_join(group_ate) %>%
    select(treatment_effect) %>%
    c() %>%
    unlist()

############################
# R Loss
############################
  # Estimate Y_hat, W_hat on test set
  Y.forest.test = regression_forest(X = as.matrix(X_test), Y = Y_test)
  Y.hat.test = predict(Y.forest.test)$predictions
  W.forest.test = regression_forest(X = as.matrix(X_test), Y = W_test)
  W.hat.test = predict(W.forest.test)$predictions
  
  # Calculate sq error loss
  sq_error_loss <- tibble(
    ate = (Y_test - Y.hat.test - (W_test - W.hat.test) * tauhat_group_ate) ^ 2,
    s_learner = (Y_test - Y.hat.test - (W_test - W.hat.test) * tauhat_s) ^ 2,
    t_learner = (Y_test - Y.hat.test - (W_test - W.hat.test) * tauhat_t) ^ 2,
    causal_forest = (Y_test - Y.hat.test - (W_test - W.hat.test) * tauhat_cf) ^ 2,
    x_learner = (Y_test - Y.hat.test - (W_test - W.hat.test) * tauhat_xl_test) ^ 2
    ) %>%
    summarize_all(sum)
  