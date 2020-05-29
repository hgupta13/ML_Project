################################################
# Date: 2020.05.29
# Course: ECON 293 - ML
# Project: Final Project
# Description: Estimate ATE in each matched pair
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
    select(-X1)

  # Estimate ATE within each matched pair of schools
  out <- df %>%
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
  saveRDS(out, file = 'Intermediate/2020.05.29 Matched Treatment Effect.RDS')
