# I use "fastDummies" package to create the dummies, and "Matrix" to handle the large sparse matrix with dummies
rm(list = ls())

#DATA
df <- import('/Users/andresfelipe/Dropbox/StanfordEconPhD/Year2/ML and Causal Inference/Final project/ML_Project/sample.csv') %>%
  mutate(gender = gender - 1) %>%
  select(-V1) %>%
  dplyr::filter(!is.na(NGO_donations),
                !is.na(student_wealth))

#CREATE MATRIX FOR CF
set.seed(321)
X <- select(df, -W, -IRT_score,-studentid) %>%
  fastDummies::dummy_cols(select_columns = c('grade','schoolid'),remove_first_dummy = TRUE) %>%
  select(-grade,-schoolid, groupid) %>% as.matrix() %>%
  as('dgCMatrix')
W <- as.matrix(df$W)
Y <- as.matrix(df$IRT_score)

#EST CF
cf <- causal_forest(X,Y,W,num.trees = 2000)
cf_cate <- predict(cf, estimate.variance=TRUE)
tauhat_cf <- cf_cate$predictions
hist(tauhat_cf)

#CHECK ATE OF CF
ate_cf_aipw = average_treatment_effect(cf, target.sample = "overlap") 
  #Larger than real ATE (0.29), problem of not including group
  #If re run including group with dummies, ATE decreases to something like 0.24, not great but lower.
    #What could be going on???




###Group analysis
group_ate <- tibble(group=df$groupid,cate=tauhat_cf) %>%
  group_by(group) %>%
  summarise(meancate=mean(cate))
hist(group_ate$meancate)

reg_gr <- paste0('IRT_score ~ ',
                 # Student fixed characteristics
                 paste('age', 'gender', 'student_wealth', 'factor(grade)', sep = ' + '),
                 ' + ',
                 # School fixed characteristics
                 paste('school_facilities', 'time_to_bank', 'enrollment_2015', 
                       'rural', 'NGO_donations', sep = ' + '))
reggroup <- lm(reg_gr,data = df) %>%
  predict(type = 'response')

df <- mutate(df, Yhat=reggroup)
df1 <- df %>% 
  dplyr::filter(W==1) %>%
  group_by(groupid) %>%
  summarise(Y1=mean(IRT_score))
df0 <- df %>% dplyr::filter(W==0) %>%
  group_by(groupid) %>%
  summarise(Y0=mean(IRT_score))

df_g1 <- dplyr::full_join(df1,df0,by='groupid') %>%
  mutate(ATE=Y1-Y0)

df1 <- df %>% 
  dplyr::filter(W==1) %>%
  group_by(groupid) %>%
  summarise(Y1=mean(Yhat))
df0 <- df %>% dplyr::filter(W==0) %>%
  group_by(groupid) %>%
  summarise(Y0=mean(Yhat))

df_g2 <- dplyr::full_join(df1,df0,by='groupid') %>%
  mutate(ATE=Y1-Y0)

hist(df_g1$ATE)
hist(df_g2$ATE)


