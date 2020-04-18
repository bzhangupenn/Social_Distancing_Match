library(sensitivitymult)

#############Primary Outcome############
data_matched<-read.csv("/Users/heng6/Desktop/Social_Distancing_Match-master/matched_pairs_with_outcome.csv")
data_matched<-data_matched[, -c(14, 15)]

#We calculate the residuals for the covariance adjusted outcomes
p_value<-NULL
beta_0=0
outcome_adjusted<-log(data_matched$average_obs_illness_outcome/100)-beta_0*data_matched$avg_social_distance_score
data_matched_temp<-cbind(data_matched, outcome_adjusted)
fit <- lm(outcome_adjusted ~ female + above_65 + median_house_income + unemployment
          + smoking + flu_vaccination + some_college
          + social_association + traffic + RUCC + poverty
          + population.density + avg_obs_ill_0316_0322
          + avg_exp_ill_0316_0322 + avg_exp_ill_0330_0405 + bl_hisp, data=data_matched_temp)
residuals_cov_adjust=fit$residuals
#The treatment indicator
t_ind<-rep(0, nrow(data_matched))
for (i in 1:(nrow(data_matched)/2)){
  t_ind[2*i]=1
}
dev<-senm(residuals_cov_adjust, t_ind, data_matched$id, alternative = "greater")$deviate
p_value=1-pchisq(dev^2, df = 1)


find_p_value<-function(beta_0, data_input, alter_direction){
  outcome_adjusted<-log(data_input$average_obs_illness_outcome/100)-beta_0*data_input$avg_social_distance_score
  data_input_temp<-cbind(data_input, outcome_adjusted)
  fit <- lm(outcome_adjusted ~ female + above_65 + median_house_income + unemployment
            + smoking + flu_vaccination + some_college
            + social_association + traffic + RUCC + poverty
            + population.density + avg_obs_ill_0316_0322
            + avg_exp_ill_0316_0322 + avg_exp_ill_0330_0405 + bl_hisp, data=data_input_temp)
  residuals_cov_adjust=fit$residuals
  t_ind<-rep(0, nrow(data_input_temp))
  for (i in 1:(nrow(data_input_temp)/2)){
    t_ind[2*i]=1
  }
  dev<-senm(residuals_cov_adjust, t_ind, data_input_temp$id, alternative = alter_direction)$deviate
  data_input_temp<-NULL
  p_value=1-pchisq(dev^2, df = 1)
  return(p_value)
}

beta_0_vec<-seq(from=-0.1, to=0.1, by=0.0001)
indicate_beta_two_sided<-rep(0, length(beta_0_vec))
for (i in 1:length(indicate_beta_two_sided)){
  p_value_two_sided=find_p_value(beta_0 = beta_0_vec[i], data_input = data_matched, alter_direction = "greater")
  indicate_beta_two_sided[i]=as.numeric(p_value_two_sided>0.05)
  if(i%%200==0){
    print(i)
  }
}

#CI
exp(max(beta_0_vec[indicate_beta_two_sided==1]))
exp(min(beta_0_vec[indicate_beta_two_sided==1]))

##############Secondary Outcome##########
data_matched<-read.csv("/Users/heng6/Desktop/Social_Distancing_Match-master/matched_pairs_with_outcome_406.csv")
ata_matched<-data_matched[, -c(14, 15)]
data_matched$average_obs_illness_outcome[data_matched$average_obs_illness_outcome==0]=0.0001

#We calculate the residuals for the covariance adjusted outcomes
p_value<-NULL
beta_0=0
outcome_adjusted<-log(data_matched$average_obs_illness_outcome/100)-beta_0*data_matched$avg_social_distance_score
data_matched_temp<-cbind(data_matched, outcome_adjusted)
fit <- lm(outcome_adjusted ~ female + above_65 + median_house_income + unemployment
          + smoking + flu_vaccination + some_college
          + social_association + traffic + RUCC + poverty
          + population.density + avg_obs_ill_0316_0322
          + avg_exp_ill_0316_0322 + avg_exp_ill_0330_0405 + bl_hisp, data=data_matched_temp)
residuals_cov_adjust=fit$residuals
#The treatment indicator
t_ind<-rep(0, nrow(data_matched))
for (i in 1:(nrow(data_matched)/2)){
  t_ind[2*i]=1
}
dev<-senm(residuals_cov_adjust, t_ind, data_matched$id, alternative = "greater")$deviate
p_value=1-pchisq(dev^2, df = 1)


find_p_value<-function(beta_0, data_input, alter_direction){
  outcome_adjusted<-log(data_input$average_obs_illness_outcome/100)-beta_0*data_input$avg_social_distance_score
  data_input_temp<-cbind(data_input, outcome_adjusted)
  fit <- lm(outcome_adjusted ~ female + above_65 + median_house_income + unemployment
            + smoking + flu_vaccination + some_college
            + social_association + traffic + RUCC + poverty
            + population.density + avg_obs_ill_0316_0322
            + avg_exp_ill_0316_0322 + avg_exp_ill_0330_0405 + bl_hisp, data=data_input_temp)
  residuals_cov_adjust=fit$residuals
  t_ind<-rep(0, nrow(data_input_temp))
  for (i in 1:(nrow(data_input_temp)/2)){
    t_ind[2*i]=1
  }
  dev<-senm(residuals_cov_adjust, t_ind, data_input_temp$id, alternative = alter_direction)$deviate
  data_input_temp<-NULL
  p_value=1-pchisq(dev^2, df = 1)
  return(p_value)
}

beta_0_vec<-seq(from=-0.1, to=0.1, by=0.0001)
indicate_beta_two_sided<-rep(0, length(beta_0_vec))
for (i in 1:length(indicate_beta_two_sided)){
  p_value_two_sided=find_p_value(beta_0 = beta_0_vec[i], data_input = data_matched, alter_direction = "greater")
  indicate_beta_two_sided[i]=as.numeric(p_value_two_sided>0.05)
  if(i%%200==0){
    print(i)
  }
}

#CI
exp(max(beta_0_vec[indicate_beta_two_sided==1]))
exp(min(beta_0_vec[indicate_beta_two_sided==1]))
