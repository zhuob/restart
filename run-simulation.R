# environment setup 


rm(list = ls())
setwd("C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code")

source('C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code/simulation-all-in-one.R', encoding = 'UTF-8')
source('C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code/correlated-orr-os-params.R')



# Efficacy parameters ------------- For stage 1
nstage1 <- 60;
n_ia   <- 30;
theta2 <- 0.25;
theta1 <- 0.30;
delta1 <- 0.95;
delta2 <- 0.95;
gamma1 <- 0.05;
gamma2 <- 0.90
orr0   <- 0.25
orr1   <- 0.45

a0 <- 1; 
b0 <- 1;

theta_null <- 0.25
theta_alt  <- 0.45

orr_fu_time <- 2


compute_alpha_power(nsbj = nstage1, a0 = a0, b0 = b0, theta2 = theta2, 
                    delta2 = delta2, theta_null = orr0, theta_alt = orr1)



# Efficacy parameters --------------   for stage 2
mos0 <- 8.5; 
mos1 <- 17
rho1 <- 0.5
nsbj <- 180 # if expansion, then up to 180 subjects in total
nstage2 <- nsbj - nstage1

enrl_rate <- 5
drop_timecut <- 12; 
drop_prob <- 0.02; 
sig_level <- 0.05


lambda0 <- find_lambda(orr0, rho1, mos0)
rho2 <- find_rho2(orr1, lambda0, rho1, mos1)
rho2 <- ifelse(rho2 > 1, 1, ifelse(rho2 < 0, 0, rho2))

mos <- log(2)/c(lambda0, rho1*lambda0, rho2*lambda0, rho1*rho2*lambda0)
names(mos) <- c("SOC/Non-responder", "SOC/Responder", "TRT/Non-Responder", "TRT/Responder")
mos

library(rpact)
d1 <- getSampleSizeSurvival(alpha = sig_level, sided = 2, beta = 0.1, lambda1 = log(2)/mos0, 
                            lambda2 = log(2)/mos1, accrualTime = c(0, nstage2/enrl_rate), accrualIntensity = 5)

nevent <- ceiling(d1$eventsFixed)

compute_alpha_power(nsbj = nstage1, a0 = a0, b0 = b0, theta2 = theta2, 
                    delta2 = delta2, theta_null = 0.25, theta_alt = 0.45)

#############
pp_distr <- pred_prob_all(n1 = n_ia, n = nstage1, a0 = a0, b0 = b0, theta1 = theta1, theta2 = theta2)


batch1 <- run_one_simu(nsbj = nsbj, nstage1 = nstage1, n_ia = n_ia, 
                       orr0 = orr0, orr1 = orr1, rho1 = rho1, rho2 = rho2, 
                       lambda0 = lambda0, enrl_timecut = 0, enrl_rate = enrl_rate, 
                       drop_timecut = drop_timecut, drop_prob = drop_prob, orr_fu_time = orr_fu_time,
                       pp_distr = pp_distr, a0 = a0, b0 = b0,
                       theta1 = theta1, theta2 = theta2, delta1 = delta1, delta2 = delta2, 
                       gamma1 = gamma1, gamma2 = gamma2, nevent = nevent, sig_level = sig_level,
                       seed = 123454)

all_result <- NULL 
for(i in 1:10){
  
  
  batchi <- run_one_simu(nsbj = nsbj, nstage1 = nstage1, n_ia = n_ia, 
                         orr0 = orr0, orr1 = orr1, rho1 = rho1, rho2 = rho2, 
                         lambda0 = lambda0, enrl_timecut = 0, enrl_rate = enrl_rate, 
                         drop_timecut = drop_timecut, drop_prob = drop_prob, orr_fu_time = orr_fu_time,
                         pp_distr = pp_distr, a0 = a0, b0 = b0,
                         theta1 = theta1, theta2 = theta2, delta1 = delta1, delta2 = delta2, 
                         gamma1 = gamma1, gamma2 = gamma2, nevent = nevent, sig_level = sig_level,
                         seed = i)
  
  all_result <- bind_rows(all_result, batchi %>% mutate(seed = i))
}

# summarize simulations

sim_summary <- function(all_result){
  
  ia_prob <- all_result %>% group_by(ia_decision) %>% summarize(ia_prob = n()/nrow(all_result))
  
  ia_summary <- all_result %>% summarize(
    ia_time = mean(ia_time), ia_n = mean(ia_n_enrl), theta = mean(ia_n_resp/ia_n_eval)
  )
  
  fa_summary <- all_result %>% mutate(
    fa_success = case_when(ia_decision == "Expand to Two-arm" & fa_p <  sig_level ~ TRUE, 
                           ia_decision == "Expand to Two-arm" & fa_p >= sig_level ~ FALSE, 
                           TRUE ~ fa_success)
  ) %>% mutate(
    success_type = case_when(ia_decision == "Continue Single Arm" & fa_success ~ "Single Arm Success", 
                             ia_decision == "Continue Single Arm" & !fa_success ~ "Single Arm Failure", 
                             ia_decision == "Expand to Two-arm" & !fa_success ~ "Two arm Expansion Failure",
                             ia_decision == "Expand to Two-arm" & fa_success ~ "Two arm Expansion Success",
                             ia_decision == "Futile"  & fa_success ~ "Futile then Success", 
                             ia_decision == "Futile"  & !fa_success ~ "Futile then Failure")
  )
  
  
  
}

names(all_result) <- c("IA - number of Evaluable Sbjs", "IA - Number of Enrolled Sbjs", "IA - Number of Responds", 
                       "IA time", "IA decision", "IA Pred.Prob-Futile", "IA Pred.Prob.Success", "FA time", 
                       "FA N", "FA Number of responds", "FA Pos. Prob.", "FA Success?", "Simulation seed", "FA P value (if Expand)", 
                       "# Events", "HR", "Lower 95", "Upper 95", "mOS SOC", "mOS TRT")

write_csv(all_result, "sim1.csv", na = "")



# varying parameters to run multiple scenarios ----------------------------------------------------
rm(list = ls())
library(tidyverse)
r_files <- c('C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code/simulation-all-in-one.R',
             'C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code/correlated-orr-os-params.R')

for(kk in r_files){
  source(kk)
}

########### fixed parameters
nsbj         <- 180
nstage1      <- 60
nstage2      <- nsbj - nstage1
drop_timecut <- 12
drop_prob    <- 0.02
enrl_timecut <- 0
orr_fu_time  <- 1
a0           <- 1
b0           <- 1
theta2 <- 0.25;
theta1 <- 0.30;
delta1 <- 0.95;
delta2 <- 0.95;
gamma1 <- 0.05;
gamma2 <- 0.90
sig_level <- 0.05

mos0 <- 8.5; 
mos1 <- 17


########### parameters to vary
n_ia <- c(20, 30, 50)
rho1 <- c(0.1, 0.5, 1) # hazard ratio between non-responder and responder
enrl_rate <- c(3, 12)

# get the sample size for Stage 2
d1 <- rpact::getSampleSizeSurvival(alpha = sig_level, sided = 2, beta = 0.1, lambda1 = log(2)/mos0, 
                                   lambda2 = log(2)/mos1, accrualTime = c(0, nstage2/enrl_rate[1]),
                                   accrualIntensity = enrl_rate[1])
nevent <- ceiling(d1$eventsFixed)

p0 <- bind_rows(tidyr::expand_grid(orr0 = 0.25, orr1 = 0.25, mos0 = 8.5,  mos1 = 8.5,   rho1 = rho1),
                tidyr::expand_grid(orr0 = 0.25, orr1 = 0.35, mos0 = 8.5,  mos1 = 12.75, rho1 = rho1),
                tidyr::expand_grid(orr0 = 0.25, orr1 = 0.45, mos0 = 8.5,  mos1 = 17,    rho1 = rho1))
                


param1 <- purrr::pmap_df(.l = as.list(p0), .f = find_lambda_rho2) %>% 
  mutate(`Relationship between ORR and OS` = rep(c("Strong", "Moderate", "Weak"), 3))

params <- tidyr::expand_grid(bind_cols(p0, param1 %>% select(-rho1)), n_ia, enrl_rate) %>% 
    mutate(scenario = 1:n()) %>% select(scenario, orr1, rho1, n_ia, enrl_rate, everything())

write_csv(params, "sim-result/params-20-30-50.csv")

combo_sim <- NULL 

for(k in 1:nrow(params)){
  
  p1 <- params %>% slice(k)
  cat("\r", k)
  pp_distr <- pred_prob_all(n1 = p1$n_ia, n = nstage1, a0 = a0, b0 = b0, theta1 = theta1, theta2 = theta2)
  
  # run RESTART
  tmp <- beverage::run_parallel_sim(ncores = 3, nsim = 1e3, core_fun = run_one_simu, seed = k,
                                    combine_method = bind_rows, parallel = FALSE, package_used = "beverage",
                                    file_to_source = r_files,
                                    nsbj = nsbj, nstage1 = nstage1, n_ia = p1$n_ia,
                                    orr0 = p1$orr0, orr1 = p1$orr1, rho1 = p1$rho1, rho2 = p1$rho2,
                                    lambda0 = p1$lambda0, enrl_timecut = enrl_timecut, enrl_rate = p1$enrl_rate,
                                    drop_timecut = drop_timecut, drop_prob = drop_prob, orr_fu_time = orr_fu_time,
                                    pp_distr = pp_distr, a0 = a0, b0 = b0,
                                    theta1 = theta1, theta2 = theta2, delta1 = delta1, delta2 = delta2,
                                    gamma1 = gamma1, gamma2 = gamma2, nevent = nevent, sig_level = sig_level)
 
  
  # run comparator
  # tmp <- beverage::run_parallel_sim(ncores = 3, nsim = 1e3,  core_fun = run_one_simu_comparator, seed = k,
  #                                   combine_method = bind_rows, parallel = FALSE, package_used = "beverage",
  #                                   file_to_source = r_files, 
  #                                   nsbj = nsbj, nstage1 = nstage1, 
  #                                   orr0 = p1$orr0, orr1 = p1$orr1, rho1 = p1$rho1, rho2 = p1$rho2, 
  #                                   lambda0 = p1$lambda0, enrl_timecut = enrl_timecut, enrl_rate = p1$enrl_rate, 
  #                                   drop_timecut = drop_timecut, drop_prob = drop_prob, orr_fu_time = orr_fu_time,
  #                                   a0 = a0, b0 = b0,
  #                                   theta2 = theta2, delta3 = 0.10, delta4 = 0.90, 
  #                                   nevent = nevent, sig_level = sig_level)
  # 
  # s1 <- run_one_simu_comparator(nsbj = 180, nstage1 = 60, orr0 = 0.25, orr1 = 0.35,
  #                               rho1 = 0.1, rho2 = 0.813, lambda0 = 0.118, 
  #                               enrl_timecut = 0, enrl_rate = 3, 
  #                               drop_timecut = 12, drop_prob = 0.02, orr_fu_time = 1,
  #                               a0 = 1, b0 = 1, theta2 = 0.25, delta3 = 0.80, delta4 = 0.98, 
  #                               nevent = 88, sig_level = 0.05,
  #                               seed = 1234)
  # 
  
  
  tmp <- tmp %>% mutate(scenario = k) %>% select(scenario, everything())
  
  combo_sim <- bind_rows(combo_sim, tmp)
  
  
  write_csv(combo_sim, "sim-result/sim-v5-start-ia-20-30-50.csv", na = "")
  
  
}

# combo_summary <- summarize_sim(combo_sim = combo_sim)
combo_summary <- summarize_sim_comparator(combo_sim = combo_sim)
s0 <- left_join(params, combo_summary, by = "scenario")

plot_sim_result(s0)

plot_sim_result(s0, var_of_interest = "fa_n", "Total Sample Size")
plot_sim_result(s0, var_of_interest = "fa_time", "Study Duration (months)")



