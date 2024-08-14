rm(list = ls())
setwd("C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code")
# source("sim-functions.R")
# source("result-summary.R")
source('C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code/simulation-all-in-one.R', encoding = 'UTF-8')
library(dplyr)
library(tidyr)
library(ggplot2)

nsbj     <- 60           # number of subjects
p0       <- 0.26         # reference control ORR
n_ia     <- 30           # number subjects at interim for expansion cohort
a0_info  <- 7.8          # to be discussed
b0_info  <- 22.2         # to be discussed
a0_noninfo <- 1          # to be discussed
b0_noninfo <- 1          # to be discussed
a1         <- 0.5        # prior for treatment
b1         <- 0.5        # prior for treatment
delta      <- 0          # 
delta_ia   <- 0.05
theta_suc  <- 0.8        # success cutoff PP(Pos(ORR > pref_suc) > theta_suc)
theta_fut  <- 0.975        # similar to success, futility cutoff
ppfut_cut  <- 0.05       # PP(Pos(ORR > pref_fut) > theta_fut) < ppfut_cut
pref_final <- 0.2       # final success criteria lower95CI > pref_fianl is success

###############################  Tune this parameter #####################
pref_fut <- 0.20
pref_suc <- 0.25

#########################################################################
ppsuc_cut0 <- c(0.7, 0.8, 0.9)
x_step <- seq(20, 50, by = 10)
p1 <- c(0.1, 0.2, 0.3, 0.4)  # treatment arm ORR scenarios

params <- expand.grid(nsbj = nsbj, theta_suc = theta_suc, 
                      theta_fut = theta_fut, 
                      pref_fut = pref_fut, pref_suc = pref_suc, 
                      a0 = 0.5, b0 = 0.5, ppsuc_cut = ppsuc_cut0, 
                      ppfut_cut = ppfut_cut,
                      eff_stop = FALSE, fut_stop = TRUE)





# power calibration
x1 <- 0:nsbj
pos1 <- pos_prob(x1, n1 = nsbj, a0 = 0.5, b0 = 0.5, p0 = pref_suc)
xdens1 <- dbinom(x1, size = nsbj, prob = 0.40)
xm <- tibble(pos1 = pos1, x1 = x1, xdens = xdens1, 
             actual_power = rev(cumsum(rev(xdens1)))) %>% mutate(success = pos1 > 0.95)

power0 <- xm %>% filter(success) %>% pull(xdens) %>% sum()

pos2 <- pos_prob(x1, n1 = nsbj, a0 = 0.5, b0 = 0.5, p0 = 0.3)
xdens2 <- dbinom(x1, size = nsbj, prob = 0.26)
xn <- tibble(pos2 = pos2, x1 = x1, xdens = xdens2, 
             actual_alpha = cumsum(xdens2)) %>% mutate(futile = pos2 > 0.95)
alpha0 <- xn %>% filter(futile) %>% pull(xdens) %>% sum()
c(alpha0, power0)


# evaluate cutoff values


result <- run_sim(params = params, x_step = c(x_step, 60), p1 = p1)

fig <- plot_status_prob(result %>% select(-futility) %>%
                          rename(`No change` = success, Futility = com_prob_fut, 
                                 `Expand` = continue, `conflicted` = conflicted), 
                        view_scenario = p1[1:4], var_start = 6, theta_name = "ppsuc_cut",  
                        status_order = c("conflicted", "Futility", "Expand", "No change"),
                        view_theta_u = ppsuc_cut0, show_text = TRUE)

fig


ggsave("/userdata/cfda/dni/simulation/onc/amg509/Simulation/Bin/state-prob.png", width = 9, height = 7)



##  flip flop simulation

tmp <- run_one_trial_ff(nsbj = nsbj, nsbj_step = x_step, 
                        a1 = a1, b1 = b1, pref_fut = pref_fut, pref_suc = pref_suc, 
                        theta_fut = theta_fut, theta_suc = theta_suc, p1 = 0.4, 
                        pref_final = pref_final, ppfut_cut = ppfut_cut, ppsuc_cut = 0.8)

nsim <- 1e4
ff_result <- NULL
for(i in 1:length(p1)){
  
  temp2 <- run_parallel_sim(nsim = nsim, ncores = 25, seed = 123 + i, parallel = TRUE,
                            core_fun = run_one_trial_ff, combine_method = bind_rows,
                            nsbj = nsbj, nsbj_step = x_step, a1 = a1, b1 = b1, 
                            pref_fut = pref_fut, pref_suc = pref_suc, 
                            theta_fut = theta_fut, theta_suc = theta_suc, p1 = p1[i], 
                            pref_final = pref_final, ppfut_cut = ppfut_cut, ppsuc_cut = 0.8) 
  
  ff_result <- bind_rows(ff_result, temp2 %>% mutate(scenario = p1[i]))
  
}

ff_data <- flip_flop(result = ff_result, nsim = nsim) %>% 
  mutate(prob_rel = prob / prob_interim) 


plot_flip_flop(ff_data %>% select(-prob) %>% rename(prob = prob_rel))
ggsave("/userdata/cfda/dni/simulation/onc/amg509/Simulation/Bin/flip-flop.png", width = 9, height = 7)
write_csv(ff_data, "/userdata/cfda/dni/simulation/onc/amg509/Simulation/Bin/flip-flop-table.csv")
plot_flip_flop2(ff_data)
ggsave("/userdata/cfda/dni/simulation/onc/amg509/Simulation/Bin/flip-flop-all-categories.png", width = 9, height = 7)


# an example

step_result <- purrr::map_df(.x = c(20, 30, 40, 50, 60), .f = step_decision,  
                             nsbj = 60, theta_fut = 0.9, theta_suc = 0.8, 
                             pref_fut = 0.1, pref_suc = 0.3, ppfut_cut = 0.05,
                             ppsuc_cut = 0.8, a0 = 0.5, b0 = 0.5,
                             eff_stop = FALSE, fut_stop = TRUE)

###################  RUN INDEPENDENT ONE TRIAL #################################
tmp <- run_one_trial_indep(nsbj = 60, n_ia = 30, a0_info = 10.14, p1 = 0.10, p0 = 0.26,
                           b0_info = 28.86, a0_noninfo = 1, b0_noninfo = 1, pref_final = 0.2,
                           pref_fut = 0.10, pref_suc = 0.30, theta_fut = 0.90, theta_suc = 0.80,
                           ppfut_cut = 0.05, ppsuc_cut = 0.8, a1 = 0.5, b1 = 0.5, delta_ia = 0.05, 
                           delta = 0, seed = 89987)


rbind_result <- NULL 
nsim <- 1e4

for (k in 1:length(p1)){
  temp2 <- run_parallel_sim(nsim = nsim, ncores = 25, seed = 123 + k, parallel = TRUE, 
                            package_used = c("beverage", "tidyverse"),
                            core_fun = run_one_trial_indep, combine_method = bind_rows,
                            # temp2 <- run_many_trial(nsim = nsim, ncores = 2, seed = 123 + k,
                            nsbj = nsbj, n_ia = n_ia, a0_info = a0_info, 
                            b0_info = b0_info, a0_noninfo = a0_noninfo, 
                            b0_noninfo = b0_noninfo, p0 = p0, p1 = p1[k], 
                            pref_final = pref_final, 
                            pref_fut = pref_fut, pref_suc = pref_suc, 
                            theta_fut = theta_fut, theta_suc = theta_suc,
                            ppfut_cut = ppfut_cut, ppsuc_cut = 0.8,
                            a1 = a1, b1 = b1, delta_ia = delta_ia, delta = delta) 
  
  
  rbind_result <- bind_rows(rbind_result, temp2 %>% mutate(scenario = p1[k], sim = 1:nsim))
  
}

powertb <- power_option1_5(rbind_result, nsim = nsim, poscut = 0.9)

plot_option1_5(powertb)
write_csv(powertb, "/userdata/cfda/dni/simulation/onc/amg509/Simulation/Bin/option1-5-power.csv")
ggsave("/userdata/cfda/dni/simulation/onc/amg509/Simulation/Bin/option1-5-power.png", width = 9, height = 7)


#####################################################################################################################
# EXPLORING RELATIONSHIP BETWEEN ORR AND OS
#####################################################################################################################
df <- data_gen_orr_os(nsbj = 180, nstage1 = 60, orr0 = 0.17, orr1 = 0.47, 
                      rho1 = 0.07, rho2 = 1, lambda0 = 0.1048986, 
                      enrl_timecut = 0, enrl_rate = 5, drop_timecut = 12, 
                      drop_prob = 0.02, seed = 1234)
table(df$arm)

df %>% group_by(arm) %>% summarize(orr = mean(resp))
df %>% group_by(resp) %>% summarize(os = median(os))

n_ia <- 30;       # number of subjects at interim
pp_distr <- pred_prob_all(n1 = n_ia, n = 60, a0 = 1, b0 = 1, theta1 = 0.3, theta2 = 0.2)
orr_fu_time <- 2; # each patient will be followed for 2 months before ORR assessment

delta1 <- 0.90; delta2 = 0.95

stage1_result <- stage1_ia(df, n_ia = 30, orr_fu_time = 2, pp_distr = pp_distr, delta1 = 0.9,
                delta2 = 0.99, gamma1 = 0.05, gamma2 = 0.9)

fa_result <- run_final_analysis(df, stage1_result, nstage1 = nstage1, a0 = 1, 
                                b0 = 1, theta2 = 0.2, delta2 = 0.95, sig_level = 0.05)


n_ia <- 30;       # number of subjects at interim
pp_distr <- pred_prob_all(n1 = n_ia, n = 60, a0 = 1, b0 = 1, theta1 = 0.3, theta2 = 0.2)

batch1 <- run_one_simu(nsbj = 180, nstage1 = 60, n_ia = 30, 
                       orr0 = 0.17, orr1 = 0.47, rho1 = 0.07, rho2 = 1, 
                       lambda0 = 0.1048986, enrl_timecut = 0, enrl_rate = 5, 
                       drop_timecut = 12, drop_prob = 0.02, orr_fu_time = 2,
                       pp_distr = pp_distr, a0 = 1, b0 = 1,
                       theta1 = 0.3, theta2 = 0.2, delta1 = 0.9, delta2 = 0.99, 
                       gamma1 = 0.05, gamma2 = 0.9, sig_level = 0.05,
                       seed = 1234)

## several questions:
# 1. given results from code in line 64-73, the probability of reaching
# inconclusive result at N = 60 is very slim, are we adjusting the boundaries?
# 2. We probably need boundary values for futility/efficacy/expansion to 2 arms
# 3. if expanded to 2-arm, how to evaluate OC when combined with stage 1

compute_alpha_power(nsbj = 60, a0 = 1, b0 = 1,  theta2 = 0.25, 
                    delta2 = 0.95, theta_null = 0.25, theta_alt = 0.45)


# sample size for survival 

library(rpact)
d1 <- getSampleSizeSurvival(alpha = 0.05, sided = 2, beta = 0.1, lambda1 = log(2)/8.5, 
                            lambda2 = log(2)/17, accrualTime = c(0, 24), accrualIntensity = 5)

nevent <- ceiling(d1$eventsFixed)

rho1 <- 0.01; mos1 <- 8.5; mos2 <- 17; orr1 <- 0.25; orr2 <- 0.45
lambda0 <- find_lambda(orr1, rho1, mos1)
rho2 <- find_rho2(orr2, lambda0, rho1, mos2)
rho2 <- ifelse(rho2 > 1, 1, ifelse(rho2 < 0, 0, rho2))

mos <- log(2)/c(lambda0, rho1*lambda0, rho2*lambda0, rho1*rho2*lambda0)
names(mos) <- c("SOC/Non-responder", "SOC/Responder", "TRT/Non-Responder", "TRT/Responder")
mos


t1 <- get_quantity(rho1, rho2, lambda0, orr1, orr2) %>% mutate(hr = hz_trt/hz_soc)

t2 <- t1 %>%  pivot_longer(cols = 2:8, names_to = "type", values_to = "value") %>% 
  mutate(category = case_when(
    stringr::str_detect(type, "hz") ~ "hazard rate", 
    stringr::str_detect(type, "resp|noresp") ~ "response status", 
    stringr::str_detect(type, "soc|trt") ~ "arm", 
    stringr::str_detect(type, "hr") ~ "hazard ratio", 
  ))

t3 <- t2 %>% split(f = t2$category)

p1 <- ggplot(data = t3$arm, aes(x = time, y = value)) + geom_line(aes(color = type)) + 
  facet_wrap(facets = "category", scale = 'free')

p2 <- p1 %+% t3$`hazard rate`
p3 <- p1 %+% t3$`hazard ratio`
p4 <- p1 %+% t3$`response status`
gridExtra::grid.arrange(p1, p2, p3, p4)




