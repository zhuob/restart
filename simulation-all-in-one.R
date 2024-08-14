#‘ given required parameters, generate correlated orr and OS
#'
#' @param nsbj number of subjects to enroll
#' @param nstage1 sample size in stage 1 where only ORR will be assessed
#' @param orr0 response rate at control arm 
#' @param orr1 response rate at treatment arm
#' @param lambda0 hazard rate for non-responder in the control arm 
#' @param enrl_timecut start time of enrollment
#' @param rho1 hazard ratio for non-responder vs responder 
#' @param rho2 additional hazard ratio for treatment effect independent of ORR
#' @param drop_timecut 
#' @param drop_prob 
#' @param seed 
#' @param enrl_rate enrollment speed
#'
#' @return
#' @export
#'
#' @examples
data_gen_orr_os <- function(nsbj, nstage1, orr0, orr1, rho1, rho2, lambda0, enrl_timecut = 0,
                            enrl_rate = 5, drop_timecut = 12, drop_prob = 0.05, seed = 12345){
  
  # rho1 <- 1; rho2 <- 0.5; 
  # # mOS SOC 8.5 months
  # # mOS treatment 17 months
  # # ORR: SOC 0.17 vs treatment 0.47
  # orr0 <- 0.17; orr1 <- 0.47
  # mos0 <- 8.5;  mos1 <- 17
  # nsbj <- 1e3
  set.seed(seed)
  library(dplyr)
  timein <- beverage::rand_inout(nsbj = nsbj, enrl_timecut = enrl_timecut, 
                                 enrl_rate = enrl_rate, drop_timecut = drop_timecut, 
                                 drop_prob = drop_prob)
  
  # treatment assignment for stage 1 and stage 2
  arms_stage1 <- rep("TRT", nstage1)
  arms_stage2 <- beverage::rand_arm(nsbj = nsbj - nstage1, ratio = c(1, 1), arm_name = c("SOC", "TRT"))
  arms <- c(arms_stage1, arms_stage2)
  
  # simulate responder or non-responder
  resp <- rep(NA, nsbj)
  soc_ind <- arms == "SOC"
  resp[soc_ind] <- rbinom(n = sum(soc_ind), 1, orr0)
  resp[!soc_ind] <- rbinom(n = sum(!soc_ind), 1, orr1)
  
  df0 <- timein %>% mutate(arm = arms, id = 1:nsbj, resp = resp) %>% select(id, arm, everything())
  
  ####### matrix 
  #       Non-Responder      Responder
  # SOC    lambda              rho1*lambda
  # COmbo  rho2*lambda         rho1*rho2*lambda
  #
  # rho1: hazard ratio for the responder and non-responder
  # rho2: additional hazard ratio treatment effect independent of ORR 
  # lambda: hazard 
  # we can specify a range of relationship between ORR and OS by varying rho1 and rho2
  
  # based on response status, simulate OS 
  lambda_soc_0  <- lambda0;      # hazard rate for non-responder in SOC
  lambda_soc_1  <- lambda0*rho1; # hazard rate for responder in SOC
  lambda_trt_0  <- lambda0*rho2  # hazard rate for non-responder in combo
  lambda_trt_1  <- lambda0*rho1*rho2 # hazard rate for responder in combo
  
  df0_soc_0 <- df0 %>% filter(arm == "SOC" & resp == 0) %>% mutate(os = rexp(n = n(), rate = lambda_soc_0))
  df0_soc_1 <- df0 %>% filter(arm == "SOC" & resp == 1) %>% mutate(os = rexp(n = n(), rate = lambda_soc_1))
  df0_trt_0 <- df0 %>% filter(arm == "TRT" & resp == 0) %>% mutate(os = rexp(n = n(), rate = lambda_trt_0))
  df0_trt_1 <- df0 %>% filter(arm == "TRT" & resp == 1) %>% mutate(os = rexp(n = n(), rate = lambda_trt_1))
  
  df0 <- bind_rows(df0_soc_0, df0_soc_1, df0_trt_0, df0_trt_1) %>% arrange(id)
  
  return(df0)
}


## calculate predictive probabilities based on 
#         https://trialdesign.org/one-page-shell.html#BEMPR

#' posterior probability of final efficacy given current data
#'
#' @param x1 number of responder in the treatment arm
#' @param n1 number of subjects in the treatment arm
#' @param a0 prior for treatment
#' @param b0 prior for treatment
#' @param p0 a fixed number the treatment arm will compare to
#'
#' @return
#' @export
#'
pos_prob <- function(x1, n1, a0, b0, p0){
  
  pbeta(p0, a0 + x1, b0 + n1 - x1, lower.tail = FALSE)
}


## Y follows a beta-binomial distribution
#' density function for beta-binomial distribution
#'
#' @param y  future number of responders
#' @param n2 future sample size
#' @param x1 past number of responder
#' @param n1 past sample size
#' @param a0 prior 
#' @param b0 prior
#'
#' @return
#' @export
#'
#' @examples

beta_binom_y <- function(y, n2, x1, n1, a0, b0){
  
  f1 <- function(theta){
    dbinom(y, n2, theta)*dbeta(theta, a0 + x1, b0 + n1 - x1)
  }
  
  integrate(f1, 0, 1)$value
}


#' based on the current data, calculate the distribution for future data
#'
#' @param x1 current number of responder
#' @param n1 current sample size evaluated
#' @param n total sample size
#' @param a0,b0 prior for beta distribution
#' @param theta1 reference response rate to be compared, to calculate P(p1>theta1)
#' @param theta2 reference response rate to be compared, to calculate P(p1>theta2)
#' 
#'
#' @return
#' @export
#'
#' @examples
pred_prob <- function(x1, n1, n, a0, b0, theta1, theta2){
  
  y <- 0:(n-n1)
  n2 <- n - n1
  pred_y_prob <- purrr::map_dbl(.x = y, .f = beta_binom_y, n2 = n2, 
                                x1 = x1, n1 = n1, a0 = a0, b0 = b0)
  
  # for futility, calculate Pos(ORR > theta1)
  pos_fut <- purrr::map_dbl(.x = x1 + y, .f = pos_prob, n1 = n, 
                            a0 = a0, b0 = b0, p0 = theta1)
  
  # for success, calculate Pos(ORR > theta2)
  pos_suc <- purrr::map_dbl(.x = x1 + y, .f = pos_prob, n1 = n, 
                            a0 = a0, b0 = b0, p0 = theta2)
  
  result <- tibble::tibble(x1 = x1, n1 = n1, y = y, pred_y = pred_y_prob, 
                           pos_fut = pos_fut, pos_suc = pos_suc)
  return(result)
}


# given sample size at interim, enumerate distribution for all possible outcome
# for future data
pred_prob_all <- function(n1, n, a0, b0, theta1, theta2){
  
  res <- NULL 
  for(i in 0:n1){
    tmp <- pred_prob(x1 = i, n1 = n1, n = n, a0 = a0, b0 = b0, theta1, theta2) 
    res <- bind_rows(res, tmp)
  }
  
  return(res)
}

# x0 <- pred_prob(x1 = 10, n1 = 40, n = 60, a0 = 1, b0 = 1, theta1 = 0.3, theta2 = 0.2)
# x1 <- pred_prob_all(n1 = 40, n = 60, a0 = 1, b0 = 1, theta1 = 0.3, theta2 = 0.2)


#' Calculate predictive probablity based on a given cutoff success criteria
#'
#' @param d_pp density of predicitive probability, resulted from \code{pred_prob}
#' @param delta1 the cutoff value for early futility P(theta > theta1) > delta1
#' @param delta2 the cutoff value for early success P(theta > theta2) > delta2
#' @return
#' @export
#'
#' @examples
pp <- function(d_pp, delta1, delta2){
  
  pp_fut <- d_pp %>% mutate(ind_u = pos_fut > delta1) %>% 
    filter(ind_u == 1) %>% pull(pred_y) %>% sum
  
  pp_suc <- d_pp %>% mutate(ind_u = pos_suc > delta2) %>% 
    filter(ind_u == 1) %>% pull(pred_y) %>% sum
  
  return(tibble(pp_fut = pp_fut, pp_suc = pp_suc))
  
}



# Perform analysis of stage 1


#' Interim analysis for stage 1, evaluating ORR and calculate predictive
#' probability, based on the PP decide whether the trial should go with
#' single-arm, expand to two-arm or stop due to futility
#' decision criteria 
#' if PP(Pos(theta > theta1) > delta1) < gamma1 then futile
#' if PP(Pos(theta > theta2) > delta2) > gamma1 then continue with single arm
#' otherwise expand to two-arm trial
#'
#' @param df simulated data from \code{data_gen_orr_os}
#' @param n_ia,nstage1 number of subjects needed for IA and total N for stage 1
#' @param orr_fu_time follow-up time for patient to be eligible for evaluation
#' @param pp_distr the predictive probability distribution from \code{pred_prob_all}
#' @param delta1,delta2 cutoff value for posterior probability to claim futility
#'   (delta1) or one-arm (delta2)
#' @param gamma1,gamma2 cutoff values for predictive probabilities
#'
#' @return
#' @export
#'
#' @examples
stage1_ia <- function(df, n_ia, nstage1, orr_fu_time, pp_distr, delta1, delta2, gamma1, gamma2){
  
  # decide the time at which interim will be performed
  ia_time <- df$timein[n_ia] + orr_fu_time
  # number of subjects enrolled at IA: if enrollment is too fast, then cap stage 1
  # sample size at nstage1
  n_ia_enrl <- min(nstage1, sum(df$timein <= ia_time))
  
  # numeer of responders at IA
  ia_snapshot <- df %>% slice(1:n_ia)
  nresp_ia <- sum(ia_snapshot$resp == 1)  
  
  # calculate predictive probability based on IA outcome
  dens_pp <- pp_distr %>% filter(x1 == nresp_ia)
  pred_prob <- pp(d_pp = dens_pp, delta1 = delta1, delta2 = delta2)
  # make interim decision 
  ia_decision <- case_when(pred_prob$pp_fut < gamma1 ~ "Futile", 
                           pred_prob$pp_suc > gamma2 ~ "Continue Single Arm", 
                           TRUE ~ "Expand to Two-arm")
  
  # report the results
  res <- tibble::tibble(ia_n_eval = n_ia, ia_n_enrl = n_ia_enrl, 
                        ia_n_resp = nresp_ia, ia_time, ia_decision, 
                        ia_pp_fut = pred_prob$pp_fut, ia_pp_suc = pred_prob$pp_suc)
  
  return(res)
  
  
}


#' Run final analysis based on interim result
#'
#' @param df the simulated data
#' @param stage1_result interim analysis result
#' @param nstage1 number of subjects planned for stage 1
#' @param theta2,delta2 calculate posterior probability, using same criteria for
#'   evaluating interim "continue", i.e., Pos(theta > theta2) > delta2
#' @param a0,b0 posterior probabilites
#' @param sig_level significance level for survival analysis
#'
#' @return
#' @export
#'
#' @examples
run_final_analysis <- function(df, stage1_result, nstage1, a0, b0, theta2, delta2, nevent, sig_level){
  
  # get interim decision
  ia_decision <- stage1_result$ia_decision
  # calculate follow-up time for ORR
  orr_fu_time <- stage1_result$ia_time - df$timein[stage1_result$ia_n_eval]
  
  if(ia_decision == "Continue Single Arm"){
    fa_time <- df$timein[nstage1] + orr_fu_time
    # number of responders at FA
    fa_snapshot <- df %>% slice(1:nstage1)
    nresp_fa <- sum(fa_snapshot$resp == 1)
    # calculate posterior probability based on the ORR outcome
    fa_pos <- pos_prob(x1 = nresp_fa, n1 = nstage1, a0 = a0, b0 = b0, p0 = theta2)
    
    res <- tibble::tibble(fa_time = fa_time, fa_n = nstage1, fa_n_resp = nresp_fa, 
                          fa_pos, fa_success = fa_pos > delta2)
    
  } else if (ia_decision == "Expand to Two-arm"){
    # use only stage 2 data to perform the analysis
    # note: in this case, the enrollment will start immediately after IA, not
    # after completion of Stage 1
    nstage2 <- nrow(df) - nstage1
    # get the updated enrollment time for stage 2
    timein_stage2 <- df %>% slice((stage1_result$ia_n_enrl + 1):(stage1_result$ia_n_enrl + nstage2)) %>% pull(timein)
    fa_sub <- df %>% slice(nstage1+1:nrow(df)) %>% 
      mutate(t_pfs = os, t_os = os, timein = timein_stage2)
    fa_snapshot <- beverage::take_snapshot_event(surv_dat = fa_sub, events = nevent, type = "os")
    fa_time <- fa_snapshot$snapshot_time
    snapshot <- fa_snapshot$snapshot
    res0 <- run_survival(time = snapshot$os, censor = snapshot$os_censor,
                         arm = snapshot$arm, control = "SOC", alpha = sig_level)
    res <- res0 %>% mutate(fa_time = fa_time) %>% 
      select(fa_n = n, fa_time, fa_p = pvalue, event, hr, lower95, upper95, median_surv_1, median_surv_2) 
    
  } else if(ia_decision == "Futile"){
  
    fa_snapshot <- df %>% slice(1:stage1_result$ia_n_enrl)
    # final analysis timing
    fa_time <- max(fa_snapshot$timein) + orr_fu_time
    fa_n_enrl <- stage1_result$ia_n_enrl
    # final analysis 
    nresp_fa <- sum(fa_snapshot$resp == 1)
    # calculate posterior probability based on the ORR outcome
    fa_pos <- pos_prob(x1 = nresp_fa, n1 = fa_n_enrl, a0 = a0, b0 = b0, p0 = theta2)
    
    res <- tibble::tibble(fa_time = fa_time, fa_n = fa_n_enrl, fa_n_resp = nresp_fa,
                          fa_pos, fa_success = fa_pos > delta2)
    
  }
  
  res1 <- bind_cols(stage1_result, res)
  
  return(res1)
}
  
  
#' Run a single trial
#'
#' @param rho1 hazard ratio for the responder vs non-responder
#' @param rho2 additional hazard ratio for treatment vs control
#' @param lambda0 hazard rate for non-responder in control arm 
#' @param enrl_timecut start time of enrollment 
#' @param enrl_rate enrollment speed 
#' @param drop_timecut time span to calculate dropout 
#' @param drop_prob drop out probability 
#' @param orr_fu_time ORR follow up time
#' @param pp_distr predictive probability distribution
#' @param sig_level significance level for stage 2 trial 
#' @param seed simulation seed for reproducibility purpose
#' @param nsbj total number of subjects (stage1 + stage2)
#' @param nstage1 N at stage 1
#' @param n_ia N at interim of Stage 1
#' @param orr0 response rate for control arm
#' @param orr1 response rate for treatment arm
#' @param a0 beta prior 
#' @param b0 beta prior
#' @param theta1,theta2,delta1,delta2,gamma1,gamma2 
#' decision criteria 
#' if PP(Pos(theta > theta1) > delta1) < gamma1 then futile
#' if PP(Pos(theta > theta2) > delta2) > gamma1 then continue with single arm
#' otherwise expand to two-arm trial
#' @param nevent number of events to trigger final analysis
#'
#' @return
#' @export
#'
#' @examples
run_one_simu <- function(nsbj, nstage1, n_ia, 
                         orr0, orr1, rho1, rho2, lambda0, 
                         enrl_timecut, enrl_rate, drop_timecut, drop_prob, orr_fu_time,
                         pp_distr, a0, b0,
                         theta1, theta2, delta1, delta2, gamma1, gamma2, nevent, sig_level,
                         seed = 1234){
  
  # generate data
  df <- data_gen_orr_os(nsbj = nsbj, nstage1 = nstage1, orr0 = orr0, orr1 = orr1, 
                        rho1 = rho1, rho2 = rho2, lambda0 = lambda0, 
                        enrl_timecut = enrl_timecut, enrl_rate = enrl_rate, drop_timecut = drop_timecut, 
                        drop_prob = drop_prob, seed = seed)
  # run stage1 analysis
  stage1_result <- stage1_ia(df, n_ia = n_ia, nstage1 = nstage1, orr_fu_time = orr_fu_time, 
                             pp_distr = pp_distr, delta1 = delta1,
                             delta2 = delta2, gamma1 = gamma1, gamma2 = gamma2)
  
  # run final analysis
  fa_result <- run_final_analysis(df, stage1_result, nstage1 = nstage1, 
                                  a0 = a0, b0 = b0, theta2 = theta2, delta2 = delta2, 
                                  nevent = nevent, sig_level = sig_level)    

  return(fa_result)
  
}   
  
  

#' evaluate type 1 error and power 
#'
#' @param nsbj number of subjects in stage 1
#' @param a0,b0 prior for beta  
#' @param theta2,delta2 parameters such that Pos(theta > theta2) > delta2 
#' @param theta_null scenario under the null
#' @param theta_alt scenario under the alternative
#'
#' @return
#' @export
#'
#' @examples
compute_alpha_power <- function(nsbj, a0, b0, theta2, delta2, theta_null, theta_alt){
  
  # power calibration
  x1 <- 0:nsbj
  pos1 <- pos_prob(x1, n1 = nsbj, a0 = a0, b0 = b0, p0 = theta2)
  xdens1 <- dbinom(x1, size = nsbj, prob = theta_alt)
  xm <- tibble(pos1 = pos1, x1 = x1, xdens = xdens1) %>% mutate(success = pos1 > delta2)
  
  power0 <- xm %>% filter(success) %>% pull(xdens) %>% sum()
  # type 1 error evaluation
  pos2 <- pos_prob(x1, n1 = nsbj, a0 = a0, b0 = b0, p0 = theta2)
  xdens2 <- dbinom(x1, size = nsbj, prob = theta_null)
  xn <- tibble(pos2 = pos2, x1 = x1, xdens = xdens2) %>% mutate(futile = pos2 > delta2)
  alpha0 <- xn %>% filter(futile) %>% pull(xdens) %>% sum()
  c(alpha0, power0)
  
}

#' For comparator analysis, stage 1 decision criteria
#'   Pos(theta > theta2) > delta4   -- success
#'   Pos(theta > theta2) < delta3   -- failure
#'   Otherwise                      -- Expand to two arms  

run_one_simu_comparator <- function(nsbj, nstage1, orr0, orr1, rho1, rho2, lambda0, 
                                    enrl_timecut, enrl_rate, drop_timecut, drop_prob, orr_fu_time,
                                    a0, b0, theta2, delta3, delta4, nevent, sig_level,
                                    seed = 1234){
  
  # generate data
  df <- data_gen_orr_os(nsbj = nsbj, nstage1 = nstage1, orr0 = orr0, orr1 = orr1, 
                        rho1 = rho1, rho2 = rho2, lambda0 = lambda0, 
                        enrl_timecut = enrl_timecut, enrl_rate = enrl_rate, drop_timecut = drop_timecut, 
                        drop_prob = drop_prob, seed = seed)
  
  # run the comparator analysis: no interim; after stage 1, use posterior
  # probability to make decision
  fa_result <- run_comparator_analysis(df, nstage1 = nstage1, a0 = a0, b0 = b0, 
                                       theta2 = theta2, delta3 = delta3, delta4 = delta4,
                                       orr_fu_time = orr_fu_time, nevent, 
                                       sig_level = sig_level)    
  
  
  return(fa_result)
  
}


#' Run a comparator design using all stage 1 data
#' For comparator analysis, stage 1 decision criteria
#'   Pos(theta > theta2) > delta4   -- success
#'   Pos(theta > theta2) < delta3   -- failure
#'   Otherwise                      -- Expand to two arms  
#' @param df the simulated data
#' @param nstage1 number of subjects needed in stage 1
#' @param a0,b0 beta prior 
#' @param theta2,delta3,delta4 as explained above 
#' @param orr_fu_time follow up time for ORR assessment
#' @param nevent Number of events needed to trigger stage 2 survival analysis
#' @param sig_level significance level for stage 2
#'
#' @return
#' @export
#'
#' @examples
run_comparator_analysis <- function(df, nstage1, a0, b0, theta2, delta3, delta4,
                                    orr_fu_time, nevent, sig_level){
  
  # run stage 1 analysis
  
  df_stage1 <- df %>% slice(1:nstage1)
  x0 <- sum(df_stage1$resp);
  pos1 <- pos_prob(x0, nstage1, a0 = a0, b0 = b0, p0 = theta2)
  
  stage1_decision <- dplyr::case_when(pos1 > delta4 ~ "Success",
                                      pos1 < delta3 ~ "Failure", 
                                      TRUE ~ "Expansion")
    
  s1_result <- tibble::tibble(nstage1 = nstage1, n_resp1 = x0, 
                              stage1_time = max(df_stage1$timein) + orr_fu_time, 
                              prob_pos = pos1, decision = stage1_decision)
  
  if (stage1_decision %in% c("Expansion")){
    
    fa_sub <- df %>% slice((nstage1+1):nrow(df)) %>% 
      mutate(t_pfs = os, t_os = os)
    fa_snapshot <- beverage::take_snapshot_event(surv_dat = fa_sub, events = nevent, type = "os")
    fa_time <- fa_snapshot$snapshot_time
    snapshot <- fa_snapshot$snapshot
    res0 <- run_survival(time = snapshot$os, censor = snapshot$os_censor,
                         arm = snapshot$arm, control = "SOC", alpha = sig_level)
    
    res <- res0 %>% mutate(fa_time = fa_time) %>% 
      select(fa_n = n, fa_time, fa_p = pvalue, event, hr, lower95, upper95, median_surv_1, median_surv_2) 
    
  } else{
    res <- NULL
  }
  
  fa1_result <- dplyr::bind_cols(s1_result, res)   
  
  return(fa1_result)
}







#' @title Run survival analysis
#' @details basic version of the function run_test, primary analysis only, no
#'   interim analysis
#' @param time time to event data
#' @param censor indicator for censor or event, 0 = censored, 1 = event
#' @param arm indicator for arm. Has to be at least two arms
#' @param alpha the significance level to claim a success for 1-sided test
#' @param control the name of the control arm
#' @import survival
#' @return a data frame containing results of sample size, p-value, decision,
#'   number of events, hazard ratio with 95% confidence interval,
#'   median survival, median follow-up, ect
#' @export
#' @seealso \code{\link[beverage]{take_snapshot_event}}
#'
#' @examples
#' rand_arm(nsbj = 1, ratio = c(1, 1))
#' time <- c(5.68, 0.34, 4.94, 1.49, 4.72, 1.32, 3.48, 3.42, 3.41, 2.93)
#' censor <- c(0, 0, 0, 1, 0, 1, 0, 0, 0, 1)
#' arm <- c("arm_1", "arm_2", "arm_1", "arm_2", "arm_2", "arm_1", "arm_2", "arm_1", "arm_1", "arm_2")
#' control <- "arm_1"
#' res1 <- run_survival(time, censor, arm, control = "arm_1", alpha = 0.025)
#'
run_survival <- function(time, censor, arm, control=NA, alpha = 0.025){
  
  library(survival)
  if(!is.na(control)){arm <- as.numeric(arm != control)}
  
  data <- tibble::tibble(time = time, censor = censor, arm = arm)
  
  # log-rank test p-value
  surv <- Surv(time,censor)~arm
  lr <- survdiff(surv,data=data)  # run log-rank test
  lr_pvalue <- pchisq(lr$chisq, length(lr$n)-1, lower.tail = FALSE) / 2 # obtain 1-sided p-value
  
  # cox HR with 95% confidence interval
  cox <- summary(coxph(surv,data=data))
  hr <- cox$coefficients[,"exp(coef)"]
  hr_lower95 <- cox$conf.int[,"lower .95"]
  hr_upper95 <- cox$conf.int[,"upper .95"]
  hr_ConfInt <- c(hr_lower95, hr_upper95)
  
  # adjust 1-sided p-value for control better than treatment case
  if (hr > 1){lr_pvalue <- 1 - lr_pvalue}
  
  # trial decision
  if (lr_pvalue < alpha & hr < 1) {decision <- 'success'} else{
    decision <- 'failure'
  }
  
  # median K-M survival
  km <- survfit(surv, type="kaplan-meier", conf.type="log", data=data)
  median_arms <- as.vector(quantile(km, 0.5)$quantile)
  
  median_diff <- median_arms[2] - median_arms[1]
  
  # subjects number
  n <- sum(km$n)
  
  # events number
  events <- sum(data$censor)
  
  # median follow-up time
  data$lfu_censor <- 1 - data$censor
  surv_rev <- Surv(time,lfu_censor)~arm
  follow_up <- survfit(surv_rev, type="kaplan-meier", conf.type="log", data=data)
  median_fu_arms <- as.vector(quantile(follow_up, 0.5)$quantile)
  median_fu_diff <- median_fu_arms[2] - median_fu_arms[1]
  
  # mean follow-up time
  temp <- survival:::survmean(follow_up, rmean=999)
  mean_fu_arms<- temp$matrix[, "rmean"]
  
  
  # output
  return(tibble::tibble(n=n, pvalue=lr_pvalue, decision=decision, event = events,
                        hr=hr, lower95=hr_lower95, upper95=hr_upper95,
                        median_surv_1=median_arms[1], median_surv_2=median_arms[2],median_surv_diff=median_diff,
                        median_fu_1=median_fu_arms[1], median_fu_2=median_fu_arms[2],median_fu_diff=median_fu_diff,
                        mean_fu_1=mean_fu_arms[1],  mean_fu_2=mean_fu_arms[2]))
}
