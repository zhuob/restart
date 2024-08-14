


summarize_sim <- function(combo_sim){
  
  result <- NULL
  
  for(k in unique(combo_sim$scenario)){
    
    tmp <- combo_sim %>% filter(scenario == k)
    
    tmp1 <- tmp %>% mutate(
      fa_success = case_when(ia_decision == "Expand to Two-arm" & fa_p <  sig_level ~ TRUE, 
                             ia_decision == "Expand to Two-arm" & fa_p >= sig_level ~ FALSE, 
                             TRUE ~ fa_success),
      fa_n = case_when(ia_decision == "Expand to Two-arm" ~ fa_n + ia_n_enrl, 
                       TRUE ~ fa_n)
    ) %>% mutate(
      success_type = case_when(ia_decision == "Continue Single Arm" & fa_success ~ "Single Arm Success", 
                               ia_decision == "Continue Single Arm" & !fa_success ~ "Single Arm Failure", 
                               ia_decision == "Expand to Two-arm" & !fa_success ~ "Two arm Expansion Failure",
                               ia_decision == "Expand to Two-arm" & fa_success ~ "Two arm Expansion Success",
                               ia_decision == "Futile"  & fa_success ~ "Futile then Success", 
                               ia_decision == "Futile"  & !fa_success ~ "Futile then Failure")
    )
    
    # calculate study duration w.r.t. IA decisions
    fa_duration <- tmp1 %>% group_by(ia_decision) %>% 
      summarize(ia_time = mean(ia_time), fa_time = mean(fa_time))
    
    fa_summary <- tmp1 %>% summarize(
      ia_time = mean(ia_time), 
      ia_n = mean(ia_n_enrl), 
      prob_ia_1arm = mean(ia_decision == "Continue Single Arm"),
      prob_ia_2arm = mean(ia_decision == "Expand to Two-arm"), 
      prob_ia_fut  = mean(ia_decision == "Futile"),
      orr_est = mean(ia_n_resp/ia_n_eval),
      fa_time = mean(fa_time), fa_n = mean(fa_n), 
      prob_fa_success = mean(fa_success),
      prob_suc_ia_1arm = mean(success_type == "Single Arm Success"), 
      prob_suc_ia_2arm = mean(success_type == "Two arm Expansion Success"), 
      prob_suc_ia_fut  = mean(success_type == "Futile then Success")
    ) %>% mutate(scenario = k, .before = 1) %>% 
      mutate(fa_time_1arm = fa_duration$fa_time[fa_duration$ia_decision == "Continue Single Arm"], 
             fa_time_2arm = fa_duration$fa_time[fa_duration$ia_decision == "Expand to Two-arm"], 
             fa_time_fut = fa_duration$fa_time[fa_duration$ia_decision == "Futile"])
    
    result <- bind_rows(result, fa_summary)
  }    
  
  return(result)
  
}


plot_prob_succ <- function(s0){
  
  relation_order <- c("Weak", "Moderate", "Strong")
  success_order <- c("Single Arm Success", "Two Arm Expansion Success", "Futile then Success")
  tmp <- s0 %>% mutate(orr1 = factor(orr1, levels = unique(s0$orr1)), 
                       `Relationship between ORR and OS` = factor(`Relationship between ORR and OS`, levels = relation_order))
  
  tmp1 <- tmp %>% select(enrl_rate, n_ia, orr1, starts_with("prob_suc_ia"), 
                         `Relationship between ORR and OS`) %>% 
    pivot_longer(cols = 4:6, names_to = "success_type", values_to = "prob_suc") %>%
    mutate(`Success Type` = case_when(success_type == "prob_suc_ia_1arm" ~ success_order[1], 
                                      success_type == "prob_suc_ia_2arm" ~ success_order[2], 
                                      success_type == "prob_suc_ia_fut"  ~ success_order[3])) %>%
    filter(success_type != "prob_suc_ia_fut") %>%
    mutate(`Success Type` = factor(`Success Type`, levels = success_order)) 
  
  
  tmp1 %>% group_by(enrl_rate, n_ia, orr1, `Relationship between ORR and OS` ) %>% mutate(cum_prob = cumsum(prob_suc)) %>% 
    ggplot(aes(orr1, cum_prob, fill = `Relationship between ORR and OS` )) + 
    geom_col(data = . %>% filter(`Success Type` == success_order[1]),  position = position_dodge(width = 0.9), alpha = 1) + 
    geom_col(data = . %>% filter(`Success Type` == success_order[2]),  position = position_dodge(width = 0.9), alpha = 0.6) + 
    geom_col(data = . %>% filter(`Success Type` == success_order[3]),  position = position_dodge(width = 0.9), alpha = 0.2) + 
    geom_tile(aes(y = NA_integer_, alpha = `Success Type`)) + 
    scale_alpha_manual(values = c(1, 0.7, 0.4)) + 
    facet_wrap(n_ia ~ enrl_rate) + 
    theme(legend.position = "bottom", legend.box = "vertical") + 
    labs(x = "Response Rate", y = "Probability of Success")  
  
  
}


plot_prob_ia <- function(s0){
  
  relation_order <- c("Weak", "Moderate", "Strong")
  ia_order <- c("Continue Single Arm", "Expand to Two-arm", "Futile")
  tmp <- s0 %>% mutate(orr1 = factor(orr1, levels = unique(s0$orr1)), 
                       `Relationship between ORR and OS` = factor(`Relationship between ORR and OS`, levels = relation_order))
  
  tmp1 <- tmp %>% select(enrl_rate, n_ia, orr1, starts_with("prob_ia"), 
                         `Relationship between ORR and OS`) %>% 
    pivot_longer(cols = 4:6, names_to = "ia_type", values_to = "probs") %>%
    mutate(`Interim Decision` = case_when(ia_type == "prob_ia_1arm" ~ ia_order[1], 
                                      ia_type == "prob_ia_2arm" ~ ia_order[2], 
                                      ia_type == "prob_ia_fut"  ~ ia_order[3])) %>%
    mutate(`Interim Decision` = factor(`Interim Decision`, levels = ia_order)) 
  
  
  tmp1 %>% group_by(enrl_rate, n_ia, orr1, `Relationship between ORR and OS` ) %>% mutate(cum_prob = cumsum(probs)) %>% 
    ggplot(aes(orr1, cum_prob, fill = `Relationship between ORR and OS` )) + 
    #theme_bw() +
    geom_col(data = . %>% filter(`Interim Decision` == ia_order[1]),  position = position_dodge(width = 0.9), alpha = 1) + 
    geom_col(data = . %>% filter(`Interim Decision` == ia_order[2]),  position = position_dodge(width = 0.9), alpha = 0.7) + 
    geom_col(data = . %>% filter(`Interim Decision` == ia_order[3]),  position = position_dodge(width = 0.9), alpha = 0.4) + 
    geom_tile(aes(y = NA_integer_, alpha = `Interim Decision`)) + 
    scale_alpha_manual(values = c(1, 0.7, 0.4)) + 
    facet_wrap(n_ia ~ enrl_rate) + 
    theme(legend.position = "bottom", legend.box = "vertical") + 
    labs(x = "Response Rate", y = "Probability")  
  
  
}


plot_duration_ia <- function(s0, eval_n_ia = 10){
  
  relation_order <- c("Weak", "Moderate", "Strong")
  ia_order <- c("Continue Single Arm", "Expand to Two-arm", "Futile")
  tmp <- s0 %>% mutate(orr1 = factor(orr1, levels = unique(s0$orr1)), 
                       `Relationship between ORR and OS` = factor(`Relationship between ORR and OS`, levels = relation_order))
  
  tmp1 <- tmp %>% select(enrl_rate, n_ia, orr1, starts_with("fa_time_"), 
                         `Relationship between ORR and OS`) %>% 
    pivot_longer(cols = 4:6, names_to = "ia_type", values_to = "duration") %>%
    mutate(`Interim Decision` = case_when(ia_type == "fa_time_1arm" ~ ia_order[1], 
                                          ia_type == "fa_time_2arm" ~ ia_order[2], 
                                          ia_type == "fa_time_fut"  ~ ia_order[3])) %>%
    mutate(`Interim Decision` = factor(`Interim Decision`, levels = ia_order)) 
  
  
  tmp2 <- tmp1 %>% filter(n_ia == eval_n_ia)
  
  ggplot(data = tmp2, aes(x = orr1, y = duration, fill = `Relationship between ORR and OS`)) + 
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(enrl_rate ~ `Interim Decision`) +
    theme(legend.position = "bottom") + 
    labs(x = "Response Rate", y = "Study Duration (months)", title = paste("IA Sample Size: ", eval_n_ia))
  
  
}



plot_sim_result <- function(s0, var_of_interest = "prob_fa_success", 
                            labely = "Probability of Success"){
  
  relation_order <- c("Weak", "Moderate", "Strong")
  tmp <- s0 %>% mutate(orr1 = factor(orr1, levels = unique(s0$orr1)), 
                       `Relationship between ORR and OS` = factor(`Relationship between ORR and OS`, levels = relation_order))
  
  ggplot(data = tmp, aes(x = orr1, y = !!rlang::sym(var_of_interest), 
                         fill = `Relationship between ORR and OS`)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(n_ia ~ enrl_rate ) + 
    theme(legend.position = "bottom") + 
    labs(x = "Response Rate", y = labely)  
  
  
}


summarize_sim_comparator <- function(combo_sim){
  
  result <- NULL
  
  for(k in unique(combo_sim$scenario)){
    t0 <- combo_sim %>% filter(scenario == k)
    t1 <- t0 %>% mutate(fa_success = case_when(
      decision == "Failure" ~ F, 
      decision == "Success" ~ T, 
      decision == "Expansion" & fa_p >= sig_level ~ FALSE,
      decision == "Expansion" & fa_p <  sig_level ~ TRUE
    )) %>% 
      mutate(fa_n = case_when(decision == "Expansion" ~ fa_n + nstage1, 
                              TRUE ~ nstage1), 
             success_type = case_when(
               decision == "Expansion" & fa_success ~ "Expansion then Success", 
               decision == "Expansion" & !fa_success ~ "Expansion then Failure", 
               TRUE  ~ decision
             ), 
             fa_time = case_when(decision  %in% c("Failure", "Success") ~ stage1_time, 
                                 TRUE ~ fa_time)) 
    
    # calculate study duration w.r.t. IA decisions
    fa_duration <- t1 %>% group_by(decision) %>% 
      summarize(duration = mean(fa_time))
    durfail <- fa_duration$duration[fa_duration$decision == "Failure"]
    
    fa_summary <- t1 %>% summarize(
      prob_1arm = mean(decision %in% c("Failure", "Success")),
      prob_1arm_suc = mean(decision == "Success"), 
      prob_1arm_fut = mean(decision == "Failure"),
      prob_2arm = mean(decision == "Expansion"), 
      fa_time = mean(fa_time), fa_n = mean(fa_n), 
      prob_fa_success = mean(fa_success),
      prob_suc_1arm = mean(success_type == "Success"), 
      prob_suc_2arm = mean(success_type == "Expansion then Success")
    ) %>% mutate(scenario = k, .before = 1) %>% 
      mutate(fa_1arm_suc = fa_duration$duration[fa_duration$decision == "Success"], 
             fa_2arm = fa_duration$duration[fa_duration$decision == "Expansion"], 
             fa_1arm_fut =ifelse(length(durfail)==0, NA, durfail))
    
    result <- bind_rows(result, fa_summary)
    
  }
  
  return(result)
}


library(ggpattern)
plot_prob_succ_comparator <- function(s0, n0_ia = 40){
  
  tmp <- s0 %>% filter(n_ia == n0_ia) %>% 
    select(orr1, enrl_rate, 15:27) %>% 
    pivot_longer(cols = 4:15, names_to = "Type", values_to = "value") %>% 
    mutate(type0 = case_when(
      Type == "prob_1arm" ~ "Single Arm", 
      Type == "prob_1arm_suc" ~ "Single Arm Success",
      Type == "prob_1arm_fut" ~ "Single Arm Failure",
      Type == "prob_2arm" ~ "Expansion",
      Type == "fa_time" ~ "Study Duration",
      Type == "fa_n" ~ "Sample Size",
      Type == "prob_fa_success" ~ "Prob of Final Success",
      Type == "prob_suc_1arm" ~ "Prob. of Success (Single Arm)",
      Type == "prob_suc_2arm" ~ "Prob. of Success (Expansion)",
      Type == "fa_1arm_suc" ~ "Duration (Single Arm Success)",
      Type == "fa_2arm" ~ "Duration (Expansion)",
      Type == "fa_1arm_fut" ~ "Duration (Single Arm Failure)",
    ),
    orr1 = factor(orr1, levels = unique(s0$orr1)), 
    `Relationship between ORR and OS` = 
      factor(`Relationship between ORR and OS`,
             levels = c("Weak", "Moderate", "Strong")))
  
  tmp0 <- tmp %>% filter(Type %in% c("fa_time"))
  ggplot(data = tmp0, aes(x = orr1, y = value, pattern = Type)) + 
    geom_col_pattern(pattern_size = 0.25) + 
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(enrl_rate ~ `Relationship between ORR and OS`) + 
    theme(legend.position = "bottom", legend.box = "vertical") + 
    labs(x = "Response Rate", y = "Duration")
  
}

