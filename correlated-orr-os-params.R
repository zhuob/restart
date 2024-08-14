####### matrix 
#       Non-Responder      Responder
# SOC    lambda              rho1*lambda
# COmbo  rho2*lambda         rho1*rho2*lambda
#
# rho1: hazard ratio for the responder and non-responder
# rho2: additional hazard ratio treatment effect independent of ORR 
# lambda: hazard 
# we can specify a range of relationship between ORR and OS by varying rho1 and rho2


# given median OS and ORR for control (1) and treatment (2), respectively, derive 
# the parameters: lambda0, rho1, rho2 

find_lambda <- function(orr0, rho1, mos0){
  
  obj_fun <- function(lambda0){
    (1-exp(-lambda0*rho1*mos0))*orr0 + (1-exp(-lambda0*mos0))*(1-orr0) - 0.5
  }
  uniroot(f = obj_fun, interval = c(0, 1), tol = 1e-5)$root
  
}


find_rho2 <- function(orr1, lambda0, rho1, mos1){
  
  obj_fun2 <- function(rho2){
    (1-exp(-lambda0*rho1*rho2*mos1))*orr1 + (1-exp(-lambda0*rho2*mos1))*(1-orr1) - 0.5
  }
  result <- uniroot(f = obj_fun2, interval = c(0, 2), tol = 1e-5)$root
  result
}


# given orr0, orr1, mos0, mos1 and one of rho1, rho2, solve the equations

find_rho1 <- function(orr0, orr1, mos0, mos1, rho2){



   rho1_lambda <- function(lambda0){
     t1 <- log(1 - (1/orr1)*(0.5 - (1-orr1)*(1-exp(-lambda0*mos1*rho2)) ))/(mos1*rho2)
     t2 <- log(1 - (1/orr0)*(0.5 - (1-orr0)*(1-exp(-lambda0*mos0))      ))/mos0

     t1 - t2
   }

   lambda0 <- uniroot(f = rho1_lambda, interval = c(0, 10), tol = 1e-5)$root
   # t1 <- log(1 - (0.5 - (1-orr1)*(1-exp(-lambda0*mos1)))/orr1)/mos1
   t2 <- log(1 - (0.5 - (1-orr0)*(1-exp(-lambda0*mos0)))/orr0)/mos0

   rho1 <-  -t2/lambda0

}




# the hazard function for control group (note, the hazard rate is not constant)
hz_fun1 <- function(lambda0, rho1, time, orr0){
  
  # density function
  ft <- (1-orr0)*lambda0*exp(-lambda0*time) + orr0*rho1*lambda0*exp(-rho1*lambda0*time)
  # survival function
  st <- (1-orr0)*exp(-lambda0*time) + orr0*exp(-rho1*lambda0*time)
  # hazard function
  ft/st
  
}

# the hazard function for treatment group
hz_fun2 <- function(lambda0, rho1, rho2, time, orr1){
  
  ft <- (1-orr1)*rho2*lambda0*exp(-rho2*lambda0*time) + orr1*rho2*rho1*lambda0*exp(-rho2*rho1*lambda0*time)
  st <- (1-orr1)*exp(-rho2*lambda0*time) + orr1*exp(-rho2*rho1*lambda0*time)
  
  ft/st
  
}

survival_prob_by_response <- function(rho1, rho2, lambda0, t, responder = "yes"){
  # suppose randonmization ratio is 1:1
  if(responder == "yes"){
    st <- 0.5*exp(-rho1*lambda0*t) + 0.5*exp(-rho1*rho2*lambda0*t)
  } else{
    st <- 0.5*exp(-lambda0*t) + 0.5*exp(-rho2*lambda0*t)
  }
  
  return(st)
}

survival_prob_by_arm <- function(rho1, rho2, lambda0, t, arm = "SOC"){
  
  # survival function
  if(arm == "TRT"){
    st <- (1-orr1)*exp(-rho2*lambda0*t) + orr1*exp(-rho2*rho1*lambda0*t)
  } else if (arm == "SOC"){
    st <- (1-orr0)*exp(-lambda0*t) + orr0*exp(-rho1*lambda0*t)
  }
  
  return(st)
}


get_quantity <- function(rho1, rho2, lambda0, orr0, orr1){
  
  times <- seq(0.1, 60, by = 0.1)
  # hazard function for SOC
  hz1 <- hz_fun1(lambda0 = lambda0, rho1, times, orr0)
  # hazard function for TRT
  hz2 <- hz_fun2(lambda0 = lambda0, rho1, rho2, times, orr1)
  
  # survival prob for responder
  surv1 <- survival_prob_by_response(rho1, rho2, lambda0, t = times, responder = "yes")
  surv2 <- survival_prob_by_response(rho1, rho2, lambda0, t = times, responder = "no")
  
  # survival prob for arms
  surv3 <- survival_prob_by_arm(rho1, rho2, lambda0, t = times, arm = "SOC")
  surv4 <- survival_prob_by_arm(rho1, rho2, lambda0, t = times, arm = "TRT")
  
  r0 <- tibble(time = times, hz_soc = hz1, hz_trt = hz2, 
               surv_resp = surv1, surv_noresp = surv2, 
               surv_soc = surv3, surv_trt = surv4)
  
  return(r0)
}



find_lambda_rho2 <- function(orr0, orr1, mos0, mos1, rho1){
  
  lambda0 <- find_lambda(orr0, rho1, mos0)
  
  rho2 <- find_rho2(orr1, lambda0, rho1, mos1)  
  rho2 <- ifelse(rho2 > 1, 1, ifelse(rho2 < 0, 0, rho2))
  
  mos <- log(2)/c(lambda0, rho1*lambda0, rho2*lambda0, rho1*rho2*lambda0)
  # names(mos) <- c("SOC/Non-responder", "SOC/Responder", "TRT/Non-Responder", "TRT/Responder")
  
  outcome_profile <- tibble(rho1, rho2, lambda0, 
                            `SOC/Non-Responder` = mos[1], 
                            `SOC/Responder`     = mos[2], 
                            `TRT/Non-Responder` = mos[3], 
                            `TRT/Responder`     = mos[4])
  
  return(outcome_profile)
  
  
}


plot_orr_survival <- function(rho1, orr0, orr1, mos0, mos1){
  
  
  res0 <- find_lambda_rho2(orr0, orr1, mos0, mos1, rho1)
  
  t1 <- get_quantity(rho1, rho2 = res0$rho2, lambda0 = res0$lambda0, orr0, orr1) %>% 
    mutate(hr = hz_trt/hz_soc)
  
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
  obj <- gridExtra::grid.arrange(p1, p2, p3, p4)

  return(obj)
  
}

