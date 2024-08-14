# organize simulation result
library(tidyverse)
source("C:/Users/zhuobin/OneDrive - Boehringer Ingelheim/personal/restart-paper/code/simulation-summary.R")

sim1 <- read_csv("sim-result/sim-v2.csv")
sim2 <- read_csv("sim-result/sim-v3-delta3-0.80-delta4-0.98.csv")
params <- read_csv("sim-result/params.csv")

sig_level <- 0.05
sim1_summary <- summarize_sim(combo_sim = sim1)
sim2_summary <- summarize_sim_comparator(combo_sim = sim2)
names(sim2_summary) <- paste("cp", names(sim2_summary), sep = "_")

sim_result <- left_join(sim1_summary, sim2_summary, by =c("scenario" = "cp_scenario"))
sim_result <- left_join(params, sim_result, by = "scenario")

col1 <- "Treatment ORR";
col2 <- "Correlation btw. ORR and OS"

# probability of interim decisions
df1 <- sim_result %>% select(orr1, n_ia, enrl_rate, `Relationship between ORR and OS`, 
                             starts_with("prob_ia"), cp_prob_1arm_suc, cp_prob_1arm_fut, cp_prob_2arm)
df1 <- df1  %>% pivot_wider(names_from = 2:3, values_from = 5:10)

df1 <- df1 %>% select(orr1, `Relationship between ORR and OS`, ends_with("_40_3")) %>% 
      select(c(1:6, 8,7))
library(xtable)
 df1_latex <- xtable(df1, caption = "Interim Decision Distribution", 
                     align = "l{3cm}l{3cm}r{1cm}r{1cm}r{1cm}|r{1cm}r{1cm}r{1cm}r{1cm}", 
                    digits = c(1, 2, 3, 3, 3, 3, 3, 3, 3))
colnames(df1_latex) <- c(col1, col2, 
                         "P(Single arm)", "P(Expansion)", "P(Futility)", 
                         "P(Success)", "P(Expansion)", "P(Failure)")
print(df1_latex, include.rownames = FALSE)


# Prob of Success
df2 <- sim_result %>% select(orr1, n_ia, enrl_rate, `Relationship between ORR and OS`, 
                             prob_fa_success, cp_prob_fa_success)
df2 <- df2 %>% pivot_wider(names_from = 2:3, values_from = 5:6)
df2 <- df2 %>% select(c(1,2,3,5,7))
df2_latex <- xtable(df2, caption = "Probability of Success (irrespective of Decision)",
                    align = "rr{2cm}r{3cm}|r{1cm}r{1cm}r{1cm}",
                    digits = c(1, 2, 3, 3, 3, 3))
colnames(df2_latex) <- c(col1, col2, "IA at 10", "IA at 40", "Tratidional")
print(df2_latex, include.rownames = FALSE)


# sample size 
df3 <- sim_result %>% select(orr1, n_ia, enrl_rate, `Relationship between ORR and OS`,
                             fa_n, cp_fa_n)
df3 <- df3 %>% pivot_wider(names_from = 2:3, values_from = 5:6) 
df3 <- df3 %>% select(c(1,2,3,5,7))
df3_latex <- xtable(df3, caption = "Sample size comparison", align = "rr{2cm}|r{3cm}rrr",
                    digits = c(1, 2, 3, 0, 0, 0))
colnames(df3_latex) <- c(col1, col2, "N(10)", "N(40)", "N")
print(df3_latex, include.rownames = FALSE)


# study duration 
df4 <- sim_result %>% select(orr1, n_ia, enrl_rate, `Relationship between ORR and OS`, 
                             fa_time, cp_fa_time)
df4 <- df4 %>% pivot_wider(names_from = 2:3, values_from = 5:6)
df4 <- df4 %>% select(c(1,2,3,5,7))
df4_latex <- xtable(df4, caption = "Study Duration comparison (months), assuming enrollment speed is 3 patients per month", 
                    align = "rr{2cm}|r{3cm}r{1cm}r{1cm}r{1cm}",
                    digits = c(1, 2, 3, 1, 1, 1))
colnames(df4_latex) <- c(col1, col2, "N(10)", "N(40)", "N")
print(df4_latex, include.rownames = FALSE)




## design parameters

p1 <- params %>% select( `Relationship between ORR and OS`, orr0, orr1, 
                        mos0, mos1, rho1, rho2, lambda0) %>% 
  distinct(orr0, orr1, `Relationship between ORR and OS`, .keep_all = TRUE)
p1_latex <- xtable(p1, caption = "Relationship between ORR and survival time", 
                   align = "rp{3cm}|r{1cm}r{1cm}r{1cm}r{1cm}r{1cm}r{1cm}r{1cm}", 
                   digits = c(1, 1, 2, 2, 1, 1, 2, 2, 3))
colnames(p1_latex) <- c(col2, "p0", "p1", "m0", "m1", 
                        "rho1", "rho2", "lambda_0")
print(p1_latex, include.rownames = FALSE)



## plot the study duration across different timing of IA
sig_level <- 0.05
params <- read_csv("sim-result/params.csv") 
params2 <- read_csv("sim-result/params-20-30-50.csv")
sim1_summary <- summarize_sim(combo_sim = read_csv("sim-result/sim-v2.csv")) # 10/40
sim2_summary <- summarize_sim(combo_sim = read_csv("sim-result/sim-v5-start-ia-20-30-50.csv")) # 20/30/50

sim3_summary <- summarize_sim_comparator(combo_sim =  read_csv("sim-result/sim-v3-delta3-0.80-delta4-0.98.csv")) # traditional, 60

sim1 <- left_join(params, sim1_summary, by = "scenario") %>%  
  select(scenario, orr1, `Relationship between ORR and OS`, n_ia, enrl_rate, ia_time, fa_time, fa_n)

sim2 <- left_join(params2, sim2_summary, by = "scenario") %>%  
  select(scenario, orr1, `Relationship between ORR and OS`, n_ia, enrl_rate, ia_time, fa_time, fa_n)

sim3 <- left_join(params, sim3_summary, by = "scenario") %>%
  select(scenario, orr1, `Relationship between ORR and OS`, enrl_rate, fa_time, fa_n) %>% 
  distinct(orr1, `Relationship between ORR and OS`, enrl_rate, .keep_all = TRUE)

sim_duration <- bind_rows(sim1, sim2, sim3 %>% mutate(n_ia = 60)) %>% 
  mutate(`Relationship between ORR and OS` = factor(`Relationship between ORR and OS`, levels = c("Strong", "Moderate", "Weak")), 
         orr1 = paste("ORR =", orr1))


dur0 <- sim_duration %>% filter(enrl_rate == 12)
dur1 <- dur0 %>% filter(n_ia != 60)
dur2 <- dur0 %>% filter(n_ia == 60)


ggplot(data = dur1, aes(x = n_ia, y = fa_time)) + 
  geom_line(aes(linetype = `Relationship between ORR and OS`), size = 0.8) + 
  geom_point(data = dur0, aes(x = n_ia, y = fa_time, shape = `Relationship between ORR and OS`), size = 2) + 
  facet_wrap(~orr1) + 
  scale_x_continuous(breaks = unique(sim_duration$n_ia), limits = c(0, max(dur0$n_ia))) + 
  scale_y_continuous(breaks = seq(0, 150, by = 10)) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 15, face = "bold"), 
        axis.text = element_text(size = 10, face = "bold"), 
        legend.key.width = unit(4, "cm"), 
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = "Number of Subjects for Analysis", y = "Study Duration (months)")

ggsave(filename = "sim-result/duration-enrl-12.png", width = 12, height = 8)

ggplot(data = dur1, aes(x = n_ia, y = fa_n)) + 
  geom_line(aes(linetype = `Relationship between ORR and OS`), size = 0.8) + 
  geom_point(data = dur0, aes(x = n_ia, y = fa_n, shape = `Relationship between ORR and OS`), size = 2) + 
  facet_wrap(~orr1) + 
  scale_x_continuous(breaks = unique(sim_duration$n_ia), limits = c(0, max(dur0$n_ia))) + 
  scale_y_continuous(breaks = seq(0, 150, by = 10)) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 15, face = "bold"), 
        axis.text = element_text(size = 10, face = "bold"), 
        legend.key.width = unit(4, "cm"), 
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = "Number of Subjects for Analysis", y = "Total Sample Size")

ggsave(filename = "sim-result/sample-size-12.png", width = 12, height = 8)
