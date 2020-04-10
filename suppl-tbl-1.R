library(tidyverse)
library(modelr)

impc_df <- read.csv("impc_df.csv")

act_wt <- impc_df %>% filter(group=="control" & !is.na(loc_act))
no_act_wt <- impc_df %>% filter(group == "control" & is.na(loc_act))

act_ko <- impc_df %>% filter(group=="experimental" & !is.na(loc_act))
no_act_ko <- impc_df %>% filter(group == "experimental" & is.na(loc_act))

act_m_wt <- act_wt %>% filter(sex == "male")
act_f_wt <- act_wt %>% filter(sex == "female")

act_m_ko <- act_ko %>% filter(sex == "male")
act_f_ko <- act_ko %>% filter(sex == "female")

no_act_m_wt <- no_act_wt %>% filter(sex == "male")
no_act_f_wt <- no_act_wt %>% filter(sex == "female")

no_act_m_ko <- no_act_ko %>% filter(sex == "male")
no_act_f_ko <- no_act_ko %>% filter(sex == "female") %>% filter(institution != "China")

act_m_model <- lm(ee ~ total_mass_1 + room_temp + loc_act, data = act_m_wt)
act_m_res <- act_m_ko %>% add_residuals(model=act_m_model) %>% 
  group_by(institution, strain, zygosity, sex) %>% summarize_if(is.numeric, mean, na.rm=T) %>%
  mutate(sd = resid/sd(.$resid, na.rm = T))

act_f_model <- lm(ee ~ total_mass_1 + room_temp + loc_act, data = act_f_wt)
act_f_res <- act_f_ko %>% add_residuals(model=act_f_model) %>% 
  group_by(institution, strain, zygosity, sex) %>% summarize_if(is.numeric, mean, na.rm=T) %>%
  mutate(sd = resid/sd(.$resid, na.rm = T))

no_act_m_model <- lm(ee ~ total_mass_1 + room_temp, data = no_act_m_wt)
no_act_m_res <- no_act_m_ko %>% add_residuals(model=no_act_m_model) %>% 
  group_by(institution, strain, zygosity, sex) %>% summarize_if(is.numeric, mean, na.rm=T) %>%
  mutate(sd = resid/sd(.$resid, na.rm = T))

no_act_f_model <- lm(ee ~ total_mass_1 + room_temp, data = no_act_f_wt)
no_act_f_res <- no_act_f_ko %>% add_residuals(model=no_act_f_model) %>% 
  group_by(institution, strain, zygosity, sex) %>% summarize_if(is.numeric, mean, na.rm=T) %>%
  mutate(sd = resid/sd(.$resid, na.rm = T))

m_res_table <- bind_rows(act_m_res, no_act_m_res) %>% arrange(desc(sd))
f_res_table <- bind_rows(act_f_res, no_act_f_res) %>% arrange(desc(sd))
