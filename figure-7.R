library(tidyverse)
library(cowplot)
library(janitor)
library(modelr)

age_groups <- read_csv("rawData_B6Aging_Auwerx.csv") %>% dplyr::select(EarTag,group,date_of_birth) %>%
  mutate(date_of_birth = as.Date(date_of_birth, "%Y.%m.%d"))
age_clams <- read_csv("CLAMS_B6Aging_LISP_RH-Auwerx.csv") %>% clean_names() %>% 
  group_by(ear_tag,experiment_started) %>% summarise_if(is.numeric, mean, na.rm=T) %>% filter(ear_tag != "EMPTY") %>% 
  left_join(age_groups, by=c("ear_tag"="EarTag")) %>% mutate(date = as.Date(experiment_started, "%d.%m.%Y"), 
                                                         age_in_weeks = round(difftime(date, date_of_birth, units = "weeks"),digits=0),
                                                         age_in_years = age_in_weeks/52,
                                                         age_in_weeks = factor(age_in_weeks), var = 1)

#7B
age_clams %>% mutate(age_in_weeks = factor(age_in_weeks,levels = c(15,55,94),labels=c("15 w","55 w","94 w")), heat_kcal_hr = heat_kcal_hr * 4.184) %>%
  ggplot(aes(x=subject_mass_g, y=heat_kcal_hr, color = age_in_weeks, fill=age_in_weeks)) +  geom_smooth(method = "lm", se=F, size = 2) +
  geom_point(pch=21,color="black",size=3) + xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)") + 
  theme(legend.position = c(0.05,0.85)) + labs(color = "age") + scale_color_manual(values = c('#a6bddb','#3690c0','#016450')) +
  scale_fill_manual(values = c('#a6bddb','#3690c0','#016450')) + guides(fill=F) +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

mass_3w <- read_csv("oneal_mass.csv") %>% gather("week","total_mass",-Subject)
ee_3w <- read_csv("oneal_ee.csv") %>% gather("week","ee",-Subject) %>% mutate(ee=as.numeric(ee), ee = ee/24) %>%
  left_join(mass_3w) %>% filter(week != "Habituation") %>% mutate(var = 1)
#7C
ee_3w %>% mutate(exercise = factor(week, labels = c("0 w exercise","1 w exercise","2 w exercise","3 w exercise")), ee = ee * 4.184) %>%
  ggplot(aes(x=total_mass,y=ee,color=exercise, fill = exercise)) + geom_smooth(method="lm",se=F,size=2) +
  geom_point(pch=21,color="black",size=3) + theme(legend.position = c(0.7,0.85)) + guides(fill=F) +
  scale_color_manual(values = c('#b3cde3','#8c96c6','#8856a7','#810f7c')) +
  scale_fill_manual(values = c('#b3cde3','#8c96c6','#8856a7','#810f7c')) +
  xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)") +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

pre_cl <- read_csv("CL_baseline.csv") %>% clean_names() %>% filter(time_of_day=="Full day") %>% mutate(CL=FALSE)
post_cl <- read_csv("CL_effect.csv") %>% clean_names() %>% filter(time_of_day=="Full day") %>% mutate(CL=TRUE)
ee_cl <- bind_rows(pre_cl, post_cl) %>% mutate(CL = factor(CL, levels = c(T,F), labels = c("CL at 30°C","30°C")), var = 1)
#7D
ee_cl %>% mutate(energy_expenditure = energy_expenditure * 4.184) %>% 
  ggplot(aes(x=total_mass, y=energy_expenditure, color=CL,fill=CL)) +
  geom_smooth(method="lm", se = F, size =2) +geom_point(pch=21,color="black",size=3) +
  xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)") + guides(fill=F) +
  theme(legend.position = c(0.05,0.85)) + scale_color_manual(values = c("#abd9e9","#d73027")) +
  scale_fill_manual(values = c("#abd9e9","#d73027")) +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

temp_bw <- read_csv("TempMassesBefore.csv") %>% clean_names()
temp_clams <- read_csv("TempMice.csv") %>% clean_names() %>% 
  filter(group=="unshaved mice", transition == F, menthol == F, olive_oil == F) %>%
  group_by(subject_id, enclosure_set_point_temperature) %>%
  summarise_if(is.numeric, mean, na.rm=T) %>% dplyr::select(-c(x1:cage)) %>% left_join(.,temp_bw, by=c("subject_id"="x1")) %>% 
  mutate(temp = factor(enclosure_set_point_temperature), var = 1) %>% ungroup()
#7E
temp_clams %>% mutate(temp=factor(temp,labels = c("6°C","10°C","14°C","18°C","22°C","25°C","28°C","30°C")), ee = ee * 4.184) %>% 
  ggplot(aes(x=total_mass, y=ee, color=temp,fill=temp)) + geom_smooth(method="lm",se=F, size = 2) + 
  labs(color = "temperature") + guides(fill=F)+ #xlim(c(22.4,28))+
  xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)") + 
  theme(legend.position = c(0.74,0.4)) + scale_x_continuous(breaks = c(23, 24, 25,26), limits = c(22.4,28)) +
  scale_color_brewer(palette="RdYlBu",direction = -1) + scale_fill_brewer(palette="RdYlBu",direction = -1) +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

age_delt <- age_clams %>% filter(age_in_weeks == 15) %>% 
  group_by(age_in_weeks) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% bind_rows(age_clams,.) %>%
  full_join(
    x = filter(., !is.na(ear_tag)),
    y = filter(., is.na(ear_tag)),
    by = "var",
    suffix = c("_exp", "_con")) %>%
  mutate(ee = heat_kcal_hr_exp - heat_kcal_hr_con) %>%
  group_by(age_in_weeks_exp) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  mutate(age_in_weeks_exp = factor(age_in_weeks_exp, levels = c(15,55,94),labels = c("15 w old","55 w old","94 w old")), exp = "age") %>%
  dplyr::select(ee, var, exp, group=age_in_weeks_exp)

ex_delt <- ee_3w %>% filter(week == "Baseline") %>% mutate(Subject = as.character(Subject)) %>%
  group_by(week) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% bind_rows(ee_3w,.) %>% 
  full_join(
    x = filter(., !is.na(Subject)),
    y = filter(., is.na(Subject)),
    by = "var",
    suffix = c("_exp", "_con")) %>%
  mutate(ee = ee_exp - ee_con) %>% 
  group_by(week_exp) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  mutate(week_exp = factor(week_exp, labels = c("0 w exercise","1 w exercise","2 w exercise","3 w exercise")), exp = "exercise") %>%
  dplyr::select(ee, var, exp, group=week_exp)

cl_delt <- pre_cl %>% mutate(var=1, CL = "30°C") %>% 
  group_by(CL) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% bind_rows(ee_cl,.) %>%
  full_join(
    x = filter(., !is.na(subject_id)),
    y = filter(., is.na(subject_id)),
    by = "var",
    suffix = c("_exp", "_con")) %>%
  mutate(ee = energy_expenditure_exp - energy_expenditure_con) %>% 
  group_by(CL_exp) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)  %>% mutate(exp = "CL") %>%
  dplyr::select(ee, var, exp, group=CL_exp)

temp_delt <- temp_clams %>% mutate(var = 1) %>% 
  mutate(subject_id = as.character(subject_id)) %>% filter(temp == 22) %>% 
  group_by(temp) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% bind_rows(temp_clams,.) %>% 
  full_join(
    x = filter(., !is.na(subject_id)),
    y = filter(., is.na(subject_id)),
    by = "var",
    suffix = c("_exp", "_con")) %>%
  mutate(ee = ee_exp - ee_con) %>%
  group_by(temp_exp) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  mutate(exp = "temperature", temp_exp =factor(temp_exp,labels = c("6°C","10°C","14°C","18°C","22°C","25°C","28°C","30°C"))) %>%
  dplyr::select(ee, var, exp, group=temp_exp)

delt_ee <- bind_rows(age_delt,ex_delt,cl_delt,temp_delt) %>%
  unite("x", exp, group, remove = F)

impc_df <- read.csv("impc_df.csv")
impc_m <- impc_df %>% filter(sex=="male") 
impc_wt_m<-impc_m %>% filter(group =="control" & ee >0.3*4.184) 
impc_ko_m<-impc_m %>% filter(group =="experimental" & strain != "Fbxl19") 
ko_model <- lm(ee ~ total_mass_1 + room_temp, data = impc_wt_m)
ko_res <- impc_ko_m %>% add_residuals(model=ko_model)
ko_res_mean <- impc_ko_m %>% add_residuals(model=ko_model) %>% 
  group_by(strain, zygosity) %>% summarise_if(is.numeric, mean, na.rm=T)

cols <- c('#a6bddb','#3690c0','#016450',"#d73027","#abd9e9",
          '#b3cde3','#8c96c6','#8856a7','#810f7c',
          '#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43','#d73027','#4575b4')

#7A
delt_ee %>% mutate(exp = factor(exp, levels = c("age","exercise","CL","temperature")), ee = ee * 4.184) %>% group_by(exp) %>% 
  ggplot(data=., aes(y = ee, x = var, fill = x),position = position_jitter(w = 0.05, h = 0)) +
  geom_hline(yintercept=0,linetype="dotted", color = "black") +
  geom_hline(yintercept=sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_hline(yintercept=-sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_hline(yintercept=2*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_hline(yintercept=-2*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_hline(yintercept=3*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_hline(yintercept=-3*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_point(pch=21,color="black",size=4) + scale_fill_manual(values=cols) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  ylab("∆ energy expenditure (kJ/hr)") + xlab(NULL) + facet_wrap(~exp, nrow = 1) +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "∆ energy expenditure (kcal/hr)"))
