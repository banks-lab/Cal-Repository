library(tidyverse)
library(cowplot)
library(relaimpo)
library(broom)
library(lsmeans)

impc_df <- read.csv("impc_df.csv")
impc_wt <- impc_df %>% filter(group=="control") %>% mutate(sex=factor(sex, levels = c("male","female")))
impc_3 <- impc_wt %>% filter(location %in% c("Germany","Oxfordshire, UK","Canada") & !is.na(lean_mass)) %>%
  mutate(temp_range = as.character(cut(room_temp, breaks=c(18,22,26,29),labels=c("18-22","22-26","26-29"))))

#4A
impc_3 %>% 
  ggplot(.,aes(y = ee, x = total_mass_1, color=sex, fill = sex)) +
  scale_fill_manual(values = c( "#00B0F6","#FF67A4")) +
  scale_color_manual(values = c( "#00B0F6","#FF67A4")) +
  theme(legend.title=element_blank(),legend.position = c(.05,.9)) +
  xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)")+xlim(15,37.5)+
  geom_point(colour="black",pch=21, size=2)+
  geom_smooth(method = "lm",se=F, aes(color=sex), size=2) + guides(fill=F) +
  scale_y_continuous(limits = c(.3*4.184,.7*4.184), sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

impc_3 %>%
  group_by(sex) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(ee ~ total_mass_1, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance)
  ) %>% 
  unnest(tidied)
sex_lm <- lm(ee ~ total_mass_1*sex, data = impc_3)
anova(sex_lm)
sex_ls <- lstrends(sex_lm,"sex", var="total_mass_1")
pairs(sex_ls)
cld(sex_ls)

#4B
impc_3 %>% mutate(location = factor(location, levels = c("Oxfordshire, UK","Germany","Canada"))) %>%
  ggplot(.,aes(y = ee, x = total_mass_1)) +
  scale_fill_manual(values = c('#b2df8a', '#1f78b4','#ff7f00')) +
  scale_color_manual(values = c('#b2df8a', '#1f78b4','#ff7f00')) +
  theme(legend.title=element_blank(),legend.position = c(.05,.9)) +
  xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)")+xlim(15,37.5)+
  geom_point(aes(fill=location),colour="black",pch=21, size=2)+
  geom_smooth(method = "lm",se=F, aes(color=location), size=2) + guides(fill=F) +
  scale_y_continuous(limits = c(.3*4.184,.7*4.184), sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

impc_3 %>%
  group_by(location) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(ee ~ total_mass_1, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance)
  ) %>% 
  unnest(tidied)
loc_lm <- lm(ee ~ total_mass_1*location, data = impc_3)
anova(loc_lm)
loc_ls <- lstrends(loc_lm,"location", var="total_mass_1")
pairs(loc_ls)
cld(loc_ls)

#4C
impc_3 %>% 
  ggplot(.,aes(y = ee, x = total_mass_1)) +
  scale_color_manual(values=c("#4575b4", "#ffffbf","#d73027"))+
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027",midpoint = 23) +
  theme(legend.position = c(.05,.8)) + 
  labs(fill = "ambient\ntemp (°C)", color = NULL) +
  xlab("body mass (g)") + ylab("energy expenditure (kJ/hr)")+xlim(15,37.5)+
  geom_point(aes(fill=room_temp),colour="black",pch=21, size=2)+
  geom_smooth(method = "lm",se=F, aes(color=temp_range), size=2) + guides(fill = guide_colorbar(reverse = TRUE), color=F) +
  scale_y_continuous(limits = c(.3*4.184,.7*4.184), sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

impc_3 %>%
  group_by(temp_range) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(ee ~ total_mass_1, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance)
  ) %>% 
  unnest(tidied)
temp_lm <- lm(ee ~ total_mass_1*temp_range, data = impc_3) 
anova(temp_lm)
temp_ls <- lstrends(temp_lm,"temp_range", var="total_mass_1")
pairs(temp_ls)
cld(temp_ls)

i3_var <- lm(ee~sex+lean_mass+fat_mass+season+room_temp, data=impc_3, na.action = na.exclude)
i3_calc <- calc.relimp(i3_var)
round(i3_calc@R2 * 100,digits=2) #r2
i3_variance <-rownames_to_column(as.data.frame(i3_calc$lmg)) %>% 'colnames<-' (c("var","pct")) %>% mutate(pct = pct*100) %>%
  arrange(desc(pct))
#4D
ggplot(i3_variance,aes(reorder(var,-pct),pct))+
  geom_bar(stat="identity",aes(fill=-pct))+ 
  xlab(NULL)+ylab("% of energy expenditure variance") + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  annotate("text",label=expression(R^{2} ~ "= 67.20%"), x = 4, y=30, size = 8) +
  scale_x_discrete(labels=c("ambient\ntemp","lean mass","sex","fat mass","season"))
#4E
data.frame("pct"=c(99.9,0.1),"var"=c("other","institution"),"total"="total") %>% 
  mutate(var = factor(var, levels=c("other","institution"))) %>%
  ggplot(aes(x=total,y=pct,fill=var))+geom_bar(position="stack",stat="identity") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.position = "none") + 
  ylab("% of residual variance\nexplained by institution") +
  scale_fill_manual(values=c("gray","#132c43"))
