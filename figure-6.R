library(tidyverse)
library(cowplot)
library(modelr)
library(viridis)
library(ggrepel)

impc_df <- read.csv("impc_df.csv")
impc_m <- impc_df %>% filter(sex=="male") 
impc_wt_m<-impc_m %>% filter(group =="control" & ee > (0.3*4.184))
impc_ko_m<-impc_m %>% filter(group =="experimental" & strain != "Fbxl19") 
ko_model <- lm(ee ~ total_mass_1 + room_temp, data = impc_wt_m)
ko_res <- impc_ko_m %>% add_residuals(model=ko_model)
ko_res_mean <- impc_ko_m %>% add_residuals(model=ko_model) %>% 
  group_by(strain, zygosity) %>% summarise_if(is.numeric, mean, na.rm=T)

#6A
ko_res_mean %>% ungroup() %>% 
  ggplot(.,aes(y = resid, x = total_mass_1)) +
  geom_hline(yintercept = 0, linetype="dotted", color="red") +
  xlim(17,42) +
  geom_point(pch=21, color="black", fill="#00B0F6") + theme(legend.position = "none") + 
  xlab("body mass of KO strain (g)") + ylab("residual of energy expenditure (kJ/hr)") +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "residual of energy expenditure (kcal/hr)"), limits = c(-0.2*4.184, 0.2*4.184))
  
#6B
ko_res_mean %>% filter(resid< 0.2*4.184 & resid > -0.2*4.184) %>% ggplot(aes(x=resid)) + 
  xlab("residual of energy expenditure (kJ/hr)") + ylab("number of KO strains") + 
  geom_vline(xintercept=sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_vline(xintercept=-sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_vline(xintercept=2*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_vline(xintercept=-2*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_vline(xintercept=3*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_vline(xintercept=-3*sd(ko_res_mean$resid, na.rm = T),linetype="dashed", color = "red") + theme(legend.position = "none") +
  geom_histogram(aes(y=..count..,fill=..count..),binwidth = 0.02*4.184)+ 
  xlim(c(-0.2*4.184,0.2*4.184)) +
  stat_function(
    fun = function(x, mean, sd, n, bw){
      dnorm(x = x, mean = mean, sd = sd) * n * bw
    }, color="red",
    args = c(mean = mean(ko_res_mean$resid, na.rm = T), sd = sd(ko_res_mean$resid, na.rm = T), n = sum(!is.na(ko_res_mean$resid)), bw = 0.02*4.184))


ko_res_zyg <- ko_res_mean %>% mutate(zyg = ifelse(zygosity=="homozygote","-/-","+/-")) %>% unite("str",c(strain, zyg), sep = " ", remove=F) %>% 
  filter((ee < 0.7*4.184) & (resid < 0.2*4.184) & (resid > -0.2*4.184))
phen_res <- ko_res_zyg %>% filter (strain %in% c("Ucp1","Pparg","Fgfr4", "Pcsk1","Cpe","Mrap2", "Gh","Acer1","Gdf15","Gatm") )
#6C
ko_res_zyg %>% 
  ggplot(., aes(x=total_mass_1, y=ee))+
  scale_fill_viridis(option="magma", direction = -1)+
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey", size = 2) +  
  geom_point(data=ko_res_zyg,size=1,color="darkgrey")+labs(fill="residual")+
  geom_point(data=phen_res,aes(fill=resid),pch=21,colour="black", size=3)+
  geom_label_repel(data=phen_res, color="black", size=3, aes(label=c(str)), box.padding = unit(1, 'lines'),force=10)+ xlim(15,42) +
  xlab("body mass (g)")+ylab("energy expenditure (kJ/hr)") + theme(legend.position = c(.78,.25)) +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))

#6D
ko_res_zyg %>% 
  ggplot(., aes(x=total_mass_1, y=resid))+
  scale_fill_viridis(option="magma", direction = -1)+
  geom_hline(yintercept=sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_hline(yintercept=-sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_hline(yintercept=2*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_hline(yintercept=-2*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_hline(yintercept=3*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_hline(yintercept=-3*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype="dotted", color = "dimgrey") +
  geom_point(data=ko_res_zyg,size=1,color="dimgrey")+
  geom_point(data=phen_res,aes(fill=resid),pch=21,colour="black", size=3)+
  geom_label_repel(data=phen_res, color="black", size=3, aes(label=str), box.padding = unit(1, 'lines'),force=10)+ xlim(15,42)+
  xlab("body mass (g)")+ylab("residual of energy expenditure (kJ/hr)") + theme(legend.position = "none") +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "residual of energy expenditure (kcal/hr)"))

GWAS_obesity_list<-c("Klf12","Pacs1","Chst8","Pald1","Tfap2b","Pax5","Pepd")
gwas_res <- ko_res_zyg %>% filter (strain %in% GWAS_obesity_list)
#6E
ko_res_zyg %>% 
  ggplot(., aes(x=total_mass_1, y=ee))+
  scale_fill_viridis(option="magma", direction = -1)+
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey", size = 2) +  
  geom_point(data=ko_res_zyg,size=1,color="darkgrey")+labs(fill="residual")+
  geom_point(data=gwas_res,aes(fill=resid),pch=21,colour="black", size=3)+
  geom_label_repel(data=gwas_res, color="black", size=3, aes(label=c(str)), box.padding = unit(1, 'lines'),force=10)+ xlim(15,42) +
  xlab("body mass (g)")+ylab("energy expenditure (kJ/hr)") + theme(legend.position = c(.78,.25)) +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))
#6F
ko_res_zyg %>% 
  ggplot(., aes(x=total_mass_1, y=resid))+
  scale_fill_viridis(option="magma", direction = -1)+
  geom_hline(yintercept=sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_hline(yintercept=-sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "green") +
  geom_hline(yintercept=2*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_hline(yintercept=-2*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "orange") +
  geom_hline(yintercept=3*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_hline(yintercept=-3*sd(ko_res_zyg$resid, na.rm = T),linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype="dotted", color = "dimgrey") +
  geom_point(data=ko_res_zyg,size=1,color="dimgrey")+
  geom_point(data=gwas_res,aes(fill=resid),pch=21,colour="black", size=3)+
  geom_label_repel(data=gwas_res, color="black", size=3, aes(label=str), box.padding = unit(1, 'lines'),force=10)+ xlim(15,42)+
  xlab("body mass (g)")+ylab("residual of energy expenditure (kJ/hr)") + theme(legend.position = "none") +
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "residual of energy expenditure (kcal/hr)"))
