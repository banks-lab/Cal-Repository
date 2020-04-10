library(tidyverse)
library(cowplot)
library(broom)
library(janitor)
library(modelr)

impc_df <- read.csv("impc_df.csv")
cdk <- read_csv("cdkal1.csv") %>% clean_names() %>% 
  dplyr::select(total_mass_1 = total_mass, ee = energy_expenditure , group, loc_act = locomotor_activity) %>%
  mutate(location = "Massachusetts, USA", sex = "male", room_temp = 23,zygosity="homozygote") %>%
  mutate(strain=ifelse(group=="WT",NA_character_,"Cdkal1"),group=ifelse(group=="WT","control","experimental"), ee = ee * 4.184)
impc_m <- impc_df %>% bind_rows(cdk) %>% filter(sex=="male") %>%
  mutate(location = factor(location, levels = c("Cambridgeshire, UK","Germany","Oxfordshire, UK","California, USA","France","China","Japan","Canada","Korea", "Texas, USA", "Massachusetts, USA")))
impc_wt_m<-impc_m %>% filter(group =="control" & ee >0.3*4.184) 
impc_ko_m<-impc_m %>% filter(group =="experimental" & strain != "Fbxl19") 
ko_model <- lm(ee ~ total_mass_1 + room_temp, data = impc_wt_m)
ko_res <- impc_ko_m %>% add_residuals(model=ko_model)
ko_res_mean <- impc_ko_m %>% add_residuals(model=ko_model) %>% 
  group_by(location, strain, zygosity) %>% summarise_if(is.numeric, mean, na.rm=T)
ko_res_mean %>% ungroup() %>% count(strain, zygosity, sort = T) %>% top_n(7)
ko_sites<-ko_res %>% filter ((strain == "Dnase1l2" & zygosity=="homozygote") | (strain == "Ap4e1" & zygosity=="homozygote") | 
                                    (strain == "Dbn1" & zygosity=="heterozygote") | (strain == "Nxn" |strain == "Prkab1"|strain == "Cdkal1" ) | 
                                    (strain=="Rnf10" & zygosity=="homozygote" & location != "Canada") )
ko_sites_mean<-ko_res_mean %>% filter ((strain == "Dnase1l2" & zygosity=="homozygote") | (strain == "Ap4e1" & zygosity=="homozygote")| 
                                              (strain == "Dbn1" & zygosity=="heterozygote") | (strain == "Nxn" |strain == "Prkab1"|strain == "Cdkal1" ) | 
                                              (strain=="Rnf10" & zygosity=="homozygote" & location != "Canada"))
all_sites<-impc_m %>% filter ((strain == "Dnase1l2" & zygosity=="homozygote") | (strain == "Ap4e1" & zygosity=="homozygote") | 
                            (strain == "Dbn1" & zygosity=="heterozygote") | (strain == "Nxn" |strain == "Prkab1"|strain == "Cdkal1" ) | 
                            (strain=="Rnf10" & zygosity=="homozygote" & location != "Canada") | is.na(strain))
impc_colors<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928')
#5A
ko_sites %>% 
  mutate(strain = factor(strain, labels = c("Ap4e1 -/-","Cdkal1 -/-","Dbn1 +/-","Dnase1l2 -/-","Nxn +/-","Prkab1 -/-","Rnf10 -/-"))) %>%
  ggplot(., aes(y=ee, x=total_mass_1)) + 
  geom_point(pch=21, color="black",size=2,aes(fill=location))+
  geom_smooth(method = "lm",se=F, size=1, aes(color=location))+ 
  geom_smooth(method = "lm",se=F, size=2, color="black")+ 
  scale_fill_manual(values=impc_colors) +
  scale_color_manual(values=impc_colors) +
  ylab("energy expenditure (kJ/hr)") + xlab("body mass (g)")+ 
  xlim(c(17, 35))+ facet_wrap(~strain, ncol=1)+
  theme(legend.position = "none") + 
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "energy expenditure (kcal/hr)"))
#5B
ko_sites_mean %>% ungroup() %>% 
  mutate(strain = factor(strain, labels = c("Ap4e1 -/-","Cdkal1 -/-","Dbn1 +/-","Dnase1l2 -/-","Nxn +/-","Prkab1 -/-","Rnf10 -/-"))) %>%
  ggplot(., aes(y=resid, x=total_mass_1, color=location)) + 
  geom_hline(yintercept=0, linetype=3, color="red")+
  geom_point(pch=21, color = "black",size=2, aes(fill=location))+
  scale_fill_manual(values=impc_colors) +
  scale_color_manual(values=impc_colors) +
  ylab("residual of energy expenditure (kJ/hr)") +  xlab("body mass (g)")+ 
  facet_wrap(~strain, ncol=1)+xlim(c(17, 35))+
  theme(legend.position = "none") + 
  scale_y_continuous(sec.axis = sec_axis(~./4.184, name = "residual of energy expenditure (kcal/hr)"))


act_sites <- all_sites %>% filter(!(location %in% c("Germany","Texas, USA","China", "California, USA")))
no_act_sites <- all_sites %>% filter(location %in% c("Germany","Texas, USA","China", "California, USA"))

strain_model <- function(ko,sites) {
  a <- act_sites %>% filter((is.na(strain) | strain == ko) & location %in% sites) %>% 
    group_by(location) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(ee ~ total_mass_1 + room_temp + loc_act + group, data = .x)),
      tidied = map(fit, tidy),
      glanced = map(fit, glance)
    ) %>% 
    unnest(tidied)
  b <- no_act_sites %>% filter((is.na(strain) | strain == ko) & location %in% sites) %>% 
    group_by(location) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(ee ~ total_mass_1 + room_temp + group, data = .x)),
      tidied = map(fit, tidy),
      glanced = map(fit, glance)
    ) %>% 
    unnest(tidied)
  bind_rows(a,b) %>% filter(term == "groupexperimental") %>% mutate(strain = ko) %>% arrange(p.value)
}

strain_model("Ap4e1",c("Cambridgeshire, UK","Germany","Oxfordshire, UK","California, USA","France","Japan","Canada","Korea", "Texas, USA"))
strain_model("Cdkal1",c("Massachusetts, USA", "Germany","Cambridgeshire, UK"))
strain_model("Dbn1",c("Cambridgeshire, UK","Germany","California, USA","France","Japan","Canada","Korea"))
all_sites %>% filter((is.na(strain) | strain =="Dnase1l2") & location == "Korea") %>% lm(ee~total_mass_1 + room_temp + group, data = .) %>%
  tidy() %>% filter(term == "groupexperimental") %>% mutate(strain = "Dnase1l2", location = "Korea") %>% 
  bind_rows(strain_model("Dnase1l2",c("Cambridgeshire, UK","Germany","California, USA","France","Japan","Canada","Texas, USA")))
strain_model("Nxn",c("Cambridgeshire, UK","Germany","California, USA","France","China","Japan","Canada","Texas, USA"))
strain_model("Rnf10",c("Korea","France","Germany","Cambridgeshire, UK"))
strain_model("Prkab1",c("Cambridgeshire, UK","Oxfordshire, UK","California, USA","France","Japan","Canada","Korea"))
