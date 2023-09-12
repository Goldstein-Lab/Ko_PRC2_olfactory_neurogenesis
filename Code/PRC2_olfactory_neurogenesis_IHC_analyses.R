library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(rstatix)
library(ggrepel)

##### import data #####
setwd("/Users/story/Box/Goldstein Lab/Tiffany/R/")
EZH2_BrdU_Normal_OBX<- read_excel("EZH2_BrdU_Normal_OBX_quantification_RC_Sep2018.xlsx")
Mx14d_EZH2 <- read_excel("Mx_14d_EZH2_quantification.xlsx")
OBx_EED <- read_excel("OBx_EED_Jarid.xlsx", sheet = "EED")
OBx_Jarid <- read_excel("OBx_EED_Jarid.xlsx", sheet = "Jarid")
EED <- read_excel("EED-KO-het-quantification-EED-H3K27Me3-CBX8.xlsx", sheet = "EED")
CBX8 <- read_excel("EED-KO-het-quantification-EED-H3K27Me3-CBX8.xlsx", sheet = "CBX8")
EED_KO_Het_Ki67 <- read_excel("EED KO-Het 9.17.22 Tiffanay.xlsx", 
                              sheet = "Ki67-tidy")
EED_KO_Het_Caspase_tom <- read_excel("EEDKOHET_Casp3_Tom_counts.xlsx")
Sox9_EEDKOHET <- read_excel("C:/Users/story/Box/Goldstein Lab/Tiffany/Images/SOX9-EEDKOHET/blinded counts Sox9 EEDKOHET.xlsx", 
                            sheet = "final decoded counts")
shape <- read_excel("EEDKOHET-cellshape/EEDKOHET_cellshape.xlsx")

##### Mx14d EZH2 #####
#Data Management
normal_EZH2 <- EZH2_BrdU_Normal_OBX %>% 
  filter(Condition == "Normal")
sum_normal_EZH2 <-normal_EZH2 %>% 
  group_by(Unique_id) %>% 
  summarise(EZH2_sum = sum(EZH2), length_sum = sum(Length)) %>% 
  mutate(EZH2_per_mm = EZH2_sum/length_sum*1000) %>% 
  mutate(condition = "Normal") %>% 
  mutate(Animal = Unique_id) %>% 
  select(-Unique_id)

sum_Mx14d_EZH2 <- Mx14d_EZH2 %>% 
  group_by(Animal) %>% 
  summarise(EZH2_sum = sum(EZH2_cells), length_sum = sum(Length)) %>% 
  mutate(EZH2_per_mm = EZH2_sum/length_sum*1000) %>% 
  mutate(condition = "14dpMx")
sum_Mx14d_EZH2$Animal <- as.character(sum_Mx14d_EZH2$Animal)

sum_Mx14d_normal_EZH2 <- bind_rows(sum_normal_EZH2, sum_Mx14d_EZH2)

mean_Mx14d_normal_EZH2 <- sum_Mx14d_normal_EZH2 %>% 
  group_by(condition) %>% 
  summarise(n=n(),
            mean_EZH2_per_mm = mean(EZH2_per_mm),
            sd = sd(EZH2_per_mm)) %>% 
  mutate(se=sd/sqrt(n))

MxEZH2_sum_mean <- left_join(mean_Mx14d_normal_EZH2, sum_Mx14d_normal_EZH2, by="condition")

MxEZH2_sum_mean$condition <- factor(MxEZH2_sum_mean$condition, 
                                    levels = c("Normal","14dpMx"))

#Normality Test
shapiro.test(sum_normal_EZH2$EZH2_per_mm)
qqnorm(sum_normal_EZH2$EZH2_per_mm, pch = 1, frame = FALSE)
qqline(sum_normal_EZH2$EZH2_per_mm, col = "steelblue", lwd = 2)
shapiro.test(sum_Mx14d_EZH2$EZH2_per_mm)
qqnorm(sum_Mx14d_EZH2$EZH2_per_mm, pch = 1, frame = FALSE)
qqline(sum_Mx14d_EZH2$EZH2_per_mm, col = "steelblue", lwd = 2)

#T-test
EZH2.stat.test <- sum_Mx14d_normal_EZH2 %>% 
  t_test(EZH2_per_mm ~ condition, var.equal = TRUE) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p < "0.05", paste("*"), paste("n.s."))))

#Data Visualization
MxEZH2plot <- MxEZH2_sum_mean %>% ggplot(aes(condition, EZH2_per_mm,
                                             group=condition, color = condition))+
  geom_bar(aes(fill = condition), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 1.5, 
                                             jitter.width = 0.2))+
  geom_errorbar(aes(x=condition, ymin=mean_EZH2_per_mm-se, ymax=mean_EZH2_per_mm+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "condition", 
                     labels = c("Normal", "14dpMx"))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,500)+
  xlab("")+
  ylab("EZH2+ cells
       (per mm epithelium)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=18, colour = "black"), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18, colour = "black"))+
  scale_x_discrete(labels = c("Normal", "14dpMx"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE), 
                     label = "p.signif", label.x.npc = "center", label.y=450,
                     size=10)
MxEZH2plot

##### OBx EED normalized to epithelial length#####
#Data Management 
sum_OBx_EED <-OBx_EED %>% 
  group_by(Unique_id) %>% 
  summarise(EED_sum = sum(EED), length_sum = sum(Length), 
            across(c("Animal", "Condition"))) %>% 
  mutate(EED_per_mm = EED_sum/length_sum*1000) 

sum_OBx_EED<-sum_OBx_EED[!duplicated(sum_OBx_EED), ]

mean_OBx_EED <- sum_OBx_EED %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_EED_per_mm = mean(EED_per_mm),
            sd = sd(EED_per_mm)) %>% 
  mutate(se=sd/sqrt(n))

OBx_EED_sum_mean <- left_join(mean_OBx_EED, sum_OBx_EED, by="Condition")

sum_OBx_EED_OBx <- sum_OBx_EED %>% 
  filter(Condition =="OBx")
sum_OBx_EED_norm <- sum_OBx_EED %>% 
  filter(Condition =="Normal")

#Normality test
shapiro.test(sum_OBx_EED_norm$EED_per_mm)
qqnorm(sum_OBx_EED_norm$EED_per_mm, pch = 1, frame = FALSE)
qqline(sum_OBx_EED_norm$EED_per_mm, col = "steelblue", lwd = 2)
shapiro.test(sum_OBx_EED_OBx$EED_per_mm) # p=0.01
qqnorm(sum_OBx_EED_OBx$EED_per_mm, pch = 1, frame = FALSE)
qqline(sum_OBx_EED_OBx$EED_per_mm, col = "steelblue", lwd = 2)

#T-test
EED.stat.test <- t.test(sum_OBx_EED_norm$EED_per_mm, sum_OBx_EED_OBx$EED_per_mm,
                        paired = TRUE)

#Data Visualization
OBx_EEDplot <- OBx_EED_sum_mean %>% ggplot(aes(Condition, EED_per_mm,
                                               group=Condition, color = Condition))+
  geom_bar(aes(fill = Condition), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = Condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 1.5, 
                                             jitter.width = 0.69))+
  geom_errorbar(aes(x=Condition, ymin=mean_EED_per_mm-se, ymax=mean_EED_per_mm+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "Condition", 
                     labels = c("Normal", "OBx"))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,100)+
  xlab("")+
  ylab("EED+ cells
       (per mm epithelium)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=16, colour = "black"), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18, colour = "black"))+
  scale_x_discrete(labels = c("Normal", "OBx"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE, paired = TRUE), 
                     label = "p.signif", label.x.npc = "center", label.y=90,
                     size=10)
OBx_EEDplot

##### OBx JARID2 normalized to epithelial length #####
#Data Management
sum_OBx_Jarid <-OBx_Jarid %>% 
  group_by(Unique_id) %>% 
  summarise(Jarid_sum = sum(JARID), length_sum = sum(Length), 
            across(c("Animal", "Condition"))) %>% 
  mutate(Jarid_per_mm = Jarid_sum/length_sum*1000) 

sum_OBx_Jarid<-sum_OBx_Jarid[!duplicated(sum_OBx_Jarid), ]

mean_OBx_Jarid <- sum_OBx_Jarid %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_Jarid_per_mm = mean(Jarid_per_mm),
            sd = sd(Jarid_per_mm)) %>% 
  mutate(se=sd/sqrt(n))

OBx_Jarid_sum_mean <- left_join(mean_OBx_Jarid, sum_OBx_Jarid, by="Condition")


sum_OBx_Jarid_OBx <- sum_OBx_Jarid %>% 
  filter(Condition =="OBx") %>% 
  select(Condition, Animal, Jarid_per_mm)
sum_OBx_Jarid_norm <- sum_OBx_Jarid %>% 
  filter(Condition =="Normal") %>% 
  select(Condition, Animal, Jarid_per_mm)

#Normality Test
shapiro.test(sum_OBx_Jarid_norm$Jarid_per_mm)
qqnorm(sum_OBx_Jarid_norm$Jarid_per_mm, pch = 1, frame = FALSE)
qqline(sum_OBx_Jarid_norm$Jarid_per_mm, col = "steelblue", lwd = 2)
shapiro.test(sum_OBx_Jarid_OBx$Jarid_per_mm) 
qqnorm(sum_OBx_Jarid_OBx$Jarid_per_mm, pch = 1, frame = FALSE)
qqline(sum_OBx_Jarid_OBx$Jarid_per_mm, col = "steelblue", lwd = 2)

#T-test
Jarid.stat.test <- t.test(sum_OBx_Jarid_norm$Jarid_per_mm, sum_OBx_Jarid_OBx$Jarid_per_mm,
                          paired = TRUE)
Jarid.stat.test

#Data Visualization
OBx_Jaridplot <- OBx_Jarid_sum_mean %>% ggplot(aes(Condition, Jarid_per_mm,
                                                   group=Condition, color = Condition))+
  geom_bar(aes(fill = Condition), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = Condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 1.5, 
                                             jitter.width = 0.5))+
  geom_errorbar(aes(x=Condition, ymin=mean_Jarid_per_mm-se, ymax=mean_Jarid_per_mm+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "Condition", 
                     labels = c("Normal", "OBx"))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,320)+
  xlab("")+
  ylab("JARID2+ cells
       (per mm epithelium)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=16, colour = "black"), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18, colour = "black"))+
  scale_x_discrete(labels = c("Normal", "OBx"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE, paired = TRUE), 
                     label = "p.signif", label.x.npc = "center", label.y=300,
                     size=10)
OBx_Jaridplot

##### OBx EZH2 normalized to epithelial length ####
#Data Management
sum_OBx_Ezh2_Brdu <-EZH2_BrdU_Normal_OBX %>% 
  group_by(Unique_id) %>% 
  summarise(EZH2_sum = sum(EZH2), 
            Brdu_sum= sum(Overlap),
            length_sum = sum(Length), 
            across(c("Animal", "Condition"))) %>% 
  mutate(EZH2_per_mm = EZH2_sum/length_sum*1000) %>% 
  mutate(Brdu_per_mm = Brdu_sum/length_sum*1000)

sum_OBx_Ezh2_Brdu<-sum_OBx_Ezh2_Brdu[!duplicated(sum_OBx_Ezh2_Brdu), ]

mean_OBx_Ezh2_Brdu<- sum_OBx_Ezh2_Brdu %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_EZH2_per_mm = mean(EZH2_per_mm),
            mean_Brdu_per_mm = mean(Brdu_per_mm),
            sd_EZH2 = sd(EZH2_per_mm),
            sd_Brdu = sd(Brdu_per_mm)) %>% 
  mutate(se_EZH2=sd_EZH2/sqrt(n),
         se_Brdu=sd_Brdu/sqrt(n))

OBx_Ezh2_Brdu_sum_mean <- left_join(mean_OBx_Ezh2_Brdu, sum_OBx_Ezh2_Brdu, by="Condition")


sum_OBx_Ezh2_Brdu_OBx <- sum_OBx_Ezh2_Brdu %>% 
  filter(Condition =="OBX") %>% 
  select(Condition, Animal, EZH2_per_mm, Brdu_per_mm)
sum_OBx_Ezh2_Brdu_norm <- sum_OBx_Ezh2_Brdu %>% 
  filter(Condition =="Normal") %>% 
  select(Condition, Animal, EZH2_per_mm, Brdu_per_mm)

#Normality Test
shapiro.test(sum_OBx_Ezh2_Brdu_norm$EZH2_per_mm)
qqnorm(sum_OBx_Ezh2_Brdu_norm$EZH2_per_mm, pch = 1, frame = FALSE)
qqline(sum_OBx_Ezh2_Brdu_norm$EZH2_per_mm, col = "steelblue", lwd = 2)
shapiro.test(sum_OBx_Ezh2_Brdu_OBx$EZH2_per_mm) 
qqnorm(sum_OBx_Ezh2_Brdu_OBx$EZH2_per_mm, pch = 1, frame = FALSE)
qqline(sum_OBx_Ezh2_Brdu_OBx$EZH2_per_mm, col = "steelblue", lwd = 2)

#T-tests
Ezh2.stat.test <- t.test(sum_OBx_Ezh2_Brdu_norm$EZH2_per_mm, sum_OBx_Ezh2_Brdu_OBx$EZH2_per_mm,
                         paired = TRUE)
EZH2.stat.test

Brdu.stat.test <- t.test(sum_OBx_Ezh2_Brdu_norm$Brdu_per_mm, sum_OBx_Ezh2_Brdu_OBx$Brdu_per_mm,
                         paired = TRUE)
Brdu.stat.test

mean_OBx_Ezh2<- sum_OBx_Ezh2_Brdu %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_pos_per_mm = mean(EZH2_per_mm),
            sd = sd(EZH2_per_mm)) %>% 
  mutate(se=sd/sqrt(n))

mean_OBx_Brdu<- sum_OBx_Ezh2_Brdu %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_pos_per_mm = mean(Brdu_per_mm),
            sd = sd(Brdu_per_mm)) %>% 
  mutate(se=sd/sqrt(n))

sum_OBx_Ezh2 <- sum_OBx_Ezh2_Brdu %>% 
  mutate(pos_per_mm = EZH2_per_mm) %>% 
  mutate(pos_sum = EZH2_sum) %>% 
  select(-c("Brdu_sum", "Brdu_per_mm", "EZH2_per_mm", "EZH2_sum")) %>% 
  mutate(count_type = "Ezh2")

sum_OBx_Brdu <- sum_OBx_Ezh2_Brdu %>% 
  mutate(pos_per_mm = Brdu_per_mm) %>% 
  mutate(pos_sum = Brdu_sum) %>% 
  select(-c("Brdu_sum", "Brdu_per_mm", "EZH2_per_mm", "EZH2_sum")) %>% 
  mutate(count_type = "Brdu")

mean_sum_OBx_Ezh2 <- left_join(mean_OBx_Ezh2, sum_OBx_Ezh2, by = "Condition")
mean_sum_OBx_Brdu <- left_join(mean_OBx_Brdu, sum_OBx_Brdu, by = "Condition")

plotdf_OBx_Ezh2_Brdu <- bind_rows(mean_sum_OBx_Brdu, mean_sum_OBx_Ezh2)
plotdf_OBx_Ezh2_Brdu$count_type <- factor(plotdf_OBx_Ezh2_Brdu$count_type,
                                          levels = c("Ezh2","Brdu"))
#Data Visualization
Ezh2_brdu.stat.test <- plotdf_OBx_Ezh2_Brdu %>% 
  group_by(count_type) %>% 
  t_test(pos_per_mm ~ Condition, var.equal = FALSE) %>% 
  add_xy_position(x = "count_type", dodge = 0.8) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p<"0.05", paste("*"),paste("n.s."))))

OBx_EZH2_Brdu_plot <- plotdf_OBx_Ezh2_Brdu %>% ggplot(aes(count_type, pos_per_mm, 
                                                          group= Condition,
                                                          color = Condition))+
  stat_summary(aes(fill = Condition),geom = "col", fun = mean, 
               position = position_dodge())+
  geom_errorbar(aes(x=count_type,group=Condition, 
                    ymin=mean_pos_per_mm-se, ymax=mean_pos_per_mm+se), 
                position = position_dodge(width = 0.9), width=0.4,
                colour="black",  size=0.5)+
  geom_point(aes(fill=Condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.97, 
                                             jitter.width = 0.36))+
  scale_color_manual(values = c("#000000","#000000"), guide = "none")+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"),name = "Condition", 
                    labels = c("EZH2+","EZH2+/BrdU+"))+
  ylim(0,450)+
  ylab("Cells per mm OE")+
  xlab ("")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.text = element_text(size=20, colour = "black"),
        axis.title = element_text(size = 24, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "none")+
  stat_pvalue_manual(Ezh2_brdu.stat.test, label = "{p.signif}", size = 6)

OBx_EZH2_Brdu_plot

##### Mx14d EZH2 normalized to DAPI#####
normal_EZH2 <- EZH2_BrdU_Normal_OBX %>% 
  filter(Condition == "Normal")
sum_normal_EZH2 <-normal_EZH2 %>% 
  group_by(Unique_id) %>% 
  summarise(EZH2_sum = sum(EZH2), dapi_sum = sum(DAPI)) %>% 
  mutate(EZH2_per_dapi = EZH2_sum/dapi_sum) %>% 
  mutate(condition = "Normal") %>% 
  mutate(Animal = Unique_id) %>% 
  select(-Unique_id)

sum_Mx14d_EZH2 <- Mx14d_EZH2 %>% 
  group_by(Animal) %>% 
  summarise(EZH2_sum = sum(EZH2_cells), dapi_sum = sum(DAPI)) %>% 
  mutate(EZH2_per_dapi = EZH2_sum/dapi_sum) %>% 
  mutate(condition = "14dpMx")
sum_Mx14d_EZH2$Animal <- as.character(sum_Mx14d_EZH2$Animal)

sum_Mx14d_normal_EZH2 <- bind_rows(sum_normal_EZH2, sum_Mx14d_EZH2)

mean_Mx14d_normal_EZH2 <- sum_Mx14d_normal_EZH2 %>% 
  group_by(condition) %>% 
  summarise(n=n(),
            mean_EZH2_per_dapi = mean(EZH2_per_dapi),
            sd = sd(EZH2_per_dapi)) %>% 
  mutate(se=sd/sqrt(n))

shapiro.test(sum_normal_EZH2$EZH2_per_dapi)
qqnorm(sum_normal_EZH2$EZH2_per_dapi, pch = 1, frame = FALSE)
qqline(sum_normal_EZH2$EZH2_per_dapi, col = "steelblue", lwd = 2)
shapiro.test(sum_Mx14d_EZH2$EZH2_per_dapi)
qqnorm(sum_Mx14d_EZH2$EZH2_per_dapi, pch = 1, frame = FALSE)
qqline(sum_Mx14d_EZH2$EZH2_per_dapi, col = "steelblue", lwd = 2)

EZH2.stat.test <- sum_Mx14d_normal_EZH2 %>% 
  t_test(EZH2_per_dapi ~ condition, var.equal = TRUE) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p < "0.05", paste("*"), paste("n.s."))))
MxEZH2_sum_mean <- left_join(mean_Mx14d_normal_EZH2, sum_Mx14d_normal_EZH2, by="condition")

MxEZH2_sum_mean$condition <- factor(MxEZH2_sum_mean$condition, 
                                    levels = c("Normal","14dpMx"))


shapiro.test(sum_normal_EZH2$EZH2_per_dapi)
qqnorm(sum_normal_EZH2$EZH2_per_dapi, pch = 1, frame = FALSE)
qqline(sum_normal_EZH2$EZH2_per_dapi, col = "steelblue", lwd = 2)
shapiro.test(sum_Mx14d_EZH2$EZH2_per_dapi)
qqnorm(sum_Mx14d_EZH2$EZH2_per_dapi, pch = 1, frame = FALSE)
qqline(sum_Mx14d_EZH2$EZH2_per_dapi, col = "steelblue", lwd = 2)

EZH2.stat.test <- sum_Mx14d_normal_EZH2 %>% 
  t_test(EZH2_per_dapi ~ condition, var.equal = TRUE) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p < "0.05", paste("*"), paste("n.s."))))

MxEZH2plot <- MxEZH2_sum_mean %>% ggplot(aes(condition, EZH2_per_dapi,
                                             group=condition, color = condition))+
  geom_bar(aes(fill = condition), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 1.5, 
                                             jitter.width = 0.5))+
  geom_errorbar(aes(x=condition, ymin=mean_EZH2_per_dapi-se, ymax=mean_EZH2_per_dapi+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "condition", 
                     labels = c("Normal", "14dpMx"))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,0.3)+
  xlab("")+
  ylab("EZH2+ cells/DAPI")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=18, colour = "black"), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18, colour = "black"))+
  scale_x_discrete(labels = c("Normal", "14dpMx"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE), 
                     label = "p.signif", label.x.npc = "center", label.y=0.27,
                     size=10)
MxEZH2plot

##### OBx EED normalized to DAPI #####
sum_OBx_EED <-OBx_EED %>% 
  group_by(Unique_id) %>% 
  summarise(EED_sum = sum(EED), dapi_sum = sum(DAPI), 
            across(c("Animal", "Condition"))) %>% 
  mutate(EED_per_dapi = EED_sum/dapi_sum) 

sum_OBx_EED<-sum_OBx_EED[!duplicated(sum_OBx_EED), ]

mean_OBx_EED <- sum_OBx_EED %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_EED_per_dapi = mean(EED_per_dapi),
            sd = sd(EED_per_dapi)) %>% 
  mutate(se=sd/sqrt(n))

OBx_EED_sum_mean <- left_join(mean_OBx_EED, sum_OBx_EED, by="Condition")

sum_OBx_EED_OBx <- sum_OBx_EED %>% 
  filter(Condition =="OBx")
sum_OBx_EED_norm <- sum_OBx_EED %>% 
  filter(Condition =="Normal")

shapiro.test(sum_OBx_EED_norm$EED_per_dapi)
qqnorm(sum_OBx_EED_norm$EED_per_dapi, pch = 1, frame = FALSE)
qqline(sum_OBx_EED_norm$EED_per_dapi, col = "steelblue", lwd = 2)
shapiro.test(sum_OBx_EED_OBx$EED_per_dapi) # p=0.01
qqnorm(sum_OBx_EED_OBx$EED_per_dapi, pch = 1, frame = FALSE)
qqline(sum_OBx_EED_OBx$EED_per_dapi, col = "steelblue", lwd = 2)

EED.stat.test <- t.test(sum_OBx_EED_norm$EED_per_dapi, sum_OBx_EED_OBx$EED_per_dapi,
                        paired = TRUE)


OBx_EEDplot <- OBx_EED_sum_mean %>% ggplot(aes(Condition, EED_per_dapi,
                                               group=Condition, color = Condition))+
  geom_bar(aes(fill = Condition), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = Condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 1.3, 
                                             jitter.width = 0.64))+
  geom_errorbar(aes(x=Condition, ymin=mean_EED_per_dapi-se, ymax=mean_EED_per_dapi+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "Condition", 
                     labels = c("Normal", "OBx"))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,0.07)+
  xlab("")+
  ylab("EED+ cells/DAPI")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=16, colour = "black"), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18, colour = "black"))+
  scale_x_discrete(labels = c("Normal", "OBx"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE, paired = TRUE), 
                     label = "p.signif", label.x.npc = "center", label.y=0.06,
                     size=10)
OBx_EEDplot

##### OBx JARID2 normalized to DAPI #####
sum_OBx_Jarid <-OBx_Jarid %>% 
  group_by(Unique_id) %>% 
  summarise(Jarid_sum = sum(JARID), dapi_sum = sum(DAPI), 
            across(c("Animal", "Condition"))) %>% 
  mutate(Jarid_per_dapi = Jarid_sum/dapi_sum) 

sum_OBx_Jarid<-sum_OBx_Jarid[!duplicated(sum_OBx_Jarid), ]

mean_OBx_Jarid <- sum_OBx_Jarid %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_Jarid_per_dapi = mean(Jarid_per_dapi),
            sd = sd(Jarid_per_dapi)) %>% 
  mutate(se=sd/sqrt(n))

OBx_Jarid_sum_mean <- left_join(mean_OBx_Jarid, sum_OBx_Jarid, by="Condition")


sum_OBx_Jarid_OBx <- sum_OBx_Jarid %>% 
  filter(Condition =="OBx") %>% 
  select(Condition, Animal, Jarid_per_dapi)

sum_OBx_Jarid_norm <- sum_OBx_Jarid %>% 
  filter(Condition =="Normal") %>% 
  select(Condition, Animal, Jarid_per_dapi)

shapiro.test(sum_OBx_Jarid_norm$Jarid_per_dapi)
qqnorm(sum_OBx_Jarid_norm$Jarid_per_dapi, pch = 1, frame = FALSE)
qqline(sum_OBx_Jarid_norm$Jarid_per_dapi, col = "steelblue", lwd = 2)
shapiro.test(sum_OBx_Jarid_OBx$Jarid_per_dapi) 
qqnorm(sum_OBx_Jarid_OBx$Jarid_per_dapi, pch = 1, frame = FALSE)
qqline(sum_OBx_Jarid_OBx$Jarid_per_dapi, col = "steelblue", lwd = 2)

Jarid.stat.test <- t.test(sum_OBx_Jarid_norm$Jarid_per_dapi, sum_OBx_Jarid_OBx$Jarid_per_dapi,
                          paired = TRUE)
Jarid.stat.test

OBx_Jaridplot <- OBx_Jarid_sum_mean %>% ggplot(aes(Condition, Jarid_per_dapi,
                                                   group=Condition, color = Condition))+
  geom_bar(aes(fill = Condition), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = Condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 1.5, 
                                             jitter.width = 0.5))+
  geom_errorbar(aes(x=Condition, ymin=mean_Jarid_per_dapi-se, ymax=mean_Jarid_per_dapi+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "Condition", 
                     labels = c("Normal", "OBx"))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,0.2)+
  xlab("")+
  ylab("JARID2+ cells/DAPI")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=16, colour = "black"), 
        axis.text.y = element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18, colour = "black"))+
  scale_x_discrete(labels = c("Normal", "OBx"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE, paired = TRUE), 
                     label = "p.signif", label.x.npc = "center", label.y=0.18,
                     size=10)
OBx_Jaridplot

##### OBx EZH2 normalized to DAPI ####
sum_OBx_Ezh2_Brdu <-EZH2_BrdU_Normal_OBX %>% 
  group_by(Unique_id) %>% 
  summarise(EZH2_sum = sum(EZH2), 
            Brdu_sum= sum(Overlap),
            dapi_sum = sum(DAPI), 
            across(c("Animal", "Condition"))) %>% 
  mutate(EZH2_per_dapi = EZH2_sum/dapi_sum) %>% 
  mutate(Brdu_per_dapi = Brdu_sum/dapi_sum)

sum_OBx_Ezh2_Brdu<-sum_OBx_Ezh2_Brdu[!duplicated(sum_OBx_Ezh2_Brdu), ]

mean_OBx_Ezh2_Brdu<- sum_OBx_Ezh2_Brdu %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_EZH2_per_dapi = mean(EZH2_per_dapi),
            mean_Brdu_per_dapi = mean(Brdu_per_dapi),
            sd_EZH2 = sd(EZH2_per_dapi),
            sd_Brdu = sd(Brdu_per_dapi)) %>% 
  mutate(se_EZH2=sd_EZH2/sqrt(n),
         se_Brdu=sd_Brdu/sqrt(n))

OBx_Ezh2_Brdu_sum_mean <- left_join(mean_OBx_Ezh2_Brdu, sum_OBx_Ezh2_Brdu, by="Condition")


sum_OBx_Ezh2_Brdu_OBx <- sum_OBx_Ezh2_Brdu %>% 
  filter(Condition =="OBX") %>% 
  select(Condition, Animal, EZH2_per_dapi, Brdu_per_dapi)
sum_OBx_Ezh2_Brdu_norm <- sum_OBx_Ezh2_Brdu %>% 
  filter(Condition =="Normal") %>% 
  select(Condition, Animal, EZH2_per_dapi, Brdu_per_dapi)

shapiro.test(sum_OBx_Ezh2_Brdu_norm$EZH2_per_dapi)
qqnorm(sum_OBx_Ezh2_Brdu_norm$EZH2_per_dapi, pch = 1, frame = FALSE)
qqline(sum_OBx_Ezh2_Brdu_norm$EZH2_per_dapi, col = "steelblue", lwd = 2)
shapiro.test(sum_OBx_Ezh2_Brdu_OBx$EZH2_per_dapi) 
qqnorm(sum_OBx_Ezh2_Brdu_OBx$EZH2_per_dapi, pch = 1, frame = FALSE)
qqline(sum_OBx_Ezh2_Brdu_OBx$EZH2_per_dapi, col = "steelblue", lwd = 2)



mean_OBx_Ezh2<- sum_OBx_Ezh2_Brdu %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_pos_per_dapi = mean(EZH2_per_dapi),
            sd = sd(EZH2_per_dapi)) %>% 
  mutate(se=sd/sqrt(n))

mean_OBx_Brdu<- sum_OBx_Ezh2_Brdu %>% 
  group_by(Condition) %>% 
  summarise(n=n(),
            mean_pos_per_dapi = mean(Brdu_per_dapi),
            sd = sd(Brdu_per_dapi)) %>% 
  mutate(se=sd/sqrt(n))

sum_OBx_Ezh2 <- sum_OBx_Ezh2_Brdu %>% 
  mutate(pos_per_dapi = EZH2_per_dapi) %>% 
  mutate(pos_sum = EZH2_sum) %>% 
  select(-c("Brdu_sum", "Brdu_per_dapi", "EZH2_per_dapi", "EZH2_sum")) %>% 
  mutate(count_type = "Ezh2")

sum_OBx_Brdu <- sum_OBx_Ezh2_Brdu %>% 
  mutate(pos_per_dapi = Brdu_per_dapi) %>% 
  mutate(pos_sum = Brdu_sum) %>% 
  select(-c("Brdu_sum", "Brdu_per_dapi", "EZH2_per_dapi", "EZH2_sum")) %>% 
  mutate(count_type = "Brdu")

mean_sum_OBx_Ezh2 <- left_join(mean_OBx_Ezh2, sum_OBx_Ezh2, by = "Condition")
mean_sum_OBx_Brdu <- left_join(mean_OBx_Brdu, sum_OBx_Brdu, by = "Condition")

plotdf_OBx_Ezh2_Brdu <- bind_rows(mean_sum_OBx_Brdu, mean_sum_OBx_Ezh2)
plotdf_OBx_Ezh2_Brdu$count_type <- factor(plotdf_OBx_Ezh2_Brdu$count_type,
                                          levels = c("Ezh2","Brdu"))

Ezh2_brdu.stat.test <- plotdf_OBx_Ezh2_Brdu %>% 
  group_by(count_type) %>% 
  t_test(pos_per_dapi ~ Condition, paired = TRUE) %>% 
  add_xy_position(x = "count_type", dodge = 0.8) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p<"0.05", paste("*"),paste("n.s."))))

OBx_EZH2_Brdu_plot <- plotdf_OBx_Ezh2_Brdu %>% ggplot(aes(count_type, pos_per_dapi, 
                                                          group= Condition,
                                                          color = Condition))+
  stat_summary(aes(fill = Condition),geom = "col", fun = mean, 
               position = position_dodge())+
  geom_errorbar(aes(x=count_type,group=Condition, 
                    ymin=mean_pos_per_dapi-se, ymax=mean_pos_per_dapi+se), 
                position = position_dodge(width = 0.9), width=0.4,
                colour="black",  size=0.5)+
  geom_point(aes(fill=Condition),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.80, 
                                             jitter.width = 0.45))+
  scale_color_manual(values = c("#000000","#000000"), guide = "none")+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"),name = "Condition", 
                    labels = c("EZH2+","EZH2+/BrdU+"))+
  ylim(0,0.5)+
  ylab("Cells/DAPI")+
  xlab ("")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.text = element_text(size=20, colour = "black"),
        axis.title = element_text(size = 24, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "none")+
  stat_pvalue_manual(Ezh2_brdu.stat.test, label = "{p.signif}", size = 6)

OBx_EZH2_Brdu_plot

##### Eed fl/fl vs fl/wt EED #####
#Data Management 
EED <- separate(data = EED, col = fileName, into = c("group", "date","stain","location","number"), sep = "-")
EED <- EED %>% 
  mutate(ID = group) 

EEDsum <- EED %>% 
  group_by(ID) %>% 
  summarise(sumNeuron = sum(neuron), sumSus = sum(sustentacular), 
            sumMV = sum(microvillar), sumBas = sum(basal), 
            sumEEDneu = sum(`EED-Neu`), sumEEDsus = sum(`EED-Sus`),
            sumEEDmv = sum(`EED-MV`), sumEEDbas = sum(`EED-Bas`)) %>% 
  mutate(EEDpercent = (sumEEDneu+sumEEDsus+sumEEDbas+sumEEDmv)/(sumNeuron + 
                                                                  sumSus+
                                                                  sumMV+
                                                                  sumBas)*100) %>% 
  mutate(EEDpercentBasal = (sumEEDneu+sumEEDbas)/(sumNeuron + sumBas)*100)
EEDsum <- EEDsum %>% 
  mutate(group = ID) 
EEDsum$group <- gsub('EED', '', EEDsum$group)
EEDsum$group <- gsub('[0-9]', '', EEDsum$group)

#Normality Test
hetEEDsum <- EEDsum %>% 
  filter(group == "HET")
koEEDsum <- EEDsum %>% 
  filter(group == "KO")

shapiro.test(hetEEDsum$EEDpercent)
shapiro.test(koEEDsum$EEDpercent)

var.test(hetEEDsum$EEDpercent, koEEDsum$EEDpercent)

#Data Visualization
EEDmeans <- EEDsum %>% 
  group_by(group) %>% 
  summarise(n=n(),meanEEDpercent = mean(EEDpercentBasal),
            sd = sd(EEDpercentBasal)) %>% 
  mutate(se=sd/sqrt(n))

EEDsum_means <- left_join(EEDmeans, EEDsum, by=c("group"))


EEDplot <- EEDsum_means %>% ggplot(aes(group, EEDpercentBasal,group=group, 
                                       fill = group))+
  stat_summary(aes(fill = group), color = "black",geom = "col", fun = mean, 
               position = position_dodge())+
  geom_errorbar(aes(group=group,
                    ymin=meanEEDpercent-se, ymax=meanEEDpercent+se), 
                position = position_dodge(width = 0.9), width=0.4,
                colour="black",  size=0.5)+
  geom_point(aes(fill = group),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.2))+
  scale_color_manual(values = c("#000000","#000000"), guide = "none")+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"),name = "Group", 
                    labels = c(bquote(EED^"fl/wt"), bquote(EED^"fl/fl")))+
  ylim(0,16)+
  xlab("")+
  ylab("EED+ lineage-traced cells (%)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=24, colour = "black"), 
        axis.text.y = element_text(size=18, colour = "black"),
        axis.title.y = element_text(size=24, colour = "black"))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE), 
                     label = "p.signif", label.x = 1.5, label.y = 14,
                     size = 20, bracket.size = 0.4)+
  scale_x_discrete(labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))
EEDplot

#Adjust p-values for EED with Benjamini-Hochberg correction
eedTtestPvals <- as.numeric(c("0.0068"))
eedTtestAdjPvals <- p.adjust(eedTtestPvals, method = "hochberg")
eed.adj.stat.test <- EEDsum %>% 
  t_test(EEDpercentBasal ~ group, var.equal = TRUE) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  mutate(p.signif = ifelse(p.adj < "0.01",paste ("**"),
                           ifelse(p.adj < "0.05", paste("*"), paste("n.s.,"))))





##### Eed fl/fl vs fl/wt CBX8 cell type #####
# Data Management
CBX8 <- separate(data = CBX8, col = fileName, into = c("group", "date","stain","location","number"), sep = "-")
CBX8 <- CBX8 %>% 
  mutate(ID = group) 
CBX8sum <- CBX8 %>% 
  group_by(ID) %>% 
  summarise(sumNeuron = sum(neuron), sumSus = sum(sustentacular), 
            sumMV = sum(microvillar), sumBas = sum(basal), 
            sumCBX8neu = sum(`CBX8-neu`)) %>% 
  mutate(CBX8percent = (sumCBX8neu)/(sumNeuron)*100)
CBX8sum <- CBX8sum %>% 
  mutate(group = ID) 
CBX8sum$group <- gsub('CBX8', '', CBX8sum$group)
CBX8sum$group <- gsub('[0-9]', '', CBX8sum$group)

cellCounts <- CBX8sum %>% 
  group_by(ID) %>% 
  mutate(sumTotal = sum(sumNeuron, sumSus, sumMV, sumBas)) %>% 
  mutate(Neuron = sumNeuron/sumTotal*100) %>% 
  mutate(Sustentacular = sumSus/sumTotal*100) %>% 
  mutate(Microvillar = sumMV/sumTotal*100) %>% 
  mutate(Basal = sumBas/sumTotal*100) %>% 
  gather("cellType", "percent", 10:13)

neuCounts <- cellCounts %>% 
  filter(cellType == "Neuron")

susCounts <- cellCounts %>% 
  filter(cellType == "Sustentacular")

basCounts <- cellCounts %>% 
  filter(cellType == "Basal")

mvCounts <- cellCounts %>% 
  filter(cellType == "Microvillar")
cellCountsAverages <- cellCounts %>% 
  group_by(group, cellType) %>% 
  summarise(n=n(),meanPercent = mean(percent), 
            sd = sd(percent))%>%
  mutate(se=sd/sqrt(n)) 

cellCounts_sum_mean <- left_join(cellCountsAverages, cellCounts, by=c("group","cellType"))

#Normality Test
hetNeurons <- cellCounts %>% 
  filter(group == "EEDHET") %>% 
  filter(cellType == "Neuron")
koNeurons <- cellCounts %>% 
  filter(group == "EEDKO") %>% 
  filter(cellType == "Neuron")
hetSus <- cellCounts %>% 
  filter(group == "EEDHET") %>% 
  filter(cellType == "Sustentacular")
koSus <- cellCounts %>% 
  filter(group == "EEDKO") %>% 
  filter(cellType == "Sustentacular")

shapiro.test(hetNeurons$percent) 
ggdensity(hetNeurons$percent)
shapiro.test(koNeurons$percent) 
ggdensity(koNeurons$percent)
shapiro.test(hetSus$percent) 
ggdensity(hetSus$percent)
shapiro.test(koSus$percent)
ggdensity(koSus$percent)
shapiro.test(cellCounts$percent) 

var.test(hetNeurons$percent, koNeurons$percent)
var.test(hetSus$percent, koSus$percent)

#T-test with Benjamini-Hochberg adjusted p-values
counts.adj.stat.test <- cellCounts %>% 
  group_by(cellType) %>% 
  t_test(percent ~ group, var.equal = TRUE) %>% 
  adjust_pvalue(method = "hochberg") %>% 
  add_xy_position(x = "cellType", dodge = 0.8) %>% 
  mutate(p.signif = ifelse(p.adj < "0.05",paste ("*"),paste("n.s.")))

#Data Visualization
cellCountsplot <- cellCounts_sum_mean %>% ggplot(aes(cellType, percent, group= group,
                                                     color = group))+
  stat_summary(aes(fill = group),geom = "col", fun = mean, 
               position = position_dodge())+
  geom_errorbar(aes(x=cellType,group=group, 
                    ymin=meanPercent-se, ymax=meanPercent+se), 
                position = position_dodge(width = 0.9), width=0.4,
                colour="black",  size=0.5)+
  geom_point(aes(fill=group),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.2))+
  scale_color_manual(values = c("#000000","#000000"), guide = "none")+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"),name = "Group", 
                    labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))+
  ylim(0,65)+
  ylab("Lineage-traced cells (%)")+
  xlab ("Cell type")+
  theme_bw()+
  theme(axis.ticks.x=element_blank(),
        axis.text = element_text(size=20, colour = "black"),
        axis.title = element_text(size = 24, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))+
  stat_pvalue_manual(counts.adj.stat.test, label = "{p.signif}", size = 6)

cellCountsplot

##### Eed fl/fl vs fl/wt SOX9 ####
#Data Management
sox9_sum <- Sox9_EEDKOHET %>% 
  group_by(Animal) %>% 
  summarise(sox9Sum = sum(`SOX9+, tdTom+`), allTomSum = sum(`all tdTom+`),
            totalEpLength = sum(`mm ep`)) %>% 
  mutate(percentSox9tom = sox9Sum/allTomSum*100) %>% 
  mutate(group = Animal) %>% 
  mutate(tdTomPerEpLength = allTomSum/totalEpLength) %>% 
  mutate(logTdTomPerEpLength = log(tdTomPerEpLength))

sox9_sum <- sox9_sum %>% 
  filter(Animal != "EEDKO2")

sox9_sum$group <- gsub('EED', '', sox9_sum$group)
sox9_sum$group <- gsub('[0-9]', '', sox9_sum$group)

meanSdSox9 <- sox9_sum %>% 
  group_by(group) %>% 
  summarise(n=n(),meanPercentSox9 = mean(percentSox9tom), 
            sd = sd(percentSox9tom))%>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ic=1.96*sd)

sox9_sum_mean <- left_join(meanSdSox9, sox9_sum, by="group")

#Normality Test
koSox9 <- sox9_sum%>% 
  filter(group == "KO")

hetSox9 <- sox9_sum %>% 
  filter(group=="HET")

shapiro.test(koSox9$percentSox9tom)
qqnorm(koSox9$percentSox9tom, pch = 1, frame = FALSE)
qqline(koSox9$percentSox9tom, col = "steelblue", lwd = 2)

shapiro.test(hetSox9$percentSox9tom)
qqnorm(hetSox9$percentSox9tom)
qqline(hetSox9$percentSox9tom, col = "steelblue", lwd = 2)

#T-test
sox9stat.test <- sox9_sum %>% 
  t_test(percentSox9tom ~ group, var.equal = TRUE) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p < "0.05", paste("*"), paste("n.s."))))

#Data Management
SOX9plot <- sox9_sum_mean %>% ggplot(aes(group, percentSox9tom,
                                         group=group, color = group))+
  geom_bar(aes(fill = group), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = group),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.3))+
  geom_errorbar(aes(x=group, ymin=meanPercentSox9-se, ymax=meanPercentSox9+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"), guide = "none")+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"),name = "group")+
  ylim(0,40)+
  xlab("")+
  ylab("SOX9+ lineage-traced cells (%)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=24, 
                                                             colour = "black"), 
        axis.text.y = element_text(size=18, colour = "black"),
        axis.title.y = element_text(size=24, colour = "black"))+
  scale_x_discrete(labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE),
                     label = "p.signif", label.x=1.3, label.y=36,
                     size=20)

SOX9plot

##### Eed fl/fl vs fl/wt Ki67 #####
#Data Management
Ki67_EEDKOHET <- EED_KO_Het_Ki67 %>% 
  filter(!is.na(Sample)) %>% 
  select(Sample,StainSection,Animal,Genotype,`total No. of mCherry Cells`, 
         `Merge (red/green)`) 
Ki67_sum <- EED_KO_Het_Ki67 %>% 
  group_by(Animal) %>% 
  summarise(Ki67Sum = sum(`Merge (red/green)`),
            allTomSum=sum(`total No. of mCherry Cells`))  %>% 
  mutate(percentKi67tom = Ki67Sum/allTomSum*100) %>% 
  mutate(group = Animal)
Ki67_sum$group <- gsub('EED', '', Ki67_sum$group)
Ki67_sum$group <- gsub('[0-9]', '', Ki67_sum$group)  
Ki67_sum$group <- gsub('-', '', Ki67_sum$group)

meanSdKi67 <- Ki67_sum %>% 
  group_by(group) %>% 
  summarise(n=n(),meanPercentKi67 = mean(percentKi67tom), 
            sd = sd(percentKi67tom))%>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ic=1.96*sd)
ki67_sum_mean <- left_join(meanSdKi67, Ki67_sum, by="group")

#Normality Test
koKi67 <- Ki67_sum %>% 
  filter(group == "KO") 
hetKi67 <- Ki67_sum %>% 
  filter(group == "Het") 

shapiro.test(koKi67$percentKi67tom)
qqnorm(koKi67$percentKi67tom, pch = 1, frame = FALSE)
qqline(koKi67$percentKi67tom, col = "steelblue", lwd = 2)
shapiro.test(hetKi67$percentKi67tom)
qqnorm(hetKi67$percentKi67tom, pch = 1, frame = FALSE)
qqline(hetKi67$percentKi67tom, col = "steelblue", lwd = 2)

#T-test
ki67stat.test <- ki67_sum_mean %>% 
  t_test(percentKi67tom ~ group, var.equal = TRUE) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p < "0.05", paste("*"), paste("n.s."))))

#Data Visualization
Ki67plot <- ki67_sum_mean %>% ggplot(aes(group, percentKi67tom,
                                         group=group, color = group))+
  geom_bar(aes(fill = group), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = group),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.3))+
  geom_errorbar(aes(x=group, ymin=meanPercentKi67-se, ymax=meanPercentKi67+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "group", 
                     labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,45)+
  xlab("")+
  ylab("Ki67+ 
       lineage-traced cells (%)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=20, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16, colour = "black"))+
  scale_x_discrete(labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE), 
                     label = "p.signif",label.x.npc = "center", label.y=44,
                     size=8)
Ki67plot
##### Eed fl/fl vs fl/wt Casp-3 #####
#Data Management
caspTom_sum <- EED_KO_Het_Caspase_tom %>% 
  group_by(Animal) %>% 
  summarise(CaspaseSum = sum(casp_tomato),
            totalEpLength=sum(length),
            totalTom = sum(total_tomato),
            neuBasTom = sum(bas_neu_tomato))  %>% 
  mutate(percentCaspTom= CaspaseSum/totalTom*100) %>% 
  mutate(group = Animal)

caspTom_sum$group <- gsub('EED', '', caspTom_sum$group)
caspTom_sum$group <- gsub('[0-9]', '', caspTom_sum$group)  
caspTom_sum$group <- gsub('-', '', caspTom_sum$group)

meanSdcaspTom <- caspTom_sum %>% 
  group_by(group) %>% 
  summarise(n=n(),
            meanPercentCaspTom = mean(percentCaspTom), 
            sd = sd(percentCaspTom))%>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ic=1.96*sd)

#Normality test
koCaspase <- caspTom_sum %>% 
  filter(group == "KO") 
hetCaspase <- caspTom_sum %>% 
  filter(group == "HET") 

shapiro.test(koCaspase$percentCaspTom)
qqnorm(koCaspase$percentCaspTom, pch = 1, frame = FALSE)
qqline(koCaspase$percentCaspTom, col = "steelblue", lwd = 2)
shapiro.test(hetCaspase$percentCaspTom)
qqnorm(hetCaspase$percentCaspTom, pch = 1, frame = FALSE)
qqline(hetCaspase$percentCaspTom, col = "steelblue", lwd = 2)

#T-test
caspTom.stat.test <- caspTom_sum %>% 
  t_test(percentCaspTom ~ group, var.equal = TRUE) %>% 
  mutate(p.signif = ifelse(p < "0.01",paste ("**"),
                           ifelse(p < "0.05", paste("*"), paste("n.s."))))

#Data Visualization
caspTom_sum_mean <- left_join(meanSdcaspTom, caspTom_sum, by="group")

caspTomPlot <- caspTom_sum_mean %>% ggplot(aes(group, percentCaspTom,
                                               group=group, color = group))+
  geom_bar(aes(fill = group), position = position_dodge(), stat = "summary", fun = "mean")+
  geom_point(aes(fill = group),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.3))+
  geom_errorbar(aes(x=group, ymin=meanPercentCaspTom-se, ymax=meanPercentCaspTom+se), 
                color = "black", width=0.3)+
  scale_color_manual(values = c("#000000","#000000"),name = "group", 
                     labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))+
  scale_fill_manual(values = c("#FFFFFF","#D3D3D3"))+
  ylim(0,5)+
  xlab("")+
  ylab("Casp3+ 
       lineage-traced cells (%)")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=20, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16, colour = "black"))+
  scale_x_discrete(labels = expression(italic(Eed^"fl/wt"),italic(Eed^"fl/fl")))+
  stat_compare_means(method = "t.test",method.args = list(var.equal = TRUE), 
                     label = "p.signif",label.x.npc = "center", label.y=4,
                     size=8)
caspTomPlot



##### Non-neuronal OE cell morphological characteristics #####
#Data Management
shape_mean_area_circ <- shape %>% 
  group_by(Animal, Celltype) %>% 
  summarise(meanArea = mean(Area),
            meanCirc = mean(Circularity))

shape_mean_apWid_proc = shape %>% 
  filter(Celltype == c("Sus","MV")) %>% 
  group_by(Animal, Celltype) %>%
  summarise(meanApWid = mean(as.numeric(Apical_width)),
            meanProc = mean(as.numeric(Processes)))
sus_apwid_proc <- shape_mean_apWid_proc %>% 
  filter(Celltype == "Sus")
mv_apwid_proc <- shape_mean_apWid_proc %>% 
  filter(Celltype == "MV")

mean_sd_area_circ <- shape_mean_area_circ %>% 
  group_by(Celltype) %>% 
  summarise(totMeanArea = mean(meanArea),
            totSdArea = sd(meanArea),
            totMeanCirc = mean(meanCirc),
            totSdCirc = sd(meanCirc))

mean_sd_apwid <- shape_mean_apWid_proc %>% 
  group_by(Celltype) %>% 
  summarise(totMeanApwid = mean(meanApWid),
            totSdApwid = sd(meanApWid))

sus_area_circ <- shape_mean_area_circ %>% 
  filter(Celltype=="Sus")
mv_area_circ <- shape_mean_area_circ %>% 
  filter(Celltype=="MV")
bas_area_circ <- shape_mean_area_circ %>% 
  filter(Celltype=="Bas")
#Apical Width Normality Test
shapiro.test(mv_apwid_proc$meanApWid)
shapiro.test(sus_apwid_proc$meanApWid)
var.test(sus_apwid_proc$meanApWid, mv_apwid_proc$meanApWid) 
var.test(sus_apwid_proc$meanApWid, sus_apwid_proc$meanApWid) 
#Apical width T-test
t.test(sus_apwid_proc$meanApWid, mv_apwid_proc$meanApWid, var.equal = TRUE) 
p.adjust(1.45e-11, method = "bonferroni", n = 7) 

#Apical Width data visualization
apWid_sum_means <- left_join(shape_mean_apWid_proc, mean_sd_apwid, by=c("Celltype"))

apWidPlot <- apWid_sum_means %>% 
  ggplot(aes(Celltype,meanApWid))+
  stat_summary(aes(fill=Celltype),geom = "col",fun = "mean", colour = "black")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=24, colour = "black"), 
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size=20, colour = "black"),
        axis.title.y = element_text(size=24, colour = "black"))+
  geom_errorbar(aes(group=Celltype,
                    ymin=totMeanApwid-totSdApwid, ymax=totMeanApwid+totSdApwid), 
                position = position_dodge(width = 0.9), width=0.4,
                colour="black",  size=0.5)+
  geom_point(aes(fill = Celltype),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.2))+
  scale_fill_manual(values = c("#CFCFCF","#6A6A6A"),name = "Group")+
  ylim(0,8)+
  xlab("Cell Type")+
  ylab(expression("Mean Apical Width"~"("*mu*m*")"))+
  geom_bracket(xmin = "MV",
               xmax = "Sus",
               y.position = 7.8,
               label = "***",
               tip.length = 0.07,
               label.size = 14,
               vjust = 0.5)
apWidPlot

#Area  Normality test
shapiro.test(sus_area_circ$meanArea)
shapiro.test(mv_area_circ$meanArea)
shapiro.test(bas_area_circ$meanArea)
var.test(sus_area_circ$meanArea, mv_area_circ$meanArea)
var.test(sus_area_circ$meanArea, bas_area_circ$meanArea)
var.test(mv_area_circ$meanArea, bas_area_circ$meanArea)

#Area T-test
t.test(sus_area_circ$meanArea, mv_area_circ$meanArea, var.equal = TRUE)
p.adjust(0.02128, method = "bonferroni", n = 7) 
t.test(sus_area_circ$meanArea, bas_area_circ$meanArea, var.equal = TRUE)
p.adjust(1.44e-08, method = "bonferroni", n=7) 
t.test(mv_area_circ$meanArea, bas_area_circ$meanArea, var.equal = TRUE)
p.adjust(1.425e-06, method = "bonferroni", n=7) 

#Area Data visualization
areaPlot <-circ_sum_means %>% 
  ggplot(aes(Celltype, meanArea))+
  stat_summary(aes(fill=Celltype), geom = "col", fun = "mean", colour = "black")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=24, colour = "black"), 
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size=20, colour = "black"),
        axis.title.y = element_text(size=24, colour = "black"))+
  geom_errorbar(aes(group=Celltype,
                    ymin=totMeanArea-totSdArea, ymax=totMeanArea+totSdArea),
                position = position_dodge(width = 0.9), width=0.4)+
  geom_point(aes(fill = Celltype),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.2))+
  scale_fill_manual(values = c("#FFFFFF","#CFCFCF","#6A6A6A"),name = "Group")+
  ylim(0, 150)+
  ylab(expression("Mean Area"~"("*mu*m^2*")"))+
  geom_bracket(xmin = "Basal",
               xmax = "MV",
               y.position = 120,
               label = "***",
               tip.length = 0.07,
               label.size = 14,
               vjust = 0.5)+
  geom_bracket(xmin = "Basal",
               xmax = "Sus",
               y.position = 145,
               label = "***",
               tip.length = 0.07,
               label.size = 14,
               vjust = 0.5)

areaPlot

#Circularity Normality test
shapiro.test(sus_area_circ$meanCirc)
shapiro.test(mv_area_circ$meanCirc)
shapiro.test(bas_area_circ$meanCirc)
var.test(sus_area_circ$meanCirc, mv_area_circ$meanCirc)
var.test(sus_area_circ$meanCirc, bas_area_circ$meanCirc)
var.test(mv_area_circ$meanCirc, bas_area_circ$meanCirc)

#Circularity T-test
t.test(sus_area_circ$meanCirc, mv_area_circ$meanCirc, var.equal = TRUE)
t.test(sus_area_circ$meanCirc, bas_area_circ$meanCirc, var.equal = FALSE)
p.adjust(9.661e-06, method = "bonferroni", n=7) 
t.test(mv_area_circ$meanCirc, bas_area_circ$meanCirc, var.equal = FALSE)
p.adjust(5.448e-07, method = "bonferroni", n=7) 

#Circularity data visualization
circ_sum_means <- left_join(shape_mean_area_circ, mean_sd_area_circ, by=c("Celltype"))
circ_sum_means$Celltype<-str_replace(circ_sum_means$Celltype, "Bas","Basal")

circPlot <- circ_sum_means %>% 
  ggplot(aes(Celltype,meanCirc))+
  stat_summary(aes(fill=Celltype),geom = "col",fun = "mean", colour = "black")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size=24, colour = "black"), 
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size=20, colour = "black"),
        axis.title.y = element_text(size=24, colour = "black"))+
  geom_errorbar(aes(group=Celltype,
                    ymin=totMeanCirc-totSdCirc, ymax=totMeanCirc+totSdCirc), 
                position = position_dodge(width = 0.9), width=0.4,
                colour="black",  size=0.5)+
  geom_point(aes(fill = Celltype),size = 3, color = "black",
             position = position_jitterdodge(dodge.width = 0.9, 
                                             jitter.width = 0.2))+
  scale_fill_manual(values = c("#FFFFFF","#CFCFCF","#6A6A6A"),name = "Celltype",
                    labels = c("Basal","MV","Sus"))+
  ylim(0,1.05)+
  xlab("Cell Type")+
  ylab("Mean Circularity")+
  geom_bracket(xmin = "Basal",
               xmax = "MV",
               y.position = 0.92,
               label = "***",
               tip.length = 0.07,
               label.size = 14,
               vjust = 0.5)+
  geom_bracket(xmin = "Basal",
               xmax = "Sus",
               y.position = 1,
               label = "***",
               tip.length = 0.07,
               label.size = 14,
               vjust = 0.5)
circPlot


