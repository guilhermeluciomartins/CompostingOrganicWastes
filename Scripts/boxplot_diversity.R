### Boxplot for diversity index ###
rm(list=ls(all=TRUE)) #apaga todo o ambiente R

### Definir diretorio de trabalho ###
getwd()
setwd("C:/Users/gui_l/OneDrive/Mestrado/Artigo 1 - De novo/Scripts")
path<-"C:/Users/gui_l/OneDrive/Mestrado/Artigo 1 - De novo/Scripts"
list.files(path)

### Carregar pacotes ###
library(multcompView)
library(tidyverse)
library(tidyr)
library(ggthemes)
library(patchwork)
library(ggplot2)
library(ggpubr)

### Definir cores do boxplot ###
corFC = c("tomato2", "tomato4")
corPL = c("seagreen2", "seagreen4")
corCM = c("steelblue2", "steelblue4")

### Abrir conjunto de dados ###
library(readxl)
shannon <- read_excel("dados-paper.xlsx",
                      sheet = "Diversity")
shannon

### Selecionar dados a serem analisados ###
shannon_FC <-shannon[48:53,]
shannon_FC

shannon_PL <-shannon[54:59,]
shannon_PL

shannon_CM <-shannon[60:65,]
shannon_CM

### Dar nome aos dados a serem analisados ###
ID_FC=shannon_FC$ID
value_FC=shannon_FC$Shannon
group_FC=shannon_FC$Group
state_FC=shannon_FC$State

ID_PL=shannon_PL$ID
value_PL=shannon_PL$Shannon
group_PL=shannon_PL$Group
state_PL=shannon_PL$State

ID_CM=shannon_CM$ID
value_CM=shannon_CM$Shannon
group_CM=shannon_CM$Group
state_CM=shannon_CM$State

### Criar dataframe com os dados ###
data_FC=data.frame(state_FC, value_FC, group_FC)
data_PL=data.frame(state_PL, value_PL, group_PL)
data_CM=data.frame(state_CM, value_CM, group_CM)

### Reordenar a apresentacao dos grupos ###
#data_FC$group_FC = factor(data_FC$group_FC, levels = c("fresh", "composted"))

### Definir quais serao as comparacoes ###
#compair = list(c("fresh", "composted"))

### Gerar bloxplot ###
FC = ggboxplot(data_FC, "state_FC", "value_FC", fill = "state_FC")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corFC)+
  labs(x = NULL, y = "Shannon-Weiner (H)")+
  ylim(2,6.5)+
  rremove("legend")+
  rremove("x.text")+
  facet_wrap("group_FC")

FC

######### 
PL = ggboxplot(data_PL, "state_PL", "value_PL", fill = "state_PL")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corPL)+
  labs(x = NULL, y = NULL)+
  ylim(2,6.5)+
  rremove("legend")+
  rremove("y.text")+
  rremove("x.text")+
  facet_wrap("group_PL")

PL

########
CM =  ggboxplot(data_CM, "state_CM", "value_CM", fill = "state_CM")+
  stat_compare_means(label = "p.signif" , method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corCM)+
  labs(x = NULL, y = NULL)+
  ylim(2,6.5)+
  rremove("legend")+
  rremove("y.text")+
  rremove("x.text")+
  facet_wrap("group_CM")
CM

### Juntar figuras ###
SHANNON <- FC + PL + CM  # exportar em .tiff (700 x 350) 
SHANNON

#  save plot with 600 dpi resolution
dev.print(tiff, "shannon_diversity.tiff", compression = "lzw", res=600, height=10, width=15, units="cm")


#######################################
        ### Simpson index ###

### Abrir conjunto de dados ###
library(readxl)
simpson <- read_excel("dados-paper.xlsx",
                      sheet = "Diversity")
simpson

### Selecionar dados a serem analisados ###
simpson_FC <-simpson[48:53,]
simpson_FC

simpson_PL <-simpson[54:59,]
simpson_PL

simpson_CM <-simpson[60:65,]
simpson_CM

### Dar nome aos dados a serem analisados ###
ID_FC=simpson_FC$ID
value_FC=simpson_FC$Simpson
group_FC=simpson_FC$Group
state_FC=simpson_FC$State

ID_PL=simpson_PL$ID
value_PL=simpson_PL$Simpson
group_PL=simpson_PL$Group
state_PL=simpson_PL$State

ID_CM=simpson_CM$ID
value_CM=simpson_CM$Simpson
group_CM=simpson_CM$Group
state_CM=simpson_CM$State

### Criar dataframe com os dados ###
data_FC=data.frame(state_FC, value_FC, group_FC)
data_PL=data.frame(state_PL, value_PL, group_PL)
data_CM=data.frame(state_CM, value_CM, group_CM)

### Gerar bloxplot ###
FC = ggboxplot(data_FC, "state_FC", "value_FC", fill = "state_FC")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("fresh", "composted")))+
  theme_bw()+
  scale_fill_manual(values = corFC)+
  labs(x = NULL, y = "Simpson (1-D)")+
  ylim(0.7,1)+
  rremove("legend")+
  facet_wrap("group_FC")

FC

######### 
PL = ggboxplot(data_PL, "state_PL", "value_PL", fill = "state_PL")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("fresh", "composted")))+
  theme_bw()+   
  scale_fill_manual(values = corPL)+
  labs(x = NULL, y = NULL)+
  ylim(0.7,1)+
  rremove("legend")+
  rremove("y.text")+
  facet_wrap("group_PL")

PL

########
CM =  ggboxplot(data_CM, "state_CM", "value_CM", fill = "state_CM")+
  stat_compare_means(label = "p.signif" , method = "t.test", comparisons = list(c("fresh", "composted")))+
  theme_bw()+
  scale_fill_manual(values = corCM)+
  labs(x = NULL, y = NULL)+
  ylim(0.7,1)+
  rremove("legend")+
  #rremove("y.text")+
  facet_wrap("group_CM")
CM

### Juntar figuras ###
FC + PL + CM  # exportar em .tiff (700 x 350)
#  save plot with 600 dpi resolution
dev.print(tiff, "simpson_diversity.tiff", compression = "lzw", res=600, height=10, width=15, units="cm")


##################################
        ### Chao1 ###

### Abrir conjunto de dados ###
library(readxl)
chao1 <- read_excel("dados-paper.xlsx",
                      sheet = "Diversity")
chao1

### Selecionar dados a serem analisados ###
chao1_FC <-chao1[48:53,]
chao1_FC

chao1_PL <-chao1[54:59,]
chao1_PL

chao1_CM <-chao1[60:65,]
chao1_CM

### Dar nome aos dados a serem analisados ###
ID_FC=chao1_FC$ID
value_FC=chao1_FC$Chao1
group_FC=chao1_FC$Group
state_FC=chao1_FC$State

ID_PL=chao1_PL$ID
value_PL=chao1_PL$Chao1
group_PL=chao1_PL$Group
state_PL=chao1_PL$State

ID_CM=chao1_CM$ID
value_CM=chao1_CM$Chao1
group_CM=chao1_CM$Group
state_CM=chao1_CM$State

### Criar dataframe com os dados ###
data_FC=data.frame(state_FC, value_FC, group_FC)
data_PL=data.frame(state_PL, value_PL, group_PL)
data_CM=data.frame(state_CM, value_CM, group_CM)

### Gerar bloxplot ###
FC2 = ggboxplot(data_FC, "state_FC", "value_FC", fill = "state_FC")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corFC)+
  labs(x = NULL, y = "Chao1 (S)")+
  ylim(0,1000)+
  rremove("legend")+
  rremove("x.text")+
  facet_null()

FC2

######### 
PL2 = ggboxplot(data_PL, "state_PL", "value_PL", fill = "state_PL")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corPL)+
  labs(x = NULL, y = NULL)+
  ylim(0,1000)+
  rremove("legend")+
  rremove("y.text")+
  rremove("x.text")+
  facet_null()

PL2

########
CM2 =  ggboxplot(data_CM, "state_CM", "value_CM", fill = "state_CM")+
  stat_compare_means(label = "p.signif" , method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corCM)+
  labs(x = NULL, y = NULL)+
  ylim(0,1000)+
  rremove("legend")+
  rremove("y.text")+
  rremove("x.text")+
  facet_null()

CM2

CHAO1 <- FC2 + PL2 + CM2
CHAO1
##################################
### Observed OTU ###

### Abrir conjunto de dados ###
library(readxl)
obs <- read_excel("dados-paper.xlsx",
                    sheet = "Diversity")
obs

### Selecionar dados a serem analisados ###
obs_FC <-obs[48:53,]
obs_FC

obs_PL <-obs[54:59,]
obs_PL

obs_CM <-obs[60:65,]
obs_CM

### Dar nome aos dados a serem analisados ###
ID_FC=obs_FC$ID
value_FC=obs_FC$Chao1
group_FC=obs_FC$Group
state_FC=obs_FC$State

ID_PL=obs_PL$ID
value_PL=obs_PL$Chao1
group_PL=obs_PL$Group
state_PL=obs_PL$State

ID_CM=obs_CM$ID
value_CM=obs_CM$Chao1
group_CM=obs_CM$Group
state_CM=obs_CM$State

### Criar dataframe com os dados ###
data_FC=data.frame(state_FC, value_FC, group_FC)
data_PL=data.frame(state_PL, value_PL, group_PL)
data_CM=data.frame(state_CM, value_CM, group_CM)

### Gerar bloxplot ###
FC3 = ggboxplot(data_FC, "state_FC", "value_FC", fill = "state_FC")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corFC)+
  labs(x = NULL, y = "Observed ASVs")+
  ylim(0,1000)+
  rremove("legend")+
  facet_null()

FC3

######### 
PL3 = ggboxplot(data_PL, "state_PL", "value_PL", fill = "state_PL")+
  stat_compare_means(label = "p.signif", method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corPL)+
  labs(x = NULL, y = NULL)+
  ylim(0,1000)+
  rremove("legend")+
  rremove("y.text")+
  facet_null()

PL3

########
CM3 =  ggboxplot(data_CM, "state_CM", "value_CM", fill = "state_CM")+
  stat_compare_means(label = "p.signif" , method = "t.test", comparisons = list(c("Fresh", "Composted")))+
  theme_bw()+
  scale_fill_manual(values = corCM)+
  labs(x = NULL, y = NULL)+
  ylim(0,1000)+
  rremove("legend")+
  rremove("y.text")+
  facet_null()

CM3

### Juntar figuras ###
OBSERVED <- FC3 + PL3 + CM3  # exportar em .tiff (700 x 350)
OBSERVED

DIVERSITY <- SHANNON / CHAO1 / OBSERVED
DIVERSITY

################################################
### Arrange plots on one page ###
ggarrange(DIVERSITY, MDS,
          labels = c("A", "B"),
          widths = c(1.4,2.1), # numerical vector of relative columns widths
          common.legend = FALSE,
          legend = "right",
          ncol = 2, nrow = 1)


#  save plot with 600 dpi resolution
dev.print(tiff, "FIG1.tiff", compression = "lzw", res=600, height=15, width=30, units="cm")
