### NMDS with Envfit for environmental variables ###

### Get and set work directory ###
rm(list=ls(all=TRUE)) #remove all created objects
getwd()
setwd("C:/Users/gui_l/OneDrive/Mestrado/Artigo 1 - De novo/Scripts")
path <- "C:/Users/gui_l/OneDrive/Mestrado/Artigo 1 - De novo/Scripts"
list.files(path)

### Install pairwiseAdonis package ###
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

### Open required packages ###
library (vegan)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(cluster)
library(dplyr)
library(tidyr)
library(ARTool)
library(lsmeans)
library(pairwiseAdonis)
library(readxl)

### Load data files ###
chem <- read_excel("dados-paper.xlsx",
                   sheet = "Chemical")
chem$Group <- factor(chem$Group, levels = unique(chem$Group))

### PERMANOVA ###
adonis2(chem[4:15] ~ Group*State, data = chem, permutations = 999, method = "gower", p.adjusted = p.adjust(p.value,method="fdr"))
#Df SumOfSqs      R2      F Pr(>F)    
#Group        2  0.31470 0.32319 5.2855  0.001 ***
#State        1  0.09736 0.09999 3.2705  0.019 *  
#Group:State  2  0.20443 0.20994 3.4335  0.003 ** 
#Residual    12  0.35725 0.36688                  
#Total       17  0.97374 1.00000

#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.adonis(chem[4:15], chem$Group, sim.function = "vegdist", sim.method = "gower", p.adjust.m = "fdr", perm = 999)
#pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#1 FC vs PL  1 0.2434621 4.630329 0.3164884   0.008      0.012   .
#2 FC vs CM  1 0.2855863 5.237366 0.3437186   0.008      0.012   .
#3 PL vs CM  1 0.0815821 1.337636 0.1179819   0.238      0.238      

### Create matrix and order as NMDS ###
matrix.nmds<-metaMDS(chem[4:15], k=2, distance = 'gower', noshare=FALSE, autotransform=FALSE) # k=2 transforma em nmds de duas dimensoes, autotransform = true para dados de sequencias
matrix.nmds$stress #0.1552347
#A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation. 

NMDS1<-matrix.nmds$points[,1]
NMDS2<-matrix.nmds$points[,2]
NMDS=data.frame(NMDS1=NMDS1,NMDS2=NMDS2,Group=chem$Group,State=chem$State)

### Plot NMDS ###
ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Group, shape=State), size=5,  alpha=0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_colour_brewer(palette = "Set1") +
  theme_gray() +
  theme(legend.position="right",
        legend.text = element_text(size=10),
        axis.text.x= element_text(angle= 0, size= 10),
        axis.text.y= element_text(size= 10),
        axis.title= element_text(size= 10, face= "bold"))

#######################################################
   ### NMDS for bacterial community with envfit ###
mds_gen <- read_excel("dados-paper.xlsx",
                   sheet = "MDS_Gen")
mds_gen$Group <- factor(mds_gen$Group, levels=unique(mds_gen$Group))

### NMDS ordenation ###
matrix.mds_gen<-metaMDS(mds_gen[4:570], k=2, distance = 'bray', noshare=FALSE, autotransform=FALSE)
matrix.mds_gen$stress #0.02777157 >0.5 = excellent representation

NMDS1<-matrix.mds_gen$points[,1]
NMDS2<-matrix.mds_gen$points[,2]
NMDS=data.frame(NMDS1=NMDS1,NMDS2=NMDS2,Group=mds_gen$Group,State=mds_gen$State)

### Aplicar vetores envfit ###
vec.sp<-envfit(matrix.mds_gen, chem, perm=999)
vec.sp
#NMDS1    NMDS2     r2 Pr(>r)    
#pH    0.03126  0.99951 0.0848  0.507    
#M.O. -0.06154  0.99810 0.0392  0.702    
#P     0.41775  0.90856 0.2207  0.170    
#K     0.77954 -0.62636 0.0969  0.461    
#Ca    0.23294  0.97249 0.2561  0.089 .  
#Mg    0.95690  0.29043 0.0191  0.869    
#H+Al  0.38212  0.92411 0.0704  0.572    
#S    -0.12395 -0.99229 0.0485  0.715    
#Fe   -0.56148 -0.82749 0.4034  0.016 *  
#Mn    0.99922  0.03942 0.7655  0.001 ***
#Cu    0.33730  0.94140 0.5184  0.006 ** 
#Zn    0.70799  0.70622 0.6724  0.001 ***

#Goodness of fit:
#r2 Pr(>r)    
#ID    0.9962  0.001 ***
#Group 0.3976  0.016 *  
#State 0.3652  0.007 ** 

vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*vec.sp$vectors$r)
vec.sp.df$species<-rownames(vec.sp.df)

### Plot NMDS ###
ggplot(NMDS, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color=Group, shape=State), size=5,  alpha=0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +  
  scale_colour_brewer(palette = "Set1") +
  geom_segment(data=vec.sp.df[c("Ca","Fe","Mn","Cu","Zn"),],aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),colour="grey") +
  geom_text(data=vec.sp.df[c("Ca","Fe","Mn","Cu","Zn"),],aes(x=NMDS1,y=NMDS2,label=species),size=5)+
  theme_gray()+
  theme(legend.position="right",
        legend.text = element_text(size=10),
        axis.text.x= element_text(angle= 0, size= 10),
        axis.text.y= element_text(size= 10),
        axis.title= element_text(size= 10, face= "bold"))

### PERMANOVA ###
adonis2(mds_gen[4:570] ~ Group*State, data = mds_gen, permutations = 999, method = "bray", p.adjusted = p.adjust(p.value,method="tukey"))
#Df SumOfSqs      R2      F Pr(>F)    
#Group        2   2.3139 0.37612 51.798  0.001 ***
#State        1   1.7021 0.27667 76.204  0.001 ***
#Group:State  2   1.8680 0.30364 41.816  0.001 ***
#Residual    12   0.2680 0.04357                  
#Total       17   6.1520 1.00000                  
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


################################
    ### Diversity index ###

diversity(mds_gen[4:570], index = "shannon", MARGIN = 1, base = exp(1))
#[1] 3.927380 3.887648 3.743499 4.424698 4.348587 4.252712 3.158315 3.230927
#[9] 3.083951 3.201866 3.123318 3.040068 1.932856 2.596174 2.270015 4.065777
#[17] 4.029217 3.976788

### Richness ###
specnumber(mds_gen[4:570], MARGIN = 1)
#[1] 169 167 140 204 180 137 128 122 115 151 126 145  48  91  70 161 165 155

### Evenness (Pielou) ###
shannon <- diversity(mds_gen[4:570])
shannon/log(specnumber(mds_gen[4:570]))
#[1] 0.7655863 0.7596039 0.7575415 0.8320042 0.8374009 0.8643757 0.6509265
#[8] 0.6725464 0.6499464 0.6381677 0.6458098 0.6108560 0.4992911 0.5755387
#[15] 0.5343103 0.8001285 0.7891225 0.7885093


#######################
  ### Statistics ###

### Diversidade ###
diversidade$Group <- as.factor(diversidade$Group); diversidade$State <- as.factor(diversidade$State)
diversity_estatistica <- art(Diversity ~ Group * State, data = diversidade, check.errors.are.factors = TRUE); summary(diversity_estatistica); anova(diversity_estatistica)
