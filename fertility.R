library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(report)
library(DescTools)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# ---------- data input ----------
f<-read.table("fertility.csv", header=T, sep=",")
f$HatchingRate <- (f$Hatched/f$Embryos)*100
f$Experiment <- paste(f$Drive, f$Crosstype)

fn <- f[c("Drive","Crosstype", "Embryos")]
(fn)

fn %>%
  group_by(Drive, Crosstype) %>%
  dplyr::summarise(across(Embryos,sum))

fn %>%
  group_by(Drive) %>%
  dplyr::summarise(across(Embryos,sum))

# ---------- plot ----------

boxp_fbar <- ggplot(data=f, aes(x=Experiment, y=Embryos, fill=Experiment,color=Experiment)) +
  stat_summary(geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               width = .3, size = 0.8) +
  geom_bar(stat = "summary", width = .8, size = .8) +
  scale_fill_manual(values=c("#E0E0E0","#606060","#606060","#606060","#606060","#606060","#606060","#606060","#606060","#606060")) +
  scale_color_manual(values=c("#404040","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")) + 
  theme_classic(base_size = 16)+
  scale_y_continuous(limits = c(0,550), expand = c(0, 0)) 
plot(boxp_fbar)

boxp_fbarH <- ggplot(data=f, aes(x=Experiment, y=HatchingRate, fill=Experiment,color=Experiment)) +
  stat_summary(geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               width = .3, size = 0.8) +
  geom_bar(stat = "summary", width = .8, size = .8) +
  scale_fill_manual(values=c("#E0E0E0","#606060","#606060","#606060","#606060","#606060","#606060","#606060","#606060","#606060")) +
  scale_color_manual(values=c("#404040","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")) + 
  theme_classic(base_size = 16)+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) 
plot(boxp_fbarH)

ggarrange(boxp_fbar,boxp_fbarH,  
          labels = c("A", "B"),legend="none",
          ncol = 2, nrow = 1, widths = c(1, 1))


#----------------------------------- Stats  ----------------

f.aov <- aov(Embryos ~ Experiment, data = f)
#summary(f.aov)
#report(f.aov)
TukeyHSD(f.aov)

fH.aov <- aov(HatchingRate ~ Experiment, data = f)
#summary(f.aov)
#report(f.aov)
TukeyHSD(fH.aov)

DunnettTest(x=f$Embryos, g=f$Experiment)
DunnettTest(x=f$HatchingRate, g=f$Experiment)

