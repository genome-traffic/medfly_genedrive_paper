library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)
library("extrafont")
loadfonts()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# ---------- data input ----------

t<-read.table("tra_homing.csv", header=T, sep=",")
t <- subset(t, t$Male_parent == "VasaGD/tra+")
t$percent <- with(t, ave(Number, Replicate, FUN = prop.table) * 100)

t_neg <- subset(t, t$Fluorescence == "DsRed negative")
t_neg$percent <- with(t_neg, ave(Number, Replicate, FUN = prop.table) * 100)
t_neg2 <- ddply(t_neg, c("Sex_Phenotype"), summarise, sum = sum(Number))
t_neg2$percent <- with(t_neg2, ave(sum, FUN = prop.table) * 100)

t_all <- ddply(t, c("Sex_Phenotype"), summarise, sum = sum(Number))
t_all$percent <- with(t_all, ave(sum, FUN = prop.table) * 100)

t_F <- ddply(t, c("Replicate","Fluorescence"), summarise, sum = sum(Number))
t_F$percent <- with(t_F, ave(sum, Replicate, FUN = prop.table) * 100)

# ---------- plot ----------

box <- ggplot(data=t, aes(x=Sex_Phenotype, y=percent, fill=Sex_Phenotype, color=Sex_Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#FF6666","#C0C0C0","#3399FF")) +
  scale_color_manual(values=c("#990000","#404040","#003366")) + 
  theme_classic(base_size = 12) +
  scale_y_continuous(n.breaks=10) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1)+
 facet_wrap(~Fluorescence) + theme(legend.position="none")
plot(box)

boxp_homing <- ggplot(data=t_F, aes(x=Fluorescence, y=percent, fill=Fluorescence, color=Fluorescence)) +
  geom_hline(yintercept = 50, size=0.25, alpha=0.5, linetype="dotted") +
  ylim(40,90) +
  geom_boxplot() +
  scale_fill_manual(values=c("#606060","#606060")) +
  scale_color_manual(values=c("#000000","#000000")) + 
  theme_classic(base_size = 12) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1)+
  theme(legend.position="none")
plot(boxp_homing)


# ---------- stats ----------

an <- aov(h$transmission ~ h$condition)
summary(an)
TukeyHSD(an)
