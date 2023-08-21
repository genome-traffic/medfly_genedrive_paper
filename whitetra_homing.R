library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# ---------- data input ----------

w<-read.table("tra_homing.csv", header=T, sep=",")
w <- subset(w, w$Male_parent != "VasaGD/tra+")
w$Female_parent <- NULL
w$Male_parent[w$Male_parent == "VasaGD;tra+/w+;tra+"] <- "male cross"
w$Male_parent[w$Male_parent == "w+;tra+/w+;tra+"] <- "female cross"
colnames(w)[1] <- "Crosstype"
w_reps <- ddply(w, c("Replicate"), summarise, Number = sum(Number))

sum(w[w$Crosstype == 'male cross',]$Number)
sum(w[w$Crosstype == 'female cross',]$Number)

# ---------- overall homing ----------

w_homing <- ddply(w, c("Replicate","Fluorescence","Crosstype"), summarise, Number = sum(Number))
w_homing$percent <- with(w_homing, ave(Number, Replicate, Crosstype ,FUN = prop.table) * 100)
w_homing <- subset(w_homing, w_homing$Fluorescence == "DsRed positive")

mean(w_homing[w_homing$Crosstype == "female cross",]$percent)
mean(w_homing[w_homing$Crosstype == "male cross",]$percent)

boxp_hom <- ggplot(data=w_homing, aes(x=Crosstype, y=percent, fill=Crosstype, color=Crosstype)) +
  geom_hline(yintercept = 50, size=0.25, alpha=0.5, linetype="dotted") +
  geom_boxplot() +
  scale_fill_manual(values=c("#606060","#606060")) +
  scale_color_manual(values=c("#000000","#000000")) + 
  ylim(40,80) +
  theme_classic(base_size = 12) + theme(legend.position="none") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1)
plot(boxp_hom)

# ---------- overall sex plots ----------

w_sex <- ddply(w, c("Replicate","Sex_Phenotype","Crosstype"), summarise, Number = sum(Number))
w_sex$percent <- with(w_sex, ave(Number, Replicate, Crosstype ,FUN = prop.table) * 100)

sex_means <- ddply(w_sex, c("Sex_Phenotype","Crosstype"), summarise,
                   N    = sum(Number),
                   mean = mean(percent),
                   median = median(percent),
                   sd   = sd(percent),
                   se   = sd / sqrt(N))

boxp_sex <- ggplot(data=w_sex, aes(x=Sex_Phenotype, y=percent, fill=Sex_Phenotype, color=Sex_Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#FF6666","#C0C0C0","#3399FF")) +
  scale_color_manual(values=c("#990000","#404040","#003366")) + 
  theme_classic(base_size = 12) + theme(legend.position="none") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1)+
  facet_wrap(~Crosstype)
plot(boxp_sex)

ggarrange(boxp_hom,boxp_sex,  
          labels = c("A", "B"),
          ncol = 2, nrow = 1, widths = c(1, 2.5))

# ---------- white-eye & fluorescence breakdown ----------

w$percent <- with(w, ave(Number, Replicate, Crosstype, FUN = prop.table) * 100)

boxp <- ggplot(data=w, aes(x=Sex_Phenotype, y=percent, fill=Fluorescence, color=Fluorescence)) +
  geom_boxplot() +
  #scale_fill_manual(values=c("#6666FF","#FF6666")) +
  #scale_color_manual(values=c("#000099","#990000")) +
  scale_fill_manual(values=c("#E0E0E0","#CC0000")) +
  scale_color_manual(values=c("#404040","#660000")) + 
  theme_classic(base_size = 12) +
  facet_grid(vars(Crosstype), vars(White_Phenotype)) + scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=exp(1)), breaks=c(0,1,2,3,4,5,10,20,40,60)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1)
plot(boxp)

sum(w[w$Crosstype == 'male cross' & w$White_Phenotype == 'Mosaic Eyes',]$percent)
sum(w[w$Crosstype == 'male cross' & w$White_Phenotype != 'Red Eyes' & w$Fluorescence =='DsRed positive',]$percent)

