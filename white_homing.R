library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(ggbeeswarm)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# ---------- data input ----------

h<-read.table("white_homing.csv", header=T, sep=",")

h$Cross <- paste(h$Male_parent, " x ", h$Female_parent)
h$Test <- grepl("w-", h$Cross, fixed = TRUE)
h$Test <- ifelse(h$Test, 'to_white','to_wt')
h$T <- grepl("Zpg", h$Cross,fixed = TRUE)
h$Promoter[h$T == "TRUE"] <- "zpg"
h$T <- grepl("Nanos", h$Cross,fixed = TRUE)
h$Promoter[h$T == "TRUE"] <- "nanos"
h$T <- grepl("Vasa", h$Cross,fixed = TRUE)
h$Promoter[h$T == "TRUE"] <- "vasa"
h$T <- NULL
h$Crosstype <- str_length(h$Male_parent) > str_length(h$Female_parent)
h$Crosstype[h$Crosstype == "TRUE"] <- "male cross"
h$Crosstype[h$Crosstype == "FALSE"] <- "female cross"

# ------------- homing stats -----------
totest <- xtabs(Number ~ Fluorescence, subset(h, c(h$Crosstype == "male cross")
                                              & h$Promoter=="zpg" 
                                              & h$Test=="to_wt"))
print(totest)
chisq.test(totest)

# ------------- homing plots -----------

h_homing <- ddply(h, c("Replicate","Fluorescence","Test","Promoter","Crosstype"), summarise, Number = sum(Number))
h_homing$percent <- with(h_homing, ave(Number, Replicate, Crosstype, Promoter, Test ,FUN = prop.table) * 100)
h_homing <- subset(h_homing, h_homing$Fluorescence == "DsRed positive")

h_homing_means <- ddply(h_homing, c("Test","Promoter","Crosstype"), summarise,
              N    = length(percent),
              mean = mean(percent),
              median = median(percent),
              sd   = sd(percent),
              se   = sd / sqrt(N))

sum(h[h$Promoter == 'nanos' & h$Crosstype == 'male cross',]$Number)
sum(h[h$Promoter == 'nanos' & h$Crosstype == 'female cross',]$Number)
sum(h[h$Promoter == 'vasa' & h$Crosstype == 'male cross',]$Number)
sum(h[h$Promoter == 'vasa' & h$Crosstype == 'female cross',]$Number)
sum(h[h$Promoter == 'zpg' & h$Crosstype == 'male cross',]$Number)
sum(h[h$Promoter == 'zpg' & h$Crosstype == 'female cross',]$Number)


boxp_hom <- ggplot(data=h_homing, aes(x=Promoter, y=percent, fill=Test, color=Test)) +
  geom_hline(yintercept = 50, size=0.25, alpha=0.5, linetype="dotted") + geom_boxplot() +
  scale_fill_manual(values=c("#E0E0E0","#606060")) +
  scale_color_manual(values=c("#404040","#000000")) + 
  theme_classic(base_size = 12) +
  facet_grid(Test ~ Crosstype) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1)
plot(boxp_hom)


# -------------------------- white-eye phenotype plots ------------------------

h_neg <- subset(h, h$Fluorescence == "DsRed negative")
h_de_neg <- ddply(h_neg, c("Replicate","White_Phenotype","Test","Promoter","Crosstype"), summarise, Number = sum(Number))
h_de_neg$percent <- with(h_de_neg, ave(Number, Replicate, Crosstype, Promoter, Test ,FUN = prop.table) * 100)
h_de_neg <- subset(h_de_neg, h_de_neg$White_Phenotype != "Red Eyes")

h_de_neg[h_de_neg$White_Phenotype == 'Mosaic Eyes' & h_de_neg$Test == 'to_wt' & h_de_neg$Crosstype == 'male cross' & h_de_neg$Promoter == 'vasa' ,]

h_act_neg_means <- ddply(h_de_neg, c("White_Phenotype","Test","Promoter","Crosstype"), summarise,
                        N    = length(percent),
                        mean = mean(percent),
                        median = median(percent),
                        sd   = sd(percent),
                        se   = sd / sqrt(N))

boxp_neg <- ggplot(data=h_de_neg, aes(x=White_Phenotype, y=percent, fill=Test, color=Test)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E0E0E0","#606060")) +
  scale_color_manual(values=c("#404040","#000000")) + 
  theme_classic(base_size = 12) +
  facet_grid(Test ~ Crosstype + Promoter) + scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=exp(1)), breaks=c(0,1,2,3,4,5,10,20,40,60,80,100)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(boxp_neg)

h_pos <- subset(h, h$Fluorescence == "DsRed positive")
h_de_pos <- ddply(h_pos, c("Replicate","White_Phenotype","Test","Promoter","Crosstype"), summarise, Number = sum(Number))
h_de_pos$percent <- with(h_de_pos, ave(Number, Replicate, Crosstype, Promoter, Test ,FUN = prop.table) * 100)
h_de_pos <- subset(h_de_pos, h_de_pos$White_Phenotype != "Red Eyes")

boxp_pos <- ggplot(data=h_de_pos, aes(x=White_Phenotype, y=percent, fill=Test, color=Test)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E0E0E0","#606060")) +
  scale_color_manual(values=c("#404040","#000000")) + 
  theme_classic(base_size = 12) +
  facet_grid(Test ~ Crosstype + Promoter) + scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=exp(1)), breaks=c(0,1,2,3,4,5,10,20,40,60,80,100)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, size=0.1)+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(boxp_pos)


