library(MASS)
library(asbio)
library(stats)
library(car)

#--Data--#

Data.2019 <- read.csv("/Users/spenc/Desktop/2019 Data/Data.2019.csv")

climate <- Data.2019[,c(28:47)] 
locational <- Data.2019[,c(25:27)]
precip <- Data.2019[,c(29,30,43:47)]
temp <- Data.2019[,c(28,31:42)]             
phys <- as.data.frame(Data.2019[,c(7:24)])   
ploidy <- as.factor(Data.2019[,4])
subspecies <- as.factor(Data.2019[,3])   
sub.cyto <- as.factor(Data.2019[,5])
Distance <- Data.2019[,48]

#--Variables--#

temp.pca1 <- prcomp(temp, center=TRUE, scale.=TRUE)
temp.PC <- as.data.frame(temp.pca1$x[,1:2])
names(temp.PC) <- c('TPC1', 'TPC2')
temp.PC1 <- as.data.frame(temp.PC[,1])
names(temp.PC1) <- c('TPC1')
temp.PC2 <- as.data.frame(temp.PC[,2])
names(temp.PC2) <- c('TPC2')

precip.pca1 <- prcomp(precip, center = TRUE, scale. = TRUE)
precip.PC <- as.data.frame(precip.pca1$x[,1:2])
names(precip.PC) <- c('PPC1', 'PPC2')
precip.PC1 <- as.data.frame(precip.PC[,1])
names(precip.PC1) <- c('PPC1')
precip.PC2 <- as.data.frame(precip.PC[,2])
names(precip.PC2) <- c('PPC2')

predictors <- cbind(sub.cyto, temp.PC, precip.PC, Distance)

#--Hypothesis Testing--#

#-Transpiration (E)-#
lmodel.E <- lm(phys$E ~ ., data = predictors)
summary(lmodel.E)

#-Photosynthesis (A)-#
lmodel.A <- lm(phys$A ~ ., data = predictors)
summary(lmodel.A)

#-Respiration (R)-#
lmodel.R <- lm(phys$R ~ ., data = predictors)
summary(lmodel.R)

#-A/R-#
lmodel.A.R <- lm(phys$A.R ~ ., data = predictors)
summary(lmodel.A.R)

#-Stomatal Conductance (gsw)-#
lmodel.gsw <- lm(phys$gsw ~ ., data = predictors)
summary(lmodel.gsw)

#-Midday Photosystem II Efficiency (PhiPS2)-#
lmodel.PhiPS2 <- lm(log(phys$PhiPS2) ~ ., data = predictors)
summary(lmodel.PhiPS2)

#-Predawn PSII Efficiency-#
lmodel.Fv.Fm <- lm(phys$Fv.Fm ~ ., data = predictors)
summary(lmodel.Fv.Fm)

#-Water Use Efficiency (WUE)-#
lmodel.WUE <- lm(phys$WUE ~ ., data = predictors)
summary(lmodel.WUE)

#-Intrinsic Water Use Efficiency (WUEi)-#
lmodel.WUEi <- lm(phys$WUEi ~ ., data = predictors)
summary(lmodel.WUEi)

#-Relative Water Content (RWC)-#
lmodel.RWC <- lm(phys$RWC ~ ., data = predictors)
summary(lmodel.RWC)

#-Specific Leaf Area (SLA)-#
lmodel.SLA <- lm(phys$SLA ~ ., data = predictors)
summary(lmodel.SLA)

#-Total Leaf Area/Specific Leaf Area-#
lmodel.LA.SA <- lm(phys$LA.SA ~ ., data = predictors)
summary(lmodel.LA.SA)

#-Maximum Branch Hydraulic Conductivity (Kbm)-#
lmodel.Kbm <- lm(phys$Kbm ~ ., data = predictors)
summary(lmodel.Kbm)

#-Sapwood Specific Maximum Branch Hydraulic Conductivity (Kbms)-#
lmodel.Kbms <- lm(phys$Kbms ~ ., data = predictors)
summary(lmodel.Kbms)

#-Leaf Area Specific Maximum Branch Hydraulic Conductivity (Kbml)-# 
lmodel.Kbml <- lm(phys$Kbml ~ ., data = predictors)
summary(lmodel.Kbml)

#-Predawn Water Potential (Ypd)-# 
lmodel.Ypd <- lm(phys$Ypd ~ ., data = predictors)
summary(lmodel.Ypd)

#-Midday Water Potential (Ymd)-#
lmodel.Ymd <- lm(phys$Ymd ~ ., data = predictors)
summary(lmodel.Ymd)

#-Maximum Shrub Height-#
lmodel.Max.Height <- lm(phys$Max.Height ~ ., data = predictors)
summary(lmodel.Max.Height)


#--Posthoc Pairwise Analysis--#

A.pair <- pairw.anova(phys$A[!is.na(phys$A)], sub.cyto[!is.na(phys$A)], MSE = 8.820, df.err = 77)
A.pair

PhiPS2.pair <- pairw.anova(phys$PhiPS2[!is.na(phys$PhiPS2)], sub.cyto[!is.na(phys$PhiPS2)], MSE = 0.08258, df.err = 72)
PhiPS2.pair

LASA.pair <- pairw.anova(phys$LA.SA[!is.na(phys$LA.SA)], sub.cyto[!is.na(phys$LA.SA)], MSE = 5355.7, df.err = 80)
LASA.pair

Kbm.pair <- pairw.anova(phys$Kbm[!is.na(phys$Kbm)], sub.cyto[!is.na(phys$Kbm)], MSE = 4407.3, df.err = 81)
Kbm.pair

Kbms.pair <- pairw.anova(phys$Kbms[!is.na(phys$Kbms)], sub.cyto[!is.na(phys$Kbms)], MSE = 485468, df.err = 80)
Kbms.pair

Kbml.pair <- pairw.anova(phys$Kbml[!is.na(phys$Kbml)], sub.cyto[!is.na(phys$Kbml)], MSE = 0.0054281, df.err = 82)
Kbml.pair

Ymd.pair <- pairw.anova(phys$Ymd[!is.na(phys$Ymd)], sub.cyto[!is.na(phys$Ymd)], MSE = 0.325, df.err = 80)
Ymd.pair

MH.pair <- pairw.anova(phys$Max.Height[!is.na(phys$Max.Height)], sub.cyto[!is.na(phys$Max.Height)], MSE = 347.8, df.err = 81)
MH.pair

#--Model Selection--#


#-Transpiration (E)-# 
steplmodel.E <- lm(phys$E ~ ., data = predictors)
stepAIC(steplmodel.E, direction = "both")

#-Photosystem (A)-# 
steplmodel.A <- lm(phys$A ~ ., data = predictors)
stepAIC(steplmodel.A, direction = "both")

#-Respiration (R)-#
steplmodel.R <- lm(phys$R ~ ., data = predictors)
stepAIC(steplmodel.R, direction = "both")

#-A/R-# 
steplmodel.A.R <- lm(phys$A.R ~ ., data = predictors)
stepAIC(steplmodel.A.R, direction = "both")

#-Stomatal Conductance (gsw)-#
steplmodel.gsw <- lm(phys$gsw ~ ., data = predictors)
stepAIC(steplmodel.gsw, direction = "both")

#-Midday PhiPS2-# 
steplmodel.PhiPS2 <- lm(phys$PhiPS2 ~ ., data = predictors)
stepAIC(steplmodel.PhiPS2, direction = "both")

#-Predawn PhiPS2-# 
steplmodel.FvFm <- lm(phys$Fv.Fm ~ ., data = predictors)
stepAIC(steplmodel.FvFm, direction = "both")

#-Water Use Efficiency (WUE)-#
steplmodel.WUE <- lm(phys$WUE ~ ., data = predictors)
stepAIC(steplmodel.WUE, direction = "both")

#-Intrinsic Water Use Efficiency (WUEi)-#
steplmodel.WUEi <- lm(phys$WUEi ~ ., data = predictors)
stepAIC(steplmodel.WUEi, direction = "both")

#-Relative Water Content (RWC)-#
steplmodel.RWC <- lm(phys$RWC ~ ., data = predictors)
stepAIC(steplmodel.RWC, direction = "both")

#-Specific Leaf Area (SLA)-#
steplmodel.SLA <- lm(phys$SLA ~ ., data = predictors)
stepAIC(steplmodel.SLA, direction = "both")

#-Total Leaf Area/Specific Leaf Area-#
steplmodel.LA.SA <- lm(phys$LA.SA ~ ., data = predictors)
stepAIC(steplmodel.LA.SA, direction = "both")

#-Maximum Branch Hydraulic Conductivity (Kbm)-#
steplmodel.Kbm <- lm(phys$Kbm ~ ., data = predictors)
stepAIC(steplmodel.Kbm, direction = "both")

#-Sapwood Specific Maximum Branch Hydraulic Conductivity (Kbms)-#
steplmodel.Kbms <- lm(phys$Kbms ~ ., data = predictors)
stepAIC(steplmodel.Kbms, direction = "both")

#-Leaf Area Specific Maximum Branch Hydraulic Conductivity (Kbml)-#
steplmodel.Kbml <- lm(phys$Kbml ~ ., data = predictors)
stepAIC(steplmodel.Kbml, direction = "both")

#-Predawn Water Potential-#
steplmodel.Ypd <- lm(phys$Ypd ~ ., data = predictors)
stepAIC(steplmodel.Ypd, direction = "both")

#-Midday Water Potential-#
steplmodel.Ymd <- lm(phys$Ymd ~ ., data = predictors)
stepAIC(steplmodel.Ymd)

#-Maximum Height-# 
steplmodel.Max.Height <- lm(phys$Max.Height ~ ., data = predictors)
stepAIC(steplmodel.Max.Height)

#--Graphs--#

library(ggplot2)
library(ggpmisc)
library(cowplot)
library(gridExtra)
library(devtools)

#-Subspecies:Cytotype-#

Apairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.007, 0.058, 0.959))
boxplot.A <- ggplot(Data.2019, aes(sub.cyto, phys$A, color = sub.cyto)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4, notch = FALSE) +
  labs(x = "", 
       y = expression("A (umol m"^-2*"s"^-1*")")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 2, y = 15, label = "B") + 
  annotate("text", x = 3, y = 15, label = "AB") + 
  annotate("text", x = 1, y = 15, label = "A") + 
  annotate(geom = 'table', 
           x = 4,
           y = 15, 
           label = list(Apairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(A)")

PhiPS2pairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.910, 0.838, 0.984))
boxplot.PhiPS2 <- ggplot(Data.2019, aes(sub.cyto, phys$PhiPS2, color = sub.cyto)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "",
       y = expression(paste("Midday", phi, "PSII"))) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 0.4, label = "A") +
  annotate("text", x = 2, y = 0.4, label = "A") + 
  annotate("text", x = 3, y = 0.4, label = "A") + 
  annotate(geom = 'table', 
           x = 4,
           y = 0.4, 
           label = list(PhiPS2pairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(B)")

LASApairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.232, 0.059, 0.685))
boxplot.LASA <- ggplot(Data.2019, aes(sub.cyto, phys$LA.SA, color = sub.cyto)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "",
       y = expression(paste("Leaf Area/Sapwood Area (cm"^2*"/mm"^2*")"))) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 400, label = "A") + 
  annotate("text", x = 2, y = 400, label = "A") + 
  annotate("text", x = 3, y = 400, label = "A") + 
  annotate(geom = 'table', 
           x = 4,
           y = 400, 
           label = list(LASApairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(C)")

Kbmpairdata <- data.frame(pair = c("T2-T4","T3-W4","T4-W4"), p = c(0.374, 0.879, 0.322))
boxplot.Kbm <- ggplot(Data.2019, aes(sub.cyto, phys$Kbm, color = sub.cyto)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "",
       y = expression(paste("K"["max"]*" (mol m"^-1* "s" ^-1* "mPa" ^-1* ")"))) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 300, label = "A") + 
  annotate("text", x = 2, y = 300, label = "A") + 
  annotate("text", x = 3, y = 300, label = "A") +
  annotate(geom = 'table',
           x = 4,
           y = 300,
           label = list(Kbmpairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(D)")

Kbmspairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.933, 0.450, 0.693))
boxplot.Kbms <- ggplot(Data.2019, aes(sub.cyto, phys$Kbms, color = sub.cyto)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "",
       y = expression(paste("Sapwood Area Specific K"["max"]*" (mol m"^-1* "s"^-1* "mPa"^-1*")"))) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 6000, label = "A") + 
  annotate("text", x = 2, y = 6000, label = "A") + 
  annotate("text", x = 3, y = 6000, label = "A") + 
  annotate("table",
           x = 4, 
           y = 6000, 
           label = list(Kbmspairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(F)")

Kbmlpairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.031, 0.305, 0.811))
boxplot.Kbml <- ggplot(Data.2019, aes(sub.cyto, phys$Kbml, color = sub.cyto)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "",
       y = expression(paste("Leaf Area Specific K"["max"]*" (mol m"^-1* "s"^-1* "mPa"^-1*")"))) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 0.55, label = "A") + 
  annotate("text", x = 2, y = 0.55, label = "B") + 
  annotate("text", x = 3, y = 0.55, label = "AB") + 
  annotate("table",
           x = 4, 
           y = 0.55, 
           label = list(Kbmlpairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(E)")

Ymdpairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.993, 0.509 ,0.640))
boxplot.Ymd <- ggplot(Data.2019, aes(sub.cyto, phys$Ymd, color = sub.cyto)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "Subspecies:Cytotype", 
       y = expression(paste(psi, " "["m"]*" (MPa)"))) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 1, label = "A") + 
  annotate("text", x = 2, y = 1, label = "A") + 
  annotate("text", x = 3, y = 1, label = "A") + 
  annotate("table",
           x = 4, 
           y = 1, 
           label = list(Ymdpairdata)) + 
  annotate("text", x = 3.75, y = -5, label = "(G)")

MHpairdata <- data.frame(pair = c("T2-T4","T2-W4","T4-W4"), p = c(0.022, "<0.000",0.0002))
boxplot.MaxHeight <- ggplot(Data.2019, aes(sub.cyto, phys$Max.Height, color = sub.cyto)) + 
  geom_boxplot(outlier.color = "red", outlier.shape = 8, outlier.size = 4, notch = FALSE) + 
  labs(x = "Subspecies:Cytotype", 
       y = "Maximum Shrub Height (cm)") + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
  theme_classic() + 
  scale_color_grey() + 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 175, label = "A") + 
  annotate("text", x = 2, y = 175, label = "B") + 
  annotate("text", x = 3, y = 175, label = "C") + 
  annotate("table",
           x = 4, 
           y = 175, 
           label = list(MHpairdata)) + 
  annotate("text", x = 3.75, y = 0, label = "(H)")

SubCyto.Graphs <- plot_grid(boxplot.A, 
                            boxplot.PhiPS2,
                            boxplot.LASA,
                            boxplot.Kbm,
                            boxplot.Kbml,
                            boxplot.Kbms,
                            boxplot.Ymd,
                            boxplot.MaxHeight,
                            ncol = 2, nrow = 4)

png("/Users/spenc/Desktop/SubCyto.Graphs.png", width = 15, height = 20, units = "in", res = 1200, pointsize = 4)
SubCyto.Graphs
dev.off()

#-Temperature Principal Components-#

TPC1.AR.Graph <- ggplot(Data.2019, aes(predictors$TPC1, phys$A.R)) + 
  geom_point() + 
  geom_segment(aes(x = -5.118, 
                   xend = 6.534, 
                   y = coef(lmodel.A.R)[1] + coef(lmodel.A.R)[4],
                   yend = coef(lmodel.A.R)[1] + coef(lmodel.A.R)[4]*6.534),
               colour = "#E41A1C") + 
  labs(x = "",
       y = expression("A:R (umol m"^2*"s"^-1*"/ umol m"^2*"s"^-1*")")) + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 1, y = 1.70, label = "p = 0.051", size = 5) + 
  annotate("text", x = 6, y = 0.1, label = "(A)")

TPC2.R.Graph <- ggplot(Data.2019, aes(predictors$TPC2, phys$R)) + 
  geom_point() + 
  geom_segment(aes(x = -3.355, 
                   xend = 3.181, 
                   y = coef(l.model.R)[1] + coef(l.model.R)[2], 
                   yend = coef(l.model.R)[1] + coef(l.model.R)[2]*3.181), 
               colour = "#E41A1C") +
  labs(x = "", 
       y = expression("R (umol m"^-2*"s"^-1*")")) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 2.5, y = 6, label = "p = 0.038", size = 5) + 
  annotate("text", x = 3, y = 4, label = "(B)")

TPC1.Ymd.Graph <- ggplot(Data.2019, aes(predictors$TPC1, phys$Ymd)) + 
  geom_point() +   
  geom_segment(aes(x = -5.118, 
                   xend = 6.534, 
                   y = coef(l.model.Ymd)[1] + coef(l.model.Ymd)[2], 
                   yend = coef(l.model.Ymd)[1] + coef(l.model.Ymd)[2]*6.534), 
               colour = "#E41A1C") +
  labs(x = "Temperature Principal Component 1", 
       y = expression(paste(psi, " "["m"]*" (mPa)"))) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)) + 
  annotate("text", x = 1, y = -3.3, label = "p = 0.035", size = 5) + 
  annotate("text", x = 6, y = -5, label = "(C)")

TPC2.Ymd.Graph <- ggplot(Data.2019, aes(predictors$TPC2, phys$Ymd)) + 
  geom_point() +   
  geom_segment(aes(x = -3.355, 
                   xend = 3.181, 
                   y = coef(l.model.Ymd)[1] + coef(l.model.Ymd)[3], 
                   yend = coef(l.model.Ymd)[1] + coef(l.model.Ymd)[3]*3.181), 
               colour = "#E41A1C") +
  labs(x = "Temperature Principal Component 2", 
       y = expression(paste(psi, " "["m"]*" (mPa)"))) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)) + 
  annotate("text", x = 2.3, y = -2.65, label = "p = 0.055", size = 5) + 
  annotate("text", x = 3, y = -5, label = "(D)")

TPC.Graphs <- plot_grid(TPC1.AR.Graph,
                        TPC2.R.Graph,
                        TPC1.Ymd.Graph,
                        TPC2.Ymd.Graph,
                        ncol = 2, nrow = 2)

png("/Users/spenc/Desktop/TPC.Graphs.png", width = 20, height = 15, units = "in", res = 1200, pointsize = 4)
TPC.Graphs
dev.off()

#-Precipitation Principal Component-#

PPC1.AR.Graph <- ggplot(Data.2019, aes(predictors$PPC1, phys$A.R)) + 
  geom_point() + 
  geom_segment(aes(x = -3.436, 
                   xend = 2.161, 
                   y = coef(l.model.AR)[1] + coef(l.model.AR)[3], 
                   yend = coef(l.model.AR)[1] + coef(l.model.AR)[3]*2.161), 
               colour = "#E41A1C") +
  labs(x = "", 
       y = expression("A:R (umol m"^2*"s"^-1*"/ umol m"^2*"s"^-1*")")) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = -1.8, y = 1.1, label = "p = 0.073", size = 5) + 
  annotate("text", x = 2.5, y = 0, label = "(A)")

PPC2.WUEi.Graph <- ggplot(Data.2019, aes(predictors$PPC2, phys$WUEi)) + 
  geom_point() + 
  geom_segment(aes(x = -2.77,
                   xend = 4.304, 
                   y = coef(l.model.WUEi)[1] + coef(l.model.WUEi)[3],
                   yend = coef(l.model.WUEi)[1] + coef(l.model.WUEi)[3]*4.304),
               colour = "#E41A1C") +
  labs(x = "", 
       y = expression("WUEi (umol m"^-2*"s"^-1*"/mol m"^-2*"s"^-1*")")) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 3, y = 175, label = "p = 0.035", size = 5) + 
  annotate("text", x = 4.5, y = 0, label = "(B)")


PPC1.PhiPS2.Graph <- ggplot(Data.2019, aes(predictors$PPC1, phys$PhiPS2)) + 
  geom_point() + 
  geom_segment(aes(x = -3.436, 
                   xend = 2.161, 
                   y = coef(l.model.PhiPS2)[1] + coef(l.model.PhiPS2)[5], 
                   yend = coef(l.model.PhiPS2)[1] + coef(l.model.PhiPS2)[5]*2.161), 
               colour = "#E41A1C") +
  labs(x = "", 
       y = expression(paste("Midday", phi, "PSII"))) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = -2, y = 0.28, label = "p = 0.049", size = 5) + 
  annotate("text", x = 2.5, y = 0, label = "(C)")

PPC2.Ymd.Graph <- ggplot(Data.2019, aes(predictors$PPC2, phys$Ymd)) + 
  geom_point() +   
  geom_segment(aes(x = -2.77,
                   xend = 4.304, 
                   y = coef(l.model.Ymd)[1] + coef(l.model.Ymd)[4],
                   yend = coef(l.model.Ymd)[1] + coef(l.model.Ymd)[4]*4.304),
               colour = "#E41A1C") +
  labs(x = "", 
       y = expression(paste(psi, " "["m"]*" (mPa)"))) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15)) + 
  annotate("text", x = 3.9, y = -2.6, label = "p = 0.041") + 
  annotate("text", x = 4.5, y = -5, label = "(D)")

PPC1.Kbm.Graph <- ggplot(Data.2019, aes(predictors$PPC1, phys$Kbm)) + 
  geom_point() +   
  geom_segment(aes(x = -3.436, 
                   xend = 2.161, 
                   y = coef(l.model.Kbm)[1] + coef(l.model.Kbm)[4], 
                   yend = coef(l.model.Kbm)[1] + coef(l.model.Kbm)[4]*2.161), 
               colour = "#E41A1C") +
  labs(x = "Precipitation Principal Component 1", 
       y = expression(paste("K"["max"]*" (mol m"^-1* "s" ^-1* "mPa" ^-1* ")"))) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)) + 
  annotate("text", x = -2, y = 88, label = "p = 0.055", size = 5) + 
  annotate("text", x = 2.5, y = 0, label = "(E)")

PPC2.MH.Graph <- ggplot(Data.2019, aes(predictors$PPC2, phys$Max.Height)) + 
  geom_point() + 
  geom_segment(aes(x = -2.77,
                   xend = 4.304, 
                   y = coef(l.model.MH)[1] + coef(l.model.MH)[4],
                   yend = coef(l.model.MH)[1] + coef(l.model.MH)[4]*4.304),
               colour = "#E41A1C") +
  labs(x = "Precipitation Principal Component 2", 
       y = expression("Maximum Shrub Height (cm"^2*")")) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)) + 
  annotate("text", x = 2.75, y = 112, label = "p = 0.032", size = 5) + 
  annotate("text", x = 4.5, y = 0, label = "(F)")

PPC.Graphs <- plot_grid(PPC1.AR.Graph,
                        PPC2.WUEi.Graph,
                        PPC1.PhiPS2.Graph,
                        PPC2.Ymd.Graph,
                        PPC1.Kbm.Graph,
                        PPC2.MH.Graph,
                        ncol = 2, nrow = 3)

png("/Users/spenc/Desktop/PPC.Graphs.png", width = 20, height = 15, units = "in", res = 1200, pointsize = 4)
PPC.Graphs
dev.off()

#-Euclidian Distance-#

Dist.PhiPS2.Graph <- ggplot(Data.2019, aes(predictors$Distance, phys$PhiPS2)) + 
  geom_point() + 
  geom_segment(aes(x = 0, 
                   xend = 2002.649, 
                   y = coef(l.model.PhiPS2)[1] + coef(l.model.PhiPS2)[6], 
                   yend = coef(l.model.PhiPS2)[1] + coef(l.model.PhiPS2)[6]*2002.649), 
               colour = "#E41A1C") +
  labs(x = "Euclidian Climate Distance", 
       y = expression(paste("Midday", phi, "PSII"))) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)) + 
  annotate("text", x = 400, y = 0.25, label = "p = 0.073", size = 5)

png("/Users/spenc/Desktop/Dist.Graphs.png", width = 10, height = 7, units = "in", res = 1200, pointsize = 4)
Dist.PhiPS2.Graph
dev.off()