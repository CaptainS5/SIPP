library(ggplot2)
library(ez)
# library(Hmisc)
# library(reshape2)
library(psychReport)
# library(lsr)
# library(bayestestR)
# library(BayesFactor)
# library(TOSTER)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R")
source("pairwise.t.test.with.t.and.df.R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R/plots/")
### modify these parameters to plot different conditions
# varName <- "response"
dataFileName <- "summaryDataDiffSubgroup.csv"
# pdfFileName <- paste("diff_", varName, "_all.pdf", sep = "")
## for catch-up saccades
# varName <- "numY"
# dataFileName <- "summaryCatchUpSaccades.csv"
# pdfFileName <- paste("saccadeBias_consistency", varName, "_all.pdf", sep = "")

# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 3

data <- read.csv(dataFileName)
subAll <- unique(data$sub)
subTotalN <- length(subAll)

#### catch-up saccade analysis
sub <- data["sub"]
rdkApertureAngle <- data["rdkApertureAngle"]
sacGroup <- data["sacGroup"]
measure <- data[varName]
measure[is.na(measure)] = 0
dataAnova <- data.frame(sub, sacGroup, rdkApertureAngle, measure)
dataAnova$sacGroup <- as.factor(dataAnova$sacGroup)
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$rdkApertureAngle <- as.factor(dataAnova$rdkApertureAngle)
colnames(dataAnova)[4] <- "measure"

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(sacGroup, rdkApertureAngle), type = 3, return_aov = TRUE, detailed = TRUE)
# print(aonvaData)
aovEffectSize(anovaData, 'pes')

dataPlot <- data.frame(sub, sacGroup, rdkApertureAngle, measure)
colnames(dataPlot)[4] <- "measure"
dataPlot$sub <- as.factor(dataPlot$sub)
dataPlot$sacGroup <- as.factor(dataPlot$sacGroup)
# saccades
ylimLow <- -1
ylimHigh <- 1

p <- ggplot(dataPlot, aes(x = rdkApertureAngle, y = measure, color = sacGroup)) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun = mean, geom = "line", width = 1) +
        stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subTotalN)), geom = 'errorbar', width = 1, size = 1) +
        stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_point(size = dotSize, shape = 1) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-9, -12), y = c(ylimLow, ylimLow), xend = c(9, -12), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
        # scale_y_continuous(name = "Bias in sum vertical amplitudes (deg)", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        scale_y_continuous(name = "Bias in number", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = varName, limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) + 
        scale_x_continuous(name = "Object motion direction (deg)", breaks=c(-9,-6,-3,0,3,6,9), limits = c(-12, 12), expand = c(0, 0)) +
        # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "90")) +
        scale_colour_manual(name = "Saccade direction", labels = c("same", "opposite"), values = c("#DBB167", "#2ADB77")) +
        theme(axis.text=element_text(colour="black"),
                      axis.ticks=element_line(colour="black", size = axisLineWidth),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size = textSize, colour = "black"),
                      legend.background = element_rect(fill="transparent"),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      legend.position="top")
        # facet_wrap(~exp)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))


#### compare different dot motion conditions
sub <- data["sub"]
rdkApertureAngle <- data["rdkApertureAngle"]
rdkInternalDir <- data["rdkInternalDir"]
measure <- data[varName]
measure[is.na(measure)] = 0
dataAnova <- data.frame(sub, rdkInternalDir, rdkApertureAngle, measure)
dataAnova$rdkInternalDir <- as.factor(dataAnova$rdkInternalDir)
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$rdkApertureAngle <- as.factor(dataAnova$rdkApertureAngle)
# dataAnova$timeBin <- as.factor(dataAnova$timeBin)
colnames(dataAnova)[4] <- "measure"
# dataAnova <- aggregate(perceptualErrorMean ~ sub * rotationSpeed * exp,
    # data = dataTemp, FUN = "mean")

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(rdkInternalDir, rdkApertureAngle), type = 3, return_aov = TRUE, detailed = TRUE)
# print(aonvaData)
aovEffectSize(anovaData, 'pes')

# ## t-test of the simple main effect of probability in the control experiment
# dataD <- dataAnova[dataAnova$exp==3,]
# # show(dataD)
# res <- pairwise.t.test.with.t.and.df(x = dataD$measure, g = dataD$prob, paired = TRUE, p.adj="none")
# show(res) # [[3]] = p value table, un adjusted
# res[[5]] # t-value
# res[[6]] # dfs
# res[[3]]
# p.adjust(res[[3]], method = "bonferroni", n = 4) 
# cohensd <- cohensD(subset(dataD, prob==50)$measure, subset(dataD, prob==90)$measure, method = 'paired')
# show(cohensd)

## interaction plot
dataPlot <- data.frame(sub, rdkApertureAngle, rdkInternalDir, measure)
colnames(dataPlot)[4] <- "measure"
dataPlot$sub <- as.factor(dataPlot$sub)
dataPlot$rdkInternalDir <- as.factor(dataPlot$rdkInternalDir)

# # pursuit bias clp
# ylimLow <- -15
# ylimHigh <- 15
# response
ylimLow <- -15
ylimHigh <- 15

p <- ggplot(dataPlot, aes(x = rdkApertureAngle, y = measure, color = rdkInternalDir)) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun = mean, geom = "line", width = 1) +
        stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subTotalN)), geom = 'errorbar', width = 1, size = 1) +
        stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_point(size = dotSize, shape = 1) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-9, -12), y = c(ylimLow, ylimLow), xend = c(9, -12), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
        scale_y_continuous(name = "Bias in perceived direction (deg)", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = "Bias in clp pursuit direction (deg)", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = varName, limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) + 
        scale_x_continuous(name = "Object motion direction (deg)", breaks=c(-9,-6,-3,0,3,6,9), limits = c(-12, 12), expand = c(0, 0)) +
        theme(axis.text=element_text(colour="black"),
                      axis.ticks=element_line(colour="black", size = axisLineWidth),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size = textSize, colour = "black"),
                      legend.background = element_rect(fill="transparent"),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      legend.position="top")
        # facet_wrap(~exp)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))

#### Comparison between perceptual bias groups
### separetly plot each group
pdfFileName <- paste("diff_", varName, "_nobiasGroup.pdf", sep = "")

dataSub <- subset(data, group==3)
subAll <- unique(dataSub["sub"])
subTotalN <- dim(subAll)[1]

sub <- dataSub["sub"]
rdkApertureAngle <- dataSub["rdkApertureAngle"]
rdkInternalDir <- dataSub["rdkInternalDir"]
measure <- dataSub[varName]
measure[is.na(measure)] = 0
dataAnova <- data.frame(sub, rdkInternalDir, rdkApertureAngle, measure)
dataAnova$rdkInternalDir <- as.factor(dataAnova$rdkInternalDir)
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$rdkApertureAngle <- as.factor(dataAnova$rdkApertureAngle)
# dataAnova$timeBin <- as.factor(dataAnova$timeBin)
colnames(dataAnova)[4] <- "measure"
# dataAnova <- aggregate(perceptualErrorMean ~ sub * rotationSpeed * exp,
    # data = dataTemp, FUN = "mean")

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(rdkInternalDir, rdkApertureAngle), type = 3, return_aov = TRUE, detailed = TRUE)
# print(aonvaData)
aovEffectSize(anovaData, 'pes')

### plot
dataPlot <- data.frame(sub, rdkApertureAngle, rdkInternalDir, measure)
colnames(dataPlot)[4] <- "measure"
dataPlot$sub <- as.factor(dataPlot$sub)
dataPlot$rdkInternalDir <- as.factor(dataPlot$rdkInternalDir)

# response
ylimLow <- -15
ylimHigh <- 15

p <- ggplot(dataPlot, aes(x = rdkApertureAngle, y = measure, color = rdkInternalDir)) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun = mean, geom = "line", width = 1) +
        stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subTotalN)), geom = 'errorbar', width = 1, size = 1) +
        # stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_point(size = dotSize, shape = 1) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-9, -12), y = c(ylimLow, ylimLow), xend = c(9, -12), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
        scale_y_continuous(name = "Bias in perceived direction (deg)", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = varName, limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) + 
        scale_x_continuous(name = "Object motion direction (deg)", breaks=c(-9,-6,-3,0,3,6,9), limits = c(-12, 12), expand = c(0, 0)) +
        # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "90")) +
        # scale_colour_discrete(name = "After reversal\ndirection", labels = c("CCW", "CW")) +
        theme(axis.text=element_text(colour="black"),
                      axis.ticks=element_line(colour="black", size = axisLineWidth),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size = textSize, colour = "black"),
                      legend.background = element_rect(fill="transparent"),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      legend.position="top")
        # facet_wrap(~exp)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))

### compare with temporal dynamics
pdfFileName <- paste("diffPursuitTemporal_subgroup.pdf", sep = "")
varNames <- c("rdkInternalDir", "group", "dirOlp", "dirClpEarly", "dirClpLate")
dataAnova <- data.frame(sub=double(60),
                 pursuitPhase=double(60),
                 group=double(60),
                 measure=double(60),
                 stringsAsFactors=FALSE)
rowN <- 1
for (subN in subAll) {
    dataT <- subset(data, sub==subN, select = varNames)
    dataT[which(dataT["rdkInternalDir"]<0), 3:5] <- -dataT[which(dataT["rdkInternalDir"]<0), 3:5]
    for (phaseN in 1:3) {
        dataAnova$sub[rowN] <- subN
        dataAnova$group[rowN] <- dataT$group[1]
        dataAnova$pursuitPhase[rowN] <- phaseN
        dataAnova$measure[rowN] <- mean(dataT[, phaseN+2])
        rowN <- rowN+1
    }
}
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$group <- as.factor(dataAnova$group)
dataAnova$pursuitPhase <- as.factor(dataAnova$pursuitPhase)

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(pursuitPhase), between = .(group), type = 3, return_aov = TRUE, detailed = TRUE)
# print(aonvaData)
aovEffectSize(anovaData, 'pes')

### plot
dataPlot <- dataAnova

# response
ylimLow <- -5
ylimHigh <- 15

p <- ggplot(dataPlot, aes(x = pursuitPhase, y = measure, color = group)) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5, position = position_dodge(width = 0.75)) +
        # stat_summary(fun = mean, geom = "line", width = 1, position = position_dodge()) +
        stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subTotalN)), geom = 'errorbar', width = 0.1, size = 1, position = position_dodge(width = 0.75)) +
        # stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_point(size = dotSize, shape = 1, position = position_jitterdodge()) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(1, 0.5), y = c(ylimLow, ylimLow), xend = c(3, 0.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
        scale_y_continuous(name = "Mean pursuit bias (deg)", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        # scale_x_continuous(name = "Object motion direction (deg)", breaks=c(-9,-6,-3,0,3,6,9), limits = c(-12, 12), expand = c(0, 0)) +
        scale_x_discrete(name = "Pursuit phase", breaks=1:3, labels = c("Initiation", "Early steady-state", "Late steady-state")) +
        scale_colour_manual(name = "Perceptual bias group", labels = c("Assimilation", "Contrast", "No-bias"), values=c("#C77CFF", "#7CAE00", "#000000")) +
        theme(axis.text=element_text(colour="black"),
                      axis.ticks=element_line(colour="black", size = axisLineWidth),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size = textSize, colour = "black"),
                      legend.background = element_rect(fill="transparent"),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      legend.position="top")
        # facet_wrap(~exp)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))

### fit individual slope for each person, and compare between groups
library(rmcorr)
library(lme4)
library(lmerTest)

varNames <- c("rdkInternalDir", "group", "dirOlp", "dirClpEarly", "dirClpLate")
dataBase <- data.frame(sub=double(60),
                 pursuitPhase=double(60),
                 group=double(60),
                 measure=double(60),
                 stringsAsFactors=FALSE)
rowN <- 1
for (subN in subAll) {
    dataT <- subset(data, sub==subN, select = varNames)
    dataT[which(dataT["rdkInternalDir"]<0), 3:5] <- -dataT[which(dataT["rdkInternalDir"]<0), 3:5]
    for (phaseN in 1:3) {
        dataBase$sub[rowN] <- subN
        dataBase$group[rowN] <- dataT$group[1]
        dataBase$pursuitPhase[rowN] <- phaseN
        dataBase$measure[rowN] <- mean(dataT[, phaseN+2])
        rowN <- rowN+1
    }
}
dataBase$sub <- as.factor(dataBase$sub)
dataBase$group <- as.factor(dataBase$group)

slopeData <- data.frame(sub=double(20),
                 group=double(20),
                 slope=double(20),
                 stringsAsFactors=FALSE)
for (subN in subAll) {
    dataTemp <- subset(dataBase, sub==subN)
    LM <- lm( measure ~ pursuitPhase, dataTemp)
    slopeData$sub[subN] <- subN
    slopeData$group[subN] <- dataTemp$group[1]
    slopeData$slope[subN] <- coef(LM)["pursuitPhase"]
}
show(slopeData)
slopeData$sub <- as.factor(slopeData$sub)
slopeData$group <- as.factor(slopeData$group)

options(contrasts=c("contr.sum","contr.poly"))
anovaData <- ezANOVA(slopeData, dv = .(slope), wid = .(sub),
    between = .(group), type = 3)
print(anovaData)

## post-hoc t-test
aovData <- aov(slope ~ group, slopeData)
res <- TukeyHSD(aovData)
print(res)

model = lm(slope ~ group, slopeData)
library(emmeans)
marginal = emmeans(model, ~ group)
pairs(marginal, adjust="tukey")


pdfFileName <- paste("slopePursuitBiasTemporal_subgroup.pdf", sep = "")
dataPlot <- slopeData
# # pursuit bias magnitude
# ylimLow <- -2
# ylimHigh <- 15
# average response
ylimLow <- -7
ylimHigh <- -2

p <- ggplot(data=dataPlot, aes(x = group, y = slope)) +
        stat_summary(aes(y = slope), fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1.96/sqrt(subTotalN)),
               geom = 'linerange', size = 1) +
        # geom_line(aes(x = group, y = slope, group = sub), size = 0.5, linetype = "dashed") +
        geom_point(aes(x = group, y = slope), size = dotSize, shape = 1, position = position_jitter(w = 0.1, h = 0)) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(1, 0.5), y = c(ylimLow, ylimLow), xend = c(3, 0.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
        scale_y_continuous(name = "Slope of the temporal change in pursuit bias", breaks = c(ylimLow, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = "Time of the start of the window after RDK onset (ms)", breaks = c(ylimLow, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        # scale_x_continuous(name = "Preceded perceived direction", breaks=c(0, 1), limits = c(-0.5, 1.5), expand = c(0, 0)) + # preceded perception
        scale_x_discrete(name = "Perceptual bias group", breaks=1:3, labels = c("Assimilation", "Contrast", "No-bias")) + 
        # scale_colour_discrete(name = "After reversal\ndirection", labels = c("CCW", "CW")) +
        theme(axis.text=element_text(colour="black"),
              axis.ticks=element_line(colour="black", size = axisLineWidth),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              text = element_text(size = textSize, colour = "black"),
              legend.background = element_rect(fill="transparent"),
              legend.key = element_rect(colour = "transparent", fill = "white"))
        # facet_wrap(~prob)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))