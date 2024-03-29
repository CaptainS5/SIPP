library(ggplot2)
library(ez)
# library(Hmisc)
# library(lsr)
# library(BayesFactor)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R")
source("pairwise.t.test.with.t.and.df.R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R/plots/")

### modify these parameters to plot different conditions
varName <- "response"
dataFileName <- "summaryData.csv"

# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 3

# source("pairwise.t.test.with.t.and.df.R")
data <- read.csv(dataFileName)
subAll <- unique(data$sub)
subTotalN <- length(subAll)



#### plot baseline perception response
pdfFileName <- paste("perceptionBias_baseline_all.pdf", sep = "")
dataB <- subset(data, rdkInternalDir==0)

sub <- dataB["sub"]
rdkApertureAngle <- dataB["rdkApertureAngle"]
measure <- dataB[varName]
dataPlot <- data.frame(sub, rdkApertureAngle, measure)
dataPlot$sub <- as.factor(dataPlot$sub)
colnames(dataPlot)[3] <- "measure"

# plotting
ylimLow <- -25
ylimHigh <- 25

p <- ggplot(data=dataPlot, aes(x = rdkApertureAngle, y = measure)) +
        stat_summary(aes(y = measure), fun = mean, geom = "point", shape = 95, size = 17.5) +
        # stat_summary(fun = mean, geom = "line", width = 1) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1.96/sqrt(subTotalN)),
               geom = 'errorbar', width = 1, size = 0.75) +
        geom_point(size = dotSize, shape = 1) +
        geom_segment(aes(x = -9, xend = 9, y = -9, yend = 9)) + # add the identity line
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-9, -12), y = c(ylimLow, ylimLow), xend = c(9, -12), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
        scale_y_continuous(name = "Bias in perceived direction (deg)", limits = c(ylimLow, ylimHigh), breaks = seq(from = ylimLow, to = ylimHigh, by = 5), expand = c(0, 0)) +
        scale_x_continuous(name = "Object motion direction (deg)", breaks=c(-9,-6,-3,0,3,6,9), limits = c(-12, 12), expand = c(0, 0)) + 
        # scale_colour_discrete(name = "Internal motion", labels = c("Up", "Down"), values = c("#67a9cf", "ef8a62")) +
        theme(axis.text=element_text(colour="black"),
              axis.ticks=element_line(colour="black", size = axisLineWidth),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              text = element_text(size = textSize, colour = "black"),
              legend.background = element_rect(fill="transparent"),
              legend.key = element_rect(colour = "transparent", fill = "white"),
              aspect.ratio = 0.85)
        # facet_wrap(~prob)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))




#### process temporal dynamics in pursuit: three variables
pdfFileName <- paste("pursuitBias_differentPhases_all.pdf", sep = "")
### get the data ready
varNames <- c("rdkInternalDir", "dirOlp", "dirClpEarly", "dirClpLate")
dataAnova <- data.frame(sub=double(60),
                 pursuitPhase=double(60),
                 measure=double(60),
                 stringsAsFactors=FALSE)
rowN <- 1
for (subN in subAll) {
    dataT <- subset(data, sub==subN, select = varNames)
    dataT[which(dataT["rdkInternalDir"]<0), 2:4] <- -dataT[which(dataT["rdkInternalDir"]<0), 2:4]
    for (phaseN in 1:3) {
        dataAnova$sub[rowN] <- subN
        dataAnova$pursuitPhase[rowN] <- phaseN
        dataAnova$measure[rowN] <- mean(dataT[, phaseN+1])
        rowN <- rowN+1
    }
}
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$pursuitPhase <- as.factor(dataAnova$pursuitPhase)

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(pursuitPhase), type = 3)
print(anovaData)

### plot
# pursuit bias
ylimLow <- -2
ylimHigh <- 15

p <- ggplot(data=dataAnova, aes(x = pursuitPhase, y = measure)) +
        stat_summary(aes(y = measure), fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1.96/sqrt(subTotalN)),
               geom = 'linerange', size = 1) +
        # geom_line(aes(x = group, y = measure, group = sub), size = 0.5, linetype = "dashed") +
        geom_point(aes(x = pursuitPhase, y = measure), size = dotSize, shape = 1, position = position_jitter(w = 0.1, h = 0)) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(1, 0.5), y = c(ylimLow, ylimLow), xend = c(3, 0.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
        scale_y_continuous(name = "Bias in pursuit direction (deg))", breaks = c(ylimLow, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = "Time of the start of the window after RDK onset (ms)", breaks = c(ylimLow, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        # scale_x_continuous(name = "Preceded perceived direction", breaks=c(0, 1), limits = c(-0.5, 1.5), expand = c(0, 0)) + # preceded perception
        scale_x_discrete(name = "Pursuit phase", breaks=1:3, labels = c("Initiation", "Early steady-state", "Late steady-state")) + 
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


#### Comparison between perceptual bias groups
### about the temporal relationship
varName <- c("wMiddle")
# pdfFileName <- paste("groupCompTemporal_", varName, ".pdf", sep = "")

sub <- data["sub"]
group <- data["group"]
measure <- data["wStart"]+data["wLength"]/2
# measure <- data[varName]
dataAnova <- data.frame(sub, group, measure)
dataAnova$group <- as.factor(dataAnova$group)
dataAnova$sub <- as.factor(dataAnova$sub)
colnames(dataAnova)[3] <- "measure"
# dataAnova <- dataAnova[-10, ]

options(contrasts=c("contr.sum","contr.poly"))
anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    between = .(group), type = 3)
print(anovaData)

## post-hoc t-test
aovData <- aov(measure ~ group, dataAnova)
res <- TukeyHSD(aovData)
print(res)

model = lm(measure ~ group, dataAnova)
library(emmeans)
marginal = emmeans(model, ~ group)
pairs(marginal, adjust="tukey")

dataPlot <- dataAnova
# windows
ylimLow <- 200
ylimHigh <- 700
# R2adjusted
# ylimLow <- 0
# ylimHigh <- 1

p <- ggplot(data=dataPlot, aes(x = group, y = measure)) +
        stat_summary(aes(y = measure), fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1.96/sqrt(subTotalN)),
               geom = 'linerange', size = 1) +
        # geom_line(aes(x = group, y = measure, group = sub), size = 0.5, linetype = "dashed") +
        geom_point(aes(x = group, y = measure), size = dotSize, shape = 1, position = position_jitter(w = 0.1, h = 0)) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(1, 0.5), y = c(ylimLow, ylimLow), xend = c(3, 0.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
        scale_y_continuous(name = "Max R^(2) adjusted (ms)", breaks = c(ylimLow, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
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

### need to average across conditions for each participant
pdfFileName <- paste("dirClpBias_subgroups.pdf", sep = "")
### get the data ready
varNames <- c("rdkInternalDir", "group", "dirClp")
dataAnova <- data.frame(sub=double(subTotalN),
                 group=double(subTotalN),
                 measure=double(subTotalN),
                 stringsAsFactors=FALSE)
rowN <- 1
for (subN in subAll) {
    dataT <- subset(data, sub==subN, select = varNames)
    dataT[which(dataT["rdkInternalDir"]<0), 3] <- -dataT[which(dataT["rdkInternalDir"]<0), 3]
    dataAnova$sub[rowN] <- subN
    dataAnova$group[rowN] <- dataT$group[1]
    dataAnova$measure[rowN] <- abs(mean(dataT[, 3]))
    rowN <- rowN+1
}
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$group <- as.factor(dataAnova$group)

options(contrasts=c("contr.sum","contr.poly"))
anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    between = .(group), type = 3)
print(anovaData)

aovData <- aov(measure ~ group, dataAnova)
print(aovData)

## post-hoc t-test
res <- TukeyHSD(aovData)
print(res)

model = lm(measure ~ group, dataAnova)
library(emmeans)
marginal = emmeans(model, ~ group)
pairs(marginal, adjust="tukey")

# dataD <- dataAnova[dataAnova$group!=3,]
# show(dataD)
# res <- pairwise.t.test.with.t.and.df(x = dataD$measure, g = dataD$group, paired = TRUE, p.adj="none")
# show(res) # [[3]] = p value table, un adjusted
# res[[5]] # t-value
# res[[6]] # dfs
# res[[3]]
# p.adjust(res[[3]], method = "bonferroni", n = 4) 
# cohensd <- cohensD(subset(dataD, group==1)$measure, subset(dataD, group==2)$measure, method = 'unequal')
# show(cohensd)

dataPlot <- dataAnova
# # pursuit bias magnitude
# ylimLow <- -2
# ylimHigh <- 15
# average response
ylimLow <- 0
ylimHigh <- 6

p <- ggplot(data=dataPlot, aes(x = group, y = measure)) +
        stat_summary(aes(y = measure), fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1.96/sqrt(subN)),
               geom = 'linerange', size = 1) +
        # geom_line(aes(x = group, y = measure, group = sub), size = 0.5, linetype = "dashed") +
        geom_point(aes(x = group, y = measure), size = dotSize, shape = 1, position = position_jitter(w = 0.1, h = 0)) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(1, 0.5), y = c(ylimLow, ylimLow), xend = c(3, 0.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
        scale_y_continuous(name = "Mean pursuit bias magnitude (deg)", breaks = c(ylimLow, 0, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
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
