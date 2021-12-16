library(ggplot2)
library(ez)
library(Hmisc)
library(lsr)
library(BayesFactor)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/micropursuit/analysis/R")
source("pairwise.t.test.with.t.and.df.R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/micropursuit/analysis/R/plots/")

### modify these parameters to plot different conditions
# varName <- "dirClp"
dataFileName <- "tempCorrDataSubgroup.csv"

# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 3

# source("pairwise.t.test.with.t.and.df.R")
data <- read.csv(dataFileName)
subAll <- unique(data$sub)
subTotalN <- length(subAll)

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


#### Comparison between perceptual bias groups
### about the temporal relationship
varName <- c("wMiddle")
pdfFileName <- paste("groupCompTemporal_", varName, ".pdf", sep = "")

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
ylimHigh <- 10

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
