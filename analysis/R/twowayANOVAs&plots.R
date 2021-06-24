library(ggplot2)
library(ez)
library(Hmisc)
library(reshape2)
library(psychReport)
library(lsr)
library(bayestestR)
library(BayesFactor)
library(TOSTER)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/micropursuit/analysis/R")
# source("pairwise.t.test.with.t.and.df.R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/micropursuit/analysis/R/plots/")
### modify these parameters to plot different conditions
varName <- "dirClpChange"
dataFileName <- "summaryDataDiffGroup1.csv"
pdfFileName <- paste("diff_", varName, "_group1.pdf", sep = "")

# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 3
# # slope
# ylimLow <- 10
# ylimHigh <- 50
# # PSE
# ylimLow <- -0.15
# ylimHigh <- 0.15
# # ASP
# ylimLow <- -1
# ylimHigh <- 5
# # # ASP gain
# # ylimLow <- -0.1
# # ylimHigh <- 0.4
# # clp gain in context trials
# ylimLow <- 0
# ylimHigh <- 1

data <- read.csv(dataFileName)
subAll <- unique(data["sub"])
subTotalN <- dim(subAll)[1]
# dataD <- read.csv(dataDFileName)
# data <- data[data.exp==3]
# # exclude bad fitting...
# data <- subset(data[which(data$sub!=8),])

## compare two experiments
# PSE anova
sub <- data["sub"]
rdkApertureAngle <- data["rdkApertureAngle"]
rdkInternalDir <- data["rdkInternalDir"]
measure <- data[varName]
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
# print(anovaData)
aovEffectSize(anovaData, 'pes')

# # Equivalence test for the differences of PSE between experiments
# dataE <- data.frame(matrix(ncol=3,nrow=dim(subAll)[1], dimnames=list(NULL, c("sub", "exp1", "exp2"))))
# for (subN in 1:subTotalN) {
#   data150 <- dataAnova[which(sub==subAll[subN, 1] & exp==1 & prob==50), ]$measure
#   data190 <- dataAnova[which(sub==subAll[subN, 1] & exp==1 & prob==90), ]$measure
#   data250 <- dataAnova[which(sub==subAll[subN, 1] & exp!=1 & prob==50), ]$measure
#   data290 <- dataAnova[which(sub==subAll[subN, 1] & exp!=1 & prob==90), ]$measure 
#   dataE["sub"][subN, 1] <- subN
#   dataE["exp1"][subN, 1] <- data190-data150
#   dataE["exp2"][subN, 1] <- data290-data250
# }
# show(dataE)

# dataTOSTpaired(data = dataE, pairs = list((c(i1="exp1",i2="exp2"))), low_eqbound = -0.36, high_eqbound = 0.36, eqbound_type = "d", alpha = 0.05, desc = TRUE, plots = TRUE)

# # compute Bayes Factor inclusion...
# bf <- anovaBF(measure ~ prob + timeBin + prob*timeBin + sub, data = dataAnova, 
#              whichRandom="sub")
# bayesfactor_inclusion(bf, match_models = TRUE)

# p <- ggplot(dataAnova, aes(x = prob, y = measure, color = exp)) +
#         stat_summary(aes(y = measure), fun.y = mean, geom = "point", shape = 95, size = 15) +
#         stat_summary(fun.data = 'mean_sdl',
#                fun.args = list(mult = 1.96/sqrt(subTotalN)),
#                geom = 'errorbar', width = .1) +
# # geom = 'smooth', se = 'TRUE') +
#         # stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
#         geom_point(aes(x = prob, y = measure), size = dotSize, shape = 1) +
#         # geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(50, 40), xend = c(90, 40), y = c(-0.1, -0.1), yend = c(-0.1, 0.15)), size = axisLineWidth) +
#         scale_y_continuous(name = "Anticipatory pursuit velocity (°/s)") + #, limits = c(-0.1, 0.55), expand = c(0, 0)) +
#         # scale_y_continuous(name = "PSE") + 
#         scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "90")) +
#         # scale_x_discrete(name = "Probability of rightward motion", breaks=c(50, 70, 90)) +
#         # scale_colour_discrete(name = "After reversal\ndirection", labels = c("CCW", "CW")) +
#         theme(axis.text=element_text(colour="black"),
#                       axis.ticks=element_line(colour="black", size = axisLineWidth),
#                       panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(),
#                       panel.border = element_blank(),
#                       panel.background = element_blank(),
#                       text = element_text(size = textSize, colour = "black"),
#                       legend.background = element_rect(fill="transparent"),
#                       legend.key = element_rect(colour = "transparent", fill = "white"))
#         # facet_wrap(~exp)
# print(p)
# ggsave(paste(plotFolder, pdfFileName, sep = ""))

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
# dataPlot$rdkApertureAngle <- as.factor(dataPlot$rdkApertureAngle)
# is.numeric(dataPlot$timeBin)
dataPlot$rdkInternalDir <- as.factor(dataPlot$rdkInternalDir)
# dataPlot <- aggregate(measure ~ rdkInternalDir+prob, data = dataPlot, FUN = "mean")
# show(dataPlot)

# response
ylimLow <- -15
ylimHigh <- 15

# # for time bin plots
# p <- ggplot(dataPlot, aes(x = timeBin, y = measure, color = prob)) +
#         stat_summary(fun.y = mean, geom = "point", shape = 95, size = 17.5) +
#         stat_summary(fun.y = mean, geom = "line", width = 1) +
#         stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subTotalN)), geom = 'errorbar', width = 1.5, size = 1) +
#         scale_x_continuous(name = "time bin of trials", breaks=c(1, 2), limits = c(0.5, 2.5), expand = c(0, 0)) +
p <- ggplot(dataPlot, aes(x = rdkApertureAngle, y = measure, color = rdkInternalDir)) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun = mean, geom = "line", width = 1) +
        stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subTotalN)), geom = 'errorbar', width = 1.5, size = 1) +
        stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_point(size = dotSize, shape = 1) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-9, -12), y = c(ylimLow, ylimLow), xend = c(9, -12), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
        scale_y_continuous(name = "Bias in pursuit direction (deg)", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = varName, limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) + 
        scale_x_continuous(name = "Aperture motion direction (deg)", breaks=c(-9,-6,-3,0,3,6,9), limits = c(-12, 12), expand = c(0, 0)) +
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

# ## t-test of the difference
# sub <- dataD["sub"]
# exp <- dataD["exp"]
# measure <- dataD["slopeDiff"]
# dataDtemp <- data.frame(sub, exp, measure)
# dataDtemp$sub <- as.factor(dataDtemp$sub)
# dataDtemp$exp <- as.factor(dataDtemp$exp)
# colnames(dataDtemp)[3] <- "measure"
# # dataDttest <- aggregate(measure ~ exp, data = dataDtemp, FUN = "mean")

# # res <- pairwise.t.test.with.t.and.df(x = dataDtemp$measure, g = dataDtemp$exp, paired = TRUE, p.adj="none")
# # show(res) # [[3]] = p value table, un adjusted
# # res[[5]] # t-value
# # res[[6]] # dfs
# # res[[3]]
# # p.adjust(res[[3]], method = "bonferroni", n = 3) 

# # # bias in PSE
# # ylimLow <- -0.05
# # ylimHigh <- 0.2
# # # bias in ASP
# # ylimLow <- 0
# # ylimHigh <- 3
# # bias in slope
# ylimLow <- -40
# ylimHigh <- 30

# p <- ggplot(dataDtemp, aes(x = exp, y = measure)) +
#         stat_summary(aes(y = measure), fun.y = mean, geom = "point", shape = 95, size = 15) +
#         stat_summary(fun.data = 'mean_sdl',
#                fun.args = list(mult = 1.96/sqrt(subTotalN)),
#                geom = 'linerange', size = 1) +
#         geom_line(aes(x = exp, y = measure, group = sub), size = 0.5, linetype = "dashed") +
#         geom_point(aes(x = exp, y = measure), size = dotSize, shape = 1) +
#         geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(0), y = c(ylimLow), xend = c(0), yend = c(ylimHigh)), size = axisLineWidth) +
#         # scale_y_continuous(name = "Bias of PSE") + #, limits = c(0, 0.15), expand = c(0, 0.01)) +
#         scale_y_continuous(name = "Bias of slope") + #, limits = c(0, 0.15), expand = c(0, 0.01)) +
#         # scale_y_continuous(name = "Bias of anticipatory pursuit velocity(deg/s)") + #, limits = c(0, 0.15), expand = c(0, 0.01)) +
#         scale_x_discrete(name = "Experiment", limits = c("1", "3"), labels = c("1" = "Exp1", "3" = "Exp3")) +
#         # scale_x_discrete(name = "Experiment", limits = c("1", "2"), labels = c("1" = "Exp1", "2" = "Exp2")) +
#         theme(axis.text=element_text(colour="black"),
#               axis.ticks=element_line(colour="black", size = axisLineWidth),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.border = element_blank(),
#               panel.background = element_blank(),
#               text = element_text(size = textSize, colour = "black"),
#               legend.background = element_rect(fill="transparent"),
#               legend.key = element_rect(colour = "transparent", fill = "white"))
#         # facet_wrap(~prob)
# print(p)
# ggsave(paste(plotFolder, pdfFileNameD, sep = ""))