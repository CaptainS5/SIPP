library(ggplot2)
library(ez)
library(Hmisc)
library(reshape2)
library(psychReport)
library(lsr)
library(bayestestR)
library(BayesFactor)
library(lme4)
library(lmerTest)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/analysis/R")
source("pairwise.t.test.with.t.and.df.R")
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/analysis/R/Exp1")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/results/manuscript/figures/rawPlots/")
### modify these parameters to plot different conditions
dataFileName <- "data_clpGainX_probeConsistency_exp1.csv"
pdfInteractionFileName <- "clpGainX_probeConsistency_exp1.pdf"
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
# clp velocity in probe trials
ylimLow <- -4
ylimHigh <- 4
# # clp gain in probe trials
ylimLow <- -0.4
ylimHigh <- 0.7

data <- read.csv(dataFileName)
subs <- unique(data$sub)
subN <- length(subs)
# dataD <- read.csv(dataDFileName)
# data <- data[data.exp==3]
# # exclude bad fitting...
# data <- subset(data[which(data$sub!=8),])

# ## linear mixed effects model
# sub <- data["sub"]
# prob <- data["prob"]
# visualDir <- data["visualDir"]
# perceivedDir <- data["dir"]
# measure <- data["measure"]
# dataCorr <- data.frame(sub, prob, visualDir, perceivedDir, measure)
# colnames(dataCorr)[4] <- "perceivedDir"
# # colnames(dataCorr)[5] <- "measure"
# # dataCorr$prob <- as.factor(dataCorr$prob)
# dataCorr$sub <- as.factor(dataCorr$sub)
# dataCorr$visualDir <- as.factor(dataCorr$visualDir)
# dataCorr$perceivedDir <- as.factor(dataCorr$perceivedDir)
# fm1 <- lmer(measure ~ visualDir + perceivedDir + visualDir*perceivedDir + prob + prob*visualDir*perceivedDir + (1 | sub), dataCorr, REML = F)
# summary(fm1)

# p <- ggplot(dataCorr, aes(x = prob, y = measure, color = interaction(visualDir, perceivedDir))) +
#         stat_summary(fun= mean, geom = "point", shape = 95, size = 17.5, position = position_dodge(width=15)) +
#         stat_summary(fun.data = 'mean_sdl',
#                fun.args = list(mult = 1.96/sqrt(subN)),
#                geom = 'linerange', size = 1, position = position_dodge(width=15)) +
#         geom_point(size = dotSize, shape = 1, position = position_jitterdodge(  jitter.width = NULL, jitter.height = 0, dodge.width = 15)) +
#         geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(50, 40), xend = c(90, 40), y = c(ylimLow, ylimLow), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
#         scale_y_continuous(name = "Steady-state pursuit gain", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) + 
#         # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "90")) +
#         scale_x_continuous(name = "Probability of rightward motion", limits = c(40, 100), breaks=c(50, 70, 90), expand = c(0, 0)) +
#         # scale_colour_discrete(name = "Motion direction", labels = c("CCW", "CW")) +
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

## ANOVA
sub <- data["sub"]
prob <- data["prob"]
congruency <- data["congruency"]
# dir <- data["dir"]
measure <- data["measure"]
dataAnova <- data.frame(sub, prob, congruency, measure)
colnames(dataAnova)[4] <- "measure"
# colnames(dataAnova)[3] <- "visualDirection"
dataAnova$prob <- as.factor(dataAnova$prob)
dataAnova$sub <- as.factor(dataAnova$sub)
dataAnova$congruency <- as.factor(dataAnova$congruency)
# dataAnova$visualDirection <- as.factor(dataAnova$visualDirection)

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(prob, congruency), type = 3, return_aov = TRUE, detailed = TRUE)
# print(anovaData)
aovEffectSize(anovaData, 'pes')

# bf <- anovaBF(measure ~ prob + dir + prob*dir + sub, data = dataAnova, 
#              whichRandom="sub")
# bayesfactor_inclusion(bf, match_models = TRUE)

# # one sample t-test to zero
# res <- t.test(dataAnova[which(dataAnova$prob==50 & dataAnova$congruency==-1),]$measure, mu=0, alternative = 'less')
# show(res)
# # p.adjust(res[[3]], method = "bonferroni", n = 3) # 
# cohensd <- cohensD(subset(dataAnova, prob==50 & congruency==-1)$measure, mu = 0)
# show(cohensd)

## interaction plot
dataPlot <- data.frame(sub, prob, congruency, measure)
colnames(dataPlot)[4] <- "measure"
# colnames(dataPlot)[3] <- "visualDirection"
dataPlot$sub <- as.factor(dataPlot$sub)
# dataPlot$prob <- as.factor(dataPlot$prob)
# is.numeric(dataPlot$timeBin)
dataPlot$congruency <- as.factor(dataPlot$congruency)
# dataPlot$visualDirection <- as.factor(dataPlot$visualDirection)
show(dataPlot[which(dataPlot$congruency==-1),]$measure)

p <- ggplot(dataPlot, aes(x = prob, y = measure, color = congruency)) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5, position = position_dodge(width=10)) +
        # stat_summary(fun.y = mean, geom = "line", width = 1) +
        stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subN)), geom = 'errorbar', width = 1.5, size = 1, position = position_dodge(width=10)) +
        # stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_point(size = dotSize, shape = 1, position = position_jitterdodge(  jitter.width = NULL, jitter.height = 0, dodge.width = 10)) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(50, 40), y = c(ylimLow, ylimLow), xend = c(90, 40), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
        # scale_y_continuous(name = "Anticipatory pursuit velocity (°/s)", breaks = seq(ylimLow, ylimHigh, 1), expand = c(0, 0)) +
        # scale_y_continuous(name = "Visually-guided pursuit velocity (°/s)", breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        scale_y_continuous(name = "Visually-guided pursuit gain", breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
        coord_cartesian(ylim=c(ylimLow, ylimHigh)) +
        # scale_y_continuous(name = "PSE", limits = c(ylimLow, ylimHigh), breaks = c(ylimLow, 0, ylimHigh), expand = c(0, 0)) + 
        scale_x_continuous(name = "Probability of rightward motion", breaks=c(50, 70, 90), limits = c(40, 100), expand = c(0, 0)) +
        # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "90")) +
        # scale_colour_discrete(name = "After reversal\ndirection", labels = c("CCW", "CW")) +
        scale_color_manual(name = "Motion\ncongruency"， values=c("#293352", "#4e84c4")) +
        theme(axis.text=element_text(colour="black"),
                      axis.ticks=element_line(colour="black", size = axisLineWidth),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size = textSize, colour = "black"),
                      legend.background = element_rect(fill="transparent"),
                      legend.key = element_rect(colour = "transparent", fill = "white"))
        # facet_wrap(~exp)
print(p)
ggsave(paste(plotFolder, pdfInteractionFileName, sep = ""))

# ## plot trial number...
# sub <- data["sub"]
# prob <- data["prob"]
# measure <- data["trialNumber"]
# congruency <- data["consistency"]
# dataTN <- data.frame(sub, prob, congruency, measure)
# colnames(dataTN)[4] <- "measure"
# colnames(dataTN)[3] <- "congruency"
# # dataTN$prob <- as.factor(dataTN$prob)
# dataTN$sub <- as.factor(dataTN$sub)
# dataTN$congruency <- as.factor(dataTN$congruency)

# ylimLow <- 0
# ylimHigh <- 160

# p <- ggplot(dataTN, aes(x = prob, y = measure, color = congruency)) +
#         stat_summary(fun = mean, geom = "point", shape = 95, size = 17.5, position = position_dodge(width=10)) +
#         # stat_summary(fun.y = mean, geom = "line", width = 1) +
#         stat_summary(fun.data = 'mean_sdl', fun.args = list(mult = 1.96/sqrt(subN)), geom = 'errorbar', width = 1.5, size = 1, position = position_dodge(width=10)) +
#         # stat_summary(aes(y = measure), fun.data = mean_se, geom = "errorbar", width = 0.1) +
#         geom_point(size = dotSize, shape = 1, position = position_jitterdodge(  jitter.width = NULL, jitter.height = 0, dodge.width = 10)) +
#         geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(50, 40), y = c(ylimLow, ylimLow), xend = c(90, 40), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = FALSE) +
#         scale_y_continuous(name = "Trial number", breaks = seq(ylimLow, ylimHigh, 20), expand = c(0, 0)) +
#         coord_cartesian(ylim=c(ylimLow, ylimHigh)) +
#         scale_x_continuous(name = "Probability of rightward motion", breaks=c(50, 70, 90), limits = c(40, 100), expand = c(0, 0)) +
#         # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "90")) +
#         scale_colour_discrete(name = "Congruency", labels = c("incongruent", "congruent")) +
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
# ggsave(paste(plotFolder, pdfInteractionFileName, sep = ""))

# ## Exp1, plot, perceived direction x visual consistency
# sub <- data["sub"]
# visualDir <- data["visualDir"]
# # visualConsistency <- data["visualConsistency"]
# perceivedDir <- data["dir"]
# prob <- data["prob"]
# measure <- data["measure"]
# dataAnova <- data.frame(sub, prob, visualDir, perceivedDir, measure)
# colnames(dataAnova)[4] <- "perceivedDir"
# colnames(dataAnova)[5] <- "measure"
# # # dataAnova$prob <- as.factor(dataAnova$prob)
# # dataAnova$sub <- as.factor(dataAnova$sub)
# # dataAnova$visualConssitency <- as.factor(dataAnova$visualConsistency)
# # dataAnova$perceivedDir <- as.factor(dataAnova$perceivedDir)

# # anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
# #     within = .(prob, exp), type = 3, return_aov = TRUE, detailed = TRUE)
# # # print(anovaData)
# # aovEffectSize(anovaData, 'pes')

# # # compute Bayes Factor inclusion...
# # bf <- anovaBF(measure ~ prob + timeBin + prob*timeBin + sub, data = dataAnova, 
# #              whichRandom="sub")
# # bayesfactor_inclusion(bf, match_models = TRUE)
# # show(dataAnova)


# ## difference between consistent and inconsistent visual by perceived direction
# sub <- data["sub"]
# prob <- data["prob"]
# # visualDir <- data["visualDir"]
# measure <- data["perceptEffect"]
# dataAnova <- data.frame(sub, prob, measure)
# colnames(dataAnova)[3] <- "measure"
# # dataAnova$prob <- as.factor(dataAnova$prob)
# dataAnova$sub <- as.factor(dataAnova$sub)
# # dataAnova$visualDir <- as.factor(dataAnova$visualDir)

# anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
#     within = .(prob), type = 3)
# print(anovaData)

# # t-test for each probability condition...
# res <- t.test(dataAnova[which(dataAnova$prob==90),]$measure, mu=0)
# res[[3]]
# p.adjust(res[[3]], method = "bonferroni", n = 3) 

# res <- pairwise.t.test.with.t.and.df(x = dataAnova$measure, g = dataAnova$sub, paired = TRUE, p.adj="none")
# show(res) # [[3]] = p value table, un adjusted
# res[[5]] # t-value
# res[[6]] # dfs
# res[[3]]
# p.adjust(res[[3]], method = "bonferroni", n = 4) 
# cohensd <- cohensD(subset(dataD, prob==50)$measure, subset(dataD, prob==90)$measure, method = 'paired')
# show(cohensd)

# # dataPlot <- aggregate(measure ~ sub * prob,
#     # data = dataAnova, FUN = "mean")

# p <- ggplot(data=dataAnova, aes(x = prob, y = measure)) +
#         stat_summary(fun.y = mean, geom = "point", shape = 95, size = 17.5) +
#         stat_summary(fun.data = 'mean_sdl',
#                fun.args = list(mult = 1.96/sqrt(subN)),
#                geom = 'linerange', size = 1) +
#         geom_line(aes(x = prob, y = measure, group = sub), size = 0.5, linetype = "dashed") +
#         geom_point(aes(x = prob, y = measure), size = dotSize, shape = 1) +
#         # geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(0, -0.5), y = c(ylimLow, ylimLow), xend = c(1, -0.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
#         scale_y_continuous(name = "Effect of perception") + #, breaks = seq(ylimLow, ylimHigh, 0.5), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
#         # scale_x_continuous(name = "Probability of rightward motion", breaks=c(50, 90), limits = c(45, 95), expand = c(0, 0)) + # Exp2&3
#         scale_x_continuous(name = "Probability of rightward motion", breaks=c(50, 70, 90), limits = c(45, 95), expand = c(0, 0)) + # Exp1
#         # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "70", "90")) +
#         # scale_colour_discrete(name = "After reversal\ndirection", labels = c("CCW", "CW")) +
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
# ggsave(paste(plotFolder, pdfFileName, sep = ""))