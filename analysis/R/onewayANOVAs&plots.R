library(ggplot2)
library(ez)
library(Hmisc)
library(lsr)
library(BayesFactor)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/analysis/R")
source("pairwise.t.test.with.t.and.df.R")
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/analysis/R/Exp1")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/results/manuscript/figures/rawPlots/")
### modify these parameters to plot different conditions
dataFileName <- "data_clpGainX_contextMerged_exp1.csv"
pdfFileName <- "clpGain_contextTrials_exp1.pdf"
# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 4
# slope
ylimLow <- 10
ylimHigh <- 50
# PSE 
ylimLow <- -0.15
ylimHigh <- 0.15
# # PSE by preceding perception
# ylimLow <- -0.2
# ylimHigh <- 0.1
# # OSE
# ylimLow <- -0.15
# ylimHigh <- 0.15
# # ASP
# ylimLow <- -1
# ylimHigh <- 5 #for all probability blocks
# # ASP by preceding perception, in the 50% block
# ylimLow <- -1
# ylimHigh <- 2.5 #for all probability blocks

# source("pairwise.t.test.with.t.and.df.R")
data <- read.csv(dataFileName)
subs <- unique(data$sub)
subN <- length(subs)
# # exclude bad fitting...
# data <- subset(data[which(data$sub!=8),])

#### perceptual illusion repeated measures ANOVA
### 2 way for perception--rotational speed x after-reversal direction
## Exp1
sub <- data["sub"]
# prePercept <- data["precededPerception"]
prob <- data["prob"]
# preDir <- data["precededPerception"]
# timeBin <- data["timeBin"]
# measure <- data["PSE"]
measure <- data["measure"]
# dataAnova <- data.frame(sub, preDir, measure)
dataAnova <- data.frame(sub, prob, measure)
# colnames(dataAnova)[2] <- "prePercept"
dataAnova$prob <- as.factor(dataAnova$prob)
dataAnova$sub <- as.factor(dataAnova$sub)
# dataAnova$prePercept <- as.factor(dataAnova$prePercept)
# dataAnova$timeBin <- as.factor(dataAnova$timeBin)
colnames(dataAnova)[3] <- "measure"
# dataAnova <- aggregate(perceptualErrorMean ~ sub * rotationSpeed * exp,
#     data = dataTemp, FUN = "mean")

anovaData <- ezANOVA(dataAnova, dv = .(measure), wid = .(sub),
    within = .(prob), type = 3)
print(anovaData)

# # compute Bayes factor
# bf10 = anovaBF(measure ~ prob + sub, data = dataAnova, 
#              whichRandom="sub")
# bf10

# dataTemp <- data.frame(sub, prob, measure)
# colnames(dataTemp)[3] <- "measure"
# # one sample t-test to zero
# t.test(dataTemp[which(dataTemp$prob==90),]$measure, mu=0)

# # dataPH <- aggregate(. ~ sub * prob, data = dataTemp, FUN = "mean")
# # # show(dataPH[which(dataPH$prob==90), ])
# dataPH <- dataAnova
# # show(dataPH[which(dataPH$prePercept==1), ])
# res <- pairwise.t.test.with.t.and.df(x = dataPH$measure, g = dataPH$prePercept, paired = TRUE, p.adj="none")
# show(res) # [[3]] = p value table, un adjusted
# res[[5]] # t-value
# res[[6]] # dfs
# res[[3]]
# cohensd <- cohensD(subset(dataPH, prePercept==0)$measure, subset(dataPH, prePercept==1)$measure, method = 'paired')
# show(cohensd)
# # p.adjust(res[[3]], method = "bonferroni", n = 10) # interaction between timeWindow & afterReversalD

# # t-test to compare PSE for different preceding perception
# res <- pairwise.t.test.with.t.and.df(x = dataAnova$measure, g = dataAnova$precededPerception, paired = TRUE, p.adj="none")
# # show(res) # [[3]] = p value table, un adjusted
# res[[5]] # t-value
# res[[6]] # dfs
# res[[3]]
# # p.adjust(res[[3]], method = "bonferroni", n = 9) # interaction between timeWindow & afterReversalD
# cohensd <- cohensD(subset(dataAnova, precededPerception==0)$measure, subset(dataAnova, precededPerception==1)$measure, method = 'paired')
# show(cohensd)

# to do the plot, do not turn prob into factor, check the values and data type!
dataPlot <- data.frame(sub, prob, measure)
colnames(dataPlot)[3] <- "measure"
dataPlot$sub <- as.factor(dataPlot$sub)

p <- ggplot(data=dataPlot, aes(x = prob, y = measure)) +
        stat_summary(aes(y = measure), fun.y = mean, geom = "point", shape = 95, size = 17.5) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1.96/sqrt(subN)),
               geom = 'linerange', size = 1) +
        geom_line(aes(x = prob, y = measure, group = sub), size = 0.5, linetype = "dashed") +
        geom_point(aes(x = prob, y = measure), size = dotSize, shape = 1) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(50, 45), y = c(ylimLow, ylimLow), xend = c(90, 45), yend = c(ylimLow, ylimHigh)), size = axisLineWidth) +
        # scale_y_continuous(name = "Anticipatory pursuit velocity (°/s)", breaks = seq(ylimLow, ylimHigh, 0.5), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        scale_y_continuous(name = "PSE", breaks = c(ylimLow, 0, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        # scale_y_continuous(name = "Slope", breaks = seq(ylimLow, ylimHigh， 10), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        # scale_x_continuous(name = "Preceded perceived direction", breaks=c(0, 1), limits = c(-0.5, 1.5), expand = c(0, 0)) + # preceded perception
        scale_x_continuous(name = "Probability of rightward motion", breaks=c(50, 90), limits = c(45, 95), expand = c(0, 0)) + # Exp2&3
        # scale_x_continuous(name = "Probability of rightward motion", breaks=c(50, 70, 90), limits = c(45, 95), expand = c(0, 0)) + # Exp1
        # scale_x_discrete(name = "Probability of rightward motion", breaks=c("50", "70", "90")) +
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