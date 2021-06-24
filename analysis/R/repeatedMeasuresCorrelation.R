library(ggplot2)
library(rmcorr)
library(lme4)
library(lmerTest)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/micropursuit/analysis/R")
# source("pairwise.t.test.with.t.and.df.R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/micropursuit/analysis/R/plots/")
### modify these parameters to plot different conditions
varName <- "dirClpEarly"
dataFileName <- "summaryDataDiffGroup1.csv"
# pdfFileName <- paste("lmeDiff_allVSperception_", varName, "_group1.pdf", sep = "")
pdfFileName <- paste("correct_lmeDiff_", varName, "VSperception_group1.pdf", sep = "")
# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 3

# load data
data <- read.csv(dataFileName)
subAll <- unique(data["sub"])
subTotalN <- dim(subAll)[1]
# dataD <- read.csv(dataDFileName)
# data <- data[data.exp==3]
# # exclude bad fitting...
# data <- subset(data[which(data$sub!=8),])

# prepare the dataset
sub <- data["sub"]
rdkApertureAngle <- data["rdkApertureAngle"]
rdkInternalDir <- data["rdkInternalDir"]
measure1 <- data["response"]
measure2 <- data["dirOlp"]
measure3 <- data["dirClpEarly"]
measure4 <- data["dirClpLate"]

## linear mixed effects model for all eye measures together
dataCorr <- data.frame(sub, rdkApertureAngle, rdkInternalDir, rdkInternalDir, measure1, measure2, measure3, measure4)
colnames(dataCorr)[5] <- "measure1"
colnames(dataCorr)[6] <- "measure2"
colnames(dataCorr)[7] <- "measure3"
colnames(dataCorr)[8] <- "measure4"
dataCorr$sub <- as.factor(dataCorr$sub)
# dataCorr$rdkApertureAngle <- as.factor(dataCorr$rdkApertureAngle)
dataCorr$rdkInternalDir <- as.factor(dataCorr$rdkInternalDir)
# show(dataCorr)

# using lme4
# fm1 <- lmer(measure1 ~ rdkInternalDir*rdkApertureAngle + (1 | rdkInternalDir:sub) + (1 | rdkApertureAngle:sub), dataCorr, REML = F)
fm1 <- lmer(measure1 ~ rdkInternalDir + measure2 + measure3 + measure4 + (1 | rdkApertureAngle) + (1 | sub), dataCorr, REML = F)
# fm1 <- lmer(measure1 ~ rdkInternalDir + measure2 + measure3 + measure4 + (1 | sub) , dataCorr, REML = F)
summary(fm1)
plot(fm1, which = 2)
plot(fm1, rdkInternalDir ~ resid(.), abline = 0 )
qqnorm(resid(fm1))

# # anovaModelRM = aov(measure1 ~ rdkInternalDir*rdkApertureAngle + Error(sub/rdkInternalDir*rdkApertureAngle), dataCorr)
# # summary(anovaModelRM)

# # # separate data into subgroups
# # dataOne <- subset(dataCorr, rdkInternalDir==-90)
# # dataTwo <- subset(dataCorr, rdkInternalDir==90)
# # newdat <- expand.grid(rdkInternalDir = unique(dataCorr$rdkInternalDir), measure1 = c(min(dataCorr$measure1), max(dataCorr$measure1)))
# # show(dataCorr[, c(5,7)])

# # # first generate the newdat frame for plotting the prediction lines
# # newdat <- expand.grid(rdkInternalDir = unique(dataCorr$rdkInternalDir), sub = unique(dataCorr$sub), rdkApertureAngle = unique(dataCorr$rdkApertureAngle), measure2 = 0, measure4 = 0) 
# # # then fill in the min and max values for each observer...
# # # not done yet... but we can use a different way to do this, see below
# # show(newdat)
# # # pred <- predict(fm1, newdat)

## plotting prediction of a single eye measure
# generate the plot data
measure2 <- data[varName]
dataPlot <- data.frame(sub, rdkApertureAngle, rdkInternalDir, rdkInternalDir, measure1, measure2)
colnames(dataPlot)[5] <- "measure1"
colnames(dataPlot)[6] <- "measure2"
dataPlot$sub <- as.factor(dataPlot$sub)
dataPlot$rdkInternalDir <- as.factor(dataPlot$rdkInternalDir)

# generate the data for plotting the prediction of one factor...
# fm2 <- lmer(measure1 ~ rdkInternalDir + measure2 + (1 + measure2 | sub), dataPlot, REML = F)
# fm2 <- lmer(measure1 ~ rdkInternalDir + measure2 + (1 | sub) + (0 + measure2 | sub), dataPlot, REML = F)
# fm2 <- lmer(measure1 ~ rdkInternalDir + measure2 + (1 | sub) + + (1 | rdkInternalDir:sub), dataPlot, REML = F)
fm2 <- lmer(measure1 ~ rdkInternalDir + measure2 + (1 | rdkApertureAngle)+ (1 | sub), dataPlot, REML = F)
# fm2 <- lmer(measure1 ~ rdkInternalDir + measure2 + (1 | sub), dataPlot, REML = F)
# fm2 <- lmer(measure1 ~ rdkInternalDir + measure2 + (1 + measure2 | sub), dataPlot, REML = F)
summary(fm2)
# plot(fm2, which=1)
# plot(fm2, rdkInternalDir ~ resid(.), abline = 0 )
# qqnorm(resid(fm2))

# p <- ggplot(dataCorr, aes(x = rdkApertureAngle, y = measure1, colour = sub, shape = rdkInternalDir)) +
p <- ggplot(dataPlot, aes(x = measure2, y = measure1, colour = sub, shape = rdkInternalDir)) +
        geom_point(size = dotSize) +
        # geom_line(data = dataCorr[, c(5,7)], aes(y = predict(fm1))) +
        geom_line(aes(y=predict(fm2))) +
        # scale_size_manual(name = "Predictions", values = c("Observer" = 0.5, "Group" = 2)) +
        # geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-1, -1.5), y = c(ylimLow, ylimLow), xend = c(5, -1.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = F) +
        scale_x_continuous(name = paste("Bias in", varName, "(deg)", sep = " ")) +
        # , breaks = c(ylimLow, 0, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
        scale_y_continuous(name = "Bias in perceived direction (deg)") +
        # , breaks = seq(-1, 5, 1), limits = c(-1.5, 5.5), expand = c(0, 0)) +
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
              aspect.ratio=1)  
        # facet_wrap(~rdkInternalDir)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))

# ## repeated-measures correlation and plot
# my.rmc <- rmcorr(participant = sub, measure1 = measure3, measure2 = measure4, dataset = dataCorr, CI.level = 0.95, CIs = c("analytic", "bootstrap"), nreps = 1000)
# my.rmc

# # plot the results
# pdf(paste("corrDiff_", varName, "VSperception_group1.pdf", sep = ""))
# p <- plot(my.rmc, overall = FALSE, lty = 2, ylab = "Bias in perceived direction (deg)", xlab = paste("Bias in", varName, "(deg)", sep = " "))
# print(p)
# dev.off()