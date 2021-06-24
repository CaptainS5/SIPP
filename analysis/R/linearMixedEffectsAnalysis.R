library(ggplot2)
library(rmcorr)
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
dataFileName <- "eyeVSpercept_exp1.csv"
pdfFileName <- "corrSlopeASPvsPSE_exp1.pdf"
# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 4
# PSE 
ylimLow <- -0.15
ylimHigh <- 0.15
# # ASP
# ylimLow <- -1
# ylimHigh <- 5 #for all probability blocks
# # ASP by preceding perception, in the 50% block
# ylimLow <- -1
# ylimHigh <- 2.5 #for all probability blocks

# source("pairwise.t.test.with.t.and.df.R")
data <- read.csv(dataFileName)

#### perceptual illusion repeated measures ANOVA
### 2 way for perception--rotational speed x after-reversal direction
## Exp1
sub <- data["sub"]
prob <- data["prob"]
measure1 <- data["asp"]
measure2 <- data["PSE"]
dataCorr <- data.frame(sub, prob, measure1, measure2)
colnames(dataCorr)[3] <- "measure1"
colnames(dataCorr)[4] <- "measure2"
# dataCorr$prob <- as.factor(dataCorr$prob)
dataCorr$sub <- as.factor(dataCorr$sub)
# show(dataCorr)
# # using rmcorr
# corrResults <- rmcorr(participant = sub, measure1 = measure1, measure2 = measure2, dataset = dataCorr, CIs =
# c("analytic", "bootstrap"), nreps = 1000, bstrap.out = T)
# # show(corrResults)
# pdf(paste(plotFolder, pdfFileName, sep = ""))
# p <- plot(corrResults, overall = F, lty = 2, xlab = "Effects on PSE", ylab = "Effects on anticipatory pursuit velocity")
# dev.off()

# using lme4
fm1 <- lmer(measure2 ~ measure1 + prob + (1 | sub), dataCorr, REML = F)
summary(fm1)

# fitting the probability-induced slopes of asp/PSE for each observer
subAll <- unique(data["sub"])
slopeData <- data.frame(matrix(ncol=3,nrow=dim(subAll[1]), dimnames=list(NULL, c("sub", "asp", "PSE"))))
for (subN in 1:dim(subAll)[1]) {
	dataTemp <- subset(data, sub==subAll[subN, 1])
	aspLM <- lm(asp ~ prob, dataTemp)
	pseLM <- lm(PSE ~ prob, dataTemp)
	slopeData["sub"][subN, 1] <- subN
	slopeData["asp"][subN, 1] <- coef(aspLM)["prob"]
	slopeData["PSE"][subN, 1] <- coef(pseLM)["prob"]
}
show(slopeData)

ylimLow <- .0005
ylimHigh <-.003

 p <- ggplot(data=slopeData, aes(x = asp, y = PSE)) +
        geom_point(size = dotSize) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(.01, .005), y = c(ylimLow, ylimLow), xend = c(.06, .005), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = F) +
        scale_y_continuous(name = "Slope of PSE", breaks = c(ylimLow, 0, ylimHigh), limits = c(ylimLow, ylimHigh+.0005), expand = c(0, 0)) +
        scale_x_continuous(name = "Slope of anticipatory pursuit velocity", breaks = seq(0.01, .06, .01), limits = c(.005, .065), expand = c(0, 0)) +
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
        # facet_wrap(~prob)
print(p)
ggsave(paste(plotFolder, pdfFileName, sep = ""))

# p <- ggplot(data=dataCorr, aes(x = measure1, y = measure2, group = prob)) +
#         geom_point(aes(shape = prob, color = sub), size = dotSize) +
#         geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-1, -1.5), y = c(ylimLow, ylimLow), xend = c(5, -1.5), yend = c(ylimLow, ylimHigh)), size = axisLineWidth, inherit.aes = F) +
#         scale_y_continuous(name = "PSE", breaks = c(ylimLow, 0, ylimHigh), limits = c(ylimLow, ylimHigh), expand = c(0, 0)) +
#         scale_x_continuous(name = "Anticipatory pursuit velocity (°/s)", breaks = seq(-1, 5, 1), limits = c(-1.5, 5.5), expand = c(0, 0)) +
#         # scale_colour_discrete(name = "After reversal\ndirection", labels = c("CCW", "CW")) +
#         theme(axis.text=element_text(colour="black"),
#               axis.ticks=element_line(colour="black", size = axisLineWidth),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.border = element_blank(),
#               panel.background = element_blank(),
#               text = element_text(size = textSize, colour = "black"),
#               legend.background = element_rect(fill="transparent"),
#               legend.key = element_rect(colour = "transparent", fill = "white"),
#               aspect.ratio=1)
#         # facet_wrap(~prob)
# print(p)
# ggsave(paste(plotFolder, pdfFileName, sep = ""))