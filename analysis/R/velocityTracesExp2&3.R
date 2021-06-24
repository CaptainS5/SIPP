library(ggplot2)
library(RColorBrewer)
library(matrixStats)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/analysis/R/Exp3")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/2ndYear/AnticipatoryPursuit/AnticipatoryPursuitMotionPerception/results/manuscript/figures/rawPlots/")
### modify these parameters to plot different conditions
dataFileName <- "velocityTrace_congruent_exp3_prob"
pdfFileName <- "velocityTrace_congruent_all_Exp3.pdf"
# dataFileName <- "slidingW_perception_"
# pdfFileName <- "slidingW_perception_STE_all_Exp3.pdf"
n <- 8 # number of participants
# for plotting
vtPlotWidth <- 16 # width for the pdf file
textSize <- 25
axisLineWidth <- 0.5
velAlpha <- 0.2
# perceptual trials, velocity traces
ylimLow <- -4.2 # range for x axis line
ylimHigh <- 4.2
ytickLow <- -4 # range of y axis line
ytickHigh <- 4
# # standard trials
# ylimLow <- -11 # range for x axis line
# ylimHigh <- 11
# ytickLow <- -10 # range of y axis line
# ytickHigh <- 10

# load and process data
# mean velocity traces in different conditions
# first get mean for each condition and the confidence intervals
probName <-c("50", "90")
timePoints <- seq(from = -200, to = 700, by = 1)
velTrace <- list()

for (probI in 1:2) {
    fileName = paste(dataFileName,probName[probI],".csv", sep = "")
    velData <- read.csv(fileName, header = FALSE, sep=",")
    velData <- as.matrix(velData)

    for (conI in 1:2) {
      idxStart <- (conI-1)*n+1
      idxEnd <- conI*n

      velDataT <- velData[idxStart:idxEnd, 1000:1900] # -200 to 700 ms
      tempM <- colMeans(velDataT, na.rm=TRUE)
      velMean <- as.matrix(tempM)
      velStd <- colSds(velDataT, na.rm=TRUE)
      error <- qt(0.975, df=n-1)*velStd/sqrt(n)
      velLower <- tempM-error
      velUpper <- tempM+error
      velTrace[[(probI-1)*2+conI]] <- data.frame(timePoints, velMean, velLower, velUpper)
    }
}
colourCode <- c("#dad9d9", "#737474", "#000000") # different shades of grey
# colourCode <- c("#c6dbef", "#4292c6", "#08306b") # all blue
###

## plot velocity traces with 95% CI
# pdf(paste(folder1, "velocityTraces.pdf", sep = ""))
p <- ggplot(velTrace[[1]], aes(x = timePoints, y = velMean)) +
        geom_line(aes(colour = "50%"), size = 1, linetype = "twodash") + 
        # geom_ribbon(aes(x = timePoints, ymin=velLower, ymax=velUpper, fill = "50%"), alpha = velAlpha) +
        geom_line(data = velTrace[[2]], aes(x = timePoints, y = velMean, colour = "50%"), size = 1, linetype = "solid") + 
        # geom_ribbon(data = velTrace[[2]], aes(x = timePoints, ymin=velLower, ymax=velUpper, fill = "50%"), alpha = velAlpha) +
        geom_line(data = velTrace[[3]], aes(x = timePoints, y = velMean, colour = "90%"), size = 1, linetype = "twodash") + 
        # geom_ribbon(data = velTrace[[3]], aes(x = timePoints, ymin=velLower, ymax=velUpper, fill = "90%"), alpha = velAlpha) +
        geom_line(data = velTrace[[4]], aes(x = timePoints, y = velMean, colour = "90%"), size = 1, linetype = "solid") + 
        # geom_ribbon(data = velTrace[[4]], aes(x = timePoints, ymin=velLower, ymax=velUpper, fill = "90%"), alpha = velAlpha) +
        geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(-200, -250), xend = c(700, -250), y = c(ylimLow, ytickLow), yend = c(ylimLow, ytickHigh)), size = axisLineWidth) +
        scale_y_continuous(name = "Horinzontal eye velocity (°/s)",limits = c(ylimLow, ylimHigh), breaks=seq(ytickLow, ytickHigh, 2), expand = c(0, 0)) +
        scale_x_continuous(name = "Time (ms)", limits = c(-250, 700), breaks=seq(-200, 700, 100), expand = c(0, 0)) +
        # coord_cartesian(ylim=c(-4, 4)) +
        scale_colour_manual("Probability of rightward motion",
                            breaks = c("50%", "90%"),
                            values = c("50%" = colourCode[1], "90%" = colourCode[3])) +
        # scale_fill_manual("Probability of rightward motion",
        #                     breaks = c("50%", "70%", "90%"),
        #                     values = c("50%" = colourCode[1], "70%" = colourCode[2], "90%" = colourCode[3])) +
        theme(axis.text=element_text(colour="black"),
              axis.ticks=element_line(colour="black", size = axisLineWidth),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              text = element_text(size = textSize),
              legend.background = element_rect(fill="transparent"),
              legend.key = element_rect(colour = "transparent", fill = "white"))
print(p)
# dev.off()
ggsave(paste(plotFolder, pdfFileName, sep = ""), width = vtPlotWidth)

# sliding window analysis
# mean velocity traces in different conditions
# first get mean for each condition and the confidence intervals

# # perceptual trials, sliding windows, adjust according to needs
# # for perception
# ylimLow <- -0.3 # range for x axis line
# ylimHigh <- 0.3
# ytickLow <- -0.3 # range of y axis line
# ytickHigh <- 0.3

# # # for perception ratio
# # ylimLow <- -0.25 # range for x axis line
# # ylimHigh <- 1.25
# # ytickLow <-  -0.25 # range of y axis line
# # ytickHigh <- 1.25

# # # for APvelX
# # ylimLow <- -0.5 # range for x axis line
# # ylimHigh <- 2
# # ytickLow <- -0.5 # range of y axis line
# # ytickHigh <- 2

# probName <-c("50", "90")
# velTrace <- list()

# for (probI in 1:2) {
#     fileName = paste(dataFileName,probName[probI],".csv", sep = "")
#     velData <- read.csv(fileName, header = FALSE, sep=",")
#     velData <- as.matrix(velData)
#     lengthD <- dim(velData)[2]
#     timePoints <- seq(from = 1, to = lengthD, by = 1)

#       tempM <- colMeans(velData, na.rm=TRUE)
#       velMean <- as.matrix(tempM)
#       velStd <- colSds(velData, na.rm=TRUE)
#       error <- velStd/sqrt(n) # ste
#       # error <- qt(0.975, df=n-1)*velStd/sqrt(n) # 95% CI
#       velLower <- tempM-error
#       velUpper <- tempM+error
#       velTrace[[probI]] <- data.frame(timePoints, velMean, velLower, velUpper)
# }
# # colourCode <- c("#ffb687", "#897eff", "#71cc64") # different hues
# # colourCode <- c("#c6dbef", "#4292c6", "#08306b") # all blue
# colourCode <- c("#d9d9d9", "#737373", "#000000") # all black
# ###

# ## plot velocity traces with 95% CI
# # pdf(paste(folder1, "velocityTraces.pdf", sep = ""))
# p <- ggplot(velTrace[[1]], aes(x = timePoints, y = velMean)) +
#         geom_line(aes(colour = "50%"), size = 1, linetype = "solid") + geom_ribbon(aes(x = timePoints, ymin=velLower, ymax=velUpper, fill = "50%"), alpha = velAlpha) +
#         geom_line(data = velTrace[[2]], aes(x = timePoints, y = velMean, colour = "90%"), size = 1, linetype = "solid") + geom_ribbon(data = velTrace[[2]], aes(x = timePoints, ymin=velLower, ymax=velUpper, fill = "90%"), alpha = velAlpha) +
#         geom_segment(aes_all(c('x', 'y', 'xend', 'yend')), data = data.frame(x = c(0, 0), xend = c(130, 0), y = c(ylimLow, ytickLow), yend = c(ylimLow, ytickHigh)), size = axisLineWidth) +
#         # scale_y_continuous(name = "Anticipatory pursuit velocity (°/s)", limits = c(ylimLow, ylimHigh), breaks=c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
#         scale_y_continuous(name = "P(perceiving right)-P(motion right)", limits = c(ylimLow, ylimHigh), breaks=c(ylimLow, 0, ylimHigh), expand = c(0, 0)) +
#         scale_x_continuous(name = "Perceptual trial number", limits = c(-0, 130), breaks=c(0, 50, 100, 130), expand = c(0, 0)) +
#         scale_colour_manual("Probability of rightward motion",
#                             breaks = c("50%", "90%"),
#                             values = c("50%" = colourCode[1], "90%" = colourCode[3])) +
#         scale_fill_manual("Probability of rightward motion",
#                             breaks = c("50%", "90%"),
#                             values = c("50%" = colourCode[1], "90%" = colourCode[3])) +
#         theme(axis.text=element_text(colour="black"),
#               axis.ticks=element_line(colour="black", size = axisLineWidth),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.border = element_blank(),
#               panel.background = element_blank(),
#               text = element_text(size = textSize),
#               legend.background = element_rect(fill="transparent"),
#               legend.key = element_rect(colour = "transparent", fill = "white"))
# print(p)
# # dev.off()
# ggsave(paste(plotFolder, pdfFileName, sep = ""), width = vtPlotWidth)