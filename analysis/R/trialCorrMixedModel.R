# load in libraries:
library(ez)
# library(plyr)
# library(dplyr)
# library(lsr)
library(lme4)
library(lmerTest)
# library(car)
# library(afex)

#### clear environment
rm(list = ls())


# set working directory: 
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R/plots/")

# Read in the data and rename the columns
data <- read.csv("trialDataCorr.csv")

# Define factorial and numerical variables:
data$sub           <- as.factor(data$sub)
data$rdkInternalCon         <- as.factor(data$rdkInternalCon)
data$rdkApertureAngle            <- as.numeric(data$rdkApertureAngle)
data$dirClp           <- as.numeric(data$dirClp)
data$response           <- as.numeric(data$response)

# # scale TTCeye and TTChand:
# data$TTCeye = scale(data$TTCeye, scale = TRUE)
# data$TTChand = scale(data$TTChand, scale = TRUE)
#data$condition = scale(data$condition, scale = TRUE)

# Run Linear Mixed Model using the afex package: https://cran.r-project.org/web/packages/afex/afex.pdf

# add contrasts to the condition factor:
# contrasts(data$condition)         <- MASS::contr.sdif(5)

fm1 <- lmer(response ~ rdkApertureAngle + rdkInternalCon + dirClp + (1+dirClp|sub), data, REML = F)
summary(fm1)

# # run model:
# full_model <- afex::mixed(response ~ rdkApertureAngle + rdkInternalCon + dirClp + rdkInternalCon*dirClp + (1+dirClp|sub), data = data, method = "S")
# summary(full_model)
# #full_model <- afex::mixed(TTChand ~ condition*TTCeye + (TTCeye|subject), data = data, method = "LRT")

# # perform full model comparison:
# afex::nice(full_model)

################################################################################