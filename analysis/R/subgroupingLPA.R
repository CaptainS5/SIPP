library(ggplot2)
library(tidyLPA)
library(dplyr)

#### clear environment
rm(list = ls())

#### load data
# on Inspiron 13
setwd("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R")
source("pairwise.t.test.with.t.and.df.R")
plotFolder <- ("C:/Users/wuxiu/Documents/PhD@UBC/Lab/3rdYear/SIPP/analysis/R/plots/")
### modify these parameters to plot different conditions
dataFileName <- "bootStrapPerceptualBias.csv"

# for plotting
textSize <- 25
axisLineWidth <- 0.5
dotSize <- 3

data <- read.csv(dataFileName)
subAll <- unique(data$sub)
subTotalN <- length(subAll)

m <- data[, 2:8] %>%
    single_imputation() %>%
    estimate_profiles(2) 

# compare_solutions(statistics = c("AIC", "BIC"))

# pdf(file = "subGrouping2.pdf", width = 4, height = 4)
p <- plot_profiles(m, size = 20)
p <- p + ylim(c(-10, 10)) + scale_colour_manual(values = c("#7cae00", "#c77cff")) + scale_size_manual(values = c(30, 30))
print(p)
ggsave(plot = p, file = "subGrouping2.pdf")
# dev.off()

get_estimates(m) %>% 
    write.csv("subGroupingEstimates2.csv")

get_data(m) %>% 
    write.csv("subGrouping2.csv")
    # tidyr::gather(Class_prob, Probability, contains("CPROB"))