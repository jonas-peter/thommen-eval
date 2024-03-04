setwd("/Users/ariel/src/unibe-msb/thommen-eval")

source("helper.R")

################################################################################
# Read data form CSV
################################################################################
# Read data of WP1 of the experimentally tested samples

data <- read.csv(paste0('data_generated/df_st.csv', sep=''))

dataE <- filter(data, group == "Element")
dataN <- filter(data, group == "Nevo")


################################################################################
# Evaluation ST
################################################################################

# all results of the tested samples
bvtvgE <- round(dataE$BVTV, 2)
bvtvgN <- round(dataN$BVTV, 2)

yE <- dataE$ST
yN <- dataN$ST


################################################################################
# Ancova Analysis
################################################################################

# IT
group <- "group"
dependent <- "ST"
independent <- "BVTV"

ancova <- ancova_analysis(data, group, dependent, independent)

dir.create("plots")
jpeg(paste0('plots/emmenas_', dependent, '.jpg'), units = "mm", width = 150, height = 130, res = 400)
ancova[10]
dev.off()


################################################################################
# # Stiffness With respect to BV/TV:
################################################################################

titel <- 'Stiffness vs. BV/TV'

dir.create("plots")
jpeg(paste0('plots/EXP_ST.jpg'), units = "mm", width = 150, height = 130, res = 400)

plot2setts(bvtvgE, yE, bvtvgN, yN,
           'BV/TV []', 'Stiffness [N/mm]',
           'red', 'blue',
           'Element', 'Nevo',
           titel)
dev.off()


################################################################################
# # Values for descriptive stat for Stiffness:
################################################################################

cat("\nDescriptive Stat Stiffness\n\n")
# mean:
meanE <- mean(yE)
cat(paste0("mean ELEMENT:\t", meanE, "\n"))
meanN <- mean(yN)
cat(paste0("mean NEVO:\t\t", meanN, "\n"))
# StD:
sdE <- sd(yE)
cvE <- sdE/meanE
cat(paste0("coefficient of variation ELEMENT:\t", cvE, "\n"))
sdN <- sd(yN)
cvN <- sdN/meanN
cat(paste0("coefficient of variation NEVO:\t\t", cvN, "\n"))

################################################################################
