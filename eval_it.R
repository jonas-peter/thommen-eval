setwd("/Users/ariel/src/unibe-msb/thommen-eval")

source("helper.R")

################################################################################
# Read data form CSV
################################################################################
# Read data of WP1 of the experimentally tested samples

data <- read.csv(paste0('data_generated/df_it.csv', sep=''))

dataE <- filter(data, group == "Element")
dataN <- filter(data, group == "Nevo")

################################################################################
# Ancova Analysis
################################################################################

# IT
group <- "group"
dependent <- "IT"
independent <- "BVTV"

ancova <- ancova_analysis(data, group, dependent, independent)

dir.create("plots")
jpeg(paste0('plots/emmenas_', dependent, '.jpg'), units = "mm", width = 150, height = 130, res = 400)
ancova[10]
dev.off()


################################################################################
# Evaluation IT
################################################################################

# all results of the tested samples
bvtvgE <- round(dataE$BVTV, 2)
bvtvgN <- round(dataN$BVTV, 2)

yE <- dataE$IT
yN <- dataN$IT


################################################################################
# # Insertion Torque With respect to BV/TV:
################################################################################

titel <- 'Insertion Torque vs. BV/TV'

dir.create("plots")
jpeg(paste0('plots/EXP_IT.jpg'), units = "mm", width = 150, height = 130, res = 400)

 plot2setts(bvtvgE, yE, bvtvgN, yN,
            'BV/TV []', 'Implantation Torque [Nmm]',
            'red', 'blue',
            'Element', 'Nevo',
            titel)
 dev.off()


################################################################################
 # # Values for descriptive stat for Insertion Torque:
################################################################################

cat("\nDescriptive Stat Insertion Torque\n\n")
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
