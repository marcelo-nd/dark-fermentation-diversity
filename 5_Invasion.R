library("dplyr")
library("tidyr")
library("ggplot2")
source("C:/Users/marce/OneDrive/DiversidadH2/1_scripts/DA_helper_functions.R")

# Read Data
####################################################################################

otu_table <- read.csv("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/otu_table.csv", row.names = 1)

invasion_diversidad <- select(otu_table, starts_with("A.i"))
invasion_diversidad <- filter_otus_by_counts_col_percent(invasion_diversidad, min_count = 20, percentage = 0.20)


# All samples separately
####################################################################################
# Ordering by col names
invasion_diversidad <- invasion_diversidad[ ,order(colnames(invasion_diversidad))]
# Correcting order of columns
invasion_diversidad <- cbind(invasion_diversidad[1], invasion_diversidad[5:12], invasion_diversidad[2:4], invasion_diversidad[13], invasion_diversidad[17:24], invasion_diversidad[14:16])

#scaled by column
invasion_diversidad_scaled <- scale(invasion_diversidad)

# Graph heatmap
heatmap(invasion_diversidad_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"))


invasion_diversidad_2 <- invasion_diversidad

invasion_diversidad_2["bacteria"] <- row.names(invasion_diversidad)

invasion_diversidad_g <- gather(invasion_diversidad_2, colnames(invasion_diversidad), key = "time", value = "counts")

invasion_diversidad_g$time <- factor(invasion_diversidad_g$time, levels = colnames(invasion_diversidad))

ggplot(invasion_diversidad_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="stack", stat="identity")

ggplot(invasion_diversidad_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity")

####################################################################################

# Averages of invasion times
####################################################################################
####################################################################################
