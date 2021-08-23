library("dplyr")
library("tidyr")
library("ggplot2")
source("C:/Users/marce/Desktop/microbiome-help/microbiome_helper_functions.R")

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
invasion_diversidad_scaled <- as.data.frame(scale(invasion_diversidad))

# Graph heatmap
heatmap(as.matrix(invasion_diversidad_scaled), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"))

# Todas las muestras de invasión

invasion_diversidad_2 <- invasion_diversidad

invasion_diversidad_2["bacteria"] <- row.names(invasion_diversidad)

invasion_diversidad_g <- gather(invasion_diversidad_2, colnames(invasion_diversidad), key = "time", value = "counts")

invasion_diversidad_g$time <- factor(invasion_diversidad_g$time, levels = colnames(invasion_diversidad))

ggplot(invasion_diversidad_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="stack", stat="identity")

ggplot(invasion_diversidad_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity")

# Muestras antes de la invasión
invasion_antes <- select(invasion_diversidad, starts_with("A.i1"))

spcs_names <- colnames(invasion_antes)

invasion_antes["bacteria"] <- row.names(invasion_diversidad)

invasion_antes_g <- gather(invasion_antes, spcs_names, key = "time", value = "counts")

invasion_antes_g$time <- factor(invasion_antes_g$time, levels = colnames(invasion_diversidad))

ggplot(invasion_antes_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))

# Muestras después de la invasión

invasion_despues <- select(invasion_diversidad, starts_with("A.i8"))

spcs_names_d <- colnames(invasion_despues)

invasion_despues["bacteria"] <- row.names(invasion_diversidad)

invasion_despues_g <- gather(invasion_despues, spcs_names_d, key = "time", value = "counts")

invasion_despues_g$time <- factor(invasion_despues_g$time, levels = colnames(invasion_diversidad))

ggplot(invasion_despues_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))


####################################################################################

# Averages of invasion times
####################################################################################
####################################################################################
