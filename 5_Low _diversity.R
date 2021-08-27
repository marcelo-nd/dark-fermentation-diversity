library("readxl")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("dplyr")

source("C:/Users/marce/Desktop/microbiome-help/microbiome_helper_functions.R")


# Diversity
####################################################################################

# Time 0
####################################################################################

#Reading OTU table
otu_table <- read_qiime_otu_table("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/Rarefaction/11_table.from_biom_w_taxonomy.txt")

# Selecting only "alta diversidad"
baja_diversidad <- select(otu_table, starts_with("B.0"), starts_with("B.4"), starts_with("B.7"), starts_with("B.1"), starts_with("B.2"), starts_with("B.6"))
# Filter only outs with at least one count in one replicate
baja_table_f <- filter_otus_by_counts_col_percent(baja_diversidad, min_count = 20, percentage = 0.2)

# Selecting only "time zero" data
time_zero_baja <- select(baja_table_f, starts_with("B.0."))

#scaled by column
#time_zero_scaled <- scale(time_zero_baja)

# Graph heatmap
heatmap(as.matrix(time_zero_baja), distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"), scale = "none")

# Euclidean distance
edist <- as.matrix(vegan::vegdist(x = t(time_zero_baja), method="euclidean", binary = TRUE))
heatmap(edist, scale = "none")

# Jaccard simmilarity
jdist <- as.matrix(vegan::vegdist(x = t(time_zero_baja), method="jaccard", binary = TRUE))
heatmap(jdist, scale = "none")

# Dendogram
plot(hclust(as.dist(edist), method="ward.D"))

plot(hclust(as.dist(jdist), method="ward.D"))
####################################################################################

### Time series
####################################################################################

otu_table <- read_qiime_otu_table("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/Rarefaction/11_table.from_biom_w_taxonomy.txt")

# Selecting only "alta diversidad"
baja_diversidad <- select(otu_table, starts_with("B.0"), starts_with("B.4"), starts_with("B.7"), starts_with("B.1"), starts_with("B.2"), starts_with("B.6"))
# Filter only outs with at least one count in one replicate
baja_table_f <- filter_otus_by_counts_col_percent(baja_diversidad, min_count = 20, percentage = 0.2)

#days <- c(4, 7 , 11, 14, 21, 23, 27, 47, 60)
times <- c("B.0.", "B.2.", "B.4.", "B.7.", "B.11.", "B.14.", "B.18.", "B.23.", "B.25.", "B.27.")

# Getting the means of each OTU.
time_mean_diversity <- data.frame()[1:2, ] # Empty dataframe
row.names(time_mean_diversity) <- row.names(baja_table_f)
for (iter_time in times) {
  print(iter_time)
  current_values <- select(baja_table_f, starts_with(iter_time))
  time_mean_diversity[iter_time] <- rowMeans(current_values)
}

time_mean_mat <- data.matrix(time_mean_diversity)

heatmap(time_mean_mat, scale = "none")

## Graph bars
time_mean_div_f_2 <- time_mean_diversity

time_mean_div_f_2["bacteria"] <- row.names(time_mean_diversity)

time_div_g <- gather(time_mean_div_f_2, "B.0.", "B.2.", "B.4.", "B.7.", "B.11.", "B.14.", "B.18.", "B.23.", "B.25.", "B.27.", key = "time", value = "counts")

# Order
time_div_g$time <- factor(time_div_g$time, levels = c("B.0.", "B.2.", "B.4.", "B.7.", "B.11.", "B.14.", "B.18.", "B.23.", "B.25.", "B.27."))

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="stack", stat="identity")

cbbPalette <- c("darkolivegreen3", "lightskyblue4")

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)


####################################################################################
