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
alta_diversidad <- select(otu_table, starts_with("A.0"), starts_with("A.4"), starts_with("A.7"), starts_with("A.1"), starts_with("A.2"), starts_with("A.6"))
# Filter only outs with at least one count in one replicate
otu_table_f <- filter_otus_by_counts_col_percent(alta_diversidad, min_count = 20, percentage = 0.2)

#days <- c(4, 7 , 11, 14, 21, 23, 27, 47, 60)
times <- c("A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60.")

# Getting the means of each OTU.
time_mean_diversity <- data.frame()[1:22, ] # Empty dataframe
row.names(time_mean_diversity) <- row.names(otu_table_f)
for (iter_time in times) {
  current_values <- select(otu_table_f, starts_with(iter_time))
  time_mean_diversity[iter_time] <- rowMeans(current_values)
}

time_mean_mat <- data.matrix(time_mean_diversity)

heatmap(time_mean_mat, scale = "row")

## Graph bars
time_mean_div_f_2 <- time_mean_diversity

time_mean_div_f_2["bacteria"] <- row.names(time_mean_diversity)

time_div_g <- gather(time_mean_div_f_2, "A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60.", key = "time", value = "counts")

# Order
time_div_g$time <- factor(time_div_g$time, levels = c("A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60."))

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="stack", stat="identity")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3", "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4", "purple4", "springgreen4", "firebrick3", "gold3", "cyan3", "plum", "mediumspringgreen")

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)


####################################################################################
