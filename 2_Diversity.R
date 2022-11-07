library("tidyr")
library("dplyr")
library("vegan")
library("ggplot2")


source("C:/Users/marce/Documents/Repos/microbiome-help/diversity_data_helper_functions.R")
source("C:/Users/marce/Documents/Repos/microbiome-help/functional_data_helper_functions.R")

# Reading data
####################################################################################
# Reading otu table from dada2 output. Aggregated by species level.
biom_file <- "C:/Users/marce/OneDrive/Sci/DiversidadH2/Análisis/2_resultados/diversidad_h2_dada/9_table_w_tax.biom"
otu_table <- get_otu_table_dada(biom_file = biom_file, starting_col = 16, level = "Species")

alta_diversidad <- select(otu_table, starts_with("A.0"), starts_with("A.4"), starts_with("A.7"), starts_with("A.1"), starts_with("A.2"), starts_with("A.6"))

alta_diversidad_f <- filter_otus_by_counts_col_percent(alta_diversidad, min_count = 20, percentage = 0.2)

# Reading otu table from deblur output.
#otu_table <- read_qiime_otu_table("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/Rarefaction2/11_table.from_biom_w_taxonomy.txt")

####################################################################################

# Rarefaction curves
####################################################################################
#Reading OTU table

alta_diversidad_t <- t(alta_diversidad_f)

raremax <- min(rowSums(alta_diversidad_t))

rarecurve(alta_diversidad_t, step = 1, sample = raremax,label = FALSE)

####################################################################################

# T-0 Analyses
####################################################################################
# Selecting only "time zero" data
time_zero <- select(alta_diversidad_f, starts_with("A.0."))

#scaled by column
time_zero_scaled <- scale(time_zero)

# Graph heatmap
heatmap(time_zero_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"))

# Euclidean distance
edist <- as.matrix(vegan::vegdist(x = t(time_zero), method="euclidean", binary = TRUE))
heatmap(edist, scale = "none")

# Jaccard simmilarity
jdist <- as.matrix(vegan::vegdist(x = t(time_zero), method="jaccard", binary = TRUE))
corrplot::corrplot(jdist, order = 'hclust', addrect = 2)
corrplot::corrplot.mixed(jdist, order = 'hclust', addrect = 2)
heatmap(jdist, scale = "none")

# Bray-Curtis index is based on abundance data, while the Sorensen index is based on presence/absence data
bdist <- as.matrix(vegan::vegdist(x = t(time_zero_f), method="bray"))
corrplot::corrplot(bdist, order = 'hclust', addrect = 2)
corrplot::corrplot.mixed(bdist, order = 'hclust', addrect = 2)
heatmap(bdist, scale = "none")

# Dendogram
plot(hclust(as.dist(edist), method="ward.D"))

plot(hclust(as.dist(jdist), method="ward.D"))

####################################################################################

# T-60 Analyses
####################################################################################
#Reading OTU table
#otu_table <- read_qiime_otu_table("C:/Users/marce/OneDrive/Sci/DiversidadH2/Análisis/2_resultados/resultados_h2/diversidad_std/11_table.from_biom_w_taxonomy.txt")

# Selecting only "alta diversidad"
alta_diversidad <- select(otu_table, starts_with("A.0"), starts_with("A.4"), starts_with("A.7"), starts_with("A.1"), starts_with("A.2"), starts_with("A.6"))
# Filter only outs with at least one count in one replicate
otu_table_f <- filter_otus_by_counts_col_percent(alta_diversidad, min_count = 20, percentage = 0.2)

# Selecting only "time 60" data
time_60 <- select(otu_table_f, starts_with("A.60."))

#scaled by column
time_60_scaled <- scale(time_60)

# Graph heatmap
heatmap(time_60_scaled, distfun = function(x) dist(x, method="euclidian"), hclustfun = function(x) hclust(x, method="ward.D"))

# Euclidean distance
edist <- as.matrix(vegan::vegdist(x = t(time_60_f), method="euclidean", binary = TRUE))
heatmap(edist, scale = "none")

# Jaccard similarity
jdist <- as.matrix(vegan::vegdist(x = t(time_60_f), method="jaccard", binary = TRUE))
corrplot::corrplot(jdist, order = 'hclust', addrect = 2)
corrplot::corrplot.mixed(jdist, order = 'hclust', addrect = 2)
heatmap(jdist, scale = "none")

# Bray-Curtis index is based on abundance data, while the Sorensen index is based on presence/absence data
bdist <- as.matrix(vegan::vegdist(x = t(time_zero_f), method="bray"))
corrplot::corrplot(bdist, order = 'hclust', addrect = 2)
corrplot::corrplot.mixed(bdist, order = 'hclust', addrect = 2)
heatmap(bdist, scale = "none")

####################################################################################

# Time series analyses
####################################################################################
#days <- c(4, 7 , 11, 14, 21, 23, 27, 47, 60)
times <- c("A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60.")

# Getting the means of each OTU.
time_mean_diversity <- data.frame()[1:20, ] # Empty dataframe
row.names(time_mean_diversity) <- row.names(alta_diversidad_f)
for (iter_time in times) {
  current_values <- select(alta_diversidad_f, starts_with(iter_time))
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

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3", "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4", "purple4", "springgreen4", "firebrick3", "gold3", "cyan3", "plum", "mediumspringgreen", "blue", "red")

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)

####################################################################################

# Time series analyses (rare species only)
####################################################################################
#Reading OTU table
otu_table <- read_qiime_otu_table("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/11_table.from_biom_w_taxonomy.txt")
#days <- c(4, 7 , 11, 14, 21, 23, 27, 47, 60)
times <- c("A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60.")

# Selecting only "alta diversidad"
alta_diversidad <- select(otu_table, starts_with("A.0"), starts_with("A.4"), starts_with("A.7"), starts_with("A.1"), starts_with("A.2"), starts_with("A.6"))
# Filter only outs with at least one count in one replicate
otu_table_f <- filter_otus_by_counts_col_percent(alta_diversidad, min_count = 20, percentage = 0.2)

time_mean_diversity <- data.frame()[1:20, ]
row.names(time_mean_diversity) <- row.names(otu_table_f)
for (iter_time in times) {
  current_values <- select(otu_table_f, starts_with(iter_time))
  time_mean_diversity[iter_time] <- rowMeans(current_values)
}

time_div_rare <- as.data.frame(t(time_mean_diversity))

time_div_rare <- as.data.frame(select(time_div_rare, -starts_with(c("Clostridium pasteurianum_1", "Ethanoligenens harbinense_1","Enterococcus dispar_1", "Clostridium butyricum_1", "Citrobacter amalonaticus_1"))))

time_div_rare_mat <- data.matrix(t(time_div_rare))

heatmap(time_div_rare_mat, scale = "row")

# pending do it with ggplot
#ggplot(time_div_rare_mat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()


time_mean_div_f_2 <- as.data.frame(t(time_div_rare))

time_mean_div_f_2["bacteria"] <- row.names(time_mean_div_f_2)

time_div_g <- gather(time_mean_div_f_2, "A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60.", key = "time", value = "counts")

# Order
time_div_g$time <- factor(time_div_g$time, levels = c("A.0.", "A.4.", "A.7.", "A.11.", "A.14.", "A.21.", "A.23.", "A.27.", "A.47.", "A.60."))

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="stack", stat="identity")

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity")

####################################################################################

# Species correlation heatmap
####################################################################################
#Reading OTU table
#library("corrplot")

otu_table <- read_qiime_otu_table("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/Rarefaction2/11_table.from_biom_w_taxonomy.txt")

alta_diversidad <- select(otu_table, starts_with("A.0"), starts_with("A.4"), starts_with("A.7"), starts_with("A.1"), starts_with("A.2"), starts_with("A.6"))

otu_table_f <- filter_otus_by_counts_col_percent(alta_diversidad, min_count = 20, percentage = 0.2)

otu_table_t <- t(otu_table_f)

otuXotu <- cor(otu_table_t, method = "pearson")

heatmap(otuXotu, scale = "none")
corrplot::corrplot(otuXotu, order = 'hclust')
corrplot::corrplot.mixed(otuXotu)

### Adding significance stars on corrplot

# Performing significance tests
testRes = cor.mtest(t(otu_table_f), conf.level = 0.95)

# Making corrplot with stars
corrplot(otuXotu, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5, 
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

####################################################################################

# Species co-occurrence analyses
####################################################################################
library(cooccur)

# Transforming abundance data to presence/abscence
otu_table_pa <- vegan::decostand(alta_diversidad_f, method = "pa")

# Infering co-ocurrences
cooccur.otus <- cooccur(otu_table_pa,
                           type = "spp_site",
                           spp_names = TRUE)

summary(cooccur.otus)
plot(cooccur.otus)

### Inferring a network from cooccur results

# Getting only the significant interactions
co <- print(cooccur(otu_table_pa, type = "spp_site", spp_names = TRUE))

# Create a data frame of the nodes in the network. 
nodes <- data.frame(id = 1:nrow(otu_table_pa),
                    label = rownames(otu_table_pa),
                    color = "#606482",
                    shadow = TRUE) 

# Create an edges dataframe from the significant pairwise co-occurrences.
edges <- data.frame(from = co$sp1, to = co$sp2,
                    color = ifelse(co$p_lt <= 0.05,
                                   "#B0B2C1", "#3C3F51"),
                                   dashes = ifelse(co$p_lt <= 0.05, TRUE, FALSE))

# Plotting network
library(visNetwork)
visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(layout = "layout_with_kk")



####################################################################################

#Export Data
write.csv(otu_table, "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/otu_table.csv")

write.csv(alta_diversidad_f, "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/altad_otu_table.csv")

write.csv(time_mean_diversity, "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/altad_means_otu_table.csv")

write.table(time_mean_diversity, "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/altad_means_otu_table.txt", sep = "\t")
