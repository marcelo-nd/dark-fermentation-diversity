library("dplyr")
library("corrplot")
library("vegan")
library("tidyr")

source("C:/Users/marce/Documents/Repos/microbiome-help/diversity_data_helper_functions.R")
source("C:/Users/marce/Documents/Repos/microbiome-help/functional_data_helper_functions.R")

# Read Data
####################################################################################
alta_diversidad <- read.csv("C:/Users/marce/OneDrive/Sci/DiversidadH2/Análisis/2_resultados/altad_otu_table.csv", row.names = 1)

metabolite_data <- read.csv("C:/Users/marce/OneDrive/Sci/DiversidadH2/Análisis/2_resultados/metabolite_data.csv", row.names = 1)

####################################################################################

# Correlation Heatmap
####################################################################################
a_div_mat <- t(alta_diversidad)
a_div_mat <- a_div_mat[order(row.names(a_div_mat)), ] # Ordering by row names


metab_mat <- t(metabolite_data)
metab_mat <- metab_mat[order(row.names(metab_mat)), ] # Ordering by row names


adivXmetab <- cor(a_div_mat, metab_mat, method = "pearson")
heatmap(adivXmetab)
corrplot::corrplot(adivXmetab)

####################################################################################

# Correlation Heatmap (with significance)
####################################################################################
library("Hmisc")

# Transforming data to matrices
a_div_mat <- t(alta_diversidad)
a_div_mat <- a_div_mat[order(row.names(a_div_mat)), ] # Ordering by row names

metab_mat <- t(metabolite_data)
metab_mat <- metab_mat[order(row.names(metab_mat)), ] # Ordering by row names

# Getting correlations

adivXmetab <- cor(a_div_mat, metab_mat, method = "pearson")
# Converting rownames to column 1
adivXmetab_df <- as.data.frame(adivXmetab)
adXm <- cbind(rownames(adivXmetab_df), data.frame(adivXmetab_df, row.names=NULL))
colnames(adXm)[1] <- "species"
# Gathering data by species
adXm_g <- gather(adXm, "Acetico", "Butirico", "Isobutirico", "Isovalerico", "Propionico", "Valerico", "Biogas", key = "compound", value = "correlation")

# Getting significance data
adxm.pval <- as.data.frame(rcorr(a_div_mat, metab_mat, type = "pearson")$P[, 21:27][1:20,])
# Converting rownames to column 1
adxm.pval <- cbind(rownames(adxm.pval), data.frame(adxm.pval, row.names=NULL))
colnames(adxm.pval)[1] <- "species"
# Gathering data by species
adXm_pvs_g <- gather(adxm.pval, "Acetico", "Butirico", "Isobutirico", "Isovalerico", "Propionico", "Valerico", "Biogas", key = "compound", value = "p_value")
# Create column of significance labels
adXm_pvs_g$stars <- cut(adXm_pvs_g$p_value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

# Adding stars column to original df
adXm_g["stars"] <- adXm_pvs_g$stars

# Plotting
ggplot(aes(x=compound, y=species, fill=correlation), data=adXm_g) +
  geom_tile() +
  scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") +
  geom_text(aes(label=stars), color="black", size=2) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))

####################################################################################

# CCA
####################################################################################
a_div_mat <- t(alta_diversidad)
a_div_mat <- a_div_mat[order(row.names(a_div_mat)), ] # Ordering by row names

metab_mat <- t(metabolite_data)
metab_mat <- metab_mat[order(row.names(metab_mat)), ] # Ordering by row names
metab_mat_scaled <- apply(metab_mat, 2, scale)

a_div_df <- as.data.frame(a_div_mat)
metab_df <- as.data.frame(metab_mat_scaled)

divXmet.cca <- cca(a_div_df ~ Biogas + Butirico + Isovalerico + Acetico + Isobutirico + Valerico + Propionico, data=metab_df)

plot(divXmet.cca, scaling = 1, type = "none")
text(divXmet.cca, "bp", col="blue", cex=1, scaling = 1)
points(divXmet.cca, "sites", pch=1, col="black", bg="yellow", cex=0.5, scaling = 1)
points(divXmet.cca, "species", pch=21, col="green", bg="yellow", cex=1, scaling = 1)
text(divXmet.cca, "species", col="red", cex=0.8, scaling = 1)

anova(divXmet.cca, perm=1000)
anova(divXmet.cca, by = "axis", perm=1000)
anova(divXmet.cca, by = "term", perm=1000)
#anova(divXmet.cca, by = "mar", perm=1000)

####################################################################################

# CCA (color coded)
####################################################################################
a_div_mat <- t(alta_diversidad)
a_div_mat <- a_div_mat[order(row.names(a_div_mat)), ] # Ordering by row names

metab_mat <- t(metabolite_data)
metab_mat <- metab_mat[order(row.names(metab_mat)), ] # Ordering by row names
metab_mat_scaled <- apply(metab_mat, 2, scale)

a_div_df <- as.data.frame(a_div_mat)
metab_df <- as.data.frame(metab_mat_scaled)

sample_names <- c()
for (sample in row.names(metab_mat)) {
  sample_splt <- strsplit(sample, ".", fixed = TRUE)
  sample_names <- c(sample_names, paste(sample_splt[[1]][1], sample_splt[[1]][2], sep = "."))
}

metab_df["sample_time"] <- sample_names

metab_df$sample_time <- factor(metab_df$sample_time, levels = c("A.0", "A.4", "A.7", "A.11", "A.14", "A.21", "A.23", "A.27", "A.47", "A.60"))

divXmet.cca <- cca(a_div_df ~ Biogas + Butirico + Isovalerico + Acetico + Isobutirico + Valerico + Propionico, data=metab_df)

#colvec <- c("turquoise2", "snow4", "mediumblue", "maroon2", "aquamarine4", "tan3", "darkblue", "slateblue1", "lightpink3", "springgreen4")

colvec <- c("firebrick3", "#cb3040", "#c6235a", "#b92573", "#a33189", "#853e9a", "#7344a1", "#5d49a5", "#434ea8", "#1b51a9")

plot(divXmet.cca, scaling = 1, type = "none")
text(divXmet.cca, "bp", col="goldenrod1", cex=1, scaling = 1)
points(divXmet.cca, "sites", pch=21, col = colvec[metab_df$sample_time], bg=colvec[metab_df$sample_time], cex=1, scaling = 1)
points(divXmet.cca, "species", pch=21, col="forestgreen", bg="forestgreen", cex=1, scaling = 1)
text(divXmet.cca, "species", col="forestgreen", cex=0.7, scaling = 1)
legend("bottomleft", legend = levels(metab_df$sample_time), bty = "n",col = colvec, pch = 21, pt.bg = colvec)

summary(divXmet.cca)
anova(divXmet.cca, perm=1000)
anova(divXmet.cca, by = "axis", perm=1000)
anova(divXmet.cca, by = "term", perm=1000)
#anova(divXmet.cca, by = "mar", perm=1000)

####################################################################################

