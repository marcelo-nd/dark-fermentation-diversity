library("readxl")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("dplyr")

source("C:/Users/marce/Desktop/microbiome-help/microbiome_helper_functions.R")

# Metabolites
####################################################################################
# Biogas
####################################################################################
biogas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "biogas", range = "A36:O64")
colnames(biogas) <- c("tiempo", "fecha", "semana", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

head(biogas)

biogas_by_cases <- gather(biogas, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                          key = "replicate", value = "biogas_ml")
head(biogas_by_cases)

ggplot(biogas_by_cases, aes(tiempo, biogas_ml)) +
  geom_boxplot(aes(group = cut_width(tiempo, 1)), outlier.color = "#BE3A34", alpha = 1) +
  geom_point(size = 1, alpha = 0.2) +
  scale_x_continuous(breaks = biogas$tiempo) +
  ylim(500, 2000) +
  xlab("Time (days)") +
  ylab("Production (biogas mL/L)") +
  ggtitle("Biogas produced (n = 12)") +
  theme_few() +
  geom_vline(xintercept = c(0, 2, 4, 7, 11, 14, 18, 23, 25, 27), linetype = "dotted", color = "#00B8DE", alpha = 0.8 , size = 1) +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=7, angle = 45),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13))

####################################################################################

# Biogas invasion
####################################################################################
biogas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "biogas", range = "A36:O64")
colnames(biogas) <- c("tiempo", "fecha", "semana", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
biogas <- select(biogas, "tiempo", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

biogas_by_cases <- gather(biogas, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                          key = "replicate", value = "biogas_ml")


biogas_invasion <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "biogas_invasion", range = "A13:N22")
biogas_invasion <- select(biogas_invasion, "tiempo", ends_with("i"))[2:9,]
head(biogas_invasion)
biogas_inv_by_cases <- gather(biogas_invasion, "3i", "4i", "5i", "6i", "8i", "11i",
                              key = "replicate", value = "biogas_ml")

biogas_by_cases <- rbind(biogas_by_cases, biogas_inv_by_cases)


head(biogas_by_cases)

ggplot(biogas_by_cases, aes(tiempo, biogas_ml)) +
  geom_boxplot(aes(group = cut_width(tiempo, 1)), outlier.color = "#BE3A34", alpha = 1) +
  geom_point(size = 1, alpha = 0.2) +
  scale_x_continuous(breaks = c(0, 2, 4, 7, 10, 11, 14, 18, 20, 23, 25, 27, 28, 30, 35)) +
  ylim(500, 2000) +
  xlab("Time (days)") +
  ylab("Production (biogas mL/L)") +
  ggtitle("Biogas produced (n = 12)") +
  theme_few() +
  geom_vline(xintercept = c(27), color = "firebrick", alpha = 0.2 , size = 2) +
  geom_vline(xintercept = c(0, 2, 4, 7, 11, 14, 18, 23, 25, 27, 28, 35), linetype = "dotted", color = "#00B8DE", alpha = 0.8 , size = 1) +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=12, angle = 45),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13))

####################################################################################

# AGVs
####################################################################################
agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "agvs_areas", range = "A6:AA12")

std_values <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "agvs_areas", range = "A130:B136")

agv_means <- get_replicate_means(replicate_data = agvs_areas, id_col_name = "tiempo", id = 0, replicates = 12, subreplicates = 2, standard = TRUE)

agv_ppm_per_day <- get_qnty_from__std(agv_means, id_col = 2, replicates = 12, std_values = std_values)

days <- c(2, 4, 7, 11, 14, 18, 23, 25, 27)

day_counter <- 1

for(file_upper_row in seq(16, 96, by = 10)){
  print(days[day_counter])
  file_lower_row <- file_upper_row + 6
  file_range <- paste("A", toString(file_upper_row), ":AA", toString(file_lower_row), sep="")
  agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "agvs_areas", range = file_range)
  agv_means <- get_replicate_means(replicate_data = agvs_areas, id_col_name = "tiempo", id = days[day_counter], replicates = 12, subreplicates = 2, standard = TRUE)
  agv_ppm <- get_qnty_from__std(agv_means, id_col = 2, replicates = 12, std_values = std_values)
  agv_ppm_per_day <- rbind(agv_ppm_per_day, agv_ppm)
  day_counter <- day_counter + 1
}

colnames(agv_ppm_per_day) <- c("Compuesto", "tiempo", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

head(agv_ppm_per_day)

agvs_by_cases <- gather(agv_ppm_per_day, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                        key = "replicate", value = "agv_ppm")
head(agvs_by_cases)

ggplot(data = agvs_by_cases, aes(factor(tiempo), agv_ppm)) +
  geom_boxplot(aes(fill = Compuesto), alpha = 0.5, outlier.shape=NA) +
  geom_point(aes(colour = Compuesto), alpha = 1, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75)) +
  ylim(0, 1200) +
  xlab("Tiempo (dias)") +
  ylab("Concentracion AGV's (PPM)") +
  ggtitle("Concentracion AGV's (PPM); 12 replicas") +
  theme_few() +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13)) +
  labs(colour = "AGV's", fill = "AGV's")
####################################################################################

# AGVs (includes invasion experiment concentration data data)
####################################################################################
# Exploraciï¿½n de datos de producciï¿½n de ï¿½cidos grasos volatiles. 12 rï¿½plicas.

# Cargando datos crudos de AGVs. Este DF tiene subrï¿½plicas de mediciï¿½n para cada reactor rï¿½plica.
# Tambiï¿½n contiene mediciones rï¿½plica del estï¿½ndard usado y que ayudarï¿½ a calcular la concentraciï¿½n en PPM.
# Esta tabla es sï¿½lo para el tiempo 0.
agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "agvs_areas", range = "A6:AA12")

std_values <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "agvs_areas", range = "A126:B132")

agv_means <- get_replicate_means(replicate_data = agvs_areas, id_col_name = "tiempo", id = 0, replicates = 12, subreplicates = 2, standard = TRUE)

agv_ppm_per_day <- get_qnty_from__std(agv_means, id_col = 2, replicates = 12, std_values = std_values)

days <- c(2, 4, 7, 11, 14, 18, 23, 25, 27)
day_counter <- 1

for(file_upper_row in seq(16, 96, by = 10)){
  file_lower_row <- file_upper_row + 6
  file_range <- paste("A", toString(file_upper_row), ":AA", toString(file_lower_row), sep="")
  agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "agvs_areas", range = file_range)
  agv_means <- get_replicate_means(replicate_data = agvs_areas, id_col_name = "tiempo", id = days[day_counter], replicates = 12, subreplicates = 2, standard = TRUE)
  agv_ppm <- get_qnty_from__std(agv_means, id_col = 2, replicates = 12, std_values = std_values)
  agv_ppm_per_day <- rbind(agv_ppm_per_day, agv_ppm)
  day_counter <- day_counter + 1
}

colnames(agv_ppm_per_day) <- c("Compuesto", "tiempo", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

head(agv_ppm_per_day)

agvs_by_cases <- gather(agv_ppm_per_day, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                        key = "replicate", value = "agv_ppm")
head(agvs_by_cases)

# Read invasion experiment AGVs concentration data
agvs_areas_inv <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "agvs_areas", range = "A106:AA112")

agv_means_inv <- get_replicate_means(replicate_data = agvs_areas_inv, id_col_name = "tiempo", id = 28, replicates = 12, subreplicates = 2, standard = TRUE)

agv_ppm_per_day_inv <- get_qnty_from__std(agv_means_inv, id_col = 2, replicates = 12, std_values = std_values)

days <- c(35)
day_counter <- 1

for(file_upper_row in seq(116, 125, by = 10)){
  file_lower_row <- file_upper_row + 6
  file_range <- paste("A", toString(file_upper_row), ":AA", toString(file_lower_row), sep="")
  agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/baja_diversidad.xlsx", sheet = "agvs_areas", range = file_range)
  agv_means <- get_replicate_means(replicate_data = agvs_areas, id_col_name = "tiempo", id = days[day_counter], replicates = 12, subreplicates = 2, standard = TRUE)
  agv_ppm <- get_qnty_from__std(agv_means, id_col = 2, replicates = 12, std_values = std_values)
  agv_ppm_per_day_inv <- rbind(agv_ppm_per_day_inv, agv_ppm)
  day_counter <- day_counter + 1
}

agv_ppm_per_day_inv <- select(agv_ppm_per_day_inv, "Compuesto", "tiempo", "R2_", "R4_","R6_", "R8_", "R10", "R12")

colnames(agv_ppm_per_day_inv) <- c("Compuesto", "tiempo", "2", "4", "6", "8", "10", "12")

agvs_by_cases_inv <- gather(agv_ppm_per_day_inv, "2", "4", "6", "8", "10", "12",
                            key = "replicate", value = "agv_ppm")
head(agvs_by_cases_inv)
# Joining initial experiment with invasion experiment data

agvs_by_cases <- rbind(agvs_by_cases, agvs_by_cases_inv)


ggplot(data = agvs_by_cases, aes(factor(tiempo), agv_ppm)) +
  geom_boxplot(aes(fill = Compuesto), alpha = 0.5, outlier.shape=NA) +
  geom_point(aes(colour = Compuesto), alpha = 1, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75)) +
  ylim(0, 1200) +
  xlab("Time (days)") +
  ylab("VFAs Concentration(PPM)") +
  ggtitle("VFAs Concentration(PPM); 12 replicates") +
  geom_vline(xintercept = c(10), color = "firebrick", alpha = 0.2 , size = 2) +
  theme_few() +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13)) +
  labs(colour = "VFAs", fill = "VFAs")
####################################################################################

# Export metabolite data
####################################################################################
# Transforming data and exporting to make it usable in other analyses.
# Preparing AGVs data.

# Se añade el nombre de las muestras correctamente
agvs_by_cases["exp"] <- "B"
agvs_by_cases <- unite(agvs_by_cases, exp, tiempo, replicate, col = "sample", sep = ".")
# We separate the DF, with samples as key and preserving agv concentration value
agvs_by_samples <- spread(agvs_by_cases, key = "sample", value = "agv_ppm")

agvs_by_samples

# Preparing Biogas Data
biogas_by_cases2 <- biogas_by_cases

# Adding temporary column 
biogas_by_cases2["Compuesto"] <- "Biogas"
# Selecting relevant columns only
biogas_by_cases2 <- select(biogas_by_cases2, Compuesto, tiempo, replicate, biogas_ml)
# Choosing appropriate samples only
biogas_by_cases2 <- filter(biogas_by_cases2, tiempo == 0 | tiempo == 2 | tiempo == 4| tiempo == 7| tiempo == 11| tiempo == 14| tiempo == 18| tiempo == 23| tiempo == 25| tiempo == 27)
# Adding the name of samples correctly
biogas_by_cases2["exp"] <- "B"
biogas_by_cases2 <- unite(biogas_by_cases2, exp, tiempo, replicate, col = "sample", sep = ".")

# We separate the DF, with samples as key and preserving biogas production value
biogas_by_samples <- spread(biogas_by_cases2, key = "sample", value = "biogas_ml")

biogas_by_samples

# Merging the two data sets
# This data set contains as first column the name of measured parameters. The rest of columns are samples. Rows are the parameters.

metabolite_data <- rbind(agvs_by_samples, biogas_by_samples)

metabolite_data

write.csv(x = metabolite_data, file = "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/metabolite_data_baja.csv", row.names = FALSE)

####################################################################################

# Diversity-Time 0
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

# Diversity-Time series
####################################################################################

otu_table <- read.csv("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/otu_table.csv", row.names = 1)

# Selecting only "alta diversidad"
baja_diversidad <- select(otu_table, starts_with("B.0"), starts_with("B.4"), starts_with("B.7"), starts_with("B.1"), starts_with("B.2"), starts_with("B.6"))
# Filter only outs with at least one count in one replicate
baja_table_f <- filter_otus_by_counts_col_percent(baja_diversidad, min_count = 20, percentage = 0.2)

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

cbbPalette <- c("#0072B2", "darkorange3")

ggplot(time_div_g, aes(x=time, y=counts, fill=bacteria)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)


####################################################################################

# Diversity X Metabolites
####################################################################################
library("dplyr")
library("corrplot")
library("vegan")
library("tidyr")
library("Hmisc")

source("C:/Users/marce/Desktop/microbiome-help/microbiome_helper_functions.R")

otu_table <- read.csv("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/otu_table.csv", row.names = 1)

metabolite_data <- read.csv("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/metabolite_data_baja.csv", row.names = 1)

baja_diversidad <- select(otu_table, starts_with("B.0"), starts_with("B.2"), starts_with("B.4"), starts_with("B.7"), starts_with("B.1"))
baja_diversidad <- filter_otus_by_counts_col_percent(baja_diversidad, min_count = 20, percentage = 0.20)

write.csv(baja_diversidad, "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/baja_diversidad_otu_table.csv", row.names =  TRUE)

# Correlation Heatmap (with significance)

# Transforming data to matrices
b_div_mat <- t(baja_table_f)
b_div_mat <- b_div_mat[order(row.names(b_div_mat)), ] # Ordering by row names

metab_mat <- t(metabolite_data)
metab_mat <- metab_mat[order(row.names(metab_mat)), ] # Ordering by row names

# Getting correlations

adivXmetab <- cor(b_div_mat, metab_mat, method = "pearson")
# Converting rownames to column 1
adivXmetab_df <- as.data.frame(adivXmetab)
adXm <- cbind(rownames(adivXmetab_df), data.frame(adivXmetab_df, row.names=NULL))
colnames(adXm)[1] <- "species"
# Gathering data by species
adXm_g <- gather(adXm, "Acetico", "Butirico", "Isobutirico", "Isovalerico", "Propionico", "Valerico", "Biogas", key = "compound", value = "correlation")

# Getting significance data
adxm.pval <- as.data.frame(rcorr(b_div_mat, metab_mat, type = "pearson")$P[, 3:9][1:2,])
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

# Invasion
####################################################################################
library("dplyr")
library("tidyr")
library("ggplot2")
source("C:/Users/marce/Desktop/microbiome-help/microbiome_helper_functions.R")

otu_table <- read.csv("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/otu_table.csv", row.names = 1)
invasion_baja <- select(otu_table, starts_with("B.i"))
invasion_baja <- filter_otus_by_counts_col_percent(invasion_baja, min_count = 20, percentage = 0.20)
# Ordering by col names
invasion_baja <- invasion_baja[ ,order(colnames(invasion_baja))]
# Correcting order of columns
invasion_baja <- cbind(invasion_baja[1], invasion_baja[5:12], invasion_baja[2:4], invasion_baja[13], invasion_baja[17:24], invasion_baja[14:16])

# Muestras antes de la invasión
#invasion_antes <- select(invasion_baja, starts_with("B.i1"))

invasion_antes <- select(invasion_baja, "B.i1.2", "B.i1.4", "B.i1.6", "B.i1.8", "B.i1.10", "B.i1.12")

spcs_names <- colnames(invasion_antes)

invasion_antes["bacteria"] <- row.names(invasion_baja)

invasion_antes_g <- gather(invasion_antes, spcs_names, key = "time", value = "counts")

invasion_antes_g$time <- factor(invasion_antes_g$time, levels = colnames(invasion_baja))

cbbPalette <- c("#0072B2", "darkorange3","lightskyblue4")

ggplot(invasion_antes_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)

# Muestras después de la invasión

invasion_despues <- select(invasion_baja, "B.i8.2", "B.i8.4", "B.i8.6", "B.i8.8", "B.i8.10", "B.i8.12")

spcs_names_d <- colnames(invasion_despues)

invasion_despues["bacteria"] <- row.names(invasion_baja)

invasion_despues_g <- gather(invasion_despues, spcs_names_d, key = "time", value = "counts")

invasion_despues_g$time <- factor(invasion_despues_g$time, levels = colnames(invasion_baja))

cbbPalette <- c("#0072B2", "darkorange3","lightskyblue4")

ggplot(invasion_despues_g, aes(x=time, y=counts, fill=bacteria)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=cbbPalette)


####################################################################################
