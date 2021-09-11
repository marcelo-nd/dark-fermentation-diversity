library("readxl")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("dplyr")

source("C:/Users/marce/Desktop/microbiome-help/microbiome_helper_functions.R")

# Biogas
####################################################################################
# Biogas data exploration.

biogas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "biogas", range = "A68:O129")
colnames(biogas) <- c("tiempo", "fecha", "semana", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

head(biogas)

biogas_by_cases <- gather(biogas, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                          key = "replicate", value = "biogas_ml")
head(biogas_by_cases)

ggplot(biogas_by_cases, aes(tiempo, biogas_ml)) +
  geom_boxplot(aes(group = cut_width(tiempo, 1)), outlier.color = "#BE3A34", alpha = 1) +
  geom_point(size = 1, alpha = 0.2) +
  scale_x_continuous(breaks = biogas$tiempo) +
  ylim(0, 1200) +
  xlab("Time (days)") +
  ylab("Production (biogas mL/L)") +
  ggtitle("Biogas produced (n = 12)") +
  theme_few() +
  #geom_vline(xintercept = c(0, 27, 60), color = "#74AA50", alpha = 0.2 , size = 2) +
  #geom_vline(xintercept = c(0, 1, 2, 4, 7, 9, 11, 14, 16, 18, 21, 23, 25, 26, 27, 47, 53, 60), linetype = "dotted", color = "#00B8DE", alpha = 0.8 , size = 1) +
  geom_vline(xintercept = c(0, 4, 7, 11, 14, 21, 23, 27, 47, 60), linetype = "dotted", color = "#00B8DE", alpha = 0.8 , size = 1) +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=12, angle = 90),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13))

####################################################################################
# AGVs
####################################################################################
# Exploraci�n de datos de producci�n de �cidos grasos volatiles. 12 r�plicas.

# Cargando datos crudos de AGVs. Este DF tiene subr�plicas de medici�n para cada reactor r�plica.
# Tambi�n contiene mediciones r�plica del est�ndard usado y que ayudar� a calcular la concentraci�n en PPM.
# Esta tabla es s�lo para el tiempo 0.
agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "agvs_areas", range = "A6:AA12")

std_values <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "agvs_areas", range = "A130:B136")

agv_means <- get_replicate_means(replicate_data = agvs_areas, id_col_name = "tiempo", id = 0, replicates = 12, subreplicates = 2, standard = TRUE)

agv_ppm_per_day <- get_qnty_from__std(agv_means, id_col = 2, replicates = 12, std_values = std_values)

days <- c(4, 7 , 11, 14, 21, 23, 27, 47, 60)
day_counter <- 1

for(file_upper_row in seq(16, 96, by = 10)){
  file_lower_row <- file_upper_row + 6
  file_range <- paste("A", toString(file_upper_row), ":AA", toString(file_lower_row), sep="")
  agvs_areas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "agvs_areas", range = file_range)
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
  xlab("Tiempo (d�as)") +
  ylab("Concentraci�n AGV's (PPM)") +
  ggtitle("Concentraci�n AGV's (PPM); 12 r�plicas") +
  theme_few() +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=13),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13)) +
  labs(colour = "AGV's", fill = "AGV's")
####################################################################################

# Biogas con invasion
####################################################################################
# Biogas data exploration.

biogas <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "biogas", range = "A68:O129")
colnames(biogas) <- c("tiempo", "fecha", "semana", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
biogas <- select(biogas, "tiempo", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

biogas_by_cases <- gather(biogas, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                          key = "replicate", value = "biogas_ml")

biogas_invasion <- read_xlsx(path = "C:/Users/marce/OneDrive/DiversidadH2/0_datos/alta_diversidad.xlsx", sheet = "biogas_invasion", range = "A13:N22")
biogas_invasion <- select(biogas_invasion, "tiempo", ends_with("i"))[2:9,]
head(biogas_invasion)
biogas_inv_by_cases <- gather(biogas_invasion, "1_1i", "1_3i", "1_4i", "1_5i", "1_6i", "1_8i",
                          key = "replicate", value = "biogas_ml")

biogas_by_cases <- rbind(biogas_by_cases, biogas_inv_by_cases)


ggplot(biogas_by_cases, aes(tiempo, biogas_ml)) +
  geom_boxplot(aes(group = cut_width(tiempo, 1)), outlier.color = "#BE3A34", alpha = 1) +
  geom_point(size = 1, alpha = 0.2) +
  #scale_x_continuous(breaks = c(biogas$tiempo, biogas_invasion$tiempo)) +
  scale_x_continuous(breaks = c(0, 4, 7, 11, 14, 20, 21, 23, 27, 40, 47, 60, 61, 68)) +
  ylim(0, 1200) +
  xlab("Time (days)") +
  ylab("Production (biogas mL/L)") +
  ggtitle("Biogas produced (n = 12)") +
  theme_few() +
  geom_vline(xintercept = c(60), color = "firebrick", alpha = 0.2 , size = 2) +
  #geom_vline(xintercept = c(0, 27, 60), color = "#74AA50", alpha = 0.2 , size = 2) +
  #geom_vline(xintercept = c(0, 1, 2, 4, 7, 9, 11, 14, 16, 18, 21, 23, 25, 26, 27, 47, 53, 60), linetype = "dotted", color = "#00B8DE", alpha = 0.8 , size = 1) +
  geom_vline(xintercept = c(0, 4, 7, 11, 14, 21, 23, 27, 47, 60, 61, 68), linetype = "dotted", color = "#00B8DE", alpha = 0.8 , size = 1) +
  theme(axis.title.x =  element_text(size=18), axis.text.x = element_text(size=12, angle = 45),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size=13))


####################################################################################

# Export data
####################################################################################
# Transforming data and exporting to make it usable in other analyses.
# Preparing AGVs data.

# Se a�ade el nombre de las muestras correctamente
agvs_by_cases["exp"] <- "A"
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
biogas_by_cases2 <- filter(biogas_by_cases2, tiempo == 0 | tiempo == 4 | tiempo == 7| tiempo == 11| tiempo == 14| tiempo == 21| tiempo == 23| tiempo == 27| tiempo == 47| tiempo == 60)
# Adding the name of samples correctly
biogas_by_cases2["exp"] <- "A"
biogas_by_cases2 <- unite(biogas_by_cases2, exp, tiempo, replicate, col = "sample", sep = ".")

# We separate the DF, with samples as key and preserving biogas production value
biogas_by_samples <- spread(biogas_by_cases2, key = "sample", value = "biogas_ml")

biogas_by_samples

# Merging the two data sets
# This data set contains as first column the name of measured parameters. The rest of columns are samples. Rows are the parameters.

metabolite_data <- rbind(agvs_by_samples, biogas_by_samples)

metabolite_data

write.csv(x = metabolite_data, file = "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/metabolite_data.csv", row.names = FALSE)
