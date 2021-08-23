source("C:/Users/marce/Desktop/update-gg-taxonomy-to-ncbi/update_taxonomy_refseq.R")


# Reading required data.

# Load taxonomy table to be updated with NCBI's taxonomy.

tax_table <- read.table("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/9_tsv/taxonomy.tsv", sep = "\t", header = TRUE)
head(tax_table)

# Load sequences to be analyzed in fasta format. Ids corresponding to taxonomy OTUs.
fasta_data <- readDNAStringSet("C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/5_vsearch_output/clustered_sequences/4d823db1-2e16-49a8-b2e1-08565db0d033/data/dna-sequences.fasta")
# Check fasta sequences
fasta_data

# Reading microbial data base
bacteria_database <- get_database(phyl_group = "bacteria", path = "C:/Users/marce/OneDrive/DiversidadH2/1_scripts")

# set NCBI's entrez api key
Sys.setenv(ENTREZ_KEY = "ed4870836e8f61529227d9176a7c4a994c07")
# Check that key variable is in path.
getkey(service = "entrez")
# Check if blast is in PATH.
Sys.which("blastn")

tax_updated <- update_taxonomy_refseq(taxonomy_table = tax_table, microbial_database = bacteria_database, data_fasta = fasta_data, level = "spcs", phyl_group = "bacteria", update_all = TRUE)

bacteria_tax_table_updated <- tax_updated[[1]]

head(bacteria_tax_table_updated)

id_gb_table <- tax_updated[[2]]

write.table(bacteria_tax_table_updated, file = "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/resultados_h2diversidad_std/9_tsv/taxonomy_updated.tsv", sep = "\t", row.names = FALSE)

write.csv(id_gb_table, file = "C:/Users/marce/OneDrive/DiversidadH2/2_resultados/id_gb_table.csv", row.names = FALSE)
