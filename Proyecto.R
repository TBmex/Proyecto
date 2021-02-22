#Libraries
library(fastbaps)
library(ape)

#Loading data OK
fasta.file.name <- "ESP.fq"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)

#Detalle OK
sparse.data <- optimise_prior(sparse.data, type = "baps")

#Running fastbaps "baps.hc" es el archivo a obtener.
baps.hc <- fast_baps(sparse.data)

# Obtener tablas

#Bayesian hierarchical partition
ESP_baps <- multi_res_baps(sparse.data, levels = 4)

#Phylogeny partition, esto trabaja con el arbol pero creo que preferimos datos crudos.
#d <- snp_dist(sparse.data)
#d <- as.dist(d/max(d))
#h <- hclust(d, method="ward.D2")
#multi4 <- multi_level_best_baps_partition(sparse.data, h, levels = 6)

#Function to perform bootstrap replicated of fastbaps
boot.result <- boot_fast_baps(sparse.data)