library(universalmotif)
library(MotifDb)
motifs <- filter_motifs(MotifDb, organism = "Athaliana")
write_meme(motifs, "../0_data/B_userData/processed_files/JASPAR_motifs.meme")
