
## When data is downloaded from CIS-BP, a TF_information.txt file will be generated
TF_Information_file <- "/Users/TedCCLeung/Documents/Projects/Photoperiod/2_analysis/0_data/A_pipelineData/CISBP/TF_Information.txt"
## When data is downloaded from CIS-BP, a directory with the positional weight matrices will be generated
PWM_dir <- "/Users/TedCCLeung/Documents/Projects/Photoperiod/2_analysis/0_data/A_pipelineData/CISBP/pwms_all_motifs/"
## Output MEME file
out_file <- "/Users/TedCCLeung/Documents/Projects/Packages/MotifAnalysis/raw-data/CISBP/CISBP_motifs.meme"

## Read in the motif data and only retrieve motifs with direct evidence
## Motifs with indirect evidence are redundant
df <- read.table(TF_Information_file, sep = "\t", header = TRUE)
df <- df[df$TF_Status == "D", ]


## Retrieve the motif files that were selected
motif_ids <- df$Motif_ID
motif_files <- paste0(PWM_dir, motif_ids, ".txt")

## There are some motif files that are empty
## First, we remove them from the files we are going to read
file_check <- sapply(motif_files, function(f){ifelse(nrow(read.table(file = f, sep = "\t")) == 1, FALSE, TRUE)}, USE.NAMES = FALSE)
motif_files_checked <- motif_files[file_check]

## Load the motifs
motifs <- lapply(motif_files_checked, universalmotif::read_cisbp) %>% universalmotif::to_df()
motifs$name <- df$Motif_ID[file_check]
motifs$altname <- df$TF_Name[file_check]
motifs$family <- df$Family_Name[file_check]
motifs$organism <- df$TF_Species[file_check]
CISBP_motifset <- motifs %>% universalmotif::to_list()

write_meme(CISBP_motifset, file = out_file, overwrite = TRUE)
use_data(CISBP_motifset)
use_data_raw()
