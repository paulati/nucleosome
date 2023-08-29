# change this value based on the genome to be analyzed
# cobfigured options: TcruziCL or Sylvio
genome_id <- "TcruziCL" 

if(genome_id == "TcruziCL") {
  config_file_path <- file.path(base_path, 'common', "config.yaml")  
} else if(genome_id == "Sylvio") {
  config_file_path <- file.path(base_path, 'common', "config_sylvio.yaml") 
} else {
  print("Error: genome_id not configured")
}

config <- yaml.load_file(config_file_path)

qc_out_base_path <- config$qc_out_base_path

alignment_summary_out_base_path <- config$alignment_summary_out_base_path

utrme_file_name <- config$utrme_file_name

cores <- config$cores

genomes <- config$genomes

heatmaps <- config$heatmaps

bigwigs <- config$bigwigs

replicates <- config$replicates

repl_config <- function(repl_id) replicates[[repl_id]]



