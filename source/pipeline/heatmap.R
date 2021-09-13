suppressPackageStartupMessages({
  library(yaml)
})

get_command_line_args_base <- function(data_file_path) {
  
  option_file <- paste0("--file=", data_file_path)
  option_type <- paste0("--type=", heatmaps$plot2DO_option_type)
  option_minLength <- paste0("--minLength=", heatmaps$plot2DO_option_minLength)
  option_maxLength <- paste0("--maxLength=", heatmaps$plot2DO_option_maxLength)
  option_genome <- paste0("--genome=", heatmaps$plot2DO_option_genome)
  option_align <- paste0("--align=", heatmaps$plot2DO_option_align)
  command_line_args_base = c(option_file, option_type, option_minLength, option_maxLength, option_genome, option_align)
  
  return(command_line_args_base)
}

init <- function(repl_id, genome_id, type, tool) {

  plot2DO_base_path <- file.path(dirname(base_path), 'tools', 'tcruzi_plot2DO')
  
  setwd(plot2DO_base_path)
  
  # Load default paths from config.yaml
  config <- yaml.load_file("config/config.yaml")
  
  #change global value:
  sourceBasePath <<- file.path(plot2DO_base_path, 'source')

  #change global value:
  outputBasePath_root <- heatmaps$plot2DO_application_paths_output 
  #config$application$paths$output
  outputBasePath <<- file.path(outputBasePath_root, repl_id, tool, type, genome_id)
  if(! dir.exists(outputBasePath)) {
    dir.create(outputBasePath, recursive = TRUE, showWarnings = FALSE)
  }
  
  #change global value:
  annotationsBasePath <<- file.path(plot2DO_base_path, 'annotations')

  configBasePath <<- file.path(plot2DO_base_path, 'config')

  main <- file.path(sourceBasePath, "plot2DO_main.R")  
  
  source(main)
  
}

plot_heatmap <- function(repl_id, genome_id, type, tool) {
  
  init(repl_id, genome_id, type, tool)
  
  paths <- get_alignmet_out_files_path(repl_id, genome_id, type, tool)
  sam_file_path <- paths$alignment
  bam_file_path <- str_replace(sam_file_path, pattern = ".sam", replacement = ".bam")
  
  command_line_args_base <- get_command_line_args_base(bam_file_path)
  sites_file_name <- utrme_file_name
  sites_arg <-  paste0("--sites=", sites_file_name)
  
  command_line_args <- c(command_line_args_base, sites_arg)
  
  Main(command_line_args)
  
}



