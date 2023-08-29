
get_reads_path <- function(id, type = 'raw') {
  
  object_config <- repl_config(id)
  
  #repl_files_names <- get_files_names(id)
  
  repl_file_path_r1  <- file.path(object_config[[type]]$data_base_path, 
                                  id, 
                                  object_config$repl_file_name$r1)
  
  repl_file_path_r2  <- file.path(object_config[[type]]$data_base_path, 
                                  id, 
                                  object_config$repl_file_name$r2)
  
  result <- list('r1'=repl_file_path_r1, 'r2' = repl_file_path_r2)

  return(result)  
}

get_md5sum <- function(id) {
  object_config <- repl_config(id)
  result <- object_config$repl_file_md5sum
  return(result)
}

get_qc_out_base_path <- function(id, type = 'raw') {
  
  result <- file.path(qc_out_base_path, id, type)
  return(result)
  
}

get_overrepresented_seqs_file_path <- function(id) {
  
  object_config <- repl_config(id)
  result <- file.path(object_config$raw$overrepresented_seqs_base_path, 
                      object_config$raw$overrepresented_seqs_file)
  return(result)
  
}

get_clean_data_base_path <- function(id) {

  object_config <- repl_config(id)  
  clean_data_base_path <- object_config$clean$data_base_path
  
  #repl_files_names <- get_files_names(id)
  
  result <- list('r1' = file.path(clean_data_base_path, id, object_config$repl_file_name$r1), 
                 'r2' = file.path(clean_data_base_path, id, object_config$repl_file_name$r2))
  
  return(result)
  
}

get_genome_file_path <- function(genome_id) {

  genomes_base_path <- genomes[[genome_id]]$genome_base_path
  result <- file.path(genomes_base_path, genomes[[genome_id]]$file_name)
  
  return(result)
  
}

get_index_base_path <- function(genome_id, tool) {

  file_name_parts <- unlist(strsplit(x = genomes[[genome_id]]$file_name, split = "\\."))
  file_name_base <- file_name_parts[1]
  
  base_path <- file.path(genomes[[genome_id]]$index_base_path, tool, file_name_base)
  
  if(! dir.exists(base_path)) {
    dir.create(base_path, showWarnings = FALSE)
  }
  
  result <- file.path(base_path, file_name_base)
  
  return(result)
  
  
}

get_alignmet_out_files_path <- function(repl_id, genome_id, type, tool) {
  
  file_name_parts <- unlist(strsplit(x = genomes[[genome_id]]$file_name, split = "\\."))
  file_name_base <- file_name_parts[1]
  
  result_base_path <-file.path(genomes[[genome_id]]$alignment_base_path, tool, type, file_name_base)
  
  if(! dir.exists(result_base_path)) {
    dir.create(result_base_path, showWarnings = FALSE)
  }
  
  object_config <- repl_config(repl_id)
  
  alignment_file_name <- paste0(object_config$repl_alignment_file_name, '.sam')
  
  summary_file_name <- paste0(object_config$repl_alignment_file_name, '_summary.txt')
  
  log_file_name <- paste0(object_config$repl_alignment_file_name, '.log')
  
  result <- list(
    'alignment' =  file.path(result_base_path, alignment_file_name),
    'summary' = file.path(result_base_path, summary_file_name),
    'log' =  file.path(result_base_path, log_file_name)
  )
  
  return(result)
  
}


