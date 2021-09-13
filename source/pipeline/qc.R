library(tools)


check_md5sum <- function(id) {
  
  type <- 'raw'
  
  r_paths <- get_reads_path(id, type)
  
  files <- c(r_paths$r1, r_paths$r2)
  
  md5_value <- md5sum(files)
  
  target_values <- get_md5sum(id)
  
  check_ok <- (md5_value[r_paths$r1] == target_values$r1) &&
    (md5_value[r_paths$r2] == target_values$r2)
    
  return(check_ok)
  
}


# types: clean or raw
quality_control <- function(id, type='raw') {
  
  r_paths <- get_reads_path(id, type)
  
  out_path <- get_qc_out_base_path(id, type)
  
  if(! dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  command <- paste("fastqc", r_paths$r1, r_paths$r2, "-o", out_path, sep =" ")
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)
  
}

clean_repetitive_sequences <- function(id) {
  
  # https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads
  
  cores_param <- paste("-j", cores, sep = " ")
  
  remove_param <- paste("-m", 10, sep = " ")
  
  overrepresented_seqs_file_path <- get_overrepresented_seqs_file_path(id)
  overrepresented_param <- paste0("-a file:", overrepresented_seqs_file_path)
  
  output_paths <- get_clean_data_base_path(id)
  
  output_paths_param <- paste("-o", output_paths$r1, "-p", output_paths$r2, sep = " ")
  
  input_paths <- get_reads_path(id)
  input_paths_param <- paste(input_paths$r1, input_paths$r2, sep = " ")
  
  command <- paste("cutadapt", 
                   cores_param, 
                   remove_param, 
                   overrepresented_param,
                   output_paths_param,
                   input_paths_param,
                   sep=" ")
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)
  
}

