build_hisat2_index <- function(genome_id) {
  
  # http://daehwankimlab.github.io/hisat2/manual/#:~:text=The%20hisat2-build%20indexer
  
  genome_file_path <- get_genome_file_path(genome_id)
  
  fasta_param <- paste('-f', genome_file_path, sep = ' ')
  out_base_param <- get_index_base_path(genome_id, 'hisat2')
  
  command <- paste('hisat2-build', fasta_param, out_base_param, sep = " ")
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)
  
}

build_bowtie2_index <- function(genome_id) {
  
  # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
  
  genome_file_path <- get_genome_file_path(genome_id)
  
  fasta_param <- paste('-f', genome_file_path, sep = ' ')
  out_base_param <- get_index_base_path(genome_id, 'bowtie2')
  
  command <- paste('bowtie2-build', fasta_param, out_base_param, sep = " ")
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)

  }

align_hisat2 <- function(repl_id, genome_id, type) {

  # http://daehwankimlab.github.io/hisat2/manual/
  # https://www.biostars.org/p/191995/
  # --summary-file *.txt --met-file *.txt
  
  index_base <- get_index_base_path(genome_id, "hisat2")
  input_paths <- get_reads_path(repl_id, type)
  
  out_paths <- get_alignmet_out_files_path(repl_id, genome_id, type, "hisat2")
  
  alignment_path <- out_paths$alignment
  summary_path <- out_paths$summary
  log_path <- out_paths$log
    
  max_seeds_param <- "--max-seeds 20"
  
  # Bowtie2: In end-to-end alignment mode, the default minimum score threshold is -0.6 + -0.6 * L, where L is the read length
  score_min_param <- "--score-min L,-0.6,-0.6"
  
  pred_quality_param <- "--phred33"
  cores_param <- paste("-p", cores, sep = " ")
  max_frag_len_valid_paired_param <- "-X 1000"
  additional_params <- "--no-spliced-alignment --no-mixed --no-discordant --no-unal"
  index_param <- paste("-x", index_base, sep = " ")
  r1_param <- paste("-1", input_paths$r1, sep = " ")
  r2_param <- paste("-2", input_paths$r2, sep = " ")
  alignment_path_param <- paste("-S", alignment_path, sep = " ")
  summary_param <- paste("--summary-file", summary_path, sep = " ")
  log_param <-  paste(">", log_path, sep = " ")
  
  
  command <- paste("hisat2", 
                   max_seeds_param,
                   score_min_param,
                   pred_quality_param,
                   cores_param,
                   max_frag_len_valid_paired_param,
                   additional_params,
                   index_param, 
                   r1_param,
                   r2_param,
                   alignment_path_param,
                   summary_param,
                   log_param,
                   "&",  #run in background
                   sep = " ")  
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)
  
}

align_bowtie2 <- function(repl_id, genome_id, type) {
  
  # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
  
  index_base <- get_index_base_path(genome_id, "bowtie2")
  input_paths <- get_reads_path(repl_id, type)
  
  out_paths <- get_alignmet_out_files_path(repl_id, genome_id, type, "bowtie2")
  
  alignment_path <- out_paths$alignment
  summary_path <- out_paths$summary
  log_path <- out_paths$log
  
  score_min_param <- "--score-min L,-0.6,-0.6"
  
  cores_param <- paste("-p", cores, sep = " ")
  
  max_frag_len_valid_paired_param <- "-X 1000"
   
  # Bowtie 2 with the --very-sensitive option is the same as running 
  # with options: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50.
  # individual preset values can be overridden by providing 
  # the specific options e.g. the configured seed length of 20 
  # in the [--very-sensitive] preset above can be changed to 25 
  # by also specifying the -L 25 parameter anywhere on the command line.
  additional_params <- "--end-to-end --very-sensitive --no-mixed --no-discordant --no-unal"
  
  index_param <- paste("-x", index_base, sep = " ")
  r1_param <- paste("-1", input_paths$r1, sep = " ")
  r2_param <- paste("-2", input_paths$r2, sep = " ")
  alignment_path_param <- paste("-S", alignment_path, sep = " ")
  log_param <-  paste(">", log_path, sep = " ")
  summary_param <- paste("2>", summary_path, sep = " ")
  
  
  command <- paste("bowtie2", 
                   cores_param,
                   score_min_param,
                   max_frag_len_valid_paired_param,
                   additional_params,
                   index_param, 
                   r1_param,
                   r2_param,
                   alignment_path_param,
                   log_param,
                   summary_param,
                   "&",  #run in background
                   sep = " ")  
  
  print(command)
  
  system(command, intern = FALSE,
        ignore.stdout = FALSE, ignore.stderr = FALSE,
        wait = TRUE, input = NULL, show.output.on.console = TRUE,
        timeout = 0)

}

sam_to_bam <- function(repl_id, genome_id, type, tool) {
  
  paths <- get_alignmet_out_files_path(repl_id, genome_id, type, tool)
  
  sam_file_path <- paths$alignment
    
  bam_file_path <- str_replace(sam_file_path, pattern = ".sam", replacement = ".bam")
  
  command <- paste("samtools view -Sb",
                   sam_file_path,
                   ">",
                   bam_file_path,
                   "&", #background
                   sep = " ")
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)

}


bam_to_bigwig <- function(repl_id, genome_id, type, tool) {
  
  script_file_path <- file.path(dirname(base_path), 'tools', "generate_BW_from_BAM.R")
  
  paths <- get_alignmet_out_files_path(repl_id, genome_id, type, tool)
  sam_file_path <- paths$alignment
  bam_file_path <- str_replace(sam_file_path, pattern = ".sam", replacement = ".bam")
  
  out_base_path <- dirname(bam_file_path)
  
  command <- paste(
    "Rscript", 
    script_file_path,
    "--file", bam_file_path,
    "--type", bigwigs$option_type,
    "--minLength", bigwigs$option_minLength,
    "--maxLength", bigwigs$option_maxLength,
    "--outputDir", out_base_path,
    #"&", #background
    sep = " ")
  
  print(command)
  
  system(command, intern = FALSE,
         ignore.stdout = FALSE, ignore.stderr = FALSE,
         wait = TRUE, input = NULL, show.output.on.console = TRUE,
         timeout = 0)
  
}