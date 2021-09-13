library(stringr)


summary_map = list('reads' = list(
                      'line_number' = 1,
                      'expression' = '[:digit:]+'),
                   'paired' = list(
                     'line_number' = 2,
                     'expression' = '[:digit:]+'),
                   'paired_perc' = list(
                     'line_number' = 2,
                     'expression' = '([:digit:]+\\.[:digit:]+\\%)'),
                   'aligned_concordantly_0' = list(
                     'line_number' = 3,
                     'expression' = '[:digit:]+'),
                   'aligned_concordantly_0_perc' = list(
                     'line_number' = 3,
                     'expression' = '([:digit:]+\\.[:digit:]+\\%)'),
                   'aligned_concordantly_1' = list(
                     'line_number' = 4,
                     'expression' = '[:digit:]+'),
                   'aligned_concordantly_1_perc' = list(
                     'line_number' = 4,
                     'expression' = '([:digit:]+\\.[:digit:]+\\%)'),
                   'aligned_concordantly_gt_1' = list(
                     'line_number' = 5,
                     'expression' = '[:digit:]+'),
                   'aligned_concordantly_gt_1_perc' = list(
                     'line_number' = 5,
                     'expression' = '([:digit:]+\\.[:digit:]+\\%)'),
                   'overall_alignment_rate' = list(
                     'line_number' = 6,
                     'expression' = '([:digit:]+\\.[:digit:]+\\%)')
                   )

parse <- function(repl_id, genome_id, type, tool) {
  
  out_paths <- get_alignmet_out_files_path(repl_id, genome_id, type, tool)
  summary_path <- out_paths$summary
  
  result <- c()
  result_names <- c()
  
  if(file.exists(summary_path)) {
    
    data <- readLines(summary_path)
    
    for(key in names(summary_map)) {
    
      line_number <- summary_map[[key]]$line_number
      pattern <- summary_map[[key]]$expression
      value <- data[line_number]
      parsed_value <- str_extract(value, pattern)
      
      #result[[key]] <- parsed_value
      result <- c(result, parsed_value)
      result_names <- c(result_names, key)
      
    }
  }
  
  names(result) <- result_names
  
  return(result)
  
}

build_paper_table <- function(repl_ids, types, genome_ids, tools) {
  
  result <- data.table(
    'tool' = character(),
    'type' = character(),
    'genome' = character(),
    stringsAsFactors = FALSE
  )
  
  result[, names(summary_map)] <- rep('', length(summary_map))
  
  for (repl_id in repl_ids) {
    for (type in types) {
      for (genome_id in genome_ids) {
        for (tool in tools) {
          
          print(paste("processing:", repl_id, type, genome_id, tool))
          
          row_values_parse <- parse(repl_id, genome_id, type, tool)
          
          row_values <- c(tool, type, genome_id, row_values_parse)
          
          names(row_values) <- colnames(result)
          
          result <- rbind(result, t(row_values))

        }
      }
    }
  }  
  

  out_file_name <- paste(
    paste(repl_ids, collapse = "-"),
    paste(types, collapse = "-"),
    paste(genome_ids, collapse = "-"),
    paste(tools, collapse = "-"),
    sep = "_")
  out_file_name <- paste0(out_file_name, ".cvs")
  
  out_file_path <- file.path(alignment_summary_out_base_path, out_file_name)

  write.table(x = result, file = out_file_path, sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)

}


