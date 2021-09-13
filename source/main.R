base_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(file.path(base_path, 'common', 'config.R'))
source(file.path(base_path, 'common', 'data.R'))
source(file.path(base_path, 'pipeline', 'qc.R'))
source(file.path(base_path, 'pipeline', 'alignment.R'))
source(file.path(base_path, 'pipeline', 'heatmap.R'))
source(file.path(base_path, 'pipeline', 'summary_parser.R'))
  

build_indexes <- function(genome_ids) {
 
  for(genome_id in genome_ids) {
    
    build_bowtie2_index(genome_id)
    build_hisat2_index(genome_id)
    
  }

}

preprocess <- function(repl_ids) {

  for(repl_id in repl_ids) {
    
    check_ok <- check_md5sum(repl_id)
    
    if(check_ok) {
      quality_control(repl_id, 'raw')
      clean_repetitive_sequences(repl_id)
      quality_control(repl_id, 'clean')
    } else {
      print("Invalid md5sum value")
    }
  }
  
}

align <- function(repl_ids, types, genome_ids, tools) {
  for (repl_id in repl_ids) {
    for (type in types) {
      for (genome_id in genome_ids) {
        for (tool in tools) {
          
          print(paste("processing:", repl_id, type, genome_id, tool))
          
          if(tool == 'hisat2') {
            align_hisat2(repl_id, genome_id, type)
          } else if(tool == 'bowtie2')  {
            align_bowtie2(repl_id, genome_id, type)  
          }
        }
      }
    }
  }
}

convert_to_bam <- function(repl_ids, types, genome_ids, tools) {

  for (repl_id in repl_ids) {
    for (type in types) {
      for (genome_id in genome_ids) {
        for (tool in tools) {
          
          print(paste("processing:", repl_id, type, genome_id, tool))
          
          sam_to_bam(repl_id, genome_id, type, tool)
          
        }
      }
    }
  }
}

convert_to_bigwig <- function(repl_ids, types, genome_ids, tools) {
  
  for (repl_id in repl_ids) {
    for (type in types) {
      for (genome_id in genome_ids) {
        for (tool in tools) {
          
          print(paste("processing:", repl_id, type, genome_id, tool))
          
          bam_to_bigwig(repl_id, genome_id, type, tool)
          
        }
      }
    }
  }  
  
}
  
plot <- function(repl_ids, types, genome_ids, tools) {
  for (repl_id in repl_ids) {
    for (type in types) {
      for (genome_id in genome_ids) {
        for (tool in tools) {
          
          print(paste("processing:", repl_id, type, genome_id, tool))
          
          plot_heatmap(repl_id, genome_id, type, tool)
        }
      }
    }
  }
}

summary_tables <- function(repl_ids, types, genome_ids, tools) {

  for (repl_id in repl_ids) {
    
    print(paste("processing:", repl_id, types, genome_ids, tools))
    
    build_paper_table(c(repl_id), types, genome_ids, tools)
    
  }  
  
}

main <- function(step_name) {
  
  print(step_name)

  genome_ids <- names(genomes)
  
  repl_ids <- names(replicates)

  types <- c('raw', 'clean')
  
  tools <- c('hisat2', 'bowtie2')
  
  arguments_values <- switch(  
    step_name,  
    "build_indexes" = list(genome_ids),  
    "preprocess" = list(repl_ids),  
    "align" = list(repl_ids, types, genome_ids, tools),  
    "convert_to_bam" = list(repl_ids, types, genome_ids, tools),
    "plot" = list(repl_ids, types, genome_ids, tools),
    "summary_tables" = list(repl_ids, types, genome_ids, tools),
    "convert_to_bigwig"= list(repl_ids, types, genome_ids, tools)
  )    
  
  do.call(step_name, arguments_values)

  
}

# available steps: 
#   build_indexes, 
#   preprocess, 
#   align, 
#   convert_to_bam,
#   plot, 
#   summary_tables, 
#   convert_to_bigwig

#main('preprocess')
main('plot')
#main('summary_tables')
#main('convert_to_bigwig')


