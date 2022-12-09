library(dplyr)
library(tidyr)
library(tximport)
library(rlang)
library(DESeq2)
library(annotables)
library(tidyverse)
library(optparse)
library(yaml)
library(data.table)


  


run_salmon_deseq = function(salmon_quant_directory, 
                            metadata_filepath, 
                            tx2gene, 
                            outputdir)
{
  #1. import data
  #Read in the comparison design from .yaml
  deseq2Config = read_yaml('config/DESeq2comparisons.yaml') 
  comparison = purrr::map(deseq2Config, as.data.table)
  
  #read in the tx2gene file
  tx2gene <- data.table::fread(tx2gene,header=FALSE)
  colnames(tx2gene) = c("TXNAME", "GENEID")
  
  #Read in the metadata
  metadata_orig = read.csv(metadata_filepath,header=TRUE)
  
  
  for (x in comparison){
  
  #Which columns to read in the metadata
  column_name = unique(x[[1]])
  baseline = x[[2]]
  contrast = x[[3]]
  
  #What folder and what  will my outputs be written to?
  baseline_name = colnames(x)[2]
  contrast_name = colnames(x)[3]
  exp = paste0(baseline_name, "_", contrast_name)
  output_path = paste0(outputdir,"/", exp)
  
  #Make a new minimal metadata file containing the sample for this comparison
  #heads up! there's some kinda advance R code below turning a variable (column_name) into 
  #quotation in order to be able to select from a data.frame using what's stored in the variable
  
  metadata = metadata_orig %>%  
    filter(is.na(exclude_sample_downstream_analysis)) %>% 
    dplyr::select(unit, sample_name, !!(column_name)) %>% #!!(column_name) - this is called 'unquoting
    mutate(comparison_condition = case_when(!!as.symbol(column_name) %in% baseline ~ 'baseline', #here we need as symbol due to the logic of case_when
                                            !!as.symbol(column_name) %in% contrast ~ 'contrast',
                                            TRUE ~ NA_character_)) %>% 
    filter(!is.na(comparison_condition))   
  
  #Generate a vector of the wanted file names.
  
  #files = unique(file.path(salmon_quant_directory,paste0(metadata$unit,"_",metadata$sample_name),"quant.sf")) 
  files = unique(file.path(salmon_quant_directory,metadata$sample_name,"quant.sf")) 
  
  names(files) = unique(metadata$sample_name)
  
  #To check if all the files exist
  if(all(file.exists(files)) == FALSE) {
    stop("It seems that I cannot find those files...Please check if your directory is correct.")
  }
  
  
  #import transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages
  #output matrix: average transcript length, weighted by sample-specific transcript abundance estimates 
  txi.tx <- tximport(files, 
                     type="salmon", 
                     tx2gene=tx2gene,
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
                     txOut = TRUE)
  
  #For downstream analysis it can sometimes
  #be helpful to have a table of all the transcript counts
  transcript_counts = as.data.frame(txi.tx$abundance) %>% 
    tibble::rownames_to_column('TXNAME') %>% 
    left_join(tx2gene)
  
  
  fwrite(transcript_counts,paste0(output_path,".transcript_counts.csv"))
  
  txi.sum <- summarizeToGene(txi.tx, tx2gene)
  
  #Run a default DeSeq2
  dds = DESeqDataSetFromTximport(txi.sum,
                                 colData = metadata,
                                 design = ~ comparison_condition) 
  dds = DESeq(dds)
  
  results_table = results(dds) %>%  
    as.data.frame() %>% 
    tibble::rownames_to_column('ensgene') %>%   
    #remove fractions in ensgene
    mutate(ensgene = gsub("\\..*","",ensgene))  %>% 
    left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol)) %>% 
    dplyr::rename(gene_name = symbol) %>% 
    unique() %>% as.data.table()
  
  # Now, write everything out
  
  output_path=paste0(outputdir,"/",exp)
  
  fwrite(results_table,paste0(output_path,".DESEQ2_results.csv"))
  
  #Also save the dds object to later capture
  #more metadata
  saveRDS(dds, file.path(paste0(output_path,".DESEQ2_object.RDS")))
  
  }
  
  return()
}



option_list = list(
    make_option(c("-q", "--salmon_quant"), type = "character", default = NULL, help = "file path for the folder containing salmon_quant output", metavar = "character"),
    make_option(c("-s", "--samplesheet"), type = "character", default = NULL, help = "file path for the sample sheet", metavar = "character"),
    make_option(c("-g", "--tx2gene"), type = "character", default = NULL, help = "file path for reference tx2gene file", metavar = "character"),
    make_option(c("-o", "--outputdir"), type = "character", default = NULL, help = "directory for DESeq2 outputs", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

deseq2_output = run_salmon_deseq(opt$salmon_quant, 
                                 opt$samplesheet, 
                                 opt$tx2gene, 
                                 opt$outputdir)


