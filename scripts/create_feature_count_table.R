
create_feature_count_table = function(feature_count_folder, suffix = ".Aligned.sorted.out.bam"){
  library(data.table)

  # feature counts gives you feature coutns and summariues, make sure that all the files end in this pattern
  # _featureCounts_results.txt
  feature_count_files = grep(list.files(path=feature_count_folder,full.names = TRUE), pattern='_featureCounts_results.txt', value=T)
  feature_count_files = grep(feature_count_files,pattern = "summary", invert = T, value = T)

#throw an erros if the length is wrong
  if(length(feature_count_files) == 0){
    print("Hey there, you don't have any feature counts files, are you sure you did the right folder?")
    return(1)
  }
  # ensure that no directories are brought in
  feature_count_files = feature_count_files[!file.info(feature_count_files)$isdir]

  # from the list, read them all
  l <- lapply(feature_count_files, fread)
  # get the file path that every folder ends with
  joint_file_path = l[[1]] %>% colnames() %>% tail(.,1) 
  #remove everything starting with this
  prefix = gsub(basename(joint_file_path),"",joint_file_path) %>% make.names()

  # make it into a data table
  wide_feature_counts = setDT(unlist(l, recursive = FALSE), check.names = TRUE)[]
  # there's only a few columns we really want from this table it's going to be
  # Geneid and the things that labeled by their file name so we do two things
  # first when we do this setDT unlist thing names become X.path.to. so we'lre gooing to take the feature counts files and replace the / with .
  selector_list = c("Geneid",names(wide_feature_counts)[grepl("^X",names(wide_feature_counts))])
  if('gene_name' %in% colnames(wide_feature_counts)){
    selector_list = c(selector_list, 'gene_name')
  }

# take only those columns
  wide_feature_counts = wide_feature_counts[,.SD, .SDcols = selector_list]
  # ----remove unnecessary information from the column names
  if(suffix != ""){
    colnames(wide_feature_counts) = gsub(suffix,"",colnames(wide_feature_counts))
  }
  if(prefix != ""){
    colnames(wide_feature_counts) = gsub(prefix,"",colnames(wide_feature_counts))
  }
  return(wide_feature_counts)
}

calc_rpkm = function(feature_counts_table,species = "human"){
  
  gene_lengths = fread("/Users/annaleigh/Documents/data/reference_genomes/gencode.v31_gene_lengths.txt")
  y <-  edgeR::DGEList(counts=feature_counts_table,genes=gene_lengths)
  y <-  edgeR::calcNormFactors(y)
  RPKM <- edgeR::rpkm(y)
  print(y$samples)
  return(RPKM)
}