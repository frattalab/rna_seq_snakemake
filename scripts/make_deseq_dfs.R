#this function takes the total rna tables produced by featureCounts, and gives a reasonable output data frame and
#metadata frame
make_deseq_dfs = function(total_table, grep_pattern = "", leave_out = "", base_grep = "", contrast_grep = ""){

  if(grep_pattern == ""){

    grep_pattern = glue::glue("{base_grep}|{contrast_grep}")

  }
  #grep pattern is being used to select small parts of this overall
  total_table = as.data.table(total_table, keep.rownames = TRUE)

  if('gene_name' %in% colnames(total_table)){
    total_table$gene_name = NULL

  }
  if(leave_out == ""){
    conv_df = as.data.frame(total_table[,round(.SD),
                                        .SDcols = grep(grep_pattern,names(total_table))])
  }else{
    conv_df = as.data.frame(total_table[,round(.SD),
                                        .SDcols = grep(grep_pattern,names(total_table))])
    drop_col = paste0(leave_out_sample,".aligned.sorted.out.bam")
    print(paste("Dropping column:",drop_col ))
    conv_df = conv_df[ , !(names(conv_df) %in% drop_col)]
  }

  if("geneid" %in% names(total_table)){
    rownames(conv_df) = total_table$geneid
  }else if("Geneid" %in% names(total_table)){
    rownames(conv_df) = total_table$Geneid
  }else if("rn" %in% names(total_table)){
    rownames(conv_df) = total_table$rn
  }
  else{
    rownames(conv_df) = total_table$gene
  }
  ###damn I really need to delte this...it's vestigal
  coldata = as.data.table(names(conv_df))
  if(base_grep == "" & contrast_grep == ""){
  }else if(base_grep != ""){
    coldata[grep(base_grep,V1), cond := "base"]
  }else if(contrast_grep != ""){
    coldata[grep(contrast_grep,V1), cond := "contrast"]
  }
  coldata[is.na(cond), cond := "contrast"]
  print("This is your metaData")
  print(coldata)
  coldata = as.data.frame(coldata[,2:ncol(coldata)])
  rownames(coldata) = names(conv_df)
  morphed = list(conv_df,coldata)
  names(morphed) = c("conv_df","coldata")


  return(morphed)
}


# this function takes a deseq metadata df, the name you want to call the baseline and contrast conditions and returns the metadata table with
# a new column called 'comparison' which is a factor with baseline as the first level
rename_relevel_for_deseq = function(coldata, baseName = "", contrastName = ""){
  coldata$comparison  = factor(ifelse(coldata$cond == "base", baseName, contrastName), levels = c(baseName, contrastName))
  return(coldata)
}

# this function filters a count table
filter_count_table = function(count_table){
  keep <- rowSums(edgeR::cpm(count_table) > 0.5) >= 2
  print("Filtered Genes by CPM greater than 0.5 in a least 2 samples")
  print(table(keep))
  return(count_table[keep, ])
}

# this function takes the big feature counts table and determines which samples are female
find_females = function(featureCountsTables, species = "human"){
  if(species == "human"){
    gene_table = fread("/Users/annaleigh/Documents/data/reference_genomes/gencode.v31.parsed_names.txt")
    xist = "XIST"
  }else if(species == "mouse"){
    gene_table = fread("/Users/annaleigh/Documents/data/reference_genomes/gencode.vM22.parsed_names.txt")
    xist = "Xist"
  }else{
    print("Species either 'human' or 'mouse")
    return(1)
  }
  # check that the right species was chosen
  if(sum(featureCountsTables$Geneid %in% gene_table$gene_id) < 10 ){
    print("are you sure you chose the right species? Nothing found ")
    return(1)
  }
  # --- find all the female samples
  female_xist = featureCountsTables %>%
    left_join(gene_table %>% dplyr::select(gene_id,gene_name) %>% unique(),by = c("Geneid" = "gene_id")) %>%
    filter(gene_name == xist) %>%
    select_if(is.numeric) %>% select_if( . > 100) %>% colnames()
  return(female_xist)
}

# this function adds a metadata column for gender
add_gender = function(coldata, femalelist){
  coldata = coldata %>% rownames_to_column("samp") %>%
    mutate(sex = ifelse(samp %in% femalelist, "F", "M")) %>%
    column_to_rownames('samp')
  return(coldata)
}
