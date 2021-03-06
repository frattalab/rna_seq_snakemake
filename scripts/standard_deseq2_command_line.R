
library("optparse")
option_list = list(
    make_option(c("-f", "--folder_of_featurecounts"), type="character", default=NULL,
                help="a folder containing featurecounts output", metavar="character"),
    make_option(c("-b", "--base_grep"), type="character", default=NULL,
        help="a folder containing featurecounts output", metavar="character"),
    make_option(c("-c", "--contrast_grep"), type="character", default=NULL,
                help="output file name", metavar="character")
    make_option(c("-o", "--out"), type="character", default="out.txt",
          help="output file name; will output 3 files", metavar="character")
    make_option(c("-o", "--out"), type="character", default="out.txt",
          help="output file name; will output 3 files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

standard_output = run_standard_deseq(opt$folder_of_featurecounts,
                              base_grep = opt$base_grep,
                              contrast_grep = opt$contrast_grep,
                              grep_pattern = opt$grep_pattern,
                              suffix = opt$suffix,
                              baseName = opt$baseName,
                              contrastName = opt$contrastName
                              )

normed_counts_standard =  as.data.frame(counts(standard_output$deseq_obj, normalized = T)) %>% rownames_to_column("sample_name")
unnormed_counts_standard =  as.data.frame(counts(standard_output$deseq_obj)) %>% rownames_to_column("sample_name")

standard_output$results_table %>% fwrite(paste0(opt$out, "results.csv"))

meta_data = as.data.frame(colData(standard_output$deseq_obj)) %>% rownames_to_column("sample_name")
meta_data %>% fwrite(paste0(opt$out, "meta_data.csv"))
unnormed_counts_standard %>% fwrite(paste0(opt$out, "unnormed_counts.csv.gz"))
normed_counts_standard %>% fwrite(paste0(opt$out, "normed_counts.csv.gz"))
