
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
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_standard_deseq = function(opt$folder_of_featurecounts,
                              opt$base_grep = "Ctrl",
                              opt$contrast_grep = "TDPKD",
                              opt$grep_pattern = "",
                              opt$suffix = ".Aligned.sorted.out.bam",
                              opt$baseName = "control",
                              opt$contrastName = 'TDPKD'
                              ){

output = janitor::clean_names(read_STAR(path = opt$STAR,reshape = TRUE))

data.table::fwrite(output, opt$out)

buratti_dzap = run_standard_deseq(folder_of_featurecounts = "/Users/annaleigh/Documents/GitHub/sinai_splice/data/buratti_dzap/",
                             base_grep = "LUC",
                             contrast_grep = 'TDP',
                             suffix = "unique_rg_fixed.bam")
