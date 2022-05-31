
# The tximport pipeline. 
# So this code has these following input that can be uptaken:
# (1) necessary: the father direcotry where all the salmon files folder locate + metadata that has sample_name column containing the name of each salmon folder.
# (2) necassary: specify the species. mouse/human will lead to the usage of a online ensembl data. users can also input "none" and use "--gtf" to input their own reference columns. But it is super slow, any suggestion...?
# (3) optional: counts from abundance. There are 4 options. read ?tximport for more information about difference between these 4.
# (4) optional: transcript or not. default is FALSE and will output the gene-level summarization from salmon files. but if users want the output table to be at transcript level (just a output of the original data), this command can achieve.
# (5) optional: save directory. the default is user's working directory.
# (6) optional: this function can also use the imported data to perform a default deseq2 test. Default is FALSE.
# (7) optional: for deseq2, users can use --setcontrol and --setcomparegroup to specify which group of samples are comparing and which one is the control that other groups are comparing with.


#make option that intake the directory etc
suppressPackageStartupMessages(require(optparse))


option_list = list(
 make_option(c("-s","--salmondirectory"),action="store",type="character",default=NA,help="the directory where all the salmon quant files locate"),
 make_option(c("-m","--metadata"),action="store",type="character",default=NA,help="the directory where the metadata file locates"),
 make_option(c("-t","--tx2gene"),action="store",type="character",default=NA,help="Tx2gene file for the GTF the data was aligned to"),
 make_option(c("-o","--outputdir"),action="store",default=getwd(), help="set the directory you want to save the results"),
 make_option(c("-g","--column_name"),action="store", default="NA",help="Column name containing the metadata of baseline and contrasts that will be used"),
 make_option(c("-b","--baseline"),action="store",default="baseline", help="What the control is called in the column name"),
 make_option(c("-c","--contrast"),action="store",default="baseline", help="What the contrast is called in the column name"),
 make_option(c("-n","--controls_name"),action="store",default="baseline", help="What the control should be called in the output files"),
 make_option(c("-t","--contrast_name"),action="store",default="baseline", help="What the contrast should be called in the output files"),

)
opt = parse_args(OptionParser(option_list=option_list))



#directory
salmon_quant_directory = opt$directory
outputdir = opt$outputdir
metadata_dir = opt$metadata
tx2gene = opt$tx2gene
column_name = opt$column_name
baseline = opt$baseline
contrast = opt$contrast
controls_name = opt$controls_name
contrast_name = opt$contrast_name
column_name = opt$column_name

output_path=paste0(controls_name,controls_name,"-",contrast_name)

# ============================ section 1: import data ===========================
library(dplyr)
library(tidyr)
library(tximport)
library(rlang)
library(DESeq2)


tx2gene <- data.table::fread(tx2gene,header=FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

#(1) First read in the metadata. if only a subset of the files are used, the opt$pattern option will be taken.
metadata = read.csv(metadata_dir,header=TRUE)

metadata = metadata %>% 
    select(sample_name, !!as.symbol(column_name)) %>% 
    mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline',
                                            !!as.symbol(column_name) == contrast ~ 'contrast',
                                            TRUE ~ NA_character_)) %>% 
    filter(!is.na(comparison_condition))

#(2) Generate a vector of the wanted file names.
files = unique(file.path(salmon_quant_directory,metadata$sample_name,"quant.sf")) 
names(files) = unique(metadata$sample_name)


#(3) To check if all the files exist
if(all(file.exists(files)) == FALSE) {
  stop("It seems that I cannot find those files...Please check if your directory is correct.")
}


# ====================== section 3: import salmon files ==============================
# files is a vector of directory where quant.sf file locates.
# just ignore the version... to make it easier for following steps.
txi.tx <- tximport(files, 
               type="salmon", 
               tx2gene=tx2gene,
               ignoreTxVersion = TRUE,
               ignoreAfterBar = TRUE,
               txOut = TRUE)

txi.sum <- summarizeToGene(txi.tx, tx2gene)

# make it csv
TPM_transcripts = as.data.frame(txi.tx$abundance) %>% 
  tibble::rownames_to_column(.,var="transcript_id")
TPM_gene = as.data.frame(txi.sum$abundance) %>% 
  tibble::rownames_to_column(.,var="gene_id")

write.csv(TPM_transcripts,file=paste0(output_path,"TPM_transcripts.csv"))
write.csv(TPM_gene,file=paste0(output_path,"TPM_gene.csv"))


# ========================================== section 4: RUN A DEFAULT DESEQ 2 (optional) =============================================================


dds = DESeqDataSetFromTximport(txi.sum,
                               colData = metadata,
                               design = ~ comparison_condition) 


# 'Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts' -- <Deseq2 vignette> (just a note - Deseq() wraps the normalization step inside)
# perform the Deseq function
dds = DESeq(dds)

# Now, extract the result and named them by their contrast group
results_table = results(dds) %>%
            as.data.frame() %>%
            tibble::rownames_to_column('Geneid')
# Now, extract the DESeq2 normed counts
normed_counts = counts(dds, normalized = TRUE) %>%
            as.data.frame() %>%
            tibble::rownames_to_column('Geneid')


write.csv(results_table,file=paste0(output_path,"DESeq2_results.csv"))
write.csv(normed_counts,file=paste0(output_path,"DESeq2_normalized_counts.csv"))