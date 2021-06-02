
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
 make_option(c("-d","--directory"),action="store",type="character",default=NA,help="the directory where all the salmon quant files locate"),
 make_option(c("-m","--metadata"),action="store",type="character",default=NA,help="the directory where the metadata file locates"),
 make_option(c("-a","--species"),action="store",type="character",default="mouse",help="what species is your data from? option: mouse/human/none. default is mouse.If you want to use your own gtf, please choose none and specify your gtf file in the --g option."),
 make_option(c("-g","--gtf"),action="store",type="character",default=NA,help="if you need to input your own GTF file to make transcript-gene annotation list, this is where you put the GTF file directory. Also, please use --species and input 'none'. Please note that this method is very VERY time-consuming."),
 make_option(c("-c","--countsfromabundance"),action="store",default="no",help="import the counts from abundance -- inputs you can make: no/scaledTPM/lengthScaledTPM/dtuScaledTPM"),
 make_option(c("-s","--savedirectory"),action="store",default=getwd(), help="set the directory you want to save your data"),
 make_option(c("-t","--deseq2ornot"), action="store_true",default=FALSE,help="please put this option if you want the deseq2 is also done"),
 make_option(c("-b","--setcontrol"),action="store",default="control", help="specify the name of the control"),
 make_option(c("-e","--setcomparegroup"),action="store",default="group",help="specify the category you want to compare. Should be a category in metadata."),
 make_option(c("-p","--pattern"),action="store", default="NA",help="if only need to include a subset of the data,use regular expression to specify the pattern, as select a subset of rows in the metadata. Because this function use information inside metadata to import data"),
 make_option(c("-x","--transcriptornot"),action="store_true",default=FALSE,help ="if you want your output table to be associated with transcript level. if yes, put this in the command line.")
 )
opt = parse_args(OptionParser(option_list=option_list))



#directory

salmon_quant_directory = opt$directory
metadata_dir = opt$metadata
gtf_dir = opt$gtf

# ============================ section 1: import data ===========================
library(tximportData)
library(tidyverse)
library(tidyr)
library(dplyr)
library(tximport)

pattern = opt$pattern

#(1) First read in the metadata. if only a subset of the files are used, the opt$pattern option will be taken.
metadata = read.csv(metadata_dir,header=TRUE)
if(!is.na(opt$pattern)){
  metadata = metadata %>% filter(test == "A" | test =="B")
}
#(2) Generate a vector of the wanted file names.
files = unique(file.path(salmon_quant_directory,metadata$sample_name,"quant.sf")) 
names(files) = unique(metadata$sample_name)


#(3) To check if all the files exist
if(all(file.exists(files)) == FALSE) {
  stop("It seems that I cannot find those files...Please check if your directory is correct.")
}

# ============================ section 2: provide the transcript | gene name table =======================================
# we need to associate the geneIDs with transcripts, because salmon output only provides the transcript ID.
# I want the next chunk can reach this effect: you can choose to use the online ensembl database, or you can input gtf to make one. But reading gtf is significantly slower than that graping online data.
# I guess this is because for the online source, biomart allow you to just extract a subset of information, but reading gtf is reading a extremely large GRange file.
# I tried MakeTxDbfromgtf(), I then give up because (1) it is very slow (2) txdb do not store information like ensembl ID, gene name etc...

species = opt$species

if(species =="mouse"){library(annotables)
  an = annotables::grcm38_tx2gene %>% 
    left_join(annotables::grch38 %>% 
                dplyr::select(c("ensgene","symbol")))
  tx2gene = an_ms %>% select(c("enstxp","ensgene"))
}

if(species == "human"){library(annotables)
  an = annotables::grch38_tx2gene %>% 
    left_join(annotables::grch38 %>% 
                dplyr::select(c("ensgene","symbol")))
  tx2gene = an_hm %>% select(c("enstxp","ensgene"))
}

if(species == "none"){library(rtracklayer)
  gtf = import.gff2(gtf_dir)
  an = as.data.frame(gtf@elementMetadata@listData)[,c("gene_id","transcript_id","gene_name")] %>%
    separate("transcript_id",c("transcript_id","transcript_version")) %>% 
    separate("gene_id",c("gene_id","gene_version")) %>%
    dplyr::select(c("transcript_id","gene_id","gene_name")) %>%
    unique(.)
  tx2gene = select(an, c("transcript_id","gene_id","gene_name"))
}

if(species != "mouse"& species!= "human"& species != "none"){
  stop("please specify your species, or provide a reference gtf")
}

# ====================== section 3: import salmon files ==============================
# files is a vector of directory where quant.sf file locates.
# just ignore the version... to make it easier for following steps.
txi = tximport(files, 
               type="salmon", 
               tx2gene=tx2gene,
               ignoreTxVersion = TRUE,
               ignoreAfterBar = TRUE,
               countsFromAbundance = opt$countsfromabundance,
               txOut = opt$transcriptornot) 

# make it csv
TPM = as.data.frame(txi$abundance) %>% 
  tibble::rownames_to_column(.,var="ensembl_gene_id")
# The annotation table contains 3 columns, gene id/transcript id/gene name. because there is possibility that user will use their gtf, so I think it is better to make the columns name with a fixed name.
an = an %>% dplyr::rename(ensembl_transcript_id= colnames(an)[1]) %>%
  rename(ensembl_gene_id = colnames(an)[2]) %>%
  rename(gene_name = colnames(an)[3])
symbol = unique(an[,c("gene_name","ensembl_gene_id")])
TPM = left_join(TPM,symbol,by="ensembl_gene_id") %>% tibble::column_to_rownames(.,var = "ensembl_gene_id")

write.csv(TPM,file=paste0(opt$directory,"TPM.csv"))


# ========================================== section 4: RUN A DEFAULT DESEQ 2 (optional) =============================================================

if (opt$deseq2ornot == TRUE) {
library(DESeq2)
deseq2_metadata = unique(metadata[,!(names(metadata) %in% c("unit","fast1","fast2","exclude_sample_downstream_analysis"))])
if(!is.na(opt$pattern)){
  deseq2_metadata = deseq2_metadata %>% filter(opt$pattern)
}

design = opt$setcomparegroup
dds = DESeqDataSetFromTximport(txi,
                               colData = deseq2_metadata,
                               design = formula(paste("~", design, collapse=" "))) 

# set the control if there are multiple groups.

dds[[design]] = relevel(dds[[design]],ref= opt$setcontrol)

# 'Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts' -- <Deseq2 vignette> (just a note - Deseq() wraps the normalization step inside)
# perform the Deseq function
dds = DESeq(dds)

# Now, extract the result and named them by their contrast group
for (disease_vs_control in resultsNames(dds)) {
  a= as.data.frame(DESeq2::results(dds,name= disease_vs_control)) %>%
    tibble::rownames_to_column("ensembl_gene_id")%>%
    full_join(symbol,by="ensembl_gene_id")
  assign(paste0(disease_vs_control,"_result"),a)%>%
    write.csv(file = paste0(opt$savedirectory,disease_vs_control,".csv"))
}
}
