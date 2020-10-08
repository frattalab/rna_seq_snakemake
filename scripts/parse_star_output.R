#!/usr/bin/env Rscript
#gratefully stolen from github
library(tidyverse)

#' Parse sample from file names
#'
#' Sample names are parsed from file names without the extension.  If the file name is not unique,
#'    the parent directory is used.
#'
#' @param files file name with path
#'
#' @return A vector
#'
#' @author Chris Stubben
#'
#' @examples
#'  # File name or parent directory should be unique
#'  extract_samples(c("align1/1355X1.counts", "align2/1355X2.counts"))
#'  extract_samples(c("align1/1355X1/Log.out", "align2/1355X2/Log.out"))
#' @export


extract_samples <- function(files){
  ## capture sample in file name  or path
  x <- lapply( strsplit(files, "/"), rev)
  samples <- sapply(x, "[", 1 )
  # remove file extension
  samples <- gsub("\\..*", "", samples)
  if( any(duplicated(samples )) ){
      # use parent directory, 13555X2/Log.final.out
      samples <- sapply(x, "[", 2 )
  }
  if( any(duplicated(samples )) ){
     stop("Sample names are not unique:  \n  ", paste(files, collapse="\n  "), call.=FALSE)
  }
  samples
}
#' Read and combine common output files
#'
#' Read output files and parse the file name to add sample IDs in the first column
#'
#' @param path the path to output files, the default corresponds to the working directory.
#' @param pattern regular expression for file name matching
#' @param delim separator in output files, default table
#' @param \dots additional options such as col_names passed to \code{read_delim}.
#'
#' @note Sample names are parsed from file names without extensions.  If the file name is not unique,
#'    the parent directory is used.   Requires tibble > version 1.2 to avoid error in add_column
#'
#' @return A list with coverage and stats data.frames
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    #FeatureCounts summary (second column name with *.bam is always unique, so skip and assign)
#'    fc <- read_sample_files(".summary$", skip=1, col_names=c("status", "count"))
#'  filter(fc, count!=0) %>%
#'    hchart("bar", x=sample, y=count, group=status) %>%
#'     hc_plotOptions(bar = list(stacking = "normal"))
#' }
#' @export

read_sample_files <- function(path=".", pattern="\\.counts$", delim="\t",  ...){
   outfiles <- list.files(path, pattern, recursive=TRUE, full.names=TRUE)

   outfiles = outfiles[!grepl("pass1",outfiles)]
   if(length(outfiles) == 0) stop("No ", pattern, " files found in ", path, call.=FALSE)
   samples <- extract_samples(outfiles)

   out1 <- vector("list", length(outfiles))
   for(i in seq_along(outfiles)){
       message("Reading ", outfiles[i])
       x <- suppressMessages( readr::read_delim(outfiles[i], delim=delim, ...) )
       ## requires tibble > 1.2 (to add single name into column)
       out1[[i]] <- tibble::add_column (x, sample= samples[i], .before=1)
   }
   dplyr::bind_rows(out1)
}

#' Read STAR log files
#'
#' Read STAR Log.final.out files and optionally reshape into wide format.
#'
#' @param path the path to STAR log files, the default corresponds to the working directory.
#' @param pattern regular expression for file name matching, default .final.out
#' @param reshape reshape percent mapping into wide format with samples in rows
#'
#' @note Reading output files requires a unique sample identifier in either the file name or parent directory
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- read_STAR( pattern=".star.out$")
#'   reads <- c("Uniquely mapped",  "Mapped to multiple loci",
#'                "Mapped to too many loci", "Unmapped reads")
#'   y <- filter(x, stat %in% reads) %>% mutate( stat= factor(stat, levels=reads ))
#'          hchart(y, "bar", x=sample ,  y=value, group=stat ) %>%
#'             hc_plotOptions(bar = list(stacking = "normal")) %>%
#'              hc_colors( c('#437bb1', '#7cb5ec', '#f7a35c', '#b1084c') ) %>%
#'                hc_yAxis(reversedStacks = FALSE)
#'  read_STAR( pattern=".star.out$", reshape=TRUE)
#' }
#' @export

read_STAR <- function( path=".", pattern, reshape=FALSE){
   if(missing(pattern))  pattern <- "\\.final.out$"
   # suppress warnings about subheaders....  Warning: 4 parsing failures.
   x <- suppressWarnings( read_sample_files(path, pattern, col_names=c("stat", "value") , skip=5, trim_ws=TRUE ))

   x <- dplyr::filter(x, !is.na(value)) %>%                    # drop subheader with NA value ... UNIQUE READS:
  dplyr::mutate( stat = gsub(" \\|$", "", stat),              # remove pipe from end of statistic
                 stat = gsub("Number of reads m", "M", stat), # shorten  Number of reads mapped to...
                 stat = gsub(" reads number", "", stat),      #  shorten Uniquely mapped reads number
                value = gsub("%", "", value),                 # drop % from value
                value = as.numeric(value))                    # change value to numeric

   # ADD unmapped reads  -
   mapped <- c("Uniquely mapped",  "Mapped to multiple loci", "Mapped to too many loci" )
   # see http://stackoverflow.com/questions/40749742/add-missing-subtotals-to-each-group-using-dplyr
   x <- dplyr::group_by(x, sample) %>%
          dplyr::summarize(  value = value[stat == 'Number of input reads'] - sum(value[stat %in% mapped ]),
                     stat = 'Unmapped reads' ) %>%
          dplyr::bind_rows(x) %>%
           dplyr::select(sample, stat, value) %>%  # back to original column order
            dplyr::arrange(sample)

   ##  split unmapped reads into too many mismatches, too short and  other?
   #  divide %unmapped too short  by total %unmapped  and muliply by Unmapped reads to get estimate

   if(reshape){
      n <- c(mapped,  "Unmapped reads", "Number of input reads" )
      x <- dplyr::filter(x, stat %in% n) %>%
            dplyr::mutate( stat = factor(stat, levels= n)) %>%
             tidyr::spread( stat, value)  %>%
              dplyr::mutate_each( dplyr::funs(as.integer ) , -1)
   }
   x
}



library("optparse")
option_list = list(
    make_option(c("-d", "--STAR"), type="character", default=NULL,
                help="a folder containing STAR.out logs", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


output = janitor::clean_names(read_STAR(path = opt$STAR,reshape = TRUE))

data.table::fwrite(output, opt$out)
