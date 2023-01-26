#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

parser <- OptionParser()

parser <- add_option(parser, 
    c("-f", "--prefix_filt_scrubb"), 
    action = "store",
    type = "character",
    help = "A string with a prefix to unequivocally detect filtered and scrubbed files in the curren directory. 
    NOT A GLOB!!")

parser <- add_option(parser, 
    c("-r", "--prefix_raw"), 
    action = "store",
    type = "character",
    help = "A string with a prefix to unequivocally detect raw files in the current directory. NOT A GLOB!!")

parser <- add_option(parser, 
    c("-o", "--out"), 
    action = "store",
    type = "character",
    help = "The output file name.")

if (length(commandArgs(TRUE)) == 0) {
  print_help(parser)
  quit(status = 0)
}

opt <- parse_args(parser)


filt_files <- list.files(pattern = opt$prefix_filt_scrubb)
raw_files <- list.files(pattern = opt$prefix_raw)
filt_files <- setNames(filt_files, gsub(opt$prefix_filt_scrubb, "", tools::file_path_sans_ext(filt_files)))
raw_files <- setNames(raw_files, gsub(opt$prefix_raw, "", tools::file_path_sans_ext(raw_files)))

lp <- lapply(list(RawReads = raw_files, Filtered = filt_files), function(x) {
	vapply(x, readLines, FUN.VALUE=NA_character_)
})

x <- do.call(cbind, lapply(lp, function(x) { x[match(names(lp[[1]]), names(x))] }) )
x <- as.data.frame(x)
x$BarCode <- rownames(x)
x <- x[, c("BarCode", "RawReads", "Filtered")]

write.table(x, file = opt$out, sep="\t", quote=F, row.names=FALSE, col.names=TRUE)