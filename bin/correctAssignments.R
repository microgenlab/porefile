#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

parser <- OptionParser()

parser <- add_option(parser, 
    c("-n", "--new"), 
    action = "store",
    type = "numeric",
    help = "The new read assignments.")

parser <- add_option(parser, 
    c("-o", "--old"), 
    action = "store",
    type = "character",
    help = "The old read assignments.")

parser <- add_option(parser, 
    c("-c", "--corrected"), 
    action = "store",
    type = "character",
    help = "The output corrected read assignment files.")

if (length(commandArgs(TRUE)) == 0) {
  print_help(parser)
  quit(status = 0)
}

opt <- parse_args(parser)

new <- opt$new
old <- opt$old
cor <- opt$corrected

ndf <- read.table(new, sep="\t", col.names=c("read", "rank", "ncbi_id", "path"))
odf <- read.table(old, sep="\t", col.names=c("read", "rank", "ncbi_id", "path"))

odf[odf$read %in% ndf$read, ] <- ndf

write.table(odf, file = cor, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
