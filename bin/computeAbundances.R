#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(optparse)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))

parser <- OptionParser()

parser <- add_option(parser, 
                     c("-t", "--lowAbundanceThreshold"), 
                     action = "store",
                     type = "numeric",
                     help = "The abundance threshold to be filtered out.")

parser <- add_option(parser, 
                     c("-i", "--in_suffix"), 
                     action = "store",
                     type = "character",
                     help = "The read info extension.")

parser <- add_option(parser, 
                     c("-s", "--synonyms"), 
                     action = "store",
                     type = "character",
                     help = "The synonyms file.")

parser <- add_option(parser, 
                     c("-v", "--out_silva"), 
                     action = "store",
                     type = "character",
                     help = "The output txt file with silva headers.")

parser <- add_option(parser, 
                     c("-o", "--out_suffix"), 
                     action = "store",
                     type = "character",
                     help = "The output txt file suffix with read headers.")

if (length(commandArgs(TRUE)) == 0) {
  print_help(parser)
  quit(status = 0)
}

opt <- parse_args(parser)

in_suffix <- opt$in_suffix
synonyms <- opt$synonyms
lat <- opt$lowAbundanceThreshold
out_silva <- opt$out_silva
out_suffix <- opt$out_suffix

read_info <- list.files(pattern = paste0("[.]", in_suffix, "$"))
syno <- read.table(synonyms, colClasses = "character")

info <- lapply(setNames(read_info, sub(paste0("[.]", in_suffix, "$"), "", read_info)), read.table, sep="\t") %>%
  mapply(function(x, nm){
    colnames(x) <- c("read_id", "rank", "ncbi_id", "path")
    x$barcode <- nm
    return(x)
  }, x = ., nm = names(.), 
  SIMPLIFY = FALSE, 
  USE.NAMES = FALSE) %>%
  do.call(rbind, .)


tbl <- info %>%
  filter(rank == "S") %>% 
  pull(ncbi_id) %>%
  table() %>%
  sort() %>%
  proportions()

high <- names(which(tbl >= lat))
high_silva <- syno$V1[syno$V2 %in% high]

cat(high_silva, sep = "\n", file = out_silva)

low <- names(which(tbl < lat))
info %>%
  filter(ncbi_id %in% low | rank != "S") %>%
  select(barcode, read_id) %>%
  split(x = .$read_id, f = .$barcode) %>%
  mapply(function(x, nm){
    file <- paste0(nm, ".", out_suffix)
    cat(x, sep = "\n", file = file)
  }, x = ., nm = sub("minimap2_", "", names(.)), 
  SIMPLIFY = FALSE, 
  USE.NAMES = FALSE)

