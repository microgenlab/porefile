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
                     c("-c", "--out_counts"), 
                     action = "store",
                     type = "character",
                     help = "The output tsv file with taxa counts.")

parser <- add_option(parser, 
                     c("-x", "--out_taxcla"), 
                     action = "store",
                     type = "character",
                     help = "The output tsv file with taxa classification.")

parser <- add_option(parser, 
                     c("-k", "--out_keys"), 
                     action = "store",
                     type = "character",
                     help = "The output tsv file with TAXA to ncbi_id.")

parser <- add_option(parser, 
                     c("-p", "--polish"), 
                     action = "store",
                     type = "character",
                     help = "If this scripts must compute files for later polishing step.")

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
out_counts <- opt$out_counts
out_taxcla <- opt$out_taxcla
out_keys <- opt$out_keys
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

# Create labels ("TAXA_xxxx") and write key
keys <- unique(info[c("ncbi_id", "path")])
keys <- keys[order(keys$ncbi_id), ]
lntx <- dim(keys)[1]
wdtx <- nchar(lntx)
keys$taxa <- paste0("TAXA_", formatC(seq_len(lntx), width=wdtx, format = "d", flag = "0"))
# write.table(keys, out_keys, sep="\t", quote=F)

# Write counts
counts <- table(info$ncbi_id, info$barcode)
rownames(counts) <- keys$taxa[match(rownames(counts), keys$ncbi_id)]
write.table(counts, out_counts, sep="\t", quote=FALSE)

# Parse taxon names
taxcla <- keys[c("taxa", "path")]
spl <- strsplit(setNames(taxcla$path, taxcla$taxa), ";")
rks <- c("[D]", "[P]", "[C]", "[O]", "[F]", "[G]", "[S]")
tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rks <- setNames(rks, tax)
tt <- lapply(spl, function(x) {
    unlist(lapply(rks, function(y){
	    grep(y, x, fixed = TRUE, value = TRUE)
    }))
})
tt <- lapply(tt, function(y) gsub("[ ]*\\[\\w{1,2}\\][ ]", "", y))

# Create TAX table and write
taxt <- matrix(NA_character_, 
			nrow = length(tt), 
			ncol = length(tax), 
			dimnames = list(names(tt), tax))
for (i in seq_along(tt)){
	taxt[names(tt)[i], names(tt[[i]])] <- tt[[i]]
}
write.table(taxt, out_taxcla, sep = "\t", quote = FALSE)


if (opt$polish == "true"){
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
}
