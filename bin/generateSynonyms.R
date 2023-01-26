#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

parser <- OptionParser()

parser <- add_option(parser, 
    c("-t", "--taxmap"), 
    action = "store",
    type = "character",
    help = "The 'tax_ncbi-species_ssu_ref_nr99_VERSION.txt' file from SILVAdb.")

parser <- add_option(parser, 
    c("-n", "--ncbisp"), 
    action = "store",
    type = "character",
    help = "The 'taxmap_slv_ssu_ref_nr_VERSION.txt' file from SILVAdb.")

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

tms <- read.csv(opt$taxmap, sep = "\t")
head(tms)

ns <- read.csv(opt$ncbisp, sep = "\t", header = FALSE)
head(ns)

# Clean taxmap organism_name
cln <- grep("[U|u]ncultured|[E|e]nvironmental|[U|u]nidentified|[M|m]etagenome", 
        tms$organism_name)
tms$organism_name[cln] <- NA_character_

# Clean taxmap path
tms$path <- gsub("uncultured;$", "", tms$path)

# Create full path
rps <- apply(tms, 1, function(x) grepl(x[["organism_name"]], x[["path"]], fixed = TRUE))
tms$organism_name[rps] <- ""
tms$organism_name[cln] <- ""
tms$path <- paste(tms$path, tms$organism_name, sep = "")

# Identify silva indexes for occurrences
idx <- split(seq_len(dim(tms)[1]), tms$path)

# Clean ncbi taxonomy
ns$V1 <- gsub(" <\\w+>", "", ns$V1)

# Create hashmap for ncbi taxonomy
ns$tip <- vapply(strsplit(ns$V1, ";"), function(x) rev(x)[1], FUN.VALUE = NA_character_)
ee <- list2env(split(ns$V2, ns$tip), hash = TRUE)

# Find taxids 
pths <- strsplit(setNames(names(idx), names(idx)), ";")
pths <- lapply(pths, rev)
txid <- lapply(pths, function(x){
    r <- NULL
    y <- x
    while(is.null(r)){
	    a <- y[1]
	    y <- y[-1]
	    r <- ee[[a]]
	    r <- r[1]
    }
    r
})

# Repeat taxids as many times as it was assigned
rp <- mapply(rep, txid, lengths(idx))

# Sort taxids by index
ultx <- unlist(rp, use.names = FALSE)
tms$map <- ultx[order(unlist(idx, use.names = FALSE))]

df <- data.frame(Accs = paste(tms$primaryAccession, tms$start, tms$stop, sep = "."),
				Syno = tms$map)

write.table(df, 
			file = opt$out, 
			sep = "\t", 
			quote = FALSE, 
			row.names = FALSE, 
			col.names = FALSE)
        