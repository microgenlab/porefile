#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

parser <- OptionParser()

parser <- add_option(parser, 
    c("-s", "--suffix"), 
    action = "store",
    type = "character",
    help = "A string with a suffix (extension) to detect info files in the current directory.")

parser <- add_option(parser, 
    c("-w", "--workflow"), 
    action = "store",
    type = "character",
    help = "A string with the workkflow selected.")

if (length(commandArgs(TRUE)) == 0) {
  print_help(parser)
  quit(status = 0)
}

opt <- parse_args(parser)
	
# Read *.info files
pattern <- paste0("[.]", opt$suffix, "$")
infos <- list.files(pattern = pattern)
lp <- lapply(setNames(infos, infos), function(x) {
a <- try(read.csv(x, sep = "\t", header = FALSE))
if (class(a) != "try-error"){
	setNames(a$V2, a$V1)
}else{
	NULL
}
})

# Get present taxa
lvls <- unique(unlist(lapply(lp, names)))

# Create new names (TAXA_*)
lnlv <- length(lvls)
wdth <- nchar(lnlv)
bc <- paste0("TAXA_", formatC(seq_along(lvls),
							width = wdth,
							format = 'd',
							flag = '0'))

# Parse tax paths for each TAXA
rr <- data.frame(otu = bc, raw = lvls)
spl <- strsplit(setNames(rr$raw, rr$otu), ";")
rks <- c("[D]", "[P]", "[C]", "[O]", "[F]", "[G]", "[S]")
tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rks <- setNames(rks, tax)
tt <- lapply(spl, function(x) {
    unlist(lapply(rks, function(y){
	    grep(y, x, fixed = TRUE, value = TRUE)
    }))
})
tt <- lapply(tt, function(y) gsub("[ ]*\\[\\w{1,2}\\][ ]", "", y))

# Create TAX table
taxt <- matrix(NA_character_, 
			nrow = length(tt), 
			ncol = length(tax), 
			dimnames = list(names(tt), tax))
for (i in seq_along(tt)){
	taxt[names(tt)[i], names(tt[[i]])] <- tt[[i]]
}


re <- as.data.frame(cbind(rr, taxt))

# Rename OTU counts
lp <- lapply(lp, function(x){
	if (length(x)){
		setNames(x, re$otu[match(names(x), re$raw)])
	}
})

# Create "OTU" (TAXA) Table
otut <- matrix(0L, nrow = length(re$otu), ncol = length(lp),
			dimnames = list(re$otu, names(lp)))
for (i in seq_along(lp)){
    otut[names(lp[[i]]), names(lp)[i]] <- lp[[i]]
}
colnames(otut) <- sub("[.]info$", "", colnames(otut))


# Write Taxa Classification Table and Counts Table
write.table(taxt, paste0(opt$workflow, "_TAXCLA.tsv"), quote = FALSE, sep = "\t")
write.table(otut, paste0(opt$workflow, "_COUNTS.tsv"), quote = FALSE, sep = "\t")