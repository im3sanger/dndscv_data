# Extracting a slide of the human genome and of the transcript table to use for 
# the buildref tutorial (dNdScv)

setwd("/warehouse/casm/team268im/im3/dNdScv/Other_species_development_area")
segment = c(178.8e6,179.3e6) # 500kb of the genome containing PIK3CA
chr = "3"
genomefile = "genome.fa"

## Extracting the genome segment and saving it as a fasta file
gr = GenomicRanges::GRanges(chr, IRanges::IRanges(segment[1], segment[2]))
seq = as.vector(Rsamtools::scanFa(genomefile, gr))

fileConn = file("chr3_segment.fa")
writeLines(c(sprintf(">3 %0.0f:%0.0f",segment[1],segment[2]), seq), fileConn)
close(fileConn)

## Saving transcript table slice
reftable = read.table("BioMart_table_human_GRCh37.p13.txt", header=1, sep="\t", stringsAsFactors=F)
colnames(reftable) = c("gene.id","gene.name","cds.id","chr","chr.coding.start","chr.coding.end","cds.start","cds.end","length","strand")
reftable = reftable[!is.na(reftable$chr.coding.start) & (reftable$chr==chr) & (reftable$chr.coding.start>(segment[1])) & (reftable$chr.coding.end<(segment[2])), ]

# Correcting the coordinates to reflect that this is a slice
reftable$chr.coding.start = reftable$chr.coding.start-segment[1]+1
reftable$chr.coding.end = reftable$chr.coding.end-segment[1]+1

write.table(reftable, file="BioMart_human_GRCh37_chr3_segment.txt", row.names=F, col.names=T, sep="\t", quote=F)

# Running the example file
library(dndscv)

path_cds_table = system.file("extdata", "BioMart_human_GRCh37_chr3_segment.txt", package = "dndscv", mustWork = TRUE)
path_genome_fasta = system.file("extdata", "chr3_segment.fa", package = "dndscv", mustWork = TRUE)

buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "kk.rda")

library(dndscv)
segment = c(178.8e6,179.3e6) # 500kb of the genome containing PIK3CA used for this example
data("dataset_simbreast", package="dndscv")
mutations = mutations[mutations$chr=="3" & mutations$pos>segment[1] & mutations$pos<segment[2], ]
mutations$pos = mutations$pos-segment[1]+1

dndsout = dndscv(mutations, refdb="RefCDS_chr3_segment.rda", cv=NULL)


# Generate these two files (fasta slice and transcript table) for this region, 
# save them in the data directory and use them in the example.
# Then explain how these tables could be generated:
# 1. Using BioMart
# 2. Restricting the table upfront to the chosen transcripts. Changing gene name to reflect multiple isoforms (e.g. CDKN2A)
# 3. In bacteria or without splicing, simply provide start and end of CDS
# 4. Or feed a fasta with one CDS per entry
# 5. Finally explain that domains or fractions of genes can be inputted as long as they are in frame and multiple of 3