#!/usr/bin/env Rscript
# Author: Mei Wu, https://github.com/projectoriented

# Logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

# Extract arguments
paf_file_name <- snakemake@input[["oriented_ava_paf"]]
out_png_name <- snakemake@output[["ava_png"]]

param_dict <- snakemake@params[[1]]
seqnames_order <- param_dict$seqnames_order
figure_title <- param_dict$figure_title


# Load the library
library(SVbyEye)
library(ggplot2)

print(param_dict)

paf_table <- readPaf(paf.file = paf_file_name, include.paf.tags = TRUE, restrict.paf.tags = "cg")

plt <- plotAVA(
            paf.table = paf_table,
            color.by = "identity",
            min.deletion.size = 50,
            min.insertion.size = 50,
            highlight.sv = 'outline',
            seqnames.order = seqnames_order,
        )

# Add Annotation
annot_bed_path <- snakemake@input[["anno_bed_shifted"]]
#, col.names=c("chrom", "start", "end", "label", "color")

if (file.size(annot_bed_path) != 0) {
    target_annot_df <- read.table(annot_bed_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    target_annot_gr <- GenomicRanges::makeGRangesFromDataFrame(target_annot_df, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
    plt <- addAnnotation(
        ggplot.obj = plt,
        annot.gr = target_annot_gr,
        coordinate.space = 'self',
        shape = "rectangle",
        offset.annotation = TRUE,
        annotation.group = 'label',
        fill.by = 'label',
        label.by = 'label',
        color.palette = target_annot_gr$color
        )
} else {
    print("empty annotation bed file, ignoring")
}

# Add title
plt <- plt + ggtitle(figure_title)

# Write out the figure
png(filename = out_png_name, width = 20, height = 10, res = 300, units = 'in')
plt
dev.off()