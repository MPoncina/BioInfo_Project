library(Seurat)
library(Matrix)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(dplyr)

options(warn = 1)  # Show warnings but don't stop for docker

matrix_file <- "matrix.mtx.gz"
features_file <- "features.tsv.gz"
barcodes_file <- "barcodes.tsv.gz"

#STEP 1
sparse_mat <- readMM(matrix_file)
features <- fread(features_file, header = FALSE)
barcodes <- fread(barcodes_file, header = FALSE)
head(features)

#features = genes (rows), barcodes = cells (columns)
rownames(sparse_mat) <- features$V1
colnames(sparse_mat) <- barcodes$V1
# Convert to full matrix
full_mat <- as.matrix(sparse_mat)
# Convert to data.table
dt_matrix <- as.data.table(full_mat, keep.rownames = "Gene")
#write.table(dt_matrix, "full_matrix.tsv", sep = "\t", row.names = TRUE) if want to save


#STEP 2
# Separate gene expression data and ATAC
gene_expression_dt <- dt_matrix[grepl("^ENSG", Gene)]
atac_peak_dt <- dt_matrix[grepl("^chr[0-9XYM]+:\\d+-\\d+", Gene)]
# missing some rows to reach 117757
# Rows not in genes or peaks
missing_rows <- dt_matrix[!grepl("^ENSG", Gene) & !grepl("^chr[0-9XYM]+:\\d+-\\d+", Gene)]
missing_rows[, .N]   # Number of missing rows
missing_rows[, Gene] # Check what these missing features are
#If I want to include them:
#atac_peak_dt <- dt_matrix[grepl("^[A-Za-z0-9._]+:\\d+-\\d+", Gene)]


#STEP 3
# Total expression per gene
total_expression <- rowSums(gene_expression_dt[, -1, with = FALSE])
# Total chromatin accessibility per peak region
total_accessibility <- rowSums(atac_peak_dt[, -1, with = FALSE])
# Attach names
names(total_expression) <- gene_expression_dt$Gene
names(total_accessibility) <- atac_peak_dt$Gene


#STEP 4
peak_features <- features[grepl("^chr[0-9XYM]+:\\d+-\\d+", features$V1)]
# Create GRanges object for accessibility
peak_gr <- GRanges(
  seqnames = peak_features$V4,
  ranges = IRanges(start = as.numeric(peak_features$V5),
                   end = as.numeric(peak_features$V6)),
  accessibility = total_accessibility
)
#peak_gr <- reduce(peak_gr, with.revmap = TRUE) it is the same lenght
#saveRDS(peak_gr, "peak_gr.rds") takes too long!

gene_features <- features[features$V3 == "Gene Expression"] #is the same as grep
#add "chr" because I have some NA in V4
gene_features$V4 <- as.character(gene_features$V4)
gene_features$V4 <- ifelse(grepl("^chr", gene_features$V4),
                           gene_features$V4, paste0("chr", gene_features$V4))
# Create GRanges object for gene expression
gene_gr <- GRanges(
  seqnames = gene_features$V4,
  ranges = IRanges(start = as.numeric(gene_features$V5),
                   end = as.numeric(gene_features$V6)),
  strand = "*",
  gene_id = gene_features$V1,
  gene_symbol = gene_features$V2,
  expression = total_expression
)

#saveRDS(gene_gr, "gene_gr.rds") too long!

#STEP 5
gtf <- import("Homo_sapiens.GRCh38.114.gtf.gz")
gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding"]
# Find overlaps (each peak may overlap 0 or 1 genes)
#missing the name chr in seqlevels on gtf so I need to add it
seqlevels(gtf_genes) <- paste0("chr", seqlevels(gtf_genes))

overlaps <- findOverlaps(peak_gr, gtf_genes)
overlaps
#Remap the ATAC-seq GenomicRanges to this object and attach the summarized peak data from step 4???
# For each gene, sum the accessibility of overlapping peaks
# Get indices to extract which peak overlaps which gene:
peak_idx <- queryHits(overlaps)
gene_idx <- subjectHits(overlaps)
# Create a data.table with gene and peak info
dt_gene_peak <- data.table(
  gene_id = gtf_genes$gene_id[gene_idx],
  accessibility = mcols(peak_gr)$accessibility[peak_idx]
)
# Sum accessibility per gene
peak_summary <- dt_gene_peak[, .(total_peak_accessibility = sum(accessibility)), by = gene_id]
peak_summary
# Match peak_summary to gtf_genes by gene_id
match_idx <- match(gtf_genes$gene_id, peak_summary$gene_id)
# Create a new metadata column
mcols(gtf_genes)$total_peak_accessibility <- peak_summary$total_peak_accessibility[match_idx]
# Optional: set NA to 0 for genes with no overlapping peaks
mcols(gtf_genes)$total_peak_accessibility[is.na(mcols(gtf_genes)$total_peak_accessibility)] <- 0


#Step 6
protein_coding_ids <- gtf_genes$gene_id
# Keep only entries in gene_gr whose gene_id is protein-coding
gene_gr_pc <- gene_gr[mcols(gene_gr)$gene_id %in% protein_coding_ids]
# Create a lookup table from GTF: gene_id -> gene_name
gtf_lookup <- data.table(
  gene_id = gtf_genes$gene_id,
  gene_symbol = gtf_genes$gene_name
)
# Match by gene_id
lookup_idx <- match(mcols(gene_gr_pc)$gene_id, gtf_lookup$gene_id)
# Add gene symbols to metadata
mcols(gene_gr_pc)$gene_symbol <- gtf_lookup$gene_symbol[lookup_idx]
gene_gr_pc

#STEP 7
peak_dt <- as.data.table(peak_gr)
peak_dt[, accessibility := mcols(peak_gr)$accessibility]
peak_dt[, chr := as.character(seqnames)]
# Boxplot of total accessibility per chromosome from atac
p <- ggplot(peak_dt, aes(x = chr, y = accessibility)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  coord_cartesian(ylim = quantile(peak_dt$accessibility, c(0.05, 0.95), na.rm = TRUE)) +
  theme_minimal() +
  labs(title = "ATAC Peak Accessibility per Chromosome",
       x = "Chromosome", y = "Accessibility Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

# Extract matrix: assume gene_expression_dt includes genes as rows, samples as columns
expr_mat <- as.matrix(gene_expression_dt[, -1, with = FALSE])  # remove gene ID column
# CPM normalization
expr_cpm <- log2((t(t(expr_mat) / colSums(expr_mat)) * 1e6) + 1)
# Assume atac_peak_dt is similar: peaks as rows, columns as cells/samples
atac_mat <- as.matrix(atac_peak_dt[, -1, with = FALSE])
# CPM normalization
atac_cpm <- log2((t(t(atac_mat) / colSums(atac_mat)) * 1e6) + 1)

# Create data.table for merging
expr_dt <- data.table(
  gene_id = mcols(gene_gr_pc)$gene_id,
  expr_cpm = mcols(gene_gr_pc)$expression,
  chr = as.character(seqnames(gene_gr_pc))
)

atac_dt <- data.table(
  gene_id = gtf_genes$gene_id,
  atac_cpm = mcols(gtf_genes)$total_peak_accessibility
)

# Merge by gene_id
merged_dt <- merge(expr_dt, atac_dt, by = "gene_id",all=TRUE)

used_peak_ids <- unique(queryHits(overlaps))
unmerged_peak_ids <- setdiff(seq_along(peak_gr), used_peak_ids)

unmerged_peaks_dt <- as.data.table(peak_gr[unmerged_peak_ids])
unmerged_summary <- unmerged_peaks_dt[, .N, by = seqnames]
unmerged_summary

p1 <-ggplot(unmerged_summary, aes(x = seqnames, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Unmerged ATAC Peaks per Chromosome", x = "Chromosome", y = "Count")

#missing the genes with no atac correlation but this gives 0 output
genes_without_atac <- merged_dt[is.na(atac_cpm)]

merged_dt <- merged_dt[!is.na(expr_cpm) & !is.na(atac_cpm)]
# Apply CPM normalization and log2 transform
merged_dt[, expr_log := log2(expr_cpm + 1)]
merged_dt[, atac_log := log2(atac_cpm + 1)]

p2<-ggplot(merged_dt, aes(x = expr_log, y = atac_log)) +
  geom_point(alpha = 0.3, color = "steelblue", size = 1) +
  theme_minimal() +
  labs(title = "Expression vs ATAC Accessibility",
       x = "log2(CPM Expression + 1)",
       y = "log2(CPM ATAC + 1)")
p2 #crowded plot to split in chrs
# Filter to 24 chromosomes only
merged_dt_chr <- merged_dt[chr %in% paste0("chr", c(1:22, "X", "Y"))]

# Order chromosomes
merged_dt_chr[, chr := factor(chr, levels = paste0("chr", c(1:22, "X", "Y")))]

# Plot with facet_wrap
p3<-ggplot(merged_dt_chr, aes(x = expr_log, y = atac_log)) +
  geom_point(alpha = 0.3, color = "steelblue", size = 0.8) +
  theme_minimal() +
  facet_wrap(~chr, scales = "free") +
  labs(title = "Expression vs ATAC by Chromosome",
       x = "log2(CPM Expression + 1)",
       y = "log2(CPM ATAC + 1)")
p3
ggsave("expression_vs_atac_by_chr.png", plot = plot, width = 10, height = 8, dpi = 300, bg = "white")
save.image(file = "my_workspace.RData")
