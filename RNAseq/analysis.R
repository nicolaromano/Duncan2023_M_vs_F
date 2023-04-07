library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)
library(genefilter)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(pheatmap)
library(goseq)
library(GOfuncR)

# This assumes bam files are in the current directory
bamfilenames <- dir(".", "bam$")
bfl <- BamFileList(bamfilenames, yieldSize = 5e6)

# Exons by gene
annot <- makeTxDbFromGFF("../Mus_musculus.GRCm38.98.gtf", format="gtf")
exonsByGene <- exonsBy(annot, by="gene");
overlaps <- summarizeOverlaps(exonsByGene, bfl, singleEnd = FALSE, ignore.strand = TRUE, fragments = TRUE)
saveRDS(overlaps, "overlaps.rds")

# Convert Ensembl ids to gene names
# We can use listDatasets() to see what datasets are available
# listDatasets(useEnsembl("ensembl"))

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
# Similarly, we can use listAttributes() to see what attributes are available to pass to getBM
# grep("mgi", listAttributes(mart)$name, value = T)

genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
               values = rownames(overlaps), mart = mart)

genes <- genes[-which(duplicated(genes$ensembl_gene_id)),]
saveRDS(genes, "genenames.rds")

# We can directly load the data from the files and start from here
# overlaps <- readRDS("overlaps.rds")
# genes <- readRDS("genenames.rds")

se <- SummarizedExperiment()
# Construct metadata
treat <- factor(rep(c("F.CTRL", "M.CTRL"), each = 3))
# If we want to use males as the reference, we can do this
treat <- relevel(treat, ref = "M.CTRL")

# Build the metadata
exp.meta <- data.frame(
  run = colnames(overlaps),
  treatment = treat,
  replicate = factor(1:3)
)

colData(overlaps) <- cbind(colData(overlaps), exp.meta)

# Plot total reads in each group
totalreads <- data.frame(
  reads = colSums(assay(overlaps)),
  treatment = colData(overlaps)$treatment
)

ggplot(totalreads, aes(treatment, log10(reads))) +
  geom_point() +
  ylim(6, 8) +
  ylab(expression(log[10]~"(total reads)")) +
  theme_bw() +
  xlab("") +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18)
        )

# Define the design
dds <- DESeqDataSet(overlaps, design = ~ replicate + treatment)

colnames(dds) <- (c(paste("Female", 1:3, sep = "_"), 
                    paste("Male", 1:3, sep = "_")))

########### FILTERING ############

nrow(dds) # ~55k genes
# Only keep genes with a count of 10 or higher in at least 2 samples
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
nrow(dds) # ~18.6k genes left

saveRDS(dds, "dds_MF.rds")

# Apply variance stabilising transformation to the data
# Note that this should NOT be used for DE analysis, but only for initial
# exploratory analyses, such as PCA
dds.vst <- vst(dds, blind = FALSE)

all_counts <- data.frame(assay(dds))

colnames(all_counts) <- paste(rep(c("Female", "Male"), each = 3), 1:3, sep = "_")
all_counts$ENSEMBL_ID <- rownames(all_counts)
all_counts$MGI_symbol <- genes$mgi_symbol[match(all_counts$ENSEMBL_ID, genes$ensembl_gene_id)]
all_counts$Description <- genes$description[match(all_counts$ENSEMBL_ID, genes$ensembl_gene_id)]
write.csv(all_counts, "MF_all_data_raw_counts.csv", row.names = FALSE)

# Calculate euclidean distances between samples
sampleDists <- dist(t(assay(dds.vst)))
# Cluster samples
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds.vst$treatment, dds.vst$replicate, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(blues9))(255)

png("sample_distances.png", width = 10, height = 10, units = "in", res = 300)
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  fontsize_row = 16
)
dev.off()

# PCA
png("RNAseq_sample_PCA.png", width = 10, height = 10, units = "in", res = 300)
plotPCA(dds.vst, intgroup = c("treatment")) +
theme(axis.title = element_text(size = 18),
      axis.text = element_text(size = 16))
dev.off()

## DIFFERENTIAL EXPRESSION ANALYSIS

# Run DESeq2 for DE gene analysis
dds <- DESeq(dds)

# Coefficient 4 is MvF (1 is intercept, 2 and 3 are between replicates so we ignore them)
# SAVE DE GENES
contrast <- resultsNames(dds)[4]
print(paste("Contrast:", contrast))
res <- results(dds, name = contrast, alpha = .05)
res <- lfcShrink(dds,
  res = res,
  coef = contrast, type = "apeglm"
)

# out of 18631 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1332, 7.1%
# LFC < 0 (down)     : 1148, 6.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 723, 3.9%
# (mean count < 8)#
summary(res)

res$ensembl.id <- rownames(res)
# Append gene names and descriptions to the table
res$gene.symbol <- genes$mgi_symbol[match(rownames(res), genes$ensembl_gene_id)]
res$gene.description <- genes$description[match(rownames(res), genes$ensembl_gene_id)]

# Filter for significant genes (adjusted p-value < 0.05)
siggenes005 <- data.frame(res) %>%
  subset(padj < 0.05) %>%
  arrange(log2FoldChange)

write.csv(siggenes005, file = "FvM_p_less_005.csv", row.names = F)

# All genes arranged by log2 fold change, whether significant or not
allgenes <- data.frame(res) %>%
  arrange(log2FoldChange)

write.csv(allgenes, file = "FvM_all_genes.csv", row.names = F)


# MA plots
# plotMA(res, ylim = c(-6, 6))

# A nicer MA plot using ggplot2
pdf("MA_plot_all_genes_padj_less_0.05.pdf", width = 15, height = 10)
p <- ggplot(allgenes, aes(log(baseMean), log2FoldChange)) +
  geom_point(aes(col = (padj > 0.05 | is.na(padj))), size = 1) +
  scale_color_manual(values = c("#b83122", "#b8b8b8")) +
  geom_hline(yintercept = 0, lty = "dotted") +
  ylim(-11, 11) +
  xlab("Average expression") +
  ylab(expression("Average" ~ log[2] ~ "fold change")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  annotate("text",
    x = -Inf, y = -Inf,
    hjust = -0.5, vjust = -2,
    label = "Enriched in males",
    size = 5
  ) +
  annotate("text",
    x = -Inf, y = Inf,
    hjust = -0.5, vjust = 2,
    label = "Enriched in females",
    size = 5
  ) +
  ggtitle("MA plot - all genes") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20)
  )
print(p)
dev.off()

# A little handy function to plot MA plots for genes associated with a specific GO term
MAplotGO <- function(genes, go, go_terms) {
  res <- lapply(go, function(x) {
    ifelse(length(grep(paste0("(", paste(go_terms, collapse = "|"), ")"), x)),
      TRUE, FALSE
    )
  })

  res <- do.call("rbind", res)
  genes$hasGO <- res

  # Filter for genes associated with the specified GO terms
  go_rows <- which(genes[match(allgenes$gene.symbol, genes$mgi_symbol), ]$hasGO == TRUE)
  go_subset <- allgenes[go_rows, ]

  gonames <- paste(get_names(go_terms)$go_name, collapse = ", ")
  goids <- paste(get_names(go_terms)$go_id, collapse = ", ")

  print(paste("GO terms:", gonames))
  print(paste("Number of genes:", nrow(go_subset)))
  print(paste("Number of genes with p < 0.05:", nrow(subset(go_subset, padj < 0.05))))
  print(paste(
    "Number of genes with p < 0.05 and abs(log2FC) > 1:",
    nrow(subset(go_subset, padj < 0.05 & abs(log2FoldChange) > 1))
  ))

  p <- ggplot(go_subset, aes(log(baseMean), log2FoldChange)) +
    geom_point(aes(col = (padj > 0.05 | is.na(padj))), size = 1.5) +
    scale_color_manual(values = c("#801F15", "#CCCCCC")) +
    geom_hline(yintercept = 0, lty = "dotted") +
    ylim(
      -max(abs(go_subset$log2FoldChange)),
      max(abs(go_subset$log2FoldChange))
    ) +
    xlab("Average expression") +
    ylab(expression("Average" ~ log[2] ~ "fold change")) +
    annotate("text",
      x = -Inf, y = -Inf,
      hjust = -0.5, vjust = -2,
      label = "Enriched in males",
      size = 5
    ) +
    annotate("text",
      x = -Inf, y = Inf,
      hjust = -0.5, vjust = 2,
      label = "Enriched in females",
      size = 5
    ) +
    geom_text_repel(
      data = subset(go_subset, padj < 0.05),
      aes(label = gene.symbol), size = 4, max.overlaps = 20
    ) +
    ggtitle(paste("MA plot", goids)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 20)
    )

  print(p)
}

go <- getgo(genes = genes$ensembl_gene_id, id = "ensGene", genome = "mm9")

pdf("MA_plot_ion_channels_padj_less_0.05.pdf", width = 15, height = 10)
# Plot MA plot for GO terms:
# GO:0005216 - ion channel activity
MAplotGO(genes, go, go_terms = c(
  "GO:0005216"
))
dev.off()

pdf("MA_plot_steroids_padj_less_0.05.pdf", width = 15, height = 10)
# Plot MA plot for GO terms:
# GO:0008207 - steroid hormone metabolic process
# GO:0035929 - steroid hormone secretion
# GO:0090030 - regulation of steroid hormone biosynthetic process
# GO:0048545 - response to steroid hormone
MAplotGO(genes, go, go_terms = c(
  "GO:0008207", "GO:0035929",
  "GO:0090030", "GO:0048545"
))
dev.off()

# Filter for ion channels
channels_rows <- which(genes[match(allgenes$gene.symbol, genes$mgi_symbol), ]$Channel == TRUE)
channels <- allgenes[channels_rows, ]

channels %>%
  subset(padj < 0.05) %>%
  write.csv("FvM_p_less_005_channels.csv", row.names = FALSE)

# Heatmap of significantly DE ion channels
# z-scored expression data
chan_005 <- channels %>%
  subset(padj < 0.05)

expr_z_score <- t(apply(counts(dds)[ensembl_chan_005$ensembl.id, ], 1, function(x) {
  scale(x, center = TRUE, scale = TRUE)
}))

# Change to MGI symbols
rownames(expr_z_score) <- chan_005$gene.symbol
colnames(expr_z_score) <- 1:6
annot <- data.frame(Sex = rep(c("F", "M"), each = 3))

pdf("heatmap_channels_padj_less_0.05.pdf", width = 8, height = 14)
pheatmap(expr_z_score, cluster_cols = FALSE, annotation_col = annot, 
  annotation_colors = list(Sex = c("M" = "#44a4e4", "F" = "#ff8fc9")),
  show_colnames = FALSE, fontsize_row = 12, cutree_rows = 2, cellheight = 12,
  treeheight_row = 100, fontsize = 15)
dev.off()

