library(vegan)
library(tidyverse)
library(dendextend)
library(phyloseq)
library(DESeq2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(grid)



# For local testing:

# runsheet <- "249_Truncsheet.csv"
# sample_info <- "unique-sample-IDs.txt"
# counts <- "counts.tsv")
# taxonomy <- "taxonomy.tsv"

# final_outputs_dir <- "/testing_amplicon/Final_Outputs/"

# For Snakemake:
args <- commandArgs(trailingOnly = TRUE)
# Assign arguments to variables based on snakemake shell line
runsheet <- paste0(args[1])
sample_info <- paste0(args[2])
counts <- paste0(args[3])
taxonomy <- paste0(args[4])

final_outputs_dir <- paste0(args[5]) 


# Make plot output directories
dendrogram_out_dir <- paste0(final_outputs_dir, "dendrogram/")
if (!dir.exists(dendrogram_out_dir)) {
  dir.create(dendrogram_out_dir)
}

pcoa_out_dir <- paste0(final_outputs_dir, "PCoA/")
if (!dir.exists(pcoa_out_dir)) {
  dir.create(pcoa_out_dir)
}

rarefaction_out_dir <- paste0(final_outputs_dir, "rarefaction/")
if (!dir.exists(rarefaction_out_dir)) {
  dir.create(rarefaction_out_dir)
}

richness_out_dir <- paste0(final_outputs_dir, "richness/")
if (!dir.exists(richness_out_dir)) {
  dir.create(richness_out_dir)
}

taxonomy_out_dir <- paste0(final_outputs_dir, "taxonomy/")
print(taxonomy_out_dir)
if (!dir.exists(taxonomy_out_dir)) {
  dir.create(taxonomy_out_dir)
}

de_out_dir <- paste0(final_outputs_dir, "de/")
if (!dir.exists(de_out_dir)) {
  dir.create(de_out_dir)
}

# 2. Reading in processed data

runsheet <- as.data.frame(read.table(file = runsheet, 
                                     header = TRUE, sep = ",", 
                                     row.names = 1))
# Use only samples listed in sample_info
sample_names <- readLines(sample_info)

# Reorder the runsheet df to match the unique_sample_ids
# Identify the matching rows by removing suffix from basename of file

remove_suffix <- function(path, suffix) {
  file_name <- basename(path)
  sub(suffix, "", file_name)
}
matching_rows <- sapply(sample_names, function(sn) {
  which(sapply(1:nrow(runsheet), function(i) {
    remove_suffix(runsheet$read1_path[i], runsheet$raw_R1_suffix[i]) == sn
  }))
})

matching_rows <- unlist(matching_rows, use.names = FALSE)
matching_rows <- unique(matching_rows)

# Subset the runsheet
runsheet <- runsheet[matching_rows, ]

# Remove the longest common prefix from the sample names (for visualizations)
longest_common_prefix <- function(strs) {
  if (length(strs) == 1) return(strs)
  
  prefix <- strs[[1]]
  for (str in strs) {
    while (substring(str, 1, nchar(prefix)) != prefix) {
      prefix <- substr(prefix, 1, nchar(prefix) - 1)
    }
  }
  
  return(prefix)
}

remove_common_prefix <- function(strs) {
  prefix <- longest_common_prefix(strs)
  sapply(strs, function(x) substr(x, nchar(prefix) + 1, nchar(x)))
}

shortened_row_names <- remove_common_prefix(rownames(runsheet))
rownames(runsheet) <- shortened_row_names

count_tab <- read.table(file = counts, 
                        header = TRUE, row.names = 1, sep = "\t")
# Convert sample names to match those in counts table
transformed_sample_names <- make.names(sample_names, unique = TRUE)
# Keep only samples used in the sample info (should be redundant step)
count_tab <- count_tab[, transformed_sample_names]

# Keep only genes with at least 1 count
count_tab <- count_tab[rowSums(count_tab) > 0, ]
# Check if every gene has a 0 in the row, add +1 pseudocount for VST
if (all(apply(count_tab, 1, any))) {
  # Add pseudocount of 1 to the entire counts data frame
  count_tab <- count_tab + 1
}


# Move shortened sample names to sample_names
sample_names <- rownames(runsheet)
# Rename counts columns with shortened names

colnames(count_tab) <- rownames(runsheet)
tax_tab <- read.table(file = taxonomy, 
                      header = TRUE, row.names = 1, sep = "\t")
deseq_counts <- DESeqDataSetFromMatrix(countData = count_tab, 
                                       colData = runsheet, 
                                       design = ~1)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)


# Add colors to runsheet
num_colors <- length(unique(runsheet$groups))
colors <- brewer.pal(num_colors, "Set1") #or colors_list <-circlize:randcolor(num_groups)
group_colors <- setNames(colors, unique(runsheet$groups))
runsheet <- runsheet %>% mutate(color = group_colors[groups])


#op <- par(no.readonly = TRUE)
#par(op)

width_in_inches <- 11.1
height_in_inches <- 8.33
dpi <- 300
width_in_pixels <- width_in_inches * dpi
height_in_pixels <- height_in_inches * dpi

# Adjust parameters using group labels

# 3A: Hierarchical Clustering

euc_dist <- dist(t(vst_trans_count_tab))
euc_dist
euc_clust <- hclust(d = euc_dist, method = "ward.D2")


# output 1: Uncolored

png(paste0(dendrogram_out_dir, "dendrogram.png"), width = width_in_pixels, height = height_in_pixels, res = dpi)
plot(euc_clust)
dev.off()

# output 2: Dendrograms colored by group
euc_dend <- as.dendrogram(euc_clust, h = .1)
sample_info_tab <- runsheet[, c('groups', 'color')]
dend_cols <- sample_info_tab$color[order.dendrogram(euc_dend)]
labels_colors(euc_dend) <- dend_cols
png(file.path(dendrogram_out_dir, paste0("dendrogram_colored", ".png")), width = width_in_pixels, height = height_in_pixels, res = dpi)
plot(euc_dend, ylab = "VST Euc. dist.")
dev.off()
#ggsave(filename = paste0(pcoa_out_dir, "phyloseq_PCoA", ".png"), plot=pcoa_plot)
########################

# 3B Ordination

# making a phyloseq object with our transformed table
vst_count_phy <- otu_table(object = vst_trans_count_tab, taxa_are_rows = TRUE)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
vst_physeq

# generating a PCoA with phyloseq
vst_pcoa <- ordinate(physeq = vst_physeq, method = "PCoA", distance = "euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

# Calculate the percentage of variance
percent_variance <- eigen_vals / sum(eigen_vals) * 100
label_PC1 <- sprintf("PC1 [%.1f%%]", percent_variance[1])
label_PC2 <- sprintf("PC2 [%.1f%%]", percent_variance[2])

pcoa_plot <- plot_ordination(vst_physeq, vst_pcoa, color = "groups") + 
  geom_point(size = 1) + 
  labs(
    col = "groups", 
    x = label_PC1,
    y = label_PC2
  ) + labs(col = "groups") + 
  geom_text(aes(label = rownames(sample_info_tab), hjust = 0.3, vjust = -0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  scale_color_manual(values = unique(sample_info_tab$color[order(sample_info_tab$groups)])) + 
  theme_bw() + theme(legend.position = "none",  text = element_text(size = 15)) + ggtitle("PCoA")
ggsave(filename = paste0(pcoa_out_dir, "PCoA", ".png"), plot=pcoa_plot, width = 11.1, height = 8.33, dpi = 300)

#4. Alpha diversity

# 4a. Rarefaction curves

png(file = paste0(rarefaction_out_dir, "rarefaction", ".png"))
rarecurve(x = t(count_tab), step = 100, col = sample_info_tab$color, 
          lwd = 2, ylab = "ASVs", label = FALSE)
dev.off()

# 4b. Richness and diversity estimates
# create a phyloseq object similar to how we did above in step 3B, only this time also including our taxonomy table:
count_tab_phy <- otu_table(count_tab, taxa_are_rows = TRUE)
tax_tab_phy <- tax_table(as.matrix(tax_tab))
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

richness_plot <- plot_richness(ASV_physeq, color = "groups", measures = c("Chao1", "Shannon")) + 
  scale_color_manual(values = unique(sample_info_tab$color)) + 
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title.align = 0.5,
    axis.text.x = element_blank(),
    legend.title = element_blank()
  )
ggsave(paste0(richness_out_dir, "richness", ".png"), plot=richness_plot, width = 11.1, height = 8.33, dpi = 300)

richness_by_group <- plot_richness(ASV_physeq, x = "groups", color = "groups", measures = c("Chao1", "Shannon")) + 
  scale_color_manual(values = unique(sample_info_tab$color)) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title.align = 0.5,
    axis.text.x = element_blank(),
    legend.title = element_blank()
  )
ggsave(filename = paste0(richness_out_dir, "richness_by_group", ".png"), plot=richness_by_group, width = 11.1, height = 8.33, dpi = 300)

# 5. Taxonomic summaries

proportions_physeq <- transform_sample_counts(ASV_physeq, function(ASV) ASV / sum(ASV))

relative_phyla <- plot_bar(proportions_physeq, x = "groups", fill = "phylum") + 
  theme_bw() + theme(text = element_text(size = 9))
ggsave(filename = paste0(taxonomy_out_dir, "relative_phyla", ".png"), plot=relative_phyla, width = 11.1, height = 8.33, dpi = 300)

relative_classes <- plot_bar(proportions_physeq, x = "groups", fill = "class") + 
  theme_bw() + theme(text = element_text(size = 9))
ggsave(filename = paste0(taxonomy_out_dir, "relative_classes", ".png"), plot=relative_classes, width = 11.1, height = 8.33, dpi = 300)

# 6 Statistically testing for differences

# statistical significance stuff
betadisper(d = euc_dist, group = sample_info_tab$groups) %>% anova()
adonis_res <- adonis2(formula = euc_dist ~ sample_info_tab$groups)
r2_value <- adonis_res$R2[1]
prf_value <- adonis_res$`Pr(>F)`[1]

label_PC1 <- sprintf("PC1 [%.1f%%]", percent_variance[1])
label_PC2 <- sprintf("PC2 [%.1f%%]", percent_variance[2])

ordination_plot <- plot_ordination(vst_physeq, vst_pcoa, color = "groups") + 
  geom_point(size = 1) + 
  labs(
    col = "groups", 
    x = label_PC1,
    y = label_PC2
  ) + labs(col = "groups") + 
  geom_text(aes(label = rownames(sample_info_tab), hjust = 0.3, vjust = -0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  scale_color_manual(values = unique(sample_info_tab$color[order(sample_info_tab$groups)])) + 
  theme_bw() + theme(legend.position = "bottom",  text = element_text(size = 15, ),
                     legend.direction = "vertical",
                     legend.justification = "center",
                     legend.box.just = "center",
                     legend.title.align = 0.5) +
  annotate("text", x = Inf, y = -Inf, label = paste("R2:", toString(round(r2_value, 3))), hjust = 1.1, vjust = -2, size = 4)+
  annotate("text", x = Inf, y = -Inf, label = paste("Pr(>F)", toString(round(prf_value,4))), hjust = 1.1, vjust = -0.5, size = 4)+ ggtitle("PCoA")
ggsave(filename=paste0(pcoa_out_dir, "PCoA_anova", ".png"), plot=ordination_plot, width = 11.1, height = 8.33, dpi = 300)


#### pairwise comparisons

unique_groups <- unique(runsheet$groups)
deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq, design = ~groups)
deseq_modeled <- DESeq(deseq_obj)

# make the volcanoplot
plot_comparison <- function(group1, group2) {
  
  deseq_res <- results(deseq_modeled, contrast = c("groups", group1, group2))
  norm_tab <- counts(deseq_modeled, normalized = TRUE) %>% data.frame()
  
  volcano_data <- as.data.frame(deseq_res)
  
  p_val <- 0.05
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj < p_val #also logfc cutoff?
  
  # ASVs promoted in space on right, reduced on left
  p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(pvalue), color=significant)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("black", "red")) +
    theme_bw() +
    labs(title="Volcano Plot",
         x=paste("Log2 Fold Change\n(",group1," vs ",group2,")"),
         y="-Log10 P-value",
         color=paste0("Significant < ", p_val)) +
    theme(legend.position="top")
  
  # label points and plot
  top_points <- volcano_data %>%
    arrange(padj) %>%
    filter(significant) %>%
    head(10)
  
  volcano_plot <- p + geom_text_repel(data=top_points, aes(label=row.names(top_points)), size=3)
  ggsave(filename=paste0(de_out_dir,"volcano_",
                         gsub(" ", "_", group1),
                         "_vs_",
                         gsub(" ", "_", group2), ".png"),
         plot=volcano_plot,
         width = 11.1, height = 8.33, dpi = 300)
}


# setting up pairwise comparisons and running
comparisons <- expand.grid(group1 = unique_groups, group2 = unique_groups)
comparisons <- subset(comparisons, group1 != group2)

apply(comparisons, 1, function(pair) plot_comparison(pair['group1'], pair['group2']))


##########
# Extract legend from richness_by_group
# Print to new file
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend 
} 
legend <- g_legend(richness_by_group)
grid.newpage()
grid.draw(legend)
legend_filename <- paste0(final_outputs_dir, "color_legend.png")
ggsave(legend_filename, plot = legend, device = "png", width = 11.1, height = 4, dpi = 300)
