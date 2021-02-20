# File name: Finalize.R
# Last Edit: 03/01/2020
#
# Zaro BW, Noh JJ, et al. 2020
# 
# https://github.com/jnoh4/PofHemat
#
# Description:
# Converts mirDB target list into a list of files for each miRNA
# with targets represented as entrezID.

if(!('ggfortify' %in% rownames(installed.packages()))) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ggfortify")
}
library(ggfortify)
source('./pathways.R')
source('./data.R')

# List of miRNAs pulled
miRDB_list <- read.table('./miRDB_mmu_list_publication.txt', header = TRUE)
mmu_miRNAs <- sapply(unlist(miRDB_list), function(x){
  file_name <- paste('./miRDB mmu_publication/', x, '.txt', sep = '')
  return(unique(read.table(file_name, header = TRUE)))
})
names(mmu_miRNAs) <- unlist(miRDB_list)
gene_set <- mmu_miRNAs
miRNAs = as.character(unlist(gene_set))

# Color and order specification
pref_order <- c('AgedHSC', 'AgedMPPa', 'AgedMPPb', 'AgedMPPc', 'HSC', 'MPPa',
                'MPPb', 'MPPc', 'CLP', 'CMP', 'GMP', 'MEP')
HSPC_order <- c('HSC', 'MPPa', 'MPPb', 'MPPc')
my_colors <- c(HSC = "#FC6565", MPPa = "#F7AD14",
               MPPb = "#27A01A", MPPc = "#4ED7B9",
               CLP = '#118EFE', CMP = '#38026F',
               GMP = '#E20E40', MEP = '#E078FC',
               AgedH = '#780677', AgedA = '#985697', 
               AgedB = '#b886b7', AgedC = '#d8a6d7')

# Raw entrez protein and mRNA
raw_df <- read.csv('3_entrez_intensity_raw.csv', row.names = 1, header = TRUE)
raw_order <- sapply(colnames(raw_df), function(x) {
  sub('[0-9].+' , '', x)
})
raw_df_r <- read.table('mRNA_entrez_norm.txt', header = TRUE, row.names = 1)
raw_r_order <- sapply(colnames(raw_df_r), function(x) {
  sub('f', 'c', sub('_.+', '', x))
})
# Norm entrez protein
norm_df <- data.frame(sapply(raw_df, function(x) {
  1000000 * x / sum(x)
}))
rownames(norm_df) <- rownames(raw_df)

# Combined entrez protein and mRNA
p_df <- read.csv('3_entrez_intensity_comb.csv', row.names = 1, header = TRUE)
p_df <- p_df[, order(match(colnames(p_df), pref_order))]
r_df <- read.table('mRNA_entrez_comb.txt', row.names = 1, header = TRUE)
colnames(r_df)[colnames(r_df) == 'MPPf'] <- 'MPPc'
r_df <- r_df[, order(match(colnames(r_df), colnames(p_df)))]

# Combined gene protein
comb_df <- read.csv('3_gene_intensity_comb.txt', header = TRUE, row.names = 2)
comb_df$X <- NULL
comb_df <- comb_df[,order(match(colnames(comb_df), pref_order))]

# Mapped entrez to protein (use for translating; also used before for mRNA to entrez)
maps <- read.table('mRNA_entrez_filter.txt', sep = ',', header = TRUE)
maps <- maps[!is.na(maps$entrez),]
rownames(maps) <- maps$gene

### Revision ###
################

yh_tf <- comb_df[, c('AgedHSC', 'HSC')] > 0

# Expressed in AgedHSC but not adult HSC
eAgedNEAdult <- yh_tf[, 'AgedHSC'] * (!yh_tf[, 'HSC'])
genes_ENE <- rownames(yh_tf)[as.logical(eAgedNEAdult)]

# Expressed in at least three other cell types besides AgedHSCs
norm_df <- read.table('3_gene_intensity_norm.txt', sep = ',', heade = TRUE, row.names = 2)
norm_df$X <- NULL
ndf_cells <- sapply(colnames(norm_df), function(x){
  return(sub('[0-9].+', '', x))
})
norm_no_HSC <- norm_df[, ndf_cells != 'HSC'] > 0
norm_aHSC <- norm_df[, ndf_cells == 'AgedHSC'] > 0
nnH_ge3 <- apply(norm_no_HSC, 1, sum) >= 3
nnH_ge3 <- rownames(norm_no_HSC)[nnH_ge3]
naH_ge2 <- apply(norm_aHSC, 1, sum) >= 3
naH_ge2 <- rownames(norm_aHSC)[naH_ge2]

overlap3 <- intersect(genes_ENE, intersect(nnH_ge3, naH_ge2))
yAgednAdult <- comb_df[overlap3,c('AgedHSC', 'HSC')]

write.table(yAgednAdult, 'yesAgednoAdult.csv', sep = ',', row.names = T, col.names = T)

## Fold change higher in aged
# Expressed in both, at least 3 times in aged
exp_3 <- comb_df[intersect(rownames(comb_df)[as.logical((comb_df$AgedHSC > 0) * (comb_df$HSC > 0))], naH_ge2),c('AgedHSC', 'HSC')]
e3fc <- exp_3$AgedHSC / exp_3$HSC
e3fc_rank <- rank(e3fc)
of_interest <- e3fc[order(-e3fc_rank)][1:floor(length(e3fc)/40)]
names(of_interest) <- rownames(exp_3)[order(-e3fc_rank)][1:floor(length(e3fc)/40)]

write.table(of_interest, 'topAgedFold.csv', sep = ',', row.names = T, col.names = T)

# Combined list for miRNA
fold_filter <- rownames(exp_3)[order(-e3fc_rank)][1:floor(length(e3fc)/40)]
fold_entrez <- maps$entrez[which(maps$gene %in% fold_filter)]
binary_filter <- overlap3
binary_entrez <- maps$entrez[which(maps$gene %in% binary_filter)]
comb_filter <- union(fold_filter, binary_filter)
comb_entrez <- maps$entrez[which(maps$gene %in% comb_filter)]

miRNAs_ranked <- sapply(mmu_miRNAs, function(x) {
  length(intersect(x, comb_entrez))
})
miRNAs_sorted <- sort(miRNAs_ranked, decreasing = TRUE)

write.table(miRNAs_sorted, 'AgedmiRNA_filtered.csv', sep = ',', row.names = T, col.names = T)

# Higher aged miRNA detect
miRNAs_detected <- sapply(fold_entrez, function(x) {
  miRNA_TF <- sapply(mmu_miRNAs, function(y){
    return(x %in% y)
  })
  return(paste(names(miRNA_TF)[miRNA_TF], collapse = ', '))
})
fold_e_g <- maps$gene[which(maps$entrez %in% fold_entrez)]
names(miRNAs_detected) <- fold_e_g
write.table(miRNAs_detected, 'fold_miRNA_list.csv', sep = ',', row.names = T, col.names = T)

# Presence aged miRNA detect
miRNAs_detected2 <- sapply(binary_entrez, function(x) {
  miRNA_TF <- sapply(mmu_miRNAs, function(y){
    return(x %in% y)
  })
  return(paste(names(miRNA_TF)[miRNA_TF], collapse = ', '))
})
binary_e_g <- maps$gene[which(maps$entrez %in% binary_entrez)]
length(binary_entrez) == length(binary_e_g)
names(miRNAs_detected2) <- binary_e_g
write.table(miRNAs_detected2, 'binary_miRNA_list.csv', sep = ',', row.names = T, col.names = T)

# miR29a-3p expression
mir29a_genes <- maps$gene[which(maps$entrez %in% mmu_miRNAs$'mmu-miR-29a-3p')]

write.table(mir29a_genes, 'miR29a-3p_targets.csv', row.names = F, col.names = F)

aged_p_overlap_29 <- intersect(rownames(comb_df)[comb_df$AgedHSC > 0], mir29a_genes)
young_p_overlap_29 <- intersect(rownames(comb_df)[comb_df$HSC > 0], mir29a_genes)
length(aged_p_overlap_29) / sum(comb_df$AgedHSC > 0)
length(young_p_overlap_29) / sum(comb_df$HSC > 0)
young_p_overlap_29_Ma <- intersect(rownames(comb_df)[comb_df$MPPa > 0], mir29a_genes)
young_p_overlap_29_GMP <- intersect(rownames(comb_df)[comb_df$GMP > 0], mir29a_genes)
values_29_Ma <- comb_df[as.character(young_p_overlap_29_Ma), 'MPPa']
values_29_GMP <- comb_df[as.character(young_p_overlap_29_GMP), 'GMP']
values_29 <- (comb_df[as.character(young_p_overlap_29), 'HSC'])

Ma_29a <- data.frame(val = values_29_Ma, gene = young_p_overlap_29_Ma, cell = 'MPPa')
GMP_29a <- data.frame(val = values_29_GMP, gene = young_p_overlap_29_GMP, cell = 'GMP')
HSC_29a <- data.frame(val = values_29, gene = young_p_overlap_29, cell = 'HSC')
cells_29a <- do.call(rbind, list(Ma_29a, GMP_29a, HSC_29a))

write.table(cells_29a, 'miR29a-3p_cells.csv', sep = ',', row.names = TRUE, col.names = TRUE)

ggplot(cells_29a, aes(x = log2(val), color = cell, fill = cell)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha = 0.5, binwidth = 0.5) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 0.35))

cells_mean <- c(mean(Ma_29a$val), mean(GMP_29a$val), mean(HSC_29a$val))
cells_median <- c(median(Ma_29a$val), median(GMP_29a$val), median(HSC_29a$val))
cells_order <- c('MPPa', 'GMP', 'HSC')

cells_summary <- data.frame(cells_mean, cells_median, cells_order)
write.table(cells_summary, 'miR29a-3p_cells_sum.csv', sep = ',', row.names = T, col.names = T)

values_29 <- data.frame(values_29)
values_29$typ <- 'miR'
colnames(values_29) <- c('val', 'typ')
values_29$genes <- as.character(young_p_overlap_29)

values_tot <- (comb_df[as.character(rownames(comb_df)[comb_df$HSC > 0]), 'HSC'])
values_tot <- data.frame(values_tot)
values_tot$typ <- 'tot'
colnames(values_tot) <- c('val', 'typ')
values_tot$genes <- as.character(rownames(comb_df)[comb_df$HSC > 0])

values_plot <- rbind(values_29, values_tot)

write.table(values_plot, 'miR29a-3p_young.csv', sep = ',', row.names = TRUE, col.names = TRUE)

ggplot(values_plot, aes(x = log2(val), color = typ, fill = typ)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha = 0.5, binwidth = 0.5) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 0.35))

col_29 <- c('miR', 'all_hsc')
means <- c(mean(values_29$val), mean(values_tot$val))
medians <- c(median(values_29$val), median(values_tot$val))
mir29.data <- data.frame(means, medians, col_29)
write.table(mir29.data, 'miR29a-3p_stats.csv', sep = ',', row.names = T, col.names = T)

################
### Revision ###
# 
# ### Analyzing data ###
# # Reduce to genes within both DATASETS
# common <- intersect(rownames(p_df), rownames(r_df))
# p_common <- p_df[which(rownames(p_df) %in% common), ]
# r_common <- r_df[which(rownames(r_df) %in% common), ]
# r_common <- r_common[order(match(rownames(r_common), rownames(p_common))),]
# 
# # Genes in both datasets for all four cell types
# all_ex <- ((rowSums((p_common[,HSPC_order] > 0)) == 4) * (rowSums((r_common[,HSPC_order] > 0)) == 4)) == 1
# p_all <- p_common[all_ex, HSPC_order]
# r_all <- r_common[all_ex, HSPC_order]
# 
# # Genes within each dataset overlapping with miRNA targets
# m_r_common <- intersect(miRNAs, rownames(r_df))
# m_p_common <- intersect(miRNAs, rownames(p_df))
# 
# # RNA expressed by each HSPC cell type that are miRNA targets
# m_targets <- sapply(HSPC_order, function(x) {
#   rownames(r_df[m_r_common,])[r_df[m_r_common, x] > 0]
# })
# m_targets_c <- sapply(HSPC_order, function(x) {
#   length(rownames(r_df[m_r_common,])[r_df[m_r_common, x] > 0])
# })
# # Protein expressed by each HSPC cell type that are miRNA targets
# p_targets <- sapply(HSPC_order, function(x) {
#   rownames(p_df[m_p_common,])[p_df[m_p_common, x] > 0]
# })
# p_targets_c <- sapply(HSPC_order, function(x) {
#   length(rownames(p_df[m_p_common,])[p_df[m_p_common, x] > 0])
# })
# # mRNA AND protein expressed by each cell type that are miRNA targets
# mp_targets <- sapply(HSPC_order, function(x) {
#   intersect(m_targets[x][[1]], p_targets[x][[1]])
# })
# mp_targets_c <- sapply(HSPC_order, function(x) {
#   length(intersect(m_targets[x][[1]], p_targets[x][[1]]))
# })
# # mRNA only for miRNA targets
# m_only <- sapply(HSPC_order, function(x) {
#   setdiff(m_targets[x][[1]], mp_targets[x][[1]])
# })
# m_only_c <- sapply(HSPC_order, function(x) {
#   length(setdiff(m_targets[x][[1]], mp_targets[x][[1]]))
# })
# # protein only for miRNA targets
# p_only <- sapply(HSPC_order, function(x) {
#   setdiff(p_targets[x][[1]], mp_targets[x][[1]])
# })
# p_only_c <- sapply(HSPC_order, function(x) {
#   length(setdiff(p_targets[x][[1]], mp_targets[x][[1]]))
# })
# # miRNA targets uniquely m_only in one cell type and mp_targets in other cells
# {
# miu_HSC <- length(intersect(setdiff(setdiff(setdiff(m_only$HSC, m_only$MPPa), m_only$MPPb), m_only$MPPc), 
#                             union(mp_targets$MPPa, union(mp_targets$MPPb, mp_targets$MPPc))))
# miu_HSC_names <- intersect(setdiff(setdiff(setdiff(m_only$HSC, m_only$MPPa), m_only$MPPb), m_only$MPPc), 
#                            union(mp_targets$MPPa, union(mp_targets$MPPb, mp_targets$MPPc)))
# miu_Ma <- length(intersect(setdiff(setdiff(setdiff(m_only$MPPa, m_only$HSC), m_only$MPPb), m_only$MPPc), 
#                            union(mp_targets$HSC, union(mp_targets$MPPb, mp_targets$MPPc))))
# miu_Ma_names <- intersect(setdiff(setdiff(setdiff(m_only$MPPa, m_only$HSC), m_only$MPPb), m_only$MPPc), 
#                           union(mp_targets$HSC, union(mp_targets$MPPb, mp_targets$MPPc)))
# miu_Mb <- length(intersect(setdiff(setdiff(setdiff(m_only$MPPb, m_only$MPPa), m_only$HSC), m_only$MPPc), 
#                            union(mp_targets$MPPa, union(mp_targets$HSC, mp_targets$MPPc))))
# miu_Mb_names <- intersect(setdiff(setdiff(setdiff(m_only$MPPb, m_only$MPPa), m_only$HSC), m_only$MPPc), 
#                           union(mp_targets$MPPa, union(mp_targets$HSC, mp_targets$MPPc)))
# miu_Mc <- length(intersect(setdiff(setdiff(setdiff(m_only$MPPc, m_only$MPPa), m_only$MPPb), m_only$HSC), 
#                            union(mp_targets$MPPa, union(mp_targets$MPPb, mp_targets$HSC))))
# miu_Mc_names <- intersect(setdiff(setdiff(setdiff(m_only$MPPc, m_only$MPPa), m_only$MPPb), m_only$HSC), 
#                           union(mp_targets$MPPa, union(mp_targets$MPPb, mp_targets$HSC)))
# }
# # miu_HSC expressed at least three times in other cell types
# {
#   raw_miu <- raw_df[miu_HSC_names, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   miu_HSC_names <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   miu_HSC_o <- miu_HSC
#   miu_HSC <- length(miu_HSC_names)
# }
# write.csv(maps$gene[which(maps$entrez %in% miu_HSC_names)], 
#           './Finalize/miuHSC.csv', row.names = FALSE, quote = FALSE)
# # miu_Ma expressed at least three times in other cell types
# {
#   raw_miu <- raw_df[miu_Ma_names, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   miu_Ma_names <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   miu_Ma_o <- miu_Ma
#   miu_Ma <- length(miu_Ma_names)
# }
# # miu_Mb expressed at least three times in other cell types
# {
#   raw_miu <- raw_df[miu_Mb_names, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   miu_Mb_names <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   miu_Mb_o <- miu_Mb
#   miu_Mb <- length(miu_Mb_names)
# }
# # miu_Mc expressed at least three times in other cell types
# {
#   raw_miu <- raw_df[miu_Mc_names, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   miu_Mc_names <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   miu_Mc_o <- miu_Mc
#   miu_Mc <- length(miu_Mc_names)
# }
# 
# # Figure 4A: Distribution of all mRNAs #
# figure20 <- function() {
#   # Expression of mRNA/protein for each cell type
#   m_w <- sapply(HSPC_order, function(x) {
#     rownames(r_df)[r_df[, x] > 0]
#   })
#   p_w <- sapply(HSPC_order, function(x) {
#     rownames(p_df)[p_df[, x] > 0]
#   })
#   
#   # mRNA and protein
#   mp_w <- sapply(HSPC_order, function(x) {
#     intersect(m_w[x][[1]], p_w[x][[1]])
#   })
#   mp_wc <- sapply(HSPC_order, function(x) {
#     length(intersect(m_w[x][[1]], p_w[x][[1]]))
#   })
#   
#   # mRNA only
#   mo_w <- sapply(HSPC_order, function(x) {
#     setdiff(m_w[x][[1]], p_w[x][[1]])
#   })
#   mo_wc <- sapply(HSPC_order, function(x) {
#     length(setdiff(m_w[x][[1]], p_w[x][[1]]))
#   })
#   
#   # protein only
#   po_w <- sapply(HSPC_order, function(x) {
#     setdiff(p_w[x][[1]], m_w[x][[1]])
#   })
#   po_wc <- sapply(HSPC_order, function(x) {
#     length(setdiff(p_w[x][[1]], m_w[x][[1]]))
#   })
#   
#   # unique to each cell type
#   hsc <- intersect(setdiff(setdiff(setdiff(mo_w$HSC, mo_w$MPPa), mo_w$MPPb), mo_w$MPPc),
#                    union(union(mp_w$MPPa, mp_w$MPPb), mp_w$MPPc))
#   mppa <- intersect(setdiff(setdiff(setdiff(mo_w$MPPa, mo_w$HSC), mo_w$MPPb), mo_w$MPPc),
#                     union(union(mp_w$HSC, mp_w$MPPb), mp_w$MPPc))
#   mppb <- intersect(setdiff(setdiff(setdiff(mo_w$MPPb, mo_w$MPPa), mo_w$HSC), mo_w$MPPc),
#                     union(union(mp_w$MPPa, mp_w$HSC), mp_w$MPPc))
#   mppc <- intersect(setdiff(setdiff(setdiff(mo_w$MPPc, mo_w$MPPa), mo_w$MPPb), mo_w$HSC),
#                     union(union(mp_w$MPPa, mp_w$MPPb), mp_w$HSC))
#   
#   # expressed at least three times in other cell types
#   # HSC
#   raw_miu <- raw_df[hsc, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   hsc <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   
#   for_file <- data.frame(gene = hsc)
#   for_file$in_miu <- for_file$gene %in% miu_HSC_names
#   for_file <- for_file[which(for_file$gene %in% maps$entrez),]
#   for_file$gene <- maps$gene[which(maps$entrez %in% for_file$gene)]
#   write.csv(for_file, './Finalize/unique_missing.csv', row.names = FALSE,
#             quote = FALSE)
#   # MPPa
#   raw_miu <- raw_df[mppa, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   mppa <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   # MPPb
#   raw_miu <- raw_df[mppb, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   mppb <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   # MPPc
#   raw_miu <- raw_df[mppc, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[HSPC_order, ]
#   length(sapply(raw_agg, max))
#   mppc <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   # Compiled
#   miu_group = c(HSC = length(hsc), MPPa = length(mppa),
#                 MPPb = length(mppb), MPPc = length(mppc))
#   
#   # Put into graphing form
#   by_count <- data.frame(list(counts = c(both = mp_wc, 
#                                                m = mo_wc - miu_group,
#                                                u = miu_group),
#                                     presence = c('both', 'both', 'both', 'both',
#                                                  'm', 'm', 'm', 'm',
#                                                  'u', 'u', 'u', 'u'),
#                                     cell = c('HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc')))
#   by_count$presence <- factor(by_count$presence, levels = c('m',
#                                                                   'u',
#                                                                   'both'))
#   
#   # Save info
#   write.csv(by_count, './Finalize/transcriptome_dist.csv')
#   
#   # Plot
#   ggplot(by_count, aes(x = cell, y = counts, fill = presence)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
#     theme_classic() +
#     ggtitle('All Distribution') +
#     annotate('text', x = HSPC_order, y = mp_wc + mo_wc + po_wc + 500,
#              label = sapply(miu_group, function(d) {
#                paste('U = ', d, sep = '')
#              }), size = 4.8) +
#     xlab('Cell') +
#     ylab('Targets Detected') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 17000)) +
#     scale_fill_manual(NULL, labels = c(both = 'Both',
#                                        m = 'mRNA Only',
#                                        u = 'mRNA Only (Unique)'),
#                       values=c(both = 'tomato',
#                                m = 'springgreen3',
#                                u = 'gold'))
# }
# figure20()
# 
# # Figure 4B, Supp 4A: mRNA Protein Log2 Plot #
# figure3 <- function() {
#   scatters <- lapply(c('HSC', 'MPPa', 'MPPb', 'MPPc'), function(x) {
#     plot_df <- data.frame(cbind(prot = p_common[, x], mRNA = r_common[, x]))
#     rownames(plot_df) <- rownames(p_common)
#     plot_df <- plot_df + 0.0001
#     plot_df <- log2(plot_df)
#     
#     ggplot(plot_df, aes(x = mRNA, y = prot)) +
#       geom_point(color = my_colors[x], size = 2.5) +
#       geom_point(color = 'black', shape = 1, size = 2.7) +
#       theme_classic() +
#       ggtitle(x) +
#       xlab('mRNA (log2)') +
#       ylab('Protein (log2)') +
#       theme(axis.text=element_text(size=14, colour = 'black'), 
#             axis.title=element_text(size=14),
#             plot.title = element_text(size = 14),
#             aspect.ratio = 1)
#   })
#   scatters
# }
# figure3()
# 
# # Figure 4C: Spearman Correlation #
# figure4 <- function() {
#   rhos <- sapply(HSPC_order, function(x){
#     bexp <- (p_common[, x] > 0) * (r_common[, x] > 0) == 1
#     column <- data.frame(cbind(p_common[bexp, x], r_common[bexp, x]))
#     cor.test( ~ X1 + X2,  data=column, method = "spearman",
#               continuity = FALSE, conf.level = 0.95)$estimate
#   })
#   names(rhos) <- sapply(names(rhos), function(x){sub('\\..+', '', x)})
#   rho_df <- round(data.frame(rhos), 3)
#   rho_df$cell <- rownames(rho_df)
#   ggplot(rho_df, aes(x = cell, y = rhos, fill = cell)) +
#     theme_classic()+
#     geom_bar(aes(fill = cell), color = '#000000', stat = 'identity', width = 0.5) + 
#     scale_fill_manual(NULL, values=c(my_colors['HSC'],
#                                      my_colors['MPPa'],
#                                      my_colors['MPPb'],
#                                      my_colors['MPPc'])) +
#     geom_text(aes(label=sapply(rhos, function(d) {
#       paste('rho = ', d, sep = '')
#     })), vjust = -0.5, size = 4.8) +
#     ggtitle('Protein vs mRNA Spearman Correlation') +
#     xlab('Cell Type') +
#     ylab('Correlation Score') + 
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.position = 'none') +
#     scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0))
# }
# figure4()
# 
# # Figure 4D: Pearson Correlation of Fold Changes #
# figure16 <- function() {
#   boxes <- sapply(HSPC_order, function(x){
#     both_ex <- ((p_all[, x] > 0) * (r_all[, x] > 0)) == 1
#     log2fold <- log2(p_all[both_ex, x] / r_all[both_ex, x])
#     return(log2fold)
#   })
#   boxes <- data.frame(boxes)
#   
#   pearson_matrix <- lapply(names(boxes), function(x) {
#     sapply(names(boxes), function(y) {
#       c(value = pearson(boxes[x][[1]], boxes[y][[1]]), X = x, Y = y)
#     })
#   })
#   pearson_matrix <- data.frame(t(do.call('cbind', pearson_matrix)))
#   pearson_matrix <- data.frame(pearson_matrix)
#   pearson_matrix$value <- as.double(as.character(pearson_matrix$value))
#   pearson_matrix$Y <- factor(pearson_matrix$Y, levels = c('HSC',
#                                                           'MPPa',
#                                                           'MPPb',
#                                                           'MPPc'))
#   pearson_matrix$X <- factor(pearson_matrix$X, levels = c('MPPc',
#                                                           'MPPb',
#                                                           'MPPa',
#                                                           'HSC'))
#   
#   ggplot(pearson_matrix, aes(x = X, y = Y)) +
#     geom_tile(aes(fill = value)) +
#     theme_classic() +
#     scale_fill_distiller(palette = "YlGnBu", direction = 1) +
#     ggtitle('Fold Change Pearson Correlation') +
#     theme(axis.line = element_blank(), axis.ticks = element_blank(),
#           axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_blank(),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.title =element_text(size = 14),
#           legend.text = element_text(size = 14)) +
#     labs(fill = 'Correlation Score')
# }
# figure16()
# 
# # Figure 4E, Supp 4C: Protein/mRNA log2 Fold Change Values for Underexpressed Detection #
# figure10 <- function() {
#   # Determine genes overexpressed in non-HSCs
#   overexp <- lapply(c('MPPa', 'MPPb', 'MPPc'), function(x) {
#     start_c <- 'HSC'
#     end_c <- x
#     rexp <- ((r_common[, start_c] > 0) * (r_common[, end_c] > 0)) == 1
#     r_exp_df <- log2(r_common[rexp,] + 0.00001)
#     p_exp_df <- log2(p_common[rexp,] + 0.00001)
#     p_zero_df <- p_exp_df == log2(0.00001)
#     p_zero_HM <- (p_zero_df[, start_c] + p_zero_df[, end_c]) > 0
#     p_nz_HM <- !p_zero_HM
#     
#     H_fc <- data.frame(p_exp_df[,start_c] - r_exp_df[,start_c])
#     rownames(H_fc) <- rownames(p_exp_df)
#     M_fc <- data.frame(p_exp_df[,end_c] - r_exp_df[,end_c])
#     rownames(M_fc) <- rownames(p_exp_df)
#     df = data.frame(H_fc,M_fc)
#     colnames(df) <- c('Start', 'End')
#     rownames(df) <- rownames(H_fc)
#     df$micro <- 'Normal'
#     
#     df_diff <- df$End - df$Start
#     upper <- percentile(df_diff[p_nz_HM], 97.5)
#     lower <- percentile(df_diff[p_nz_HM], 2.5)
#     
#     overexp <- ((df_diff > upper) * (p_nz_HM)) == 1
#     underexp <- ((df_diff <= lower) * (p_nz_HM)) == 1
#     return(rownames(df)[overexp])
#   })
#   overexp <- intersect(intersect(overexp[[1]], overexp[[2]]), overexp[[3]])
#   
#   # Narrow down to at least three in other cell types
#   raw_miu <- raw_df[overexp, ]
#   raw_df_t <- t(raw_miu)
#   raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#     sum(x > 0)
#   })
#   rownames(raw_agg) <- raw_agg$Group.1
#   raw_agg$Group.1 <- NULL
#   raw_agg <- raw_agg[c('MPPa', 'MPPb', 'MPPc'), ]
#   length(sapply(raw_agg, max))
#   overexp <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#   overexp_c <- length(overexp)
#   overexp <- data.frame(overexp)
#   colnames(overexp) <- c('Gene')
#   overexp$miRNA <- sapply(overexp$Gene, function(x) x %in% miRNAs)
#   overexp$Gene <- maps$gene[which(maps$entrez %in% overexp$Gene)]
#   write.csv(overexp, './Finalize/Only HSC Underexpress.csv', quote = FALSE,
#             row.names = FALSE)
#   
#   lapply(c('MPPa', 'MPPb', 'MPPc'), function(x) {
#     start_c <- 'HSC'
#     end_c <- x
#     rexp <- ((r_common[, start_c] > 0) * (r_common[, end_c] > 0)) == 1
#     r_exp_df <- log2(r_common[rexp,] + 0.00001)
#     p_exp_df <- log2(p_common[rexp,] + 0.00001)
#     p_zero_df <- p_exp_df == log2(0.00001)
#     p_zero_HM <- (p_zero_df[, start_c] + p_zero_df[, end_c]) > 0
#     p_nz_HM <- !p_zero_HM
#     
#     H_fc <- data.frame(p_exp_df[,start_c] - r_exp_df[,start_c])
#     rownames(H_fc) <- rownames(p_exp_df)
#     M_fc <- data.frame(p_exp_df[,end_c] - r_exp_df[,end_c])
#     rownames(M_fc) <- rownames(p_exp_df)
#     df = data.frame(H_fc,M_fc)
#     colnames(df) <- c('Start', 'End')
#     rownames(df) <- rownames(H_fc)
#     df$micro <- 'Normal'
#     
#     df_diff <- df$End - df$Start
#     diff_mu <- mean(df_diff[p_nz_HM])
#     diff_sig <- sd(df_diff[p_nz_HM])
#     upper <- percentile(df_diff[p_nz_HM], 97.5)
#     lower <- percentile(df_diff[p_nz_HM], 2.5)
#     
#     overexp <- ((df_diff > upper) * (p_nz_HM)) == 1
#     underexp <- ((df_diff <= lower) * (p_nz_HM)) == 1
#     
#     of_choice <- intersect(gene_set$TTTGCAC_MIR19A_MIR19B, rownames(df))
#     zeros <- rownames(df)[p_zero_HM]
#     
#     df[overexp, 'micro'] <- 'Overexp'
#     # df[of_choice, 'micro'] <- 'Micro'
#     df$micro <- factor(df$micro, levels = c(#'Micro',
#       'Overexp',
#       'Normal'))
#     df <- df[p_nz_HM,]
#     ggplot(df, aes(x = Start, y = End, colour = micro)) +
#       xlab(paste(start_c, 'Log2 Fold Change')) +
#       ylab(paste(end_c, 'Log2 Fold Change')) +
#       geom_point(data = subset(df, micro == 'Normal'),
#                  aes(x = Start, y = End, color = 'Normal')) +
#       geom_point(data = subset(df, micro == 'Overexp'),
#                  aes(x = Start, y = End, color = 'Overexp')) +
#       # geom_point(data = subset(df, micro == 'Micro'),
#       #            aes(x = Start, y = End, color = 'Micro')) +
#       theme_classic() +
#       # geom_segment(aes(x = -8, y = -8, xend = 10, yend = 10), colour = 'blue') + 
#       # geom_segment(aes(x = -8, y = -8 + upper, xend = 10, 
#       #                  yend = 10 + upper), colour = 'red') + 
#       # geom_segment(aes(x = -8, y = -8 + lower, xend = 10, 
#       #                  yend = 10 + lower), colour = 'red') + 
#       scale_color_manual(NULL, labels = c(Normal = "Bottom 97.5%", 
#                                           # Micro = 'Mir19A & Mir19B Targets',
#                                           Overexp = 'Top 2.5%'),
#                          values=c(Normal = '#d3d3d3',
#                                   # Micro = '#ff6347',
#                                   Overexp = 'gold')) +
#       theme(axis.text=element_text(size=14, colour = 'black'), 
#             axis.title=element_text(size=14),
#             plot.title = element_text(size = 14),
#             aspect.ratio = 1, 
#             legend.title = element_blank(), legend.text = element_text(size = 14)) +
#       ggtitle('Protein/mRNA log2 Fold Change Values')
#   })
# }
# figure10()
# 
# # Figure 5C: miRNA Target Distribution all throughout #
# figure5 <- function() {
#   miu_group = c(HSC = miu_HSC, MPPa = miu_Ma,
#                 MPPb = miu_Mb, MPPc = miu_Mc)
#   
#   to_write <- data.frame(list(counts = c(both = mp_targets_c, 
#                                          p =p_only_c, 
#                                          m = m_only_c - miu_group,
#                                          u = miu_group),
#                               presence = c('both', 'both', 'both', 'both',
#                                            'p', 'p', 'p', 'p',
#                                            'm', 'm', 'm', 'm',
#                                            'u', 'u', 'u', 'u'),
#                               cell = c('HSC', 'MPPa', 'MPPb', 'MPPc',
#                                        'HSC', 'MPPa', 'MPPb', 'MPPc',
#                                        'HSC', 'MPPa', 'MPPb', 'MPPc',
#                                        'HSC', 'MPPa', 'MPPb', 'MPPc')))
#   
#   to_write$presence <- factor(to_write$presence, levels = c('m',
#                                                             'p',
#                                                             'u',
#                                                             'both'))
#   
#   miRNA_by_count <- data.frame(list(counts = c(both = mp_targets_c, 
#                                                m = m_only_c - miu_group,
#                                                u = miu_group),
#                                     presence = c('both', 'both', 'both', 'both',
#                                                  'm', 'm', 'm', 'm',
#                                                  'u', 'u', 'u', 'u'),
#                                     cell = c('HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc')))
#   miRNA_by_count$presence <- factor(miRNA_by_count$presence, levels = c('m',
#                                                                         'u',
#                                                                         'both'))
#   write.csv(to_write, './Finalize/miRNA_by_count.csv', 
#             quote = FALSE, row.names = TRUE)
#   ggplot(miRNA_by_count, aes(x = cell, y = counts, fill = presence)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
#     theme_classic() +
#     ggtitle('miRNA Target Distribution') +
#     annotate('text', x = HSPC_order, y = mp_targets_c + 
#                p_only_c + m_only_c + 450,
#              label = sapply(miu_group, function(d) {
#                paste('U = ', d, sep = '')
#              }), size = 4.8) +
#     xlab('Cell') +
#     ylab('Targets Detected') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 15000)) +
#     scale_fill_manual(NULL, labels = c(both = 'Both',
#                                        m = 'mRNA Only',
#                                        u = 'mRNA Only (Unique)'),
#                       values=c(both = 'tomato',
#                                m = 'springgreen3',
#                                u = 'gold'))
# }
# figure5()
# 
# # Figure 5D: Count of gene overlap with miRNA gene set #
# figure8 <- function() {
#   miu_HSC_overlap_n <- sapply(gene_set, function(x) {
#     length(intersect(miu_HSC_names, x))
#   })
#   miu_HSC_overlap <- sapply(gene_set, function(x) {
#     intersect(miu_HSC_names, x)
#   })
#   miu_HSC_overlap_prop <- sapply(gene_set, function(x) {
#     length(intersect(miu_HSC_names, x)) / length(x)
#   })
#   mp_exp <- rownames(p_common)[
#     as.logical((p_common['HSC'] > 0) * (r_common['HSC'] > 0))
#     ]
#   all_HSC_overlap_prop <- sapply(gene_set, function(x) {
#     length(intersect(mp_exp, x)) / length(x)
#   })
#   mh_df <- data.frame(cbind(miu_HSC_overlap_n, miu_HSC_overlap_prop,
#                             all_HSC_overlap_prop))
#   mh_df <- mh_df[order(mh_df$miu_HSC_overlap_n),]
#   colnames(mh_df) <- c('count', 'prop', 'all.prop')
#   mh_overlap_df <- data.frame(sort(miu_HSC_overlap_n, decreasing = FALSE))
#   colnames(mh_overlap_df) <- 'count'
#   write.csv(mh_df, './Finalize/miRNA_target_overlap.csv',
#             quote = FALSE)
#   mh_overlap_df$entrez <- rownames(mh_overlap_df)
#   mh_overlap_df$entrez <- factor(mh_overlap_df$entrez, levels = mh_overlap_df$entrez)
#   
#   mh_overlap_df <- mh_overlap_df[mh_overlap_df$count > 0, ]
#   
#   p3 <- ggplot(mh_overlap_df, aes(x = entrez, y = count)) +
#     geom_bar(stat = 'identity', width = 1, fill = 'springgreen3', color = 'black') +
#     theme_classic() +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, axis.text.x = element_blank(),
#           axis.ticks.x = element_blank()) +
#     ggtitle('HSC mRNA only overlap with miRNA target lists') +
#     xlab('miRNA Target Gene Set') +
#     ylab('Count') +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 250))
#   return(list(p3))
# }
# figure8()
# 
# # Figure 5E: Expression of miu_HSC proteins in aged HSCs #
# figure22 <- function() {
#   aexp <- rownames(p_df)[p_df$AgedHSC > 0]
#   df <- data.frame(c(length(intersect(aexp, miu_HSC_names)), 
#                      length(miu_HSC_names) - length(intersect(aexp, miu_HSC_names))),
#                    c('expressed', 'not'))
#   colnames(df) <- c('count', 'type')
#   write.table(df, './Finalize/agedExpressed.csv', row.names = TRUE, quote = FALSE,
#               sep = ',')
#   ggplot(df, aes(x = '', y = count, fill = count)) +
#     geom_bar(stat = 'identity', color = 'white', size = 1.5) +
#     coord_polar('y', start = 0) +
#     theme_classic() +
#     ggtitle('Percent unique expressed in aged HSCs') +
#     theme(axis.text=element_blank(),
#           axis.title=element_blank(),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.line = element_blank(), legend.title = element_blank(),
#           legend.position = 'none') +
#     scale_fill_gradient2(low = 'white', high = 'steelblue') +
#     guides(fill = guide_legend(reverse = TRUE))
# }
# figure22()
# 
# # Supp 4B: PCA with genes as features #
# figure18 <- function() {
#   boxes <- sapply(HSPC_order, function(x){
#     both_ex <- ((p_all[, x] > 0) * (r_all[, x] > 0)) == 1
#     log2fold <- log2(p_all[both_ex, x] / r_all[both_ex, x])
#     return(log2fold)
#   })
#   boxes <- data.frame(boxes)
#   
#   pca <- prcomp(boxes)
#   
#   autoplot(pca, loadings = TRUE) +
#     theme_classic() +
#     ggtitle('PCA by Protein/mRNA Fold Change') +
#     theme(axis.text=element_text(size=14, colour = 'black'),
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14))
# }
# figure18()
# 
# # Supp 4D: Null Distribution of all mRNAs #
# figure21 <- function() {
#   m_w <- sapply(HSPC_order, function(x) {
#     rownames(r_df)[r_df[, x] > 0]
#   })
#   p_w <- sapply(HSPC_order, function(x) {
#     rownames(p_df)[p_df[, x] > 0]
#   })
#   
#   mp_w <- sapply(HSPC_order, function(x) {
#     intersect(m_w[x][[1]], p_w[x][[1]])
#   })
#   mp_wc <- sapply(HSPC_order, function(x) {
#     length(intersect(m_w[x][[1]], p_w[x][[1]]))
#   })
#   
#   mo_w <- sapply(HSPC_order, function(x) {
#     setdiff(m_w[x][[1]], p_w[x][[1]])
#   })
#   mo_wc <- sapply(HSPC_order, function(x) {
#     length(setdiff(m_w[x][[1]], p_w[x][[1]]))
#   })
#   
#   po_w <- sapply(HSPC_order, function(x) {
#     setdiff(p_w[x][[1]], m_w[x][[1]])
#   })
#   po_wc <- sapply(HSPC_order, function(x) {
#     length(setdiff(p_w[x][[1]], m_w[x][[1]]))
#   })
#   
#   HSC <- data.frame(m_w['HSC'])
#   colnames(HSC) <- 'entrez'
#   HSC$cate <- c(rep('mp', mp_wc['HSC']), rep('m', mo_wc['HSC']))
#   Ma <- data.frame(m_w['MPPa'])
#   colnames(Ma) <- 'entrez'
#   Ma$cate <- c(rep('mp', mp_wc['MPPa']), rep('m', mo_wc['MPPa']))
#   Mb <- data.frame(m_w['MPPb'])
#   colnames(Mb) <- 'entrez'
#   Mb$cate <- c(rep('mp', mp_wc['MPPb']), rep('m', mo_wc['MPPb']))
#   Mc <- data.frame(m_w['MPPc'])
#   colnames(Mc) <- 'entrez'
#   Mc$cate <- c(rep('mp', mp_wc['MPPc']), rep('m', mo_wc['MPPc']))
#   
#   randomized <- lapply(c(1:1000), function(z) {
#     set.seed((z-1) * 4 + 1)
#     HSC$cate <- sample(HSC$cate)
#     set.seed((z-1) * 4 + 2)
#     Ma$cate <- sample(Ma$cate)
#     set.seed((z-1) * 4 + 3)
#     Mb$cate <- sample(Mb$cate)
#     set.seed((z-1) * 4 + 4)
#     Mc$cate <- sample(Mc$cate)
#     
#     comb <- list(HSC = HSC,Ma = Ma, Mb = Mb, Mc = Mc)
#     
#     m_t <- sapply(names(comb), function(x) {
#       comb[x][[1]][comb[x][[1]][, 'cate'] == 'm', 'entrez']
#     })
#     m_t_c <- sapply(names(comb), function(x) {
#       length(comb[x][[1]][comb[x][[1]][, 'cate'] == 'm', 'entrez'])
#     })
#     mp_t <- sapply(names(comb), function(x) {
#       comb[x][[1]][comb[x][[1]][, 'cate'] == 'mp', 'entrez']
#     })
#     mp_t_c <- sapply(names(comb), function(x) {
#       length(comb[x][[1]][comb[x][[1]][, 'cate'] == 'mp', 'entrez'])
#     })
#     h_m_t <- intersect(setdiff(setdiff(setdiff(m_t$HSC, m_t$Ma), m_t$Mb), m_t$Mc),
#                        union(union(mp_t$Ma, mp_t$Mb), mp_t$Mc))
#     h_m_t_c <- length(h_m_t)
#     ma_m_t <- intersect(setdiff(setdiff(setdiff(m_t$Ma, m_t$HSC), m_t$Mb), m_t$Mc),
#                         union(union(mp_t$HSC, mp_t$Mb), mp_t$Mc))
#     ma_m_t_c <- length(ma_m_t)
#     mb_m_t <- intersect(setdiff(setdiff(setdiff(m_t$Mb, m_t$HSC), m_t$Ma), m_t$Mc),
#                         union(union(mp_t$HSC, mp_t$Ma), mp_t$Mc))
#     mb_m_t_c <- length(mb_m_t)
#     mc_m_t <- intersect(setdiff(setdiff(setdiff(m_t$Mc, m_t$HSC), m_t$Ma), m_t$Mb),
#                         union(union(mp_t$HSC, mp_t$Ma), mp_t$Mb))
#     mc_m_t_c <- length(mc_m_t)
#     
#     raw_miu <- raw_df[h_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     h_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     h_m_t_c <- length(h_m_t)
#     
#     raw_miu <- raw_df[ma_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     ma_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     ma_m_t_c <- length(ma_m_t)
#     
#     raw_miu <- raw_df[mb_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     mb_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     mb_m_t_c <- length(mb_m_t)
#     
#     raw_miu <- raw_df[mc_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     mc_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     mc_m_t_c <- length(mc_m_t)
#     
#     m_t_group <- c(HSC = h_m_t_c, Ma = ma_m_t_c, Mb = mb_m_t_c, Mc = mc_m_t_c)
#     return(m_t_group)
#   })
#   rand_df <- data.frame(do.call(rbind, randomized))
#   m_t_group <- sapply(rand_df, mean)
#   m_t_group
#   m_t_sd <- sapply(rand_df, sd)
#   lower <- round(m_t_group - 2 * m_t_sd)
#   upper <- round(m_t_group + 2 * m_t_sd)
#   
#   miu_group = c(HSC = miu_HSC, MPPa = miu_Ma,
#                 MPPb = miu_Mb, MPPc = miu_Mc)
#   
#   miRNA_by_count <- data.frame(list(counts = c(both = mp_wc, 
#                                                m = mo_wc - m_t_group,
#                                                u = m_t_group),
#                                     presence = c('both', 'both', 'both', 'both',
#                                                  'm', 'm', 'm', 'm',
#                                                  'u', 'u', 'u', 'u'),
#                                     cell = c('HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc')))
#   miRNA_by_count$presence <- factor(miRNA_by_count$presence, levels = c('m',
#                                                                         'u',
#                                                                         'both'))
#   ggplot(miRNA_by_count, aes(x = cell, y = counts, fill = presence)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
#     theme_classic() +
#     annotate('text', x = HSPC_order, y = mp_wc + mo_wc + po_wc + 500,
#              label = c(paste('U = ', '[',lower['HSC'],', ', upper['HSC'],']', sep = ''),
#                        paste('U = ', '[',lower['Ma'],', ', upper['Ma'],']', sep = ''), 
#                        paste('U = ', '[',lower['Mb'],', ', upper['Mb'],']', sep = ''), 
#                        paste('U = ', '[',lower['Mc'],', ', upper['Mc'],']', sep = '')),
#              size = 3.5) +
#     ggtitle('All NULL Distribution') +
#     xlab('Cell') +
#     ylab('Targets Detected') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 17000)) +
#     scale_fill_manual(NULL, labels = c(both = 'Both',
#                                        m = 'mRNA Only',
#                                        u = 'mRNA Only (Unique)'),
#                       values=c(both = 'lightpink2',
#                                m = 'darkseagreen2',
#                                u = 'peachpuff2'))
# }
# figure21()
# 
# # Supp 5A: Ribosomal protein pie chart #
# figure15 <- function() {
#   ribo_genes <- rownames(comb_df)[sapply(rownames(comb_df), function(x) grepl('Rp[sl][0-9]+[^k]', x))]
#   min_dist <- sapply(colnames(comb_df), function(x) {
#     COI = x
#     ribo_zero_HSC <- ribo_genes[apply(comb_df[ribo_genes,], 1, function(x) {
#       if(sum(x == 0) == 1) {
#         if(x[colnames(comb_df) == COI] == 0) {
#           return(TRUE)
#         } else {
#           return(FALSE)
#         }
#       } else {
#         FALSE
#       }
#     })]
#     ribo_min_HSC <- ribo_genes[apply(comb_df[ribo_genes,], 1, function(x) {
#       if(sum(x == 0) >= 1) {
#         return(FALSE)
#       } else {
#         if(min(x) == x[colnames(comb_df) == COI]) {
#           if(sum(x == min(x)) == 1) {
#             return(TRUE)
#           } else{
#             return(FALSE)
#           }
#         } else{
#           return(FALSE)
#         }
#       }
#     })]
#     length(ribo_min_HSC) + length(ribo_zero_HSC)
#   })
#   min_dist <- data.frame(min_dist)
#   colnames(min_dist) <- 'value'
#   min_dist$cell <- rownames(min_dist)
#   min_dist <- rbind(min_dist, other = c(length(ribo_genes) - sum(min_dist$value), 'other'))
#   min_dist$value <- as.integer(min_dist$value)
#   min_dist <- min_dist[order(as.integer(min_dist$value), decreasing = FALSE), ]
#   min_dist$cell <- factor(min_dist$cell, levels = min_dist$cell)
#   write.csv(min_dist, './Finalize/Ribosomal Numbers.csv', 
#             quote = FALSE, row.names = TRUE)
#   ggplot(min_dist, aes(x = '', y = as.integer(value), fill = as.integer(value))) +
#     geom_bar(stat = 'identity', color = 'white', size = 1.5) +
#     coord_polar('y', start = 0) +
#     theme_classic() +
#     ggtitle('Unique minimum/zero values for ribosomal proteins') +
#     theme(axis.text=element_blank(),
#           axis.title=element_blank(),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.line = element_blank(), legend.title = element_blank(),
#           legend.position = 'none') +
#     scale_fill_gradient2(low = 'white', high = 'steelblue') +
#     guides(fill = guide_legend(reverse = TRUE))
# }
# figure15()
# 
# # Supp 5B: Null Distribution of Unique miRNA Target Enrichment #
# figure6 <- function() {
#   HSC <- data.frame(m_targets$HSC)
#   colnames(HSC) <- 'entrez'
#   HSC$cate <- c(rep('mp', mp_targets_c['HSC']), rep('m', m_only_c['HSC']))
#   Ma <- data.frame(m_targets$MPPa)
#   colnames(Ma) <- 'entrez'
#   Ma$cate <- c(rep('mp', mp_targets_c['MPPa']), rep('m', m_only_c['MPPa']))
#   Mb <- data.frame(m_targets$MPPb)
#   colnames(Mb) <- 'entrez'
#   Mb$cate <- c(rep('mp', mp_targets_c['MPPb']), rep('m', m_only_c['MPPb']))
#   Mc <- data.frame(m_targets$MPPc)
#   colnames(Mc) <- 'entrez'
#   Mc$cate <- c(rep('mp', mp_targets_c['MPPc']), rep('m', m_only_c['MPPc']))
#   
#   randomized <- lapply(c(1:1000), function(z) {
#     set.seed((z-1) * 4 + 1)
#     HSC$cate <- sample(HSC$cate)
#     set.seed((z-1) * 4 + 2)
#     Ma$cate <- sample(Ma$cate)
#     set.seed((z-1) * 4 + 3)
#     Mb$cate <- sample(Mb$cate)
#     set.seed((z-1) * 4 + 4)
#     Mc$cate <- sample(Mc$cate)
#     
#     comb <- list(HSC = HSC,Ma = Ma, Mb = Mb, Mc = Mc)
#     
#     m_t <- sapply(names(comb), function(x) {
#       comb[x][[1]][comb[x][[1]][, 'cate'] == 'm', 'entrez']
#     })
#     m_t_c <- sapply(names(comb), function(x) {
#       length(comb[x][[1]][comb[x][[1]][, 'cate'] == 'm', 'entrez'])
#     })
#     mp_t <- sapply(names(comb), function(x) {
#       comb[x][[1]][comb[x][[1]][, 'cate'] == 'mp', 'entrez']
#     })
#     mp_t_c <- sapply(names(comb), function(x) {
#       length(comb[x][[1]][comb[x][[1]][, 'cate'] == 'mp', 'entrez'])
#     })
#     h_m_t <- intersect(setdiff(setdiff(setdiff(m_t$HSC, m_t$Ma), m_t$Mb), m_t$Mc),
#                        union(union(mp_t$Ma, mp_t$Mb), mp_t$Mc))
#     h_m_t_c <- length(h_m_t)
#     ma_m_t <- intersect(setdiff(setdiff(setdiff(m_t$Ma, m_t$HSC), m_t$Mb), m_t$Mc),
#                         union(union(mp_t$HSC, mp_t$Mb), mp_t$Mc))
#     ma_m_t_c <- length(ma_m_t)
#     mb_m_t <- intersect(setdiff(setdiff(setdiff(m_t$Mb, m_t$HSC), m_t$Ma), m_t$Mc),
#                         union(union(mp_t$HSC, mp_t$Ma), mp_t$Mc))
#     mb_m_t_c <- length(mb_m_t)
#     mc_m_t <- intersect(setdiff(setdiff(setdiff(m_t$Mc, m_t$HSC), m_t$Ma), m_t$Mb),
#                         union(union(mp_t$HSC, mp_t$Ma), mp_t$Mb))
#     mc_m_t_c <- length(mc_m_t)
#     
#     # h_m_t_c <- h_m_t_c * miu_HSC / miu_HSC_o
#     # ma_m_t_c <- ma_m_t_c * miu_Ma / miu_Ma_o
#     # mb_m_t_c <- mb_m_t_c * miu_Mb / miu_Mb_o
#     # mc_m_t_c <- mc_m_t_c * miu_Mc / miu_Mc_o
#     
#     raw_miu <- raw_df[h_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     h_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     h_m_t_c <- length(h_m_t)
#     
#     raw_miu <- raw_df[ma_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     ma_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     ma_m_t_c <- length(ma_m_t)
#     
#     raw_miu <- raw_df[mb_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     mb_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     mb_m_t_c <- length(mb_m_t)
#     
#     raw_miu <- raw_df[mc_m_t, ]
#     raw_df_t <- t(raw_miu)
#     raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#       sum(x > 0)
#     })
#     rownames(raw_agg) <- raw_agg$Group.1
#     raw_agg$Group.1 <- NULL
#     raw_agg <- raw_agg[HSPC_order, ]
#     length(sapply(raw_agg, max))
#     mc_m_t <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
#     mc_m_t_c <- length(mc_m_t)
#     
#     m_t_group <- c(HSC = h_m_t_c, Ma = ma_m_t_c, Mb = mb_m_t_c, Mc = mc_m_t_c)
#     return(m_t_group)
#   })
#   rand_df <- data.frame(do.call(rbind, randomized))
#   m_t_group <- sapply(rand_df, mean)
#   m_t_group
#   m_t_sd <- sapply(rand_df, sd)
#   lower <- round(m_t_group - 2 * m_t_sd)
#   upper <- round(m_t_group + 2 * m_t_sd)
#   
#   miu_group = c(HSC = miu_HSC, MPPa = miu_Ma,
#                 MPPb = miu_Mb, MPPc = miu_Mc)
#   
#   miRNA_by_count <- data.frame(list(counts = c(both = mp_targets_c, 
#                                                m = m_only_c - m_t_group,
#                                                u = m_t_group),
#                                     presence = c('both', 'both', 'both', 'both',
#                                                  'm', 'm', 'm', 'm',
#                                                  'u', 'u', 'u', 'u'),
#                                     cell = c('HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc',
#                                              'HSC', 'MPPa', 'MPPb', 'MPPc')))
#   miRNA_by_count$presence <- factor(miRNA_by_count$presence, levels = c('m',
#                                                                         'u',
#                                                                         'both'))
#   ggplot(miRNA_by_count, aes(x = cell, y = counts, fill = presence)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
#     theme_classic() +
#     annotate('text', x = HSPC_order, y = m_only_c + mp_targets_c + 450,
#              label = c(paste('U = ', '[',lower['HSC'],', ', upper['HSC'],']', sep = ''),
#                        paste('U = ', '[',lower['Ma'],', ', upper['Ma'],']', sep = ''), 
#                        paste('U = ', '[',lower['Mb'],', ', upper['Mb'],']', sep = ''), 
#                        paste('U = ', '[',lower['Mc'],', ', upper['Mc'],']', sep = '')),
#              size = 3.5) +
#     ggtitle('miRNA Target NULL Distribution') +
#     xlab('Cell') +
#     ylab('Targets Detected') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 15000)) +
#     scale_fill_manual(NULL, labels = c(both = 'Both',
#                                        m = 'mRNA Only',
#                                        u = 'mRNA Only (Unique)'),
#                       values=c(both = 'lightpink2',
#                                m = 'darkseagreen2',
#                                u = 'peachpuff2'))
# }
# figure6()
# 
# # Supp 5C: RNA expression with miuHSC overlayed #
# figure17 <- function() {
#   plots <- lapply(c('MPPa', 'MPPb', 'MPPc'), function(x) {
#     rcom <- (r_df[, 'HSC'] > 0) * (r_df[, x] > 0)
#     rcom <- rcom == 1
#     rcom <- log2(r_df[rcom, c('HSC', x)])
#     colnames(rcom) <- c('HSC', 'MPP')
#     rcom$type <- 'Normal'
#     rcom[miu_HSC_names, 'type'] <- 'miRNA'
#     ggplot(rcom, aes(x = HSC, y = MPP, color = tye)) +
#       theme_classic() +
#       geom_point(data = subset(rcom, type == 'Normal'),
#                  aes(x = HSC, y = MPP, color = 'Normal')) +
#       geom_point(data = subset(rcom, type == 'miRNA'),
#                  aes(x = HSC, y = MPP, color = 'miRNA')) +
#       theme(axis.text=element_text(size=14, colour = 'black'), 
#             axis.title=element_text(size=14),
#             plot.title = element_text(size = 14),
#             aspect.ratio = 1, legend.text=element_text(size=14),
#             axis.text.x = element_text(angle = 45, hjust = 1)) +
#       ggtitle('Unique HSC miRNA targets') +
#       xlab('HSC mRNA (log2)') +
#       ylab(paste(x, 'mRNA (log2)')) +
#       scale_color_manual(NULL, labels = c(miRNA = 'miRNA', Normal = 'Others'),
#                          values = c(Normal = '#c3c3c3', miRNA = 'tomato'))
#   })
#   plots
# }
# figure17()
# 
# # Supp 5D: Expression of mRNA/protein in aged HSCs in relation to unique proteins #
# figure23 <- function() {
#   miu_HSC_overlap <- sapply(gene_set, function(x) {
#     length(intersect(miu_HSC_names, x))
#   })
#   miu_HSC_overlap_n <- sapply(gene_set, function(x) {
#     intersect(miu_HSC_names, x)
#   })
#   age_exp <- rownames(p_df)[p_df$AgedHSC > 0]
#   my_aged_overlap <- sapply(gene_set, function(x) {
#     100 * length(intersect(intersect(miu_HSC_names, x), age_exp)) / length(intersect(miu_HSC_names, x))
#   })
#   my_aged_overlap_n <- sapply(gene_set, function(x) {
#     intersect(intersect(miu_HSC_names, x), age_exp)
#   })
#   mh_df <- data.frame(cbind(miu_HSC_overlap, my_aged_overlap))
#   colnames(mh_df) <- c('count', 'percent')
#   mh_df_no_nan <- mh_df[!is.nan(mh_df$percent),]
#   ggplot(mh_df_no_nan, aes(x = count, y = percent)) +
#     geom_point(color = 'springgreen3')+
#     theme_classic() +
#     ggtitle('% untranslated miRNA targets expressed in agedHSC') +
#     xlab('# miRNA target genes uniquely mRNA only in HSCs') +
#     ylab('% detected in agedHSC') + 
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = 'none')
# }
# figure23()
# 
# #################################################
# 
# # 
# comb_df <- read.csv('3_gene_intensity_comb.txt', header = TRUE, row.names = 1)
# comb_df <- comb_df[,order(match(colnames(comb_df), pref_order))]
# exp_df <- comb_df > 0
# comb_cells <- colnames(comb_df)
# 
# # HSC and Ma
# miu_HMa <- length(intersect(setdiff(setdiff(intersect(m_only$HSC, m_only$MPPa), m_only$MPPb), m_only$MPPc), 
#                             union(mp_targets$MPPb, mp_targets$MPPc)))
# miu_HMa_names <- intersect(setdiff(setdiff(intersect(m_only$HSC, m_only$MPPa), m_only$MPPb), m_only$MPPc), 
#                            union(mp_targets$MPPb, mp_targets$MPPc))
# 
# miu_HMb <- length(intersect(setdiff(intersect(setdiff(m_only$HSC, m_only$MPPa), m_only$MPPb), m_only$MPPc), 
#                             union(mp_targets$MPPa, mp_targets$MPPc)))
# miu_HMb_names <- intersect(setdiff(intersect(setdiff(m_only$HSC, m_only$MPPa), m_only$MPPb), m_only$MPPc), 
#                            union(mp_targets$MPPa, mp_targets$MPPc))
# 
# raw_miu <- raw_df[miu_HMa_names, ]
# raw_df_t <- t(raw_miu)
# raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#   sum(x > 0)
# })
# rownames(raw_agg) <- raw_agg$Group.1
# raw_agg$Group.1 <- NULL
# raw_agg <- raw_agg[HSPC_order, ]
# length(sapply(raw_agg, max))
# miu_HMa_names <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
# miu_HMa_o <- miu_HMa
# miu_HMa <- length(miu_HMa_names)
# 
# raw_miu <- raw_df[miu_HMb_names, ]
# raw_df_t <- t(raw_miu)
# raw_agg <- aggregate(raw_df_t, list(raw_order), function(x){
#   sum(x > 0)
# })
# rownames(raw_agg) <- raw_agg$Group.1
# raw_agg$Group.1 <- NULL
# raw_agg <- raw_agg[HSPC_order, ]
# length(sapply(raw_agg, max))
# miu_HMb_names <- colnames(raw_agg)[sapply(raw_agg, max) >= 3]
# miu_HMb_o <- miu_HMb
# miu_HMb <- length(miu_HMb_names)
# 
# 
# 
# 
# 
# 
# # Figure 1 : Unique proteins within each cell type #
# figure1 <- function() {
#   by_gene <- rowSums(exp_df)
#   names(by_gene) <- rownames(exp_df)
#   uniques = sapply(comb_cells, function(x) {
#     singly <- by_gene == 1
#     of_cell <- ((singly * exp_df[, x]) == 1)
#     return(names(by_gene)[of_cell])
#   })
#   uniques_l <- sapply(pref_order, function(y){
#     length(uniques[y][[1]])
#   })
#   sapply(pref_order, function(y) {
#     write.table(uniques[y][[1]], paste('./Finalize/', y, '_Unique_Expression.txt', 
#                                        sep = ''), quote = FALSE, row.names = FALSE)
#   })
#   a_uniques = names(by_gene)[(((exp_df[, 'HSC'] + exp_df[, 'MPPa'] +
#                                   exp_df[, 'MPPb'] + exp_df[, 'MPPc'] + 
#                                   exp_df[, 'CLP'] + exp_df[, 'CMP'] +
#                                   exp_df[, 'GMP'] + exp_df[, 'MEP']) == 0) 
#                               * ((exp_df[, 'AgedHSC'] + exp_df[, 'AgedMPPa'] +
#                                     exp_df[, 'AgedMPPb'] + exp_df[, 'AgedMPPc']) > 0)) == 1]
#   a_uniques_l <- sapply(c('AgedHSC', 'AgedMPPa', 'AgedMPPb', 'AgedMPPc'), function(y) {
#     length(intersect(rownames(exp_df)[exp_df[, y]], a_uniques))
#   })
#   sapply(c('AgedHSC', 'AgedMPPa', 'AgedMPPb', 'AgedMPPc'), function(y) {
#     write.table(a_uniques[y][[1]], paste('./Finalize/', y, '_Group_Expression.txt', 
#                                        sep = ''), quote = FALSE, row.names = FALSE)
#   })
#   h_uniques = names(by_gene)[(((exp_df[, 'AgedHSC'] + exp_df[, 'AgedMPPa'] +
#                                   exp_df[, 'AgedMPPb'] + exp_df[, 'AgedMPPc'] + 
#                                   exp_df[, 'CLP'] + exp_df[, 'CMP'] +
#                                   exp_df[, 'GMP'] + exp_df[, 'MEP']) == 0) 
#                               * ((exp_df[, 'HSC'] + exp_df[, 'MPPa'] +
#                                     exp_df[, 'MPPb'] + exp_df[, 'MPPc']) > 0)) == 1]
#   h_uniques_l <- sapply(c('HSC', 'MPPa', 'MPPb', 'MPPc'), function(y) {
#     length(intersect(rownames(exp_df)[exp_df[, y]], h_uniques))
#   })
#   sapply(c('HSC', 'MPPa', 'MPPb', 'MPPc'), function(y) {
#     write.table(h_uniques[y][[1]], paste('./Finalize/', y, '_Group_Expression.txt', 
#                                          sep = ''), quote = FALSE, row.names = FALSE)
#   })
#   o_uniques = names(by_gene)[(((exp_df[, 'AgedHSC'] + exp_df[, 'AgedMPPa'] +
#                                   exp_df[, 'AgedMPPb'] + exp_df[, 'AgedMPPc'] + 
#                                   exp_df[, 'HSC'] + exp_df[, 'MPPa'] +
#                                   exp_df[, 'MPPb'] + exp_df[, 'MPPc']) == 0) 
#                               * ((exp_df[, 'CLP'] + exp_df[, 'CMP'] +
#                                     exp_df[, 'GMP'] + exp_df[, 'MEP']) > 0)) == 1]
#   o_uniques_l <- sapply(c('CLP', 'CMP', 'GMP', 'MEP'), function(y) {
#     length(intersect(rownames(exp_df)[exp_df[, y]], o_uniques))
#   })
#   sapply(c('CLP', 'CMP', 'GMP', 'MEP'), function(y) {
#     write.table(o_uniques[y][[1]], paste('./Finalize/', y, '_Group_Expression.txt', 
#                                          sep = ''), quote = FALSE, row.names = FALSE)
#   })
#   group_unique <- c(a_uniques_l, h_uniques_l, o_uniques_l)
#   non_unique <- colSums(exp_df) - uniques_l - group_unique
#   
#   unique_df <- data.frame(value = c(uniques_l, group_unique, non_unique),
#                           cell = c(pref_order, pref_order, pref_order),
#                           type = c(rep('unique', 12), rep('group', 12), rep('non', 12)))
#   unique_df$type <- factor(unique_df$type, levels = c('unique',
#                                                       'group',
#                                                       'non'))
#   ggplot(unique_df, aes(x = cell, y = value, fill = type)) +
#     geom_bar(stat = 'identity', color = 'black') +
#     theme_classic() +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 6500)) +
#     ggtitle('Uniquely Detected Proteins') +
#     xlab('Cell') +
#     ylab('Count Protein') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_fill_manual(NULL, labels = c(unique = 'Unique to cell',
#                                        group = 'Unique within group',
#                                        non = 'Not unique'),
#                       values=c(unique = 'tomato',
#                                group = 'springgreen3',
#                                non = 'skyblue'))
# }
# figure1()
# 
# # Gene Expression Profile #
# figure2 <- function() {
#   exp_profile = data.frame(sapply(HSPC_order, function(x) {
#     r_exp = sum(r_df[, x] > 0)
#     p_exp = sum(p_df[, x] > 0)
#     rp_exp <- sum((r_common[, x] > 0) * (p_common[, x] > 0) == 1)
#     c(rna = r_exp - rp_exp, prot = p_exp - rp_exp, both = rp_exp)
#   }))
#   total_profile <- data.frame(colSums(exp_profile))
#   colnames(total_profile) <- 'number'
#   total_profile$cell <- rownames(total_profile)
#   HSC <- data.frame(exp_profile$HSC)
#   HSC$cell <- 'HSC'
#   colnames(HSC) <- c('value', 'cell')
#   Ma <- data.frame(exp_profile$MPPa)
#   Ma$cell <- 'MPPa'
#   colnames(Ma) <- c('value', 'cell')
#   Mb <- data.frame(exp_profile$MPPb)
#   Mb$cell <- 'MPPb'
#   colnames(Mb) <- c('value', 'cell')
#   Mc <- data.frame(exp_profile$MPPc)
#   Mc$cell <- 'MPPc'
#   colnames(Mc) <- c('value', 'cell')
#   comb <- do.call(rbind, list(HSC, Ma, Mb, Mc))
#   comb$type <- rep(c('rna', 'prot', 'both'), 4)
#   comb$type <- factor(comb$type, levels = c('rna',
#                                             'prot',
#                                             'both'))
#   write.csv(comb, './Finalize/gene_exp_profile.csv', row.names = FALSE)
#   ggplot(comb, aes(x = cell, y = value, fill = type)) +
#     theme_classic() +
#     geom_bar(color = '#000000', stat = 'identity', width = 0.5) + 
#     annotate("text", x = total_profile$cell, 
#              y= total_profile$number + 600, label = sapply(total_profile$number, function(d) {
#                paste('T = ', d)
#              }), size = 4.8) +
#     ggtitle('Gene Expression Profile') +
#     xlab('Cell') +
#     ylab('Count') + 
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text = element_text(size = 14)) +
#     scale_y_continuous(limits = c(0, 17000), expand = c(0, 0)) +
#     scale_fill_manual(NULL, labels = c(rna = 'RNA Only',
#                                        prot = 'Protein Only',
#                                        both = 'Both'),
#                       values=c(rna = 'springgreen3',
#                                prot = 'skyblue',
#                                both = 'tomato'))
# }
# figure2()
# 
# 
# 
# # Log2 mRNA values for unique HSC miRNA Targets #
# figure7 <- function() {
#   rcom <- (r_df[, 'HSC'] > 0) * (r_df[, 'MPPa'] > 0)
#   rcom <- rcom == 1
#   rcom <- log2(r_df[rcom, c('HSC', 'MPPa')])
#   p1 <- ggplot(rcom[miu_HSC_names,], aes(x = HSC, y = MPPa)) +
#     geom_point() +
#     theme_classic() +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     xlim(c(-7, 14)) + ylim(c(-7, 14)) +
#     ggtitle('Unique HSC miRNA targets') +
#     xlab('HSC mRNA (log2)') +
#     ylab('MPPa mRNA (log2)')
#   
#   rcom <- (r_df[, 'HSC'] > 0) * (r_df[, 'MPPb'] > 0)
#   rcom <- rcom == 1
#   rcom <- log2(r_df[rcom, c('HSC', 'MPPb')])
#   p2 <- ggplot(rcom[miu_HSC_names,], aes(x = HSC, y = MPPb)) +
#     geom_point() +
#     theme_classic() +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     xlim(c(-7, 14)) + ylim(c(-7, 14)) +
#     ggtitle('Unique HSC miRNA targets') +
#     xlab('HSC mRNA (log2)') +
#     ylab('MPPb mRNA (log2)')
#   
#   rcom <- (r_df[, 'HSC'] > 0) * (r_df[, 'MPPc'] > 0)
#   rcom <- rcom == 1
#   rcom <- log2(r_df[rcom, c('HSC', 'MPPc')])
#   p3 <- ggplot(rcom[miu_HSC_names,], aes(x = HSC, y = MPPc)) +
#     geom_point() +
#     theme_classic() +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     xlim(c(-7, 14)) + ylim(c(-7, 14)) +
#     ggtitle('Unique HSC miRNA targets') +
#     xlab('HSC mRNA (log2)') +
#     ylab('MPPc mRNA (log2)')
#   list(p1, p2, p3)
# }
# figure7()
# 
# 
# 
# # Filter miRNAs for microarray #
# figure9 <- function() {
#   r_both <- (r_common[, 'HSC'] > 0) * (r_common[, 'MPPa'] > 0) == 1
#   p_both <- p_common[r_both, c('HSC', 'MPPa')]
#   r_both <- r_common[r_both, c('HSC', 'MPPa')]
#   
#   rc_exp <- data.frame(sapply(r_both, function(x) rownames(r_both)[x > 0]))
#   rc_exp_c <- data.frame(sapply(r_both, function(x) sum(x > 0)))
#   
#   rpc_exp <- sapply(colnames(r_both), function(x) {
#     rownames(r_both)[((r_both[, x] > 0) * (p_both[, x] > 0)) == 1]
#   })
#   rpc_exp_c <- sapply(colnames(r_both), function(x) {
#     sum(((r_both[, x] > 0) * (p_both[, x] > 0)) == 1)
#   })
#   
#   roc_exp <- sapply(colnames(r_both), function(x) {
#     setdiff(rc_exp[x][[1]], rpc_exp[x][[1]])
#   })
#   roc_exp_c <- sapply(colnames(r_both), function(x) {
#     length(setdiff(rc_exp[x][[1]], rpc_exp[x][[1]]))
#   })
#   hpmr_names <- intersect(rpc_exp$HSC, roc_exp$MPPa)
#   hrmp_names <- intersect(roc_exp$HSC, rpc_exp$MPPa)
#   prob <- roc_exp_c / (roc_exp_c + rpc_exp_c)
#   hpmr <- (1-prob['HSC']) * prob['MPPa']
#   hrmp <- (1-prob['MPPa']) * prob['HSC']
#   
#   p <- length(hrmp_names) / (length(hrmp_names) + length(hpmr_names))
#   n = (length(hpmr_names) + length(hrmp_names))
#   
#   miu_HSC_oTF <- sapply(names(miu_HSC_overlap), function(x) {
#     x <- gene_set[x][[1]]
#     hp <- length(intersect(hpmr_names, x))
#     hr <- length(intersect(hrmp_names, x))
#     n <- hp + hr
#     n >= 2
#   })
#   miu_HSC_oTF <- miu_HSC_overlap[miu_HSC_oTF]
#   gs_prediction <- sapply(names(miu_HSC_oTF), function(x){
#     x <- gene_set[x][[1]]
#     hp <- length(intersect(hpmr_names, x))
#     hr <- length(intersect(hrmp_names, x))
#     n <- hp + hr
#     new_p <- hr / n
#     sig <- sqrt(p * (1 - p) / n)
#     return((new_p - p) / sig)
#   })
#   
#   miu_df <- data.frame(cbind(miu_HSC_oTF, gs_prediction))
#   colnames(miu_df) <- c('count', 'z')
#   miu_df$p.value <- 2 * pnorm(-abs(miu_df$z))
#   miu_df$prev <- 'YES'
#   miu_agg <- aggregate(miu_df$p.value, list(miu_df$count), max)
#   count_co <- min(miu_agg$Group.1[(miu_agg$x < 0.05)])
#   
#   both_crit <- (miu_df$p.value < 0.05) == 1
#   
#   miu_df$label <- rownames(miu_df) %in% c('mmu.miR.19a.3p', 
#                                           'mmu.miR.19a.5p',
#                                           'mmu.miR.125a.5p',
#                                           'mmu.miR.130a.3p')
#   
#   miu_df[!both_crit, 'prev'] <- 'NOT'
#   miu_df$prev[miu_df$label] <- 'COOL'
#   miu_df$prev <- factor(miu_df$prev, levels = c('NOT', 'YES', 'COOL'))
#   miu_df$p.log <- log2(miu_df$p.value)
#   miu_df$label[!both_crit] <- FALSE
#   miu_df$name <- sapply(rownames(miu_df), function(x){
#     gsub('mmu.', '', x)
#   })
#   
#   for_write <- miu_df
#   for_write$name <- NULL
#   for_write$label <- NULL
#   for_write$prev[for_write$prev == 'COOL'] <- 'YES'
#   colnames(for_write)[colnames(for_write) == 'prev'] <- 'Significant'
#   
#   write.csv(for_write, './Finalize/miuHSC_count_z-score.csv', 
#             row.names = TRUE, quote = FALSE)
#   
#   ggplot(miu_df, aes(x = count, y = z, color = prev)) +
#     geom_point(data = subset(miu_df, prev == 'NOT'),
#                aes(x = count, y = z, color = 'NOT')) +
#     geom_point(data = subset(miu_df, prev == 'YES'),
#                aes(x = count, y = z, color = 'YES')) +
#     geom_point(data = subset(miu_df, prev == 'COOL'),
#                aes(x = count, y = z, color = 'COOL')) +
#     geom_text(data = subset(miu_df, prev == 'COOL'),
#               aes(x = count, y = z, label = name), size = 4.5,
#               color = 'black', nudge_x = 40) +
#     theme_classic() +
#     geom_segment(x = 0, y = log2(0.05), xend = 250, yend = log2(0.05), color = 'red') + 
#     geom_segment(x = count_co, y = -60, xend = count_co, yend = 5, color = 'blue') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, 
#           legend.title = element_blank(), legend.text = element_text(size = 14)) +
#     ggtitle('HSC mRNA only enrichment per miRNA target gene set') +
#     xlab('Count genes overlap') +
#     ylab('Log Adjusted p-value') + 
#     scale_color_manual(NULL, labels = c(NOT = 'Not significant',
#                                        YES = 'Significant',
#                                        COOL = 'Implicated'),
#                       values=c(NOT = '#a3a3a3',
#                                YES = 'springgreen3',
#                                COOL = 'tomato3'))
# }
# figure9()
# 
# 
# 
# # Generates csv for specific genes #
# figure11 <- function(genes, name = TRUE) {
#   if(name == FALSE) {
#     gene <- maps$gene[maps$entrez == gene]
#   }
#   entrez <- maps$entrez[which(maps$gene %in% genes)]
#   RNA <- raw_df_r[as.character(entrez), which(raw_r_order %in% HSPC_order)]
#   colnames(RNA) <- paste(raw_r_order, 'mRNA')
#   rownames(RNA) <- sapply(rownames(RNA), function(x){
#     maps$gene[maps$entrez == x]
#   })
#   PROT <- norm_df[as.character(entrez), which(raw_order %in% HSPC_order)]
#   colnames(PROT) <- paste(sapply(colnames(PROT), function(x) sub('[0-9].+', '', x))
#                           , 'Protein')
#   rownames(PROT) <- sapply(rownames(PROT), function(x){
#     maps$gene[maps$entrez == x]
#   })
#   bound <- cbind(RNA, PROT)
#   write.csv(bound, './Finalize/Expression_Information.csv', row.names = TRUE, quote = FALSE)
# }
# figure11(c('Adnp', 'Dnmt3a', 'Hprt'))
# 
# # mRNA Expression Profile of expressed proteins that are miRNA targets #
# figure12 <- function() {
#   prot<-p_df[, HSPC_order]
#   rna <-r_df[, HSPC_order]
#   
#   prot_exp_name <- sapply(HSPC_order, function(x) {
#     rownames(prot)[prot[, x] > 0]
#   })
#   prot_exp_name_c <- sapply(HSPC_order, function(x) {
#     length(rownames(prot)[prot[, x] > 0])
#   })
#   rna_exp_name <- sapply(HSPC_order, function(x) {
#     rownames(rna)[rna[, x] > 0]
#   })
#   rna_exp_name_c <- sapply(HSPC_order, function(x) {
#     length(rownames(rna)[rna[, x] > 0])
#   })
#   pr_exp_name <- sapply(HSPC_order, function(x) {
#     intersect(prot_exp_name[x][[1]], rna_exp_name[x][[1]])
#   })
#   pr_exp_name_c <- sapply(HSPC_order, function(x) {
#     length(pr_exp_name[x][[1]])
#   })
#   po_exp_name <- sapply(HSPC_order, function(x) {
#     setdiff(prot_exp_name[x][[1]], rna_exp_name[x][[1]])
#   })
#   po_exp_name_c <- sapply(HSPC_order, function(x) {
#     length(setdiff(prot_exp_name[x][[1]], rna_exp_name[x][[1]]))
#   })
#   
#   miu_over <- sapply(pr_exp_name, function(x){
#     intersect(x, miu_HSC_names)
#   })
#   miu_over_c <- data.frame(sapply(miu_over, length))
#   colnames(miu_over_c) <- 'value'
#   miu_over_c$cell <- rownames(miu_over_c)
#   miu_over_c$type <- 'both'
#   miu_over_p <- sapply(po_exp_name, function(x) {
#     intersect(x, miu_HSC_names)
#   })
#   miu_over_pc <- data.frame(sapply(miu_over_p, length))
#   colnames(miu_over_pc) <- 'value'
#   miu_over_pc$cell <- rownames(miu_over_pc)
#   miu_over_pc$type <- 'prot'
#   
#   df <- rbind(miu_over_c, miu_over_pc)
#   df$type <- factor(df$type, levels = c('prot', 'both'))
#   
#   ggplot(df, aes(x = cell, y = value, fill = type)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') +
#     scale_fill_manual(NULL, labels = c(prot = 'Protein Only', both = 'Protein & mRNA'),
#                       values = c(both = 'tomato1', prot = 'springgreen3')) +
#     theme_classic() +
#     annotate('text', x = miu_over_c$cell, y = miu_over_c$value + 20, 
#              label = paste(miu_over_pc$value, ' / ', 
#                            miu_over_pc$value + miu_over_c$value , sep = ''),
#              size = 4.8) +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, 
#           legend.title = element_blank(), legend.text = element_text(size = 14)) +
#     ggtitle('Expression Profile of miRNA targets with no proteins uniquely in HSCs') +
#     xlab('Cell') +
#     ylab('Count') +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 900))
# }
# figure12()
# 
# # miRNA Target Distribution in Aged #
# figure13 <- function() {
#   p_int <- p_df[which(rownames(p_df) %in% miRNAs), 
#                 c('AgedHSC', 'AgedMPPa', 'AgedMPPb', 'AgedMPPc',
#                   'HSC', 'MPPa', 'MPPb', 'MPPc')]
#   p_exp <- sapply(p_int, function(x) {
#     rownames(p_int)[x > 0]
#   })
#   p_expc <- sapply(p_int, function(x) {
#     length(rownames(p_int)[x > 0])
#   })
#   p_h_exp <- sapply(p_int, function(x) {
#     intersect(rownames(p_int)[x > 0], miu_HSC_names)
#   })
#   p_h_expc <- sapply(p_int, function(x) {
#     length(intersect(rownames(p_int)[x > 0], miu_HSC_names))
#   })
#   p_n_expc <- p_expc - p_h_expc
#   
#   p_h_expc <- data.frame(p_h_expc)
#   colnames(p_h_expc) <- 'value'
#   p_h_expc$cell <- rownames(p_h_expc)
#   p_h_expc$type <- 'nHSC'
#   p_n_expc <- data.frame(p_n_expc)
#   colnames(p_n_expc) <- 'value'
#   p_n_expc$cell <- rownames(p_n_expc)
#   p_n_expc$type <- 'Others'
#   
#   comb_df <- data.frame(rbind(p_h_expc, p_n_expc))
# 
#   nAnO <- setdiff(miu_HSC_names, p_exp$AgedHSC)
#   nAnO <- maps$gene[which(maps$entrez %in% nAnO)]
#   
#   nAyO <- intersect(miu_HSC_names, p_exp$AgedHSC)
#   nAyO <- maps$gene[which(maps$entrez %in% nAyO)]
#   
#   ggplot(comb_df, aes(x = cell, y = value, fill = type)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
#     theme_classic() +
#     ggtitle('miRNA Target Distribution') +
#     annotate('text', x = p_h_expc$cell, y = p_h_expc$value + p_n_expc$value + 250,
#              label = p_h_expc$value, size = 4.8) +
#     xlab('Cell') +
#     ylab('Targets Detected') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 5000)) +
#     scale_fill_manual(NULL, labels = c(nHSC = 'Uniquely not in HSC',
#                                        Others = 'Others'),
#                       values=c(nHSC = 'gold',
#                                Others = 'tomato'))
# }
# figure13()
# 
# # miRNA Target Distribution Expanded #
# figure14 <- function() {
#   p_int <- p_df[which(rownames(p_df) %in% miRNAs), pref_order]
#   p_exp <- sapply(p_int, function(x) {
#     rownames(p_int)[x > 0]
#   })
#   p_expc <- sapply(p_int, function(x) {
#     length(rownames(p_int)[x > 0])
#   })
#   p_h_exp <- sapply(p_int, function(x) {
#     intersect(rownames(p_int)[x > 0], miu_HSC_names)
#   })
#   p_h_expc <- sapply(p_int, function(x) {
#     length(intersect(rownames(p_int)[x > 0], miu_HSC_names))
#   })
#   p_hma_exp <- sapply(p_int, function(x) {
#     intersect(rownames(p_int)[x > 0], miu_HMa_names)
#   })
#   p_hma_expc <- sapply(p_int, function(x) {
#     length(intersect(rownames(p_int)[x > 0], miu_HMa_names))
#   })
#   p_hmb_exp <- sapply(p_int, function(x) {
#     intersect(rownames(p_int)[x > 0], miu_HMb_names)
#   })
#   p_hmb_expc <- sapply(p_int, function(x) {
#     length(intersect(rownames(p_int)[x > 0], miu_HMb_names))
#   })
#   p_n_expc <- p_expc - p_h_expc - p_hma_expc - p_hmb_expc
#   
#   p_h_expc <- data.frame(p_h_expc)
#   colnames(p_h_expc) <- 'value'
#   p_h_expc$cell <- rownames(p_h_expc)
#   p_h_expc$type <- 'nHSC'
#   p_hma_expc <- data.frame(p_hma_expc)
#   colnames(p_hma_expc) <- 'value'
#   p_hma_expc$cell <- rownames(p_hma_expc)
#   p_hma_expc$type <- 'nHMa'
#   p_hmb_expc <- data.frame(p_hmb_expc)
#   colnames(p_hmb_expc) <- 'value'
#   p_hmb_expc$cell <- rownames(p_hmb_expc)
#   p_hmb_expc$type <- 'nHMb'
#   p_n_expc <- data.frame(p_n_expc)
#   colnames(p_n_expc) <- 'value'
#   p_n_expc$cell <- rownames(p_n_expc)
#   p_n_expc$type <- 'Others'
#   
#   comb_df <- data.frame(rbind(rbind(rbind(p_h_expc,p_hma_expc), p_hmb_expc), p_n_expc))
#   
#   comb_df$type <- factor(comb_df$type, levels = c('nHMb',
#                                                   'nHMa',
#                                                   'nHSC',
#                                                   'Others'))
#   
#   comb_df$cell <- factor(comb_df$cell, levels = pref_order)
#   
#   ggplot(comb_df, aes(x = cell, y = value, fill = type)) +
#     geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
#     theme_classic() +
#     ggtitle('miRNA Target Distribution') +
#     xlab('Cell') +
#     ylab('Targets Detected') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_y_continuous(expand = c(0, 0), lim = c(0, 6000)) +
#     scale_fill_manual(NULL, labels = c(nHSC = 'Not in HSC Only',
#                                        Others = 'mRNA Only',
#                                        nHMa = 'Not in HSC & MPPa Only',
#                                        nHMb = 'Not in HSC & MPPb Only'),
#                       values=c(nHSC = 'gold',
#                                Others = 'tomato',
#                                nHMa = 'springgreen2',
#                                nHMb = 'skyblue'))
# }
# figure14()
# 
# 
# 
# 
# 
# 
# 
# 
# # Z-score of miRNA target gene enrichment within unique compartments #
# figure19 <- function() {
#   p = length(intersect(miRNAs, rownames(r_df))) / length(rownames(r_df))
#   n = length(rownames(r_df))
#   mu = n * p
#   sig = sqrt(n * p * (1 - p))
#   
#   m_w <- sapply(HSPC_order, function(x) {
#     rownames(r_df)[r_df[, x] > 0]
#   })
#   p_w <- sapply(HSPC_order, function(x) {
#     rownames(p_df)[p_df[, x] > 0]
#   })
#   
#   mp_w <- sapply(HSPC_order, function(x) {
#     intersect(m_w[x][[1]], p_w[x][[1]])
#   })
#   mp_wc <- sapply(HSPC_order, function(x) {
#     length(intersect(m_w[x][[1]], p_w[x][[1]]))
#   })
#   
#   mo_w <- sapply(HSPC_order, function(x) {
#     setdiff(m_w[x][[1]], p_w[x][[1]])
#   })
#   mo_wc <- sapply(HSPC_order, function(x) {
#     length(setdiff(m_w[x][[1]], p_w[x][[1]]))
#   })
#   
#   po_w <- sapply(HSPC_order, function(x) {
#     setdiff(p_w[x][[1]], m_w[x][[1]])
#   })
#   po_wc <- sapply(HSPC_order, function(x) {
#     length(setdiff(p_w[x][[1]], m_w[x][[1]]))
#   })
#   
#   hsc <- intersect(setdiff(setdiff(setdiff(mo_w$HSC, mo_w$MPPa), mo_w$MPPb), mo_w$MPPc),
#                    union(union(mp_w$MPPa, mp_w$MPPb), mp_w$MPPc))
#   mppa <- intersect(setdiff(setdiff(setdiff(mo_w$MPPa, mo_w$HSC), mo_w$MPPb), mo_w$MPPc),
#                     union(union(mp_w$HSC, mp_w$MPPb), mp_w$MPPc))
#   mppb <- intersect(setdiff(setdiff(setdiff(mo_w$MPPb, mo_w$MPPa), mo_w$HSC), mo_w$MPPc),
#                     union(union(mp_w$MPPa, mp_w$HSC), mp_w$MPPc))
#   mppc <- intersect(setdiff(setdiff(setdiff(mo_w$MPPc, mo_w$MPPa), mo_w$MPPb), mo_w$HSC),
#                     union(union(mp_w$MPPa, mp_w$MPPb), mp_w$HSC))
#   
#   hsc_exp <- length(hsc) * p
#   hsc_act <- length(intersect(hsc, miRNAs))
#   zhsc <- sqrt(length(hsc)) * (hsc_act - hsc_exp) / sqrt(length(hsc) * p * (1 - p))
#   ma_exp <- length(mppa) * p
#   ma_act <- length(intersect(mppa, miRNAs))
#   zma <- sqrt(length(mppa)) * (ma_act - ma_exp) / sqrt(length(mppa) * p * (1 - p))
#   mb_exp <- length(mppb) * p
#   mb_act <- length(intersect(mppb, miRNAs))
#   zmb <- sqrt(length(mppb)) * (mb_act - mb_exp) / sqrt(length(mppb) * p * (1 - p))
#   mc_exp <- length(mppc) * p
#   mc_act <- length(intersect(mppc, miRNAs))
#   zmc <- sqrt(length(mppc)) * (mc_act - mc_exp) / sqrt(length(mppc) * p * (1 - p))
#   
#   zs <- data.frame(c(zhsc, zma, zmb, zmc), HSPC_order)
#   colnames(zs) <- c('zscore', 'cell')
#   
#   ggplot(zs, aes(x = cell, y = zscore)) +
#     geom_bar(stat = 'identity', width = 0.5) +
#     annotate('text', x = zs$cell, y = zs$zscore + 2, label = round(zs$zscore, 3), size = 4.8) +
#     theme_classic() +
#     ggtitle('z-score of miRNA target enrichment') +
#     xlab('z-score') +
#     ylab('Cell') +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, legend.text=element_text(size=14),
#           axis.text.x = element_text(angle = 45, hjust = 1)) + 
#     scale_y_continuous(expand = c(0, 0), lim = c(-0.006, 65))
# }
# figure19()
# 
# 
# 
# 
# 
# 
# 
# 
# # Proportion of total mRNA that are miRNA targets #
# figure24 <- function() {
#   expressed <- sapply(r_df, function(x) {
#     rownames(r_df)[x > 0]
#   })
#   intersected <- sapply(expressed, function(x) {
#     intersect(x, unique(miRNAs))
#   })
#   
#   proportions <- sapply(c('HSC', 'MPPa', 'MPPb', 'MPPc'), function(x){
#     return(length(intersected[x][[1]]) / length(expressed[x][[1]]))
#   })
#   
#   write.csv(proportions, './Finalize/miRNA_proportion.csv', row.names = TRUE,
#             quote = FALSE)
# }
# figure24()
# 
# figure25 <- function() {
#   idf <- read.csv('./Finalize/Image Quant/Formatted.csv')
#   idf$X <- factor(idf$X, levels = c('KLS', 'HSC', 'cHSC'))
#   ggplot(idf, aes(x = Ki67.DAPI, y = Dnmt3a.DAPI, color = X)) +
#     geom_point() +
#     theme_classic() +
#     xlab('KI67_CTCF / DAPI_CTCF') +
#     ylab('DNMT3a_CTCF / DAPI_CTCF') +
#     scale_color_manual(NULL, labels = c(cHSC = "Cultred HSC", 
#                                         HSC = 'HSC',
#                                         KLS = 'KLS'), 
#                        values = c(cHSC = 'springgreen3', 
#                                   HSC = 'tomato3',
#                                   KLS = 'steelblue')) +
#     theme(axis.text=element_text(size=14, colour = 'black'), 
#           axis.title=element_text(size=14),
#           plot.title = element_text(size = 14),
#           aspect.ratio = 1, 
#           legend.title = element_text('Cell Type'), legend.text = element_text(size = 14)) +
#     ggtitle('DNMT3A vs KI67 Fluorescence')
# }
# figure25()
# 
# # Explore other two files!
# # Write methods once figures are all set
# # Explore the unique gene set stuff
# # miRNA separate list! With count & rank!
# 
# # After Balyn, color green table
# 
# # Validation List : Single Protein Detection Cases
# # Underexpressed
# # Non expressed
# # I think Notch1 was one of them?
# # Check every step for potential validation
# 
# # Why are there... Different mappings & different gene names by count
# 
