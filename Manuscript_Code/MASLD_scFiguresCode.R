## Liver single cell RNAseq
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library('ggplot2')
library('RColorBrewer')
library(readxl)  
library(grid)
library(patchwork)
library(gridExtra)

multiplesheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # print data frame
  return(data_frame)
}

theme = theme(legend.position = "bottom") + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank())

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

## SET PARAMETERS ----
options(stringsAsFactors = FALSE)
set.seed(1212)
INDIR = "/path/to/workdir"
if(!dir.exists(file.path(INDIR, 'Output'))){  dir.create(file.path(INDIR, 'Output')) }
OUTDIR = file.path(INDIR, 'Output')
date = Sys.Date()
filter <- dplyr::filter
select <- dplyr::select



########################################
##### Proteomics and metadata
########################################
## read coefficients and sig proteins 
ravi_proteomics_coef <- read.csv(file.path(INDIR, 'nafld_lasso_coefs.csv'))

# specifying the path name
xlsxfiles <- file.path(INDIR, 'HR3LAMEAN_Cardia_lm_protein.xlsx')
ravi_proteomics_fdr <- multiplesheets(xlsxfiles)

hu_meta_file <- file.path(INDIR, 'nflad_hu_metadata.xlsx')
hu_metadata <- multiplesheets(hu_meta_file)[[1]]
hu_metadata <- na.omit(hu_metadata)

## calculate expression of significant genes
sig_proteomics_genes <- ravi_proteomics_fdr$ageSexRaceBMI %>% 
  filter(der_fdr < 0.05 & val_fdr < 0.05)
sig_proteomics_genes <- sig_proteomics_genes[order(sig_proteomics_genes$der_fdr),]


########################################
##### sc RNA seq object - post process
########################################

## scale data for protein of interests
seuratObj <- readRDS(file.path(INDIR, 'seuratObj_humanRNAseq_V5COMPATIBLE.rds'))
protemics_genes <- ravi_proteomics_fdr$ageSexRaceBMI$EntrezGeneSymbol %>% unique()
protemics_genes <- protemics_genes[protemics_genes %in% rownames(seuratObj)]
protemics_genes <- unique(c(protemics_genes, 'OPA1', 'ULK1'))
seuratObj <- ScaleData(seuratObj, feature = protemics_genes)


## check which targets exsist in liver sc data
sig_proteomics_genes_scexpressed <- sig_proteomics_genes[ which(sig_proteomics_genes$EntrezGeneSymbol %in% rownames(seuratObj)), ]
sig_proteomics_genes_scexpressed_list <-  sig_proteomics_genes_scexpressed$EntrezGeneSymbol %>% unique()

########################################
##### Visium
########################################
##### Load seuratObj
seuratObjMerge<-readRDS(file=file.path(INDIR, 'seuratObj_humanVisium.rds'))
seuratObjMerge <- ScaleData(seuratObjMerge, feature = protemics_genes)
dim(seuratObjMerge)
# 16873  4333

##### Load metadata
annotHumanVisium<-read.table(file=file.path(INDIR, 'annot_humanVisium.csv'), sep = ',',
                             header = T, stringsAsFactors = F)
dim(annotHumanVisium)
# 4333    8
rownames(annotHumanVisium)<-annotHumanVisium$spot
##### Merge metaData
seuratObjMerge@meta.data$zonation<-annotHumanVisium[rownames(seuratObjMerge@meta.data),'zonation']
seuratObjMerge@meta.data$zonationGroup<-annotHumanVisium[rownames(seuratObjMerge@meta.data),'zonationGroup']
seuratObjMerge@meta.data$diet<-annotHumanVisium[rownames(seuratObjMerge@meta.data),'diet']
seuratObjMerge@meta.data$diet <- factor(seuratObjMerge@meta.data$diet, levels = c('healthy', 'fatty'))


## get model gene overlapping with visium gex ----
sig_proteomics_genes_vis.expressed <- sig_proteomics_genes[ which(sig_proteomics_genes$EntrezGeneSymbol %in% rownames(seuratObjMerge)), ]
sig_proteomics_genes_visexpressed_list <-  sig_proteomics_genes_vis.expressed$EntrezGeneSymbol %>% unique()
## final gene list visium, model sig, scRNAseq
sig_proteomics_genes_liver_atlas_expressed <- intersect(sig_proteomics_genes_scexpressed_list, sig_proteomics_genes_visexpressed_list)
length(sig_proteomics_genes_liver_atlas_expressed)
# 198

## calculated expression score from genes of interest## check which targets exsist in liver sc data
goi <- sig_proteomics_genes_liver_atlas_expressed
which(goi %in% rownames(seuratObjMerge)) %>% length()
seuratObjMerge <- AddModuleScore(seuratObjMerge,
                                 features = list(goi), slot = 'scale.data',  
                                 name="Liver_proteomics_enriched")

seuratObj <- AddModuleScore(seuratObj,assay = 'RNA', 
                            features = list(goi),slot = 'scale.data', 
                            name="Liver_proteomics_enriched")

# FIGURE 1A.1 ----
## visualize score across clusters
p1 <- FeaturePlot(seuratObj,
                  features = "Liver_proteomics_enriched1", label = TRUE, repel = TRUE,  max.cutoff = 'q99',
                  raster = TRUE, raster.dpi = c(2048,2048), pt.size =  2, label.size = 5) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

Figure1a <- p1 + 
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow()),
        axis.title = element_text(hjust = 0)) +
  theme(legend.position = "right") + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text = element_text(size = 14),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '')


p2 <- DimPlot(seuratObjMerge, label=T)
p3 <- DimPlot(seuratObjMerge, label=F, group.by = 'diet') + scale_color_manual(values=colorsDiet)
p4 <- DimPlot(seuratObjMerge, label=F, group.by = 'zonationGroup') + scale_color_manual(values=colorsZonationGroups)
## visualize score across clusters ----
p5 <- FeaturePlot(seuratObjMerge,
                  features = "Liver_proteomics_enriched1", label = F, repel = TRUE,  max.cutoff = 'q95', min.cutoff = 'q5', label.size = 5) +
  viridis::scale_color_viridis(option = 'A', direction = 1)

# FIGURE 1B.1 ----
Figure1b.1 <- p3 + 
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow()),
        axis.title = element_text(hjust = 0)) +
  theme(legend.position = "right") + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text = element_text(size = 14),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '')

# FIGURE 1B.2 ----
Figure1b.2 <- p5 + 
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow()),
        axis.title = element_text(hjust = 0)) +
  theme(legend.position = "right") + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text = element_text(size = 14),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.505))


# visualize enrichment score on tissue ----
sl1.empty <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                                images = 'slice1', pt.size.factor = 0) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.52)) +
  ggtitle(NULL) +
  labs(fill = "Expression score") +
  theme(legend.position="none")
sl2.empty <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                                images = 'slice1.1', pt.size.factor = 0) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.52)) +
  ggtitle(NULL) +
  labs(fill = "Expression score") +
  theme(legend.position="none")
sl3.empty <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                                images = 'slice1.2', pt.size.factor = 0) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.52)) +
  ggtitle(NULL)+
  labs(fill = "Expression score") +
  theme(legend.position="none")
sl4.empty <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                                images = 'slice1.3', pt.size.factor = 0) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.52)) +
  ggtitle(NULL)+
  labs(fill = "Expression score") +
  theme(legend.position="none")
sl5.empty <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                                images = 'slice1.4', pt.size.factor = 0) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.52)) +
  ggtitle(NULL)+
  labs(fill = "Expression score") +
  theme(legend.position="none")


sl1 <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                          images = 'slice1', pt.size.factor = 2) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.605)) +
  ggtitle('H35 (steatotic)') +
  labs(fill = "Expression score") +
  theme(legend.position="none",
        plot.title = element_text(size=10)) 
sl2 <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                          images = 'slice1.1', pt.size.factor = 3) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.605)) +
  ggtitle('H35 (steatotic)') +
  labs(fill = "Expression score") +
  theme(legend.position="none",
        plot.title = element_text(size=10))
sl3 <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                          images = 'slice1.2', pt.size.factor = 2) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.605)) +
  ggtitle('H36 (healthy)')+
  labs(fill = "Expression score") +
  theme(legend.position="none",
        plot.title = element_text(size=10)) 
sl4 <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                          images = 'slice1.3', pt.size.factor = 2) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.605)) +
  ggtitle('H37 (steatotic)')+
  labs(fill = "Expression score") +
  theme(legend.position="none",
        plot.title = element_text(size=10)) 
sl5 <- SpatialFeaturePlot(seuratObjMerge, features = c("Liver_proteomics_enriched1"),  
                          images = 'slice1.4', pt.size.factor = 2) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(.3,.605)) +
  ggtitle('H38 (healthy)')+
  labs(fill = "Expression score") +
  theme(legend.position="none",
        plot.title = element_text(size=10)) 
sl_plist <- list(sl3, sl5, sl1, sl2, sl4,
                 sl3.empty, sl5.empty, sl1.empty, sl2.empty, sl4.empty)
# FIGURE 1C.1 ----
Figure1c.1.supp <- arrangeGrob(grobs = sl_plist, ncol = 5)
sl_plist <- list(sl3, sl4,
                 sl3.empty, sl4.empty)
Figure1c.1 <- arrangeGrob(grobs = sl_plist, ncol = 5)


# FIGURE 1C.4 ----
Figure1c.4 <- seuratObjMerge@meta.data %>% 
  ggplot() +
  # geom_boxplot(aes(x = diet, y = Liver_proteomics_enriched1)) +
  ggforce::geom_sina(aes(x = diet, y = Liver_proteomics_enriched1, color = diet),
                     alpha = 0.25, scale = 'area', bins = 500, size = 0.5) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines")) +
  theme(text = element_text(size = 16),
        axis.ticks=element_blank(), 
        axis.title.x = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  scale_color_manual(values = c('healthy'='darkorchid4','fatty'='orange')) + 
  labs(y = 'Expression score')

wilcox.test(seuratObjMerge@meta.data$Liver_proteomics_enriched1 ~ seuratObjMerge@meta.data$diet)
# Wilcoxon rank sum test with continuity correction
# 
# data:  seuratObjMerge@meta.data$Liver_proteomics_enriched1 by seuratObjMerge@meta.data$diet
# W = 1037362, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


## plot protein set score by zonation group ----
zones <- unique(seuratObjMerge@meta.data$zonationGroup)
zone_coords <- c()
for(z in zones){
  zmin <- seuratObjMerge@meta.data %>% 
    filter(zonationGroup == z) %>% 
    pull(zonation) %>% min()
  zmax <- seuratObjMerge@meta.data %>% 
    filter(zonationGroup == z) %>% 
    pull(zonation) %>% max()
  zmed <- seuratObjMerge@meta.data %>% 
    filter(zonationGroup == z) %>% 
    pull(zonation) %>% median()
  tmp <- data_frame(zone = z, start = zmin, end = zmax, median = zmed)
  zone_coords <- rbind(zone_coords, tmp)
}
zone_coords <- zone_coords[order(zone_coords$end, decreasing = F),]
# FIGURE 1C.2 ----
Figure1c.2 <- seuratObjMerge@meta.data %>% 
  ggplot(aes(x = zonation, y = Liver_proteomics_enriched1, group = diet, color = diet)) +
  geom_smooth() + 
  theme_bw() +
  scale_x_continuous(
    breaks = sort(c(0,zone_coords$end))) +
  theme(panel.background=element_blank(),
        plot.background=element_blank(),
        text = element_text(size = 14),
        axis.ticks.length.x =  unit(2, 'mm'),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(x = '',
       y = 'Expression score', 
       color = NULL) +
  annotate(
    ymin = 0.335, ymax = 0.3435,
    xmin = zone_coords$start,
    xmax = zone_coords$end,
    geom = "rect",
    fill = colorsZonationGroups
  ) + 
  annotate("text",x= zone_coords$median[1]*0.75,y=.34,label=zone_coords$zone[1], size = 4.5)+
  annotate("text",x= zone_coords$median[2]*0.9,y=.34,label=zone_coords$zone[2], size = 4.5)+
  annotate("text",x= zone_coords$median[3],y=.34,label=zone_coords$zone[3], size = 4.5)+
  annotate("text",x= zone_coords$median[4]*0.99,y=.34,label=zone_coords$zone[4], size = 4.5)+ 
  coord_cartesian(ylim=c(0.35, 0.48),clip="off") + 
  scale_color_manual(values=colorsDiet)

grid.arrange(grobs = c(plist, list(Figure1c.2)), ncol = 3)



########## GET DE GENES ----
Idents(seuratObjMerge) <- 'diet'

## TO PLOT AS IS USE ident.1 = "healthy", ident.2 = "fatty",
all_disease_markers <- FindMarkers(seuratObjMerge, ident.1 = "healthy", ident.2 = "fatty", 
                                   logfc.threshold = 0, min.diff.pct = 0, min.pct = 0, only.pos = F, test.use = 'negbinom', latent.vars = 'orig.ident')
all_disease_markers$diff <- all_disease_markers$pct.1 - all_disease_markers$pct.2
all_disease_markers$symbol <- rownames(all_disease_markers)

### write list of de targets ----
# tmp <- all_disease_markers[which(all_disease_markers$symbol %in% sig_proteomics_genes_liver_atlas_expressed), ]
# tmp <- tmp[,c(7, 2, 1, 5, 3, 4, 6)]
# colnames(tmp) <- c('Symbol', 'log2FC', 'pvalue', 'fdr', 'pct.expr.healthy', 'pct.expr.fatty', 'pct.expr.diff')
# write.csv(tmp, file.path(OUTDIR, 'NAFLD_target_genes_de.csv'), row.names = F)
# tmp$is.DE <- ifelse(tmp$Symbol %in% disease_markers_l, 'yes', 'no')

disease_markers <- FindMarkers(seuratObjMerge, ident.1 = "healthy", ident.2 = "fatty", 
                               logfc.threshold = 0.25, min.diff.pct = 0.1, min.pct = 0.1, only.pos = F, test.use = 'negbinom', latent.vars = 'orig.ident')
disease_markers$diff <- disease_markers$pct.1 - disease_markers$pct.2
disease_markers$symbol <- rownames(disease_markers)
disease_markers <- disease_markers %>% filter(p_val_adj < 0.05)
disease_markers_l <- disease_markers$symbol[which(disease_markers$symbol %in% sig_proteomics_genes_scexpressed$EntrezGeneSymbol)]


# prepare DE for spatial visualization of top 5 genes (sc and proteomics) ----
head(sig_proteomics_genes)
sig_proteomics_genes_moi <- sig_proteomics_genes
sig_proteomics_genes_moi$Visium_l2fc <- all_disease_markers$avg_log2FC[match(sig_proteomics_genes$EntrezGeneSymbol, all_disease_markers$symbol)]
sig_proteomics_genes_moi$Visium_nucdif <- all_disease_markers$diff[match(sig_proteomics_genes$EntrezGeneSymbol, all_disease_markers$symbol)]
sig_proteomics_genes_moi$Visium_fdr <- all_disease_markers$p_val_adj[match(sig_proteomics_genes$EntrezGeneSymbol, all_disease_markers$symbol)]



## plot differentially expressed genes from Visium with model sig genes ----
pzone <- DotPlot(seuratObjMerge, features = disease_markers_l, group.by =  'zonationGroup')+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank()) 

pdiet <- DotPlot(seuratObjMerge, features = disease_markers_l, group.by =  'diet')+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank())

lay <- rbind(c(1,1),
             c(1,1),
             c(1,1),
             c(NA, NA),
             c(2,2),
             c(2,2),
             c(2,2),
             c(2,2))
pvis <- grid.arrange(grobs = list(pdiet, pzone), ncol = 1, layout_matrix = lay)
# FIGURE 1D.1 ----
Figure1d <- arrangeGrob(grobs = list(pdiet, pzone), ncol = 1, layout_matrix = lay)
pcelltype <- DotPlot(seuratObj, features = disease_markers_l, group.by =  'annot')+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank())
# FIGURE 1D.2 ----
Figure1d.2 <- pcelltype



# Plot all model genes ----
'%!in%' <- function(x,y)!('%in%'(x,y))
sig_proteomics_genes_moi$col <- 'lightgrey'
sig_proteomics_genes_moi$col[sig_proteomics_genes_moi$EntrezGeneSymbol %in% disease_markers_l] = '#ad1da4'
sig_proteomics_genes_moi$col[sig_proteomics_genes_moi$EntrezGeneSymbol %!in% disease_markers_l] = 'lightgrey'


der_beta_genes <- sig_proteomics_genes_moi %>% 
  filter(!is.na(Visium_fdr)) %>% 
  pull(EntrezGeneSymbol) %>% unique()

sig_proteomics_genes_moi_forplot <- c()
for(g in der_beta_genes){
  tmp <- subset(sig_proteomics_genes_moi, EntrezGeneSymbol == g)
  if(nrow(tmp) > 1){
    tmp$der_beta <- mean(tmp$der_beta)
    sig_proteomics_genes_moi_forplot <- rbind(sig_proteomics_genes_moi_forplot, tmp[1,])
  } else {
    sig_proteomics_genes_moi_forplot <- rbind(sig_proteomics_genes_moi_forplot, tmp) 
  }
}

# FIGURE 1E ----
Figure1e <- ggplot(sig_proteomics_genes_moi_forplot, aes(x = der_beta, y = Visium_l2fc)) +
  geom_point(color = sig_proteomics_genes_moi_forplot$col, size = 2) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.25) + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.25) +
  xlim(-0.3, 0.3) +
  ylim(-4.5,4.5) + 
  ggrepel::geom_text_repel(data = subset(sig_proteomics_genes_moi_forplot, col != 'lightgrey'),
                           aes(x = der_beta, y = Visium_l2fc, label = EntrezGeneSymbol), size = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = 'CARDIA model regression coefficient', 
       y = 'Log2 fold change (fatty/healthy)') 



# top 3 target - Visium disease markers intersect ----
## spatial ----
top3_targets_visium <- sig_proteomics_genes_moi[order(sig_proteomics_genes_moi$der_fdr), ] %>% 
  .[!is.na(.$Visium_fdr),] %>% 
  pull(EntrezGeneSymbol) %>% head(3)
top3_targets_visium <- c(top3_targets_visium, 'OPA1', 'ULK1')
comb_gene_plist <- list()
for(i in seq_along(top3_targets_visium)){
  gene <- top3_targets_visium[i]
  plot_limits <- quantile(seuratObjMerge@assays$SCT$data[gene,], c(0.01,0.99))
  voi <- c('slice1.2', 'slice1.4', 'slice1', 'slice1.1','slice1.3')
  gene_plist <- lapply(voi, function(var){
    p <- SpatialFeaturePlot(seuratObjMerge, features = gene, images = var, pt.size.factor = 2) + 
      #scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = plot_limits)
      scale_fill_gradientn(colours = brewer.pal(n = 11, name = "Purples"), limits = plot_limits)
    return(p)
  })
  comb_gene_plist[[i]] <- gene_plist
}
names(comb_gene_plist) <- top3_targets_visium
comb_gene_plist <- do.call(c, comb_gene_plist)
lay <- rbind(c(1,2,3,4,5),
             c(6,7,8,9,10),
             c(11,12,13,14,15))
geneplot <- grid.arrange(grobs = comb_gene_plist[1:15], layout_matrix = lay)
# FIGURE 2 ----
Figure2 <- arrangeGrob(grobs = comb_gene_plist[1:15], layout_matrix = lay)

## spatial violin ----
comb_vln_plist <- list()
spot_expression <- seuratObjMerge@meta.data
tagets <- c('IGFBP2', 'HMGCS1', 'ADH1A', 'AKR1C4', 'IGF1')
for(i in seq_along(tagets)){
  
  gene <- tagets[i]
  spot_expression$expression <- seuratObjMerge@assays$SCT$scale.data[gene,]
  p <- spot_expression %>% 
    ggplot() +
    # geom_boxplot(aes(x = diet, y = Liver_proteomics_enriched1)) +
    ggforce::geom_sina(aes(x = diet, y = expression, color = diet),
                       alpha = 0.25, scale = 'area', bins = 500, size = 0.5) +
    theme_bw() +
    theme(panel.spacing = unit(0, "lines")) +
    theme(text = element_text(size = 16),
          axis.ticks=element_blank(),
          axis.title.x = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position = "none") +
    scale_color_manual(values = c('healthy'='darkorchid4','fatty'='orange')) + 
    labs(y = paste(gene, 'Expression'))
  comb_vln_plist[[i]] <- p
}
names(comb_vln_plist) <- tagets
lay <- rbind(c(1,2,3,4,5))
geneplot <- grid.arrange(grobs = comb_vln_plist[1:5], layout_matrix = lay)
# Figure 2D ----
Figure2d <- arrangeGrob(grobs = comb_vln_plist[1:5], layout_matrix = lay)


## Save figures ----
Figure1a
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1a.pdf'), Figure1a, width = 6, height = 4.75, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1a.svg'), Figure1a, width = 6, height = 4.75, units = 'in', dpi = 300)
Figure1b.1
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1b_1.pdf'), Figure1b.1, width = 6, height = 4.75, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1b_1.svg'), Figure1b.1, width = 6, height = 4.75, units = 'in', dpi = 300)
Figure1b.2
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1b_2.pdf'), Figure1b.2, width = 6, height = 4.75, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1b_2.svg'), Figure1b.2, width = 6, height = 4.75, units = 'in', dpi = 300)

Figure1b <- Figure1b.1 + Figure1b.2
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1b.pdf'), Figure1b, width = 12, height = 4.75, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1b.svg'), Figure1b, width = 12, height = 4.75, units = 'in', dpi = 300)

plot(Figure1c.1)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_1.pdf'), Figure1c.1, width = 20, height = 4.75, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_1.svg'), Figure1c.1, width = 12, height = 18, units = 'in', dpi = 600)
plot(Figure1c.1.supp)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_1.supp.pdf'), Figure1c.1.supp, width = 20, height = 4.75, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_1.supp.svg'), Figure1c.1.supp, width = 30, height = 18, units = 'in', dpi = 600)
Figure1c.2
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_2.pdf'), Figure1c.2, width = 5, height = 4, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_2.svg'), Figure1c.2, width = 5, height = 4, units = 'in', dpi = 300)
Figure1c.4
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_4.pdf'), Figure1c.4, width = 4, height = 4, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1C_4.svg'), Figure1c.4, width = 4, height = 4, units = 'in', dpi = 300)
plot(Figure1d)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1D_1.pdf'), Figure1d, width = 14, height = 6, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1D_1.svg'), Figure1d, width = 14, height = 6, units = 'in', dpi = 300)
Figure1d.2
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1D_2.pdf'), Figure1d.2, width = 14, height = 6, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1D_2.svg'), Figure1d.2, width = 14, height = 6, units = 'in', dpi = 300)
Figure1e
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1E_flipped.pdf'), Figure1e, width = 7, height = 6, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure1E_flipped.svg'), Figure1e, width = 7, height = 6, units = 'in', dpi = 300)
plot(Figure2)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure2.pdf'), Figure2, width = 20, height = 13.5, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure2.svg'), Figure2, width = 20, height = 13.5, units = 'in', dpi = 300)
Figure1b.1 + Figure1b.2
plot(Figure2d)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure2d.pdf'), Figure2d, width = 17.5, height = 4, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'Figure2d.svg'), Figure2d, width = 17.5, height = 4, units = 'in', dpi = 300)



# Heatmap Supp ----
mean_express <- AverageExpression(object = seuratObj, group.by = "annot", assay = 'RNA')
mean_express <- mean_express$RNA %>% as.data.frame()
mean_express[1:5,1:5]
mean_express <- mean_express[which(rownames(mean_express) %in% goi),]

library(circlize)
library(EnrichedHeatmap)

rmeans <- rowMeans(mean_express)
rsd <- apply(mean_express, 1, function(x) sd(x))
mean_express_z <- (mean_express + rmeans)/rsd

col_fun = colorRampPalette(brewer.pal(8, "Purples"))(500)
h1 <- Heatmap(mean_express_z, column_names_rot = 45, column_names_side = "top", column_names_centered = F, 
              clustering_method_rows = 'ward.D2', clustering_distance_rows = 'euclidean',
              name = "expression (z-score)",  show_row_names = T, 
              show_column_names = T, row_names_side = "right",
              cluster_columns = T, col = col_fun, show_column_dend = F,
              ##what determines the color of the columns, tell it where you want to annotations --> left of right ,
              column_title_gp = gpar(fontsize = 12), row_title_rot = 0,row_title_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 5),
              show_row_dend = T, use_raster = TRUE, cluster_rows = T,
              heatmap_legend_param = list(direction = "horizontal", position = "bottom"))



svg(file.path(OUTDIR, "all_model_gene_expr.svg"), width = 7, height = 15)
draw(h1)
dev.off()






















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































