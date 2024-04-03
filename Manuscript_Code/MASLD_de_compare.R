library(Seurat)
library(dplyr)
library(Libra)
'%!in%' <- function(x,y)!('%in%'(x,y))


## SET PARAMETERS ----
options(stringsAsFactors = FALSE)
options(digits=4)
set.seed(1212)

seurat_de <- read.csv('/scratch/nhadad/NAFLD/Output/visium_differential_gene_expression_significant.csv')
seurat_digde <- seurat_de$EntrezGeneSymbol %>% unique()

## DE Pseudobulk ----
# run differential expression analysis using libra (pseudobulk)
seuratObjMerge@meta.data$cell_type = 'all'
psDE = run_de(seuratObjMerge, 
            de_family = 'mixedmodel', de_method = 'negbinom_offset', de_type = 'LRT',
            replicate_col = "orig.ident", label_col = 'diet', cell_type_col = 'cell_type', min_cells = 50)
de_sig <- psDE %>% filter(p_val_adj < 0.05)
which(de_sig$gene %in% seurat_digde)

length(de_sig)
print(de_sig)
# View(de_sig)

which(de_sig$gene %in% seurat_digde)



## DE nbmm ----
library(doParallel)
library(foreach)
options(stringsAsFactors = FALSE)
set.seed(1212)


'%!in%' <- function(x,y)!('%in%'(x,y))
get_test_stats <- function(glmm_model) {
  tmp_statistics <- parameters::model_parameters(glmm_model) %>% as.data.frame()
  tmp_std_es <- effectsize::standardize_parameters(glmm_model) %>% as.data.frame()
  tmp_statistics$std_Coefficient <- tmp_std_es$Std_Coefficient[match(tmp_statistics$Parameter, tmp_std_es$Parameter)]
  tmp_statistics$std_CI_low <- tmp_std_es$CI_low[match(tmp_statistics$Parameter, tmp_std_es$Parameter)]
  tmp_statistics$std_CI_high <- tmp_std_es$CI_high[match(tmp_statistics$Parameter, tmp_std_es$Parameter)]
  return(tmp_statistics)
}


gc()
nCores <- 20 # set number of cores for de

DefaultAssay(seuratObjMerge) <- "SCT"
meta <- seuratObjMerge@meta.data
counts <- seuratObjMerge@assays$SCT@counts

gene_totals <- (rowSums(counts == 0)/ncol(counts)) * 100
genes_to_include <- names(gene_totals[gene_totals < 90])

counts <- counts[which(rownames(counts) %in% sig_proteomics_genes_liver_atlas_expressed),]


log <- c()
library(doParallel)
registerDoParallel(cores=nCores)
# nrow(counts)
x <- foreach (i=1:nrow(counts), .combine = rbind) %dopar% {
  
  tryCatch({
    require(glmmTMB); require("bbmle"); require("dplyr")
    options(stringsAsFactors = FALSE)
    g <- rownames(counts)[i]
    
    meta$expression <- counts[g,] %>% unname() %>% unlist()
    ## ~~~~~~ modify here ~~~~~~ ##
    # options(contrasts = c("contr.sum","contr.poly"))
    mod.tbm0 <- glmmTMB::glmmTMB(expression ~ diet + (1|orig.ident), data = meta, family=glmmTMB::nbinom2(link = "log"), ziformula= ~0)

    
    library(parameters)
    library(effectsize)
    options(es.use_symbols = TRUE)
    
    ## get standarize effect size
    mod0_test_sts <- get_test_stats(mod.tbm0)
    mod0_test_sts$gene <- g
    
    return(mod0_test_sts)
    
  },  error = function(e){})
}
stopImplicitCluster()
x <- x %>% filter(Parameter == 'diet1')
x$Coefficient[x$Coefficient > 2] = NA
View(x)
# write.csv(x, file.path(OUTDIR, paste0("ChrY_", ct_out ,'_cell-agnostic_new_difexp90pct.csv')))


tmp <- sig_proteomics_genes_moi_forplot
tmp$glmm_coef <- x$Coefficient[match(tmp$EntrezGeneSymbol, x$gene)]
tmp$glmm_unadj.p <- x$p[match(tmp$EntrezGeneSymbol, x$gene)]
View(tmp)
tmp$VisiumDE <- ifelse(tmp$EntrezGeneSymbol %in% seurat_digde, 'Yes', 'No')

LOC_targets <- c('HMGCS1', 'SERPINE1', 'HSPA1B', 'ENO3', 'HSPA1A', 'CDA', 'PYGL', 'IL1RAP', 'SHBG', 'IGFBP2', 'ME1', 'CTSZ', 'DEFB1')
tmp$LOC_targets <- ifelse(tmp$EntrezGeneSymbol %in% LOC_targets, 'LOC', 'No')
write.csv(tmp, file.path(OUTDIR, 'CARDIA-Targets_SpatialTx_DE_wLOC_targets.csv'), row.names = F)

top_DE <- c('IL1RAP', 'SHBG', 'IGFBP2', 'HMGCS1', 'SERPINE1', 'HSPA1B', 'ENO3', 'HSPA1A', 'AKR1C3', 'GHR', 'PCOLCE2', 'AKR1D1')
tmp$top_DE <- ifelse(tmp$EntrezGeneSymbol %in% top_DE, 'LOC', 'No')

p1 <- tmp %>% 
  ggplot(aes(x = Visium_l2fc, y = glmm_coef)) + 
  geom_point() + 
  geom_point(data = subset(tmp, VisiumDE != 'No'),
             aes(x = Visium_l2fc, y = glmm_coef), color = 'purple') + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.65, color = 'red') + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.65, color = 'red') +
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
  labs(x = 'Log2 fold change (fatty/healthy)', 
       y = 'mixed-model coefficient')  +
  ggrepel::geom_text_repel(data = subset(tmp, top_DE != 'No'),
                           aes(x = Visium_l2fc, y = glmm_coef, label = EntrezGeneSymbol), size = 4) 

ggplot2::ggsave(file = file.path(OUTDIR, 'DE_effectsize_comp.pdf'), p1, width = 7, height = 6, units = 'in', dpi = 300)
ggplot2::ggsave(file = file.path(OUTDIR, 'DE_effectsize_comp.svg'), p1, width = 7, height = 6, units = 'in', dpi = 300)
#fc2403

tmp1 <- tmp[which(tmp$EntrezGeneSymbol %in% seurat_digde),]
p2 <- tmp1 %>% 
  ggplot(aes(x = Visium_l2fc, y = glmm_coef)) +
  geom_point() + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.65, color = 'red') + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.65, color = 'red') +
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
  labs(x = 'Log2 fold change (fatty/healthy)', 
       y = 'mixed-model coefficient') + 
  ggrepel::geom_text_repel(data = subset(tmp1, col != 'lightgrey'),
                           aes(x = Visium_l2fc, y = glmm_coef, label = EntrezGeneSymbol), size = 4) 

p1 + p2

cor.test(tmp$Visium_l2fc, tmp$glmm_coef)
cor.test(tmp1$Visium_l2fc, tmp1$glmm_coef)
