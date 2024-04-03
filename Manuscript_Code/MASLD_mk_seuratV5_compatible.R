library(Seurat)
library(SeuratObject)
library(SeuratDisk)

## SET PARAMETERS ----
options(stringsAsFactors = FALSE)
set.seed(1212)
INDIR = "/scratch/nhadad/NAFLD"
if(!dir.exists(file.path(INDIR, 'Output'))){  dir.create(file.path(INDIR, 'Output')) }
OUTDIR = file.path(INDIR, 'Output')
date = Sys.Date()


## Make seurat object V5 compatible ----

### Load seuratObj ----
seuratObj<-readRDS(file=file.path(INDIR, 'seuratObj_humanRNAseq.rds'))
dim(seuratObj@assays$RNA@data)
dim(seuratObj@assays$RNA@counts) ## count data contains ~2k genes more than data so it throws an error; code below should correct that

tmp_data <- seuratObj@assays$RNA@data
tmp_count <- seuratObj@assays$RNA@counts

tmp_data <- GetAssayData(seuratObj, assay = 'RNA')
tmp_count <- seuratObj@assays$RNA@counts

dim(tmp_data)
dim(tmp_count)
genes_tokeep <- intersect(rownames(tmp_count), rownames(tmp_data))
tmp_count <- tmp_count[genes_tokeep,]
dim(tmp_data)
dim(tmp_count)
mean( rownames(tmp_count) == rownames(tmp_data))
mean( colnames(tmp_count) == colnames(tmp_data))
# add new count assay to exisitng object (will replace current one)
DefaultAssay(seuratObj)<-'RNA'
seuratObj <- SetAssayData(seuratObj, slot = 'counts', new.data = tmp_count)

## add meta data to seurat object
seuratObj$Age <- hu_metadata$Age[match(seuratObj$patient, hu_metadata$`Sample Number Flow`)]
seuratObj$Gender <- hu_metadata$Gender[match(seuratObj$patient, hu_metadata$`Sample Number Flow`)]
seuratObj$Weight <- hu_metadata$Weight[match(seuratObj$patient, hu_metadata$`Sample Number Flow`)]
seuratObj$Steatosis <- hu_metadata$Steatosis[match(seuratObj$patient, hu_metadata$`Sample Number Flow`)]
seuratObj$Steatosis_over10pct <- hu_metadata$Steatosis_over10pct[match(seuratObj$patient, hu_metadata$`Sample Number Flow`)]


## scale data for protein of interests
seuratObj <- readRDS(file.path(INDIR, 'seuratObj_humanRNAseq_V5COMPATIBLE.rds'))
protemics_genes <- ravi_proteomics_fdr$ageSexRaceBMI$EntrezGeneSymbol
protemics_genes <- protemics_genes[protemics_genes %in% genes_tokeep]
seuratObj <- ScaleData(seuratObj, feature = protemics_genes)

## calculate expression of significant genes
sig_proteomics_genes <- ravi_proteomics_fdr$ageSexRaceBMI %>%
  filter(der_fdr < 0.05 & val_fdr < 0.05)
sig_proteomics_genes <- sig_proteomics_genes[order(sig_proteomics_genes$der_fdr),]

## check which targets exsist in liver sc data
sig_proteomics_genes_moi <- sig_proteomics_genes[ which(sig_proteomics_genes$EntrezGeneSymbol %in% rownames(seuratObj@assays$RNA$scale.data)), ]
## calculated expression score from genes of interest
goi <- sig_proteomics_genes_moi$EntrezGeneSymbol
seuratObj <- AddModuleScore(seuratObj,
                            features = list(goi),
                            name="Liver_proteomics_enriched")

## save output
saveRDS(seuratObj, file.path(INDIR, 'seuratObj_humanRNAseq_V5COMPATIBLE.rds'))