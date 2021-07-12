.libPaths( "/home/icb/anna.danese/miniconda3/envs/scmoib-seuratv4/lib/R/library")
library(Seurat)
library(SeuratDisk)
library(Signac)

# Load the peak ATAC
file = "/home/icb/anna.danese/project_anna/scmoib/brain_peaks_filtered_not_normalised.h5ad"
file_seurat = "/home/icb/anna.danese/project_anna/scmoib/brain_peaks_filtered_not_normalised.h5seurat"
Convert(file, dest = "h5seurat", overwrite = TRUE, assay='ATAC')
atac <- LoadH5Seurat(file_seurat)
brain.atac <- CreateSeuratObject(counts = atac[['ATAC']], assay = "ATAC", project = "SHAREseq_ATAC")


# RNA
file = "/home/icb/anna.danese/project_anna/scmoib/processed_data/brain_data/brain_preprocessed_rna_full_features.h5ad"
file_seurat = "/home/icb/anna.danese/project_anna/scmoib/processed_data/brain_data/brain_preprocessed_rna_full_features.h5seurat"
Convert(file, dest = "h5seurat", overwrite = TRUE, assay='RNA')
rna <- LoadH5Seurat(file_seurat)
brain.rna <- CreateSeuratObject(counts = rna[['RNA']], assay = "RNA", project = "SHAREseq_RNA")


# Perform standard analysis of each modality independently RNA analysis
brain.rna <- NormalizeData(brain.rna)
brain.rna <- FindVariableFeatures(brain.rna)
brain.rna <- ScaleData(brain.rna)
brain.rna <- RunPCA(brain.rna)
brain.rna <- RunUMAP(brain.rna, dims = 1:30)

#create the raw geneactivity matrix
gtf_file = "/home/icb/anna.danese/project_anna/scmoib/processed_data/gencode.vM26.annotation.gtf"
chromosome <- paste0('chr', c(1:19, "X", "Y"))
chromosome
activity.matrix <- GeneActivity(peak.matrix = atac[['ATAC']]@counts, annotation.file = gtf_file, 
                                            seq.levels = chromosome, upstream = 2000, verbose = TRUE)
