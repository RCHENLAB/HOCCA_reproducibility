# Human Ocular Surface Cell Atlas

This repository contains code and analyses for the Human Ocular Surface Cell Atlas project within the HCA framework, with an emphasis on the [Eye Biological Network](https://www.humancellatlas.org/biological-networks/). Here, the ocular surface includes the cornea, limbus, and anterior sclera. 

## System requirements
The pipeline is expected to run on standard Linux-based high-performance computing environments with R and Python installed. The full workflow was tested on Linux Rocky Linux 8.8 and above. No required non-standard hardware
- R version: R 4.3.2 
- Python version: Python v3.10 and above

## Software and code
The code used for atlas construction follows the established pipeline developed for the Human Retina Cell Atlas (HRCA) and is publicly available through [HRCA_reproducibility](https://github.com/RCHENLAB/HRCA_reproducibility). The analyses listed below were performed using the pipelines in the [HRCA_reproducibility](https://github.com/RCHENLAB/HRCA_reproducibility) with the parameters specified in the manuscript:
  - Integration, clustering, and annotation of single-cell and single-nuclei RNA-sequencing data (scRNA-seq and snRNA-seq)
  - Quality control and annotation of single-nuclei ATAC-sequencing (snACTA-seq) cells
  - Co-embedding of snATAC-seq and scRNA-seq cells
  - Identification of regulon from snATAC-seq data
  - Differential gene expression and pathway enrichment during aging
  - GWAS trait enrichment and fine-mapping of associated variants

For Data preprocessing and quality control:
  - Preprocessing: [10x Cell Ranger v7](https://www.10xgenomics.com/support/software/cell-ranger/7.2)
  - Quality control: [CellQC](https://github.com/lijinbio/cellqc)
    - Ambient RNA contamination correction: SoupX
    - Low-quality droplets identification: dropkick
    - Putative doulbets detected and removal: DoulbetFinder

The codes specifically for the downstream analysis for this atlas can be found in the following aspects:
- [Epithelial cell trajectory and module identification](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Epithelial_trajectory):
  - [RNA velocity and trajectory inference](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Epithelial_trajectory/RNA_velocity_trajectory_inference)
  - [Gene module analysis](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Epithelial_trajectory/Gene_module_analysis)
- [Spatial data analysis](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Spatial_analysis):
  - Spatial data preprocessing: [Xenium Ranger v3](https://www.10xgenomics.com/support/software/xenium-ranger/latest)
  - [Preprocessing, integration, and annotation](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Spatial_analysis/preprocessing_annotation)
  - [Neighborhood enrichment and spatial niche analysis](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Spatial_analysis/Neighborhood_niche_analysis)
  - [Cornea-limbus axis analysis](https://github.com/RCHENLAB/HOCCA_reproducibility/tree/main/Spatial_analysis/corneal_limbal_axis_analysis), including:
      - [Distance between Melanocytes and epithelial cells with different states](https://github.com/RCHENLAB/HOCCA_reproducibility/blob/main/Spatial_analysis/corneal_limbal_axis_analysis/cal_dist_to_melanocytes.R)
      - [Cell proportion of the epithelial cells in different states across the cornea-limbus axis](https://github.com/RCHENLAB/HOCCA_reproducibility/blob/main/Spatial_analysis/corneal_limbal_axis_analysis/epi_cell_proportion_across_bins.R)
 
### Versions of the software
CellRanger: v7.0.0, v7.0.1, v7.1.0, and v7.2.0; 
SoupX: v1.6.2;
dropkick: v1.2.8;
DoubletFinder: v2.0.4;
scvi-tools: v1.2.0;
NSForest: v4.0;
Scanpy: v1.10.3;
Velocyto: v0.17.17;
scVelo: v0.3.3;
escape: v2.7.3;
hotspotsc: v1.1.3;
hypeR: v2.10.0;
edgeR: v4.4.0;
variancePartition: v1.36.2;
CellChat: v2.1;
Enrichr: v3.4;
cellranger-atac-2.0.0;
ArchR package: v1.0.3;
scGLUE: v0.3.2;
SCENIC+ package: v1.0a2;
Xenium Ranger: v3.1.0
Seurat: v5.0.0 for scRNA/snRNA-related analyses and v5.3.0 for spatial-related analyses;
Harmony: v1.2.3;
spacexr: v2.2.1;
Squidpy:  v1.6.2;
MAGMA.Celltype:v2.0.8;
LDSC: v1.0.1;
EWCE: v1.16.0;
susieR: v0.12.45.

## Interactive Browsers

The atlas can be accessed through [CELLxGENE](https://cellxgene.cziscience.com/collections/0f7d022a-46c7-4e64-be4c-e34adbb78089) and [CAP](https://celltype.info/project/565).

## Questions

If you have any questions, please submit an issue to this repository.
