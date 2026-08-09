# Bioinformatics Pipelines & Analysis Scripts

Reproducible Nextflow pipelines and R analysis scripts for bulk and single-cell 
RNA-seq, spatial transcriptomics, and related omics data, built for execution 
on HPC/SLURM clusters. Developed as part of cancer research at Cedars-Sinai 
Medical Center.

## Repository Structure

### `nextflow/`
Modular DSL2 Nextflow pipelines. Each processing step is a separate module 
under `nextflow/modules/`, imported into one of two top-level workflows under 
`nextflow/workflows/`.

#### `workflows/rnaseq.nf` — Bulk RNA-seq pipeline
End-to-end bulk RNA-seq analysis, from raw FASTQs to differential expression, 
pathway analysis, and transcription factor activity:

1. FASTQ renaming, MD5 manifest generation, and input validation
2. Reference retrieval (Ensembl) and indexing (STAR, Salmon), with xenograft 
   support (human/mouse read deconvolution via Xengsort)
3. Quantification via Salmon (alignment-free) and STAR (splice-aware alignment)
4. Alignment QC (RSeQC) and aggregated reporting (MultiQC)
5. Differential expression (DESeq2), with QC plots and result tables
6. Pathway enrichment analysis and transcription factor / regulon activity 
   analysis, each with accompanying visualization
7. A final "deep audit" step that cross-validates alignment/quantification 
   stats against the DESeq2 objects for QC assurance

Supports flexible run modes (`bulk`, `de`, `bulk+de`) and can resume from any 
completed step via Nextflow's caching.

#### `workflows/scrnaseq.nf` — Single-cell RNA-seq pipeline
FASTQ validation and QC through CellRanger quantification and aggregated 
MultiQC reporting, structured to feed into downstream single-cell analysis 
(clustering, integration, cell-type annotation — see `R/scRNASeq_workflow.R`).

#### `modules/`
Individual, reusable process definitions called by both workflows — reference 
handling, alignment/quantification (STAR, Salmon, CellRanger), QC (FastQC, 
RSeQC, MultiQC), and downstream statistics (DESeq2, pathway analysis, 
transcription factor analysis).

**Output structure** (both pipelines) is organized into numbered stage 
directories (e.g., `01.FastQ/`, `02.FastQC/`, `03.Salmon/`, `04.STAR/`, 
`08.DESeq2/`, `09.Pathways/`) with a final `MultiQC` HTML report and full 
execution logs/trace files for reproducibility.

### `R/`
Wrapper and utility scripts for analysis and visualization across omics types:
- `RNASeq_workflow.R`, `scRNASeq_workflow.R` — bulk/single-cell analysis workflows
- `GeoMx_Wrapper.R`, `Spatial_Wrapper.R`, `Xenium_Wrapper.R` — spatial transcriptomics
- `Microarray_Wrapper.R`, `Metabolomics_wrapper.R` — additional omics support
- `Custom_Functions.R` — shared helper functions
- `access_GEO.R` — GEO database retrieval
- `compare_dds.R` — comparison utilities for DESeq2 objects

## Tools & Technologies
Nextflow (DSL2), R/Bioconductor (DESeq2, tximport), STAR, Salmon, CellRanger, 
RSeQC, MultiQC, SLURM/HPC, Singularity

## Related Publication
Analysis approaches in this repository supported:
Mani, S.K.K., et al. "Single-Cell Profiling of Murine Bladder Cancer Identifies 
Sex-Specific Transcriptional Programs with Prognostic Significance" *iScience*, 
2023. DOI: [10.1016/j.isci.2023.107703](https://doi.org/10.1016/j.isci.2023.107703)

## Contact
Saravana Kumar Kailasam Mani
[LinkedIn](https://linkedin.com/in/skailasammani) | skailasammani@gmail.com

