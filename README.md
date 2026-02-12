# Gut microbiome and hyperammonemia in rats

## Overview

This repository contains the full analysis pipeline and supporting data for the manuscript:

**"Hyperammonemia reduces beneficial lactobacilli and bifidobacteria populations disrupting the metabolic balance of the gut microbiome in rats".**

The study explores the impact of chronic HA on the GM in rats and analyzes the taxonomic and likely metabolic changes induced by HA.

## Abstract

> The gut microbiome (GM) plays a critical role in metabolic and neurological health and is implicated in hepatic encephalopathy (HE). Chronic hyperammonemia (HA), a major contributor to cognitive and motor impairment in HE, may influence GM structure and function, yet its specific effects in GM remain unclear. Here, it was investigated how chronic HA alters the GM using a rat model fed an ammonia-enriched diet for four weeks. Fecal microbiota profiles obtained by 16S rRNA gene sequencing revealed marked taxonomic shifts in HA rats, with beta-diversity showing clear separation from controls. Genera within the Lachnospiraceae family and ~Alistipes~ genus were enriched in HA rats, while lactic acid–producing and xylanolytic Firmicutes were reduced. Network analysis identified ~Alistipes~ as a central node in the HA microbiome. Predicted metabolic functions were significantly altered, showing negative associations between HA and pathways related to the pyruvate dehydrogenase complex, sucrose and urea degradation, and 4-aminobutyrate (GABA) degradation. Consistent with these predictions, fecal short-chain fatty acid (SCFA) analysis revealed reduced acetic and butyric acid, alongside increased valeric and isobutyric acid levels. The predicted GABA levels increasement by GM would activate GABA receptors in immune cells and would also contribute to peripheral inflammation and, eventually, neuroinflammation. Together, these findings demonstrate that chronic HA reshapes GM composition, disrupts key metabolic pathways, and alters SCFA profiles, providing mechanistic insight into how HA-associated dysbiosis may contribute to the metabolic, immune, and neurological dysfunction characteristic of HE. 


## Repository Structure

```
Rats_HE_HiperAmmon_model/
├── R/                              # R Markdown scripts for each analysis module
│   ├── 1_Dada2_pipeline_rats_hiperamm.Rmd             # Raw sequence processing and ASV generation
│   ├── 2_phyloseq_obj_hiperammon.Rmd                  # Phyloseq object creation
│   ├── 3_decontam_hiperammon.Rmd                      # Contaminant filtering
│   ├── 4_raref_anal_and_filter_hiperammon.Rmd         # Rarefaction and filtering
│   ├── 5_Alpha_DIv_hiperammon.Rmd                     # Alpha diversity analysis
│   ├── 6_Batch_eff_remv_and_BetaDiv_hipperammon.Rmd   # Beta diversity and statistical modeling
│   ├── 7_metadeconfoundR.Rmd                          # Metadata confounding analysis
│   ├── 8_Picrust2_analysis_modules.Rmd                # Picrust2
│   ├── 9_Network_analysis.Rmd                         # Microbial network analysis
│   └── ...                                     # Additional scripts for figures, tables, and reviewer responses
├── README.md                       # Project description and usage instructions
```

All scripts are written in R and require packages such as phyloseq, decontam, ggplot2, metadeconfoundR, and picrust2.
Run analysis: Follow the numbered R Markdown files in the R/ folder to reproduce the pipeline step-by-step.
Figures and tables: Scripts generate publication-ready plots and tables.

