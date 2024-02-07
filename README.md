# Project README

## Data Processing
- Transcriptome data: CEL-Seq
- Quality control: FastQC and MultiQC
- Gene expression levels: Bowtie2 Aligner
- Normalisation: Conversion to TPM

## Developmental Timepoint Prediction
- Tool: RAPToR
- Validation: Comparison with BLIND tool

## TAI calculaiton
- Performed using MyTAI

## DEGs
- Performed using DESeq2
- Identification criteria: adjusted p-value ≤ 0.05, |log2FoldChange| ≥ 1

## KEGG Pathway Analysis
- Tools: Orthofinder and Pathview
- Goal: Identify homologous protein families and pathway coverage

## WGCNA
- Tool: WGCNA package in R
- Goal: Identify gene modules and perform functional term annotations 

Detailed methods and codes can be found in the corresponding GitHub repository.
