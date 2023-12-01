# Whole Genome Sequencing NGS Data Analysis in the Cloud with WDL/Cromwell and Docker

I wrote this WDL script to perform whole genome sequencing (WGS) variant calling analysis on Google Cloud Platform (GCP) cloud computing resources using GATK (Genome Analysis Toolkit) and other tools in Docker containers. Here's a breakdown of its components and functionalities:

### Overall Functionality:
The workflow fetches sequencing data from NCBI's Sequence Read Archive (SRA) database, processes the data to generate an unaligned BAM file, aligns reads to a reference genome, performs quality control, marks duplicates, recalibrates base quality scores, performs variant calling using GATK HaplotypeCaller, filters variants using a Convolutional Neural Network (CNN), and finally, produces various outputs for downstream analysis.

### Input Files Required to Run:
- Human reference genome files (fasta, index, dict, alt, bwt, sa, amb, ann, pac)
- Contamination sites files (ud, bed, mu)
- SNP database (dbSNP_vcf, dbSNP_vcf_index)
- Known indel sites (VCFs and their indices)
- WGS calling interval list
- Resource VCFs and their indices (e.g., mills, gnomAD)
- Other associated files for alignment and variant calling

### Required Bioinformatics Tools:
- SRA Toolkit for fetching data from the SRA database
- BWA for read alignment
- GATK tools for alignment, variant calling, quality control, etc
- SAMtools, Picard, and bedtools (implicitly used within the GATK tools)

### Outputs Produced:
The workflow generates various outputs, including:
- Unaligned BAM file from SRA data retrieval
- Aligned BAM files 
- Marked duplicates metrics file
- Recalibration reports and files
- CRAM file from the final recalibrated BAM file
- VCF files: raw variants, CNN-filtered variants 
- Validation reports
- Self-sample contamination information 

I then used the resulting variant call files to build machine learning models to detect artifactual variant calls in Python with scikit-learn and LightGBM and compared their performance to GATK's 1D CNN.
