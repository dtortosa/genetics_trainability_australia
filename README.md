# Genome-Wide Association Study:  Genetic Predictors of Cardiorespiratory Fitness Response to Military Training

## Overview

This project presents a large-scale genomic analysis investigating genetic factors associated with cardiorespiratory fitness (CRF) adaptation in Australian Army recruits undergoing Basic Military Training.  As the data scientist on this project, I designed and implemented a comprehensive bioinformatics pipeline handling genomic data from 1,422 participants, performing rigorous quality control, and developing predictive polygenic score models.

## Project Highlights

### Large-Scale Data Management
- **1,422 participants** with complete genomic and phenotypic data
- **~540,000 genetic markers** initially genotyped (Infinium CoreExome-24 Kit)
- **~3.2 million genetic variants** after quality control and imputation
- Multiple phenotypes tracked: VO2max changes, Multistage Fitness Test (MSFT) scores, and total distance covered
- Longitudinal fitness assessments (Week 2 and Week 9 of training)

### Technical Stack

**Programming & Data Processing:**
- Genomic data format conversion and management
- Large-scale statistical computing
- Iterative quality control workflows
- Batch processing and pipeline automation

**Specific languages and tools:**
- Python, Bash and AWK
- Pyspark and Multiprocessing for distributed and parallel computing
- Singularity containers
- Cloud computing (University of Arizona's High Performance Computing servers)

### Role & Contributions

As the **data scientist** on this project, I was responsible for:
- Designing and implementing the complete quality control pipeline
- Managing and processing large-scale genomic datasets
- Conducting population structure and relatedness analyses
- Developing and testing multiple polygenic score models
- Implementing robust cross-validation strategies
- Generating comprehensive visualizations (Manhattan plots, PGS correlation plots, forest plots)
- Statistical analysis and interpretation of results

### Rigorous Quality Control Pipeline

Implemented a multi-stage quality control workflow following GWAS best practices:

#### **Variant-Level Filtering**
- Genotyping rate threshold (>0.99)
- Minor allele frequency filtering (MAF >0.05)
- Hardy-Weinberg Equilibrium testing (p < 1e-25)
- Imputation quality control (INFO/R2 > 0.95)

#### **Sample-Level Quality Control**
- Missing genotype rate filtering (<1% missing)
- LogR value validation (Illumina's threshold:  <0.3)
- **Relatedness detection**:  Implemented dual-method approach using both IBD (Identity by Descent) and KING algorithms to identify and remove related individuals (35 samples removed)
- **Population structure analysis**:  Iterative principal component analysis (PCA) using SMARTPCA/EIGENSOFT to detect and remove ancestry outliers (171 samples removed)
- Sex concordance verification between self-reported and genotypic data (37 discordant samples removed)

#### **Iterative Filtering Strategy**
- Sequential loops of variant and sample filtering to ensure no threshold violations
- Multiple rounds of MAF and missingness filters after each sample removal stage
- Final Hardy-Weinberg Equilibrium filter post-sample QC to detect variant calling errors

#### **Batch Effect Assessment**
- Verified no systematic differences between genotyping batches
- Average MAF similarity across batches (0.245 vs 0.246)
- Tested association between SNPs and batch membership (average p-value: 0.497)

### Advanced Imputation
- Utilized **TOPMed imputation server** for high-quality genotype imputation
- Prepared data with allele flip resolution and strand alignment
- Achieved high correlation between study data and imputation reference panel
- Stringent post-imputation filtering (INFO/R2 > 0.95)

### Sophisticated Statistical Modeling

#### **Single Variant Analysis**
- Genome-wide association testing for all ~3.2 million variants
- Linear models adjusting for relevant covariates (age, sex, weight, baseline fitness, 20 genetic PCs)
- Manhattan plots generated for visualization of association signals
- Bonferroni correction applied for multiple testing

#### **Polygenic Score (PGS) Development**
Implemented multiple PGS calculation strategies using LDAK software:

**Classical Thresholding + Clumping Approach:**
- Tested across 6 p-value thresholds (1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
- LD-based clumping to handle correlated variants

**Advanced Modeling:**
- **Elastic net approach** accounting for MAF-dependent heritability contributions
- Superior predictive performance compared to classical methods

**Robust Validation Framework:**
- **100 iterations** of 75/25 train-test splits
- Calculated 95% confidence intervals across iterations for correlation between PGS and phenotypes
- Generated **4,200 total polygenic scores** across phenotypes, models, covariate selection approaches, and data partitions
- Tested both reduced (phenotype-associated) and full covariate sets

## Key Findings

- Successfully processed and cleaned genomic data from initial 1,422 samples to final analysis dataset of 1,204 high-quality samples
- No significant batch effects detected in final dataset
- Elastic net and unrestricted classical PGS models showed highest predictive performance
- Median correlations between PGS and phenotypes ranged from 0.056 to 0.085
- Clear trend of increased fitness gains in higher PGS quartiles

## Data Quality Metrics

- **Final sample retention rate**: 84.7% (1,204/1,422)
- **Final variant count**: ~3.2 million high-quality variants
- **Genotype call rate**: >99%
- **Imputation quality**: INFO/R2 > 0.95
- **Population structure**: No ancestry outliers remaining (confirmed via PCA)
- **Relatedness**: No second-degree or closer relatives in final dataset

## Methodology Strengths

✓ **Comprehensive QC**:  Multi-stage quality control following established GWAS protocols  
✓ **Large Sample Size**:  Exceeds recommended threshold (>1,000) for complex trait association studies  
✓ **Robust Validation**: 100-fold train-test splits prevent overfitting  
✓ **Multiple Approaches**: Compared classical and state-of-the-art PGS methods  
✓ **Proper Covariate Adjustment**: Controlled for population structure, demographics, and baseline fitness  
✓ **Reproducible Pipeline**: Followed best practices from genomic analysis literature  

## References

Quality control and analysis procedures followed established guidelines: 
- Choi et al. (2020) - Polygenic risk score analysis tutorial
- Marees et al. (2018) - GWAS quality control and statistical analysis
- Truong et al. (2022) - Quality control procedures for GWAS
- Berrandou et al. (2023) - LDAK methodology for PGS calculation

---
