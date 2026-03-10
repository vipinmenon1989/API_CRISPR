# API_CRISPR

An automated pipeline for collecting and curating **CRISPR screen datasets** (Cas12 and Cas13) from scientific literature using AI-powered data extraction.

## Overview

This repository provides tools to automatically search PubMed for CRISPR screen publications, retrieve full-text articles, and use Google Gemini AI to extract structured experimental metadata. The pipeline supports both Cas12 and Cas13 CRISPR systems and is optimized for high-performance computing (HPC) environments.

## Features

- **Automated PubMed Search**: Queries NCBI Entrez API with configurable search terms targeting CRISPR screen publications in high-impact journals
- **AI-Powered Data Extraction**: Uses Google Gemini 2.0 Flash to extract structured metadata (lab, method, sample size N) from paper text
- **Multi-Strategy Text Retrieval**: Attempts full-text retrieval via PMC, DOI scraping, and Playwright browser automation as fallbacks
- **Parallel Processing**: Uses Python ProcessPoolExecutor for concurrent paper processing
- **HPC-Ready**: Includes SLURM job scripts for cluster execution
- **Dual CRISPR Support**: Separate pipelines for Cas12 and Cas13 datasets

## Repository Structure

```
API_CRISPR/
├── human_api_approach.py          # Main pipeline: PubMed search + Gemini AI extraction
├── human_api.slurm                # SLURM job script for HPC cluster
├── API_Cas12_Cas13.yml            # Conda environment specification
├── python_cas12_cas13_pure_api/   # Extended pipeline with Playwright support
│   ├── dataloader.py              # Cas12 dataloader (with browser scraping fallback)
│   ├── dataloader_cas13.py        # Cas13 dataloader
│   ├── recovery_cas13.py          # Recovery script for failed Cas13 downloads
│   ├── dataloader.slurm           # SLURM script for Cas12 dataloader
│   ├── dataloader_cas13.slurm     # SLURM script for Cas13 dataloader
│   └── recovery_cas13.slurm       # SLURM script for recovery
```

## Prerequisites

### Environment Setup

```bash
conda env create -f API_Cas12_Cas13.yml
conda activate CRISPR
```

### API Keys

Create a `.env` file in the project root:

```env
GEMINI_API_KEY=your_google_gemini_api_key
ENTREZ_API_KEY=your_ncbi_entrez_api_key
ENTREZ_EMAIL=your_email@example.com
UMD_COOKIES=your_umd_cookie_string
```

## Usage

### Main Pipeline

```bash
python human_api_approach.py
```

Steps performed:
1. Search PubMed using configured queries and forced PMIDs
2. Retrieve paper metadata (journal, title) via Entrez
3. Fetch full text (PMC first, then DOI scraping, then abstract fallback)
4. Use Gemini 2.0 Flash AI to extract lab, method, and sample size (N)
5. Save results to big_boy_wei_li_final.csv

### HPC Cluster (SLURM)

```bash
sbatch human_api.slurm
```

### Extended Cas12/Cas13 Pipeline

```bash
sbatch python_cas12_cas13_pure_api/dataloader.slurm
sbatch python_cas12_cas13_pure_api/dataloader_cas13.slurm
```

## Output

The pipeline generates a CSV with these columns:

| Column  | Description |
|---------|-------------|
| PMID    | PubMed ID |
| Journal | Journal name |
| Lab     | PI/Lab name (AI-extracted) |
| Method  | CRISPR screen method |
| N       | Sample size / number of guides |
| URL     | PubMed link |

## Key Configuration Parameters

In `human_api_approach.py`:

- `FORCE_PMIDS`: Specific PMIDs to always include (e.g., DeepCas12a, DeepCas13 benchmark papers)
- `SEARCH_QUERIES`: PubMed search queries for literature discovery
- `ROYAL_JOURNALS`: Target journals (Nature, Science, Cell, Genome Biology, etc.)
- `max_workers=5`: Parallelism — reduce if hitting API rate limits

## Dependencies

Key Python packages (see API_Cas12_Cas13.yml for full list):

- biopython: NCBI Entrez API
- google-generativeai: Google Gemini AI
- requests: HTTP requests and DOI scraping
- playwright: Browser automation (extended pipeline)
- tenacity: Retry logic with exponential backoff
- pandas: Data output

## Author

**Vipin Menon**
PhD Candidate, Computational Biology (Genome Editing)
Hanyang University, Seoul, South Korea
