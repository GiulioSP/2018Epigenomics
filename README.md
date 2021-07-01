### 2018Epigenomics
Late upload of work done in 2018 for the Hirst lab at UBC.

## Data types:
# Genomic data
- WGS: Whole genome sequencing raw data.
- VCF: Variant Call Format. Lists all variations found between a WGS file and a reference human genome, like hg38 or hg19
- 

# Methylomic data 
- WGBS: Whole genome bisulfite sequencing. The bisulfite treatment modifies Cytosines that were unmethylated, so one can extrapolate methylation data from WGBS by comparing with WGS data from the same sample.
- Compare samples to identify differentially methylated regions (DMRs)

# Hystone Modification data
- ChIP-seq 
- Compare samples to find differential histone modifications (DHMs)


# Gene Expression data
- RNA-seq and miRNA-seq: 
- Compare samples to find differentially expressed regions (DEs)


## Data source:
Canadian Epigenomes dataset

http://www.epigenomes.ca/site-data
In order to download data, select desired datasets and click "Download Tracks", then open each relevant track as a tab with your browser, right click and select "Save Page As ...".

Colorectal samples:
- Tumor
- Tumor Adjacent
- Normal

![image](https://github.com/GiulioSP/2018Epigenomics/blob/main/new_files/data_breakdown.png)

Grouping samples by type and patient:

| Patient code | Tumor sample | Tumor adjacent sample | Normal sample |
| :----: | :----: | :----: | :----: |
| AKCC46 | CEMT0062 | CEMT0052 | CEMT0034 | 
| AKCC52 | CEMT0063 | CEMT0053 | CEMT0033 | 
| AKCC58 | CEMT0064 | CEMT0055 | CEMT0054 | 
| AKCC63 | CEMT0065 | CEMT0057 | CEMT0056 | 
| AKCC70 | CEMT0066 | CEMT0059 | CEMT0058 | 
| AKCC55 | CEMT0067 | CEMT0061 | CEMT0060 | 

Datasets from patient AKCC46 are saved in repo as example data.


## Glossary:
- hg38 -> reference sequence for the human genome, build version 38
- 
