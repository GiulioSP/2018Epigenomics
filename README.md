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

## Scripts:
- [._count_analysis_P6.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/._count_analysis_P6.Rmd):
- [._count_analysis_tum.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/._count_analysis_tum.Rmd):
- [._cpgToDmrByRegion.3files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/._cpgToDmrByRegion.3files.sh):
- [._cpgToDmrByRegion.6files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/._cpgToDmrByRegion.6files.sh):
- [._finderToChipByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/._finderToChipByRegion.sh):
- [._vcfToVcfByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/._vcfToVcfByRegion.sh):
- [cgiToCgiByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cgiToCgiByRegion.sh):
- [count_analysis_P6.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/count_analysis_P6.Rmd): R code that aggregates different types of data by region, then plots comparisons and exports results to a .csv file. Set up for P6 (patient 6 in my analysis, AKCC55). Tweak file names for a different patient's datasets. Used to gleam insights into cancer evolution, by noting epigenetic and genetic factors that are different between normal, tumor adjacent and tumor samples.
- [count_analysis_tum.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/count_analysis_tum.Rmd): R code that aggregates all data from tumor samples by region, then plots comparisons and exports results to a .csv file. Used to gleam insights into common epigenetic mechanisms in cancer, often associated with deactivation of tumor suppressor genes.
- [cpgToDmrByRegion.3files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cpgToDmrByRegion.3files.sh): Creates .dmrByRegion file, which lists differentialy methylated regions across 3 files, parsing the genome into same sized chunks ("size"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: "chr:positionStart:positionEnd", for each file "adjustedResidual/max(adjustedResidual)", and "chiSquare". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: ".5mC.CpG.gz". Will create ".methyByRegion.txt" files for each input file, as well as a ".methyByRegion.join.txt". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the "chi" vector. 
- [cpgToDmrByRegion.6files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cpgToDmrByRegion.6files.sh): Creates .dmrByRegion file, which lists differentialy methylated regions across 6 files, parsing the genome into same sized chunks ("size"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: "chr:positionStart:positionEnd", for each file "adjustedResidual/max(adjustedResidual)", and "chiSquare". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: ".5mC.CpG.gz". Will create ".methyByRegion.txt" files for each input file, as well as a ".methyByRegion.join.txt". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the "chi" vector. 
- [finderToChipByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/finderToChipByRegion.sh):
- [vcfContrast.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/vcfContrast.sh): Contrasts VCF files to isolate variants unique to the main sample (the first in fileList). Can apply filters for higher specificity, at the cost of sensitivity, when -specif is added.
- [vcfToVcfByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/vcfToVcfByRegion.sh): Creates vcfByRegion file, which parses genome into same sized chunks ("size"), and reports the chunks where one or more variant calls are present, in the format: "chr:positionStart:positionEnd",  "adjustedResidual/max(adjustedResidual)", and "chiSquare". Requires enrichment BED file: FindER.bed.gz.3col.sorted.


## Analyses:
- countAnalysis.plots.xlsx 
- Lab Meeting - Aug 13th.pptx
- Lab Meeting - Aug 22.pptx
- Lab Meeting - July 4th.pptx


## Glossary:
- hg38 -> reference sequence for the human genome, build version 38
- 
