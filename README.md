### 2018Epigenomics
Late upload of work done in 2018 for the Hirst lab at UBC.

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

## Data types:
# Genomic data
- WGS: Whole genome sequencing raw data.
- VCF: Variant Call Format. Lists all variations found between a WGS file and a reference human genome, like hg38 or hg19. Derived from WGS data with BCFtools package.
- VcfByRegion: Original data format that lists variant calls and their attributes separated by regions, for easy comparison with epigenetic data sets that are otherwise impossible to synthesize. Generates with script vcfToVcfByRegion.sh.

# Methylomic data 
- WGBS: Whole genome bisulfite sequencing. The bisulfite treatment modifies Cytosines that were unmethylated, so one can extrapolate methylation data from WGBS by comparing with WGS data from the same sample.
- CPG: Pre-processed WGBS data that lists all methylated Cytosines along with position. 
- dmrByRegion: Original data format that splits the genome into small regions and counts the amount of methylated Cytosines in each region. Then, compares between samples to identify differentially methylated regions (DMRs) according to a Chi Squared test. Done with scripts cpgToDmrByRegion.3files.sh and cpgToDmrByRegion.6files.sh.  

# Hystone Modification data
- ChIP-seq: Sequencing data after a step of chromatin immunoprecipitation (ChIP) that should repress genetic signal of regular DNA and amplify signal near a targeted protein. In this case, the targeted proteins are a set of common histone modifications. 
- FindER: Absolute measurement of histone modification enrichment detected. Generated from ChIP-Seq data by the CEEHRC's FindER package.  
- chipByRegion: Original data format that splits the genome into regions and aggregates histone modifications associated with each region. Generated with script finderToChipByRegion.


## Scripts:

# Analyses in R
- [count_analysis_P6.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/count_analysis_P6.Rmd): R code that aggregates different types of data by region, then plots comparisons and exports results to a .csv file. Set up for P6 (patient 6 in my analysis, AKCC55). Tweak file names for a different patient's datasets. Used to gleam insights into cancer evolution, by noting epigenetic and genetic factors that are different between normal, tumor adjacent and tumor samples.
- [count_analysis_tum.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/count_analysis_tum.Rmd): R code that aggregates all data from tumor samples by region, then plots comparisons and exports results to a .csv file. Used to gleam insights into common epigenetic mechanisms in cancer, often associated with deactivation of tumor suppressor genes.

# Preprocessing Data with Bash
- [cgiToCgiByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cgiToCgiByRegion.sh): Creates .cgiByRegion file, which parses genome into same sized chunks ("size"), and reports regions where most bps are part of a CGI. Each region is reported as: "chr:positionStart:positionEnd" and "1". Requires CGI position BED file, in this case "hg38_CGIs".
- [cpgToDmrByRegion.3files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cpgToDmrByRegion.3files.sh): Creates .dmrByRegion file, which lists differentialy methylated regions across 3 files, parsing the genome into same sized chunks ("size"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: "chr:positionStart:positionEnd", for each file "adjustedResidual/max(adjustedResidual)", and "chiSquare". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: ".5mC.CpG.gz". Will create ".methyByRegion.txt" files for each input file, as well as a ".methyByRegion.join.txt". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the "chi" vector. 
- [cpgToDmrByRegion.6files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cpgToDmrByRegion.6files.sh): Creates .dmrByRegion file, which lists differentialy methylated regions across 6 files, parsing the genome into same sized chunks ("size"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: "chr:positionStart:positionEnd", for each file "adjustedResidual/max(adjustedResidual)", and "chiSquare". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: ".5mC.CpG.gz". Will create ".methyByRegion.txt" files for each input file, as well as a ".methyByRegion.join.txt". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the "chi" vector. 
- [finderToChipByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/finderToChipByRegion.sh): Creates chipByRegion file, which parses genome into same sized chunks ("size"), and reports the chunks where most bps are considered enriched by this histone modification, in the format "chr:positionStart:positionEnd" and "1". Requires enrichment BED file: "FindER.bed.gz.3col.sorted". The postfix will be ".chipByRegion".
- [vcfContrast.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/vcfContrast.sh): Contrasts VCF files to isolate variants unique to the main sample (the first in fileList). Can apply filters for higher specificity, at the cost of sensitivity, when -specif is added.
- [vcfToVcfByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/vcfToVcfByRegion.sh): Creates vcfByRegion file, which parses genome into same sized chunks ("size"), and reports the chunks where one or more variant calls are present, in the format: "chr:positionStart:positionEnd",  "adjustedResidual/max(adjustedResidual)", and "chiSquare". Requires enrichment BED file: FindER.bed.gz.3col.sorted.

## External Packages
- [FindER](http://www.epigenomes.ca/tools-and-software/finder) processes ChIP-Seq data, identifying enrichment from histone modifications.  
- [BCFtools](https://samtools.github.io/bcftools/bcftools.html) processes WGS data into a list of variant calls (VCF). It also processes VCF files in specific ways, such as contrasting two VCF files to list statistically significant differences in genetic variants.


## Analyses:
- countAnalysis.plots.xlsx 
- Lab Meeting - Aug 13th.pptx
- Lab Meeting - Aug 22.pptx
- Lab Meeting - July 4th.pptx


## Glossary:
- hg38 -> reference sequence for the human genome, build version 38
- 
