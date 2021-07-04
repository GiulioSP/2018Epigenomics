# 2018Epigenomics
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
### Genomic data
- WGS: Whole genome sequencing raw data.
- VCF: Variant Call Format. Lists all variations found between a WGS file and a reference human genome, like hg38 or hg19. Derived from WGS data with BCFtools package.
- VcfByRegion: Original data format that lists variant calls and their attributes separated by regions, for easy comparison with epigenetic data sets that are otherwise impossible to synthesize. Generates with script vcfToVcfByRegion.sh.

### Methylomic data 
- WGBS: Whole genome bisulfite sequencing. The bisulfite treatment modifies Cytosines that were unmethylated, so one can extrapolate methylation data from WGBS by comparing with WGS data from the same sample.
- CPG: Pre-processed WGBS data that lists all methylated Cytosines along with position. 
- dmrByRegion: Original data format that splits the genome into small regions and counts the amount of methylated Cytosines in each region. Then, compares between samples to identify differentially methylated regions (DMRs) according to a Chi Squared test. Done with scripts cpgToDmrByRegion.3files.sh and cpgToDmrByRegion.6files.sh.  

### GpG Island data
- CGI: A list of regions considered to be CpG islands in a given genome build (in this case, hg38). These regions are rich in the sequence "CG", which is the most common sequence for a methylated Cytosine; they are usually associated with the upstream of genes and serve a regulatory function, whereas fully methylated CGIs suppress expression (or at least transcription factor binding), and unmethylated CGIs promote expression (or at least allow transcription factor binding).    
- cgiByREgion: Original data format that splits the genome into small regions and notes if a region is considered to be a CGI.

### Histone Modification data
- ChIP-seq: Sequencing data after a step of chromatin immunoprecipitation (ChIP) that should repress genetic signal of regular DNA and amplify signal near a targeted protein. In this case, the targeted proteins are a set of common histone modifications. 
- FindER: Absolute measurement of histone modification enrichment detected. Generated from ChIP-Seq data by the CEEHRC's FindER package.  
- chipByRegion: Original data format that splits the genome into regions and aggregates histone modifications associated with each region. Generated with script finderToChipByRegion.

## Scripts:
### Analyses in R
- [count_analysis_P6.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/count_analysis_P6.Rmd): R code that aggregates different types of data by region, then plots comparisons and exports results to a .csv file. Set up for P6 (patient 6 in my analysis, AKCC55). Tweak file names for a different patient's datasets. Used to gleam insights into cancer evolution, by noting epigenetic and genetic factors that are different between normal, tumor adjacent and tumor samples.
- [count_analysis_tum.Rmd](https://github.com/GiulioSP/2018Epigenomics/blob/main/count_analysis_tum.Rmd): R code that aggregates all data from tumor samples by region, then plots comparisons and exports results to a .csv file. Used to gleam insights into common epigenetic mechanisms in cancer, often associated with deactivation of tumor suppressor genes.

### Preprocessing Data with Bash
- [cgiToCgiByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cgiToCgiByRegion.sh): Creates .cgiByRegion file, which parses genome into same sized chunks ("size"), and reports regions where most bps are part of a CGI. Each region is reported as: "chr:positionStart:positionEnd" and "1". Requires CGI position BED file, in this case "hg38_CGIs".
- [cpgToDmrByRegion.3files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cpgToDmrByRegion.3files.sh): Creates .dmrByRegion file, which lists differentialy methylated regions across 3 files, parsing the genome into same sized chunks ("size"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: "chr:positionStart:positionEnd", for each file "adjustedResidual/max(adjustedResidual)", and "chiSquare". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: ".5mC.CpG.gz". Will create ".methyByRegion.txt" files for each input file, as well as a ".methyByRegion.join.txt". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the "chi" vector. 
- [cpgToDmrByRegion.6files.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/cpgToDmrByRegion.6files.sh): Creates .dmrByRegion file, which lists differentialy methylated regions across 6 files, parsing the genome into same sized chunks ("size"), and reports regions where the methylated CpG count of the files inputed is sufficiently different to pass a chi squared test. Each differentially methylated region is reported as: "chr:positionStart:positionEnd", for each file "adjustedResidual/max(adjustedResidual)", and "chiSquare". Where adjustedResidual=(obs–exp)*total/√[exp*(total-rowTotal)*(total-columnTotal)]. Requires methylation counts BED file: ".5mC.CpG.gz". Will create ".methyByRegion.txt" files for each input file, as well as a ".methyByRegion.join.txt". This script can compare at most 28 files at once; to increase this number, edit the code to add more chi squared values to the "chi" vector. 
- [finderToChipByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/finderToChipByRegion.sh): Creates chipByRegion file, which parses genome into same sized chunks ("size"), and reports the chunks where most bps are considered enriched by this histone modification, in the format "chr:positionStart:positionEnd" and "1". Requires enrichment BED file: "FindER.bed.gz.3col.sorted". The postfix will be ".chipByRegion".
- [vcfContrast.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/vcfContrast.sh): Contrasts VCF files to isolate variants unique to the main sample (the first in fileList). Can apply filters for higher specificity, at the cost of sensitivity, when -specif is added.
- [vcfToVcfByRegion.sh](https://github.com/GiulioSP/2018Epigenomics/blob/main/vcfToVcfByRegion.sh): Creates vcfByRegion file, which parses genome into same sized chunks ("size"), and reports the chunks where one or more variant calls are present, in the format: "chr:positionStart:positionEnd",  "adjustedResidual/max(adjustedResidual)", and "chiSquare". Requires enrichment BED file: FindER.bed.gz.3col.sorted.

## External Packages
- [FindER](http://www.epigenomes.ca/tools-and-software/finder) processes ChIP-Seq data, identifying enrichment from histone modifications.  
- [BCFtools](https://samtools.github.io/bcftools/bcftools.html) processes WGS data into a list of variant calls (VCF). It also processes VCF files in specific ways, such as contrasting two VCF files to list statistically significant differences in genetic variants.

## Analysis Files
Analysis, background and results presented in 2018 are under folder [Analyses](https://github.com/GiulioSP/2018Epigenomics/blob/main/Analyses)
- countAnalysis.plots.xlsx 
- Lab Meeting - Aug 13th.pptx
- Lab Meeting - Aug 22nd.pptx
- Lab Meeting - July 4th.pptx

## Background
Select pages from document "Lab Meeting - Aug 22nd.pptx"
<kbd> ![ppt_what](/new_files/ppt_what.PNG) </kbd>
<kbd> ![ppt_how](/new_files/ppt_how.PNG) </kbd>
<kbd> ![ppt_models](/new_files/ppt_models.PNG) </kbd>
<kbd> ![ppt_data_tables](/new_files/ppt_data_tables.PNG) </kbd>
<kbd> ![ppt_analysis](/new_files/ppt_analysis.PNG) </kbd> 

## Results and Discussion

### Sanity checks and Validation results

Tumors of patients  5 and 6 are considered "hypomutated", as they have significantly less variant calls than would be expected, especially considering the biological baseline (the amount of variant calls that are unique between normal and tumor adjacent samples, which may be real mutations or false positives). 
Hypomutated tumors are worth studying, but their data cannot be compared directly with other samples.
Number of mutations (variant calls) from each sample:
<kbd> ![ppt_hypomutated](/new_files/ppt_hypomutated.PNG) </kbd>

As expected, the amount of unique mutations (not in any list of variants within BCFTools package) is higher for some tumor samples. 
<kbd> ![ppt_tummutation](/new_files/ppt_tummutation.PNG) </kbd>

As expected, most DMRs (differentially methylated regions) are unique to tumors, whereas normal and tumor adjacent samples are most similar. Further, the hypo or hypermethylation profile of DMRs is different for regions considered to be CGIs, as they carry the methylation signal differently. 
Furthermore, tumors 5 and 6 behave differently as they are undermutated. 

Tumor from P6 contains less DMRs than expected, this indicates support to models 1 and 2, where a less mutated sample would have less changes in gene expression and thus less changes in methylation profile (model 1) or less visible genic signals for methylation or demethylation (model 2). 

Tumor from P5 contains less DMRs, but also a unique distribution between normal, tumor adjacent and tumor samples. This indicates support for model 4, where mutations in genes coding for epigenetic machinery cause global dysruption to the methylome.  
|[]() | |
| :----: | :----: | 
|<kbd> ![ppt_tumdmr1](/new_files/ppt_tumdmr1.PNG) </kbd> | <kbd> ![ppt_tumdmr2](/new_files/ppt_tumdmr2.PNG) </kbd> |
|<kbd> ![ppt_tumdmr3](/new_files/ppt_tumdmr3.PNG) </kbd> | <kbd> ![ppt_tumdmr4](/new_files/ppt_tumdmr4.PNG) </kbd> |
|<kbd> ![ppt_tumdmr5](/new_files/ppt_tumdmr5.PNG) </kbd> | <kbd> ![ppt_tumdmr6](/new_files/ppt_tumdmr6.PNG) </kbd> |

### Colocalization of Mutations (VC) and Differential Methylation (DMRs)
These is indeed strong colocalization of epigenetic and genetic variants. These plots consider the DMRs as variations between different samples of the same patient, so they should include few DMRs between normal and tumor adjacent, and consider many regions where the tumor has dysregulated the methylome. This indicates support for model 4, where mutations in genes coding for epigenetic machinery cause global dysruption to the methylome, particularly for causing hypermethylation of the genome.

Note in data from patients P1 and P6 that both normal and tumor adjacent samples happen to have a similar amount of variant calls tagged as "modifier" in all sample types, but the amount of those that colocalize with DMRs is still small, indicating support for models 1 and 2 in healthy cells, and model 4 in cancer cells.  
|[]() | |
| :----: | :----: | 
| <kbd> ![P1.VC+DMR](/new_files/P1.VC+DMR.png) </kbd> | <kbd> ![P2.VC+DMR](/new_files/P2.VC+DMR.png) </kbd> |
| <kbd> ![P3.VC+DMR](/new_files/P3.VC+DMR.png) </kbd> | <kbd> ![P4.VC+DMR](/new_files/P4.VC+DMR.png) </kbd> |
| <kbd> ![P5.VC+DMR](/new_files/P5.VC+DMR.png) </kbd> | <kbd> ![P6.VC+DMR](/new_files/P6.VC+DMR.png) </kbd> |

### Colocalization of Histone Modifications (HM) and Differential Methylation (DMRs) 
The histone modifications are less common, or at least harder to detect with ChIP-Seq compared to methylation signals, but they still show patterns between tumors and other samples. 

Firstly, the absolute amount of histone markers seems to be comparable across sample types for any given patient. 
<kbd> ![ppt_hist](/new_files/ppt_hist.PNG) </kbd>

In particular, the combination of hypermethylation and repressive histone marks is particularly rare in normal samples and very common in tumor samples, perhaps because a hypermethylated region is most commonly a CGI in normal samples that is hypermethylated in the tumor cell thus deactiating the gene. Furthermore, the most common association of hypermethylation with repressive histone marks and hypomethylation with active histone marks is normal, so at least some of the epigenetic machinery seems to be preserved. This suggests that the histone markers may be responding to the methylation profile.    
|[]() | |
| :----: | :----: | 
| <kbd> ![P1.HM+DMR](/new_files/P1.HM+DMR.png) </kbd> | <kbd> ![P2.HM+DMR](/new_files/P2.HM+DMR.png) </kbd> |
| <kbd> ![P3.HM+DMR](/new_files/P3.HM+DMR.png) </kbd> | <kbd> ![P4.HM+DMR](/new_files/P4.HM+DMR.png) </kbd> |
| <kbd> ![P6.HM+DMR](/new_files/P6.HM+DMR.png) </kbd> | <kbd> ![P6.HM+DMR](/new_files/P6.HM+DMR.png) </kbd> |


### The Cancer's Epigenome
The amount of DMRs varies by grouping and patient, but the amount of variant regions between cancer samples is small, indicating support to model 4, where mutations in genes coding for epigenetic machinery cause global dysruption to the methylome.  
<kbd> ![ppt_dmr](/new_files/ppt_dmr.PNG) </kbd>

The profile of dysregulation of the methylome seems to vary between tumors, and may be an indicator of cancer evolution alongside genetic factors. Note how P5 and P6 in particular seem to have fewer hypermethylated regions (expected for hypomutated tumors) but P5 has an incomparable amount of hypomethylated regions. This demonstrates that epimutations may be a missing element to understanding of how hypomutated tumors can survive.
<kbd> ![Tum.DMR](/new_files/Tum.DMR.png) </kbd>

Furthermore, the variation is unique for CGIs. Note how tumor P1 seems to have a unique profile of CGI hypomethylation, whereas P5 and P6 display CGI hypermethylation. Both observations indicate a dysfunctional epigenetic machinery in these tumors, thus lending support to model 4.
<kbd> ![Tum.DMR.CGI](/new_files/Tum.DMR.CGI.png) </kbd>

Consider also that variant calls colocalize with DMRs in comparable amounts across samples, which is impressive when considering the drastic hypomutation of samples P5 and P6. Perhaps the interplay between mutations and epimutations may be a more accurate monitor of cancer evolution than only considering genetic mutations.
<kbd> ![Tum.VC+DMR](/new_files/Tum.VC+DMR.png) </kbd>

Histone markers seem to be comparable across all tumor samples. 
<kbd> ![Tum.HM](/new_files/Tum.HM.png) </kbd>

However, when observed in colocalizatin with DMRs, we can see that they are not equally distributed as active or repressive marks, and may be responding to specific methylation signals, thus being most visible in the hypomethylation in P5 and hypermethylation in P1. 
Still, the epigenetic interplay is between active marks with hypomethylation and repressive marks with hypermethylation, which is still normal behavior. 
<kbd> ![Tum.HM+DMR](/new_files/Tum.HM+DMR.png) </kbd>













