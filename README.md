# nematode-hts-toxicology
This repository holds the data and code to reproduce the analyses from the the manuscript:

> Crombie TA, *et al.* (2025). **High-Throughput Toxicity Screening with C. elegans: Current Platforms, Key Advantages, and Future Directions**.

## Usage
1. Clone the repo with `git clone https://github.com/Crombie-Lab/nematode-hts-toxicology.git`
2. Download the 5.8 Gb `INVITRODB_V4_1_SUMMARY.zip` file and move the required sourcefiles into the `/data/raw/toxcast_data/` directory (see [TOXCAST](###9.-TOXCAST)). 
2. Open the `01_data_processing.R` script in Rstudio.
    * Run the code to see how the raw data are cleaned, formatted, and joined.
3. Open and run the remaining `.R` files sequentially to reproduce the analyses in the paper.
4.  The `functions.R` code contains all the orthogonal regression functions to compare across species and platforms.

## Data Sources
### 1. Widmayer *et al.* 2022
    * Raw file path: `data/raw/Data_Andersen_All.xlsx`
    * Link: [View Internal data from Andersen Lab](https://github.com/Crombie-Lab/nematode-hts-toxicology/blob/main/data/raw/Data_Andersen_All.xlsx)
    * Paper DOI: 
### 2.  Boyd *et al.* 2016
    * Raw file path: `/data/raw/Boyd/ehp.1409645.s002.xlsx`
    * Link: [Download Boyd *et al.* 2016 supplemental data](https://ehp.niehs.nih.gov/doi/suppl/10.1289/ehp.1409645/suppl_file/ehp.1409645.s002.acco.zip)
    * Paper DOI:
### 3. EPA CompTox
    * Raw file path: `/data/raw/Comptox/20240717_comptox_download.xlsx`
    * Description: CompTox chemical dashboard - Batch Search
    * Link:  [View CompTox dashboard](https://comptox.epa.gov/dashboard/batch-search)
### 4. EnviroTox DB
    * Raw file path: `/data/raw/envirotox_20240729124104.xlsx`
    * Description: EnviroTox Database Download
    * Instructions: Click link below, then click on "Advanced" tab, then select "TEST" from “Select a Field” dropdown. Click “Search”, then click “Download as Excel File”.
    * Link: [View EnviroTox Database](https://envirotoxdatabase.org/)
### 5. NIEHS_ICE
    * Raw file path: `/data/raw/NIEHS_ICE/Acute_Oral_Toxicity.xlsx` 
    * Link: [Download file from ICE](https://ice.ntp.niehs.nih.gov/downloads/DataonICE/acute_oral.xlsx)
### 6. Karmaus *et al.* 2022
    * Raw file path: `/data/raw/Karmaus/toxsci-21-0357-File010.xlsx`
    * Link: [Download Karmaus *et al.* 2022 supplemental data](https://pmc.ncbi.nlm.nih.gov/articles/instance/9237992/bin/kfac042_supplementary_data.zip)
    * Paper DOI:
### 7. Scholz *et al.* 2016
    * Raw file path: `/data/raw/Zebrafish/Scholz et al 2016/annex2_fet_en.xlsx` 
    * Link: [Download Scholz *et al.* 2016 supplemental data](https://www.sciencedirect.com/science/article/pii/S0045653516311055?via%3Dihub#appsec1)
    * Paper DOI:
### 8. Su *et al.* 2021
    * Raw file path: `/data/raw/Zebrafish/Su et al 2021/1-s2.0-S0048969721027765-mmc2.xls` 
    * Link: [Download Su *et al.* 2021 supplemental data](https://ars.els-cdn.com/content/image/1-s2.0-S0048969721027765-mmc2.xls)
    * Paper DOI:
### 9. TOXCAST
    * Raw file paths: 
        - `/data/raw/toxcast_data/mc5-6_winning_model_fits-flags_invitrodb_v4_1_SEPT2023.Rdata`
        - `/data/raw/toxcast_data/assay_annotations_invitrodb_v4_1_SEPT2023.xlsx`
    * Description: Screening data from the US EPA's ToxCast program. These files were downloaded within the `INVITRODB_V4_1_SUMMARY.zip` file. This zip file is 5.8 Gb, therefore we do not host it locally in this repository. Please download the file at the link below to recreate the full analysis from RAW data.
    * Link: [View INVITRODB_V4_1_SUMMARY.zip download site](https://clowder.edap-cluster.com/files/64bfdb62e4b08a6b5a434d48)
