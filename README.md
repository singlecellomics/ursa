![Ursa logo](https://user-images.githubusercontent.com/5945741/165857896-912bfe07-f290-483c-bb96-d5ff21db1ab6.png)

# Ursa: an automated multi-omics package for single-cell analysis

__Ursa__ is an R package consisting of seven single-cell omics automated analysis workflows. One-liner command for each omics to run a full post-quantification analysis for the omics.

Six single-cell (sc) omics and one bulk omics include:

1. scRNA-sequencing (sc)
2. scATAC-sequencing (sc)
3. scImmune profiling (sc)
4. scCNV (sc)
5. CyTOF (sc)
6. Flow cytometry (sc)
7. Spatial transcriptomics (bulk)

## Installation

Ursa can be installed in R via the command:
```sh
install.packages("devtools")
devtools::install_github("eudoraleer/Ursa")
```
Please download the example sample files and their meta files from the following [__link__](https://www.dropbox.com/sh/zdi0554bf07spoo/AAAZNk_jsrFa53tg4CsGfU2ua?dl=0) with respect to the omics you will be running. Original file sources can be found below. Multiple samples are supported if information in the meta data is corrected provided.

## Running single-cell analysis with Ursa
### 1. scRNA-sequencing*
#### (1) Download example dataset from original source [__10X__](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high) or from the following [__link__](https://www.dropbox.com/sh/6q75ik2egtfai7q/AABkXelU7Iyz_cWmbdtSlpUMa?dl=0)
The following input file(s) from the example data are needed in the input directory before running the analysis:
- filtered gene matrix .h5 file: Feature / cell matrix HDF5 (filtered)
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195846978-3091c9a7-c5c6-4217-a39f-1450c1c3a55e.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195845913-84d8b84f-49fd-4b50-9fd6-03622eb49958.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
scRNASEQPip(project_name = 'My_scRNASeq', pheno_file = 'Ursa_scRNA_Seq_Metadata_Example.csv')
```
#### (4) Example output files for project My_scRNASeq: [__link__](https://www.dropbox.com/sh/triv03adukw2pp3/AAAYLKlcfy2zuhHSezYJ_Voca?dl=0)

### 2. scATAC-sequencing*
#### (1) Download example dataset from original source [__10X__](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard) or from the following [__link__](https://www.dropbox.com/sh/o5qx7coly4mp7l2/AABMSlfK2I6sIsdtkqM6Vkvja?dl=0)
For this omics, running this workflow on a computer with memory >=16GB is recommended due to large input file size
The following input file(s) from the example data are needed in the input directory before running the analysis:
- filtered peak matrix .h5 file: Peak by cell matrix HDF5 (filtered)
- fragment file and its index file: Fragments (TSV), Fragments index (TBI)
- single cell file: Per Barcode metrics (CSV, optional)
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195842755-a8512786-e757-45de-8a16-f439bbdfd232.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195843616-03e607ec-4979-4f7a-a168-fc5341ad7576.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
scATACPip(project_name = 'My_scATAC', pheno_file = 'Ursa_scATAC_Seq_Metadata_Example.csv')
```
#### (4) Example output files for project My_scATAC: [__link__](https://www.dropbox.com/sh/uwtb2gmw1vob94b/AAC4wDoYMqboF6z78roqvAr7a?dl=0)

### 3. scImmune profiling*
#### Download example dataset from original source [__10X__](https://www.10xgenomics.com/resources/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0) or from the following [__link__](https://www.dropbox.com/sh/03q8kpp5fmzcqf5/AAAGoGxEX9Ma4EGUs762i7B6a?dl=0)
The following input file(s) from the example data are needed in the input directory before running the analysis:
- BCR or/and TCR contig CSV file: VDJ Ig - All contig annotations (CSV)
- filtered gene matrix .h5 file (optional, only for multi-modal analysis): Gene Expression - Feature / cell matrix .h5 file (filtered)
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195844324-4956e9db-4d93-4c4e-be2c-667ab2b57309.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195845640-0a013558-6b42-4e5c-8e0e-58d7ef6198a4.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
scImmunePip(project_name = 'My_scImmune', pheno_file = 'Ursa_scImmune_Profiling_Metadata_Example.csv')
```
#### (4) Example output files for project My_scImmune: [__link__](https://www.dropbox.com/sh/u2cg56duniwr890/AADNnSK4rvbdgRm4f3IUU1FYa?dl=0)

### 4. scCNV*
#### Download example dataset from original source [__10X__](https://www.10xgenomics.com/resources/datasets/breast-tissue-nuclei-section-a-2000-cells-1-standard-1-1-0) or from the following [__link__](https://www.dropbox.com/sh/jp3gc0sigvt849g/AABQnEmxfdxJidwWdCxf7pz3a?dl=0)
The following input file(s) from the example data are needed in the input directory before running the analysis:
- mappable regions BED file: Mappable regions (BED)
- CNV calls: CNV calls (BED)
- per cell summary metrics: Per-cell summary metrics (CSV)
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195843861-b8672fc2-3b95-467b-b06e-b998dee084b9.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195844194-52d05ef9-daef-4641-89a5-fe3e6b4a1521.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
scCNVPip(project_name = 'My_scCNV', pheno_file = 'Ursa_scCNV_Metadata_Example.csv')
```
#### (4) Example output files for project My_scCNV: [__link__](https://www.dropbox.com/sh/aqlc10ami53fn85/AAAnWUx0Ic4uXOx46v5-EFRga?dl=0)

### 5. CyTOF
#### Download example dataset from original source [__Nowicka, M., et al. (2017)__](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_fcs_files.zip) or from the following [__link__](https://www.dropbox.com/sh/wfn4vhauj8s8zm5/AADlEbxJ_quTyQd10cLadqQBa?dl=0)
The following input file(s) from the example data are needed in the input directory before running the analysis:
- .fcs input files
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195842654-eaa061b0-adde-47ea-b5e1-28092ed10adc.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195840736-ee101304-4803-49e6-97e3-42cd3e78ebb1.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
CyTOFPip(project_name = 'My_CyTOF', pheno_file = 'Ursa_CyTOF_Metadata_Example.csv')
```
#### (4) Example output files for project My_CyTOF: [__link__](https://www.dropbox.com/sh/f3ip2znr9enmloa/AACw4GROCndSQwuxCpnNjaTUa?dl=0)

### 6. Flow Cytometry
#### Download example dataset from original source [__Dillon Hammill,2021__](https://github.com/DillonHammill/CytoExploreRData/tree/master/inst/extdata/Activation) or from the following [__link__](https://www.dropbox.com/sh/wlypurz70knlb32/AACK-s8SjwBBispS5Y0Ylopta?dl=0)
The following input file(s) from the example data are needed in the input directory before running the analysis:
- .fcs input files
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195842509-1229430f-9acd-4a11-b8dd-0e1983b85848.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195842219-d09218b5-c7d8-4709-a7ce-7fb8f8de0eec.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
FlowPip(project_name = 'My_Flow', pheno_file = 'Ursa_Flow_Cytometry_Metadata_Example.csv')
```
#### (4) Example output files for project My_Flow: [__link__](https://www.dropbox.com/sh/pwy395cl4f4tncm/AADwMWt0_tVoNbre9Ge0xld7a?dl=0)

### 7. Spatial Transcriptomics
#### Download example dataset from original source [__10X__](https://www.10xgenomics.com/resources/datasets/human-cervical-cancer-1-standard) or from the following [__link__](https://www.dropbox.com/sh/h02jr6l0f2ox9wd/AAAYQZ681WIcI39NKkKt4hbJa?dl=0)
The following input file(s) from the example data are needed in the input directory before running the analysis:
- filtered gene matrix .h5 file: Feature / barcode matrix HDF5 (filtered)
- spatial imaging data: Spatial imaging data (please make sure the imaging data for each sample is placed in their corresponding folder with the .h5 file, with imaging data folder named 'spatial')
- sample meta file (in .csv format) with the following file content:
![image](https://user-images.githubusercontent.com/5945741/195847522-69d5aa07-aeaa-43e7-8317-fe4d83dad42e.png)
#### (2) Set the downloaded file folder as working directory in R/RStudio:
![image](https://user-images.githubusercontent.com/5945741/195847129-63e042e9-9fab-4a47-baf0-2586fe2630d1.png)
#### (3) Run the analysis with the following commands:
```sh
library("Ursa")
SpatialPip(project_name = 'My_Spatial', pheno_file = 'Ursa_Spatial_Metadata_Example.csv')
```
#### (4) Example output files for project My_Spatial: [__link__](https://www.dropbox.com/sh/i6320yizw2uo81c/AACD7zftdCswTkfY_JAON0iVa?dl=0)

###### *Registration is needed for downloading the data for the first time on 10X website. Subsequent download would no longer require any registration.
