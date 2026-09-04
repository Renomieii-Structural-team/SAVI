# ProtSAVI

<p align="center">
  <img src="ProtSAVI.png" alt="SAVI logo" width="35%">
</p>

**ProtSAVI** (Protein Structure-Assisted Variant Interpretation) is an R-based dashboard that offers a transparent, user-friendly framework for variant interpretation and classification. The tool integrates ACMG/AMP criteria, literature data, functional annotations, in silico predictions, and structural analyses to interpret and reclassify genetic variants.

The increasing number of Variants of Uncertain Significance (VUS) poses a major challenge for clinicians, researchers, and patients. Re-interpretation of these variants can substantially influence diagnostic outcomes for individuals with inconclusive results. By integrating structural biophysical analyses, AI-driven predictions, and functional data on target gene-encoded proteins, we seek to provide actionable insights for variant interpretation.

🔗 **Access ProtSAVI:** https://paulafranklin.shinyapps.io/ProtSAVI/

To explore the platform, please use the input examples provided in this repository. For each gene available in the dashboard, we include variants classified as VUS according to ClinVar.

## Running ProtSAVI locally

To run **ProtSAVI** locally, you need to have **R** and **RStudio** installed, along with the package versions listed below.

1. Clone this repository or download it as a ZIP file.
2. Open the project in **RStudio**.
3. Install the required packages (see table below).
4. Run the `ProtSAVI.R` script to launch the application.

### Required package versions

ProtSAVI was developed and tested using the following package versions:

| Package | Version |
|----------|---------|
| here | 1.0.1 |
| rio | 1.2.3 |
| tidyverse | 2.0.0 |
| readr | 2.1.4 |
| stringr | 1.5.1 |
| colorspace | 2.1.0 |
| shiny | 1.8.1.1 |
| shinyWidgets | 0.8.6 |
| DT | 0.28 |
| NGLVieweR | 1.3.4 |
| plotly | 4.12.0 |
| shinyjs | 2.1.0 |
| bslib | 0.7.0 |
| ggplot2 | 3.5.1 |

