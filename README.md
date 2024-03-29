Shiny app: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MercaderLabAnatomy/PUB_Botos_et_al_2022_shinyapp_binder/HEAD?urlpath=shiny/bus-dashboard/)

### PUB_Botos_et_al_2022
#### Transcriptomic data meta-analysis reveals common and injury model specific gene expression changes in the regenerating zebrafish heart.

### **Meta-analysis code structured in different analysis steps.**

* **Download data with the `raw` folder**
* **Pre-process the data following the sequential steps:`multiqc` to `fastp` to optional rechech with `multiqc` to indexing the reference genome and prepare mapping with `STAR` within `indexing` to `mapping`**
* **Processing of the data using `R` code can be found in folder `r_processing`**
* **Pub Med query the genes of interest with `pubmed_query_processing`**
* **Web Application follows and requires analysis with Cytoscape and further R integration**
* **Shiny app code is in the folder `shiny_app`**

P.S - The app usually takes 2-3 minutes to load, if it takes longer please reload the link

