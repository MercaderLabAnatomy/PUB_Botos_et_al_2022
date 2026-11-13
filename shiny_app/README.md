# **Transcriptomic data meta-analysis reveals common and injury model specific gene expression changes in the regenerating zebrafish heart**

**Marius Botos<sup>1,2</sup>, Prateek Arora<sup>1,2</sup>, Panagiotis Chouvadras<sup>2,3</sup>, Nadia Mercader<sup>1,2,4</sup>**

*1 Institute of Anatomy, University of Bern, 3012 Bern, Switzerland*

*2 Department for Biomedical Research, University of Bern, 3012, Bern, Switzerland*

*3 Department of Urology, Inselspital, Bern University Hospital, Bern, Switzerland*

*4 Centro Nacional de Investigaciones Cardiovasculares CNIC, 28029 Madrid, Spain*

### **Abstract**

##### Zebrafish have the capacity to fully regenerate the heart after an injury, which lies in sharp contrast to the irreversible loss of cardiomyocytes after a myocardial infarction in humans. Transcriptomics analysis has contributed to dissect underlying signaling pathways and gene regulatory networks. Heart regeneration can be studied in response to different types of injuries. The three more broadly used techniques are ventricular resection, ventricular cryoinjury, and genetic ablation of cardiomyocytes. To which extent these different types of injuries elicit an equivalent regenerative response is still unclear. Here, we present a meta-analysis of transcriptomic data of regenerating zebrafish hearts in response to these three injury models. We reanalyzed and batch corrected 36 samples distributed in seven datasets and analyzed the differentially expressed genes (DEG) followed by downstream analysis of Gene Ontology Biological Processes (GO\:BP). We found that the three injury models share a common core of DEG encompassing genes involved in cell proliferation and the Wnt signaling pathway. We also found injury-specific gene signatures for resection and genetic ablation, and to a lower extent the cryoinjury model. Interestingly, the DEG analysis also led to specific GO\:BP. Our results identify a gene expression signature common to all analyzed models but also highlight the importance to consider injury-specific gene regulatory networks when interpreting results on cardiac regeneration in the zebrafish.

## Workflow implemented in the meta-analysis:


<img src="./www/Figure1_corrected_Shape_and_Index_v4_calibri-2.png" alt="Workflow implemented in the meta-analysis" style="width:720px;"/>

*PS:* **the genes found differentially expressed in *Danio rerio* were converted to the correspondent orthologues in *Mus musculus* and used for the whole analysis in the Over-representation analysis.**

## How to use:

-   **1:** In the left side of the browser window we can find the title and three horizontal lines(menu) when **clicking** it will show or hide the navigation bar.
-   **2:** In the Navigation Bar we can find multiple tabs (Home which we are here now)
-   **3:** Followed by the comparisons analyzed with the clustered Gene Ontology Biological Processes.
-   **4:** Followed by the Biological Processes involved in the three injury models named Heart Core Regeneration.
    -   **4.1:** In the left corner a selection of clusters that we obtained from the network display after using Cytoscape and the tools available for clustering the GO\:BP is available.
    -   **4.2:** Select a cluster to be specifically shown in the Network representation of the clustered GO\:BP.
    -   **4.3:** By clicking in a specific GO\:BP of the cluster, the information within this GO\:BP will be displayed and a Heatmap of the genes shown for this pairwise comparison.
-   **5:** Followed by the PubMed query realized using the genes found in the Heart Core Regeneration.
-   **6:** Followed by the Heatmap of the Core Regeneration genes in the three injury models and their respective up or down regulation when compared to sham.
-   **7:** Followed by Gene Seeker, where minimum 2 genes shall be selected in order to see a Heatmap and a data table of expression among all samples, and the respective description of each gene.

### 
