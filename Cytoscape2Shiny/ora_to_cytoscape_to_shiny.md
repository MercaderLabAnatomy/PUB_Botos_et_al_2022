### ORA_2_Cytoscape_2_R_2_Shiny (OCRS)
##### Mini tutorial on how to export a over-representation-analyis result to Cytoscape EnrichmentMap format use AutoAnnotate on it and then export it to RStudio again and create a network in Shiny displaying almost all the information.


Starting from a gene list which have been already filtered for the specified tresholds and converted to Entrez IDs, it is possible to perform and enrichment of those selected genes using multiples databases such as GO,KEGG,GSEA,etc.

* Perform enrichment of the genes selected:
	+ `enrichGO(d1,  				#gene dataset 
		 OrgDb = orgdb,                         #organism used such as MusMusculus or Danio rerio  
                 ont = "BP",			        #here we are using GO:BP
                 pAdjustMethod = "BH",                  #correction method Benajamini-Hochberg
                 pvalueCutoff  = pvalCutOff,		#define a pvalue cutt off (0.05)
                 readable = T)				#readable so translated entrez ids back to gene symbols`

* With the result create another dataframe which contains only the statistical significant if interested following a the next structure:
	+ Symbol\tDescription\tSample1\tSample2\tSample3\tSamplen(counts values for the samples)
	+ GO:ID\tDescription\tp.Val\tFDR\tPhenotype(optional)\tGenes(involved in the GO:ID)
* With the previous table we can upload the data to Cytoscape after having installed the EnrichmentMap plugin.
* Next step in our situation as we downloaded only the GO:BP annotation we need to retrieve all the GO:IDs and from their respective database, follow baderlab enrichmentmap gene sets respective date GO:BP in our case without infered GO:BP [link](http://download.baderlab.org/EM_Genesets/August_03_2022/Mouse/symbol/).
* With this information we could create our first network in EnrichmentMap.

* Next point would be having installed AutoAnnotate in the Cytoscape apps.
* After installing it,we can select the network and use it with different clustering algorithms to perform a groupping of our go processes and give us a label of the most important words in that group.
	+ For this click in Networks to the left panes and select a network by clicking on it.
	+ Next upp in the main bar click in Apps and select AutoAnnotate and select New Annotation Set
	+ Next select the cluster algorithm (MCL)
	+ Select edge weight using similarity coefficient
	+ Click in the little square for Layout network to prevent cluster overlap
	+ Then in the AutoAnnotate pane go the top where + refresh and triple bar buttons shall be and click on the tripple bar
	+ Select Layout Clusters and select again CoSE Layout to open the clusters and make it more readable.
	+ Furthr down to the right corner HeatMap can be selected and tunned up to show the best HeatMap possible based on the selected cluster or genes in the network.


To recrete this in to Shiny the network shall be downloaded, at the moment of writing this September 2022 the best format for me was using Export file as GraphML

Start with the wrangling using R,Cytoscape and optional create a table in Excel and associate it to the graph which we will load into R.

* From the Cytoscape exported GraphML, with RStudio we can load it using igraph.
	+ Start with `library("igraph")`
	+ Read the graph in R: `graph_object <- read_graph(file = "./data/yourGraphPath.graphml",format = "graphml")`

* Next we need to associate each AutoAnnotate generated cluster with the GO:IDs inside so we are able to graph and select them later.
* For this a column with the AutoAnnotate Cluster name shall be created\t and a column with the respective GO:IDs involved
* Add these dataframe to the graph, code is available in shiny_app folder.
* Examples from this steps:
	+ `cmc <- c("GO:0000075", "GO:0007088", "GO:0007091", "GO:0007093", "GO:0007094", "GO:0007096", "GO:0010639", "GO:0010948", "GO:0010965", "GO:0030071", "GO:0031577", "GO:0033044", "GO:0033045", "GO:0033046", "GO:0033047", "GO:0033048", "GO:0044770", "GO:0044772", "GO:0044784", "GO:0045786", "GO:0045839", "GO:0045841", "GO:0051304", "GO:0051306", "GO:0051783", "GO:0051784", "GO:0051983", "GO:0051985", "GO:0071173", "GO:0071174", "GO:0090231", "GO:0090266", "GO:1901976", "GO:1901987", "GO:1901988", "GO:1901990", "GO:1901991", "GO:1902099", "GO:1902100", "GO:1903504", "GO:1905818", "GO:1905819", "GO:2000816", "GO:2001251")`
	+ `vertex_attr(core)$AA <- case_when(vertex_attr(core)$name %in% cmc ~ "checkpoint mitotic cycle"...`
	+ `long_list_to_fill_core <- lapply(seq_along(unique(vertex_attr(core)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(core)$AA),function(n) {
  if_else(is.na(vertex_attr(core)$AA[[n]]) == TRUE,
          long_list_to_fill_core[n] <- vertex_attr(core)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_core[n] <- vertex_attr(core)$AA[[n]])
})
vertex_attr(core)$Clusters <- unlist(AA_Clusters)`


* And with this shall be ready to plot using igraph and visnetwork.
	+ `library("visNetwork")`
	+ `output$amputated_net = renderVisNetwork({
    visNetwork::visIgraph(amp) |>
      visNodes(color = "#CC6699") |>
      visOptions(selectedBy = "Clusters") |>
      visInteraction(navigationButtons = TRUE) |>
      visNetwork::visLayout(randomSeed = 123456) |>
      visEvents(select = "function(nodes) {Shiny.onInputChange('current_node_id', nodes.nodes);;}")
  })`
* This last code chunck contains the last line which is and observat event in shiny waiting and collecting the information of which node was selected to select it later as shown in the shiny_app code to select the datatable and the heatmap interactively with the clicks of the user.
* Further future a function to integrate all of this in R and display in Shiny might come `(TODO)`

*DONE*

