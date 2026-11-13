library("rlang")
library("sodium")
library("devtools")
library("shinymanager")
library("DT")
library("dplyr")
library("fontawesome")
library("shinydashboard")
library("shinycssloaders")
library("igraph")
library("rmarkdown")
library("markdown")
library("ggplot2")
library("visNetwork")
library("heatmaply")
library("plotly")
library("brew")
library("RColorBrewer")
library("curl")
library("wesanderson")
library("shinyauthr")

# print("Starting server.R execution...")

rm(list=ls())
gc()
# options(repos = BiocManager::repositories())
#.libPaths("/home/mbotos/R/x86_64-pc-linux-gnu-library/4.2/")


inactivity <- "function idleTimer() {
var t = setTimeout(logout, 120000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 120000);  // time is in milliseconds (1000 is 1 second)
}

idleTimer();"
  # cols_core <- s4c
  # cols_core$ID <- gsub('"', '', s4c[[1]])
  # cols_core <- cols_core |> dplyr::filter(ID %in% colnames(core_mat))
  # cols_core <- s4c
  # cols_core$ID <- gsub('"', '', s4c[[1]])
  # cols_core <- cols_core |> dplyr::filter(ID %in% colnames(core_mat)) 
# data.frame with credentials info
credentials <- data.frame(
  user = c("marius","nadia","panos","prateek", ""),
  password = c("marius","nadia","panos","prateek", ""),
  start = c("2022-04-01","2022-04-01","2022-04-01","2022-04-01", "2022-04-01"),
  expire = c(NA,NA,NA,NA, NA),
  admin = c(TRUE,FALSE,FALSE,TRUE, TRUE),
  comment = "Simple and secure authentification mechanism
  for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)




#################################
#                               #
#   Data Load                   #
#                               #     
#################################
#     Ablated                   #
#################################
##:%s/^\(.*\)$/"\1",/
#cat NN | tr "\n" " "
#excel to paste and make the final command
#setwd("/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/ResultsCytoscape/PossibleShiny/RShinyCode/ShinnyNetworksGOBP")
# print("Loading Ablated graph...")
abl <- read_graph(file = "./data/Ablated.graphml",format = "graphml")

# Function to get cluster members from a seed
get_cluster_members <- function(graph, seed) {
  cluster_attr <- vertex_attr(graph)$`__mclCluster`
  idx <- match(seed, vertex_attr(graph)$name)
  if (is.na(idx)) return(character(0))
  seed_cluster <- cluster_attr[idx]
  vertex_attr(graph)$name[which(cluster_attr == seed_cluster)]
}

apply_cluster_mapping <- function(graph, csv_path) {
  # print(paste("Applying cluster mapping from:", csv_path))
  cluster_map <- read.csv(csv_path)
  vertex_attr(graph)$AA <- rep(NA, vcount(graph))
  
  for(i in 1:nrow(cluster_map)) {
    seed <- cluster_map$GO_ID[i]
    c_name <- cluster_map$Cluster_Name[i]
    
    members <- get_cluster_members(graph, seed)
    
    if(length(members) > 0) {
       idx <- match(members, vertex_attr(graph)$name)
       idx <- idx[!is.na(idx)]
       vertex_attr(graph)$AA[idx] <- c_name
    }
  }
  # print("Cluster mapping applied.")
  return(graph)
}

abl <- apply_cluster_mapping(abl, "./data/Ablated_cluster_names_with_GO_IDs.csv")


long_list_to_fill_abl <- lapply(seq_along(unique(vertex_attr(abl)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(abl)$AA),function(n) {
  if_else(is.na(vertex_attr(abl)$AA[[n]]) == TRUE,
          long_list_to_fill_abl[n] <- vertex_attr(abl)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_abl[n] <- vertex_attr(abl)$AA[[n]])
})

vertex_attr(abl)$Clusters <- unlist(AA_Clusters)

#################################
#  Resection                    #
#################################
# print("Loading Resection graph...")
amp <- read_graph(file = "./data/Resection.graphml",format = "graphml")
amp <- apply_cluster_mapping(amp, "./data/resection_cluster_names_with_GO_IDs.csv")

long_list_to_fill_amp <- lapply(seq_along(unique(vertex_attr(amp)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(amp)$AA),function(n) {
  #unj_2 <- list()
  #print(n)
  if_else(is.na(vertex_attr(amp)$AA[[n]]) == TRUE,
          long_list_to_fill_amp[n] <- vertex_attr(amp)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_amp[n] <- vertex_attr(amp)$AA[[n]])
})


vertex_attr(amp)$Clusters <- unlist(AA_Clusters)


#################################
#   Uninjured                   #
#################################
# print("Loading Uninjured graph...")
unj <- read_graph(file = "./data/Uninjured.graphml",format = "graphml")
unj <- apply_cluster_mapping(unj, "./data/uninjured_cluster_names_with_GO_IDs.csv")

long_list_to_fill_unj <- lapply(seq_along(unique(vertex_attr(unj)$AA)),function(i) i)
AA_Clusters <- lapply(seq_along(vertex_attr(unj)$AA),function(n) {
  if_else(is.na(vertex_attr(unj)$AA[[n]]) == TRUE,
          long_list_to_fill_unj[n] <- vertex_attr(unj)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_unj[n] <- vertex_attr(unj)$AA[[n]])
})

vertex_attr(unj)$Clusters <- unlist(AA_Clusters)



#################################
#     Core regeneration         #
#################################
# print("Loading Core graph...")
core <- read_graph(file = "./data/core.graphml",format = "graphml")
core <- apply_cluster_mapping(core, "./data/core_cluster_names_with_GO_IDs.csv")





long_list_to_fill_core <- lapply(seq_along(unique(vertex_attr(core)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(core)$AA),function(n) {
  if_else(is.na(vertex_attr(core)$AA[[n]]) == TRUE,
          long_list_to_fill_core[n] <- vertex_attr(core)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_core[n] <- vertex_attr(core)$AA[[n]])
})
vertex_attr(core)$Clusters <- unlist(AA_Clusters)


###############################
#       Load counts           #
###############################

# print("Loading counts data...")
unj_counts <- read.delim(file = "./data/Platform_DEGs_only_Uninjured_vs_Sham_expression.txt",header = TRUE,sep = "\t")
unj_counts[nrow(unj_counts)+1,] <- 0
unj_counts[nrow(unj_counts)+1,] <- 0
unj_counts[nrow(unj_counts),"MGI_Symbol"] <- "Select_a_node"
unj_counts[nrow(unj_counts)-1,"MGI_Symbol"] <- "Select_a_node"


abl_counts <- read.delim(file = "./data/Platform_DEGs_only_Ablated_vs_Sham_expression.txt",header = TRUE,sep = "\t")
abl_counts[nrow(abl_counts)+1,] <- 0
abl_counts[nrow(abl_counts)+1,] <- 0
abl_counts[nrow(abl_counts),"MGI_Symbol"] <- "Select_a_node"
abl_counts[nrow(abl_counts)-1,"MGI_Symbol"] <- "Select_a_node"

amp_counts <- read.delim(file = "./data/Platform_DEGs_only_Resection_vs_Sham_expression.txt",header = TRUE,sep = "\t")
amp_counts[nrow(amp_counts)+1,] <- 0
amp_counts[nrow(amp_counts)+1,] <- 0
amp_counts[nrow(amp_counts),"MGI_Symbol"] <- "Select_a_node"
amp_counts[nrow(amp_counts)-1,"MGI_Symbol"] <- "Select_a_node"

cry_counts <- read.delim(file = "./data/Platform_DEGs_only_Cryoinjury_vs_Sham_expression.txt",header = TRUE,sep = "\t")
cry_counts[nrow(cry_counts)+1,] <- 0
cry_counts[nrow(cry_counts)+1,] <- 0
cry_counts[nrow(cry_counts),"Symbols"] <- "Select_a_node"
cry_counts[nrow(cry_counts)-1,"Symbols"] <- "Select_a_node"


s4c <- read.delim(file = "./data/colData_s4c.txt",header = TRUE,sep = "\t")
s4c$ID <- rownames(s4c)

#sham_counts <- cry_counts[,c(rownames(s4c %>% filter(Condition=="Sham")),"MGI_Symbol")]



core_counts <- read.delim("./data/core_regeneration_genes.txt",sep = "\t",header = TRUE)
# print(paste("core_regeneration_genes rows/cols:", nrow(core_counts), ncol(core_counts)))
core_counts[nrow(core_counts)+1,] <- 0
core_counts[nrow(core_counts)+1,] <- 0
core_counts[nrow(core_counts),"MGI_Symbol"] <- "Select_a_node"
core_counts[nrow(core_counts)-1,"MGI_Symbol"] <- "Select_a_node"
# print("core_regeneration_genes added Select_a_node rows")


go_abl <- read.delim("./data/Platform_unique_Ablated_vs_Sham_GO.txt", header=TRUE, sep="\t")
go_abl$Group <- "Ablated"
go_abl$geneID <- go_abl$Genes

go_amp <- read.delim("./data/Platform_unique_Resection_vs_Sham_GO.txt", header=TRUE, sep="\t")
go_amp$Group <- "Resection"
go_amp$geneID <- go_amp$Genes

go_unj <- read.delim("./data/Platform_unique_Uninjured_vs_Sham_GO.txt", header=TRUE, sep="\t")
go_unj$Group <- "Uninjured"
go_unj$geneID <- go_unj$Genes

enrich_go_bo_df <- rbind(go_abl, go_amp, go_unj)


core_reg <- read.delim(file = "./data/core_regeneration_GOBP.txt",sep = "\t",header = TRUE)
# print(paste("core_regeneration_GOBP rows/cols:", nrow(core_reg), ncol(core_reg)))


#gene seeker
abl_fil <- read.delim(file = "./data/Platform_DEGs_only_Ablated_vs_Sham_expression.txt",header = TRUE,sep = "\t")
amp_fil <- read.delim(file = "./data/Platform_DEGs_only_Resection_vs_Sham_expression.txt",header = TRUE,sep = "\t")
cryo_fil <- read.delim(file = "./data/Platform_DEGs_only_Cryoinjury_vs_Sham_expression.txt",header = TRUE,sep = "\t")
unj_fil <- read.delim(file = "./data/Platform_DEGs_only_Uninjured_vs_Sham_expression.txt",header = TRUE,sep = "\t")

unj_fil_f <- unj_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)
amp_fil_f <- amp_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)
abl_fil_f <- abl_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)
cryo_fil_f <- cryo_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)



gene_seeker_table <- dplyr::bind_rows(abl_fil_f, amp_fil_f, unj_fil_f, cryo_fil_f)



gene_seeker_table <- gene_seeker_table |> dplyr::group_by(MGI_Symbol) |> distinct(MGI_Symbol,.keep_all = TRUE)




#log2fc heatmap of the core regeneration genes
files <- list.files(all.files = TRUE,
                    path = "./data",
                    pattern = "Platform_DEGs_only.*_expression\\.txt$",
                    recursive = TRUE, 
                    full.names = TRUE)

tables_counts_stats <- lapply(files,function(n) read.delim(n,header = TRUE,sep = "\t"))

name_files <- gsub(x = files,
                   pattern = c("./data/Platform_DEGs_only_"),
                   replacement = "")
name_files <- gsub(x = name_files,
                   pattern = c("_expression.txt"),
                   replacement = "")

names(tables_counts_stats) <- name_files
# print("tables_counts_stats names:")
# print(names(tables_counts_stats))

#create df
cryo_df <- tables_counts_stats$Cryoinjury_vs_Sham %>% dplyr::filter(MGI_Symbol %in% core_counts$MGI_Symbol)
# print(paste("cryo_df dim:", paste(dim(cryo_df), collapse = "x")))
cryo_df <- cryo_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))
# print(paste("cryo_df after select dim:", paste(dim(cryo_df), collapse = "x")))
#dim(cryo_df)


amp_df <- tables_counts_stats$Resection_vs_Sham %>% filter(MGI_Symbol %in% core_counts$MGI_Symbol)
# print(paste("amp_df dim:", paste(dim(amp_df), collapse = "x")))
amp_df <- amp_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))
# print(paste("amp_df after select dim:", paste(dim(amp_df), collapse = "x")))
# print(dim(amp_df))

abl_df <- tables_counts_stats$Ablated_vs_Sham %>% filter(MGI_Symbol %in% core_counts$MGI_Symbol)
# print(paste("abl_df dim:", paste(dim(abl_df), collapse = "x")))
abl_df <- abl_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))
# print(paste("abl_df after select dim:", paste(dim(abl_df), collapse = "x")))
#dim(abl_df)


#abl_df %>% group_by(MGI_Symbol) %>% count(MGI_Symbol)  %>% filter(n > 1)
#abl_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol) %>% count(MGI_Symbol) %>% filter(n > 1)

abl_df_2 <- abl_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)
amp_df_2 <- amp_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)
cryo_df_2 <- cryo_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)


#rm(counts_148_selected_without_doubles)
all.equal(cryo_df_2$MGI_Symbol,abl_df_2$MGI_Symbol)
counts_148_selected_without_doubles <- as.data.frame(amp_df_2$log2FoldChange)
counts_148_selected_without_doubles <- cbind.data.frame(counts_148_selected_without_doubles,abl_df_2$log2FoldChange)
counts_148_selected_without_doubles <- cbind.data.frame(counts_148_selected_without_doubles,cryo_df_2$log2FoldChange)
# print(paste("counts_148_selected_without_doubles dim:", paste(dim(counts_148_selected_without_doubles), collapse = "x")))

rownames(counts_148_selected_without_doubles) <- amp_df_2$MGI_Symbol
colnames(counts_148_selected_without_doubles) <- c("Resection","Ablation","Cryoinjury")




server <- function(input, output, session) {
  
  result_auth <- secure_server(check_credentials = check_credentials(credentials))
  
  output$res_auth <- renderPrint({
    reactiveValuesToList(result_auth)
    
  })
  
  output$resected_net = renderVisNetwork({
    visNetwork::visIgraph(amp) |>
      visNodes(color = "#CC6699") |>
      visOptions(selectedBy = "Clusters") |>
      visInteraction(navigationButtons = TRUE) |>
      visNetwork::visLayout(randomSeed = 123456) |>
      visEvents(select = "function(nodes) {Shiny.onInputChange('current_node_id', nodes.nodes);;}")
  })
  
  output$goTable_amp <- renderDataTable({
    
    enrich_go_bo_df %>%
      group_by(Group) %>%
      filter(Group == "Resection") %>%
      filter(ID %in% vertex_attr(amp)$name[vertex_attr(amp)$name == input$current_node_id])
    
    
    
    
  })
  
  observeEvent(input$current_node_id, {
    visNetworkProxy("resected_net") %>%
      visGetNodes()
  })
  
  output$verb_amp <- renderText({ paste0("Selected GO:BP is : ",
                                         vertex_attr(amp)$`EnrichmentMap::GS_DESCR`[vertex_attr(amp)$name == input$current_node_id])
  })
  
  output$genes_go_amp <- renderValueBox(
    valueBox(
      value = tags$p(enrich_go_bo_df %>% group_by(Group) %>% filter(Group == "Resection") %>% filter(ID %in% vertex_attr(amp)$name[vertex_attr(amp)$name == input$current_node_id]) %>% pull(geneID) %>% gsub(pattern = "/",replacement=","),style = "font-size: 35%;", style ="color : black"),
      color="fuchsia",
      subtitle = paste0("Genes Involved in the ",vertex_attr(amp)$`EnrichmentMap::GS_DESCR`[vertex_attr(amp)$name == input$current_node_id]))
  )
  
  
  genes_amp <-  "Select_a_node"
  # print(genes_amp)
  
  # Select columns based on names in amp_counts
  target_cols_amp <- grep("^SRR", colnames(amp_counts), value = TRUE)
  amp_counts2 <- amp_counts[, c("MGI_Symbol", target_cols_amp), drop = FALSE]
  
  amp_counts_genes <- amp_counts2 %>% filter(MGI_Symbol %in% genes_amp) %>% distinct(MGI_Symbol,.keep_all = TRUE)
  
  amp_mat <- as.data.frame(amp_counts_genes[, -1])
  if(nrow(amp_mat) > 0) {
    rownames(amp_mat) <- amp_counts_genes$MGI_Symbol
  }
  
  # Create annotation dataframe
  cols_amp <- data.frame(ID = colnames(amp_mat))
  cols_amp$Condition <- ifelse(grepl("Sham", cols_amp$ID), "Sham", "Resection")
  cols_amp <- cols_amp[order(cols_amp$Condition, decreasing = TRUE),]
  # print(unique(cols_amp$Condition))
  if(nrow(amp_mat)<=2){cluster=TRUE}else {cluster=FALSE}
  my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
  
  output$amphm <- renderPlotly({
    heatmaply::heatmaply(x = amp_mat |> dplyr::select(cols_amp$ID[order(cols_amp$Condition,decreasing = TRUE)]),
                         
                         Rowv = FALSE,
                         Colv = FALSE,
                         
                         trace="none") 
  }) 
  
  
  observeEvent(input$current_node_id, {
    
    genes_amp <- enrich_go_bo_df %>% group_by(Group) %>% filter(Group == "Resection") %>%
      filter(ID %in% vertex_attr(amp)$name[vertex_attr(amp)$name == input$current_node_id]) %>%
      pull(geneID) %>% gsub(pattern = "/",replacement=",")  %>% strsplit(split = ",") %>% unlist()
    # print(paste(length(genes_amp), "something selected"))
    target_cols_amp <- grep("^SRR", colnames(amp_counts), value = TRUE)
    amp_counts2 <- amp_counts[, c("MGI_Symbol", target_cols_amp), drop = FALSE]
    
    amp_counts_genes <- amp_counts2 %>% filter(MGI_Symbol %in% genes_amp) %>% distinct(MGI_Symbol,.keep_all = TRUE)
    
    amp_mat <- as.data.frame(amp_counts_genes[, -1])
    if(nrow(amp_mat) > 0) {
      rownames(amp_mat) <- amp_counts_genes$MGI_Symbol
    }
    
    cols_amp <- data.frame(ID = colnames(amp_mat))
    cols_amp$Condition <- ifelse(grepl("Sham", cols_amp$ID), "Sham", "Resection")
    cols_amp <- cols_amp[order(cols_amp$Condition, decreasing = TRUE),]
    # print(unique(cols_amp$Condition))
    
    my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
    output$amphm <- renderPlotly({
      heatmaply::heatmaply(x = amp_mat |> dplyr::select(cols_amp$ID[order(cols_amp$Condition,decreasing = TRUE)]),
                           # colors = rev(colorRampPalette(brewer.pal(3, "RdBu"))(99)),
                           colors = my_palette,
                           col_side_colors = tibble::tibble("<b>Column Names</b>"= cols_amp$Condition),
                           col_side_palette = c("Sham" = "#999999",
                                                "Resection" = "magenta"),
                           hide_colorbar = FALSE,
                           showticklabels=TRUE,
                           Rowv = TRUE,
                           Colv = FALSE,
                           plot_method = "plotly",
                           scale = "none",
                           key.title = "Expression by\nGene Values",
                           subplot_widths=c(0.2, 0.03),
                           subplot_heights=c(0.03, 0.96),
                           hclust_method = "mcquitty",
                           grid_color = "white",
                           grid_width = 0.21)
      #heatmap_layers = theme(axis.line=element_blank()))
    }) 
    
    
  })
  
  
  output$ablated_net = renderVisNetwork({
    visNetwork::visIgraph(abl) |>
      visNodes(color = "#E69F00") |>
      visOptions(selectedBy = "Clusters") |>
      visInteraction(navigationButtons = TRUE) |>
      visNetwork::visLayout(randomSeed = 123456) |> 
      visEvents(select = "function(nodes) {Shiny.onInputChange('current_node_id', nodes.nodes);;}")
  })
  
  
  output$goTable_abl <- renderDataTable({
    enrich_go_bo_df %>%
      group_by(Group) %>%
      filter(Group == "Ablated") %>%
      filter(ID %in% vertex_attr(abl)$name[vertex_attr(abl)$name == input$current_node_id])
  })
  
  observeEvent(input$current_node_id, {
    visNetworkProxy("ablated_net") %>%
      visGetNodes()
  })
  
  output$verb_abl <- renderText({ paste0("Selected GO:BP is: ",vertex_attr(abl)$`EnrichmentMap::GS_DESCR`[vertex_attr(abl)$name == input$current_node_id])
  })
  
  output$genes_go_abl <- renderValueBox(valueBox(value = tags$p(enrich_go_bo_df %>% group_by(Group) %>% filter(Group == "Ablated") %>% filter(ID %in% vertex_attr(abl)$name[vertex_attr(abl)$name == input$current_node_id]) %>% pull(geneID) %>% gsub(pattern = "/",replacement=","),style = "font-size: 35%;", style ="color : black"), color="orange", subtitle = paste0("Genes Involved in the ",vertex_attr(abl)$`EnrichmentMap::GS_DESCR`[vertex_attr(abl)$name == input$current_node_id])))
  
  
  genes_abl <-  "Select_a_node"
  # print(genes_abl)
  
  target_cols_abl <- grep("^SRR", colnames(abl_counts), value = TRUE)
  abl_counts2 <- abl_counts[, c("MGI_Symbol", target_cols_abl), drop = FALSE]
  
  abl_counts_genes <- abl_counts2 %>% filter(MGI_Symbol %in% genes_abl) %>% distinct(MGI_Symbol,.keep_all = TRUE)
  
  abl_mat <- as.data.frame(abl_counts_genes[, -1])
  if(nrow(abl_mat) > 0) {
    rownames(abl_mat) <- abl_counts_genes$MGI_Symbol
  }
  
  cols_abl <- data.frame(ID = colnames(abl_mat))
  cols_abl$Condition <- ifelse(grepl("Sham", cols_abl$ID), "Sham", "Ablated")
  cols_abl <- cols_abl[order(cols_abl$Condition, decreasing = TRUE),]
  # print(unique(cols_abl$Condition))
  if(nrow(abl_mat)<=2){cluster=TRUE}else {cluster=FALSE}
  my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
  
  output$ablhm <- renderPlotly({
    heatmaply::heatmaply(x = abl_mat |> dplyr::select(cols_abl$ID[order(cols_abl$Condition,decreasing = TRUE)]),
                         
                         Rowv = FALSE,
                         Colv = FALSE,
                         
                         trace="none") 
  }) 
  
  observeEvent(input$current_node_id, {
    genes_abl <- enrich_go_bo_df %>% group_by(Group) %>% filter(Group == "Ablated") %>% filter(ID %in% vertex_attr(abl)$name[vertex_attr(abl)$name == input$current_node_id]) %>% pull(geneID) %>% gsub(pattern = "/",replacement=",")  %>% strsplit(split = ",") %>% unlist()
    
    target_cols_abl <- grep("^SRR", colnames(abl_counts), value = TRUE)
    abl_counts2 <- abl_counts[, c("MGI_Symbol", target_cols_abl), drop = FALSE]
    
    abl_counts_genes <- abl_counts2 %>% filter(MGI_Symbol %in% genes_abl) %>% distinct(MGI_Symbol,.keep_all = TRUE)
    
    abl_mat <- as.data.frame(abl_counts_genes[, -1])
    if(nrow(abl_mat) > 0) {
      rownames(abl_mat) <- abl_counts_genes$MGI_Symbol
    }
    
    cols_abl <- data.frame(ID = colnames(abl_mat))
    cols_abl$Condition <- ifelse(grepl("Sham", cols_abl$ID), "Sham", "Ablated")
    cols_abl <- cols_abl[order(cols_abl$Condition, decreasing = TRUE),]
    
    output$ablhm <-   
      renderPlotly({
        heatmaply::heatmaply(x = abl_mat |> dplyr::select(cols_abl$ID[order(cols_abl$Condition,decreasing = TRUE)]),
                             # colors = rev(colorRampPalette(brewer.pal(3, "RdBu"))(99)),
                             colors = my_palette,
                             col_side_colors = tibble::tibble("<b>Column Names</b>" = cols_abl$Condition),
                             col_side_palette = c("Sham" = "#999999",
                                                  "Ablated"= "orange"),
                             hide_colorbar = FALSE,
                             showticklabels=TRUE,
                             Rowv = TRUE,
                             Colv = FALSE,
                             plot_method = "plotly",
                             scale = "none",
                             key.title = "Expression by\nGene Values",
                             subplot_widths=c(0.2, 0.03),
                             subplot_heights=c(0.03, 0.965),
                             hclust_method = "mcquitty",
                             grid_color = "white",
                             grid_width = 0.21)
        #heatmap_layers = theme(axis.line=element_blank()))
      })
  })
  
  
  
  output$uninjured_net = renderVisNetwork({
    visNetwork::visIgraph(unj) |>
      visNodes(color = "#059E74") |>
      visOptions(selectedBy = "Clusters") |>
      visInteraction(navigationButtons = TRUE) |>
      visNetwork::visLayout(randomSeed = 123456) |>
      visEvents(select = "function(nodes) {Shiny.onInputChange('current_node_id', nodes.nodes);;}")
  })
  
  output$goTable_unj <- renderDataTable({
    enrich_go_bo_df %>%
      group_by(Group) %>%
      filter(Group == "Uninjured") %>%
      filter(ID %in% vertex_attr(unj)$name[vertex_attr(unj)$name == input$current_node_id])
  })
  
  observeEvent(input$current_node_id, {
    visNetworkProxy("uninjured_net") %>%
      visGetNodes()
  })
  
  output$verb_unj <- renderText({ paste0("Selected GO:BP is: ",vertex_attr(unj)$`EnrichmentMap::GS_DESCR`[vertex_attr(unj)$name == input$current_node_id])
  })
  
  output$genes_go_unj <- renderValueBox(valueBox(value = tags$p(enrich_go_bo_df %>% group_by(Group) %>% filter(Group == "Uninjured") %>% filter(ID %in% vertex_attr(unj)$name[vertex_attr(unj)$name == input$current_node_id]) %>% pull(geneID) %>% gsub(pattern = "/",replacement=","),style = "font-size: 35%;", style ="color : black"), color="green", subtitle =paste0("Genes Involved in the ",vertex_attr(unj)$`EnrichmentMap::GS_DESCR`[vertex_attr(unj)$name == input$current_node_id])))
  
  genes_unj <-  "Select_a_node"
  # print(genes_unj)
  
  target_cols_unj <- grep("^SRR", colnames(unj_counts), value = TRUE)
  unj_counts2 <- unj_counts[, c("MGI_Symbol", target_cols_unj), drop = FALSE]
  
  unj_counts_genes <- unj_counts2 %>% filter(MGI_Symbol %in% genes_unj) %>% distinct(MGI_Symbol,.keep_all = TRUE)
  
  unj_mat <- as.data.frame(unj_counts_genes[, -1])
  if(nrow(unj_mat) > 0) {
    rownames(unj_mat) <- unj_counts_genes$MGI_Symbol
  }
  
  cols_unj <- data.frame(ID = colnames(unj_mat))
  cols_unj$Condition <- ifelse(grepl("Sham", cols_unj$ID), "Sham", "Uninjured")
  cols_unj <- cols_unj[order(cols_unj$Condition, decreasing = TRUE),]
  # print(unique(cols_unj$Condition))
  if(nrow(unj_mat)<=2){cluster=TRUE}else {cluster=FALSE}
  my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
  
  output$unjhm <- renderPlotly({
    heatmaply::heatmaply(x = unj_mat |> dplyr::select(cols_unj$ID[order(cols_unj$Condition,decreasing = TRUE)]),
                         
                         Rowv = FALSE,
                         Colv = FALSE,
                         
                         trace="none") 
  }) 
  
  observeEvent(input$current_node_id, {
    genes_unj <- enrich_go_bo_df %>% group_by(Group) %>% filter(Group == "Uninjured") %>% filter(ID %in% vertex_attr(unj)$name[vertex_attr(unj)$name == input$current_node_id]) %>% pull(geneID) %>% gsub(pattern = "/",replacement=",")  %>% strsplit(split = ",") %>% unlist()
    
    target_cols_unj <- grep("^SRR", colnames(unj_counts), value = TRUE)
    unj_counts2 <- unj_counts[, c("MGI_Symbol", target_cols_unj), drop = FALSE]
    
    unj_counts_genes <- unj_counts2 %>% filter(MGI_Symbol %in% genes_unj) %>% distinct(MGI_Symbol,.keep_all = TRUE)
    
    unj_mat <- as.data.frame(unj_counts_genes[, -1])
    if(nrow(unj_mat) > 0) {
      rownames(unj_mat) <- unj_counts_genes$MGI_Symbol
    }
    
    cols_unj <- data.frame(ID = colnames(unj_mat))
    cols_unj$Condition <- ifelse(grepl("Sham", cols_unj$ID), "Sham", "Uninjured")
    cols_unj <- cols_unj[order(cols_unj$Condition, decreasing = FALSE),]
    
    output$unjhm <- renderPlotly({
      heatmaply::heatmaply(x = unj_mat |> dplyr::select(cols_unj$ID[order(cols_unj$Condition,decreasing = FALSE)]),
                           # colors = rev(colorRampPalette(brewer.pal(3, "RdBu"))(99)),
                           colors = my_palette,
                           col_side_colors = tibble::tibble("<b>Column Names</b>" = cols_unj$Condition),
                           col_side_palette = c("Sham" = "#999999",
                                                "Uninjured"= "#059E74"),
                           hide_colorbar = FALSE,
                           showticklabels=TRUE,
                           Rowv = TRUE,
                           Colv = FALSE,
                           scale = "none",
                           subplot_widths=c(0.2, 0.03),
                           subplot_heights=c(0.03, 0.96),
                           hclust_method = "mcquitty",
                           plot_method = "plotly",
                           key.title = "Expression by\nGene Values",
                           grid_color = "white",
                           grid_width = 0.21)
    })
  })
  
  
  output$heartcore_net = renderVisNetwork({
    visNetwork::visIgraph(core) |>
      visNodes(color = "#DCD159") |>
      visOptions(selectedBy = "Clusters") |>
      visInteraction(navigationButtons = TRUE) |>
      visNetwork::visLayout(randomSeed = 123456) |> 
      #visEvents(select = "function(nodes) {Shiny.onInputChange('current_node_id', nodes.nodes);;}")
      visEvents(select = "function(nodes) {Shiny.onInputChange('current_node_id', nodes.nodes);;}")
  })
  
  output$goTable_core <- renderDataTable({
    core_reg %>%
      group_by(ID) %>%
      filter(ID %in% vertex_attr(core)$name[vertex_attr(core)$name == input$current_node_id])
  })
  
  observeEvent(input$current_node_id, {
    visNetworkProxy("heartcore_net") %>%
      visGetNodes()
  })
  
  output$verb_core <- renderText({ paste0("Selected GO:BP is: ",vertex_attr(core)$`EnrichmentMap::GS_DESCR`[vertex_attr(core)$name == input$current_node_id])
  })

  
  output$genes_go_core <- renderValueBox(valueBox(value = tags$p(core_reg %>% group_by(ID) %>% filter(ID %in% vertex_attr(core)$name[vertex_attr(core)$name == input$current_node_id]) %>% pull(Genes)%>% gsub(pattern = "/",replacement=","),style = "font-size: 35%;", style ="color : black"), color="yellow", subtitle = paste0("Genes Involved in the ",vertex_attr(core)$`EnrichmentMap::GS_DESCR`[vertex_attr(core)$name == input$current_node_id])))
  
  
  genes_core <-  "Select_a_node"
  # print(genes_core)
  # print("head of s4c")
  # print(head(s4c))
  sample_ids <- gsub('"', '', s4c[,"ID"])
  # print(paste("colData_s4c sample_ids count:", length(sample_ids)))
  # print(sample_ids)
  keep_ids <- intersect(sample_ids, colnames(core_counts))
  # print(paste("core_counts sample columns kept:", length(keep_ids)))
  core_counts2 <- as.data.frame(core_counts[,c(keep_ids,"MGI_Symbol") , drop = FALSE])
  # print(paste("core_counts2 dim:", paste(dim(core_counts2), collapse = "x")))
  
  core_counts_genes <- core_counts2 %>% filter(MGI_Symbol %in% genes_core) %>% distinct(MGI_Symbol,.keep_all = TRUE)
  # print(paste("core_counts_genes head:"))
  # print(head(core_counts_genes))


  core_mat <- as.data.frame(core_counts_genes[,-c(length(core_counts_genes)) , drop = FALSE])
  rownames(core_mat) <- core_counts_genes$MGI_Symbol
  
  # print(paste("core_mat dim:", paste(dim(core_mat), collapse = "x")))
  # print("s4c head:")
  # print(head(s4c))
  # print(unique(s4c$Condition))
  cols_core <- s4c
  cols_core$ID <- rownames(cols_core)
  cols_core <- cols_core |> dplyr::filter(ID %in% colnames(core_mat))
  # print(unique(cols_core$Condition))
  
  # print(unique(cols_core$Condition))
  cols_core <- cols_core[order(cols_core$Condition,decreasing = TRUE),]
  # print(unique(cols_core$Condition))
  # print(head(core_mat))
  # print("test")
  # print(head(core_mat |> dplyr::select(cols_core$ID[order(cols_core$Condition,decreasing = TRUE)])))
  
  if(nrow(core_mat)<=2){cluster=TRUE}else {cluster=FALSE}
  my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
  
  output$hchm <- renderPlotly({
    heatmaply::heatmaply(x = core_mat |> dplyr::select(cols_core$ID[order(cols_core$Condition,decreasing = TRUE)]),
                         
                         Rowv = FALSE,
                         Colv = FALSE,
                         
                         trace="none")
  }) 
  
  
  #"GO:0071621"
  observeEvent(input$current_node_id, {
    genes_core <- core_reg %>% group_by(ID) %>% filter(ID %in% vertex_attr(core)$name[vertex_attr(core)$name == input$current_node_id]) %>% pull(Genes) %>% strsplit(split = ",") %>% unlist()
    #names(core_counts)
    core_counts2 <- as.data.frame(core_counts[,c(rownames(s4c %>% filter(Condition %in% c("Sham","Cryoinjury","Ablated","Resection"))),"MGI_Symbol")])
    core_counts_genes <- core_counts2 %>% filter(MGI_Symbol %in% genes_core) %>% distinct(MGI_Symbol,.keep_all = TRUE)
    core_mat <- as.data.frame(core_counts_genes[,-c(length(core_counts_genes))])
    rownames(core_mat) <- core_counts_genes$MGI_Symbol
    
    #col_annotation_core <- data.frame("Conditions" = s4c |> filter(rownames(s4c) %in% colnames(core_mat)) |> pull(Condition))
    
    # core_mat reaarrange by name and not by clustering in all data frames heatmaps
    #s4c |> dplyr::filter(rownames(s4c) %in% colnames(core_mat)) |> dplyr::pull(order(Condition))
    cols_core <- s4c |> dplyr::filter(rownames(s4c) %in% colnames(core_mat)) 
    cols_core$ID <- row.names(cols_core)
    cols_core <- cols_core[order(cols_core$Condition,decreasing = TRUE),]
    
    output$hchm <- renderPlotly({
      # heatmaply::heatmaply(x= core_mat[cols_core$ID[order(cols_core$Condition,decreasing = TRUE)]],
      heatmaply::heatmaply(x= core_mat |> dplyr::select(cols_core$ID[order(cols_core$Condition,decreasing = TRUE)]),
                           colors = my_palette,
                           col_side_colors = tibble::tibble("<b>Column Names</b>" = cols_core$Condition),
                           col_side_palette = c("Sham" = "#999999",
                                                "Ablated"= "#E69F00",
                                                "Resection" = "#CC6699",
                                                "Cryoinjury" = "#0673B3"),
                           hide_colorbar = FALSE,
                           showticklabels=TRUE,
                           Rowv = TRUE,
                           Colv = FALSE,
                           scale = "none",
                           subplot_widths=c(0.4, 0.03),
                           # subplot_heights=c(0.1,0.03, 0.8),
                           subplot_heights=c(0.04,0.95),
                           hclust_method = "mcquitty",
                           plot_method = "plotly",
                           key.title = "Expression by\nGene Values",
                           grid_color = "white",
                           grid_width = 0.21)#heatmap_layers = theme(axis.line=element_blank()))
    })
    
  })
  
  #Create Gene symbol list to select and display
  list_input_data_inputize <- gene_seeker_table |> group_by(MGI_Symbol) |> pull(MGI_Symbol)
  list_input_data_inputize <- list(list_input_data_inputize)
  updateSelectizeInput(session, 'gene', choices = list_input_data_inputize, server = TRUE)
  #col_annotation_core <- data.frame("Conditions" = s4c |> filter(rownames(s4c) %in% colnames(gene_seeker_table)) |> pull(Condition))
  
  
  observeEvent(input$gene,{
    #data table
    output$geneseeker_table <- DT::renderDataTable({
      a_new_table <- gene_seeker_table |> dplyr::group_by(MGI_Symbol) |>  filter(MGI_Symbol %in% input$gene)#c("Foxm1","Sypl2"))
      datatable(a_new_table,options(scrollX = TRUE))
      #gene_seeker_table |> dplyr::group_by(MGI_Symbol) |> filter(MGI_Symbol %in% "Foxm1") filter(MGI_Symbol %in% input$gene) |> distinct(.keep_all = TRUE)
    })
  })
  
  
  
  
  observeEvent(input$gene,{
    
    #genes_selected <- gene_seeker_table |> dplyr::group_by(MGI_Symbol) |>  filter(MGI_Symbol %in% input$gene)
    #names(core_counts)
    #dim(core_counts2)
    core_counts2 <- as.data.frame(gene_seeker_table[,c(rownames(s4c %>% filter(Condition %in% c("Sham","Cryoinjury","Ablated","Resection","Uninjured"))),"MGI_Symbol")])
    core_counts_genes <- core_counts2 %>% filter(MGI_Symbol %in% input$gene) %>% distinct(MGI_Symbol,.keep_all = TRUE)
    core_mat_input <- as.data.frame(core_counts_genes[,-c(length(core_counts_genes))])
    rownames(core_mat_input) <- core_counts_genes$MGI_Symbol
    #col_annotation_core_input <- data.frame("Conditions" = s4c |> filter(rownames(s4c) %in% colnames(core_mat_input)) |> pull(Condition))
    
    #labels_custom_hover <- genes_148__desc_mmus |> dplyr::filter(mgi_symbol %in% c("Adra2b","Bub1")) |> dplyr::pull(Description) 
    
    
    cols_core2 <- s4c |> dplyr::filter(rownames(s4c) %in% colnames(core_mat_input)) 
    cols_core2$ID <- row.names(cols_core2)
    cols_core2 <- cols_core2[order(cols_core2$Condition,decreasing = TRUE),]
    
    
    my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
    output$geneseeker_plot <- renderPlotly({
      #heatmaply::heatmaply(x = normalize(core_mat_input |> dplyr::select(cols_core2$ID[order(cols_core2$Condition,decreasing = TRUE)])),                   
      heatmaply::heatmaply(x = core_mat_input |> dplyr::select(cols_core2$ID[order(cols_core2$Condition,decreasing = TRUE)]),
                           colors = my_palette,
                           col_side_colors = tibble::tibble("<b>Column Names</b>" = cols_core2$Condition),
                           col_side_palette = c("Ablated" = "#E69F00",
                                                # "Uninjured"= "#059E74",
                                                "Sham" = "#999999",
                                                "Cryoinjury" = "#0673B3",
                                                "Resection" = "#CC6699"),
                           hide_colorbar = FALSE,
                           showticklabels = TRUE,
                           Rowv = TRUE,
                           Colv = FALSE,
                           lot_method = c("plotly"),
                           scale = "none",
                           subplot_widths=c(0.6, 0.3),
                           subplot_heights=c(0.04,0.95),
                           hclust_method = "mcquitty",
                           key.title = "Expression by\nGene Values",
                           fontsize_col = 10,
                           grid_color = "white",
                           #custom_hovertext = labels_custom_hover,
                           grid_width = 0.21)#heatmap_layers = theme(axis.line=element_blank()))
    })
  })
  
  my_palette_v <- rev(brewer.pal(n = 11, name = "RdYlBu"))
  
  
  output$core_l2fc_hm <- renderPlotly({
    heatmaply::heatmaply(x= counts_148_selected_without_doubles |> dplyr::select(order(colnames(counts_148_selected_without_doubles))),
                         scale_fill_gradient_fun = scale_fill_gradientn(
                           colours = c("darkblue", "blue", "white", "yellow","orange","red","darkred"),
                           limits = c(-5,15),
                           values = scales::rescale(c(-5,-2,0,2,5,10,15))),
                         hide_colorbar = FALSE,
                         showticklabels = TRUE,
                         Rowv = TRUE,
                         Colv = FALSE,
                         #plot_method = "plotly",
                         scale = "none",
                         key.title = "Expression by\nLog 2 Fold Change",
                         #grid_color = "white",
                         # fontsize_col = 12,
                         fontsize_row = 7,
                         subplot_widths=c(0.3, 0.2),
                         subplot_heights=c(1),
                         hclust_method = "mcquitty",
                         # showlegened = TRUE,
                         # theme(
                         #   axis.line = element_line(color = "black"),
                         #   panel.grid.major = element_blank(),
                         #   panel.grid.minor = element_blank(),
                         #   panel.border = element_blank(),
                         #   panel.background = element_blank()
                         # ),
                         main = "Core Regeneration Genes Heatmap of log2-Fold-Change values")
  })
  
  #Pubmed query results
  unj_fil <- read.delim(file = "./data/Platform_DEGs_only_Uninjured_vs_Sham_expression.txt",header = TRUE,sep = "\t")
  
  pubmed_positive <- read.csv(file = "./data/gene_hit_summary.csv",header = TRUE) %>% dplyr::filter(Has_Hits == TRUE) %>% pull(Gene)
  pubmed_positive
  pubmed_negative <- read.csv(file = "./data/gene_hit_summary.csv",header = TRUE) %>% dplyr::filter(Has_Hits == FALSE) %>% pull(Gene)
  pubmed_negative
  
  pubmed_positive_f <- core_counts |>  dplyr::filter(MGI_Symbol %in% pubmed_positive) |>  dplyr::distinct(MGI_Symbol,.keep_all = TRUE)
  rownames(pubmed_positive_f) <- pubmed_positive_f$MGI_Symbol
  pubmed_positive_f <- pubmed_positive_f[,-c(1,2)]
  # print("pubmed positive head")
  # print(head(pubmed_positive_f))
  
  pmp <- s4c |> dplyr::filter(rownames(s4c) %in% colnames(pubmed_positive_f)) 
  pmp$ID <- row.names(pmp)
  pmp <- pmp[order(pmp$Condition,decreasing = TRUE),]
  
  
  
  pubmed_negative_f <- core_counts |>  dplyr::filter(MGI_Symbol %in% pubmed_negative) |>  dplyr::distinct(MGI_Symbol,.keep_all = TRUE)

  rownames(pubmed_negative_f) <- pubmed_negative_f$MGI_Symbol
  pubmed_negative_f <- pubmed_negative_f[,-c(1,2)]
  #   print("pubmed negative head")
  # print(head(pubmed_negative_f))
  
  pmn <- s4c |> dplyr::filter(rownames(s4c) %in% colnames(pubmed_negative_f)) 
  pmn$ID <- row.names(pmn)
  pmn <- pmn[order(pmn$Condition,decreasing = TRUE),]
  
  
  
  my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))(256)
  
output$known <- renderPlotly({
    # 1. Filter metadata
    pmp <- s4c |> dplyr::filter(rownames(s4c) %in% colnames(pubmed_positive_f)) 
    pmp$ID <- row.names(pmp)
    
    # 2. CRITICAL FIX: Convert to character to remove ghost factor levels
    pmp$Condition <- as.character(pmp$Condition)
    
    # 3. FIX MISMATCH: Rename 'Amputation' to 'Resection' to match the palette
    pmp$Condition[pmp$Condition == "Amputation"] <- "Resection"
    
    # 4. Sort
    pmp <- pmp[order(pmp$Condition, decreasing = TRUE),]
    
    # 5. Debug
    # print("Conditions present in 'Known' dataset:")
    # print(unique(pmp$Condition))
    
    # 6. Prepare numeric matrix
    pubmed_positive_numeric <- pubmed_positive_f %>% 
      dplyr::select(starts_with("SRR"))
    
    heatmaply::heatmaply(
      x = pubmed_positive_numeric[, pmp$ID], # pmp is sorted
      colors = my_palette,
      col_side_colors = tibble::tibble("<b>Column Names</b>" = pmp$Condition),
      col_side_palette = c(
        "Sham" = "#999999",
        "Ablated" = "#E69F00",
        "Resection" = "#CC79A7",
        "Cryoinjury" = "#0072B2",
        "Uninjured" = "#009E73"
      ),
      hide_colorbar = FALSE,
      showticklabels = TRUE,
      Rowv = TRUE,
      Colv = FALSE,
      plot_method = "plotly",
      scale = "none",
      subplot_widths = c(0.6, 0.3),
      subplot_heights = c(0.04, 0.95),
      key.title = "Expression by\nGene Values",
      grid_color = "white",
      grid_width = 0.71,
      hclust_method = "mcquitty",
      main = "Previously Regeneration Described Genes"
    )
  })

  output$unknown <- renderPlotly({
    # 1. Filter metadata
    pmn <- s4c |> dplyr::filter(rownames(s4c) %in% colnames(pubmed_negative_f)) 
    pmn$ID <- row.names(pmn)
    
    # 2. CRITICAL FIX: Convert to character
    pmn$Condition <- as.character(pmn$Condition)
    
    # 3. FIX MISMATCH: Rename 'Amputation' to 'Resection'
    pmn$Condition[pmn$Condition == "Amputation"] <- "Resection"
    
    # 4. Sort
    pmn <- pmn[order(pmn$Condition, decreasing = TRUE),]
    
    # 5. Debug
    # print("Conditions present in 'Unknown' dataset:")
    # print(unique(pmn$Condition))
    
    # 6. Prepare numeric matrix
    pubmed_negative_numeric <- pubmed_negative_f %>% 
      dplyr::select(starts_with("SRR"))
      
    heatmaply::heatmaply(
      x = pubmed_negative_numeric[, pmn$ID], # pmn is sorted
      colors = my_palette,
      col_side_colors = tibble::tibble("<b>Column Names</b>" = pmn$Condition),
      col_side_palette = c(
        "Sham" = "#999999",
        "Ablated" = "#E69F00",
        "Resection" = "#CC79A7",
        "Cryoinjury" = "#0072B2",
        "Uninjured" = "#009E73"
      ),
      hide_colorbar = FALSE,
      showticklabels = TRUE,
      Rowv = TRUE,
      Colv = FALSE,
      plot_method = "plotly",
      scale = "none",
      subplot_widths = c(0.6, 0.3),
      subplot_heights = c(0.04, 0.95),
      key.title = "Expression by\nGene Values",
      grid_color = "white",
      grid_width = 0.71,
      hclust_method = "mcquitty",
      main = "Novel Regeneration Genes"
    )
  })

  
  
  
  
  
}