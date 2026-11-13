# List of required packages
required_packages <- c(
  "rlang",
  "sodium",
  "devtools",
  "shinymanager",
  "DT",
  "dplyr",
  "fontawesome",
  "shinydashboard",
  "shinycssloaders",
  "igraph",
  "rmarkdown",
  "markdown",
  "ggplot2",
  "visNetwork",
  "heatmaply",
  "plotly",
  "brew",
  "RColorBrewer",
  "curl",
  "wesanderson",
  "shinyauthr"
)

# Check which packages are not installed
packages_to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

# Install missing packages
if (length(packages_to_install) > 0) {
  cat("Installing missing packages:", paste(packages_to_install, collapse = ", "), "\n")
  install.packages(packages_to_install, dependencies = TRUE)
  cat("Installation complete!\n")
} else {
  cat("All required packages are already installed.\n")
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

cat("All packages loaded successfully.\n")

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

# print("Starting ui.R execution...")

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
}
idleTimer();"





myCol <- c("#CC79A7","#0072B2","#E69F00","#999999")
myCol_3 <- c("#CC79A7","#0072B2","#009E73","#999999")

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




core_counts <- read.delim("./data/core_regeneration_genes.txt",sep = "\t",header = TRUE)
# print(paste("[ui] core_regeneration_genes rows/cols:", nrow(core_counts), ncol(core_counts)))
core_counts[nrow(core_counts)+1,] <- 0
core_counts[nrow(core_counts)+1,] <- 0
core_counts[nrow(core_counts),"MGI_Symbol"] <- "Select_a_node"
core_counts[nrow(core_counts)-1,"MGI_Symbol"] <- "Select_a_node"
# print("[ui] core_regeneration_genes added Select_a_node rows")


enrich_go_bo_df <- read.delim("./data/Platform_unique_Ablated_vs_Sham_GO.txt", sep = "\t", header = TRUE)
go_amp <- read.delim("./data/Platform_unique_Resection_vs_Sham_GO.txt", sep = "\t", header = TRUE)
go_unj <- read.delim("./data/Platform_unique_Uninjured_vs_Sham_GO.txt", sep = "\t", header = TRUE)

enrich_go_bo_df$Group <- "Ablated"
go_amp$Group <- "Resection"
go_unj$Group <- "Uninjured"

enrich_go_bo_df <- rbind(enrich_go_bo_df, go_amp, go_unj)


core_reg <- read.delim(file = "./data/core_regeneration_GOBP.txt",sep = "\t",header = TRUE)
# print(paste("[ui] core_regeneration_GOBP rows/cols:", nrow(core_reg), ncol(core_reg)))


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


#create df
# print("Processing Cryoinjury...")
cryo_df <- tables_counts_stats$Cryoinjury_vs_Sham 
# print(paste("Cryo table dim:", paste(dim(cryo_df), collapse="x")))
cryo_df <- cryo_df %>% dplyr::filter(MGI_Symbol %in% core_counts$MGI_Symbol)
# print(paste("After filter:", paste(dim(cryo_df), collapse="x")))
cryo_df <- cryo_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))

# print("Processing Resection...")
amp_df <- tables_counts_stats$Resection_vs_Sham %>% dplyr::filter(MGI_Symbol %in% core_counts$MGI_Symbol)
amp_df <- amp_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))

# print("Processing Ablated...")
abl_df <- tables_counts_stats$Ablated_vs_Sham %>% dplyr::filter(MGI_Symbol %in% core_counts$MGI_Symbol)
abl_df <- abl_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))



abl_df_2 <- abl_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)
amp_df_2 <- amp_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)
cryo_df_2 <- cryo_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)




# How many log2FC values remain
# print("Number of log2FoldChange values in each dataset:")
# sapply(list(amp=amp_df_2$log2FoldChange,
            # abl=abl_df_2$log2FoldChange,
            # cryo=cryo_df_2$log2FoldChange), length) %>% print()

# Are log2FoldChange columns present?
# print("Checking presence of log2FoldChange column in each dataset...")
# sapply(list(Resection=amp_df, Ablated=abl_df, Cryoinjury=cryo_df), function(x) "log2FoldChange" %in% names(x)) %>% print()


#rm(counts_148_selected_without_doubles)
all.equal(cryo_df_2$MGI_Symbol,abl_df_2$MGI_Symbol)
counts_148_selected_without_doubles <- as.data.frame(amp_df_2$log2FoldChange)
counts_148_selected_without_doubles <- cbind.data.frame(counts_148_selected_without_doubles,abl_df_2$log2FoldChange)
counts_148_selected_without_doubles <- cbind.data.frame(counts_148_selected_without_doubles,cryo_df_2$log2FoldChange)

rownames(counts_148_selected_without_doubles) <- amp_df_2$MGI_Symbol
colnames(counts_148_selected_without_doubles) <- c("Resection","Ablation","Cryoinjury")


# ui <- secure_app(head_auth = tags$script(inactivity),
ui <-   dashboardPage(skin = "purple",
                               dashboardHeader(title = "Heart Injuries"),
                               dashboardSidebar(width = 300,
                                                sidebarMenu(
                                                  menuItem("Home",tabName = "home",icon = icon("home")),
                                                  menuItem("Resection vs Sham",tabName = "resected", icon = icon("scissors",lib = "glyphicon")),#,icon = icon("calendar",lib = "glyphicon")),
                                                  menuItem("Ablation vs Sham",tabName = "ablated",icon = icon("tag",lib = "glyphicon")),
                                                  menuItem("Uninjured vs Sham",tabName = "uninjured",icon = icon("ok",lib = "glyphicon")),
                                                  menuItem("Heart Core Regeneration Common Genes ",tabName = "heartcore",icon = icon("heart-empty",lib = "glyphicon")),
                                                  menuItem("PubMed Query",tabName = "kuk",icon = icon("list",lib = "glyphicon")),
                                                  menuItem("Core Regeneration Log2FC",tabName = "core_l2fc",icon = icon("bookmark",lib = "glyphicon")),
                                                  menuItem("Genes Seeker",tabName = "gs",icon = icon("search",lib = "glyphicon"))
                                                )
                                                ),
                               dashboardBody(
                                 tabItems(
                                   tabItem(tabName = "home",
                                           #home section and markdown
                                           includeMarkdown("./www/markdown_home.md")),
                                   # #next tab
                                   
                                   #https://stackoverflow.com/questions/70689513/how-to-have-shiny-dashboard-box-fill-the-entire-width-of-my-dashbaord
                                   #next tab
                                   tabItem(tabName = "resected",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_amp")),#tags$head(tags$style("#text1{color: red;font-size: 20px;font-style: italic;}")) %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "resected_net",width = "auto", height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_amp") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = "100%",outputId = "genes_go_amp") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "amphm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "ablated",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_abl")),#tags$head(tags$style("#text1{color: red;font-size: 20px;font-style: italic;}")) %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "ablated_net",width = "auto",height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_abl") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = 12,outputId = "genes_go_abl") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "ablhm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "uninjured",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_unj")),
                                             column(width = 12, visNetworkOutput(outputId = "uninjured_net",width = "auto",height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_unj") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = 12,outputId = "genes_go_unj") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "unjhm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "heartcore",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_core") %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "heartcore_net",width = "auto",height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_core") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = "auto",outputId = "genes_go_core") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "hchm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "kuk",
                                           fluidRow(
                                             column(width = 12,plotlyOutput(outputId = "known", width = "100%", height = "1024px") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "unknown", width = "100%", height = "1024px") %>% withSpinner())         
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "core_l2fc",
                                           fluidRow(
                                             column(width = 12,plotlyOutput(outputId = "core_l2fc_hm", width = "100%", height = "2048px") %>% withSpinner())         
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "gs",
                                           fluidRow(
                                             column(width=4,selectizeInput(inputId = "gene",
                                                                           choices = NULL,
                                                                           label = "Select at least 2 Gene Symbol",
                                                                           selected = c("Cst8", "Ms4a15"),
                                                                           multiple = TRUE)),
                                             column(width = 12, dataTableOutput(outputId = "geneseeker_table") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "geneseeker_plot",width = "100%",height = "100%") %>% withSpinner()),
                                           )
                                   )
                                 )# end tabItems
                                 # Add download button 
                               )# end dashboardBody
                 )# end dashboardPage
# )#end all
