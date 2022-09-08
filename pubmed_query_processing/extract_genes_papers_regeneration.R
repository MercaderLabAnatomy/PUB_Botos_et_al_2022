#install.packages("easyPubMed")
library("easyPubMed")
library("dplyr")
#Easy pub med
#https://cran.r-project.org/web/packages/easyPubMed/vignettes/getting_started_with_easyPubMed.html


#xy <- read.table("C:/Users/marius/Documents/Phd/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/R/core_reg_abl_vs_sham_genes.txt",sep="\t",header = TRUE)
core_counts <- read.table("/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/ResultsCytoscape/PossibleShiny/RShinyCode/ShinnyNetworksGOBP/data/core_selected_counts_148_Table.txt",sep = "\t",header = TRUE)
s4c <- read.table("/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/ResultsCytoscape/PossibleShiny/RShinyCode/ShinnyNetworksGOBP/data/colData_s4c.txt",sep = "\t",header = TRUE)
xy <- read.table("/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/coreRegeneration/core_reg_abl_vs_sham_genes.txt",sep="\t",header = TRUE)
length(xy$MGI_Symbol)
# length(xy$MGI_Symbol)
# my_gene_list <- list()
# 
# my_gene_list <- lapply(seq_along(1:length(xy$MGI_Symbol)),function(i) i)
# for (j in length(xy$MGI_Symbol)){
#   my_gene_list[[j]] <- xy$MGI_Symbol[j]}




#[TIAB] = title and abstract

#my_gene_list <- lapply(seq_along(1:length(xy$MGI_Symbol)),function(i) paste0("((",xy$MGI_Symbol[i],"[TIAB]) AND (regeneration[TIAB])) OR ((",xy$MGI_Symbol[i],"[TIAB]) AND (zebrafish[TIAB]))"))
my_gene_list <- lapply(seq_along(1:length(xy$MGI_Symbol)),function(i) paste0("((",xy$MGI_Symbol[i],"[TIAB]) AND (regeneration[TIAB]))"))
names(my_gene_list) <- xy$MGI_Symbol
my_gene_list



#Query pubmed for the string up used
out_pub_meds <- lapply(seq_along(1:length(my_gene_list)),function(n) {
  out_a <- batch_pubmed_download(pubmed_query_string = my_gene_list[[n]],dest_file_prefix = names(my_gene_list)[n], 
                                 encoding = "UTF-8",batch_size = 100000,format = "xml")})



#Clean out genes that do not have any return with the words combination used
out_b <- out_pub_meds[!unlist(lapply(out_pub_meds,is.null))]
names(out_b) <- out_b



#Retrieve the information for the ones positive
xx <- lapply(seq_along(1:length(out_b)),function(i) {
  table_articles_byAuth(pubmed_data = out_b[[i]], autofill = TRUE,included_authors = "all",encoding = "UTF-8")})

names(xx) <- out_b
#lapply(xx,head)


#Create a data list with the lenght of papers found per gene, the papers id and the gene used for it.
genes_df_pmid <- lapply(seq_along(1:length(out_b)),function(i) i)
#genes_df_pmid <- list()

for (i in seq_along(1:length(xx))){
  genes_df_pmid[[i]][["Count_of_PMIDS"]] <- length(unique(xx[[i]][["pmid"]]))
  genes_df_pmid[[i]][["PMIDs"]] <- paste(unique(xx[[i]][["pmid"]]),collapse = "_")
  genes_df_pmid[[i]][["Gene"]] <- gsub(x = names(xx[i]),pattern = "01.txt",replacement = "")
}


#Format the data.
lapply(genes_df_pmid,head)

test <- data.frame()
test <- do.call(cbind,lapply(genes_df_pmid,as.data.frame))

write.table(x = t(test),
            file = "/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/genes_counts_pubmed_in_regeneration.txt",
            #file = "C:/Users/marius/Documents/Phd/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/R/genes_in_reg.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

genes_returned_result <- list()
for (i in seq_along(1:length(test))){
  genes_returned_result[length(genes_returned_result) + 1] <- (test[[i]][[4]])}

genes_returned_result_df <- as.data.frame(unlist(genes_returned_result))
colnames(genes_returned_result_df) <- "Genes"
xy_genes <- as.data.frame(unlist(xy$MGI_Symbol))
colnames(xy_genes) <- "MainGenes"
xy_genes %>% filter(MainGenes %in%  genes_returned_result_df$Genes) %>% write.table("/media/marius/Samsung_T5/Phd/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/genes_in_core_reg_positive_result_query_pubmed.txt",col.names = TRUE,row.names = FALSE,sep = "\t")
xy_genes %>% filter(!MainGenes %in%  genes_returned_result_df$Genes) %>% write.table("/media/marius/Samsung_T5/Phd/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/genes_in_core_reg_no_result_query_pubmed.txt",col.names = TRUE,row.names = FALSE,sep = "\t")



#add description to the genes not known in regeneration
library("biomaRt")
library("dplyr")
library("curl")

#biomaRt::listEnsembl()
#Use ensembl
ensembl <- useEnsembl("ENSEMBL_MART_ENSEMBL")
#biomaRt::listAttributes(ensembl)
#biomaRt::listDatasets(ensembl)

#dre = useEnsembl(biomart = "ensembl", dataset = "drerio_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")#mirror = "useast") # query organism (from this organism)
mmus = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")#mirror = "useast") # target organism (convert genes to this organism)

#danio_ensembl_id_attributes <-c("ensembl_gene_id", "mmusculus_homolog_perc_id")#,# attributes to query from danio database (query organism) 
#"mmusculus_homolog_goc_score") #orthology score to mouse from danio side

#danio_zfin_id_attributes <-c("ensembl_gene_id", "entrezgene_id", "zfin_id_symbol","description") # attributes to convert ensembl id to zfin id 


mmus_attributes <- c("mgi_symbol", "description", "entrezgene_id") # attributes that you want from mouse database

#make a converted dataframe with ensemble ids and mouse orthology scores 

#make dataframe for dre ensembl id and zfin
m <- xy_genes %>% filter(!MainGenes %in%  genes_returned_result_df$Genes)
dim(m)
genes_desc_entrez_mmus_not_known = getBM(attributes = mmus_attributes,
                                         filters = "mgi_symbol",
                                         values = m$MainGenes,
                                         mart = mmus)
dim(genes_desc_entrez_mmus_not_known)

write.table(x = genes_desc_entrez_mmus_not_known,file = "/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/pubmed_core_reg_unknown_genes_description_entrez.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)



#Pie chart

library(ggplot2)

neg_res_pub <- xy_genes %>% filter(!MainGenes %in%  genes_returned_result_df$Genes)
pos_res_pub <- xy_genes %>% filter(MainGenes %in%  genes_returned_result_df$Genes)

df_pc = data.frame(x <- c(length(neg_res_pub$MainGenes),
                          length(pos_res_pub$MainGenes)),
                   labels <- c('Unknown','Known'),
                   value=c(length(neg_res_pub$MainGenes),
                           length(pos_res_pub$MainGenes)))

pc_res <- ggplot(df_pc, aes(x="", y=x, fill=labels)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5)) +
  labs(fill = ("Studied genes")) +
  theme_void()

#pc_res
ggsave(plot = pc_res,
       filename = "/media/marius/Samsung_T5/Phd/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/genes_in_core_reg_result_query_pubmed.pdf",height = 8,width = 12,dpi = 300)




#heatmaply

library("RColorBrewer")
library("heatmaply")
heatmap_it_core <- function(result_df,core_regeneration_counts_df=core_counts,name) {
  #Do a heatmap of the genes found
  a_df <- core_regeneration_counts_df %>% filter(MGI_Symbol %in% result_df$MainGenes) %>% distinct(MGI_Symbol,.keep_all = TRUE)
  #dim(a_df)
  col_annotation_core <- data.frame("Conditions" = s4c |> filter(rownames(s4c) %in% colnames(a_df)) |> pull(Condition))
  rownames(a_df) <- a_df$MGI_Symbol
  a_df <- a_df[,-c(1)]
  
  heatmaply::heatmaply(x=a_df,
                       colors = rev(colorRampPalette(brewer.pal(3, "RdBu"))(2560)),
                       col_side_colors = col_annotation_core$Conditions,
                       col_side_palette = c("Sham" = "#999999","Ablated"= "#E69F00","Amputation" = "#CC6699","Cryoinjury" = "#0673B3"),
                       hide_colorbar = FALSE,
                       showticklabels = TRUE,
                       Rowv = TRUE,
                       Colv = TRUE,
                       plot_method = "plotly",
                       scale = "column",
                       key.title = "Expression by\nColumn",
                       grid_color = "white",
                       grid_width = 0.71,
                       hclust_method = "mcquitty",
                       main = name[1],
                       file = paste0("/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/heatmaply/",name,"genes_regeneration_heatmaply.html"))
}

dim(neg_res_pub)
neg_res_pub_hm <- heatmap_it_core(result_df = neg_res_pub,core_regeneration_counts_df = core_counts,name = "Unknown")
neg_res_pub_hm
#htmlwidgets::saveWidget(widget = neg_res_pub_hm,file = "/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/PubMed_Query/heatmaply/heatmaply_Unknow_plot.html")


dim(pos_res_pub)
post_res_pub_hm <- heatmap_it_core(result_df = pos_res_pub,core_regeneration_counts_df = core_counts,name = "Known")
post_res_pub_hm
