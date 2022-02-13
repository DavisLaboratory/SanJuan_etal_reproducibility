# last used with R 3.6.0
# ensembl <- read.delim("christine_data/ensembl_results.txt", header = TRUE)

ensembl <- read.delim("ensembl_results.txt", header = TRUE)

# top 5 transcription factors
TFs <- ensembl[ensembl$isTF == 1,]
TFs[order(TFs$adj.P.Val, decreasing = FALSE),]
TFs[grep("FOX", TFs$Symbol),]

# regulators are the ones that are differentially expressed or is a TF picked by CCA and pamr as additional (supporting) evidence
target_regulators <- ensembl[ensembl$P.Value < 0.05,]
target_regulators <- ensembl[(ensembl$adj.P.Val < 0.05) | (ensembl$isTF ==1 & (ensembl$CCA ==1 & ensembl$pamr == 1)),]

ensembl2 <- target_regulators[,grep("enrichedMotif|score", colnames(target_regulators), invert = TRUE)]
ensembl2$isDE <- ifelse(ensembl2$P.Value < 0.05, 1,0)


# add information about druggability ---- 
# Re-visit ranking based on DrugBank and PharmaGKB --------
# genome_path <- "/stornext/General/data/academic/lab_davis/genomes/PharmGKB_062019"
# pharma_relationships <- read.delim(paste0(genome_path,"/relationships/relationships.tsv"), sep="\t",header = TRUE)
# pharma_genes <- read.delim(paste0(genome_path,"/genes/genes.tsv"), sep="\t",header = TRUE)


# DrugBank ----
library(dbparser)
library(dplyr)
library(ggplot2)
library(XML)

## parse data from XML and save it to memory
get_xml_db_rows("drugbank_full_database.xml")

## load drugs data
# drugs <- parse_drug()
# 
# ## load drug groups data
# drug_groups <- parse_drug_groups()
# 
# ## load drug targets actions data
# drug_targets_actions <- parse_drug_targets_actions()

drug_targets <- parse_drug_targets()
drug_targets[grep("Androgen receptor", drug_targets$name),]


all_target_ids <- read.csv("all.csv")

all_target_ids <- all_target_ids[match(ensembl2$Symbol,
                                       all_target_ids$Gene.Name),]



all_target_ids <- all_target_ids[!is.na(all_target_ids$Gene.Name),]

drugbank_dbparser_target_ids <- gsub("BE0+([1-9]+)","\\1", drug_targets$id)

table(all_target_ids %in% drugbank_dbparser_target_ids)

drug_targets$target_id <- gsub("BE0+([1-9]+)","\\1", drug_targets$id)

all_targets <- read.csv("all.csv")
drug_targets$Symbol <- all_targets$Gene.Name[match(drug_targets$target_id, all_targets$ID)]


sum(tolower(ensembl2$description) %in% tolower(drug_targets$name))
#ensembl2[ensembl2$description %in% tolower(drug_targets$name),]

ensembl2$druggable <- 0
ensembl2$druggable[ensembl2$Symbol %in% drug_targets$Symbol] <- 1

# sum(ensembl2$GeneID %in% pharma_genes$NCBI.Gene.ID)

# ensembl2_pharmGKB_id <- match(pharma_genes$PharmGKB.Accession.Id[match(ensembl2$GeneID, pharma_genes$NCBI.Gene.ID)], 
#                               pharma_relationships$Entity1_id)


# isInPharma <- which(!is.na(ensembl2_pharmGKB_id)) 
# ensembl2$PharmaGKB_association <- NA
# ensembl2$PharmaGKB_association[isInPharma] <- as.character(pharma_relationships$Association[ensembl2_pharmGKB_id[!is.na(ensembl2_pharmGKB_id)]])
# ensembl2$druggability <- ifelse(ensembl2$PharmaGKB_association %in% "associated",1,0)

#ensembl2$score_2 <- rowSums(ensembl2[,c("isDE","pamr","CCA","isTF","druggability")])
ensembl2$score <- rowSums(ensembl2[,c("isDE","pamr","CCA","isTF","druggable")])
ensembl2 <- ensembl2[order(ensembl2$score, ensembl2$logFC, decreasing = TRUE), ]
which(ensembl2$Symbol=="AR")
which(ensembl2$Symbol=="GATA1")
ensembl2[1:55,]



# ensembl3 <- ensembl2[order(ensembl2$score, ensembl2$t, decreasing = TRUE), ]
# which(ensembl3$Symbol=="AR")
# which(ensembl3$Symbol=="GATA1")
# 
# 
# ensembl4 <- ensembl2[order(ensembl2$druggable, ensembl2$t, decreasing = TRUE), ]
# which(ensembl4$Symbol=="AR")



library(dplyr)
library(networkD3)
library(tidyr)


# group_start_idx <- which(!duplicated(ensembl2$score)) 
# ggdat <- ensembl2[c(group_start_idx[1]:(group_start_idx[1] + 4),
#                      group_start_idx[2]:(group_start_idx[2] + 4),
#                      group_start_idx[3]:(group_start_idx[3] + 4),
#                      group_start_idx[4]:(group_start_idx[4] + 4)),]


ggdat <- ensembl2[1:50,]
ggdat <- ggdat[,c("Symbol", "isDE", "pamr", "CCA", "isTF", "druggable")]
ggdat <- reshape2::melt(ggdat)
ggdat <- ggdat[ggdat$value ==1,]
ggdat <- ggdat[,c("Symbol","variable")]
colnames(ggdat) <- c("source","target")

# ggdat$source <- droplevels(ggdat$source)
# ggdat$source <- factor(ggdat$source, levels = score_dat$source)


n <- length(unique(ggdat$source)) + length(unique(ggdat$target))

# need to re-order nodes by biological relevance/importance here
nodeID <- c(0:(n-1))[as.factor(c(ggdat$source, ggdat$target))]

names(nodeID) <- c(as.character(ggdat$source),
                   as.character(ggdat$target))





# replace the arbitrary values in "value" by the actual statistic from the method
pamr <- read.csv("shrunkenNN_classification.csv")

library(PMA)
load("peerm.out.RData")
print(peerm.out)

CCA <- xlsx::read.xlsx("CCA_basal_luminals.xlsx", sheetIndex = 3)


# all measurements should be in range 0-1
# ----------------------------------------

#  <<< unit length normalisation >>>
norm_pamr_tic <- abs(pamr$TIC.score)/sqrt(sum(pamr$TIC.score^2))
# norm_pamr_tic <- abs(pamr$nonTIC.score)/sqrt(sum(pamr$nonTIC.score^2))
wCCA <- 1/CCA$AveRankBD
wCCA <- abs(wCCA)/sqrt(sum(wCCA^2))

pval <- ensembl2$adj.P.Val[match(ggdat$source[ggdat$target %in% "isDE"], ensembl2$Symbol)]
pval <-  abs(pval)/sqrt(sum(pval^2))

ggdat$value <- NA
ggdat$value[ggdat$target %in% "isDE"] <- pval
ggdat$value[ggdat$target %in% "pamr"] <- norm_pamr_tic[match(ggdat$source[ggdat$target %in% "pamr"], pamr$Symbol)]
ggdat$value[ggdat$target %in% "CCA"] <- wCCA[match(ggdat$source[ggdat$target %in% "CCA"], CCA$Symbol)]
ggdat$value[ggdat$target %in% "isTF"] <- 1


# store raw scores -----
norm_pamr_tic <- pamr$TIC.score
wCCA <- 1/CCA$AveRankBD
pval <- ensembl2$adj.P.Val[match(ggdat$source[ggdat$target %in% "isDE"], ensembl2$Symbol)]

gg_raw <- ggdat
gg_raw$value <- NA
gg_raw$value[gg_raw$target %in% "isDE"] <- pval
gg_raw$value[gg_raw$target %in% "pamr"] <- norm_pamr_tic[match(gg_raw$source[gg_raw$target %in% "pamr"], pamr$Symbol)]
gg_raw$value[gg_raw$target %in% "CCA"] <- wCCA[match(gg_raw$source[gg_raw$target %in% "CCA"], CCA$Symbol)]
# -------------------------

count_of_drugs <- sapply(ggdat$source[ggdat$target %in% "druggable"], FUN=function(x) {
  length(unique(drug_targets$parent_key[drug_targets$Symbol %in% x]))
  } )

count_of_drugs <- count_of_drugs/sqrt(sum(count_of_drugs^2))


#ggdat$value[ggdat$target %in% "druggable"] <- count_of_drugs
ggdat$value[ggdat$target %in% "druggable"] <- 0.1



summary(ggdat$value)


# need to obtain a summary statistics of the values in each category to ensure the values are comparable. 
ggdat$value[ggdat$target %in% "isTF"] <- 0.15


## prepare data including raw and normalised scores to output per method scores ------
score_dat <- aggregate(value ~ source, FUN= sum, data = ggdat)
score_dat <- score_dat[order(score_dat$value, decreasing = TRUE),]


overall_scores <- score_dat
colnames(overall_scores) <- gsub("value","overall_L2norm_score", colnames(overall_scores))

normalised_scores_per_method <- tidyr::spread(ggdat, "target", "value")
normalised_scores_per_method <- normalised_scores_per_method[match(overall_scores$source, 
                                                                   normalised_scores_per_method$source),]



colnames(normalised_scores_per_method)[2:4] <- paste0("L2norm_", colnames(normalised_scores_per_method)[2:4]) 





raw_scores_per_method <- tidyr::spread(gg_raw, "target", "value")
raw_scores_per_method <- raw_scores_per_method[match(overall_scores$source, 
                                                     raw_scores_per_method$source),]



colnames(raw_scores_per_method)[2:4] <- paste0("raw_", colnames(raw_scores_per_method)[2:4]) 

raw_scores_per_method <- raw_scores_per_method[, grep("drug|TF", 
                                                      colnames(raw_scores_per_method), 
                                                      invert=TRUE)]



scores_all <- cbind(raw_scores_per_method,
                    normalised_scores_per_method[,-1],
                    overall_scores[,-1, drop=FALSE])


ggdat2 <- matrix(c(0:(length(unique(c(ggdat$source, ggdat$target)))))[as.factor(c(ggdat$source, ggdat$target))], byrow = FALSE, ncol = 2)
colnames(ggdat2) <- c("source", "target")
ggdat2 <-  data.frame(ggdat2)

# ggdat2$value <- rep(10,204)
# #ggdat2$value[grep("druggable", ggdat$target)] <- 40
# ggdat2$value[grepl("ZEB1|^AR|LEF1|RFX8|TBX18", ggdat$source)] <- 80

ggdat2$value <- ggdat$value









library(ggsci)


jcolors::jcolors("pal7")
jcolors::jcolors("pal4")
jcolors::jcolors("pal6")
jcolors::jcolors("pal5")
RColorBrewer::brewer.pal(8, "Set2")



links <- ggdat2
nodes <- data.frame(name= as.character(unique(names(nodeID)[order(nodeID)])))


drugable_targets <- unique(ggdat$source[ggdat$target == "druggable"])


links$group=ifelse(grepl("^AR", ggdat$source),"AR","other")
links$group[grepl("^ZEB1", ggdat$source)] = "ZEB1"
nodes$group= ifelse(grepl("^AR", nodes$name), "AR","node_color")
nodes$group[grepl("^ZEB1", nodes$name)] = "ZEB1"



nodes$group[match(c("pamr", "CCA","isDE","isTF","druggable"), nodes$name)] <- paste0("method_node",1:5)
nodes$group <- as.factor(nodes$group)
nodes$name <- gsub("pamr","Shrunken Nearest Centroid", nodes$name)
# nodes$name <- gsub("SNC","Shrunken Nearest Centroid", nodes$name)
nodes$name <- gsub("CCA","Canonical Correlation Analysis", nodes$name)
nodes$name <- gsub("isDE","Differentially Expressed", nodes$name)
nodes$name <- gsub("isTF","Transcription Factor", nodes$name)
nodes$name <- gsub("druggable","Druggable", nodes$name)

my_color <- 'd3.scaleOrdinal() .domain(["ZEB1","AR", "other", "node_color", "method_node1",
"method_node2", "method_node3", "method_node4", "method_node5"]) .range(["#628395", "#D33B44", "#B3B3B3", "grey" ,"#7F8E39",
"#42858C","#570D32", "#E5C616", "#1d3554", "#DFE07C"])'


isTF <- unique(links$target[ggdat$target %in% "isTF"])
druggable <- unique(links$target[ggdat$target %in% "druggable"])


links$target[ggdat$target == "isTF"] <- links$source[ggdat$target == "isTF"]
links$source[ggdat$target == "isTF"] <- rep(isTF, sum(ggdat$target == "isTF"))

links$target[ggdat$target == "druggable"] <- links$source[ggdat$target == "druggable"]
links$source[ggdat$target == "druggable"] <- rep(druggable, sum(ggdat$target == "druggable"))

# have to re-arrange nodes based on importance/relevance score



networkD3::sankeyNetwork(Links = links,
                         Nodes = nodes,
                         Source = 'target',
                         Target = 'source',
                         Value = 'value',
                        NodeID = 'name',
                        fontSize = 16,
                        nodeWidth = 30,
                        nodePadding = 5,
                        #colourScale = my_color,
                        #NodeGroup = 'group',
                        #LinkGroup = 'group',
                        # width = 500,
                        fontFamily = "Helvetica",
                        #fontFamily = "Calibri",
                        
                        # original size dimensions -----
                        # width = 1100,
                        # height = 1300
                        
                        width = 600,
                        height = 900,
                        
                        iterations = 0
                        )



# Another idea is to set the distance to be proportional to the rank of the gene in our ranking scheme.


# test the color scheme ------
# links <- data.frame(
#   source=c("group_A","group_A", "group_B", "group_C", "group_C", "group_E"), 
#   target=c("group_C","group_D", "group_E", "group_F", "group_G", "group_H"), 
#   value=c(2,3, 2, 3, 1, 3)
# )
# 
# # From these flows we need to create a node data frame: it lists every entities involved in the flow
# nodes <- data.frame(
#   name=c(as.character(links$source), as.character(links$target)) %>% 
#     unique()
# )
# 
# # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
# links$IDsource <- match(links$source, nodes$name)-1 
# links$IDtarget <- match(links$target, nodes$name)-1
# 
# # prepare color scale: I give one specific color for each node.
# my_color <- 'd3.scaleOrdinal() .domain(["group_A", "group_B","group_C", "group_D", "group_E", "group_F", "group_G", "group_H"]) .range(["blue", "blue" , "blue", "red", "red", "yellow", "purple", "purple"])'
# 
# # Make the Network. I call my colour scale with the colourScale argument
# p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
#                    Value = "value", NodeID = "name", colourScale=my_color)
# p
# 
# 
# nodes$group <- as.factor(c("a","a","a","a","a","b","b","b"))
# 
# # Give a color for each group:
# my_color <- 'd3.scaleOrdinal() .domain(["a", "b"]) .range(["#69b3a2", "steelblue"])'
# 
# # Make the Network
# p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
#                    Value = "value", NodeID = "name", 
#                    colourScale=my_color, NodeGroup="group")
# p
# 
# 
# links$group <- as.factor(c("type_a","type_a","type_a","type_b","type_b","type_b"))
# 
# # Add a 'group' column to each node. Here I decide to put all of them in the same group to make them grey
# nodes$group <- as.factor(c("my_unique_group"))
# 
# # Give a color for each group:
# my_color <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "my_unique_group"]) .range(["#69b3a2", "steelblue", "grey"])'
# 
# # Make the Network
# p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
#                    Value = "value", NodeID = "name", 
#                    colourScale=my_color, LinkGroup="group", NodeGroup="group")
# 
# p

