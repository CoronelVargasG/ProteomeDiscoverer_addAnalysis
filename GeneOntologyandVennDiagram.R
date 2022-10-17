#1) export an excel file from PD using your proteins of interest, and import it in your environment
library(readxl)
EXPERIMENT1 <- read_excel("PD_myproteinsfile_EXPERIMENT1.xlsx")
EXPERIMENT2 <- read_excel("PD_myproteinsfileEXPERIMENT2.xlsx")
#2) IF NECESSARY Create a venn diagram to detect the common proteins between different experiments

library(VennDiagram)
venn.diagram(x= list (EXPERIMENT1$Accession, EXPERIMENT1$Accession),
             fill = c("orange", "blue"),
             category.names = c("Experiment1" , "Experiment2"),
             filename = "Venn_diagram_1.png",
             output = TRUE)

#create the graphic
library(ggvenn)
library(RColorBrewer)
x <- list(Experiment1=EXPERIMENT1$Accession, Experiment2=EXPERIMENT2$Accession)

mydiagram <- ggvenn(x, 
                    show_elements = F,
                    show_percentage = FALSE,
                    label_sep = "\n",
                    fill_color = brewer.pal(name="Set2",n=3),
                    set_name_size = 4,
                    text_size = 3)

#Extract all common UniprotIDs
all_common <- intersect(EXPERIMENT1$Accession, 
                        EXPERIMENT2$Accession)

#Convert the common proteins UniprotIDs to GeneSymbol
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

cols <- c("SYMBOL", "GENENAME", "ENTREZID", "UNIPROT")
my_genes <- select(org.Hs.eg.db, keys=all_common, columns=cols, keytype="UNIPROT")

write.table(my_genes, file = "common_deregulated_genes.txt", sep = "\n")
my_genes$ENTREZID
my_genes$SYMBOL

######################GO-BP###########
library(miRBaseConverter)
library(miRNAtap)
library(miRNAtap.db)
library(topGO)
library(org.Hs.eg.db)
library(GOplot)

allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL,
                         mapping="org.Hs.eg.db", ID = "entrez")
#Prepare Your data
ForRank_EXPERIMENT1<- na.omit(data.frame(Accession=EXPERIMENT1$`Gene Symbol`,
                                   AdjPvalue=EXPERIMENT1$`Abundance Ratio Adj. P-Value: (TREATMENT) / (CONTROL)`))

cols2 <- c("SYMBOL", "GENENAME", "ENTREZID")
ForRank_EXPERIMENT1_GeneID <- select(org.Hs.eg.db, keys=ForRank_EXPERIMENT1$Accession,
                                     columns=cols2, keytype="SYMBOL")
ForRank_EXPERIMENT1_GeneID_adjPvalue <- data.frame(ForRank_EXPERIMENT1_GeneID$ENTREZID,
                                             ForRank_EXPERIMENT1$AdjPvalue)

ForRankFinal_EXPERIMENT1 <- tibble::deframe(ForRank_EXPERIMENT1_GeneID_adjPvalue)
selection_EXPERIMENT1	=	function(x) TRUE
#Do the GO analysis
#EXPERIMENT1
GOdata_EXPERIMENT1	=	new('topGOdata', ontology = 'BP',
                         allGenes = ForRankFinal_EXPERIMENT1,
                         annot = annFUN.GO2genes, GO2genes = allGO2genes,
                         geneSel = selection_EXPERIMENT1	, nodeSize=20)

results.ks_EXPERIMENT1	=	runTest(GOdata_EXPERIMENT1, algorithm = "classic", statistic = "ks")

allRes_EXPERIMENT1	=	GenTable(GOdata_EXPERIMENT1, KS = results.ks_EXPERIMENT1,
                              orderBy = "KS", topNodes = 20)

allRes_data_EXPERIMENT1	=	allRes_EXPERIMENT1	[,c('GO.ID','Term','KS')]

write.csv(allRes_data_EXPERIMENT1, file='allRes_data_EXPERIMENT1.cvs')

