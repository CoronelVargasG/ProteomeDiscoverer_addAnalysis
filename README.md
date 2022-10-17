# ProteomeDiscoverer_addAnalysis
Venn diagram and GO-BP analysis using proteins from ProteomeDiscoverer differential analysis

This file is the pipeline to do a GO-BP analysis using an excel exported file from Thermofisher Proteome discoverer.
Before exporting remember to check your up and down regulated proteins from the differential analysis.


1) Create a Venn diagram if you want to find the common proteins between two or more Experiments (Example: you want to 
see if there are common deregulated proteins in different tissues from a mice after treatment with a drug)

2) Convert from "UNIPROT" Accession to "SYMBOL", and export the list of common proteins

GO-BP analysis
1) create the object allGO2genes <- annFUN.org(whichOnto='BP', feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez")
2) Convert your accessions to gene "SYMBOL"
3) Create a new deta.frame using only the columnsENTEREZID just converted and adding the column AdjPvalue from "ForRank_EXPERIMENT1" element
4) create a named vector element using the tibble::deframe function
5) Run the GO analysis
6) Export a file with your GOs and the KS values (best limit <0.05). You can use Revigo tools (http://revigo.irb.hr/) to create immages as treemaps, scatterplots, 
and Cytoscape-compatible (https://cytoscape.org/) plots.
