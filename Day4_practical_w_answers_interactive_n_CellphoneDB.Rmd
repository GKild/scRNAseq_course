---
title: "Practical_day4_interactive_n_CellphoneDB"
output: html_document
---

# In this practical you will:
- work interactive tool to get familiar with representations of single-cell RNA-seq data
- touch on biological interpretation via gene set enrichemnt analysis and cell-cell interaction analysis


1. Exploring publicly available single-cell portal.
It is possible for authors to create a browseable portal for your single cell analysis. For example, a cellxgene format can be used: https://chanzuckerberg.github.io/cellxgene/

We browse one such dataset in a portal that combines a variety of datasets.
Go to https://www.covid19cellatlas.org/, choose "Healthy Donors"
Observe available datasets.
Select Deprez et al "Total airway" cell portal for the Respiratory System and open Interactive viewer.

Task1: Explore the datset: color by cell types, Donor, Location, Method and other metadata. Explain the biases, comment on the batch effect(s).

Answer: Donors nicely integrated, location and method biases are expected due to biological differences.

Task: Search for the genes ACE2 (the receptor for cellular entry of the Coronavirus) and TMPRSS2 (a protease needed for entry). Select double-positive cells and explore which cell types they are found in. What does it tell you? Why are there so few cells? Try again with only ACE2 used for selection.

Answer: 235 cells double positive. Mainly in multiciliated, suprabasal and secretory cells, higher in N (nasal) than trachea (regular). Indicates higher contagiousness as closer to the outside. So few cells as double-positive due to the drop-out effect.

Task3: Perform differential expression between Multiciliated (trachea) and Multiciliated N (nasal) cell types.

Answer: Select 1st cell type by clicking only Multiciliated cells in the CellTypes metadata field. Select 2nd cell type by clicking on Multiciliated N. Show 10 most Differentially Expressed genes by clicking on the button with Venn Diagram.


2. Exploring Single Cell Expression Atlas by EMBL-EBI.
This is an automated pipeline outcome where analysis is standardised on all submitted data.
Go to https://www.ebi.ac.uk/gxa/sc/home

Task 1. Choose Organism part "Lung" from the dropdown menu on the left, choose the 1st entry "Ischaemic sensitivity of human tissue by single cell RNA seq"
What can you learn about the ACE2 gene? Try out different metadata by scrolling up the selection in "Color plot by:" dropdown menu. Comment on the batch effect and integration. Can you find where to download the raw counts matrix? Where can you find other single cell experiments?

Answer: There are clusters with detected ACE2. The cells cluster by metadata features (including batches), indicating that the clusters might represent batch effects as well as the cell types. For raw counts and other data, click the "Downloads" button. For more datasets click "Browse experiments" on the top.

Conclusion: We looked at two different ways of sharing your data. One is tailored approach by the author in choosing cells, thresholds, batch correctin etc manually, the other is common approach for all. Based on these two samples, can you comment on the benefits/shortcomings of both approaches?

3. Annotating cell types by Azimuth tool
In the Day 3 practical you used a code to predict cell type labels. For immune cells, there is an online reference tool Azimuth that uses it's own CITE-seq  human peripheral blood mononuclear cells reference for mapping queries: https://satijalab.org/azimuth/ (https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1). 

Task 1. Use the Azimuth tool to explore the covid dataset annotations.
Open Azimuth tool (https://satijalab.org/azimuth/) -> Visit the App
Upload the dataset used in practical (/home/training/practical_day3/reference.rds)
Wait for the data to be mapped, look at the results. Choose to show a table with prediced annotations and the author-provided annotations (annot). How many annotations match/don't match? Do the nong-matching annotations make sense? Is the use of Azimuth tool justified in this case?

Answer: Some mismatch is expected, as single cells often are not so good quality map every single one correctly. Also, some discrepancy between difficult-to-distinguish populations is very normal. Different levels of annotations will map higher category in one dataset between multiple subpopulations in more finely annotated dataset. In lack of exact population, the next bes should be chosen, which explains how macrophages mapped as monocytes. Mismatches therefore make sense. In current case the covid data is from lungs and the reference is from blood. The quite good match of populations increases confidence in our annotations, but should not be used for annotating this data as the only option.

####  Break ####


4. Exploring own dataset further
In yesterday's practical you looked at the immune cells in covid infected vs healthy lung bonchiolveolar lavage. Let's explore the files that you saved.

Task 1. Look for marker genes in Human Protein Atlas (HPA) for tissue localisation https://www.proteinatlas.org/ 
Open the file with highly expressed genes in healthy macrophages in LibreOffice:
/home/training/practical_day3/healthy_macs.tsv # Hint: if there was problem with creating the file, they are also found in /home/training/practical_day3/ELO

These are the genes upregulated in healthy macrophages, thus should be possible to find in the HPA. Search the first gene, FABP4 in HPA Tissue atlas. Choose Lung tissue, subsection Lung, open images. Can you see where the macrophages are? Why are the images so different? Another regular macrophage marker is in the top list - C1QA - check this in the HPA in a similar way. Do you notice anything in particular in the gene lists?

Answer: Macrophages are on the surface of the very thin epithelium where the oxygen exchange happens. This is only found on the best quaity images. The images can vary depending on the location of origin, donor, antibody and other technical & biological factors. The gene list in general seems to have a lot of ribosomal genes in the top. This is often considered technical.

Task 2. Explore gene set enrichment in the covid-affected macrophages.
Open the file with highly expressed genes in covid macrophages in LibreOffice:
/home/training/practical_day3/covid_macs.tsv
Select top 200 DE genes, and perform gene set enrichment analysis using gProfileR online tool: https://biit.cs.ut.ee/gprofiler/gost
Observe significant enrichments. Do they make sense? Repeat the same with top 20 genes. How consistent are the results?

Answer: There are many enrichments that are immune (response) related. There is also COVID-19 related enrichment that makes most sense. The categories are roughly similar with top 20 genes, but with lower significance.

5. Perform receptor-ligand analysis via the CellPhoneDB tool https://www.cellphonedb.org/
CellPhoneDB tool explores receptor and ligand expression in all the cell type pairs in your dataset. This allowes to infer cell signaling and infer biologicical function.
For smaller datasets, the interactive tool allows for analysis online: https://www.cellphonedb.org/explore-sc-rna-seq
The interactive tool is often busy, therefore we use a command line version in the practical. 

Package download and instructions for use (This has been installed for you):
https://github.com/Teichlab/cellphonedb

Installing was done for you before starting:
Python3 -m venv cpdb-venv


## The collowing commands are executed in the Terminal.
# Open Terminal (by clicking on the black squared icon on the left side of screen).

# For starters, let's all get to the same directory:
cd /home/training/
# Activate the environment of the CellPhoneDB
source cpdb-venv/bin/activate
# Create a directory for receptor-ligand analysis:
mkdir cellphoneDB

# Move to the newly created directory:
cd cellphoneDB

## You previously saved counts from both COVID positive and COVID negative cells. Followingly, we ask the CellPhoneDB algorithm to perform analysis of receptor-ligand expression on paired cell type type of manner.

# Perform calculation of means and p-value on the COVID positive patients:
cellphonedb method statistical_analysis ../practical_day3/cellphonedb_meta_disease.txt ../practical_day3/cellphonedb_count_disease.txt --counts-data=ensembl --iterations=200 --threshold=0.5 --project-name=COVID_pos --threads=16

# Perform calculation of means and p-value on the COVID negative patients:
cellphonedb method statistical_analysis ../practical_day3/cellphonedb_meta_healthy.txt ../practical_day3/cellphonedb_count_healthy.txt --counts-data=ensembl --iterations=200 --threshold=0.5 --project-name=COVID_neg --threads=16


# Move to the COVID positive directory and plot the number of significant interactions between cell types:
cd /home/training/cellphoneDB/out/COVID_pos/
cellphonedb plot heatmap_plot ../../../practical_day3/cellphonedb_meta_disease.txt --pvalues-path=pvalues.txt
#File: /home/training/cellphoneDB/out/COVID_neg/out/heatmap_count.pdf

# Move to the COVID negative directory and plot the number of significant interactions between cell types:
cd /home/training/cellphoneDB/out/COVID_neg/
cellphonedb plot heatmap_plot ../../../practical_day3/cellphonedb_meta_healthy.txt --pvalues-path=pvalues.txt
#File: /home/training/cellphoneDB/out/COVID_neg/out/heatmap_count.pdf

#Open the plots and look at them. Observe differences between covid and healthy. What do you observe?

#Answer: A lot more interactions in the covid data. The largest amount of interactions in the macrophages could be confounded by the largest number of cells in.
#Ready-made plots are also found in /home/training/CPDB/out/COVID_neg/out/ and /home/training/CPDB/out/COVID_pos/out/


# The next section is for plotting the significant interactions as a dot plot.
# This is done in R. Open R(studio) and follow the commands: 

```{r}
#Read in cellphoneDB analysis output for COVID positive patients:
pval=read.delim('cellphoneDB/out/COVID_pos/pvalues.txt')
dim(pval)
colnames(pval)
head(pval[,1:10])

#Every row is receptor-ligand ineraction pair, columns starting from the 11 are cell type - cell type pairs.
#Choose all the receptor-ligand pairs with significant p-value.

pval$minPval=apply(pval[,11:ncol(pval)], 1, min) # Find minimum p-value across cell type - cell type pairs
table(pval$minPval) #Observe distribution of the minimum p values
rows=pval[pval$minPval!=1, 2] #Choose the interactions that have significant p-value in at least one cell type - cell type pair.
rows
write.table(rows, "/home/training/cellphoneDB/out/rows.tsv", row.names = F, col.names = T, sep="\t", quote=F) #Save this list of receptor-ligand pairs for use in plotting
```

## Plot the dot plots with selected receptor-ligand combinations in the Terminal:
cd /home/training/cellphoneDB/out/COVID_pos
cellphonedb plot dot_plot --means-path=means.txt --pvalues-path=pvalues.txt --rows ../rows.tsv --output-name plot.pdf
#File: /home/training/cellphoneDB/out/COVID_pos/out/plot.pdf

cd /home/training/cellphoneDB/out/COVID_neg
cellphonedb plot dot_plot --means-path=means.txt --pvalues-path=pvalues.txt --rows ../rows.tsv --output-name plot.pdf
#File: /home/training/cellphoneDB/out/COVID_neg/out/plot.pdf

#Open the plots and comment on the differences between covid and healthy samples. What do you observe?

For which purpose can this method be used? What are the benefits and drawbacks of cellPhoneDB?

#Answer: More significant interactions in the covid samples, which was also evident from the previous heatmap plots. The method can be used to look at any cell-cell interactions , although the immune interactions are a big part of the data. This can really help interpretation of the single cell data (to infer function on the new cell types for example). For more comfortable usage, this tool could be integrated with the Seurat workflow.
#Ready-made plots are also found in /home/training/CPDB/out/COVID_neg/out/plot2.pdf and /home/training/CPDB/out/COVID_pos/out/plot2.pdf


## Final words
Thank you and good-bye


