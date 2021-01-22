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


Task: Search for the genes ACE2 (the receptor for cellular entry of the Coronavirus) and TMPRSS2 (a protease needed for entry). Select double-positive cells and explore which cell types they are found in. What does it tell you? Why are there so few cells? Try again with only ACE2 used for selection.


Task3: Perform differential expression between Multiciliated (trachea) and Multiciliated N (nasal) cell types.



2. Exploring Single Cell Expression Atlas by EMBL-EBI.
This is an automated pipeline outcome where analysis is standardised on all submitted data.
Go to https://www.ebi.ac.uk/gxa/sc/home

Task 1. Choose Organism part "Lung" from the dropdown menu on the left, choose the 1st entry "Ischaemic sensitivity of human tissue by single cell RNA seq"
What can you learn about the ACE2 gene? Try out different metadata by scrolling up the selection in "Color plot by:" dropdown menu. Comment on the batch effect and integration. Can you find where to download the raw counts matrix? Where can you find other single cell experiments?


Conclusion: We looked at two different ways of sharing your data. One is tailored approach by the author in choosing cells, thresholds, batch correctin etc manually, the other is common approach for all. Based on these two samples, can you comment on the benefits/shortcomings of both approaches?


####  Break ####


3. Exploring own dataset further
In yesterday's practical you looked at the immune cells in covid infected vs healthy lung bonchiolveolar lavage. Let's explore the files that you saved.

Task 1. Look for marker genes in Human Protein Atlas (HPA) for tissue localisation https://www.proteinatlas.org/ 
Open the file with highly expressed genes in healthy macrophages in LibreOffice:
/home/training/practical_day3/healthy_macs.tsv

These are the genes upregulated in healthy macrophages, thus should be possible to find in the HPA. Search the first gene, FABP4 in HPA Tissue atlas. Choose Lung tissue, subsection Lung, open images. Can you see where the macrophages are? Why are the images so different? Another regular macrophage marker is in the top list - C1QA - check this in the HPA in a similar way. Do you notice anything in particular in the gene lists?


Task 2. Explore gene set enrichment in the covid-affected macrophages.
Open the file with highly expressed genes in covid macrophages in LibreOffice:
/home/training/practical_day3/covid_macs.tsv
Select top 200 DE genes, and perform gene set enrichment analysis using gProfileR online tool: https://biit.cs.ut.ee/gprofiler/gost
Observe significant enrichments. Do they make sense? Repeat the same with top 20 genes. How consistent are the results?


4. Perform receptor-ligand analysis via the CellPhoneDB tool https://www.cellphonedb.org/
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


## Final words
Thank you and good-bye!

