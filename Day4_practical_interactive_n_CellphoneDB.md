---
title: "Practical_day4_interactive_n_CellphoneDB"
output: html_document
---

# In this practical you will:
- work with interactive tools to get familiar with representations of single-cell RNA-seq data
- touch on biological interpretation via gene set enrichment analysis and cell-cell interaction analysis

## Exercise 1: Publicly available single-cell portals - Exploring covid19cellatlas
Authors can create a browseable portal for your single cell analysis. For example, a cellxgene format can be used: https://chanzuckerberg.github.io/cellxgene/

We will explore one such dataset in a portal that combines a variety of datasets.
1. Go to https://www.covid19cellatlas.org/
2. Choose "Healthy Donors"
3. Observe the available datasets.
4. Select **Deprez et al "Total airway"** cell portal for the Respiratory System and open **Interactive viewer**.

### Task 1: Explore the covid19 atlas dataset
- Color by cell types, Donor, Location, Method and other metadata. Explain the biases and comment on the batch effect(s).

### Task 2: Search for genes
5. Search for the genes **ACE2** (the receptor for cellular entry of the Coronavirus) and **TMPRSS2** (a protease needed for entry). 
6. Select double-positive cells.
 - Explore which cell types the two genes are found in. What does it tell you? Why are there so few cells? Try again with using only ACE2 for selection.

### Task 3: Differential expression 
7. Perform differential expression analysis between Multiciliated (trachea) and Multiciliated N (nasal) cell types. 
 - What does this tell you?

## Exercise 2: Publicly available single-cell portals - Exploring Single Cell Expression Atlas by EMBL-EBI.
This is an automated pipeline outcome where analysis is standardised on all submitted data.

8. Go to https://www.ebi.ac.uk/gxa/sc/home
9. Search for **ACE2** in the *Gene ID or gene gene symbol* search bar.
10. Choose **Organism Part: Lung** at the left drop-down box.
11. Choose the 1st *Homo sapiens* entry: **Ischaemic sensitivity of human tissue by single cell RNA seq**

### Task 4: Explore the dataset. What can you learn about the ACE2 gene? 
 - Try out different metadata by scrolling up the selection in "Color plot by:" dropdown menu. 
 - Comment on the batch effect and integration
 - Can you find where to download the raw counts matrix? 
 - Where can you find other single cell experiments?


## Conclusion
We looked at two different ways of sharing your data. One is tailored approach by the author in choosing cells, thresholds, batch correction etc manually, while the other uses a common approach for all datasets. Based on these two samples, can you comment on the benefits/shortcomings of both approaches?


## Exercise 3: Exploring the practical dataset further
In yesterday's practical, you looked at the immune cells in covid infected vs healthy lung bronchioalveolar lavage. Let's explore the files that you saved.

### Task 5. Look for marker genes in Human Protein Atlas (HPA) for tissue localisation https://www.proteinatlas.org/ 
12. Open the file with highly expressed genes in healthy macrophages in LibreOffice:
/home/training/practical_day3/healthy_macs.tsv

These are the genes upregulated in healthy macrophages, thus they should be possible to find in the HPA. 

13. Search for the first gene, FABP4, in the HPA Tissue atlas. 
14. Choose Lung tissue, subsection Lung, and open the images. 
 - Can you see where the macrophages are? 
 - Why are the images so different? 
15. Another regular macrophage marker is in the top list - **C1QA**. Check this in the HPA in a similar way. 
 - Do you notice anything in particular in the gene lists?


### Task 6. Explore gene set enrichment in the covid-affected macrophages.
16. Open the file with highly expressed genes in covid macrophages in LibreOffice:
/home/training/practical_day3/covid_macs.tsv
17. Select the top 200 DE genes, and perform gene set enrichment analysis using the gProfileR online tool: https://biit.cs.ut.ee/gprofiler/gost
18. Observe significant enrichments. 
 - Do they make sense? Repeat the same with top 20 genes. How consistent are the results?


## Exercise 4: Perform receptor-ligand analysis via the CellPhoneDB tool https://www.cellphonedb.org/
The CellPhoneDB tool explores receptor and ligand expression in all the cell type pairs in your dataset. This let's you infer both cell signaling and biologicical function.
For smaller datasets, the interactive tool allows for analysis online: https://www.cellphonedb.org/explore-sc-rna-seq
The interactive tool is often busy, therefore we use a command line version in the practical. 

Package download and instructions for use (This has been installed for you):
https://github.com/Teichlab/cellphonedb

Installing was done for you before starting (**DO NOT RUN THIS INSTALLATION CODE**) <br>
```python
Python3 -m venv cpdb-venv
``` 
Installation was done for you in advance, but the rest is up to you!

### Task 7: Create and move to a directory for cellphoneDB
19. Open Terminal (by clicking on the black squared icon on the left side of screen).

20. For starters, let's all get to the same directory
```shell
cd /home/training/
```

21. Activate the environment of the CellPhoneDB
```shell
source cpdb-venv/bin/activate
```
22. Create a directory for receptor-ligand analysis
```shell
mkdir cellphoneDB
```
23. Move to the newly created directory
```shell
cd cellphoneDB
```

### Task 8: Perform cellphonedb calculations
You previously saved counts from both COVID positive and COVID negative cells. We now ask the CellPhoneDB algorithm to perform analysis of receptor-ligand expression on paired cell types.

24. Perform calculation of means and p-value on the COVID positive patients
```shell
cellphonedb method statistical_analysis ../practical_day3/cellphonedb_meta_disease.txt ../practical_day3/cellphonedb_count_disease.txt --counts-data=ensembl --iterations=200 --threshold=0.5 --project-name=COVID_pos --threads=16
```
25. Perform calculation of means and p-value on the COVID negative patients
```shell
cellphonedb method statistical_analysis ../practical_day3/cellphonedb_meta_healthy.txt ../practical_day3/cellphonedb_count_healthy.txt --counts-data=ensembl --iterations=200 --threshold=0.5 --project-name=COVID_neg --threads=16
```

26. Move to the COVID positive directory and plot the number of significant interactions between cell types
```shell
cd /home/training/cellphoneDB/out/COVID_pos/
cellphonedb plot heatmap_plot ../../../practical_day3/cellphonedb_meta_disease.txt --pvalues-path=pvalues.txt
```

File: /home/training/cellphoneDB/out/COVID_neg/out/heatmap_count.pdf

27. Move to the COVID negative directory and plot the number of significant interactions between cell types:
```shell
cd /home/training/cellphoneDB/out/COVID_neg/
cellphonedb plot heatmap_plot ../../../practical_day3/cellphonedb_meta_healthy.txt --pvalues-path=pvalues.txt
```

File: /home/training/cellphoneDB/out/COVID_neg/out/heatmap_count.pdf

28. Open the plots and look at them. 
 - Observe differences between covid and healthy. 
 - What do you observe?


### Task 9: Plot the significant interactions as a dot plot in R
29. Open R(studio) and input the following code
30. Read in cellphoneDB analysis output for COVID positive patients

```R
pval=read.delim('cellphoneDB/out/COVID_pos/pvalues.txt')
dim(pval)
colnames(pval)
head(pval[,1:10])
``` 

Every row is a receptor-ligand interaction pair. Columns starting from the 11 are cell type - cell type pairs. You'll now choose all the receptor-ligand pairs with significant p-value. 

31. Start by finding the minimum p-value across cell type - cell type pairs
```R
pval$minPval=apply(pval[,11:ncol(pval)], 1, min) 
```
32. Observe distribution of the minimum p values
```R
table(pval$minPval) 
```
33. Choose the interactions that have significant p-value in at least one cell type - cell type pair.
```R
rows=pval[pval$minPval!=1, 2]
rows
```
34. Save this list of receptor-ligand pairs for use in plotting
```R
write.table(rows, "/home/training/cellphoneDB/out/rows.tsv", row.names = F, col.names = T, sep="\t", quote=F) 
```

###  Task 10: Plot the dot plots with selected receptor-ligand combinations in the Terminal
35. Type the following into your terminal
```bash
cd /home/training/cellphoneDB/out/COVID_pos
cellphonedb plot dot_plot --means-path=means.txt --pvalues-path=pvalues.txt --rows ../rows.tsv --output-name plot.pdf
```
This creates a File: /home/training/cellphoneDB/out/COVID_pos/out/plot.pdf
```bash
cd /home/training/cellphoneDB/out/COVID_neg
cellphonedb plot dot_plot --means-path=means.txt --pvalues-path=pvalues.txt --rows ../rows.tsv --output-name plot.pdf
```
This creates a File: /home/training/cellphoneDB/out/COVID_neg/out/plot.pdf

36. Open the plots and comment on the differences between covid and healthy samples. 
 - What do you observe?
 - For which purpose can this method be used? What are the benefits and drawbacks of cellPhoneDB?

## Final words
Thank you and good-bye!

