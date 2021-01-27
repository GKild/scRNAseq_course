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

*Task 1 Answer: Donors are nicely integrated, location and method biases are expected due to biological differences.*

### Task 2: Search for genes
5. Search for the genes **ACE2** (the receptor for cellular entry of the Coronavirus) and **TMPRSS2** (a protease needed for entry). 
6. Select double-positive cells.
 - Explore which cell types the two genes are found in. What does it tell you? Why are there so few cells? Try again with using only ACE2 for selection.

*Task 2 Answer: 235 cells are double positive. Mainly in multiciliated, suprabasal and secretory cells, higher in N (nasal) than trachea (regular). Indicates higher contagiousness as closer to the outside. So few cells as double-positive due to the drop-out effect.*

### Task 3: Differential expression 
7. Perform differential expression analysis between Multiciliated (trachea) and Multiciliated N (nasal) cell types. 
 - What does this tell you?

*Task 3 Answer: Select 1st cell type by clicking only Multiciliated cells in the CellTypes metadata field. Select 2nd cell type by clicking on Multiciliated N. Show 10 most Differentially Expressed genes by clicking on the button with Venn Diagram.*

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

*Task 4 Answer: There are clusters with detected ACE2. The cells cluster by metadata features (including batches), indicating that the clusters might represent batch effects as well as the cell types. For raw counts and other data, click the "Downloads" button. For more datasets click "Browse experiments" on the top.*

## Conclusion
We looked at two different ways of sharing your data. One is tailored approach by the author in choosing cells, thresholds, batch correction etc manually, while the other uses a common approach for all datasets. Based on these two samples, can you comment on the benefits/shortcomings of both approaches?

## Exercise 3. Annotating cell types by Azimuth tool
In the Day 3 practical, you used a code to predict cell type labels. For immune cells, there is an online reference tool Azimuth that uses it's own CITE-seq  human peripheral blood mononuclear cells reference for mapping queries: https://satijalab.org/azimuth/ (https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1). 

### Task 5: Use the Azimuth tool to explore the covid dataset annotations.
12. Open the Azimuth tool (https://satijalab.org/azimuth/) -> Visit the App
13. Upload the dataset used in the practical (/home/training/practical_day3/reference.rds). Wait for the data to be mapped and then look at the results. 
14. Show a table with prediced annotations and the author-provided annotations (annot). 
 - How many annotations match/don't match? 
 - Do the non-matching annotations make sense? 
 - Is the use of Azimuth tool justified in this case?

*Task 5 Answer: Some mismatch is expected, as single cells often are not so good quality that every single cell maps correctly. Also, some discrepancy between difficult-to-distinguish populations is very normal. Different levels of annotations will map to a higher category in one dataset and between multiple subpopulations in more finely annotated dataset. In lieu of an exact population, the next best should be chosen, which explains how macrophages mapped as monocytes. Mismatches therefore make sense. In the current case, the covid data is from lungs and the reference is from blood. The quite good match of populations increases the confidence in our annotations, but should not be used for annotating this data as the only option.*

## Exercise 4: Exploring the practical dataset further
In yesterday's practical, you looked at the immune cells in covid infected vs healthy lung bronchioalveolar lavage. Let's explore the files that you saved.

### Task 6. Look for marker genes in Human Protein Atlas (HPA) for tissue localisation https://www.proteinatlas.org/ 
15. Open the file with highly expressed genes in healthy macrophages in LibreOffice:
/home/training/practical_day3/healthy_macs.tsv

*Hint: if there was problem with creating the file, they are also found in /home/training/practical_day3/ELO*

These are the genes upregulated in healthy macrophages, thus they should be possible to find in the HPA. 

16. Search for the first gene, FABP4, in the HPA Tissue atlas. 
17. Choose Lung tissue, subsection Lung, and open the images. 
 - Can you see where the macrophages are? 
 - Why are the images so different? 
18. Another regular macrophage marker is in the top list - **C1QA**. Check this in the HPA in a similar way. 
 - Do you notice anything in particular in the gene lists?

*Task 6 Answer: Macrophages are on the surface of the very thin epithelium where the oxygen exchange happens. This is only found on the best quaity images. The images can vary depending on the location of origin, donor, antibody and other technical & biological factors. The gene list in general seems to have a lot of ribosomal genes in the top. This is often considered technical.*

### Task 7. Explore gene set enrichment in the covid-affected macrophages.
19. Open the file with highly expressed genes in covid macrophages in LibreOffice:
/home/training/practical_day3/covid_macs.tsv
20. Select the top 200 DE genes, and perform gene set enrichment analysis using the gProfileR online tool: https://biit.cs.ut.ee/gprofiler/gost
21. Observe significant enrichments. 
 - Do they make sense? Repeat the same with top 20 genes. How consistent are the results?

*Task 7 Answer: There are many enrichments that are immune (response) related. There is also COVID-19 related enrichment that makes most sense. The categories are roughly similar with top 20 genes, but with lower significance.*

## Exercise 5: Perform receptor-ligand analysis via the CellPhoneDB tool https://www.cellphonedb.org/
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

### Task 8: Create and move to a directory for cellphoneDB
22. Open Terminal (by clicking on the black squared icon on the left side of screen).

23. For starters, let's all get to the same directory
```shell
cd /home/training/
```

24. Activate the environment of the CellPhoneDB
```shell
source cpdb-venv/bin/activate
```
25. Create a directory for receptor-ligand analysis
```shell
mkdir cellphoneDB
```
26. Move to the newly created directory
```shell
cd cellphoneDB
```

### Task 9: Perform cellphonedb calculations
You previously saved counts from both COVID positive and COVID negative cells. We now ask the CellPhoneDB algorithm to perform analysis of receptor-ligand expression on paired cell types.

27. Perform calculation of means and p-value on the COVID positive patients
```shell
cellphonedb method statistical_analysis ../practical_day3/cellphonedb_meta_disease.txt ../practical_day3/cellphonedb_count_disease.txt --counts-data=ensembl --iterations=200 --threshold=0.5 --project-name=COVID_pos --threads=16
```
28. Perform calculation of means and p-value on the COVID negative patients
```shell
cellphonedb method statistical_analysis ../practical_day3/cellphonedb_meta_healthy.txt ../practical_day3/cellphonedb_count_healthy.txt --counts-data=ensembl --iterations=200 --threshold=0.5 --project-name=COVID_neg --threads=16
```

29. Move to the COVID positive directory and plot the number of significant interactions between cell types
```shell
cd /home/training/cellphoneDB/out/COVID_pos/
cellphonedb plot heatmap_plot ../../../practical_day3/cellphonedb_meta_disease.txt --pvalues-path=pvalues.txt
```

File: /home/training/cellphoneDB/out/COVID_neg/out/heatmap_count.pdf

30. Move to the COVID negative directory and plot the number of significant interactions between cell types:
```shell
cd /home/training/cellphoneDB/out/COVID_neg/
cellphonedb plot heatmap_plot ../../../practical_day3/cellphonedb_meta_healthy.txt --pvalues-path=pvalues.txt
```

File: /home/training/cellphoneDB/out/COVID_neg/out/heatmap_count.pdf

31. Open the plots and look at them. 
 - Observe differences between covid and healthy. 
 - What do you observe?

*Task 9 Answer: A lot more interactions in the covid data. The largest amount of interactions in the macrophages could be confounded by the largest number of cells in.*

*Hint: Ready-made plots are also found in /home/training/CPDB/out/COVID_neg/out/ and /home/training/CPDB/out/COVID_pos/out/*

### Task 10: Plot the significant interactions as a dot plot in R
32. Open R(studio) and input the following code
33. Read in cellphoneDB analysis output for COVID positive patients

```R
pval=read.delim('cellphoneDB/out/COVID_pos/pvalues.txt')
dim(pval)
colnames(pval)
head(pval[,1:10])
``` 

Every row is a receptor-ligand interaction pair. Columns starting from the 11 are cell type - cell type pairs. You'll now choose all the receptor-ligand pairs with significant p-value. 

34. Start by finding the minimum p-value across cell type - cell type pairs
```R
pval$minPval=apply(pval[,11:ncol(pval)], 1, min) 
```
35. Observe distribution of the minimum p values
```R
table(pval$minPval) 
```
36. Choose the interactions that have significant p-value in at least one cell type - cell type pair.
```R
rows=pval[pval$minPval!=1, 2]
rows
```
37. Save this list of receptor-ligand pairs for use in plotting
```R
write.table(rows, "/home/training/cellphoneDB/out/rows.tsv", row.names = F, col.names = T, sep="\t", quote=F) 
```

###  Task 11: Plot the dot plots with selected receptor-ligand combinations in the Terminal
38. Type the following into your terminal
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

39. Open the plots and comment on the differences between covid and healthy samples. 
 - What do you observe?
 - For which purpose can this method be used? What are the benefits and drawbacks of cellPhoneDB?
 
*Task 11 Answer: More significant interactions in the covid samples, which was also evident from the previous heatmap plots. The method can be used to look at any cell-cell interactions , although the immune interactions are a big part of the data. This can really help interpretation of the single cell data (to infer function on the new cell types for example). For more comfortable usage, this tool could be integrated with the Seurat workflow.*

*Hint: Ready-made plots are also found in /home/training/CPDB/out/COVID_neg/out/plot2.pdf and /home/training/CPDB/out/COVID_pos/out/plot2.pdf*

## Final words
Thank you and good-bye


