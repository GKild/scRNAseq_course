---
title: "CellPhoneDB"
output: html_document
---


CellPhoneDB allows for receptor-ligand analysis between the annotated cell populations. 

The interactive tool is sometimse busy, therefore we use a command line version in the practical. For smaller datasets, the interactive tool allows for analysis online: https://www.cellphonedb.org/explore-sc-rna-seq

Package download and instructions for use:
https://github.com/Teichlab/cellphonedb



```{r eval=FALSE}
## The collowing commands executed in the terminal:

#Activate the environment of the cellPhoneDB
source cpdb-venv/bin/activate
cd cellphoneDB

#Perform calculation of means and p-value on the COVID positive cells:
cellphonedb method statistical_analysis covid_pos_meta.csv covid_pos_counts.csv --counts-data=hgnc_symbol --iterations=200 --threshold=0.05 --project-name=COVID_pos  --threads=16

#Perform calculation of means and p-value on the COVID positive cells:
cellphonedb method statistical_analysis covid_neg_meta.csv covid_neg_counts.csv --counts-data=hgnc_symbol --iterations=200 --threshold=0.05 --project-name=COVID_neg  --threads=16


#Plot the results for COVID positive data:
cd out/COVID_pos/
cellphonedb plot dot_plot --means-path=means.txt --pvalues-path=pvalues.txt
cellphonedb plot heatmap_plot ../../covid_pos_meta.csv --pvalues-path=pvalues.txt

#Plot the results for COVID negative data:
cd ../COVID_neg/
cellphonedb plot dot_plot --means-path=means.txt --pvalues-path=pvalues.txt
cellphonedb plot heatmap_plot ../../covid_neg_meta.csv --pvalues-path=pvalues.txt


```


