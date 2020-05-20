# plotgistic

Scripts to create nice gistic visualisations.

Please excuse the bad coding, and light documentation!

##Plotting the output (@naveedishaque)

First run gistic2.0 (this is tested in gistic 2.0.23)

###Preprocess data

First create a table that will be used for the plotting script:

```
module load perl
perl gistic_2_dataframe_for_plotting.pl scores.gistic del_genes.conf_20.txt amp_genes.conf_20.txt > myProject_gistic_table.tsv
```

###Plot

Create the plot:

```
module load R
R -f plot_gistic_table.r --args myProject_gistic_table.tsv
display myProject_gistic_table.tsv.pdf
```

###Plot with custom genes

Create the plot, but adding in other genes that you are interested in (in a 4 column BED file, no chr prefix, and 4th column is the name)

```
module load R
R -f plot_gistic_table.r --args myProject_gistic_table.tsv my_favourite_genes.bed
```

If you want to change the genes reported by gistic, find the lines with "Gene_amp" and "Gene_del", and modify the column with the gene names
