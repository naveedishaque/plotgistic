# plotgistic

Scripts to create nice gistic visualisations.

Please excuse the bad coding, and light documentation!

## Prerequisites

- Developed using perl v5.26.1 built for x86_64-linux-gnu-thread-multi
- Developed using R v3.6.1 built for x86_64-conda_cos6-linux-gnu
- AFAIK no special libraries are required so this should work on any perl v5.x.x and R v3.x.x 

## Usage

First run gistic2.0 (this is tested on the output of gistic 2.0.23)

### Preprocess data

First create a table that will be used for the plotting script:

```
perl gistic_2_dataframe_for_plotting.pl scores.gistic del_genes.conf_20.txt amp_genes.conf_20.txt > myProject_gistic_table.tsv
```

### Plot

Create the plot:

```
R -f plot_gistic_table.r --args myProject_gistic_table.tsv
display myProject_gistic_table.tsv.pdf
```

### Plot with custom genes

Create the plot, but adding in other genes that you are interested in (in a 4 column BED file, no chr prefix, and 4th column is the name)

```
R -f plot_gistic_table.r --args myProject_gistic_table.tsv my_favourite_genes.bed
```

If you want to change the genes reported by gistic, find the lines with "Gene_amp" and "Gene_del", and modify the column with the gene names

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

Naveed Ishaque

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Dorett I Odoni
