# Liver zonation and plasticity in MASLD

Hackathon project, team work. It uses single cell data of liver cells (Gribben et al., GSE202379) from healthy livers and from every stage of fatty liver disease (MASLD).

For every hepatocyte it computes:

* a score for the cell's zone (periportal or pericentral) from marker genes
* the entropy of that score, which shows how much the cell lost its zone identity
* a plasticity score, which is high when the cell expresses both hepatocyte and cholangiocyte genes

Then it averages them per donor and per disease stage, checks which one changes first across the stages, runs the statistical tests and saves the plots.

## How to run

1. export the marker gene counts from the Seurat object once in R (the R code is at the top of `zonation_plasticity.py`)
2. put `gribben_metadata.csv` and `marker_counts.csv` next to the script
3. run `python zonation_plasticity.py`
