# Revisiting mutational resistance to ampicillin and cefotaxime in *Haemophilus influenzae*

This repository contains the snakemake workflow used to perform the
genome-wide association study described in the publication. Once the
pipeline has been successfully executed, the figures can be recreated
using the Jupyter notebook supplied.

Diricks, M., Petersen, S., Bartels, L. et al. Revisiting mutational resistance to ampicillin and cefotaxime in Haemophilus influenzae. Genome Med 16, 140 (2024).
 https://doi.org/10.1186/s13073-024-01406-4


# Installation
Please install [miniconda](https://docs.anaconda.com/miniconda/) first. Afterwards run
```
conda create -n hinf --file environment.yml
```
to install all the packages required to run the pipeline.

# Running the pipeline
The pipeline was implemented using the workflow manager [snakemake](https://snakemake.readthedocs.io/en/stable/).

Run the following command to execute the pipeline:
```
snakemake --cores <number of cores> --use-conda
```
After the execution is finished, all results can be found in the subfolder named results.