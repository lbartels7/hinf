# Revisiting mutational resistance to ampicillin and cefotaxime in *Haemophilus influenzae*
Placeholder for small explanations and citation right after publication.

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