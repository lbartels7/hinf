import sgkit as sg
from sgkit.io.vcf import vcf_to_zarr

import numpy as np
import pandas as pd
import xarray as xr

input_zarr = snakemake.input[0] # type: ignore
metadata = snakemake.input[1] # type: ignore
amended_positions = pd.read_csv(snakemake.input[2], header=None).iloc[:,0]
minimum_allele_count = snakemake.params["minimum_allele_count"] # type: ignore

output_zarr = snakemake.output[0] # type: ignore

def remove_heterozgous_variants(ds : xr.Dataset) -> xr.Dataset:
    """Remove all heterozgous variants from a diploid dataset


    Args:
        ds (xr.Dataset): Input dataset from which the heterozgous variants are removed

    Returns:
        xr.Dataset: Output dataset without any heterozgous variants
    """
    
    return(
        ds
        .pipe(lambda ds: ds.sel(variants=(ds.call_genotype[:,:,0] == ds.call_genotype[:,:,1]).all(dim='samples').compute()))
        )

def filter_minimum_AC(ds: xr.Dataset, minimum_AC: int = 20) -> xr.Dataset:
    """Remove all variants with a minimum allele count below a certain threshold from a dataset.

    Args:
        ds (xr.Dataset): Input dataset
        minimum_AC (int, optional): Minimum allele count for a variant to remain in the data set. Defaults to 20.

    Returns:
        xr.Dataset: Dataset where all variant have a minimum allele count above the given threshold.
    """

    variables_to_drop = ['variant_n_called', 'variant_call_rate',
                          'variant_n_het', 'variant_n_hom_ref', 'variant_n_hom_alt', 'variant_n_non_ref',
                          'variant_allele_count', 'variant_allele_total', 'variant_allele_frequency']

    return(
        ds
        .pipe(sg.variant_stats)
        .pipe(lambda ds: ds.sel(variants=(ds.variant_allele_count[:,:2].min(dim='alleles') >= minimum_AC).compute()))
        .drop(variables_to_drop)
        )

def filter_maximum_hetero(ds: xr.Dataset, maximum_hetero: int = 10_000) -> xr.Dataset:
    """The maximum number of heterozgous sites of a sample to keep it in the dataset

    Args:
        ds (xr.Dataset): Input dataset
        maximum_hetero (int, optional): Threshold value above which all samples are removed if they contain more heterozygous sites. Defaults to 10_000.

    Returns:
        xr.Dataset: Output data set in which there are no samples with more than the threshold value for heterozygous variants.
    """

    variables_to_drop = ['sample_n_called', 'sample_call_rate', 'sample_n_het', 'sample_n_hom_ref', 'sample_n_hom_alt',
                          'sample_n_non_ref']

    return(
        ds
        .pipe(sg.sample_stats)
        .pipe(lambda ds: ds.sel(samples=(ds.sample_n_het < maximum_hetero).compute()))
        .drop(variables_to_drop)
        )

def remove_indels(ds: xr.Dataset) -> xr.Dataset:
    """Remove all indels from the dataset. Drops each variant if either the reference or the alternative allele
    is longer than one base.

    Args:
        ds (xr.Dataset): Input dataset that may contain indels

    Returns:
        xr.Dataset: Output dataset wherer all indels are removed.
    """
    ref = ds['variant_allele'][:,0].str.len() == 1
    alt = ds['variant_allele'][:,1].str.len() == 1
    return(
        ds
        .pipe(lambda ds: ds.sel(variants=(ref & alt).compute()))
        )

def filter_position(ds: xr.Dataset, positions: pd.Series) -> xr.Dataset:
    """Remove all variants from the dataset, if their position in not present in positions

    Args:
        ds (xr.Dataset): Input dataset
        positions (pd.Series): Collection of positions that should be kept in the dataset

    Returns:
        xr.Dataset: Outputdataset with all positions removed if they are not in positions
    """
    isin_mask = ds['variant_position'].isin(positions)

    return(
        ds
        .pipe(lambda ds: ds.sel(variants=(isin_mask).compute()))
        )

# Drop all variables with the ploidy dimension + variant_id_mask and reduce call_genotype to haploid
ds = (sg.load_dataset(input_zarr)
      .pipe(remove_indels)
      .pipe(filter_position, amended_positions) 
      .drop_vars(['call_genotype_mask', 'call_genotype_phased', 'variant_id_mask'])
      .assign(call_genotype = lambda ds: ds.call_genotype.max(dim='ploidy')))

# Recreate the ploidy dimension in call_genotype, because downstream method calls expect it to be present
# Same for call_genotype_mask
ds['call_genotype'] = ds.call_genotype.expand_dims({'ploidy':1},2)
ds['call_genotype_mask'] = xr.zeros_like(ds.call_genotype).astype(bool)


# Load samples with beta_lactamase status negative and set relative mic values to their lower limit ('>8' -> 8)
df_metadata = (pd.read_csv(metadata, sep='\t')
               .assign(AMP_MIC= lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
)

# Convert pandas dataframe to xarray dataArray, so we can merge in the next step
ds_metadata = (df_metadata
               .to_xarray()
               .rename({"ID":"samples"})
               .swap_dims({'index': 'samples'})
               .drop(labels='index')
)
# Merge the genotype dataset with the phenotype dataArray (resistance status, MIC values)
ds = ds.set_index({"samples": "sample_id"})
ds = ds.merge(ds_metadata, join="left")
ds = ds.reset_index('samples').reset_coords('samples').rename_vars({'samples': 'sample_id'})

# Dropping variants that are too rare or too common
# Create an ID (ID_{position}_{alt-allele}) for later use
ds = (
    ds
    .pipe(filter_minimum_AC, minimum_allele_count)
    .assign(variant_id=lambda ds: ( 'ID_' + ds.variant_position.to_series().astype(str) + '_' + ds.variant_allele[:,1].to_series()))
)


# Save to zarr
sg.save_dataset(ds, output_zarr, auto_rechunk=True) # type: ignore




