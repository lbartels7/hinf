import sgkit
from sgkit.io.vcf import vcf_to_zarr
        
vcf_to_zarr(snakemake.input[0], snakemake.output[0])