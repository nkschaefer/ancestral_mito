# ancestral_mito
Create a mitochondrial genome for an imputed ancestor of two species, given an outgroup

### How does this work?
It aligns the mitochondrial genomes of the two species and an outgroup. For each aligned position, if the two species match, it chooses their base. If not, it chooses the outgroup base.

To account for the potential of mitochondrial genomes being linearized by breaking the circle in different places, all sequences have their last 1000 bp appended to the beginning and first 1000 bp appended to the end before alignment. The number of bases (1000 by default) is controlled by the `padding` parameter, which you can set in the parameter file.

## Install
```
conda env create --file=ancestral_mito.yml
```

## Prepare
Copy `example.yml` to your working directory and edit it to point to your files.

If you are running on a cluster, create a `nextflow.config` in your working directory and add [cluster-specific settings](https://www.nextflow.io/docs/latest/config.html).

## Run
```
conda activate ancestral_mito
nextflow [/path/to/this_repository/]ancestral_mito.nf -params-file [params].yml
```
