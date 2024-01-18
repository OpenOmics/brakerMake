# brakeMake
A repository for the Braker Snakemake pipeline, designed to be run on the completed assemblies, after assembly & scaffolding.

![](pipeline_rulegraph_new.svg)

### Before running the pipeline
Edit the config.yaml as follows:

 * Update 'input_dir' to be the path where all assemblies you want to annotate are stored. Make sure that all assemblies end in '.fasta'.
 * Update 'results_dir' to the path where all results will be stored.
 * Change 'rna_dir' to the location where bulk RNA sequencing files are stored to be used for transcriptome generation.
 * Change 'rna_list' to the list of IDs for files within the rna_dir. Can also include SRA IDs which are not in the rna_dir, these will be downloaded and used for the annotation.

### Test running the pipeline

The pipeline can be test run with the following command in an interactive session:

```sh pipeline_ctrl.sh npr $PWD```

Assuming that you are in the directory from this repository.

### Running the pipeline

The pipeline can be fully run with the following command:

```sbatch --time=10-00:00:00 pipeline_ctrl.sh process $PWD```

Assuming that you are in the directory from this repository.
