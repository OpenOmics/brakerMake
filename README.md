# annotateMake
A repository for the Braker Snakemake pipeline, designed to be run on the completed assemblies, after assembly & scaffolding.


### Before running the pipeline
Edit the config.yaml as follows:

 * Update 'input_dir' to be the path where all assemblies you want to annotate are stored. Make sure that all assemblies end in '.fasta'.
 * Update 'results_dir' to the path where all results will be stored.
 * Change 'lineage_name' to the Busco odb10 database which will be used in this analysis, chosen from https://busco-archive.ezlab.org/data/lineages/.
 * Change 'augustus_name' to the most appropriate, taken from the list from running ```funannotate species```.
 * Change 'transcript_file' and 'protein_file' to the path of known transcripts/proteins for this species, or a closely related species.
 * Change 'species' to be the species name of the assembly being annotated, or a unique identifier.
 * If repeat sequences for this genome are known, change 'repeat_file' to the path with the fasta of these repeat sequences. If repeat sequences are not known, uncomment the repeatmodeler rule/output & change the input repeats for repeatmasker/fun_mask to be the repeatmodeler output.
 * Change 'funannotate_dir' to be the path to this downloaded github directory. This is where busco odb10 databases will be stored and where the parameter files & custom scripts are stored & accessed.

### Test running the pipeline

The pipeline can be test run with the following command in an interactive session:

```sh pipeline_ctrl.sh npr $PWD```

Assuming that you are in the directory from this repository.

### Running the pipeline

The pipeline can be fully run with the following command:

```sbatch --time=10-00:00:00 pipeline_ctrl.sh process $PWD```

Assuming that you are in the directory from this repository.
