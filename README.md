# brakeMake
A repository for the Braker Snakemake pipeline, designed to be run on the completed assemblies, after assembly & scaffolding.

![](pipeline_rulegraph_new.svg)

### Before running the pipeline
Edit the config.yaml as follows:

 * Update 'input_dir' to be the path where all assemblies you want to annotate are stored. Make sure that all assemblies end in '.fasta'.
 * Update 'results_dir' to the path where all results will be stored.
 * Change 'rna_dir' to the location where bulk RNA sequencing files are stored to be used for transcriptome generation. Even when no RNAseq is used this needs to be an existing directory.
 * Change 'rna_list' to the list of IDs for files within the rna_dir. Can also include SRA IDs which are not in the rna_dir, these will be downloaded and used for the annotation. If no RNAseq is used, leave this section empty.
 * If you do not wish to use the uniprot protein DB, 'protein_file' needs to be updated to the path to a fasta file of amino acid sequences to use. This fasta file then needs to be used to generate a blastp database, this protein DB will be used within the script. If you wish to use other protein fasta files alongside uniprot, they should be put in the input_dir with the assemblies to be annotated. Make sure that all protein datasets end in '.faa'.

### Test running the pipeline

The pipeline can be test run with the following command in an interactive session:

```sh pipeline_ctrl.sh npr $PWD```

Assuming that you are in the directory from this repository.

### Running the pipeline

The pipeline can be fully run with the following command:

```sbatch --time=10-00:00:00 pipeline_ctrl.sh process $PWD```

Assuming that you are in the directory from this repository.
