###########################################################################
# Genome annotation pipeline using braker in Snakemake
# Snakemake/5.13.0
###########################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards
from Bio import SeqIO

result_dir = config["result_dir"]
input_dir = config["input_dir"]
rna_dir = config["rna_dir"]
script_dir = config["script_dir"]

protein_file = config["protein_file"]
rna_list = config["rna_list"]

SAMPLE = list(glob_wildcards(join(input_dir, "{ids}.fasta")))[0]

print(SAMPLE)

rule All:
    input:
        # Repeats
        expand(join(result_dir,"{samples}-families.fa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.softMasked.fasta"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.fasta.masked"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.fasta.out.gff"),samples=SAMPLE),

        # braker files
        expand(join(result_dir,"{samples}_braker/braker.gff3"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker/braker.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker/braker.aa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker/braker.codingseq"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker/braker.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker/braker.clean.gtf"),samples=SAMPLE),

        # renamed file
        expand(join(result_dir,"{samples}_braker.gff3"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker.aa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker.cds"),samples=SAMPLE),

        # functional file
        expand(join(result_dir,"{samples}_braker.functional.gff3"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker.functional.aa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}_braker.functional.cds"),samples=SAMPLE),

rule braker:
    input:
        fa=join(result_dir, "{samples}.softMasked.fasta"),
        prot=protein_file,
    output:
        gff=join(result_dir, "{samples}_braker/braker.gff3"),
        aa=join(result_dir, "{samples}_braker/braker.aa"),
        cds=join(result_dir, "{samples}_braker/braker.codingseq"),
    params:
        rname="braker",
        species_id="{samples}",
        out_dir=join(result_dir,"{samples}_braker"),
        rna_dir=rna_dir,
        rna_list=rna_list,
    shell:
        """
        module load braker
        mkdir -p {params.out_dir}
        braker.pl --genome={input.fa} --useexisting --species={params.species_id} \
        --prot_seq={input.prot} --workingdir={params.out_dir} \
        --gff3 --threads=8 --rnaseq_sets_ids={params.rna_list}  \
        --rnaseq_sets_dir={params.rna_dir}
        """

rule gff_rename:
    input:
        gff=join(result_dir, "{samples}_braker/braker.gff3"),
        aa=join(result_dir, "{samples}_braker/braker.aa"),
        cds=join(result_dir, "{samples}_braker/braker.codingseq"),
    output:
        map=temp(join(result_dir, "{samples}.map")),
        gff=join(result_dir, "{samples}_braker.gff3"),
        aa=join(result_dir, "{samples}_braker.aa"),
        cds=join(result_dir, "{samples}_braker.cds"),
    params:
        species_id="{samples}",
        rname="gff_rename",
    shell:
        """
        module load maker
        maker_map_ids --prefix {params.species_id} --justify 5  {input.gff} > {output.map}
        cp {input.gff} {output.gff}
        map_gff_ids {output.map} {output.gff}
        cp {input.aa} {output.aa}
        cp {input.cds} {output.cds}
        map_fasta_ids {output.map} {output.aa}
        map_fasta_ids {output.map} {output.cds}
        """

rule uniprot:
    input:
        fa=expand(join(input_dir, "{samples}.fasta"),samples=SAMPLE),
    output:
        uniprot=temp("uniprot_sprot.fasta"),
    params:
        rname="uniprot",
    shell:
        """
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip uniprot_sprot.fasta.gz
        """

rule gff_annot:
    input:
        gff=join(result_dir, "{samples}_braker.gff3"),
        prot=join(result_dir, "{samples}_braker.aa"),
        cds=join(result_dir, "{samples}_braker.cds"),
        uniprot="uniprot_sprot.fasta",
    output:
        gff=join(result_dir, "{samples}_braker.functional.gff3"),
        prot=join(result_dir, "{samples}_braker.functional.aa"),
        cds=join(result_dir, "{samples}_braker.functional.cds"),
        blast=temp(join(result_dir,"{samples}.blast")),
    params:
        rname="gff_annot",
        threads=8,
    shell:
        """
        module load blast maker
        makeblastdb -dbtype prot -in {input.uniprot}
        blastp -query {input.prot} -db {input.uniprot} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output.blast} -num_threads {params.threads}
        maker_functional_gff {input.uniprot} {output.blast} {input.gff} > {output.gff}
        maker_functional_fasta {input.uniprot} {output.blast} {input.prot} > {output.prot}
        maker_functional_fasta {input.uniprot} {output.blast} {input.cds} > {output.cds}
        """

rule RepeatModeler:
  input:
    fa=join(input_dir, "{samples}.fasta"),
  output:
    fa=join(result_dir, "{samples}.fasta"),
    rep=join(result_dir,"{samples}-families.fa"),
  params:
    rname="RepeatModeler",
    dir=result_dir,
    id="{samples}",
  threads:
    48
  shell:
    """
    cd {params.dir}
    module load repeatmodeler
    ln -s {input.fa} {output.fa}
    BuildDatabase -name {params.id} {output.fa}
    RepeatModeler -database {params.id} -pa {threads} -LTRStruct >& {params.id}.out
    """

rule RepeatMasker:
  input:
    fa=join(result_dir, "{samples}.fasta"),
    rep=join(result_dir,"{samples}-families.fa"),
  output:
    fa=join(result_dir,"{samples}.fasta.masked"),
    gff=join(result_dir,"{samples}.fasta.out.gff"),
  params:
    rname="RepeatMasker",
    dir=result_dir,
    threads="48",
  shell:
    """
    cd {params.dir}
    module load repeatmasker
    RepeatMasker -u -s -poly -engine rmblast -pa {params.threads} -gff -no_is -gccalc -norna -lib {input.rep} {input.fa}
    """

rule softMask:
    input:
        fa=join(result_dir, "{samples}.fasta"),
        gff=join(result_dir,"{samples}.fasta.out.gff"),
    output:
        fa=join(result_dir, "{samples}.softMasked.fasta"),
    params:
        rname="softMask",
    shell:
        """
        module load bedtools
        bedtools maskfasta -fullHeader -soft -fi {input.fa} -bed {input.gff} -fo {output.fa}
        """

rule gff2gtf:
    input:
        gff=join(result_dir,"{samples}_braker/braker.gff3"),
    output:
        gtf=join(result_dir,"{samples}_braker/braker.gtf"),
        clean=join(result_dir,"{samples}_braker/braker.clean.gtf"),
    params:
        rname="gff2gtf",
        script_dir=script_dir,
    shell:
        """
        module load agat python
        agat_convert_sp_gff2gtf.pl --gff {input.gff} -o {output.gtf}
        {params.script_dir}/clean_gtf.py {output.gtf} > {output.clean}
        """