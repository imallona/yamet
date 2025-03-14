"""
To generate a reference file from mm10 assembly and download corresponding annotation files
"""

CHRS = [str(i) for i in range(1, 20)] + ["X", "Y"]
MM10_BASE = "mm10"


rule mm10_chr_ref:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(MM10_BASE, "{chr}.{meth_pat}.ref")),
    params:
        fa="Mus_musculus.GRCm38.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/",
    script:
        "src/get_chr_ref.sh"


rule mm10_ref:
    input:
        expand(op.join(MM10_BASE, "{chr}.{{meth_pat}}.ref"), chr=CHRS),
    params:
        base=MM10_BASE,
    output:
        op.join(MM10_BASE, "ref.{meth_pat}.gz"),
    script:
        "src/make_ref.sh"


rule get_mm10_promoters:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(MM10_BASE, "promoters.bed.gz"),
    script:
        "src/download_mm10_promoters.sh"


rule get_mm10_genes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(MM10_BASE, "genes.bed.gz"),
    script:
        "src/download_mm10_genes.sh"


ENCODE_MAP = {
    "h3k4me3": "ENCFF160SCR",
    "h3k9me3": "ENCFF658QTP",
    "h3k4me1": "ENCFF937JHP",
    "h3k27me3": "ENCFF827BBC",
    "h3k27ac": "ENCFF442GIT",
}


rule get_mm10_encode:
    output:
        op.join(MM10_BASE, "{ann}.bed.gz"),
    params:
        accessor=lambda wildcards: f"{ENCODE_MAP[wildcards.ann]}",
    shell:
        """
            curl -L https://www.encodeproject.org/files/{params.accessor}/@@download/{params.accessor}.bed.gz -o {output}
        """


rule get_mm10_genome_sizes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(MM10_BASE, "genome.sizes"),
    params:
        asm="mm10",
    script:
        "src/download_genome_sizes.sh"
