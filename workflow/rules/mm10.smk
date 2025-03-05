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
