"""
To generate a reference file from hg19 assembly and download corresponding annotation files
"""

CHRS = [str(i) for i in range(1, 23)] + ["X", "Y"]
HG19_BASE = "hg19"


rule hg19_chr_ref:
    conda:
        op.join("..", "envs", "yamet.yml")
    output:
        temp(op.join(HG19_BASE, "{chr}.ref")),
    params:
        fa="Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/",
    script:
        "src/get_chr_ref.sh"


rule hg19_ref:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        expand(op.join(HG19_BASE, "{chr}.ref"), chr=CHRS),
    params:
        base=HG19_BASE,
    output:
        op.join(HG19_BASE, "ref.gz"),
    script:
        "src/make_ref.sh"
