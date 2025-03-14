"""
To generate a reference file from hg19 assembly and download corresponding annotation files
"""

CHRS = [str(i) for i in range(1, 23)] + ["X", "Y"]
HG19_BASE = "hg19"


rule hg19_chr_ref:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(HG19_BASE, "{chr}.{meth_pat}.ref")),
    params:
        fa="Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/",
    script:
        "src/get_chr_ref.sh"


rule hg19_ref:
    input:
        expand(op.join(HG19_BASE, "{chr}.{{meth_pat}}.ref"), chr=CHRS),
    params:
        base=HG19_BASE,
    output:
        op.join(HG19_BASE, "ref.{meth_pat}.gz"),
    script:
        "src/make_ref.sh"


rule hg19_get_genes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "genes.bed.gz"),
    script:
        "src/download_hg19_genes.sh"


rule hg19_genes:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(HG19_BASE, "genes.bed.gz"),
    output:
        op.join(HG19_BASE, "genes.genes.bed"),
    shell:
        """
            gunzip -c {input} > {output}
        """


rule hg19_get_pmds:
    output:
        op.join(HG19_BASE, "pmd.bed.gz"),
    params:
        loc="https://zhouserver.research.chop.edu/GenomeAnnotation/hg19/PMD_coordinates_hg19.bed.gz",
    shell:
        """
        curl -o {output[0]} {params.loc}
        """


PMD_MAP = {"pmds": "commonPMD", "hmds": "commonHMD"}


rule hg19_pmds:
    input:
        op.join(HG19_BASE, "pmd.bed.gz"),
    output:
        temp(op.join(HG19_BASE, "{md}.pmd.bed")),
    params:
        filter=lambda wildcards: PMD_MAP[wildcards.md],
    shell:
        """
            gunzip -c {input} |
                grep "{params.filter}$" >{output}
        """


rule hg19_get_hmm:
    output:
        op.join(HG19_BASE, "hmm.bed.gz"),
    params:
        loc="https://www.encodeproject.org/files/ENCFF526MRN/@@download/ENCFF526MRN.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output[0]}
        """


rule hg19_hmm:
    input:
        op.join(HG19_BASE, "hmm.bed.gz"),
    output:
        temp(op.join(HG19_BASE, "{ann}.hmm.bed")),
    shell:
        """
            gunzip -c {input} | grep "{wildcards.ann}" > {output}
        """


CHIP_MAP = {
    "H3K27me3": "ENCFF255ARD",
    "H3K9me3": "ENCFF354CNQ",
    "H3K4me3": "ENCFF893TVK",
}


rule hg19_get_chip:
    output:
        op.join(HG19_BASE, "{chip}.bed.gz"),
    params:
        loc=lambda wildcards: f"https://www.encodeproject.org/files/{CHIP_MAP[wildcards.chip]}/@@download/{CHIP_MAP[wildcards.chip]}.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output[0]}
        """


rule hg19_chip:
    input:
        op.join(HG19_BASE, "{chip}.bed.gz"),
    output:
        temp(op.join(HG19_BASE, "{chip}.chip.bed")),
    shell:
        """
            gunzip -c {input} > {output}
        """


rule hg19_get_lamin:
    output:
        op.join(HG19_BASE, "laminb1.bed.gz"),
    params:
        loc="https://github.com/jernst98/ChromHMM/raw/refs/heads/master/COORDS/hg19/laminB1lads.hg19.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output[0]}
        """


rule hg19_lamin:
    input:
        op.join(HG19_BASE, "{lamin}.bed.gz"),
    output:
        temp(op.join(HG19_BASE, "{lamin}.lad.bed")),
    shell:
        """
            gunzip -c {input} > {output}
        """


rule get_hg19_genome_sizes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "genome.sizes"),
    params:
        asm="hg19",
    script:
        "src/download_genome_sizes.sh"


rule hg19_windows:
    input:
        op.join(HG19_BASE, "genome.sizes"),
    output:
        op.join(HG19_BASE, "bookended.custom.bed"),
    shell:
        """
            bedtools makewindows -g {input[0]} -w 10000 |
                sort -k1,1 -k2,2n >{output[0]}
        """
