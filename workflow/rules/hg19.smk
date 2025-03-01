"""
To generate a reference file from hg19 assembly and download corresponding annotation files
"""

CHRS = [str(i) for i in range(1, 23)] + ["X", "Y"]
HG19_BASE = "hg19"


rule hg19_chr_ref:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(HG19_BASE, "{chr}.ref")),
    params:
        fa="Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/",
    script:
        "src/get_chr_ref.sh"


rule hg19_ref:
    input:
        expand(op.join(HG19_BASE, "{chr}.ref"), chr=CHRS),
    params:
        base=HG19_BASE,
    output:
        op.join(HG19_BASE, "ref.gz"),
    script:
        "src/make_ref.sh"


rule get_hg19_pmds:
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
