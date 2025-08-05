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


rule hg19_genes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "genes.genes.bed.gz"),
    script:
        "src/download_hg19_genes.sh"


rule hg19_get_pmds:
    output:
        temp(op.join(HG19_BASE, "pmd.bed.gz")),
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
        op.join(HG19_BASE, "{md}.pmd.bed.gz"),
    params:
        filter=lambda wildcards: PMD_MAP[wildcards.md],
    shell:
        """
            gunzip -c {input} |
                grep "{params.filter}$" |
                gzip -c >{output}
        """


rule hg19_get_hmm:
    output:
        temp(op.join(HG19_BASE, "hmm.bed")),
    params:
        loc="https://www.encodeproject.org/files/ENCFF526MRN/@@download/ENCFF526MRN.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n >{output[0]}
        """


rule hg19_hmm:
    input:
        op.join(HG19_BASE, "hmm.bed"),
    output:
        op.join(HG19_BASE, "{ann}.hmm.bed.gz"),
    shell:
        """
            grep "{wildcards.ann}" {input} | gzip -c >{output}
        """


CHIP_MAP = {
    "H3K27me3": "ENCFF255ARD",
    "H3K9me3": "ENCFF354CNQ",
    "H3K4me3": "ENCFF893TVK",
}


rule hg19_chip:
    output:
        op.join(HG19_BASE, "{chip}.chip.bed.gz"),
    params:
        loc=lambda wildcards: f"https://www.encodeproject.org/files/{CHIP_MAP[wildcards.chip]}/@@download/{CHIP_MAP[wildcards.chip]}.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c >{output[0]}
        """


rule hg19_lamin:
    output:
        op.join(HG19_BASE, "laminb1.lad.bed.gz"),
    params:
        loc="https://github.com/jernst98/ChromHMM/raw/refs/heads/master/COORDS/hg19/laminB1lads.hg19.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output[0]}
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
        op.join(HG19_BASE, "bookended.custom.bed.gz"),
    shell:
        """
            bedtools makewindows -g {input[0]} -w 10000 |
                sort -k1,1 -k2,2n | gzip -c >{output[0]}
        """


rule hg19_windows_250k:
    input:
        op.join(HG19_BASE, "genome.sizes"),
    output:
        op.join(HG19_BASE, "bookended_250k.custom.bed.gz"),
    shell:
        """
            bedtools makewindows -g {input[0]} -w 250000 |
                sort -k1,1 -k2,2n | gzip -c >{output[0]}
        """


rule hg19_single_annotation_coverage:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        windows=op.join(HG19_BASE, "bookended.custom.bed.gz"),
        annotation=op.join(HG19_BASE, "{subcat}.{cat}.bed.gz"),
    output:
        annotated_windows=temp(op.join(HG19_BASE, "{subcat}.{cat}.annotation.frac")),
    shell:
        """
            echo "{wildcards.subcat}_{wildcards.cat}" >{output.annotated_windows}
            bedtools coverage -a {input.windows} \
                -b {input.annotation} | cut -f7 >>{output.annotated_windows}
        """


ANN_MAP = {
    "pmd": ["pmds", "hmds"],
    "hmm": [
        "0_Enhancer",
        "2_Enhancer",
        "11_Promoter",
        "12_Promoter",
        "1_Transcribed",
        "4_Transcribed",
        "5_RegPermissive",
        "7_RegPermissive",
        "6_LowConfidence",
        "3_Quiescent",
        "8_Quiescent",
        "10_Quiescent",
        "9_ConstitutiveHet",
        "13_ConstitutiveHet",
    ],
    "chip": ["H3K27me3", "H3K9me3", "H3K4me3"],
    "lad": ["laminb1"],
    "custom": ["bookended"],
}


def list_annotated_windows():
    res = []
    for cat in ANN_MAP:
        for subcat in ANN_MAP[cat]:
            res.append(f"{subcat}.{cat}.annotation.frac")
    return [op.join("hg19", item) for item in res]


rule combine_annotated_windows:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        annotated_windows=list_annotated_windows(),
    output:
        op.join(HG19_BASE, "bookend_annotation.gz"),
    shell:
        """
            paste {input.annotated_windows} | gzip -c >{output}
        """
