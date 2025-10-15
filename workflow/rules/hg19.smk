"""
To generate a reference file from hg19 assembly and download corresponding annotation files
"""

CHRS = [str(i) for i in range(1, 23)] + ["X", "Y"]
HG19_BASE = "hg19"

rule get_cpgislandext_hg9:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "cpgIslandExt.cpgIslandExt.bed")
    shell:
        """
        wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz \\
        | gunzip -c \\
        | awk 'BEGIN{{ OFS="\\t"; }}{{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }}' \\
        > {output}
        """
        
        
rule build_hg19_chr_per_chr:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(HG19_BASE, "{chr}.{meth_pat}.ref"))
    params:
        fa="Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/",
    script:
        "src/get_chr_ref.sh"

rule build_genome_hg19_ref:
    input:
        expand(op.join(HG19_BASE, "{chr}.{{meth_pat}}.ref"), chr=CHRS),
    params:
        base=HG19_BASE,
    output:
        op.join(HG19_BASE, "ref.{meth_pat}.gz"),
    script:
        "src/make_ref.sh"


rule get_genes_hg19:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "genes.bed.gz"),
    script:
        "src/download_hg19_genes.sh"


rule uncompress_hg19_genes:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(HG19_BASE, "genes.bed.gz")
    output:
        temp(op.join(HG19_BASE, "genes.genes.bed"))
    shell:
        """
            gunzip -c {input} > {output}
        """


rule get_pmds_hg19:
    output:
        op.join(HG19_BASE, "pmd.bed.gz"),
    params:
        loc="https://zhouserver.research.chop.edu/GenomeAnnotation/hg19/PMD_coordinates_hg19.bed.gz",
    shell:
        """
        curl -o {output[0]} {params.loc}
        """


PMD_MAP = {"pmds": "commonPMD", "hmds": "commonHMD"}

rule uncompress_pmds_hg19:
    input:
        op.join(HG19_BASE, "pmd.bed.gz")
    output:
        temp(op.join(HG19_BASE, "{md}.pmd.bed"))
    params:
        filter=lambda wildcards: PMD_MAP[wildcards.md],
    shell:
        """
            gunzip -c {input} |
                grep "{params.filter}$" >{output}
        """


# https://www.encodeproject.org/files/ENCFF526MRN/
# https://www.encodeproject.org/annotations/ENCSR814YSQ/
# Segway annotation of BC_COLON_H12817N
rule get_hmm_segmentation_colon_hg19:
    output:
        op.join(HG19_BASE, "hmm.bed.gz"),
    params:
        loc="https://www.encodeproject.org/files/ENCFF526MRN/@@download/ENCFF526MRN.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output}
        """


rule uncompress_and_filter_hmm_hg19:
    input:
        op.join(HG19_BASE, "hmm.bed.gz")
    output:
        temp(op.join(HG19_BASE, "{ann}.hmm.bed"))
    shell:
        """
            gunzip -c {input} | grep "{wildcards.ann}" > {output}
        """


CHIP_MAP = {
    "H3K27me3": "ENCFF255ARD",
    "H3K9me3": "ENCFF354CNQ",
    "H3K4me3": "ENCFF893TVK",
}


rule get_encode_chip_data_hg19:
    output:
        op.join(HG19_BASE, "{chip}.bed.gz")
    params:
        loc=lambda wildcards: f"https://www.encodeproject.org/files/{CHIP_MAP[wildcards.chip]}/@@download/{CHIP_MAP[wildcards.chip]}.bed.gz"
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output}
        """


rule uncompress_hgencode_chip_data_hg19:
    input:
        op.join(HG19_BASE, "{chip}.bed.gz")
    output:
        temp(op.join(HG19_BASE, "{chip}.chip.bed"))
    shell:
        """
            gunzip -c {input} > {output}
        """


rule get_lads_hg19:
    output:
        op.join(HG19_BASE, "laminb1.bed.gz"),
    params:
        loc="https://github.com/jernst98/ChromHMM/raw/refs/heads/master/COORDS/hg19/laminB1lads.hg19.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c > {output}
        """


rule uncompress_lads_hg19:
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
        op.join(HG19_BASE, "genome.sizes")
    params:
        asm="hg19",
    script:
        "src/download_genome_sizes.sh"
