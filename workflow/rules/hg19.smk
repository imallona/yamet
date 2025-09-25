"""
To generate a reference file from hg19 assembly and download corresponding annotation files
"""

CHRS = [str(i) for i in range(1, 23)] + ["X", "Y"]
HG19_BASE = "hg19"

## bedfiles
ANNOTATIONS = {
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
}


rule build_hg19_chr_per_chr:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(HG19_BASE, "{chr}.{meth_pat}.ref")),
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


rule get_sizes_hg19:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "hg19.sizes"),
    params:
        asm="hg19",
    script:
        "src/download_genome_sizes.sh"


rule make_windows_hg19:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(HG19_BASE, "hg19.sizes"),
    output:
        op.join(HG19_BASE, "windows_{win_size}_nt.bed.gz"),
    shell:
        """
        bedtools makewindows -g {input} -w {wildcards.win_size} |
            sort -k1,1 -k2,2n | gzip -c >{output}
        """


rule get_single_annotation_coverage_per_window:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        windows=op.join(HG19_BASE, "windows_{win_size}_nt.bed.gz"),
        annotation=op.join(HG19_BASE, "{subcat}.{cat}.bed.gz"),
    output:
        temp(op.join(HG19_BASE, "windows_{win_size}_nt_{subcat}_{cat}_annotation.frac")),
    shell:
        """
        echo "{wildcards.subcat}_{wildcards.cat}" >{output}
        bedtools coverage \
            -a <(gunzip -c {input.windows}) \
            -b <(gunzip -c {input.annotation}) \
            | cut -f7 >>{output}
        """


rule get_amplification_coverage_per_window:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        windows=op.join(HG19_BASE, "windows_{win_size}_nt.bed.gz"),
        annotation=op.join(HG19_BASE, "scna.bed.gz"),
    output:
        temp(op.join(HG19_BASE, "windows_{win_size}_nt_scna_annotation.frac")),
    shell:
        """
        bedtools intersect \
            -a <(gunzip -c {input.windows}) \
            -b <(gunzip -c {input.annotation}) -wao |
            awk 'BEGIN{{OFS="\t"; print "scna_frac","scna_status"}}
                {{
                    win_len = $3 - $2
                    overlap = $NF
                    frac = (win_len > 0 ? overlap / win_len : 0)
                    status = $(NF-1)
                    if (status == ".") status = "NA"
                    print frac, status
                }}' >{output}
        """


def list_annotated_windows():
    res = []
    for cat in ANNOTATIONS:
        for subcat in ANNOTATIONS[cat]:
            res.append(f"windows_{{win_size}}_nt_{subcat}_{cat}_annotation.frac")
    return [op.join("hg19", item) for item in res]


rule combine_annotated_windows:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        annotated_windows=list_annotated_windows(),
        scna=op.join("hg19", "windows_{win_size}_nt_scna_annotation.frac"),
    output:
        op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz"),
    shell:
        """
        paste {input.annotated_windows} {input.scna} | gzip -c >{output}
        """


rule get_scna_hg19:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        op.join("src", "scna_hg19.xlsx"),
    output:
        op.join(HG19_BASE, "scna.bed.gz"),
    script:
        "src/parse_scna.R"


rule get_genes_hg19:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(HG19_BASE, "genes.genes.bed.gz"),
    script:
        "src/download_hg19_genes.sh"


rule get_pmds_file_hg19:
    output:
        temp(op.join(HG19_BASE, "pmd.bed.gz")),
    params:
        loc="https://zhouserver.research.chop.edu/GenomeAnnotation/hg19/PMD_coordinates_hg19.bed.gz",
    shell:
        """
        curl -o {output[0]} {params.loc}
        """


PMD_MAP = {"pmds": "commonPMD", "hmds": "commonHMD"}


rule get_pmds_hmds_hg19:
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


rule get_hmm_file_hg19:
    output:
        temp(op.join(HG19_BASE, "hmm.bed")),
    params:
        loc="https://www.encodeproject.org/files/ENCFF526MRN/@@download/ENCFF526MRN.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n >{output}
        """


rule get_hmm_hg19:
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


rule get_encode_chip_data_hg19:
    output:
        op.join(HG19_BASE, "{chip}.chip.bed.gz"),
    params:
        loc=lambda wildcards: f"https://www.encodeproject.org/files/{CHIP_MAP[wildcards.chip]}/@@download/{CHIP_MAP[wildcards.chip]}.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c >{output}
        """


rule get_lads_hg19:
    output:
        op.join(HG19_BASE, "{lamin}.lad.bed.gz"),
    params:
        loc="https://github.com/jernst98/ChromHMM/raw/refs/heads/master/COORDS/hg19/laminB1lads.hg19.bed.gz",
    shell:
        """
            curl -L {params.loc} | gunzip -c | sort -k1,1 -k2,2n | gzip -c >{output}
        """
