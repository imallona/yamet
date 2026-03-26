"""
Reference file generation and mm10-specific annotation downloads.
"""

CHRS = [str(i) for i in range(1, 20)] + ["X", "Y"]
MM10_BASE = "mm10"

## Chromosomes used when building the CG reference.
## ecker.smk overrides this to ["10"] when ECKER_CHR10_ONLY is True so only
## chr10 is downloaded and aggregated, avoiding a full-genome sort.
MM10_CG_CHRS = CHRS

## Chromosomes to keep when decompressing annotation BED files.
## ecker.smk overrides this to ["10"] when ECKER_CHR10_ONLY is True.
MM10_BED_CHRS = CHRS


rule mm10_per_chr_ref:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(MM10_BASE, "{chr}.{meth_pat}.ref")),
    params:
        fa="Mus_musculus.GRCm38.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/",
        chr_prefix="",
    script:
        "src/build_chr_cpg_ref.sh"


rule mm10_aggregate_ref:
    input:
        lambda wildcards: expand(
            op.join(MM10_BASE, "{chr}.{{meth_pat}}.ref"),
            chr=(MM10_CG_CHRS if wildcards.meth_pat == "CG" else CHRS),
        ),
    output:
        op.join(MM10_BASE, "ref.{meth_pat}.gz"),
    shell:
        "sort -m -k1,1 -k2,2n {input} | gzip > {output}"


rule get_mm10_promoters:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(MM10_BASE, "promoters.bed.gz"),
    shell:
        """
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
              'SELECT chrom, min(txStart), max(txEnd), name2, strand
               FROM mm10.wgEncodeGencodeBasicVM25
               GROUP BY name2
               ORDER by chrom, min(txStart);' |
              awk '{{OFS=FS="\\t"; {{print $1, ($2+1-2000 < 0 ? 0 : $2+1-2000), $2+2000, $4, ".", $5}}}}' |
              gzip -c >{output}
        """


ENCODE_MAP = {
    "h3k4me3": "ENCFF160SCR",
    "h3k9me3": "ENCFF658QTP",
    "h3k4me1": "ENCFF937JHP",
    "h3k27me3": "ENCFF827BBC",
    "h3k27ac": "ENCFF442GIT",
}


rule get_mm10_encode:
    wildcard_constraints:
        ann="|".join(ENCODE_MAP.keys()),
    output:
        op.join(MM10_BASE, "{ann}.bed.gz"),
    params:
        accessor=lambda wildcards: ENCODE_MAP[wildcards.ann],
    shell:
        """
        curl -L https://www.encodeproject.org/files/{params.accessor}/@@download/{params.accessor}.bed.gz -o {output}
        """


rule decompress_mm10_annotation:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(MM10_BASE, "{annotation}.bed.gz"),
    output:
        temp(op.join(MM10_BASE, "{annotation}.bed")),
    params:
        chrs=" ".join(MM10_BED_CHRS),
    shell:
        """zcat {input} | sed 's/^chr//' | awk 'BEGIN{{n=split("{params.chrs}",a); for(i=1;i<=n;i++) k[a[i]]=1}} $1 in k' | sort -k1,1 -k2,2n | bedtools merge -i - > {output}"""
