"""
Reference file generation and mm10-specific annotation downloads.
"""

CHRS = [str(i) for i in range(1, 20)] + ["X", "Y"]
MM10_BASE = "mm10"


rule mm10_per_chr_ref:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        temp(op.join(MM10_BASE, "{chr}.{meth_pat}.ref")),
    log:
        op.join("logs", "mm10_per_chr_ref_{chr}_{meth_pat}.log"),
    params:
        fa="Mus_musculus.GRCm38.dna.chromosome.{chr}.fa",
        base="https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/",
        chr_prefix="",
    script:
        "src/build_chr_cpg_ref.sh"


rule mm10_aggregate_ref:
    wildcard_constraints:
        meth_pat="[^.]+",
    input:
        lambda wildcards: expand(
            op.join(MM10_BASE, "{chr}.{{meth_pat}}.ref"),
            chr=CHRS,
        ),
    output:
        op.join(MM10_BASE, "ref.{meth_pat}.gz"),
    shell:
        "sort -m -k1,1 -k2,2n {input} | gzip > {output}"


rule mm10_chr10_aggregate_ref:
    input:
        op.join(MM10_BASE, "10.{meth_pat}.ref"),
    output:
        op.join(MM10_BASE, "ref.{meth_pat}.chr10.gz"),
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


rule get_mm10_sizes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        op.join(MM10_BASE, "mm10.sizes")
    shell:
        """
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
          'SELECT chrom,size FROM mm10.chromInfo' |
          sed 's/^chr//' |
          grep -v '_' > {output}
        """


rule make_windows_mm10:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        sizes=op.join(MM10_BASE, "mm10.sizes"),
    output:
        op.join(MM10_BASE, "windows_{win_size}_nt.bed")
    shell:
        "bedtools makewindows -g {input.sizes} -w {wildcards.win_size} | sort -k1,1 -k2,2n > {output}"


rule make_windows_mm10_chr10:
    input:
        op.join(MM10_BASE, "windows_{win_size}_nt.bed")
    output:
        op.join(MM10_BASE, "windows_{win_size}_nt.chr10.bed")
    shell:
        "grep '^10\t' {input} > {output}"


rule decompress_mm10_annotation:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(MM10_BASE, "{annotation}.bed.gz"),
    output:
        temp(op.join(MM10_BASE, "{annotation}.bed")),
    shell:
        "zcat {input} | sed 's/^chr//' | sort -k1,1 -k2,2n | bedtools merge -i - > {output}"
