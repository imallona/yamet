"""
Assembly-agnostic annotation download rules, shared across hg19 and mm10.
"""

## Gencode table names by assembly, used in get_genes.
GENCODE_TABLE = {
    "mm10": "wgEncodeGencodeBasicVM25",
    "hg19": "wgEncodeGencodeBasicV19",
}

## Exclude genes longer than 100 kb to remove telomeric/superlong entries.
GENE_LENGTH_FILTER = "HAVING MAX(txEnd) - MIN(txStart) < 100000"


rule get_genome_sizes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        "{assembly}/genome.sizes",
    shell:
        """
        mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
              "SELECT chrom,size FROM {wildcards.assembly}.chromInfo" |
              grep -E -w "chr[0-9XY]{{1,2}}" |
              sort -k1,1 >{output}
        """


rule get_genes:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        "{assembly}/genes.bed.gz",
    params:
        table=lambda wildcards: f"{wildcards.assembly}.{GENCODE_TABLE[wildcards.assembly]}",
    shell:
        """
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
          "SELECT chrom, MIN(txStart), MAX(txEnd), name2, strand
           FROM {params.table}
           GROUP BY name2
           {GENE_LENGTH_FILTER}
           ORDER BY chrom, MIN(txStart);" |
            awk '{{OFS=FS="\\t"; {{print $1, $2, $3, $4, ".", $5}}}}' |
            sort -k1,1 -k2,2n |
            bedtools merge -i - |
            gzip -c >{output}
        """


rule get_rmsk_lines:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        "{assembly}/lines.bed.gz",
    shell:
        """
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
          "SELECT genoName, genoStart, genoEnd, repName, repClass, repFamily
             FROM {wildcards.assembly}.rmsk
             WHERE repClass = 'LINE';" |
            awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $4, $5, $6}}' |
            sort -k1,1 -k2,2n |
            bedtools merge -i - |
            gzip -c >{output}
        """


rule get_rmsk_sines:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        "{assembly}/sines.bed.gz",
    shell:
        """
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
          "SELECT genoName, genoStart, genoEnd, repName, repClass, repFamily
             FROM {wildcards.assembly}.rmsk
             WHERE repClass = 'SINE';" |
            awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $4, $5, $6}}' |
            sort -k1,1 -k2,2n |
            bedtools merge -i - |
            gzip -c >{output}
        """
