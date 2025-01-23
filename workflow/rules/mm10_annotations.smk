"""
Retrieves mm10 annotations (for Ecker)
"""


rule download_mm10_annotation_for_ecker:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        promoters = op.join('annotation', 'mm10', 'promoters.bed.gz'),
        txs = op.join('annotation', 'mm10', 'genes.bed.gz'),
        encode = op.join('annotation', 'mm10', 'h3k27ac.bed.gz')
    output:
        flag = op.join('annotation', 'mm10', 'done.flag')
    params:
        path = op.join('annotation', 'mm10')
    threads:
        1
    shell:
        """
        date > {output.flag}
        """
    
rule download_mm10_promoters:
    conda:
        op.join("..","envs", "yamet.yml")
    output:
        promoters = op.join('annotation', 'mm10', 'promoters.bed.gz')    
    params:
        path = op.join('annotation', 'mm10')
    threads:
        1
    shell:
        """
        bash src/download_mm10_promoters.sh | gzip -c > {output.promoters}
        """

rule download_mm10_genes:
    conda:
        op.join("..","envs", "yamet.yml")
    output:
        txs = op.join('annotation', 'mm10', 'genes.bed.gz')
    params:
        path = op.join('annotation', 'mm10')
    threads:
        1
    shell:
        """
        bash src/download_mm10_genes.sh | gzip -c > {output.txs}
        """
        
rule download_mm10_encode:
    conda:
        op.join("..","envs", "yamet.yml")
    output:
        last = op.join('annotation', 'mm10', 'h3k27ac.bed.gz')
    params:
        path = op.join('annotation', 'mm10')
    threads:
        1
    shell:
        """
        bash src/download_encode_mm10_epigenomics.sh
        """

rule get_mm10_genome_sizes:
    conda:
        op.join("..", "envs", "yamet.yml")
    output:
        genome = op.join('annotation', 'mm10','genome.sizes')
    params:
        path = op.join('annotation', 'mm10')
    threads:
        1
    shell:
        """
        bash src/download_mm10_genome_sizes.sh > {output.genome}
        """
