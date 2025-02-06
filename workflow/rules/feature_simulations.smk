"""
Simulates cells with pre-defined background entropy statuses with interspersed, known genomic regions of 
  differential entropy.
"""

rule simulate_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    output:
        cell = op.join('feature_simulations', 'output', 'simulations.tsv'),
        reference = op.join('feature_simulations', 'output', 'reference.tsv'),
        regions = op.join('feature_simulations', 'output', 'regions.bed')
    params:
        path = op.join('feature_simulations', 'output')        
    shell:
        """
        mkdir -p {params.path}
        Rscript src/feature_simulation.R
        mv regions.bed {params.path}
        mv simulations.tsv {params.path}
        mv reference.tsv {params.path}
        """
