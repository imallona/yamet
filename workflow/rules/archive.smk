"""
Archive rules: pack yamet output directories into zip files for Zenodo upload
or for restarting the workflow from intermediate outputs (mid-way start).

Files in output/ are already gzip-compressed; zip is used without re-compression
(-0) to avoid wasting time on already-compressed data.

Usage (from workflow/):
    snakemake --cores 1 snapshots/ecker_output.zip snapshots/crc_output.zip
    snakemake --cores 1 snapshots/argelaguet_output.zip
    snakemake --cores 1 snapshots/simulation_output.zip
"""

SNAPSHOT_DIR = config.get("snapshot_dir", "snapshots")


rule archive_ecker_outputs:
    """Zip all ecker yamet outputs (the bare minimum to re-run the report)."""
    input:
        list_ecker_yamet_outputs,
    output:
        op.join(SNAPSHOT_DIR, "ecker_output.zip"),
    shell:
        """
        mkdir -p {SNAPSHOT_DIR}
        cd ecker && zip -0 $OLDPWD/{output} output/*.gz
        """


rule archive_crc_outputs:
    """Zip all CRC yamet outputs (the bare minimum to re-run the report)."""
    input:
        list_relevant_yamet_outputs(),
    output:
        op.join(SNAPSHOT_DIR, "crc_output.zip"),
    shell:
        """
        mkdir -p {SNAPSHOT_DIR}
        cd data/crc && zip -0 $OLDPWD/{output} output/*.gz
        """


rule archive_argelaguet_outputs:
    """Zip argelaguet yamet outputs and metadata for quick report rendering from intermediate results from zenodo."""
    input:
        list_argelaguet_yamet_outputs,
        meta=op.join(ARGELAGUET_BASE, "meta.tsv.gz"),
    output:
        op.join(SNAPSHOT_DIR, "argelaguet_output.zip"),
    shell:
        """
        mkdir -p {SNAPSHOT_DIR}
        cd {ARGELAGUET_BASE} && zip -0 $OLDPWD/{output} output/*.gz meta.tsv.gz
        """


rule archive_simulation_outputs:
    """Zip simulation intermediate outputs for quick report rendering from intermediate results from zenodo."""
    input:
        op.join(SIM_RESULTS, "simulation_figure2.html"),
        op.join(SIM_RESULTS, "simulation_08_combined_figure_adj.html"),
    output:
        op.join(SNAPSHOT_DIR, "simulation_output.zip"),
    shell:
        """
        mkdir -p {SNAPSHOT_DIR}
        cd {SIM_BASE} && zip -0r $OLDPWD/{output} output/ results/ schemes/ rds/
        """
