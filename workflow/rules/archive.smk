"""
Archive rules: pack yamet output directories into zip files for Zenodo upload.

Files in output/ are already gzip-compressed; zip is used without re-compression
(-0) to avoid wasting time on already-compressed data.

Usage (from workflow/):
    snakemake --cores 1 snapshots/ecker_output.zip snapshots/crc_output.zip
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
