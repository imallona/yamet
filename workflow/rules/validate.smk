"""
Pre-validation rules: check that harmonized cell coordinates match the yamet
reference before any yamet run.  One rule per dataset; each writes a flag that
the corresponding yamet rule depends on.
"""

_VALIDATE_SCRIPT = op.join("src", "validate_yamet_coords.sh")


rule validate_argelaguet_coords:
    input:
        flag=op.join(ARGELAGUET_HARMONIZED, "done.flag"),
        ref=op.join(MM10_BASE, "ref.CG.gz"),
    output:
        flag=touch(op.join(ARGELAGUET_HARMONIZED, "coords_validated.flag")),
    params:
        harmonized=ARGELAGUET_HARMONIZED,
    shell:
        "bash {_VALIDATE_SCRIPT} {params.harmonized} {input.ref}"


rule validate_ecker_coords:
    input:
        flag=op.join(ECKER_HARMONIZED, "done.flag"),
        ref=_ECKER_REF,
    output:
        flag=touch(op.join(ECKER_HARMONIZED, "coords_validated.flag")),
    params:
        harmonized=ECKER_HARMONIZED,
    shell:
        "bash {_VALIDATE_SCRIPT} {params.harmonized} {input.ref}"


rule validate_crc_coords:
    input:
        flag=op.join(CRC, "download.flag"),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
    output:
        flag=touch(op.join(CRC_HARMONIZED, "coords_validated.flag")),
    params:
        harmonized=CRC_HARMONIZED,
    shell:
        "bash {_VALIDATE_SCRIPT} {params.harmonized} {input.ref}"
