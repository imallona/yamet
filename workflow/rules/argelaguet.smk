"""
Argelaguet et al. 2019 mouse gastrulation scNMT-seq: methylation entropy by feature.

Multi-omics profiling of mouse gastrulation at single cell resolution
Ricard Argelaguet, Stephen J Clark, Hisham Mohammed, et al.
https://pmc.ncbi.nlm.nih.gov/articles/PMC6924995
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121708

Assembly: GRCm38 (mm10). Reuses mm10 annotations from mm10.smk plus
gastrulation-specific ChIP-seq features bundled in the scnmt_gastrulation
tarball (E7.5 H3K27ac enhancers from GSE125318, ESC marks from ENCODE).
"""

ARGELAGUET_BASE          = "argelaguet"
ARGELAGUET_HARMONIZED    = op.join(ARGELAGUET_BASE, "harmonized")
ARGELAGUET_OUTPUT        = op.join(ARGELAGUET_BASE, "output")
ARGELAGUET_MAX_CELLS     = 20
ARGELAGUET_DOWNSAMPLE_SEED = 42

## Set True to restrict to chr10 for fast testing; False for full genome.
ARGELAGUET_CHR10_ONLY    = True

## GRCm38 (mm10) chromosomes. Defined explicitly rather than reusing the shared
## CHRS variable, which hg19.smk also defines (with human chromosomes 1-22).
_MM10_CHRS    = [str(i) for i in range(1, 20)] + ["X", "Y"]
_ARG_BED_CHRS = ["10"] if ARGELAGUET_CHR10_ONLY else _MM10_CHRS

## Generic mm10 annotations shared with ecker (all from mm10.smk rules).
_ARGELAGUET_MM10_ANNOTATIONS = {
    "chip":      ["h3k4me3", "h3k9me3", "h3k27me3", "h3k4me1", "h3k27ac"],
    "genes":     ["genes"],
    "lines":     ["lines"],
    "sines":     ["sines"],
    "promoters": ["promoters"],
}

## Gastrulation-specific features bundled in the scnmt_gastrulation tarball
## (features/genomic_contexts/). Keys are sanitized names used as the
## {annotation} wildcard (no underscores, no dots). Values are the source
## filenames inside the tarball.
##
## H3K27ac distal peaks: GSE125318 (companion ChIP-seq, see Methods).
## ESC marks: ENCODE.
_ARGELAGUET_GASTRO_BEDS = {
    "enh-E75-Ect":    "H3K27ac_distal_E7.5_Ect_intersect12.bed",
    "enh-E75-End":    "H3K27ac_distal_E7.5_End_intersect12.bed",
    "enh-E75-Mes":    "H3K27ac_distal_E7.5_Mes_intersect12.bed",
    "enh-E75-union":  "H3K27ac_distal_E7.5_union_intersect12.bed",
    "h3k4me3-E75-Ect":    "H3K4me3_E7.5_Ect.bed",
    "h3k4me3-E75-End":    "H3K4me3_E7.5_End.bed",
    "h3k4me3-E75-Mes":    "H3K4me3_E7.5_Mes.bed",
    "h3k4me3-E75-common": "H3K4me3_E7.5_common.bed",
    "esc-p300":        "ESC_p300.bed",
    "esc-dhs":         "ESC_DHS.bed",
}

_ARGELAGUET_GASTRO_ANNOTATIONS = {
    "enh_gastro":   [k for k in _ARGELAGUET_GASTRO_BEDS if k.startswith("enh-")],
    "h3k4me3_E75":  [k for k in _ARGELAGUET_GASTRO_BEDS if k.startswith("h3k4me3-")],
    "esc":          [k for k in _ARGELAGUET_GASTRO_BEDS if k.startswith("esc-")],
}

ARGELAGUET_ANNOTATIONS = {
    **_ARGELAGUET_MM10_ANNOTATIONS,
    **_ARGELAGUET_GASTRO_ANNOTATIONS,
}

## All annotation names that come from mm10.smk (no underscores, all lowercase).
_MM10_ANN_NAMES = [
    ann for cat in _ARGELAGUET_MM10_ANNOTATIONS
    for ann in _ARGELAGUET_MM10_ANNOTATIONS[cat]
]

## Plates without a stage prefix in sample_metadata.txt correspond to TET TKO
## cells from a companion experiment. Exclude them.
_ARGELAGUET_TET_KO_PLATES = [f"Plate{i}" for i in range(11, 17)]

ARGELAGUET_STRATIFY_BY = ["stage", "lineage10x"]


def _arg_sanitize(s):
    """Replace characters that conflict with filename delimiters."""
    return str(s).replace(" ", "-").replace(".", "-").replace("_", "-")


def get_argelaguet_groups():
    meta_fn = op.join(ARGELAGUET_BASE, "meta.tsv.gz")
    if not op.exists(meta_fn):
        return []
    meta = pd.read_csv(meta_fn, sep="\t", compression="gzip")
    available = [c for c in ARGELAGUET_STRATIFY_BY if c in meta.columns]
    if not available:
        return []
    combos = meta[available].dropna().drop_duplicates()
    return [tuple(_arg_sanitize(v) for v in row) for _, row in combos.iterrows()]


ARGELAGUET_GROUPS = get_argelaguet_groups()


rule download_argelaguet:
    output:
        flag=touch(op.join(ARGELAGUET_BASE, "downloaded.flag")),
    params:
        url="ftp://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/scnmt_gastrulation.tar.gz",
        base=ARGELAGUET_BASE,
    shell:
        """
        mkdir -p {params.base}
        curl {params.url} | tar -xz -C {params.base}/
        """


rule filter_argelaguet_metadata:
    input:
        flag=op.join(ARGELAGUET_BASE, "downloaded.flag"),
    output:
        meta=op.join(ARGELAGUET_BASE, "meta.tsv.gz"),
    params:
        raw_meta=op.join(ARGELAGUET_BASE, "sample_metadata.txt"),
        tet_ko_plates=_ARGELAGUET_TET_KO_PLATES,
    run:
        meta = pd.read_csv(params.raw_meta, sep="\t")
        meta = meta[meta["pass_metQC"] == True]
        meta = meta[~meta["plate"].isin(params.tet_ko_plates)]
        meta = meta.dropna(subset=["id_met", "stage", "lineage10x"])
        meta.to_csv(output.meta, sep="\t", index=False, compression="gzip")


checkpoint harmonize_argelaguet_cells:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        flag=op.join(ARGELAGUET_BASE, "downloaded.flag"),
        meta=op.join(ARGELAGUET_BASE, "meta.tsv.gz"),
    output:
        flag=touch(op.join(ARGELAGUET_BASE, "harmonized.flag")),
    params:
        raw=op.join(ARGELAGUET_BASE, "met", "cpg_level"),
        harmonized=ARGELAGUET_HARMONIZED,
        chr10_only=ARGELAGUET_CHR10_ONLY,
        threads=lambda wildcards, threads: threads,
    threads: 8
    script:
        "src/argelaguet_to_yamet.sh"


def get_argelaguet_harmonized_files(stage, lineage):
    import random
    checkpoints.harmonize_argelaguet_cells.get()
    meta = pd.read_csv(
        op.join(ARGELAGUET_BASE, "meta.tsv.gz"), sep="\t", compression="gzip"
    )
    mask = pd.Series([True] * len(meta), index=meta.index)
    for col, val in zip(ARGELAGUET_STRATIFY_BY, [stage, lineage]):
        if col in meta.columns:
            mask &= meta[col].astype(str).apply(_arg_sanitize) == val
    cell_ids = meta[mask]["id_met"].dropna().tolist()
    cells = [
        op.join(ARGELAGUET_HARMONIZED, f"{c}.gz")
        for c in cell_ids
        if op.exists(op.join(ARGELAGUET_HARMONIZED, f"{c}.gz"))
    ]
    if len(cells) > ARGELAGUET_MAX_CELLS:
        rng = random.Random(ARGELAGUET_DOWNSAMPLE_SEED)
        cells = rng.sample(cells, ARGELAGUET_MAX_CELLS)
    return cells


_ARG_BED_CHR_GREP = "|".join(f"^{c}\t" for c in _ARG_BED_CHRS)


rule argelaguet_filter_mm10_bed:
    """Subset a mm10.smk-generated annotation BED to the chromosomes in use."""
    conda:
        op.join("..", "envs", "processing.yml")
    wildcard_constraints:
        annotation="|".join(_MM10_ANN_NAMES),
    input:
        op.join(MM10_BASE, "{annotation}.bed"),
    output:
        temp(op.join(ARGELAGUET_BASE, "beds", "{annotation}.bed")),
    params:
        pattern=_ARG_BED_CHR_GREP,
    shell:
        "grep -E '{params.pattern}' {input} > {output}"


rule argelaguet_prep_gastro_bed:
    """Cut, sort, and merge a gastrulation-specific BED from the tarball."""
    conda:
        op.join("..", "envs", "processing.yml")
    wildcard_constraints:
        gastro_ann="|".join(_ARGELAGUET_GASTRO_BEDS.keys()),
    input:
        lambda wildcards: op.join(
            ARGELAGUET_BASE, "features", "genomic_contexts",
            _ARGELAGUET_GASTRO_BEDS[wildcards.gastro_ann],
        ),
    output:
        temp(op.join(ARGELAGUET_BASE, "beds", "{gastro_ann}.bed")),
    params:
        pattern=_ARG_BED_CHR_GREP,
        flag=op.join(ARGELAGUET_BASE, "downloaded.flag"),
    shell:
        """
        cut -f1-3 {input} |
          grep -E '{params.pattern}' |
          sort -k1,1 -k2,2n |
          bedtools merge -i - > {output}
        """


rule run_yamet_on_argelaguet_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        cells=lambda wildcards: get_argelaguet_harmonized_files(
            wildcards.stage, wildcards.lineage
        ),
        ref=op.join(MM10_BASE, "ref.CG.gz"),
        bed=op.join(ARGELAGUET_BASE, "beds", "{annotation}.bed"),
    output:
        simple_uncomp=temp(
            op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.out")
        ),
        det_uncomp=temp(
            op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.det.out")
        ),
        norm_det_uncomp=temp(
            op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.norm.det.out")
        ),
        meth_uncomp=temp(
            op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.meth.out")
        ),
        simple=op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.out.gz"),
        det=op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.det.out.gz"),
        meth=op.join(ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.meth.out.gz"),
        norm_det=op.join(
            ARGELAGUET_OUTPUT, "{annotation}_{stage}_{lineage}.norm.det.out.gz"
        ),
    log:
        op.join("logs", "yamet_argelaguet_{annotation}_{stage}_{lineage}.log"),
    params:
        path=ARGELAGUET_OUTPUT,
    threads: max(8, workflow.cores // 8)
    shell:
        """
        mkdir -p {params.path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.bed} \
         --cores {threads} \
         --no-print-sampens \
         --out {output.simple_uncomp} \
         --det-out {output.det_uncomp} \
         --meth-out {output.meth_uncomp} \
         --norm-det-out {output.norm_det_uncomp} &> {log}

        gzip --keep -f \
          {params.path}/{wildcards.annotation}_{wildcards.stage}_{wildcards.lineage}*out \
          &>> {log}
        """


def list_argelaguet_yamet_outputs(wildcards):
    checkpoints.harmonize_argelaguet_cells.get()
    res = []
    for cat in ARGELAGUET_ANNOTATIONS:
        for ann in ARGELAGUET_ANNOTATIONS[cat]:
            for stage, lineage in ARGELAGUET_GROUPS:
                if get_argelaguet_harmonized_files(stage, lineage):
                    res.append(f"{ann}_{stage}_{lineage}.det.out.gz")
    return [op.join(ARGELAGUET_OUTPUT, item) for item in res]


rule render_argelaguet_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        list_argelaguet_yamet_outputs,
        meta=op.join(ARGELAGUET_BASE, "meta.tsv.gz"),
    params:
        output_path=ARGELAGUET_OUTPUT,
        chr10_only=ARGELAGUET_CHR10_ONLY,
        meta_path=op.join(ARGELAGUET_BASE, "meta.tsv.gz"),
    threads:
        round(workflow.cores / 2)
    output:
        op.join(ARGELAGUET_BASE, "results", "argelaguet.html"),
    log:
        log=op.join("logs", "render_argelaguet.log"),
    script:
        "src/argelaguet.Rmd"
