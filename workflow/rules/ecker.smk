"""
Ecker's related rules

DNA methylation atlas of the mouse brain at single-cell resolution

Hanqing Liu, Jingtian Zhou, Wei Tian, Chongyuan Luo, Anna Bartlett, Andrew Aldridge, Jacinta Lucero, Julia K. Osteen, Joseph R. Nery, Huaming Chen, Angeline Rivkin, Rosa G. Castanon, Ben Clock, Yang Eric Li, Xiaomeng Hou, Olivier B. Poirion, Sebastian Preissl, Antonio Pinto-Duarte, Carolyn O’Connor, Lara Boggeman, Conor Fitzpatrick, Michael Nunn, Eran A. Mukamel, Zhuzhu Zhang, …Joseph R. Ecker

https://www.nature.com/articles/s41586-020-03182-8

this is mm10
"""

ECKER_BASE = "ecker"
ECKER_HARMONIZED = op.join(ECKER_BASE, "harmonized")
ECKER_OUTPUT = op.join(ECKER_BASE, "output")
ECKER_WINDOWS_OUTPUT = op.join(ECKER_BASE, "windows_output")

## maximum cells per (major_region, cell_class) group passed to yamet
ECKER_MAX_CELLS = 20
ECKER_DOWNSAMPLE_SEED = 42

## set to True to restrict to chr10 for speed; False for full genome
ECKER_CHR10_ONLY = False

## chromosomes to retain in annotation BED files for Ecker runs
_ECKER_BED_CHRS = ["10"] if ECKER_CHR10_ONLY else CHRS

_ECKER_REF = (
    op.join(MM10_BASE, "ref.CG.chr10.gz")
    if ECKER_CHR10_ONLY
    else op.join(MM10_BASE, "ref.CG.gz")
)

ECKER_ANNOTATIONS = {
    "chip": ["h3k4me3", "h3k9me3", "h3k27me3", "h3k4me1", "h3k27ac"],
    "genes": ["genes"],
    "lines": ["lines"],
    "sines": ["sines"],
    "promoters": ["promoters"],
}

## Columns to stratify cells by (intersection of all listed columns).
## Change this list to stratify by other metadata column combinations;
## update the wildcard names in run_yamet_on_ecker_features to match.
ECKER_STRATIFY_BY = ["SubRegion", "SubType"]


def _sanitize(s):
    return str(s).replace(" ", "-")


def get_ecker_groups():
    meta_fn = op.join(ECKER_BASE, "meta.tsv.gz")
    if not op.exists(meta_fn):
        return []
    meta = pd.read_csv(meta_fn, sep="\t", compression="gzip")
    available = [c for c in ECKER_STRATIFY_BY if c in meta.columns]
    if not available:
        return []
    combos = meta[available].dropna().drop_duplicates()
    return [tuple(_sanitize(v) for v in row) for _, row in combos.iterrows()]


## list of (sub_region, sub_type) tuples
ECKER_GROUPS = get_ecker_groups()


rule download_nemo_ecker_metadata:
    output:
        meta=temp(op.join(ECKER_BASE, "nemo_meta.tsv.gz")),
    params:
        loc="https://data.nemoarchive.org/biccn/grant/u19_cemba/cemba/epigenome/sncell/mCseq/mouse/processed/analysis/EckerRen_Mouse_MOp_methylation_ATAC/metadata/mc/MOp_Metadata.tsv.gz",
    shell:
        """
            curl {params.loc} -o {output.meta}
        """


rule download_ecker_paper_metadata:
    conda:
        "../envs/yamet.yml"
    output:
        meta=temp(op.join(ECKER_BASE, "paper_meta.xlsx")),
    params:
        loc="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-03182-8/MediaObjects/41586_2020_3182_MOESM9_ESM.xlsx",
    shell:
        """
            curl {params.loc} -o {output.meta}
        """


rule harmonize_ecker_metadata:
    conda:
        "../envs/r.yml"
    input:
        nemo=op.join(ECKER_BASE, "nemo_meta.tsv.gz"),
        paper=op.join(ECKER_BASE, "paper_meta.xlsx"),
    output:
        metadata=op.join(ECKER_BASE, "meta.tsv.gz"),
    script:
        "src/harmonize_ecker_metadata.R"


rule ecker_urls:
    input:
        op.join(ECKER_BASE, "meta.tsv.gz"),
    output:
        temp(op.join(ECKER_BASE, "urls")),
    params:
        base_url="https://data.nemoarchive.org/biccn/grant/u19_cemba/cemba/epigenome/sncell/mCseq/mouse/processed/counts/",
    shell:
        """
            zcat {input[0]} | cut -f29 | grep -v AllcPath >ecker_raw_urls
            awk -F'/allc/' '{{gsub(/"/, "", $0); split($(NF-1), a, "/"); sub(/\\.tsv\\.gz$/, ".tsv.tar", $NF); print "{params[0]}" a[length(a)-1]"/"a[length(a)]"/"$NF}}' ecker_raw_urls > {output[0]}
            rm ecker_raw_urls
        """


rule download_ecker:
    conda:
        "../envs/processing.yml"
    input:
        op.join(ECKER_BASE, "urls"),
    output:
        flag=touch(op.join(ECKER_BASE, "downloaded.flag")),
    params:
        raw=op.join(ECKER_BASE, "raw"),
    threads: 5
    script:
        "src/download_ecker_tars.sh"


checkpoint harmonize_ecker_cells:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        download_flag=op.join(ECKER_BASE, "downloaded.flag"),
        meta=op.join(ECKER_BASE, "meta.tsv.gz"),
    output:
        flag=touch(op.join(ECKER_HARMONIZED, "done.flag")),
    params:
        raw=op.join(ECKER_BASE, "raw"),
        harmonized=ECKER_HARMONIZED,
        chr10_only=ECKER_CHR10_ONLY,
        threads=lambda wildcards, threads: threads,
    threads: 8
    script:
        "src/ecker_allc_to_yamet.sh"


def get_ecker_harmonized_files(sub_region, sub_type):
    import random
    checkpoints.harmonize_ecker_cells.get()
    meta = pd.read_csv(
        op.join(ECKER_BASE, "meta.tsv.gz"), sep="\t", compression="gzip"
    )
    mask = pd.Series([True] * len(meta), index=meta.index)
    for col, val in zip(ECKER_STRATIFY_BY, [sub_region, sub_type]):
        if col in meta.columns:
            mask &= meta[col].astype(str).apply(_sanitize) == val
    cell_basenames = meta[mask]["basename"].dropna().tolist()
    cells = [
        op.join(ECKER_HARMONIZED, c)
        for c in cell_basenames
        if op.exists(op.join(ECKER_HARMONIZED, c))
    ]
    if len(cells) > ECKER_MAX_CELLS:
        rng = random.Random(ECKER_DOWNSAMPLE_SEED)
        cells = rng.sample(cells, ECKER_MAX_CELLS)
    return cells


_ECKER_BED_CHR_GREP = "|".join(f"^{c}\t" for c in _ECKER_BED_CHRS)


rule ecker_filter_mm10_bed:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(MM10_BASE, "{annotation}.bed"),
    output:
        temp(op.join(ECKER_BASE, "beds", "{annotation}.bed")),
    params:
        pattern=_ECKER_BED_CHR_GREP,
    shell:
        "grep -E '{params.pattern}' {input} > {output}"


rule run_yamet_on_ecker_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        cells=lambda wildcards: get_ecker_harmonized_files(wildcards.sub_region, wildcards.sub_type),
        validation=ancient(op.join(ECKER_HARMONIZED, "coords_validated.flag")),
        ref=_ECKER_REF,
        bed=op.join(ECKER_BASE, "beds", "{annotation}.bed"),
    output:
        simple_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.out")),
        det_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.det.out")),
        norm_det_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.norm.det.out")),
        meth_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.meth.out")),
        simple=op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.out.gz"),
        det=op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.det.out.gz"),
        meth=op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.meth.out.gz"),
        norm_det=op.join(ECKER_OUTPUT, "{annotation}_{sub_region}_{sub_type}.norm.det.out.gz"),
    log:
        op.join("logs", "yamet_ecker_{annotation}_{sub_region}_{sub_type}.log"),
    params:
        path=ECKER_OUTPUT,
    threads: min(workflow.cores, 8)
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

        gzip --keep -f {params.path}/{wildcards.annotation}_{wildcards.sub_region}_{wildcards.sub_type}*out &>> {log}
        """


def list_ecker_yamet_outputs(wildcards):
    checkpoints.harmonize_ecker_cells.get()
    meta = pd.read_csv(
        op.join(ECKER_BASE, "meta.tsv.gz"), sep="\t", compression="gzip"
    )
    available = [c for c in ECKER_STRATIFY_BY if c in meta.columns]
    combos = [
        tuple(_sanitize(v) for v in row)
        for _, row in meta[available].dropna().drop_duplicates().iterrows()
    ]
    res = []
    for cat in ECKER_ANNOTATIONS:
        for ann in ECKER_ANNOTATIONS[cat]:
            for sub_region, sub_type in combos:
                if get_ecker_harmonized_files(sub_region, sub_type):
                    res.append(f"{ann}_{sub_region}_{sub_type}.det.out.gz")
    return [op.join(ECKER_OUTPUT, item) for item in res]


def get_ecker_downsampled_cells(wildcards):
    """Return the union of per-group downsampled cells (same as feature runs)."""
    checkpoints.harmonize_ecker_cells.get()
    meta = pd.read_csv(
        op.join(ECKER_BASE, "meta.tsv.gz"), sep="\t", compression="gzip"
    )
    available = [c for c in ECKER_STRATIFY_BY if c in meta.columns]
    combos = [
        tuple(_sanitize(v) for v in row)
        for _, row in meta[available].dropna().drop_duplicates().iterrows()
    ]
    cells = []
    for sub_region, sub_type in combos:
        cells.extend(get_ecker_harmonized_files(sub_region, sub_type))
    return cells


rule run_yamet_on_ecker_windows:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        cells=get_ecker_downsampled_cells,
        validation=ancient(op.join(ECKER_HARMONIZED, "coords_validated.flag")),
        ref=op.join(MM10_BASE, "ref.CG.gz"),
        windows=op.join(MM10_BASE, "windows_{win_size}_nt.bed"),
    output:
        det_tmp=temp(op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.det.out")),
        norm_det_tmp=temp(op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.norm.det.out")),
        meth_tmp=temp(op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.meth.out")),
        simple_tmp=temp(op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.out")),
        det=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.det.out.gz"),
        norm_det=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.norm.det.out.gz"),
        meth=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.meth.out.gz"),
        simple=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.out.gz"),
    params:
        path=ECKER_WINDOWS_OUTPUT,
        prefix=lambda wildcards: op.join(ECKER_WINDOWS_OUTPUT, f"{wildcards.win_size}_all"),
    log:
        op.join("logs", "yamet_ecker_windows_{win_size}.log"),
    threads: workflow.cores
    shell:
        """
        mkdir -p {params.path}
        yamet \
         --cell {input.cells} \
         --reference {input.ref} \
         --intervals {input.windows} \
         --cores {threads} \
         --no-print-sampens \
         --out {params.prefix}.out \
         --det-out {params.prefix}.det.out \
         --meth-out {params.prefix}.meth.out \
         --norm-det-out {params.prefix}.norm.det.out &> {log}
        gzip --keep -f \
          {params.prefix}.out \
          {params.prefix}.det.out \
          {params.prefix}.meth.out \
          {params.prefix}.norm.det.out &>> {log}
        """


rule render_ecker_windows_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        det=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.det.out.gz"),
        norm_det=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.norm.det.out.gz"),
        meth=op.join(ECKER_WINDOWS_OUTPUT, "{win_size}_all.meth.out.gz"),
        meta=op.join(ECKER_BASE, "meta.tsv.gz"),
    output:
        op.join(ECKER_BASE, "results", "ecker_windows_{win_size}.html"),
    params:
        corrected_sce=op.join(
            ECKER_BASE, "results", "sce_windows_ecker_{win_size}_corrected.rds"
        ),
    threads:
        workflow.cores
    log:
        log=op.join("logs", "render_ecker_windows_{win_size}.log"),
    script:
        "src/ecker_windows.Rmd"


rule render_ecker_embeddings_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        windows_html=op.join(
            ECKER_BASE, "results", "ecker_windows_{win_size}.html"
        ),
    output:
        op.join(ECKER_BASE, "results", "ecker_embeddings_{win_size}.html"),
    params:
        corrected_sce=op.join(
            ECKER_BASE, "results", "sce_windows_ecker_{win_size}_corrected.rds"
        ),
    threads:
        workflow.cores
    log:
        log=op.join("logs", "render_ecker_embeddings_{win_size}.log"),
    script:
        "src/ecker_embeddings.Rmd"


rule render_ecker_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        list_ecker_yamet_outputs,
        meta=op.join(ECKER_BASE, "meta.tsv.gz"),
    params:
        output_path=ECKER_OUTPUT,
        chr10_only=ECKER_CHR10_ONLY,
        meta_path=op.join(ECKER_BASE, "meta.tsv.gz"),
    threads:
        round(workflow.cores / 2)
    output:
        op.join(ECKER_BASE, "results", "ecker.html"),
    log:
        log=op.join("logs", "render_ecker.log"),
    script:
        "src/ecker.Rmd"
