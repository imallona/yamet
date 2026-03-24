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

## set to True to restrict to chr10 for speed; False for full genome
ECKER_CHR10_ONLY = True

## tell mm10.smk which chromosomes to include in the CG reference
MM10_CG_CHRS = ["10"] if ECKER_CHR10_ONLY else CHRS

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
ECKER_STRATIFY_BY = ["MajorRegion", "CellClass"]


def get_ecker_groups():
    meta_fn = op.join(ECKER_BASE, "meta.tsv.gz")
    if not op.exists(meta_fn):
        return []
    meta = pd.read_csv(meta_fn, sep="\t", compression="gzip")
    available = [c for c in ECKER_STRATIFY_BY if c in meta.columns]
    if not available:
        return []
    combos = meta[available].dropna().drop_duplicates()
    return [tuple(str(v) for v in row) for _, row in combos.iterrows()]


## list of (major_region, cell_class) tuples
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


# ## reports the AllcPath (basename) of cells matching the harmonized metadata
# ##   'column' equals 'value'
# def slice_eckers_metadata(column, value):
#     meta_fn = op.join('ecker_data','harmonized_ecker_metadata.tsv.gz')

#     if not op.exists(meta_fn):
#         raise Exception("No metadata found.")

#     meta = pd.read_csv(filepath_or_buffer = meta_fn,
#                        sep='\t',
#                        compression='gzip', header=0, quotechar='"')

#     if set([column]).issubset(meta.columns):
#         return [i for i in meta[meta[column] == value]['basename']]
#     else:
#         return None

# ## this is a long one!
# print(slice_eckers_metadata('SubType', 'IT-L4 Shc3'))


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
        "src/get_ecker_meth_files.sh"


checkpoint harmonize_ecker_cells:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        download_flag=op.join(ECKER_BASE, "downloaded.flag"),
        meta=op.join(ECKER_BASE, "meta.tsv.gz"),
    output:
        flag=touch(op.join(ECKER_BASE, "harmonized.flag")),
    params:
        raw=op.join(ECKER_BASE, "raw"),
        harmonized=ECKER_HARMONIZED,
        chr10_only=ECKER_CHR10_ONLY,
        threads=lambda wildcards, threads: threads,
    threads: 8
    script:
        "src/harmonize_ecker_cells.sh"


def get_ecker_harmonized_files(major_region, cell_class):
    checkpoints.harmonize_ecker_cells.get()
    meta = pd.read_csv(
        op.join(ECKER_BASE, "meta.tsv.gz"), sep="\t", compression="gzip"
    )
    mask = pd.Series([True] * len(meta), index=meta.index)
    for col, val in zip(ECKER_STRATIFY_BY, [major_region, cell_class]):
        if col in meta.columns:
            mask &= meta[col].astype(str) == val
    cell_basenames = meta[mask]["basename"].dropna().tolist()
    return [
        op.join(ECKER_HARMONIZED, c)
        for c in cell_basenames
        if op.exists(op.join(ECKER_HARMONIZED, c))
    ]


rule run_yamet_on_ecker_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        cells=lambda wildcards: get_ecker_harmonized_files(wildcards.major_region, wildcards.cell_class),
        ref=op.join(MM10_BASE, "ref.CG.gz"),
        bed=op.join(MM10_BASE, "{annotation}.bed"),
    output:
        simple_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.out")),
        det_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.det.out")),
        norm_det_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.norm.det.out")),
        meth_uncomp=temp(op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.meth.out")),
        simple=op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.out.gz"),
        det=op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.det.out.gz"),
        meth=op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.meth.out.gz"),
        norm_det=op.join(ECKER_OUTPUT, "{annotation}_{major_region}_{cell_class}.norm.det.out.gz"),
    log:
        op.join("logs", "yamet_ecker_{annotation}_{major_region}_{cell_class}.log"),
    params:
        path=ECKER_OUTPUT,
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

        gzip --keep -f {params.path}/{wildcards.annotation}_{wildcards.major_region}_{wildcards.cell_class}*out &>> {log}
        """


def list_ecker_yamet_outputs():
    res = []
    for cat in ECKER_ANNOTATIONS:
        for ann in ECKER_ANNOTATIONS[cat]:
            for major_region, cell_class in ECKER_GROUPS:
                res.append(f"{ann}_{major_region}_{cell_class}.det.out.gz")
    return [op.join(ECKER_OUTPUT, item) for item in res]


rule render_ecker_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        list_ecker_yamet_outputs(),
    params:
        output_path=ECKER_OUTPUT,
        chr10_only=ECKER_CHR10_ONLY,
    threads:
        round(workflow.cores / 2)
    output:
        op.join(ECKER_BASE, "results", "ecker.html"),
    log:
        log=op.join("logs", "render_ecker.log"),
    script:
        "src/ecker.Rmd"
