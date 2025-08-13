"""
Handles CRC data

Single-cell multiomics sequencing and analyses of human colorectal cancer

Shuhui Bian https://orcid.org/0000-0002-9662-113X, Yu Hou https://orcid.org/0000-0001-7875-7087, Xin Zhou https://orcid.org/0000-0002-4048-4017, Xianlong Li https://orcid.org/0000-0001-7745-2695, Jun Yong https://orcid.org/0000-0002-3770-2108, Yicheng Wang https://orcid.org/0000-0002-3901-1482, Wendong Wang https://orcid.org/0000-0002-8631-3109, Jia Yan https://orcid.org/0000-0002-6203-0016, Boqiang Hu https://orcid.org/0000-0002-2045-3619, Hongshan Guo https://orcid.org/0000-0003-2799-2989, Jilian Wang https://orcid.org/0000-0002-3951-6931, Shuai Gao https://orcid.org/0000-0002-8208-9197, Yunuo Mao https://orcid.org/0000-0003-1428-270X, Ji Dong, Ping Zhu, Dianrong Xiu, Liying Yan https://orcid.org/0000-0001-9572-9440, Lu Wen, Jie Qiao https://orcid.org/0000-0003-2126-1376 , Fuchou Tang https://orcid.org/0000-0002-8625-7717 , and Wei Fu https://orcid.org/0000-0001-5248-7891

https://www.science.org/doi/10.1126/science.aao3791?url_ver=Z39.88-2003

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693

This is hg19

"""

from glob import glob

CRC_RAW = op.join("crc", "raw")            # bismark files
CRC_DATA = op.join("crc", "data")          # awk-ed to suit yamet
CRC_YAMET = op.join("output", "crc_yamet") # yamet outputs


rule download_crc_accessors:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        gsm=op.join("crc", "bisulfites_gsm.txt"),
    message:
        """
        CRC Accessors download
        """
    script:
        "src/get_crc_accessors.sh"

def list_raw_files_from_metadata():
    input_list_path = op.join("src", "sctrioseq_bismark_files.txt")
    with open(input_list_path, "r") as f:
        files = [line.strip() for line in f if line.strip()]
    return [op.join(CRC_RAW, file) for file in files]


rule download_crc_bismark_report:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        gsm=op.join("crc", "bisulfites_gsm.txt"),
    output:
        protected(list_raw_files_from_metadata())
    params:
        raw=CRC_RAW,
    message:
        """
        CRC download
        """
    script:
        "src/get_crc_meth_files.sh"


rule harmonize_single_crc_bismark_report:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(CRC_RAW, "download.flag"),
        list_raw_files_from_metadata()
    output:
        protected(op.join(CRC_DATA, "{file}")),
    threads:
        max(8, workflow.cores/8) # to reduce IO pressure
    params:
        raw=CRC_RAW,
    script:
        "src/parse_crc_meth_file.sh"

def get_raw_files(patient, stage):
    filepaths = glob(op.join(CRC_RAW, f"G*_{patient}_{stage}*.txt.gz"))
    return [op.basename(file) for file in filepaths]


rule run_yamet_cg:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        yamet=op.join("build", "yamet"),
        cells=lambda wildcards: expand(
            op.join(CRC_DATA, "{file}"),
            file=get_raw_files(wildcards.patient, wildcards.stage),
        ),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        intervals=op.join(HG19_BASE, "{subcat}.{cat}.bed.gz"),
    output:
        out=op.join(CRC_YAMET, "{subcat}.{cat}.{patient}.{stage}.out.gz"),
        det_out=op.join(CRC_YAMET, "{subcat}.{cat}.{patient}.{stage}.det.out.gz"),
        norm_det_out=op.join(CRC_YAMET, "{subcat}.{cat}.{patient}.{stage}.norm.out.gz"),
        meth_out=op.join(CRC_YAMET, "{subcat}.{cat}.{patient}.{stage}.meth.out.gz"),
    group:
        "yamet"
    threads: max(8, workflow.cores/8) # to reduce IO pressure
    params:
        base=CRC_YAMET,
    script:
        "src/yamet.sh"

SAMPLES = {
    "CRC01": ["NC", "PT", "LN", "ML", "MP"],
    "CRC02": ["NC", "PT", "ML", "PT"],
    "CRC04": ["NC", "PT"],
    "CRC10": ["NC", "PT", "LN"],
    "CRC11": ["NC", "PT", "LN"],
    "CRC13": ["NC", "PT", "LN"],
    "CRC15": ["NC", "PT", "LN", "ML", "MO"],
}


def crc_yamet_outputs():
    res = []
    for ann in ANN_MAP:
        for subann in ANN_MAP[ann]:
            for patient in SAMPLES:
                for sample in SAMPLES[patient]:
                    res.append(f"{subann}.{ann}.{patient}.{sample}.det.out.gz")
    return [op.join(CRC_YAMET, item) for item in res]


rule crc_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        crc_yamet_outputs(),
        annotation=op.join(HG19_BASE, "bookend_annotation.gz"),
    params:
        yamet=CRC_YAMET,
    threads: 16
    output:
        op.join("output", "crc.html"),
    script:
        "src/crc.Rmd"


rule crc_stats_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        crc_yamet_outputs(),
        annotation=op.join(HG19_BASE, "bookend_annotation.gz"),
    params:
        yamet=CRC_YAMET,
    threads: 16
    output:
        op.join("output", "crc_stats.html"),
    script:
        "src/crc_stats.Rmd"


rule crc_deletions_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        expand(
            op.join(CRC_YAMET, "bookended_250k.custom.{patsamp}.det.out.gz"),
            patsamp=[
                f"{patient}.{sample}"
                for patient, samples in SAMPLES.items()
                for sample in samples
            ],
        ),
    params:
        yamet=CRC_YAMET,
    threads: 16
    output:
        op.join("output", "crc_deletions.html"),
    script:
        "src/crc_deletions.Rmd"
