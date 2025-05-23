"""
Handles CRC data

Single-cell multiomics sequencing and analyses of human colorectal cancer

Shuhui Bian https://orcid.org/0000-0002-9662-113X, Yu Hou https://orcid.org/0000-0001-7875-7087, Xin Zhou https://orcid.org/0000-0002-4048-4017, Xianlong Li https://orcid.org/0000-0001-7745-2695, Jun Yong https://orcid.org/0000-0002-3770-2108, Yicheng Wang https://orcid.org/0000-0002-3901-1482, Wendong Wang https://orcid.org/0000-0002-8631-3109, Jia Yan https://orcid.org/0000-0002-6203-0016, Boqiang Hu https://orcid.org/0000-0002-2045-3619, Hongshan Guo https://orcid.org/0000-0003-2799-2989, Jilian Wang https://orcid.org/0000-0002-3951-6931, Shuai Gao https://orcid.org/0000-0002-8208-9197, Yunuo Mao https://orcid.org/0000-0003-1428-270X, Ji Dong, Ping Zhu, Dianrong Xiu, Liying Yan https://orcid.org/0000-0001-9572-9440, Lu Wen, Jie Qiao https://orcid.org/0000-0003-2126-1376 , Fuchou Tang https://orcid.org/0000-0002-8625-7717 , and Wei Fu https://orcid.org/0000-0001-5248-7891

https://www.science.org/doi/10.1126/science.aao3791?url_ver=Z39.88-2003

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693

This is hg19

"""

from glob import glob

CRC_RAW = op.join("crc", "raw")
CRC_DATA = op.join("crc", "data")
CRC_YAMET = op.join("crc", "yamet")


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


rule download_crc:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        gsm=op.join("crc", "bisulfites_gsm.txt"),
    output:
        touch(op.join("crc", "download.flag")),
    params:
        raw=CRC_RAW,
    message:
        """
        CRC download
        """
    script:
        "src/get_crc_meth_files.sh"


rule parse_single_crc:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join("crc", "download.flag"),
    output:
        op.join(CRC_DATA, "{file}"),
    params:
        raw=CRC_RAW,
    script:
        "src/parse_crc_meth_file.sh"


def get_raw_files(patient, stage):
    filepaths = glob(op.join(CRC_RAW, f"G*_{patient}_{stage}*.txt.gz"))
    return [op.basename(file) for file in filepaths]


rule yamet_crc_cg:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        yamet=op.join("build", "yamet"),
        cells=lambda wildcards: expand(
            op.join(CRC_DATA, "{file}"),
            file=get_raw_files(wildcards.patient, wildcards.stage),
        ),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        intervals=op.join(HG19_BASE, "{subcat}.{cat}.bed"),
    output:
        out=op.join(CRC_YAMET, "{subcat}.{cat}.{patient}.{stage}.out"),
        det_out=op.join(CRC_YAMET, "{subcat}.{cat}.{patient}.{stage}.det.out"),
    group:
        "yamet"
    threads: 16
    params:
        base=CRC_YAMET,
    script:
        "src/yamet.sh"


CAT_MAP = {
    "pmd": ["pmds", "hmds"],
    "hmm": [
        "0_Enhancer",
        "2_Enhancer",
        "11_Promoter",
        "12_Promoter",
        "1_Transcribed",
        "4_Transcribed",
        "5_RegPermissive",
        "7_RegPermissive",
        "6_LowConfidence",
        "3_Quiescent",
        "8_Quiescent",
        "10_Quiescent",
        "9_ConstitutiveHet",
        "13_ConstitutiveHet",
    ],
    "chip": ["H3K27me3", "H3K9me3", "H3K4me3"],
    "lad": ["laminb1"],
}
STAGES = ["NC", "PT"]
PATIENTS = ["CRC01", "CRC02", "CRC04", "CRC10", "CRC11", "CRC13", "CRC15"]


rule crc_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        expand(
            op.join(CRC_YAMET, "{jnt}.{patient}.{stage}.out"),
            jnt=[f"{subcat}.{cat}" for cat in CAT_MAP for subcat in CAT_MAP[cat]],
            stage=STAGES,
            patient=PATIENTS,
        ),
    params:
        yamet=CRC_YAMET,
    output:
        op.join("crc", "crc.html"),
    script:
        "src/crc.Rmd"
