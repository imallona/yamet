"""
Handles CRC data

Single-cell multiomics sequencing and analyses of human colorectal cancer

Shuhui Bian https://orcid.org/0000-0002-9662-113X, Yu Hou https://orcid.org/0000-0001-7875-7087, Xin Zhou https://orcid.org/0000-0002-4048-4017, Xianlong Li https://orcid.org/0000-0001-7745-2695, Jun Yong https://orcid.org/0000-0002-3770-2108, Yicheng Wang https://orcid.org/0000-0002-3901-1482, Wendong Wang https://orcid.org/0000-0002-8631-3109, Jia Yan https://orcid.org/0000-0002-6203-0016, Boqiang Hu https://orcid.org/0000-0002-2045-3619, Hongshan Guo https://orcid.org/0000-0003-2799-2989, Jilian Wang https://orcid.org/0000-0002-3951-6931, Shuai Gao https://orcid.org/0000-0002-8208-9197, Yunuo Mao https://orcid.org/0000-0003-1428-270X, Ji Dong, Ping Zhu, Dianrong Xiu, Liying Yan https://orcid.org/0000-0001-9572-9440, Lu Wen, Jie Qiao https://orcid.org/0000-0003-2126-1376 , Fuchou Tang https://orcid.org/0000-0002-8625-7717 , and Wei Fu https://orcid.org/0000-0001-5248-7891

https://www.science.org/doi/10.1126/science.aao3791?url_ver=Z39.88-2003

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693

This is hg19

Izaskun Mallona
GPLv3
"""

import os.path as op
from glob import glob

CRC = op.join("data", "crc")
CRC_RAW = op.join(CRC, "raw")  ## raw, as in GEO and from GEO
CRC_HARMONIZED = op.join(CRC, "harmonized")  ## ingestable by yamet
CRC_OUTPUT = op.join(CRC, "output")  ## output for features (genes, promoters etc)
CRC_WINDOWS_OUTPUT = op.join(CRC, "windows_output")  ## output for tiles/genomic windows


## patient -> biopsy sites mapping
SAMPLES = {
    "CRC01": ["NC", "PT", "LN", "ML", "MP"],
    "CRC02": ["NC", "PT", "ML", "PT"],
    "CRC04": ["NC", "PT"],
    "CRC10": ["NC", "PT", "LN"],
    "CRC11": ["NC", "PT", "LN"],
    "CRC13": ["NC", "PT", "LN"],
    "CRC15": ["NC", "PT", "LN", "ML", "MO"],
}  # What is MO, a typo?


def list_raw_files_from_metadata():
    """
    Reports a hardcoded list of bismark meth reports, to make sure download is complete
    """

    input_list_path = op.join("src", "sctrioseq_bismark_files.txt")
    with open(input_list_path, "r") as f:
        files = [line.strip() for line in f if line.strip()]
    return [op.join(CRC_RAW, file) for file in files]


rule download_crc_accessors:
    conda:
        op.join("..", "envs", "processing.yml")
    output:
        gsm=op.join(CRC, "bisulfites_gsm.txt"),
    message:
        """
        CRC Accessors download
        """
    params:
        raw=CRC_RAW,
    shell:
        """
        mkdir -p {params.raw}
        
        esearch -db sra -query PRJNA382695 |
          efetch -format runinfo |
          cut -f11,13,14,15,30 -d"," |
          grep Bis |
          cut -f5 -d"," > {output.gsm}
        """


checkpoint download_crc_bismarks:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        gsm=op.join(CRC, "bisulfites_gsm.txt"),
    output:
        urls=temp(op.join(CRC_RAW, "tmp.urls")),
        download_flag=op.join(CRC, "download.flag"),  #,
        # files = list_raw_files_from_metadata() ## commented out so snmk does not remove all these files before re-running the rule, so the finding/skipping algo makes sense
    log:
        op.join("logs", "crc_download.log"),
    params:
        raw=CRC_RAW,
    threads: 1
    message:
        """
        CRC download
        """
    shell:  ## this is stupid, because it wipes out the folder before exec
        """
        exec &> {log}
        
        set -euo pipefail
        mkdir -p "{params.raw}"
        
        while read -r gsm
          do
            echo $gsm
            short="$(echo $gsm | cut -c1-7)"
            url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"$short"nnn/"$gsm"/suppl/
        
            # wget -q -e robots=off -r -k -A gz -nd -P {params.raw} $url
            matches=$(find "{params.raw}" \
                   -maxdepth 1 \
                   -type f \
                   -name "${{gsm}}*" \
                   -print)
            echo "$matches"
 
            if [[ -n "$matches" ]]; then
               echo "Found existing report for $gsm:"
               printf "  %s Skipping download\\n\\n" "${{matches[@]}}"
            else
               printf "Downloading $gsm from $url\\n\\n"

               wget \
                  --quiet \
                  --execute=robots=off \
                  --recursive \
                  --convert-links \
                  --accept=gz \
                  --no-directories \
                  --directory-prefix={params.raw} \
                  --timestamping \
                  "$url"
            
                echo "$url" >> {output.urls}
           fi
          
        done < {input.gsm}
        
        cp {output.urls} {output.download_flag}
        """


rule harmonize_cell_report_for_yamet:
    conda:
        op.join("..", "envs", "processing.yml")
    input:
        op.join(CRC, "download.flag"),
        # list_raw_files_from_metadata(),
        op.join(CRC_RAW, "{file}"),
        # op.join(CRC_RAW, "{gsm}_{patient}_{location}_{cellid}.singleC.txt.gz")
    output:
        op.join(CRC_HARMONIZED, "{file}"),  ## these should tmp
    params:
        raw=CRC_RAW,
        harmonized=CRC_HARMONIZED,
    threads: max(4, workflow.cores / 4)
    log:
        op.join("logs", "{file}_harmonization.log"),
    shell:
        """
        exec &> {log}

        zcat {params.raw}/{wildcards.file} |
          grep "CpG$" |
        awk '
            {{OFS="\\t";}} {{
                if ($4 == "-") $2 = $2 - 2;
                else if ($4 == "+") $2 = $2 - 1;
                print $1, $2, $2+1, "", "", $4, $6, $5;
            }}
        ' |
        bedtools merge -c 7,8 -o sum |
        awk '
             {{OFS=FS="\\t";
                bin = ($4/$5 > 0.1) ? 1 : 0;
                print $1,$2,$4,$5,bin}}
        ' |
        sort -k1,1 -k2,2n |
        gzip -c > {output}
        """


## location is whether is a NC, PT, LN and so on
## patient is whether CRC01, CRC01 and so on
def get_harmonized_files(patient, location):
    ckpt = checkpoints.download_crc_bismarks.get()

    filepaths = glob(op.join(CRC_RAW, f"G*_{patient}_{location}*.txt.gz"))
    if len(filepaths) == 0:
        raise Exception("No cytosine reports matching the specs.")
    return [op.join(CRC_HARMONIZED, op.basename(file)) for file in filepaths]


rule run_yamet_on_separate_features:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        yamet=op.join("build", "yamet"),
        cells=lambda wildcards: expand(
            "{file}",
            file=get_harmonized_files(wildcards.patient, wildcards.location),
        ),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        intervals=op.join(HG19_BASE, "{subcat}.{cat}.bed.gz"),
    output:
        out=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.out.gz"),
        det_out=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.det.out.gz"),
        norm_det_out=op.join(
            CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.norm.out.gz"
        ),
        meth_out=op.join(CRC_OUTPUT, "{subcat}_{cat}_{patient}_{location}.meth.out.gz"),
    log:
        op.join("logs", "yamet_{subcat}_{cat}_{patient}_{location}.log"),
    group:
        "yamet"
    params:
        base=CRC_OUTPUT,
    threads: max(8, workflow.cores / 8)
    script:
        "src/yamet.sh"


def list_relevant_yamet_outputs():
    res = []
    for annot in ANNOTATIONS:
        for subannot in ANNOTATIONS[annot]:
            for patient in SAMPLES.keys():
                for sampling in SAMPLES[patient]:
                    res.append(f"{subannot}_{annot}_{patient}_{sampling}.det.out.gz")
    return [op.join(CRC_OUTPUT, item) for item in res]


rule render_crc_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        list_relevant_yamet_outputs(),
    params:
        output_path=CRC_OUTPUT,
    output:
        op.join("results", "crc.html"),
    log:
        log=op.join("logs", "render_crc.log"),
    threads: 8
    script:
        "src/crc.Rmd"


#####################################################################################################

## windows-based stuff ##############################################################################

#####################################################################################################

""" the idea is to get a common set of features, genomic tiles, to run
    some data mining techniques afterwards
    largey las in prototyping/crc_prototype.sh

    issue is, we need to permute backgrounds... and/or correct analytically
    favouring the analytical aproach now, skip permuting
"""


rule run_yamet_on_windows:
    conda:
        op.join("..", "envs", "yamet.yml")
    input:
        yamet=op.join("build", "yamet"),
        cells=lambda wildcards: expand(
            "{file}",
            file=get_harmonized_files(wildcards.patient, wildcards.location),
        ),
        ref=op.join(HG19_BASE, "ref.CG.gz"),
        intervals=op.join(HG19_BASE, "windows_{win_size}_nt.bed.gz"),
    output:
        out=op.join(CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.out.gz"),
        det_out=op.join(
            CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.det.out.gz"
        ),
        norm_det_out=op.join(
            CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.norm.out.gz"
        ),
        meth_out=op.join(
            CRC_WINDOWS_OUTPUT, "{win_size}_{patient}_{location}.meth.out.gz"
        ),
    log:
        op.join("logs", "yamet_{win_size}_{patient}_{location}.log"),
    group:
        "yamet"
    params:
        base=CRC_WINDOWS_OUTPUT,
    threads: max(8, workflow.cores / 8)
    script:
        "src/yamet.sh"


def list_relevant_yamet_windows_outputs():
    res = []
    for patient in SAMPLES.keys():
        for location in SAMPLES[patient]:
            res.append(f"{{win_size}}_{patient}_{location}.det.out.gz")
    return [op.join(CRC_WINDOWS_OUTPUT, item) for item in res]


rule render_crc_windows_report:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        yamet_dets=list_relevant_yamet_windows_outputs(),
        annotations=op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz"),
    params:
        output_path=CRC_WINDOWS_OUTPUT,
    output:
        op.join("results", "crc_windows_{win_size}_nt.html"),
    log:
        log=op.join("logs", "render_crc_windows_{win_size}.log"),
    threads: 8
    script:
        "src/crc_windows.Rmd"


rule crc_stats_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        list_relevant_yamet_windows_outputs(),
        annotation=op.join(HG19_BASE, "windows_{win_size}_nt_annotation.gz"),
    output:
        op.join("results", "crc_stats_{win_size}.html"),
    params:
        output_path=CRC_WINDOWS_OUTPUT,
    log:
        log=op.join("logs", "render_crc_stats_{win_size}.log"),
    threads: 16
    script:
        "src/crc_stats.Rmd"


rule crc_deletions_doc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        yamet_dets=list_relevant_yamet_windows_outputs(),
    output:
        op.join("results", "crc_deletions_{win_size}.html"),
    params:
        output_path=CRC_WINDOWS_OUTPUT,
    log:
        log=op.join("logs", "render_crc_deletions_{win_size}.log"),
    threads: 16
    script:
        "src/crc_deletions.Rmd"
