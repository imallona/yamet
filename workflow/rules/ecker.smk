"""
Ecker's related rules

DNA methylation atlas of the mouse brain at single-cell resolution

Hanqing Liu, Jingtian Zhou, Wei Tian, Chongyuan Luo, Anna Bartlett, Andrew Aldridge, Jacinta Lucero, Julia K. Osteen, Joseph R. Nery, Huaming Chen, Angeline Rivkin, Rosa G. Castanon, Ben Clock, Yang Eric Li, Xiaomeng Hou, Olivier B. Poirion, Sebastian Preissl, Antonio Pinto-Duarte, Carolyn O’Connor, Lara Boggeman, Conor Fitzpatrick, Michael Nunn, Eran A. Mukamel, Zhuzhu Zhang, …Joseph R. Ecker

https://www.nature.com/articles/s41586-020-03182-8

this is mm10
"""

ECKER_BASE = "ecker"


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
