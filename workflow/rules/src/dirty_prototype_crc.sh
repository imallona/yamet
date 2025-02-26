#!/bin/bash

##
## Processes the CRC data quick and dirty
##

# 1. automate the process substitutions for all normals
# 2. run yamet per feature
# 3. sort features per avg methylation
# 4. repeat for primary tumors
# 5. compare whether ranks are inverted
# 6. evaluate whether there is a within-cell-type or across-cell-types change in meth variability

## Started 25th Feb 2025
##
## Izaskun Mallona
## GPLv3

DATAPATH=$HOME/src/yamet/crc_data
WD=$HOME/tmp/yamet_crc

normals=$(find $DATAPATH -name "*CRC01_NC*gz" | sort | paste -sd " ")
primaries=$(find $DATAPATH -name "*CRC01_PT*gz" |  sort| paste -sd " ")

# download pmds

mkdir -p $WD
cd $WD

wget https://zhouserver.research.chop.edu/GenomeAnnotation/hg19/PMD_coordinates_hg19.bed.gz

gunzip PMD*gz

# download / installed yamet

# wget https://github.com/imallona/yamet/releases/download/v1.1.0-rc.2/yamet-Linux-x86_64-1.1.0-rc.2.deb
# dpkg-deb -e yamet-Linux-x86_64-1.1.0-rc.2.deb yamet

# well, compiling it, caution cmake 3.31.6, from the same branch

cd method
bash build.sh
export PATH=/home/imallona/src/yamet/method/build:$PATH


# yametize cell meth files

mkdir -p nc pt

for normal in $normals
do
    bn=$(basename $normal .txt.gz)
    zcat $normal |\
        grep CpG | \
        awk '{OFS="\t";} {
              if($4=="+")
                print $1, $2,   $2+1, "", "",$4,$6,$5
              else
                print $1, $2-1, $2,   "", "",$4,$6,$5
              }' | \
        bedtools merge -c 7,8 -o sum  | \
        awk '{OFS=FS="\t"; print $1,$2,$4,$5,$4/$5}' | \
        sort -k1,1 -k2,2n -k3,3n | gzip -c > nc/"$bn".yametized.gz
done

for primary in $primaries
do
    bn=$(basename $primary .txt.gz)
    zcat $primary |\
        grep CpG | \
        awk '{OFS="\t";} {
              if($4=="+")
                print $1, $2,   $2+1, "", "",$4,$6,$5
              else
                print $1, $2-1, $2, "",   "",$4,$6,$5
              }' | \
        bedtools merge -c 7,8 -o sum  | \
        awk '{OFS=FS="\t"; print $1,$2,$4,$5,$4/$5}' | \
        sort -k1,1 -k2,2n -k3,3n | gzip -c > pt/"$bn".yametized.gz
done

## get the reference CpG file. Again we fake it getting one of the large normal files and using
##   it as a reference

zcat nc/GSM2697701_scTrioSeq2Met_CRC01_NC_529.singleC.yametized.gz |\
    awk '{OFS=FS="\t"; print $1,$2,$2+1}'|  gzip -c > custom.ref.gz


## generate genomic bins file - 10k each window, book ended

mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
      'SELECT chrom,size FROM hg19.chromInfo' > hg19.sizes

bedtools makewindows -g hg19.sizes -w 10000 | sort -k1,1 -k2,2n -k3,3n > hg19_windows.bed

# run yamet on windows

# instead of using process substitutions, automate yamet-zation
# these should be temporary files within the snmk

yamet --cell nc/*NC*yametized.gz \
      --reference custom.ref.gz \
      --intervals hg19_windows.bed \
      --out normals.yamet.out \
      --cores 30 \
      --det-out normals.yamet.det.out \
      --print-sampens F 
      

yamet --cell pt/*PT*yametized.gz \
      --reference custom.ref.gz \
      --intervals hg19_windows.bed \
      --out primaries.yamet.out \
      --cores 30 \
      --det-out primaries.yamet.det.out \
      --print-sampens F 

# run yamet on hmm segmentations

wget https://www.encodeproject.org/files/ENCFF526MRN/@@download/ENCFF526MRN.bed.gz
gunzip ENCFF526MRN.bed.gz

sort -k1,1 -k2,2n -k3,3n ENCFF526MRN.bed > ENCFF526MRN.bed.sorted
mv ENCFF526MRN.bed.sorted ENCFF526MRN.bed

## now look for associations with features; what to do with the -1s?

## managed to make this segfault once with 30 cores (?) but works with 50
for color in 0_Enhancer 10_Quiescent 11_Promoter 12_Promoter 13_ConstitutiveHet 1_Transcribed \
             2_Enhancer 3_Quiescent 4_Transcribed 5_RegPermissive 6_LowConfidence \
             7_RegPermissive 8_Quiescent 9_ConstitutiveHet
do
    echo $color
    grep -w $color ENCFF526MRN.bed > "$color".bed

    echo nc
    yamet --cell nc/*NC*yametized.gz \
      --reference custom.ref.gz \
      --intervals "$color.bed" \
      --out "$color"_normals.yamet.out \
      --cores 50 \
      --det-out "$color"_normals.yamet.det.out \
      --print-sampens F 
      
    echo pt
    yamet --cell pt/*PT*yametized.gz \
      --reference custom.ref.gz \
      --intervals "$color".bed \
      --out "$color"_primaries.yamet.out \
      --cores 50 \
      --det-out "$color"_primaries.yamet.det.out \
      --print-sampens F    
done


echo 'quick and dirty, perhaps deeptools this? or intervene it?'
# conda install -c bioconda -c conda-forge -c default \
#          bioconda::intervene cmake==3.31.6 bioconda::pybedtools==0.11.0
## pff, not solvable, let's skip this and do manual processing in R

# assuming the det.out are bed-like...
# intervene pairwise -i *det.out --filenames --compute frac --htype color

echo 'knit the dirty_prototype_crc.Rmd'

echo 'unimplemented thoughts'

# rather do that by 10kbp-long windows, get all genomic annotations, and evaluate enrichment via multiBigwigSummary
      # https://deeptools.readthedocs.io/en/latest/content/feature/plotFingerprint_QC_metrics.html
      # https://www.nature.com/articles/s41467-021-21707-1 for chromhmms for cancer

      # https://github.com/jernst98/ChromHMM/tree/master/COORDS/hg19 for lamina -lads
      # chromhmm for hg19 https://www.encodeproject.org/annotations/ENCSR814YSQ/
:<<EOF

easiest would be to chop an hmm segmentation into its components, so the whole genome is there - for H1 and the roadmap healthy colon https://pmc.ncbi.nlm.nih.gov/articles/PMC4530010/
also add one for sequence conservation!

correlations heatmap
plotfingerprint
and plotting via scaling to see 5' or 3' enrich of entropy

don't forget superenhancers! as in https://www.nature.com/articles/s41467-021-21707-1#data-availability ?
EOF


:<<EOF
what about adding to yamet a param to specify the col with the total, the meth, the unmeth, and btw provide total or not?
EOF
