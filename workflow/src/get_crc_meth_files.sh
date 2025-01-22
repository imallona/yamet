#!/bin/bash
##
## Retrieves only meth cytosine reports from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693
## Single-cell Multi-omics Sequencing and Analyses of Human Colorectal Cancer
##
## Started 20thDec 2024


esearch -db gds -query GSE121708 | efetch -format docsum | xtract -pattern DocumentSummary -element FTPLink > urls.txt

esearch -db sra -query PRJNA382695 | efetch -format runinfo > runinfo.txt

cut -f11,13,14,15,30 -d"," runinfo.txt | grep Bis > bisulfites_accessions.txt

cut -f5 -d"," bisulfites_accessions.txt > bisulfites_gsm.txt

# grep -f bisulfites_gsm.txt urls.txt > bisulfites_urls.txt

# esearch -db gds -query GSM2575588 | efetch -format docsum | xtract -pattern DocumentSummary -element FTPLink | grep samples

while read -r gsm
do
    short="$(echo $gsm | cut -c1-7)"
    url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"$short"nnn/"$gsm"/suppl/
    echo "Downloading $gsm"
    wget -e robots=off -r -k -A gz $url
    echo "Done"
done < bisulfites_gsm.txt

touch 'crc_done.flag'
