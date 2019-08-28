#!/bin/bash
##
## Downloads ega data using the ega demo client
## 
## Izaskun Mallona
## 22nd July 2019

echo 'Run in taupo'

EGACLIENT=~/soft/egaclient/EGA_download_client_2.2.2/EgaDemoClient.jar

cd ~/mnt/nfs/ega_secure

java -jar $EGACLIENT 

# izaskun.mallona@dmmd.uzh.ch
# ega pass, documented elsewhere

testbandwidth 8

datasets

size dataset EGAD00001004074
## 166.3 GB

files dataset EGAD00001004074

request dataset EGAD00001004074 first_attempt request_EGAD00001004074

## with eight cores
download request_EGAD00001004074 8

## be careful, you see the password if `ps`-ing
