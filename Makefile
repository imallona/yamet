## cg_shadows
##
## Izaskun Mallona
## 12th nov 2018
## does not work


R := ~/soft/R/R-3.5.1/bin/R --no-restore --no-save
WD := ~/cg_shadows
METILENE := ~/soft/metilene/metilene_v0.2-7/metilene
METILENE_INPUT := ~/soft/metilene/metilene_v0.2-7/metilene_input.pl
METILENE_OUTPUT := ~/soft/metilene/metilene_v0.2-7/metilene_output.pl
BEDTOOLS := ~/soft/bedtools/bin/bedtools
METHTUPLE_ACTIVATE := ~/virtenvs/methtuple/bin/activate


.DEFAULT: all
.PHONY: all clean

venv:
	. $(METHTUPLE_ACTIVATE)

bam := data/sorjuela.bam
methtuple: data/sorjuela_CG.2.tsv
data/sorjuela_CG.2.tsv:
	venv ; methtuple -m 2 --methylation-type CG $(bam)
	deactivate

all: methtuple
