#### ctcf peaks were predicted using CTCF chip-seq data in Asyn ES cells. put all SMARCE1 peaks with 3 types of annotation ,Asyn_only , Mit_only and Shared together .

BETA minus -p SMARCE1.CnR.peaks.bed -g mm10 -n minus.ctcf.100000 -o peaks.minus --pn 150000 -d 100000 --bl --bf CTCF.6C.bed