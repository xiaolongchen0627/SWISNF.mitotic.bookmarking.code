# cat narrowPeak from MACS2 for each replicate together in one single file eg Asyn.SOX2.CnR.narrowPeaks.cat.bed

for i in *narrowPeaks.cat.bed ;do echo $i ; awk '$8>=9&&$7>5' $i sortBed | bedops -m - > ${i/.narrowPeaks.cat.bed/.merged1bp.peaks.bed}  ; done
