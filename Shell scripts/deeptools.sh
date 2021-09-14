# put all required information into a meta file multiBigwigSummary.meta.file.txt
# Assay_Name	peakfile_Name	bigwig_file_separated_by_blank
# Asyn_SOX2_ChIP	SOX2_ChIP.3types.peaks.bed	Asyn_SOX2_ChIP.rep1.mm10.bw Asyn_SOX2_ChIP.rep2.mm10.bw Mit_SOX2_ChIP.rep1.mm10.bw
while read -r sample peak bw ;
do
multiBigwigSummary BED-file -b $bw -o multiBWSummary/$sample.peaks.npz --BED $peak --outRawCounts multiBWSummary/$sample.peaks.fpkm --smartLabels -p 1   ;
done < multiBigwigSummary.meta.file.txt



# put all required information into a meta file computeMatrix.meta.file.txt, preferably ,Asyn_only , shared and Mit_only peaks separately in the meta file. 
# Assay_Name	peakfile_Name	bigwig_file_separated_by_blank
# Asyn_SOX2_ChIP	SOX2_ChIP.3types.peaks.bed	Asyn_SOX2_ChIP.rep1.mm10.bw Asyn_SOX2_ChIP.rep2.mm10.bw 
# Asyn_SOX2_ChIP	SOX2_ChIP.Asyn_only.peaks.bed	Asyn_SOX2_ChIP.rep1.mm10.bw Asyn_SOX2_ChIP.rep2.mm10.bw 
# Asyn_SOX2_ChIP	SOX2_ChIP.overlap.peaks.bed	Asyn_SOX2_ChIP.rep1.mm10.bw Asyn_SOX2_ChIP.rep2.mm10.bw 
# Asyn_SOX2_ChIP	SOX2_ChIP.Mit_only.peaks.bed	Asyn_SOX2_ChIP.rep1.mm10.bw Asyn_SOX2_ChIP.rep2.mm10.bw 

while read -r sample peak bw ;
do
computematrix.out computeMatrix reference-point -R $peak -S $bw -o computeMatrix/$sample.$peak.plot --outFileNameMatrix computeMatrix/$sample.$peak.matrix --outFileSortedRegions computeMatrix/$sample.$peak.sorted.bed --sortRegions descend --sortUsingSamples 1 --smartLabels --missingDataAsZero -p 4 -a 3000 -b 3000 -bs 50  --referencePoint center;
done < computeMatrix.meta.file.txt

# plot heatmap 
for i in computeMatrix/*.plot
plotHeatmap -m $i --dpi 300 --heatmapWidth 7 --whatToShow "heatmap and colorbar" --colorMap Reds --outFileName $i.pdf

############## 
##multiBigwigSummary for ESRRB , EZH2 and SOX2 ChIP in SMARCE1-MD cells 
##############

cd bookmarkers.SMARCE1.MD

multiBigwigSummary BED-file -b -o multibw/ESRRB.npz --BED ESRRB.3types.peaks.bed --outRawCounts multibw/ESRRB.fpkm --smartLabels -p1
multiBigwigSummary BED-file -b -o multibw/EZH2.npz --BED EZH2.3types.peaks.bed --outRawCounts multibw/EZH2.fpkm --smartLabels -p1
multiBigwigSummary BED-file -b -o multibw/SOX2.npz --BED SOX2.3types.peaks.bed --outRawCounts multibw/SOX2.fpkm --smartLabels -p1