# put all Histone modification bam files in a folder /chromHMM.input.bam
# config file : chromHMM.HM.txt
# BinarizeBam 

java -Xmx20g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar BinarizeBam -paired -strictthresh /home/xchen2/SWISNF/Xiaolong/mm10.genome /home/xchen2/SWISNF/Xiaolong/summary/chromHMM.input.bam /home/xchen2/SWISNF/Xiaolong/summary/chromHMM.HM.txt /home/xchen2/SWISNF/Xiaolong/summary/chromHMM.HM.strict.reps

# LearnModel try from 8 clusters to 25 clusters.

for i in {8..25}; do 
bsub -R'rusage[mem=15000]' -o chromhmm.learnmodel.stri.$i.out -e chromhmm.learnmodel.stri.$i.err \
java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar LearnModel chromHMM.HM.strict.reps "chromHMM.HM.strict.reps.learnmodel."$i $i mm10 ;
done 

# reorder the 15 states of the model

java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar reorder -f /home/xchen2/SWISNF/Xiaolong/summary/chromhmm.reorder.hmnames.txt -o /home/xchen2/SWISNF/Xiaolong/summary/chromhmm.reorder.states.txt /home/xchen2/SWISNF/Xiaolong/summary/chromHMM.HM.strict.reps.learnmodel.15/model_15.txt chromhmm.reorder

# make segmentation for the reordered model

java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar MakeSegmentation  chromhmm.reorder/model_15.txt ~/SWISNF/Xiaolong/summary/chromHMM.HM.strict.reps chromhmm.reorder

# overlap with annotation in chromHMM 

java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar OverlapEnrichment chromhmm.reorder/Asyn_15_segments.bed /home/xchen2/software/ChromHMM/COORDS/mm10/ chromhmm.reorder/Asyn_overlap

java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar OverlapEnrichment chromhmm.reorder/Mit_15_segments.bed /home/xchen2/software/ChromHMM/COORDS/mm10/ chromhmm.reorder/Mit_overlap

# overlap with SWISNF peaks predicted by MACS2 in Asyn and Mit cells
# put peaks separately in two folders 

java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar OverlapEnrichment chromhmm.reorder/Asyn_15_segments.bed chromhmm.reorder/cnr/asyn.peaks/ chromhmm.overlap/Asyn.CnR.overlap
java -Xmx8g -jar /home/xchen2/software/ChromHMM/ChromHMM.jar OverlapEnrichment chromhmm.reorder/Mit_15_segments.bed chromhmm.reorder/cnr/mit.peaks/ chromhmm.overlap/Mit.CnR.overlap