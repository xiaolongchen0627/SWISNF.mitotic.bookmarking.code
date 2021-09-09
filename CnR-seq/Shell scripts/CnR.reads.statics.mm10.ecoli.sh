WORKDIR="/home/xchen2/SWISNF/Xiaolong/CnR/"

for i in *.Ecoli.bam ;do 
sampleName=${i/.Ecoli.bam}; 
allmm10=`ls $sampleName".bam.flagstat"` ;  ######### get total reads 
mm10q1=`ls $sampleName".bam.q1.flagstat"`; ########## get mapped reads to mm10  with mapq > 1 
EColiq1=`ls $sampleName".Ecoli.bam.q1.flagstat"`;  ############get mapped reads to Ecoli with mapq > 1

totalreads=`head -n 1 $allmm10 |awk '{print $1}' ` ; 
mappedreads2mm10=`grep "mapped" $mm10q1|head -n 1 |awk '{print $1}'  ` ; 
mappedreasd2Ecoli=`awk 'NR==5{print $1}' $EColiq1 `; 

echo $sampleName $totalreads $mappedreads2mm10 $mappedreasd2Ecoli | awk 'BEGIN{FS=" ";OFS="\t"}{print $0}' >> CnR.mappedreads.statics.txt ;done
awk 'BEGIN{print "sampleName\ttotalReads\treads2mm10\treads2Ecoli"}{print $0} ' CnR.mappedreads.statics.txt >tmp ; mv tmp CnR.mappedreads.statics.txt 