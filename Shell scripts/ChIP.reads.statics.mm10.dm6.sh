WORKDIR="/home/xchen2/SWISNF/Xiaolong/ChIP/"

for i in *.dm6.bam ;do 
sampleName=${i/.dm6.bam}; 
allmm10=`ls $sampleName".bam.flagstat"` ;  ######### get total reads 
mm10q1=`ls $sampleName".bam.q1.flagstat"`; ########## get mapped reads to mm10  with mapq > 1 
EColiq1=`ls $sampleName".dm6.bam.q1.flagstat"`;  ############get mapped reads to dm6 with mapq > 1

totalreads=`head -n 1 $allmm10 |awk '{print $1}' ` ; 
mappedreads2mm10=`grep "mapped" $mm10q1|head -n 1 |awk '{print $1}'  ` ; 
mappedreasd2dm6=`awk 'NR==5{print $1}' $EColiq1 `; 

echo $sampleName $totalreads $mappedreads2mm10 $mappedreasd2dm6 | awk 'BEGIN{FS=" ";OFS="\t"}{print $0}' >> ChIP.mappedreads.statics.txt ;done
awk 'BEGIN{print "sampleName\ttotalReads\treads2mm10\treads2dm6"}{print $0} ' ChIP.mappedreads.statics.txt >tmp ; mv tmp ChIP.mappedreads.statics.txt 