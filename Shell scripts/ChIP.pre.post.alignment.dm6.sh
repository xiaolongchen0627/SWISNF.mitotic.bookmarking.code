
WORKDIR="/home/xchen2/SWISNF/Xiaolong/ChIP/"
FASTQINPUT="/home/xchen2/SWISNF/Xiaolong/ChIP/FASTQ/"


####### Trim adaptor
TRIMMEDFQ="/home/xchen2/SWISNF/Xiaolong/ChIP/FASTQ/TRIMMED"
for i in `ls $FASTQINPUT/*R1*.gz ` ;do bsub -J trimedfq trim_galore --paired --fastqc --gzip $i ${i/R1/R2} ; done 

####### Map to DM6 genome to calculated spike-in 
cd $TRIMMEDFQ

for R1 in *.R1.fq.gz ;do 
sampleName=${i/.R1.fq.gz};
R2=${R1/R1/R2};
echo "
bowtie2 -x /rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Drosophila_melanogaster/Ensembl/r97/bowtie2-index/v2.3.5.1/Drosophila_melanogaster.BDGP6.22.dna.toplevel -1 $R1 -2 $R2 -S $sampleName.dm6.sam ; samtools view -Sb -F 4 $sampleName.dm6.sam > $sampleName.dm6.bam 
samtools sort $sampleName.dm6.bam $sampleName.dm6.sorted
samtools index $sampleName.dm6.sorted
samtools view -q 1 -b $sampleName.dm6.sorted |samtools flagstat - > $sampleName.dm6.flagstat
"  >$sample.sh ; done 