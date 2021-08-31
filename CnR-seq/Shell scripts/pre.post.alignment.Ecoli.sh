
WORKDIR="/home/xchen2/SWISNF/Xiaolong/CnR/"
FASTQINPUT="/home/xchen2/SWISNF/CnR/FASTQ/"


####### Trim adaptor
TRIMMEDFQ="/home/xchen2/SWISNF/CnR/FASTQ/TRIMMED"
for i in `ls $FASTQINPUT/*R1*.gz ` ;do bsub -J trimedfq trim_galore --paired --fastqc --gzip $i ${i/R1/R2} ; done 

####### Map to DM6 genome to calculated spike-in 
cd $TRIMMEDFQ

for R1 in *.R1.fq.gz ;do 
sampleName=${i/.R1.fq.gz};
R2=${R1/R1/R2};
mm10=$(cat /home/xchen2/SWISNF/Xiaolong/chrs.mm10.lst |tr "\n" " "); 
echo "
bowtie2 -x /rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Drosophila_melanogaster/Ensembl/r97/bowtie2-index/v2.3.5.1/Drosophila_melanogaster.BDGP6.22.dna.toplevel -1 $R1 -2 $R2 -S $sampleName.Ecoli.sam ; samtools view -Sb -F 4 $sampleName.Ecoli.sam > $sampleName.Ecoli.bam 
samtools sort $sampleName.Ecoli.bam $sampleName.Ecoli.sorted
samtools index $sampleName.Ecoli.sorted
samtools view -q 1 -b $sampleName.Ecoli.sorted |samtools flagstat - > $sampleName.Ecoli.flagstat
"  >$sample.sh ; done 