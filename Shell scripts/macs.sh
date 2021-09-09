treat=$1
input=$2

macs2 callpeak -t $treat"_mm10.rmdupq1.bam" -c $input"_mm10.rmdupq1.bam" -g mm -n $treat --outdir macs2.peaks