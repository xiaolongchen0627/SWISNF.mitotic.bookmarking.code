treat=$1

samtools view -F 1804 -b  $treat"_mm10.rmdupq1.bam" |bamsort SO=queryname > $treat"_mm10.rmdupq1.bedpe.bam"
bedtools bamtobed -bedpe -i CnR/bams/$treat"_mm10.rmdupq1.bedpe.bam" > CnR/bams/$treat"_mm10.rmdupq1.bedpe"
nread_2k=$(cat CnR/bams/$treat"_mm10.rmdupq1.bedpe" |awk '$NF== "-" && $6 - $2 <= 2000{print $0}' |wc -l | awk '{printf $1}' );
scalefactor=$(echo $nread_2k |awk '{print 20000000/$1}') ; ####### normalized to 20M .

cat $treat"_mm10.rmdupq1.bedpe" | awk '$NF == "-" && $6 - $2 < 2000 && $2+$6 > 80{mid=int(($2+$6)/2); OFS="\t"; print $1,mid-40,mid+41,$7,$8,"+"}' | sort-bed --max-mem 8G - | genomeCoverageBed -bg -i - -g /home/xchen2/SWISNF/Xiaolong/mm10.genome -scale $scalefactor |sort -k1,1 -k2,2n  > $treat"_mm10.bedGraph";

bedGraphToBigWig $treat"_mm10.bedGraph" /home/xchen2/SWISNF/Xiaolong/mm10.genome $treat"_mm10.bw"