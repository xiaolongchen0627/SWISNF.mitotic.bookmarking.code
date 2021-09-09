
# calculate overlap between Asyn and Mitotic cells 

asyn=$1 # merged peak for Asyn cell
mit=$2 # merged peak for Mitotic cell

assayName=$3 #### sample name eg SOX2 , ESRRB
bedmap --echo --echo-map --delim "\t" $asyn $mit > intersect/Asyn.$assayName.overlap.bed
bedmap --echo --echo-map --delim "\t" $mit $asyn > intersect/Mit.$assayName.overlap.bed

awk -v Name=$assayName '{if(NF==3){overlap="Asyn_only"}if(NF>3){overlap="overlap"} print $1"\t"$2"\t"$3"\t"Name"\t"overlap}' intersect/Asyn.$assayName.overlap.bed  > intersect/Asyn.$assayName.tmp.bed
awk -v Name=$assayName 'NF==3{print $1"\t"$2"\t"$3"\t"Name"\tMit_only"}' intersect/Mit.$assayName.overlap.bed  > intersect/Mit.$assayName.tmp.bed

cat intersect/Asyn.$assayName.tmp.bed intersect/Mit.$assayName.tmp.bed > intersect/$assayName.3types.peaks.bed
mkdir intersect/$assayName
awk 'BEGIN{OFS="\t"}{print $0 >> "intersect/"t4"/"$5".peaks.bed"}' intersect/$assayName.3types.peaks.bed 