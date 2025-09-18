#/bin/bash
for fq in $(cat fqfiles.txt)
do
fastp -i ${fq}_1.fq.gz -o ${fq}_fil_R1.fq.gz -I ${fq}_2.fq.gz -O ${fq}_fil_R2.fq.gz --detect_adapter_for_pe 2>${fq}.fastp.log
bowtie2 -q -p 20 --very-sensitive -x /data/liuxinhua/Reference/Mouse/mm10/mm10index/mm10 -1 ${fq}_fil_R1.fq.gz -2 ${fq}_fil_R2.fq.gz 2>${fq}.bowtie.log | samtools view -@ 20 -u -O BAM - | samtools sort -@ 20 -O BAM -o ${fq}.srt.bam -
samtools index ${fq}.srt.bam
samtools view -h -b -q 5 ${fq}.srt.bam > ${fq}.srt.rmMulti.bam
samtools index ${fq}.srt.rmMulti.bam
samtools sort -@ 20 -n -O BAM -o ${fq}.sbn.rmMulti.bam ${fq}.srt.rmMulti.bam
bedtools bamtobed -bedpe -i ${fq}.sbn.rmMulti.bam > ${fq}.sbn.rmMulti.bed
awk '$1==$4 && $6-$2 <= 2000 {print $0}' ${fq}.sbn.rmMulti.bed > ${fq}.sbn.rmMulti.clean.bed
cut -f 1,2,6 ${fq}.sbn.rmMulti.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${fq}.sbn.rmMulti.frag.bed
awk '$1 ~ /^(chr[1-9XY]+)$/ {print}' ${fq}.sbn.rmMulti.frag.bed > ${fq}.sbn.rmMulti.clean.frag.bed
bedtools genomecov -bg -i ${fq}.sbn.rmMulti.clean.frag.bed -g /data/liuxinhua/Reference/Mouse/mm10/mm10.chrsize > ${fq}.sbn.rmMulti.bedGraph
/home/liuxinhua/anaconda3/envs/python38/bin/bamCoverage -p 20 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 -b ${fq}.srt.rmMulti.bam --extendReads -o ${fq}.bw
done