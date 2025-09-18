for fq in $(cat fqfiles.txt)
do
fastp -i ${fq}_R1.fq.gz -o ${fq}_R1_fil.fq.gz -I ${fq}_R2.fq.gz -O ${fq}_R2_fil.fq.gz --detect_adapter_for_pe 2>${fq}.fastp.log
STAR --runThreadN 20 --genomeDir /data/liuxinhua/Reference/Mouse/mm10/STAR-mouse --readFilesIn ${fq}_R1_fil.fq.gz ${fq}_R2_fil.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix ${fq} 2>${fq}.STAR.log
samtools view -@ 20 -bS -F 4 ${fq}Aligned.out.sam > ${fq}_rmMulti.bam
samtools sort -@ 20 -n ${fq}_rmMulti.bam -O BAM -o ${fq}_sbnrmMulti.bam
samtools index ${fq}_sbnrmMulti.bam
htseq-count -f bam ${fq}_sbnrmMulti.bam /data/liuxinhua/Reference/Mouse/mm10/mm10.ncbiRefSeq.gtf > ${fq}.count
done