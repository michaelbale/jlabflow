/*
 * Remove chrM and forbidden list (merge ENCODE + CnR) from deduped bam
 * Scope: DNA
 * Input: Forbidden list (from params)
 * Input: Deduped bam with sample ID
 * Emits: final bam file for downstream apps
 * Feeds: RGPC bigwig
 */
 
process FINALFILTER {
   tag "Removing chrM and Forbidden List from ${id}"
   label 'big_mem'
   publishDir "$params.outdir/finalBam", mode: 'copy', pattern: "${id}_final.bam"
   
   input:
   path(filterList)
   tuple val(id), path(bam)
    
   output:
   tuple val(id), file("${id}_final.bam"), emit: final_bams
   path("${id}_filterCount.txt"), optional: true, emit: forbid_list_count
   path("${id}_finalCount.txt"), emit: final_count

   script:
   def filterARGS = params.SE ? '' : '-f 3 -F 8'
   """
   samtools index ${bam}
   export CHROMOSOMES=\$(samtools view -H ${bam} \
	 | \
	 grep '^@SQ' | \
	 cut -f 2 | \
	 grep -v -e _ -e chrM -e 'VN:' | \
	 sed 's/SN://' | \
	 xargs echo)
   samtools view -b -h $filterARGS -F 3332 -q 30 \
     ${bam} \$CHROMOSOMES > tmp.bam
   echo $id \$(bedtools intersect -a $bam -b ${filterList} -ubam -u | samtools view -c -@ task.cpus) > ${id}_filterCount.txt
   bedtools subtract -A -a tmp.bam -b ${filterList} | samtools sort -@ $task.cpus - > ${id}_final.bam
   echo ${id} "\$(samtools view -@ $task.cpus -c ${id}_final.bam)" > ${id}_finalCount.txt
   """
}
