process REFORMAT_SAMPLE {
	tag "Reformating $sample.name using $task.cpus CPUs $task.memory"
	label "s_cpu"
	label "xs_mem"

	input:
	val sample

	output:
	tuple val(sample.name), val(sample)

	""" """ //this is not an error
}

process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $name using $task.cpus CPUs and $task.memory memory"
	label "s_cpu"
	label "xs_mem"

	input:
	tuple val(name), val(sample)

	output:
	tuple val(name), val(sample), path("*R{1}*.fastq.gz"). path("*R{2}*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $name
	cp  /mnt/shared/MedGen/sequencing_results/primary_data/*${sample.run}/raw_fastq/${name}*R{1,2}* ./
	"""
} 

process UMI_extract {
	tag "trimming 1 on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/${name}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), path(fwd), path(rev)

	output:
	tuple val(name), path("umi*.fastq.gz")

	script:
	""" 
	echo $fwfq
	echo $revfq
	source activate umitools
	umi_tools extract -I ${fwd} --bc-pattern=NNNNNNNN --read2-in=${rev} --stdout=${name}.umi.R1.fastq.gz --read2-out=${name}.umi.R2.fastq.gz
	"""
}

process TRIMMING_1 {
	tag "trimming 1 on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/${name}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("*.fastq.gz")

	script:
	""" 
	cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
		-o ${name}.trimmed1.R1.fastq.gz -p ${name}.trimmed1.R2.fastq.gz $reads
	"""
}


process TRIMMING_2 {
	tag "trimming 2 on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/${name}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("*.fastq.gz")

	script:
	""" 
	cutadapt -a CAAGGGGGACTGTAGATGGG...TAGGATCTGACTGCGGCTCC \
		-A GGAGCCGCAGTCAGATCCTA...CCCATCTACAGTCCCCCTTG \
		-a ACAACGTTCTGGTAAGGACAX -A TGTCCTTACCAGAACGTTGTX --overlap 4 \
		-o ${name}.trimmed2.R1.fastq.gz -p ${name}.trimmed2.R2.fastq.gz $reads
	"""
}

process FIRST_ALIGN_BAM {
	tag "first align on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/mapped/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
        tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4 ${params.refindex} $reads \
	| samtools view -Sb -o - -| samtools sort -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai	
	"""
}

process PILE_UP {
	tag "PILE_UP on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/vcf/", pattern: '*.md.ba*', mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	tuple val(name), path("*.mpileup")
	
	script:
	"""
	#
	samtools mpileup -x -B -Q 25 -d 999999 -L 999999 -F 0.0005 -f ${params.ref}.fa $bam > ${name}.mpileup
	"""
}



process VARSCAN {
	tag "VARSCAN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/vcf/", mode:'copy'
	
	input:
	tuple val(name), path(mpileup)

	output:
	tuple val(name), path('*.snv.vcf')
	tuple val(name), path('*.indel.vcf')
	
	script:
	"""
	varscan mpileup2snp $mpileup \
		--strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005 \
		--output-vcf | sed 's/Sample1/varscan_SNV/g' > ${name}.varscan.snv.vcf

	varscan mpileup2indel $mpileup \
		 --strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005 \
		--output-vcf | sed 's/Sample1/varscan_INDEL/g' > ${name}.varscan.indel.vcf
	"""
}



process VARDICT {
	tag "VARDICT on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/vcf/", mode:'copy'
	
	input:
	tuple val(name), path(bam)
	tuple val(name), path(bai)


	output:
	tuple val(name), path ("*.vcf")

	script:
	"""
	vardict -G ${params.ref}.fa -f 0.0005 -b ${bam} \
		-c 1 -S 2 -E 3 -r 8 -Q 1 -q 25 -P 2 -m 8 \
		${params.varbed} | Rscript --vanilla ${params.teststrandbias} | perl ${params.var2vcf_valid} -f 0.0005 -d 50 -c 5 -p 2 -q 25 -Q 1 -v 8 -m 8 -N vardict - > ${name}.vardict.vcf
	"""
}



process NORMALIZE_VARIANTS {
	tag "Normalizing variants on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/vcf/", mode:'copy'
	
	input:
	tuple val(name), path (varscan_snv_path)
	tuple val(name), path (varscan_indel_path)
	tuple val(name), path (vardict_path)

	output:
	tuple val(name), path ("*varscan.snv.norm.vcf")
	tuple val(name), path ("*varscan.indel.norm.vcf")
	tuple val(name), path ("*.vardict.norm.vcf")
	
	script:
	"""
	bcftools -v
	bcftools norm -f ${params.ref}.fa $varscan_snv_path -o ${name}.varscan.snv.norm.vcf
	bcftools norm -f ${params.ref}.fa $varscan_indel_path -o ${name}.varscan.indel.norm.vcf
	bcftools norm -f ${params.ref}.fa $vardict_path -o ${name}.vardict.norm.vcf
	
	"""
}


process MERGE_VARIANTS {
	tag "Merging variants on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/vcf/", mode:'copy'
	
	input:
	tuple val(name), path (varscan_snv_norm_path)
	tuple val(name), path (varscan_indel_norm_path)
	tuple val(name), path (vardict_norm_path)

	output:
	tuple val(name), path ("*allcallers.merged.vcf")
	
	script:
	"""
	java -jar ${params.gatk36} -T CombineVariants --variant:vardict ${vardict_norm_path} \
		--variant:varscan_SNV ${varscan_snv_norm_path} \
		--variant:varscan_INDEL ${varscan_indel_norm_path} \
		-R ${params.ref}.fa -genotypeMergeOptions UNSORTED \
		--disable_auto_index_creation_and_locking_when_reading_rods -o ${name}.allcallers.merged.vcf
	"""
}



process NORMALIZE_MERGED_VARIANTS {
	tag "Normalizing merged variants on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/vcf/", mode:'copy'
	
	input:
	tuple val(name), path (merged_vcf)

	output:
	tuple val(name), path ("*allmerged.norm.vcf")
	
	script:
	"""
	bcftools norm -f ${params.ref}.fa $merged_vcf -o ${name}.allmerged.norm.vcf
	"""
}


process ANNOTATE {
	tag "Annotating variants on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/annotate/", mode:'copy'
	
	input:
	tuple val(name), path (merged_normed_vcf)

	output:
	tuple val(name), path ("*allmerged.norm.annot.vcf")
	
	script:
	"""
	vep -i $merged_normed_vcf --cache --cache_version 95 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --offline --vcf --everything -o ${name}.allmerged.norm.annot.vcf
  
	"""
}



// process NORMALIZE_VEP {
// 	tag "Normalizing annotated variants on $name using $task.cpus CPUs and $task.memory memory"
// 	publishDir "${launchDir}/${name}/annotate/", mode:'copy'
	
// 	input:
// 	tuple val(name), path(annotated)

// 	output:
// 	tuple val(name), path("*allmerged.norm.annot.norm.vcf")
	
// 	script:
// 	"""
// 	biopet tool VepNormalizer -I $annotated -O ${name}.allmerged.norm.annot.norm.vcf -m explode  
// 	"""
// }


// process CREATE_TXT {
// 	tag "Simplify annotated table on $name using $task.cpus CPUs and $task.memory memory"
// 	publishDir "${launchDir}/annotate/", mode:'copy'
	
// 	input:
// 	tuple val(name), path(annotated_normed)

// 	output:
// 	tuple val(name), path("*.norm.merged.annot.normVEP.txt")
	
// 	script:
// 	"""
// 	python ${params.vcf_simplify} SimplifyVCF -toType table -inVCF $annotated_normed -out ${name}.norm.merged.annot.normVEP.txt  
// 	"""
// }


// process CREATE_FINAL_TABLE {
// 	tag "Creating final table on $name using $task.cpus CPUs and $task.memory memory"
// 	publishDir "${launchDir}/annotate/", mode:'copy'
	
// 	input:
// 	tuple val(name), path(annotated_normed)
	
// 	script:
// 	"""
// 	Rscript --vanilla ${params.create_table} ${launchDir}/annotate $run  
// 	"""
// }


process CREATE_TXT {
	tag "Simplify annotated table on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/annotate/", mode:'copy'
	
	input:
	tuple val(name), path(annotated_normed)

	output:
	tuple val(name), path("*.norm.merged.annot.normVEP.txt")
	
	script:
	"""
	bcftools norm -m-both $annotated_normed > ${name}.allmerged.norm.vcf
	python ${params.vcf2table_simple_mode} SimplifyVCF -toType table -inVCF $annotated_normed -out ${name}.norm.merged.annot.normVEP.txt  
	"""
}

// process CREATE_FINAL_TABLE {
// 	tag "Creating final table on $name using $task.cpus CPUs and $task.memory memory"
// 	publishDir "${launchDir}/annotate/", mode:'copy'
	
// 	input:
// 	tuple val(name), path(annotated_normed)
	
// 	script:
// 	"""
// 	Rscript --vanilla ${params.create_table} ${launchDir}/annotate $run  
// 	"""
// }

process COVERAGE {
	tag "Creating coverage on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	tuple val(name), path("*.txt")
	
	script:
	"""
	bedtools coverage -abam ${params.covbed} -b $bam -d > ${name}.PBcoverage.txt   
	"""
}

process COVERAGE_STATS {
	tag "Creating coverage stats on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(whatever_navaznost)
	
	script:
	"""
	Rscript --vanilla ${params.coverstat} ${launchDir}/coverage $run  
	"""
}


process MULTIQC {
	tag "first QC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/${name}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	path "report.html"

	script:
	"""
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
        picard BedToIntervalList -I ${params.covbedpicard} -O ${name}.interval_list -SD ${params.ref}.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}.fa -O ${name}.aln_metrics
	multiqc . -n report.html
	"""

}

 
workflow {
	runlist = channel.fromList(params.samples)
	reformatedInput	= REFORMAT_SAMPLE(runlist)
	rawfastq = channel.fromFilePairs("${params.datain}/CLL*R{1,2}*", checkIfExists: true)
	rawFQs =COLLECT_BASECALLED()
	UMI_extract(rawFQs).view()
	// trimmed1	= TRIMMING_1(rawfastq)
	// trimmed2	= TRIMMING_2(trimmed1)	
	// sortedbam	= FIRST_ALIGN_BAM(trimmed2)
	// pileup		= PILE_UP(sortedbam[0])
	// varscanned	= VARSCAN(pileup)
	// vardicted	= VARDICT(sortedbam[0],sortedbam[1])
	// normalized	= NORMALIZE_VARIANTS(varscanned,vardicted)
	// merged		= MERGE_VARIANTS(normalized)
	// norm_merged	= NORMALIZE_MERGED_VARIANTS(merged)
	// annotated	= ANNOTATE(norm_merged)
	// annot_norm	= NORMALIZE_VEP(annotated)
	// txt		= CREATE_TXT(annot_norm)
	// final_table	= CREATE_FINAL_TABLE(txt)

	// covered		= COVERAGE(sortedbam[0])
	// COVERAGE_STATS(covered)					
	// MULTIQC(sortedbam[0])	
}
