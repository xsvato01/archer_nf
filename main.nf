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
	tuple val(name), val(sample), path("*R1*.fastq.gz"), path("*R2*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $name
	cp  /mnt/shared/MedGen/sequencing_results/primary_data/*${sample.run}/raw_fastq/${name}* ./
	"""
} 

process UMI_extract {
	tag "UMI_extract on $name using $task.cpus CPUs and $task.memory memory"
	// publishDir  "${launchDir}/${name}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
	tuple val(name), val(sample), path("${name}.umi*.fastq.gz")

	script:
	""" 
	echo UMI_extract $name
	source activate umi_tools
	umi_tools extract -I ${fwd} --bc-pattern=NNNNNNNN --read2-in=${rev} --stdout=${name}.umi.R1.fastq.gz --read2-out=${name}.umi.R2.fastq.gz
	"""
}

process TRIMMING {
	tag "TRIMMING on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/${name}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("*.fastq.gz")

	script:
	"""
	echo TRIMMING $name
	source activate cutadapt
	cutadapt -g AACCGCCAGGAGT -m 50 -o ${name}.trimmed1.R1.fastq.gz -p ${name}.trimmed1.R2.fastq.gz $reads
	"""
}


process FIRST_ALIGN_BAM {
	tag "FIRST_ALIGN_BAM on $name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${launchDir}/${name}/mapped/", mode:'copy'
	label "m_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(reads)

	output:
    tuple val(name), val(sample), path("${name}.bam")
	//, path("${name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo FIRST_ALIGN_BAM $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads > ${name}.sam
	samtools view -Sb ${name}.sam -o ${name}.bam
	"""
}

process SORT_INDEX {
	tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/${name}/mapped/", mode:'copy'
	label "m_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam)

	output:
	tuple val(name), val(sample), path("${name}.sorted.bam"), path("${name}.sorted.bai")


	script:
	"""
	echo SORT_INDEX $name
	source activate samtools
	samtools sort $bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}


process DEDUP {
	tag "DEDUP on $name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${launchDir}/${name}/mapped/", mode:'copy'
	label "m_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
    tuple val(name), val(sample), path("${name}.deduped.bam"), path("${name}.deduped.bai")

	script:
	"""
	echo DEDUP $name
	source activate umi_tools
	umi_tools dedup -I $bam --paired -S ${name}.deduped.bam
	samtools index ${name}.deduped.bam ${name}.deduped.bai
	"""
}

process MUTECT2 {
	tag "MUTECT2 on $name using $task.cpus CPUs and $task.memory memory"
	label "m_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam), path(bai)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.vcf")

	script:
	"""
	echo MUTECT2 $name
	source activate gatk4
	gatk Mutect2 --reference ${params.ref}.fa --input ${bam} --tumor-sample $name --annotation StrandArtifact --min-base-quality-score 20 --intervals $params.ivl --output ${sample.name}.mutect.vcf
	"""
}

process FILTER_MUTECT {
	tag "FILTER_MUTECT on $name using $task.cpus CPUs and $task.memory memory"
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.filt.vcf")

	script:
	"""
	echo FILTER_MUTECT $name
	source activate gatk4
	gatk FilterMutectCalls -V $vcf_input -O ${sample.name}.mutect.filt.vcf
	"""
}

process NORMALIZE_MUTECT {
	tag "NORMALIZE_MUTECT on $name using $task.cpus CPUs $task.memory"
	label "s_mem"
	label "s_cpu"
	container "staphb/bcftools:1.10.2"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.filt.norm.vcf")

	script:
	"""
	echo NORMALIZE_MUTECT $name
	bcftools norm -m-both $vcf_input > ${sample.name}.mutect.filt.norm.vcf
	"""
}

process ANNOTATE_MUTECT {
	tag "ANNOTATE_MUTECT on $name using $task.cpus CPUs $task.memory"
	// publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${sample.name}.mutect2.filt.norm.vep.vcf")

	script:
	"""
	echo ANNOTATE_MUTECT $name
	source activate vep
	vep -i $vcf_input --cache --cache_version 95 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --offline --vcf --everything -o ${sample.name}.mutect2.filt.norm.vep.vcf
	"""	
}

process FILTER_VCF {
	tag "FILTER_VCF on $name using $task.cpus CPUs $task.memory"
	container "staphb/bcftools:1.10.2"
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${sample.name}.mutect2.filt.norm.vep.filt.vcf")

	script:
	"""
	echo FILTER_VCF $name
	bcftools view -f 'PASS,clustered_events,multiallelic' $vcf_input > ${sample.name}.mutect2.filt.norm.vep.filt.vcf
	"""	
}

process VCF2CSV {
	tag "VCF2CSV on $name using $task.cpus CPUs $task.memory"
	// publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${name}.csv")

	script:
	"""
	echo VCF2CSV $name
	source activate vcf2csv
	python $params.vcf2csv simple --build GRCh37 -i $vcf_input -o ${name}.csv
	"""	
}


process FLT3 {
	tag "FLT3 on $name using $task.cpus CPUs $task.memory"
	// publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "s_mem"
	label "s_cpu"
	// errorStrategy 'ignore'
	debug true

	input:
	tuple val(name), val(sample), path(bam), path(bai)
	
	output:
	tuple val(name), val(sample), path("${name}.deduped_FLT3_ITD_summary.txt")

	script:
	"""
	echo FLT3 $name
	source activate perl
    tar -C /tmp -xf $params.flt3tar
	perl /tmp/FLT3/FLT3_ITD_ext/FLT3_ITD_ext.pl --bam $bam --output ./ --ngstype amplicon --genome hg19 --fgbiojar $params.fgbio --picardjar $params.picard --refindex /tmp/FLT3/FLT3_bwaindex/FLT3_dna_e1415
	ls -alh
	"""	
}


process BAMQC {
	tag "BAMQC on $name using $task.cpus CPUs and $task.memory memory"
	container 'registry.gitlab.ics.muni.cz:443/450402/qc_cmbg:26'
	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	path "*"

	script:
	"""
	echo BAMQC $name
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=$params.ivl TARGET_INTERVALS=$params.ivl R=${params.ref}.fa O=${name}.hs_metrics
	picard CollectAlignmentSummaryMetrics I=$bam R=${params.ref}.fa O=${name}.aln_metrics
	"""
}

process MULTIQC {
	tag "MultiQC on all samples using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outdir}/multiqc/", mode:'copy'
	label "m_mem"
	label "s_cpu"

	input:
	path('*')

	output:
	path "*"

	script:
	"""
	multiqc . -n MultiQC.html
	"""

}



workflow {
	runlist = channel.fromList(params.samples)
	reformatedInput	= REFORMAT_SAMPLE(runlist)
	rawFQs = COLLECT_BASECALLED(reformatedInput)
	umiFQs = UMI_extract(rawFQs)
	trimmedFQs = TRIMMING(umiFQs)
	firstBAM = FIRST_ALIGN_BAM(trimmedFQs)
	sortedBamBai = SORT_INDEX(firstBAM)
	dedupedBam = DEDUP(sortedBamBai)
	BAMQC(dedupedBam)

	rawVCFs = MUTECT2(dedupedBam)
	filteredVCFs = FILTER_MUTECT(rawVCFs)
	normedVCFs = NORMALIZE_MUTECT(filteredVCFs)
	annotatedVCFs = ANNOTATE_MUTECT(normedVCFs)
	filteredAnnotatedVCFs = FILTER_VCF(annotatedVCFs)
	csv = VCF2CSV(filteredAnnotatedVCFs)

	FLT3(dedupedBam).view()

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
