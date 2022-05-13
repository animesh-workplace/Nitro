rule samtools_sort_v1_15_1:
	input:
		rules.bwa_mem_v0_7_17.output.align,
	output:
		sort = "output/{sample}/{sample}_{file_id}_trimmed_aligned_sorted.bam",
	params: 
	threads: 1
	log: 
		"output/{sample}/log/{sample}_{file_id}_sort.log"
	conda:
		"env.yaml"
	wrapper: 
		"file:rules/samtools-1.15.1/sort_wrapper.py"