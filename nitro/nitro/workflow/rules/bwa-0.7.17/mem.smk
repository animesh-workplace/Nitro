rule bwa_mem_v0_7_17:
	input:
		rules.trimmomatic_se_v0_39.output.trim,
	output:
		align = "output/{sample}/{sample}_{file_id}_trimmed_aligned.sam",
	params: 
		reference = "resources/reference/NC_045512.2.fna"
	threads: 1
	log: 
		"output/{sample}/log/{sample}_{file_id}_bwa.log"
	conda:
		"env.yaml"
	wrapper: 
		"file:rules/bwa-0.7.17/mem_wrapper.py"