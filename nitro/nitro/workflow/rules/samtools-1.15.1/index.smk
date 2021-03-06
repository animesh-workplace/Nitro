rule samtools_index_v1_15_1:
	input:
		rules.samtools_merge_v1_15_1.output.merge,
	output:
		index = config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}_merged.bam.bai",
	params: 
	threads: 1
	log: 
		config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_index.log"
	conda:
		"env.yaml"
	wrapper: 
		"file:rules/samtools-1.15.1/index_wrapper.py"
