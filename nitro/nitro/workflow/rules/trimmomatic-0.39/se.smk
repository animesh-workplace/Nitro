rule trimmomatic_se_v0_39:
	input:
		"data/{sample}_{file_id}.fastq.gz",
	output:
		trim = config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}_{file_id}_trimmomatic_trimmed.fastq.gz",
	params: 
		adapters = "resources/adapters/Adapters.fa"
	threads: 1
	log: 
		config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_{file_id}_trimmomatic.log"
	conda:
		"env.yaml"
	wrapper: 
		"file:rules/trimmomatic-0.39/se_wrapper.py"
