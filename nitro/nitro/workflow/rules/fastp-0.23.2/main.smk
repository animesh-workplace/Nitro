rule fastp_main_v0_23_2:
	input:
		"data/{sample}_{file_id}.fastq.gz",
	output:
		trim = "output/{sample}/{sample}_{file_id}_fastp_trimmed.fastq.gz",
		json = "output/{sample}/fastp_report/{sample}{file_id}_fastp_trimmed_report.json",
		html = "output/{sample}/fastp_report/{sample}{file_id}_fastp_trimmed_report.html",
	params: 
		"-q 15 -u 40 -l 25 --cut_right --cut_window_size 20 --cut_mean_quality 30"
	threads: 1
	log: 
		"output/{sample}/log/{sample}_{file_id}_fastp.log"
	conda:
		"env.yaml"
	wrapper: 
		"file:rules/fastp-0.23.2/main_wrapper.py"