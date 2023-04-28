#!/bin/python3.6
#Version QC-data-table-script-V1
#Script to collate QC-Metrics for analysis to use in building a report
import sys
import pandas as pd
import glob
import subprocess

#Function that test if a string is numeric
def is_numeric(the_str):
	try:
		float_conv = float(the_str)
		return True
	except (ValueError, TypeError):
		return False

if __name__ == "__main__":
	
	#Initializing Global Variables
	
	#1 Path to databases and reference files
	ref_path =sys.argv[1] 
	#Making sure to have "/" at the end of the workspace absolute path
	if ((ref_path[len(ref_path)-1]) != "/"):
		ref_path = "{0}/".format(ref_path)
	print("Path to databases and reference files is : '{0}'".format(ref_path))
	
	#2 Path to analysis outputs that will be the input for metrics
	analysis_input =sys.argv[2] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((analysis_input[len(analysis_input)-1]) != "/"):
		analysis_input = "{0}/".format(analysis_input)
	print("Analysis output directory to be used in input is : '{0}'".format(analysis_input))
	
	#3 Input text file with all samples to be included to create the variable for list of samples 
	samples_list_file = ref_path + sys.argv[3]
	
	#4 Path to QC outputs that will be the input for metrics
	qc_input =sys.argv[4] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((qc_input[len(qc_input)-1]) != "/"):
		qc_input = "{0}/".format(qc_input)
	print("QC pipeline output directory to be used in input is : '{0}'".format(qc_input))

	#5 Output of script directory path
	script_output =sys.argv[5] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((script_output[len(script_output)-1]) != "/"):
		script_output = "{0}/".format(script_output)
	print("Files from this script will be located at : '{0}'".format(script_output))

	#6 Arctic insert bed
	insertbed =sys.argv[6]

	#Now loop through every sample and perform the collection of QC-data for the ten categories
	data = pd.read_csv(samples_list_file, header = None, sep='\t')

    
	# Step 1: Get the number of reads in the sample starting out
	num_of_reads = []
	# Step 2: Get sample %GC content
	perc_GC = []
	# Step 3: Find the percentage of reads align to the human genome
	perc_human = []
	# Step 4: Find the percentage of SARS-CoV-2 Reads
	perc_sars = []
	# Step 5: Breadth of coverage of genome (percentage)
	genome_breadth = []
	# Step 6: Mean read depth of genome
	genome_depth = []
	# Step 7: Mean read depth of 103 amplicons
	amplicon_depth = []
	# Step 8: Number of amplicons with mean depth >= 100x
	amplicon_100xdepth = []
	# Step 9: Mean read depth of 16 ORF (orf1a, orf1b, S, ORF3a, ORF3b, ORF3c, E, M, ORF6, ORF7a, ORF7b, ORF8, N, ORF9c, ORF10, 3'UTR)
	ORF_depth = []
	# Step 10: Number of ORF with mean depth >= 100x
	ORF_100xdepth = []
	# Step 11: Number of signficant mutations
	num_sig_mutations = []
	# Step 12: Mean Depth of All Signficant mutations
	avg_depth_mutations = []
	# Step 13: Mean Depth of spike protein
	spike_depth = []
	# Step 14: Breadth of spike coverage
	spike_breadth = []


	for index, row in data.iterrows():
		sample = row[0]
		print(sample)
		if glob.glob("{0}{1}.fastp.json".format(qc_input, sample)):
			num_of_reads.append(subprocess.getoutput("cat {0}{1}.fastp.json | jq .summary.before_filtering.total_reads".format(qc_input, sample)))
			perc_GC.append(subprocess.getoutput("cat {0}{1}.fastp.json | jq .summary.before_filtering.gc_content".format(qc_input,sample)))
		else:
			num_of_reads.append("NA")
			perc_GC.append("NA")
		if glob.glob(qc_input + sample + ".initial.kraken2.report.txt"):
			step3_cmd = "awk '{if (($4 == \"S\") && ($5 == \"9606\")) {print} }' " + qc_input + sample + ".initial.kraken2.report.txt | cut -f1"
			if subprocess.getoutput(step3_cmd).strip() == "":
				perc_human.append("NA")
			else:
				perc_human.append(float(subprocess.getoutput(step3_cmd).strip()))
			step4_cmd = "awk '{if (($4 == \"S\") && ($5 == \"694009\")) {print} }' " + qc_input + sample + ".initial.kraken2.report.txt | cut -f1"
			if subprocess.getoutput(step4_cmd).strip() == "":
				perc_sars.append("NA")
			else:
				perc_sars.append(float(subprocess.getoutput(step4_cmd).strip()))
		else:
			perc_human.append("NA")
			perc_sars.append("NA")
		if glob.glob("{0}{1}_ivartrim_sorted.bam".format(analysis_input, sample)):
			genome_breadth.append(subprocess.getoutput("bedtools coverage -a {0}genome.bed -b {1}{2}_ivartrim_sorted.bam | cut -f7 ".format(ref_path, analysis_input, sample)))
			genome_depth.append(subprocess.getoutput("bedtools coverage -mean -a {0}genome.bed -b {1}{2}_ivartrim_sorted.bam | cut -f4 ".format(ref_path, analysis_input, sample)))
			step7_cmd = "bedtools coverage -mean -a " + ref_path + insertbed + " -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '{ sum += $5 } END { print sum /103}'"
			amplicon_depth.append(float(subprocess.getoutput(step7_cmd).strip()))
			step8_cmd = "bedtools coverage -mean -a " + ref_path + insertbed + " -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '($5>= 100) {++n} END {print n+0}'"
			amplicon_100xdepth.append(float(subprocess.getoutput(step8_cmd).strip()))
			step9_cmd = "bedtools coverage -mean -a " + ref_path + "SARS-CoV-2_ORF_full.bed -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '{ sum += $5 } END { print sum /16}'"
			ORF_depth.append(float(subprocess.getoutput(step9_cmd).strip()))
			step10_cmd = "bedtools coverage -mean -a " + ref_path + "SARS-CoV-2_ORF_full.bed -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '($5>= 100) {++n} END {print n+0}'"
			ORF_100xdepth.append(float(subprocess.getoutput(step10_cmd).strip()))
			spike_breadth.append(subprocess.getoutput("bedtools coverage -a {0}spike.bed -b {1}{2}_ivartrim_sorted.bam | cut -f7 ".format(ref_path, analysis_input, sample)))
			spike_depth .append(subprocess.getoutput("bedtools coverage -mean -a {0}spike.bed -b {1}{2}_ivartrim_sorted.bam | cut -f4 ".format(ref_path, analysis_input, sample)))
		else:
			genome_breadth.append("NA")
			genome_depth.append("NA")
			spike_breadth.append("NA")
			spike_depth.append("NA")
			amplicon_depth.append("NA")
			amplicon_100xdepth.append("NA")
			ORF_100xdepth.append("NA")
			ORF_depth.append("NA")
		if glob.glob(analysis_input + sample + "_variantout.tsv"):
			step11_cmd = "awk '($13 <= 0.05) {++n} END {print n+0}' " + analysis_input + sample + "_variantout.tsv"
			num_sig_mutations.append(float(subprocess.getoutput(step11_cmd).strip()))
			step12_cmd = " awk '$13 <= 0.05' " + analysis_input + sample + "_variantout.tsv | awk '{ sum += $12; n++ } END { if (n >0) print sum / n; else print 0; }'"
			avg_depth_mutations.append(float(subprocess.getoutput(step12_cmd).strip()))
		else:
			num_sig_mutations.append("NA")
			avg_depth_mutations.append("NA")


	data['Number of Reads'] = num_of_reads
	data['Percentage GC'] = perc_GC
	data['Percentage of Human Reads'] = perc_human
	data['Percentage of SARS-CoV-2 reads'] = perc_sars
	data['Breadth of coverage (Genome)'] = genome_breadth
	data['Mean depth (Genome)'] = genome_depth
	data['Mean depth (Amplicon)'] = amplicon_depth
	data['Number of Amplicons mean depth >=100x'] = amplicon_100xdepth
	data['Mean depth (ORFs)'] = ORF_depth 
	data['Number of ORFs mean depth >=100x'] = ORF_100xdepth
	data['Breadth of coverage (Spike)'] = spike_breadth
	data['Mean depth (Spike)'] = spike_depth
	data['Number of signficant SNVs'] = num_sig_mutations
	data['Mean depth of all significant mutations'] = avg_depth_mutations
	output_file_name = script_output + 'QC-Summary-Report.tsv'
	data.to_csv(output_file_name, sep='\t')

