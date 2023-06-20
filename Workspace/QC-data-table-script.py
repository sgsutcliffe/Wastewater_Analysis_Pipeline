#!/bin/python3.6
#Version QC-data-table-script-V1
#Script to collate QC-Metrics for analysis to use in building a report
#Version2
#Make the script handle samples in parallel

import sys
import pandas as pd
import glob
import subprocess
import os
import multiprocessing as mp
from multiprocessing import Pool, Manager
import time
from functools import partial

#Function that test if a string is numeric
def is_numeric(the_str):
	try:
		float_conv = float(the_str)
		return True
	except (ValueError, TypeError):
		return False
def pull_qc_info(lst, sample):
	globals()["row_"+str(sample)] = []
	globals()["row_"+str(sample)].append(sample)
	if glob.glob("{0}{1}.fastp.json".format(qc_input, sample)):
		globals()["row_"+str(sample)].append(subprocess.getoutput("cat {0}{1}.fastp.json | jq .summary.before_filtering.total_reads".format(qc_input, sample)))
		globals()["row_"+str(sample)].append(subprocess.getoutput("cat {0}{1}.fastp.json | jq .summary.before_filtering.gc_content".format(qc_input,sample)))
	else:
		globals()["row_"+str(sample)]("NA")
		globals()["row_"+str(sample)].append("NA")
	if glob.glob(qc_input + sample + ".initial.kraken2.report.txt"):
		step3_cmd = "awk '{if (($4 == \"S\") && ($5 == \"9606\")) {print} }' " + qc_input + sample + ".initial.kraken2.report.txt | cut -f1"
		if subprocess.getoutput(step3_cmd).strip() == "":
			globals()["row_"+str(sample)].append("NA")
		else:
			globals()["row_"+str(sample)].append(float(subprocess.getoutput(step3_cmd).strip()))
		step4_cmd = "awk '{if (($4 == \"S\") && ($5 == \"694009\")) {print} }' " + qc_input + sample + ".initial.kraken2.report.txt | cut -f1"
		if subprocess.getoutput(step4_cmd).strip() == "":
			globals()["row_"+str(sample)].append("NA")
		else:
			globals()["row_"+str(sample)].append(float(subprocess.getoutput(step4_cmd).strip()))
	else:
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
	if glob.glob("{0}{1}_ivartrim_sorted.bam".format(analysis_input, sample)):
		globals()["row_"+str(sample)].append(subprocess.getoutput("bedtools coverage -a {0}genome.bed -b {1}{2}_ivartrim_sorted.bam | cut -f7 ".format(ref_path, analysis_input, sample)))
		globals()["row_"+str(sample)].append(subprocess.getoutput("bedtools coverage -mean -a {0}genome.bed -b {1}{2}_ivartrim_sorted.bam | cut -f4 ".format(ref_path, analysis_input, sample)))
		step7_cmd = "bedtools coverage -mean -a " + insertbed + " -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '{ sum += $7; ++n} END {print sum/n+0}'"
		globals()["row_"+str(sample)].append(float(subprocess.getoutput(step7_cmd).strip()))
		step8_cmd = "bedtools coverage -mean -a " + insertbed + " -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '($7>= 100) {++n} END {print n+0}'"
		globals()["row_"+str(sample)].append(float(subprocess.getoutput(step8_cmd).strip()))
		step9_cmd = "bedtools coverage -mean -a " + ref_path + "SARS-CoV-2_ORF_full.bed -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '{ sum += $5; ++n} END { print sum /n}'"
		globals()["row_"+str(sample)].append(float(subprocess.getoutput(step9_cmd).strip()))
		step10_cmd = "bedtools coverage -mean -a " + ref_path + "SARS-CoV-2_ORF_full.bed -b " + analysis_input + sample + "_ivartrim_sorted.bam | awk '($5>= 100) {++n} END {print n+0}'"
		globals()["row_"+str(sample)].append(float(subprocess.getoutput(step10_cmd).strip()))
		globals()["row_"+str(sample)].append(subprocess.getoutput("bedtools coverage -a {0}spike.bed -b {1}{2}_ivartrim_sorted.bam | cut -f8 ".format(ref_path, analysis_input, sample)))
		globals()["row_"+str(sample)].append(subprocess.getoutput("bedtools coverage -mean -a {0}spike.bed -b {1}{2}_ivartrim_sorted.bam | cut -f5 ".format(ref_path, analysis_input, sample)))
	else:
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
	if glob.glob(analysis_input + sample + "_variantout.tsv"):
		step11_cmd = "awk '($13 <= 0.05) {++n} END {print n+0}' " + analysis_input + sample + "_variantout.tsv"
		globals()["row_"+str(sample)].append(float(subprocess.getoutput(step11_cmd).strip()))
		step12_cmd = " awk '$13 <= 0.05' " + analysis_input + sample + "_variantout.tsv | awk '{ sum += $12; n++ } END { if (n >0) print sum / n; else print 0; }'"
		globals()["row_"+str(sample)].append(float(subprocess.getoutput(step12_cmd).strip()))
	else:
		globals()["row_"+str(sample)].append("NA")
		globals()["row_"+str(sample)].append("NA")
	
	lst.append(pd.DataFrame(data=[globals()["row_"+str(sample)]]))

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
	
	lst_samples=[]
	with open(samples_list_file) as f:
		lst_samples = f.readlines()
	the_idx = 0
	for sample_name in lst_samples:
		lst_samples[the_idx] = sample_name.replace("\n", "")
		the_idx = the_idx + 1
	num_samples = len(lst_samples)

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

	#7 CPUs requested, we need 1 cpu per sample analyzed in parallel to run efficiently with bedtools
	requested_cpu = int(sys.argv[7])

	if requested_cpu >= os.cpu_count():
		nb_sim_process = os.cpu_count() - 1
	else:
		nb_sim_process = requested_cpu() -1

	#Set up a dataframe to receive the output 
	column_names = [
		'Sample Name'
		'Number of Reads',
		'Percentage GC',
		'Percentage of Human Reads',
		'Percentage of SARS-CoV-2 reads',
		'Breadth of coverage (Genome)',
		'Mean depth (Genome)',
		'Mean depth (Amplicon)',
		'Number of Amplicons mean depth >=100x',
		'Mean depth (ORFs)',
		'Number of ORFs mean depth >=100x',
		'Breadth of coverage (Spike)',
		'Mean depth (Spike)',
		'Number of signficant SNVs',
		'Mean depth of all significant mutations'
		]
	
	df = pd.DataFrame(columns=column_names)
	
	#With 1 cpu needed per sample to be analyzed in parallel

	dfs_list = Manager().list()
	pool = Pool(processes=nb_sim_process)
	samples = lst_samples
	res = pool.map_async(partial(pull_qc_info, dfs_list), samples)
	res.wait()
	dfs = pd.concat(dfs_list, ignore_index = True)
	output_file_name = script_output + 'QC-Summary-Report.tsv'
	dfs.to_csv(output_file_name, sep='\t')

