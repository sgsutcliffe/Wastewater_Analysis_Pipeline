#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File-Checking Version3
Version3:
Runs samples in parallel
#Version 2: Integrates into Freyja parallelization pipeline
@author: Steven Sutcliffe
"""

import sys
import pandas as pd
import glob
import subprocess
import os
import multiprocessing as mp
from multiprocessing import Pool, Manager
import time
from functools import partial

def file_check(lst, sample):
	path_to_results = path_to_output
	globals()["row_"+str(sample)] = []
	globals()["row_"+str(sample)].append(sample)
	#1 Check the FASTQ files (Forward & Reverse)
	if glob.glob(path_to_samples + sample + "_R2_001.fastq.gz"):
		if glob.glob(path_to_samples + sample + "_R1_001.fastq.gz"):
			globals()["row_"+str(sample)].append("True")
		else:
			globals()["row_"+str(sample)].append('False')
	elif glob.glob(path_to_samples + sample + "_R1_001.fastq.gz"):
		if glob.glob(path_to_samples + sample + "_R2_001.fastq.gz"):
			globals()["row_"+str(sample)].append("True")
		else:
			globals()["row_"+str(sample)].append('False')
	else:
		globals()["row_"+str(sample)].append('False')
	#2 Check for fastp output trimmed files (forward & Reverse)
	if glob.glob(path_to_results + sample + "_R1_001_trimmed_1.fastq") or glob.glob(path_to_results + sample + "_R1_001_trimmed_1.fastq.gz"):
		if glob.glob(path_to_results + sample + "_R2_001_trimmed_2.fastq*") or glob.glob(path_to_results + sample + "_R1_001_trimmed_1.fastq.gz"):
			globals()["row_"+str(sample)].append("True")
		else:
			globals()["row_"+str(sample)].append("False")
	else:
		globals()["row_"+str(sample)].append('False')
	#3 Check for decontaminated human reads
	if glob.glob(path_to_results + sample + "*_R1_001_decon_1.fastq") or glob.glob(path_to_results + sample + "*R1_001_decon_1.fastq.gz"):
		if glob.glob(path_to_results + sample + "*_R2_001_decon_2.fastq*") or glob.glob(path_to_results + sample + "*_R2_001_decon_2.fastq.gz"):
			globals()["row_"+str(sample)].append("True")
		else:
			globals()["row_"+str(sample)].append("False")
	else:
		globals()["row_"+str(sample)].append('False')
	#4 [pre iVar trim] Check to see if any remaining reads can align to reference SARS-CoV-2 genome
	if glob.glob(path_to_results + sample + "*_preprocessed_sorted.bam"):
		bampath = "ls " + path_to_results + sample + "*_preprocessed_sorted.bam"
		bamfile = subprocess.getoutput(bampath)
		samcommand = "samtools view -c " + bamfile
		stdout = subprocess.getoutput(samcommand)
		num_reads = int(stdout.splitlines()[-1])
		print(bamfile)
		print('Number of reads is ' + str(num_reads))
		num_reads = int(num_reads)

		if num_reads > 0:
			globals()["row_"+str(sample)].append("True")
		else:
			globals()["row_"+str(sample)].append('False')
	else:
		globals()["row_"+str(sample)].append('False')
	#5 [post iVar trim] Check to see if any remaining reads can align to reference SARS-CoV-2 genome
	if glob.glob(path_to_results + sample + "*_ivartrim_sorted.bam"):
		bampath = "ls " + path_to_results + sample + "*_ivartrim_sorted.bam"
		bamfile = subprocess.getoutput(bampath)
		samcommand = "samtools view -c " + bamfile
		stdout = subprocess.getoutput(samcommand)
		num_reads = int(stdout.splitlines()[-1])
		print(bamfile)
		print('Number of reads is ' + str(num_reads))
		num_reads = int(num_reads)

		if num_reads > 0:
			globals()["row_"+str(sample)].append("True")

		else:
			globals()["row_"+str(sample)].append('False')

	else:
		globals()["row_"+str(sample)].append('False')
	#6 Check Freyja variant output
	depthout = "ls " + path_to_results + sample + "*_depthout"
	depthout = subprocess.getoutput(depthout)
	variantout = "ls " + path_to_results + sample + "*_variantout.tsv"
	variantout = subprocess.getoutput(variantout)
	if glob.glob(depthout) and glob.glob(variantout):
		if (os.path.getsize(depthout) == 0) or (os.path.getsize(variantout) == 0):
			globals()["row_"+str(sample)].append("False")
		else:
			globals()["row_"+str(sample)].append('True')
	else:
		globals()["row_"+str(sample)].append('False')

	#7 Check Freyja demix output exists
	if glob.glob(path_to_results + "results/" + sample + "*_output"):
		outputfile = "ls " + path_to_results + "results/" + sample + "*_output"
		outputfile = subprocess.getoutput(outputfile)
		if (os.path.getsize(outputfile) == 0):
			globals()["row_"+str(sample)].append("False")
			globals()["row_"+str(sample)].append("False")
		else: 
			globals()["row_"+str(sample)].append("True")
			freyja_out = pd.read_csv(outputfile, header = None, sep='\t', skiprows=1)
			print(freyja_out)
			cov = float(freyja_out.iloc[4][1])
			if cov <= 0:
				globals()["row_"+str(sample)].append("False")
			else:
				globals()["row_"+str(sample)].append("True")
				
			
		
	else:
		globals()["row_"+str(sample)].append("False")
		globals()["row_"+str(sample)].append("False")

	lst.append(pd.DataFrame(data=[globals()["row_"+str(sample)]]))
		
if __name__ == "__main__":
 
	#Initializing Global Variables

	#Using the the sample list file as the names of files to check for
	path_to_file = sys.argv[1] #Path to sample list file 
	#Output folder of all the analysis that is being checked
	path_to_output = sys.argv[2]
	#Original FASTQ files that are the input
	path_to_samples = sys.argv[3]
	filename = path_to_output + 'File_Check_Output.txt'
	data = pd.read_csv(path_to_file, header = None, sep='\t')

	#Step 1: Convert the sample-list into samples so they can be run in parallel
	lst_samples=[]
	with open(path_to_file) as f:
		lst_samples = f.readlines()
	the_idx = 0
	for sample_name in lst_samples:
		lst_samples[the_idx] = sample_name.replace("\n", "")
		the_idx = the_idx + 1
	num_samples = len(lst_samples)

	#7 CPUs requested, we need 1 cpu per sample analyzed in parallel to run efficiently with bedtools
	requested_cpu = int(sys.argv[4])

	if requested_cpu >= os.cpu_count():
		nb_sim_process = os.cpu_count() - 1
	else:
		nb_sim_process = requested_cpu -1

	#Set up a dataframe to receive the output 
	column_names = [
		'sample'
		'FASTQ'
		'fastp',
		'bwa_decontaminate',
		'preiVar_sars_alignmt',
		'postiVar_sars_alignmt',
		'freyja_variant_outputs',
		'freya_final_output_exists',
		'freya_output_meets_mincriteria'
		]
	
	df = pd.DataFrame(columns=column_names)

	#With 1 cpu needed per sample to be analyzed in parallel

	dfs_list = Manager().list()
	pool = Pool(processes=nb_sim_process)
	samples = lst_samples
	res = pool.map_async(partial(file_check, dfs_list), samples)
	res.wait()
	dfs = pd.concat(dfs_list, ignore_index = True)
	dfs.to_csv(filename, sep='\t')