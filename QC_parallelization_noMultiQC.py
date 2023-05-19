#!/bin/python3.6
#Version QC_parallelization_V4
#Script for parrallelization of QC steps on WW analysis
import sys
import time
import multiprocessing as mp
import os
#import subprocess
import argparse
from argparse import ArgumentParser

#First iteration of QC pipeline to be run after analysis
#Based largely on https://nf-co.re/viralrecon

#In version 2 of QC_parallelization 
# Corrected some things, inlcuding kraken2 for final report
# For now I will be using the Kraken2 database used by GenPipes (C3G)

#In version 3 of QC_parallelization
# Addressed issue with picard, which created a new one. If only small sample number it should work
# Changed the samtools depth (lines 46/53) to bedtools coverage
# Calculated mean coverage per amplicon/ORF
# Changed the way samples are handled in parallel, can now be modified with option -n (default 8)
# Cleaned up parameters with argparse
# Generates a quantative report file (calling QC-data-table.script-V1.py)

#In version 3.1
# Naming output of already analyzed files option
# Removes Picard step
# Updates the input for QC-table script
#In version 4
#Adds parallelization to QC-data-table-script-V2.py

#time feedback
start_time = time.time()

#A function that execute the recurring step of this pipeline 
def run_QC(wp_path, current_sample, nb_t_profiling, smpl_path, output, input_folder, kraken_db, analyzed_list, insertbed):
	# #	Step 1: Raw Read Quality via fastqc 0.11.9
	os.system("fastqc -q --threads {1} {2}{0}_R1_001.fastq.gz {2}{0}_R2_001.fastq.gz --outdir {3}qc_results/".format(current_sample, nb_t_profiling, smpl_path, output))
	# Generates a fastqc.zip and fastqc.html file per each sequence (forward and reverse)
	#	Step 2: If fastp files do not exist already from Freyja parallel step
			# Read Clean Status via fastp 0.23.2 (move the JSON file to qc_results)
	if not os.path.exists("{0}{1}.fastp.json".format(input_folder, current_sample)):
		os.system("fastp -i {2}{0}_R1_001.fastq.gz -I {2}{0}_R2_001.fastq.gz -o {3}{0}_R1_001_trimmed_1.fastq.gz -O {3}{0}_R2_001_trimmed_2.fastq.gz -l 70 -x --cut_tail --cut_tail_mean_quality 20 --detect_adapter_for_pe --thread {1} --json {3}qc_results/{0}.fastp.json".format(current_sample, nb_t_profiling, smpl_path, output))
		    # Else move it from the Freyja parallel output location to avoid duplicate files
	else:
		os.system("mv {0}{1}.fastp.json {2}qc_results/{1}.fastp.json".format(input_folder, current_sample,output))
	#	Step 2b: Check Read Clean Status post fastp with fastqc
	os.system("fastqc -q --threads {1} {4}{0}_R1_001_trimmed_1.fastq.gz {4}{0}_R2_001_trimmed_2.fastq.gz --outdir {2}qc_results/".format(current_sample, nb_t_profiling, output, wp_path, input_folder))
	#	Step 3: Check for contaminants (initial files vs cleaned slides)
	# 	Note: Step 3 takes a lot of RAM so I might want to run this step seperately
	#	Step 3a: Initial files
	os.system("kraken2 --db {4} --threads {1} --report {3}qc_results/{0}.initial.kraken2.report.txt --paired {2}{0}_R1_001.fastq.gz {2}{0}_R2_001.fastq.gz --output {3}qc_results/{0}initial.kraken2.output.txt".format(current_sample, nb_t_profiling, smpl_path, output, kraken_db))
	# 	Step 3b: Files cleaned (human reads removed) and used in analysis
	os.system("kraken2 --db {4} --threads {1} --report {2}qc_results/{0}.final.kraken2.report.txt --paired {3}{0}_R1_001_decon_1.fastq.gz {3}{0}_R2_001_decon_2.fastq.gz --output {3}qc_results/{0}initial.kraken2.output.txt".format(current_sample, nb_t_profiling, output, input_folder,kraken_db))
	#	Step 4: Alignment-level QC (both with bwa of human decon reads and iVar)
	#	Step4a: Alingment-level QC on alignment of human decontaminated reads to reference genome (preiVar)
	os.system("samtools flagstat --threads {2} {3}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.flagstat".format(output, current_sample, nb_t_profiling,input_folder))
	os.system("samtools idxstats {3}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.idxstats".format(output, current_sample, nb_t_profiling, input_folder))
	os.system("samtools stats --threads {2} --ref-seq {3}MN908947_3.fa {4}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.stats".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("bedtools coverage -d -a {3}genome.bed -b {4}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.depth.bed".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("bedtools coverage -mean -a {3}{5} -b {4}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.amplicon.depth.bed".format(output, current_sample, nb_t_profiling, wp_path, input_folder,insertbed))
	os.system("bedtools coverage -mean -a {3}SARS-CoV-2_ORF_full.bed -b {4}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.ORF.depth.bed".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	#os.system("picard -Xmx5g CollectMultipleMetrics --INPUT {4}{1}_preprocessed_sorted.bam --OUTPUT {0}qc_results/{1}.preiVar.CollectMultipleMetrics --REFERENCE_SEQUENCE {3}MN908947_3.fa".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	#	Step4b: Alingment-level QC on alignment of human decontaminated reads to reference genome (postiVar)
	os.system("samtools flagstat --threads {2} {3}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.flagstat".format(output, current_sample, nb_t_profiling, input_folder))
	os.system("samtools idxstats {3}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.idxstats".format(output, current_sample, nb_t_profiling, input_folder))
	os.system("samtools stats --threads {2} --ref-seq {3}MN908947_3.fa {4}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.stats".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("bedtools coverage -d -a {3}genome.bed -b {4}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.depth.bed".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("bedtools coverage -mean -a {3}{5} -b {4}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.amplicon.depth.bed".format(output, current_sample, nb_t_profiling, wp_path, input_folder,insertbed))
	os.system("bedtools coverage -mean -a {3}SARS-CoV-2_ORF_full.bed -b {4}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.ORF.depth.bed".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	#os.system("picard -Xmx5g CollectMultipleMetrics --INPUT {4}{1}_ivartrim_sorted.bam --OUTPUT {0}qc_results/{1}.postiVar.CollectMultipleMetrics --REFERENCE_SEQUENCE {3}MN908947_3.fa".format(output, current_sample, nb_t_profiling, wp_path, input_folder))
	#Removed picard as it breaks MultiQC by making it too large!
	if (args.single):
		os.system("multiqc {0}qc_results/{1}* -o {0}qc_results/{1}_multiqc_data".format(output, current_sample))
	#	Final Step: check for completion
	print("echo {1} >> {0}{2}".format(wp_path, current_sample, analyzed_list))
	os.system("echo {1} >> {0}{2}".format(wp_path, current_sample, analyzed_list))  

#create a function that test if a string is numeric
def is_numeric(the_str):
	try:
		float_conv = float(the_str)
		return True
	except (ValueError, TypeError):
		return False

if __name__ == "__main__":
		
	#Initializing parameters and options
	parser = argparse.ArgumentParser(prog= 'Quality Control Analysis on results from Wastewater bioinformatics pipeline',
				  description= "This script does a few things, one it outputs a QC-report file, two it generates a series of files that can be visualized using MultiQC. Helpful for understanding how the run worked but not mandatory")
	parser.add_argument("workspace_path",
			 help= "The path to directory where: reference files, databases, sample list, usherbarcodes are.") #Required
	parser.add_argument("samples_list_file",
			 help = "Text file with samples to be run. Used to identify input and output files. The pattern is <sample_name>_R[1-2]_001.fastq.gz") #Required
	parser.add_argument("-s", "--sample_path", default = '.',
			 help = "Default: Current working directory. The location of raw FASTQ files") #Optional
	parser.add_argument("-i", "--input_path", default = '.',
			 help = "Default: Current working directory. The ouput directory location of previous Freyja pipeline") #Optional
	parser.add_argument("-t", "--threads", default = 4,
			 help= "Default: 4 Threads available should be > threads multiplied by samples run in parallel") #Optional
	parser.add_argument("-o", "--output", default= '.',
			 help="Default: Current working directory. Will make directory specified if not found") #Optional
	parser.add_argument("-n", "--parallel", default = 4,
			 help = "Default: 4 Number of samples run in parallel.") #Optional
	parser.add_argument("-k", "--kraken", default = "Kraken2wViruses",
			 help = "Default: Kraken2wViruses, Directory of Kraken database. If running on Digital Alliance servers see option -mugqic") #Optional
	#parser.add_argument("-m", "--mugqic", action="store_true",
			 #help = "Only for Digital Alliance Users: Use the mugqic kraken2 database version of files") #Optional
	parser.add_argument("-x", "--single", action="store_true",
			 help = "When you have too many files for one MultiQC report run each seperately") #Optional
	parser.add_argument("-a", "--analyzed", default="QC_analyzed_samples.txt",
			 help = "Output file for listing samples that have been completed") #Optional
	parser.add_argument("-V", "--arctic", default = "V4.1",
		     help = "Specify the Arctic primer scheme used: V1, V2, V3, V4, V4.1") #Optional
	args = parser.parse_args()

	#### ARGUMENTS && MISC FILE CHECKING ####
	
	#Workspace (Path) [required]
	#Workspace is the directory where all files should be located (references, databases, scripts, and barcodes)
	#This is the default location of all the files unless stated
	workspace_path = args.workspace_path 
	#Making sure to have "/" at the end of the workspace absolute path
	if ((workspace_path[len(workspace_path)-1]) != "/"):
		workspace_path = "{0}/".format(workspace_path)
	if not os.path.isdir(workspace_path):
		print("Workspace directory does not exist")
		raise SystemExit(1)
	print("Chosen Workspace is : '{0}'".format(workspace_path))
	
	#FASTQ file directory of samples [optional]: Default current directory
	#If sequence data is stored in a different location
	#Needed to do QC on the sequencing.
	sample_path = args.sample_path  #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((sample_path[len(sample_path)-1]) != "/"):
		sample_path = "{0}/".format(sample_path)
	if not os.path.isdir(sample_path):
		print("Sample directory does not exist")
		raise SystemExit(1)
	print("Sample Path is : '{0}'".format(sample_path))

	#Sample list file (same as Freyja file) should be a simple text file with each sample on a new line and located in the workspace path
	#We are only using paired reads at this point so there should be a forward and reverse read follow the pattern
	# sample + _R2_001.fastq.gz and sample + _R1_001.fastq.gz
	samples_list_file = args.samples_list_file
	print("Looking for sample list and paired read FASTQ files")
	if not os.path.exists("{0}{1}".format(workspace_path,samples_list_file)):
		print("Sample list file could not be found in chosen workspace directory " + workspace_path)
		raise SystemExit(1)
	lst_samples=[]
	with open("{0}{1}".format(workspace_path,samples_list_file)) as f:
		lst_samples = f.readlines()
	the_idx = 0
	for sample_name in lst_samples:
		lst_samples[the_idx] = sample_name.replace("\n", "")
		the_idx = the_idx + 1
	num_samples = len(lst_samples)
	for sample_name in lst_samples:
		if not os.path.isfile("{0}{1}_R1_001.fastq.gz".format(sample_path,sample_name)):
			print("{0}{1}_R1_001.fastq.gz Forward read Fastq file could be found. Make sure files are in the specified directory using options -s or --sample_path option. Also that forward reads have the format sample + _R1_001.fastq.gz".format(sample_path,sample_name))
			raise SystemExit(1)
		if not os.path.isfile("{0}{1}_R2_001.fastq.gz".format(sample_path,sample_name)):
			print("{0}{1}_R1_001.fastq.gz Reverse read Fastq file could be found. Only handles paired reads at the moment.".format(sample_path,sample_name))
			raise SystemExit(1)
	
	#Path for directory output of pipeline
	output = args.output
	if not(args.output):
		output = os.getcwd()
	#Making sure to have "/" at the end of the workspace absolute path
	if ((output[len(output)-1]) != "/"):
		output = "{0}/".format(output)
	print("Output Path is : '{0}'".format(output))
	
	#Make directory for output called qc_results if does not exist
	if not os.path.exists(output + "qc_results"):
		os.system("mkdir -p {0}qc_results".format(output))

	#Input folder, which is the output of the Freyja pipeline (for which we will run QC on)
	input_folder = args.input_path
	if not(args.input_path):
		input_folder = os.getcwd()
	#making sure to have "/" at the end of the workspace absolute path
	if ((input_folder[len(input_folder)-1]) != "/"):
		input_folder = "{0}/".format(input_folder)
	print("Input Path is : '{0}'".format(input_folder))

	#Path to users Kraken2 database or by default Kraken2wViruses in the Workspace
	kraken_db = args.kraken
	#If the flag mugqic flag is raised then negate the Kraken2 databse and use the one provided on Digital Alliance by C3G MUGQIQ
	#if args.mugqic:
		#kraken_db = "/cvmfs/soft.mugqic/CentOS6/software/kraken2/kraken2-2.1.0/db"
	if args.kraken:
		if not os.path.exists("{0}".format(kraken_db)):
			print("The Kraken database you have listed could not be found. You listed " + kraken_db)
			raise SystemExit(1)
	else:
		if not os.path.exists("{0}".format(kraken_db)):
			print("The default Kraken database has not been built/found:" + kraken_db + "If running analysis on Digital Alliance server try the mugqic option")
			raise SystemExit(1)
	
	#Making that Artic insert bed file is downloaded
	if not os.path.exists("{0}SARS-CoV-2.{1}.insert.bed".format(workspace_path, args.arctic)):
		os.system("wget -O {0}SARS-CoV-2.{1}.insert.bed https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/{1}/SARS-CoV-2.insert.bed".format(workspace_path, args.arctic))
	insertbed = "SARS-CoV-2.{0}.insert.bed".format(args.arctic)

	#Number of samples being analysed in parallel (Default 4 samples at a time)
	nb_sim_process = int(args.parallel)
	print('Running ' + str(nb_sim_process) + ' of samples simultaneously')
	
	#Number of threads to use for multithread tasks (fastp and bwa): Default 4
	threads = int(args.threads)
	print('Using ' + str(threads) + 'threads per ' + str(nb_sim_process) + ' samples analyzed simultaneously')
	
	#Theoretically the number of CPUs available should be threads * nb_sim_process which is what we can request for QC table script
	total_CPU = int(nb_sim_process) * int(threads)

	#Ceate file that records samples that have already been analyzed for QC
	analyzed_list=args.analyzed
	os.system("/bin/touch {0}{1}".format(workspace_path, analyzed_list))

	print("File checking complete, running analysis")

	
	#Analyze samples in parrallel (4 by default or whichever number selected by --parallel)
	liste_index_sample = range(0, num_samples, nb_sim_process)
	for i in liste_index_sample:
		# create a list that will record the samples that have already been analysed
		already_analysed_samples_list=[]
		with open("{0}{1}".format(workspace_path, analyzed_list)) as f:
			already_analyzed_samples_list = f.readlines()
		index_line = 0
		for line in already_analyzed_samples_list:
			already_analyzed_samples_list[index_line] = line.replace("\n", "")
			index_line = index_line + 1

		processes = []
		if (i + nb_sim_process) < num_samples:
			sub_list_index = range(i, (i+nb_sim_process), 1)
			for t in sub_list_index:
				if (not(lst_samples[t] in already_analyzed_samples_list)):
					print("Starting QC-analysis on sample : " + lst_samples[t])
					globals()["p"+str(t)] = mp.Process(target=run_QC, args=(workspace_path, lst_samples[t], threads, sample_path, output, input_folder, kraken_db, analyzed_list, insertbed,))
					processes.append(globals()["p"+str(t)])
					globals()["p"+str(t)].start()
				else:
					print(lst_samples[t] + "has already been analyzed, if re-analyzing, remove from " + analyzed_list)
		else:
			sub_list_index = range(i,num_samples,1)
			for t in sub_list_index:
				if (not(lst_samples[t] in already_analyzed_samples_list)):
					print("Starting QC-analysis on sample : " + lst_samples[t])
					globals()["p"+str(t)] = mp.Process(target=run_QC, args=(workspace_path, lst_samples[t], threads, sample_path, output, input_folder, kraken_db,analyzed_list, insertbed))
					processes.append(globals()["p"+str(t)])
					globals()["p"+str(t)].start()
				else:
					print(lst_samples[t] + "has already been analyzed, if re-analyzing, remove from " + analyzed_list)
		if (len(processes) != 0):
			for p in processes:
				p.join()
			for p in processes:
				while (p.is_alive()):
					time.sleep(5)

	#Step 6: Generate the report with MultiQC for ALL files if not flagged for single-files
	#if not (args.single):
		#os.system("multiqc {0}qc_results/ -o {0}qc_results/multiqc_data".format(output))

	#Step 7: Generate a quantative report
	os.system("python3 {0}QC-data-table-script.py {0} {1} {2} {3}qc_results/ {3}qc_results/ {4} {5}".format(workspace_path, input_folder, samples_list_file, output, insertbed, total_CPU))
	
	
