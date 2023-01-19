#!/bin/python3.6

#Script for parrallelization of QC steps on WW analysis
import sys
import time
import multiprocessing as mp
import os
#import subprocess

#First iteration of QC pipeline to be run after analysis
#Based largely on https://nf-co.re/viralrecon

#In version 2 of QC_parallelization 
# Corrected some things, inlcuding kraken2 for final report
# For now I will be using the Kraken2 database used by GenPipes (C3G)


#time feedback
start_time = time.time()

#A function that execute the recurring step of this pipeline 
def run_QC(wp_path, current_sample, nb_t_profiling, smpl_path, output_folder, input_folder):
	# #	Step 1: Raw Read Quality via fastqc 0.11.9
	os.system("fastqc -q --threads {1} {2}{0}_R1_001.fastq.gz {2}{0}_R2_001.fastq.gz --outdir {3}qc_results/".format(current_sample, nb_t_profiling, smpl_path, output_folder))
	# Generates a fastqc.zip and fastqc.html file per each sequence (forward and reverse)
	#	Step 2: Read Clean Status via fastp 0.23.2 (move the JSON file to qc_results)
	os.system("fastp -i {2}{0}_R1_001.fastq.gz -I {2}{0}_R2_001.fastq.gz -o {3}{0}_R1_001_trimmed_1.fastq.gz -O {3}{0}_R2_001_trimmed_2.fastq.gz -l 70 -x --cut_tail --cut_tail_mean_quality 20 --detect_adapter_for_pe --thread {1} --json {3}{0}.fastp.json".format(current_sample, nb_t_profiling, smpl_path, output_folder))

	#	Step 2b: Check Read Clean Status post fastp with fastqc
	os.system("fastqc -q --threads {1} {4}{0}_R1_001_trimmed_1.fastq.gz {4}{0}_R2_001_trimmed_2.fastq.gz --outdir {3}qc_results/".format(current_sample, nb_t_profiling, output_folder, wp_path, input_folder))
	#	Step 3: Check for contaminants (initial files vs cleaned slides)
	# 	Note: Step 3 takes a lot of RAM so I might want to run this step seperately
    #	Step 3a: Initial files
	os.system("kraken2 --db {4}/Kraken2wViruses --threads {1} --report {3}qc_results/{0}.initial.kraken2.report.txt --paired {2}{0}_R1_001.fastq.gz {2}{0}_R2_001.fastq.gz ".format(current_sample, nb_t_profiling, smpl_path, output_folder, wp_path))
    # 	Step 3b: Files cleaned (human reads removed) and used in analysis
	os.system("kraken2 --db {4}/Kraken2wViruses --threads {1} --report {3}qc_results/{0}.final.kraken2.report.txt --paired {5}{0}_R1_001_decon_1.fastq.gz {5}{0}_R2_001_decon_2.fastq.gz".format(current_sample, nb_t_profiling, smpl_path, output_folder, wp_path, input_folder))
    #	Step 4: Alignment-level QC (both with bwa of human decon reads and iVar)
	#	Step4a: Alingment-level QC on alignment of human decontaminated reads to reference genome (preiVar)
	os.system("samtools flagstat --threads {2} {3}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.flagstat".format(output_folder, current_sample, nb_t_profiling,input_folder))
	os.system("samtools idxstats {3}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.idxstats".format(output_folder, current_sample, nb_t_profiling, input_folder))
	os.system("samtools stats --threads {2} --ref-seq {3}MN908947_3.fa {4}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.stats".format(output_folder, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("samtools depth {4}{1}_preprocessed_sorted.bam > {0}qc_results/{1}.preiVar.common_depth_report".format(output_folder, current_sample, nb_t_profiling, wp_path, input_folder))
	# Picard seems to be failing also due heavy memory use might re-run this step seperately too
	os.system("java -jar /opt/picard.jar -Xmx5g CollectMultipleMetrics --INPUT {4}{1}_preprocessed_sorted.bam --OUTPUT {0}qc_results/{1}.preiVar.CollectMultipleMetrics --REFERENCE_SEQUENCE {3}MN908947_3.fa".format(output_folder, current_sample, nb_t_profiling, wp_path, input_folder))
	#	Step4b: Alingment-level QC on alignment of human decontaminated reads to reference genome (postiVar)
	os.system("samtools flagstat --threads {2} {3}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.flagstat".format(output_folder, current_sample, nb_t_profiling, input_folder))
	os.system("samtools idxstats {3}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.idxstats".format(output_folder, current_sample, nb_t_profiling, input_folder))
	os.system("samtools stats --threads {2} --ref-seq {3}MN908947_3.fa {4}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.stats".format(output_folder, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("samtools depth {4}{1}_ivartrim_sorted.bam > {0}qc_results/{1}.postiVar.common_depth_report".format(output_folder, current_sample, nb_t_profiling, wp_path, input_folder))
	os.system("java -jar /opt/picard.jar -Xmx5g CollectMultipleMetrics --INPUT {4}{1}_ivartrim_sorted.bam --OUTPUT {0}qc_results/{1}.postiVar.CollectMultipleMetrics --REFERENCE_SEQUENCE {3}MN908947_3.fa".format(output_folder, current_sample, nb_t_profiling, wp_path, input_folder))
    
    #	Final Step: check for completion
	print("echo {1} >> {0}QC_analyzed_samples.txt".format(wp_path, current_sample))
	os.system("echo {1} >> {0}QC_analyzed_samples.txt".format(wp_path, current_sample))  

#create a function that test if a string is numeric
def is_numeric(the_str):
	try:
		float_conv = float(the_str)
		return True
	except (ValueError, TypeError):
		return False

if __name__ == "__main__":
	#Initializing global variables
	
	#Path of the Workspace (where all reference & database files are stored)
	workspace_path =sys.argv[1] 
	#Making sure to have "/" at the end of the workspace absolute path
	if ((workspace_path[len(workspace_path)-1]) != "/"):
		workspace_path = "{0}/".format(workspace_path)
	print("Chosen Workspace is : '{0}'".format(workspace_path))
	
	#Path of the sample location (Raw FASTQ files)
	sample_path =sys.argv[2] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((sample_path[len(sample_path)-1]) != "/"):
		sample_path = "{0}/".format(sample_path)
	print("Sample Path is : '{0}'".format(sample_path))

	#Variable for list of samples to be included in QC report
	samples_list_file = sys.argv[3]
	lst_samples=[]
	with open("{0}{1}".format(workspace_path,samples_list_file)) as f:
		lst_samples = f.readlines()

	the_idx = 0
	for sample_name in lst_samples:
		lst_samples[the_idx] = sample_name.replace("\n", "")
		the_idx = the_idx + 1
	#Number of samples being analysed simultaneously
	nb_sim_process = 8

	#Number of threads to use for multithread tasks (fastp and bwa)
	nb_threads = int(sys.argv[4])
	
	#Output folder where all the quality control files will go
	output_folder = sys.argv[5]
	#making sure to have "/" at the end of the workspace absolute path
	if ((output_folder[len(output_folder)-1]) != "/"):
		output_folder = "{0}/".format(output_folder)
	print("Output Path is : '{0}'".format(output_folder))

	#Input folder of anlysis for which the intermediate files are stored (for which we will run QC on)
	input_folder = sys.argv[6]
	#making sure to have "/" at the end of the workspace absolute path
	if ((input_folder[len(input_folder)-1]) != "/"):
		input_folder = "{0}/".format(input_folder)
	print("Input Path is : '{0}'".format(input_folder))
	
    #Quality control results will be stored within the output directory in a folder named 'qc_results'
	os.system("mkdir -p {0}qc_results/".format(output_folder))


	#Ceate file that records samples that have already been analyzed for QC
	os.system("/bin/touch {0}QC_analyzed_samples.txt".format(workspace_path))
	

	#execute nb_sim_process samples at a time
	liste_index_sample = range(0,len(lst_samples),nb_sim_process)
	for i in liste_index_sample:
		# create a list that will record the samples that have already been analysed
		already_analysed_samples_list=[]
		with open("{0}QC_analyzed_samples.txt".format(workspace_path)) as f:
			already_analyzed_samples_list = f.readlines()
		index_line = 0
		for line in already_analyzed_samples_list:
			already_analyzed_samples_list[index_line] = line.replace("\n", "")
			index_line = index_line + 1

		processes = []
		#define all processes that will be executed parallelly (run only those that have not been executed yet)
		if (((i)<len(lst_samples)) and (not (lst_samples[i] in already_analyzed_samples_list))):
			p1 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p1)
			p1.start()
		if (((i+1)<len(lst_samples)) and (not (lst_samples[i+1] in already_analyzed_samples_list))):
			p2 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+1], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p2)
			p2.start()
		if (((i+2)<len(lst_samples)) and (not (lst_samples[i+2] in already_analyzed_samples_list))):
			p3 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+2], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p3)
			p3.start()
		if (((i+3)<len(lst_samples)) and (not (lst_samples[i+3] in already_analyzed_samples_list))):
			p4 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+3], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p4)
			p4.start()
		if (((i+4)<len(lst_samples)) and (not (lst_samples[i+4] in already_analyzed_samples_list))):
			p5 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+4], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p5)
			p5.start()
		if (((i+5)<len(lst_samples)) and (not (lst_samples[i+5] in already_analyzed_samples_list))):
			p6 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+5], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p6)
			p6.start()
		if (((i+6)<len(lst_samples)) and (not (lst_samples[i+6] in already_analyzed_samples_list))):
			p7 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+6], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p7)
			p7.start()
		if (((i+7)<len(lst_samples)) and (not (lst_samples[i+7] in already_analyzed_samples_list))):
			p8 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+7], nb_threads,sample_path, output_folder,input_folder,))
			processes.append(p8)
			p8.start()
		#if (((i+8)<len(lst_samples)) and (not (lst_samples[i+8] in already_analyzed_samples_list))):
		 	#p9 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+8], nb_threads,))
		 	#processes.append(p9)
		 	#p9.start()
		#if (((i+9)<len(lst_samples)) and (not (lst_samples[i+9] in already_analyzed_samples_list))):
		 	#p10 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+9], nb_threads,))
		 	#processes.append(p10)
		 	#p10.start()
		#if (((i+10)<len(lst_samples)) and (not (lst_samples[i+10] in already_analyzed_samples_list))):
		 	#p11 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+10], nb_threads,))
		 	#processes.append(p11)
		 	#p11.start()
		#if (((i+11)<len(lst_samples)) and (not (lst_samples[i+11] in already_analyzed_samples_list))):
		 	#p12 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+11], nb_threads,))
		 	#processes.append(p12)
			#p12.start()
		#if (((i+12)<len(lst_samples)) and (not (lst_samples[i+12] in already_analyzed_samples_list))):
		 	#p13 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+12], nb_threads,))
		 	#processes.append(p13)
			#p13.start()
		#if (((i+13)<len(lst_samples)) and (not (lst_samples[i+13] in already_analyzed_samples_list))):
		 	#p14 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+13], nb_threads,))
		 	#processes.append(p14)
		 	#p14.start()
		#if (((i+14)<len(lst_samples)) and (not (lst_samples[i+14] in already_analyzed_samples_list))):
		 	#p15 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+14], nb_threads,))
		 	#processes.append(p15)
		 	#p15.start()
		#if (((i+15)<len(lst_samples)) and (not (lst_samples[i+15] in already_analyzed_samples_list))):
		 	#p16 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+15], nb_threads,))
		 	#processes.append(p16)
		 	#p16.start()
		# if (((i+16)<len(lst_samples)) and (not (lst_samples[i+16] in already_analyzed_samples_list))):
		# 	p17 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+16], nb_threads,))
		# 	processes.append(p17)
		# 	p17.start()
		# if (((i+17)<len(lst_samples)) and (not (lst_samples[i+17] in already_analyzed_samples_list))):
		# 	p18 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+17], nb_threads,))
		# 	processes.append(p18)
		# 	p18.start()
		# if (((i+18)<len(lst_samples)) and (not (lst_samples[i+18] in already_analyzed_samples_list))):
		# 	p19 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+18], nb_threads,))
		# 	processes.append(p19)
		# 	p19.start()
		# if (((i+19)<len(lst_samples)) and (not (lst_samples[i+19] in already_analyzed_samples_list))):
		# 	p20 = mp.Process(target=run_QC, args=(workspace_path, lst_samples[i+19], nb_threads,))
		# 	processes.append(p20)
		# 	p20.start()
		#wait for all processes to finish
		if (len(processes) != 0):
			for p in processes:
				p.join()
			for p in processes:
				while (p.is_alive()):
					time.sleep(5)
    
    #Step 6: Generate the report with MultiQC
	os.system("multiqc {0}qc_results/ -o {0}qc_results/multiqc_data".format(output_folder))
	
	