#!/bin/python3.6

#Script for parrallelization of Freyja
import sys
import time
import multiprocessing as mp
import os
import subprocess


#time feedback
start_time = time.time()

#A function that execute the recurring step of this pipeline 
def run_core_iPMVC(wp_path,current_sample, nb_t_profiling, smpl_path, usherbarcodes, output_folder):
	# Quality Control: trimming, automatic adapters removal and low complexity reads removal
	os.system("fastp -i {2}{0}_R1_001.fastq.gz -I {2}{0}_R2_001.fastq.gz -o {3}{0}_R1_001_trimmed_1.fastq -O {3}{0}_R2_001_trimmed_2.fastq -l 70 -x --cut_tail --cut_tail_mean_quality 20 --detect_adapter_for_pe --thread {1}".format(current_sample, nb_t_profiling, smpl_path, output_folder))
	
	#Decontaminate human reads 
	os.system("bwa mem {3}Homo_sapiens.GRCh38.fa {0}{1}_R1_001_trimmed_1.fastq {0}{1}_R2_001_trimmed_2.fastq -t {2} > {0}{1}.sam".format(output_folder, current_sample, nb_t_profiling, wp_path))
	os.system("samtools view -bS {0}{1}.sam > {0}{1}.bam ".format(output_folder, current_sample))
	os.system("samtools view -b -f 12 -F 256 {0}{1}.bam > {0}{1}_human_decon.bam ".format(output_folder, current_sample))
	os.system("samtools sort -n -m 5G -@ {2} {0}{1}_human_decon.bam -o {0}{1}_human_decon_sorted.bam ".format(output_folder, current_sample, nb_t_profiling))
	os.system("samtools fastq -@ {2} {0}{1}_human_decon_sorted.bam -1 {0}{1}_R1_001_decon_1.fastq -2 {0}{1}_R2_001_decon_2.fastq".format(output_folder, current_sample, nb_t_profiling))
		
	#Align reads to reference genome 
	os.system("bwa mem {0}MN908947_3.fa {3}{1}_R1_001_decon_1.fastq {3}{1}_R2_001_decon_2.fastq -t {2} > {3}{1}_preprocessed.sam".format(wp_path, current_sample, nb_t_profiling, output_folder))
	
	#Compress, sort and filter alignment
	
	os.system("samtools view -bS {0}{1}_preprocessed.sam > {0}{1}_preprocessed.bam".format(output_folder, current_sample))
	os.system("samtools sort {0}{1}_preprocessed.bam -o {0}{1}_preprocessed_sorted.bam".format(output_folder, current_sample))
	
	os.system("ivar trim -i {2}{1}_preprocessed_sorted.bam -b {0}SARS-CoV-2.primer.bed -p {2}{1}_ivartrim".format(wp_path, current_sample,output_folder))
	# Freyja commnad
	os.system("samtools sort -o {0}{1}_ivartrim_sorted.bam {0}{1}_ivartrim.bam".format(output_folder, current_sample))

	os.system("freyja variants {2}{1}_ivartrim_sorted.bam --variants {2}{1}_variantout --depths {2}{1}_depthout --ref {0}MN908947_3.fasta".format(wp_path, current_sample, output_folder))

	os.system("freyja demix {3}{1}_variantout.tsv {3}{1}_depthout --barcodes {0}{2} --output {3}results/{1}_output".format(wp_path, current_sample, usherbarcodes, output_folder))
	print("echo {1} >> {0}Already_analyzed_samples.txt".format(wp_path, current_sample))
	os.system("echo {1} >> {0}Already_analyzed_samples.txt".format(wp_path, current_sample))  

#create a function that test if a string is numeric
def is_numeric(the_str):
	try:
		float_conv = float(the_str)
		return True
	except (ValueError, TypeError):
		return False

if __name__ == "__main__":
	#Initializing global variables
	
	#Path of the Workspace where there are all analysis will go
	workspace_path =sys.argv[1] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((workspace_path[len(workspace_path)-1]) != "/"):
		workspace_path = "{0}/".format(workspace_path)
	print("Chosen Workspace is : '{0}'".format(workspace_path))
	
	#Path of the sample location
	sample_path =sys.argv[2] #commandline argument
	#making sure to have "/" at the end of the workspace absolute path
	if ((sample_path[len(sample_path)-1]) != "/"):
		sample_path = "{0}/".format(sample_path)
	print("Sample Path is : '{0}'".format(sample_path))

	#Create the variable for list of samples 
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

	#Check for most up to date usher barcodes
	usherbarcodes = sys.argv[5]
	down_date = subprocess.check_output("curl http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt | awk -F'[()]' '{print $2}'", shell=True)
	down_date = down_date.decode("utf-8").split("\n")
	down_date = down_date[0]
	down_date =  "usher_barcodes_{0}.csv".format(down_date)
	print("Usher barcode used {0}".format(usherbarcodes))
	print("Lastest barcode should be {0}".format(down_date))
	if usherbarcodes != down_date:
		print("Updating Usherbarcordes")
		print("Output new usher barcodes in {0}".format(workspace_path))
		os.system("python3 {0}/usher_update_V2.py {0}".format(workspace_path))
		usherbarcodes = down_date
	
	#Use this for naming the output folder
	output_folder = sys.argv[6]
	#making sure to have "/" at the end of the workspace absolute path
	if ((output_folder[len(output_folder)-1]) != "/"):
		output_folder = "{0}/".format(output_folder)
	print("Output Path is : '{0}'".format(output_folder))
	
	#Make directory output if does not exist
	if not os.path.exists(output_folder):
		os.system("mkdir -p {0}results/".format(output_folder))

	#Making sure that Human refrence genome BWA index is downloaded
	if not os.path.exists("{0}Homo_sapiens.GRCh38.fa.bwt".format(workspace_path)):
		os.system("wget --quiet -O {0}Homo_sapiens.GRCh38.dna.alt.fa.gz http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz".format(workspace_path))
		os.system("gzip -d {0}Homo_sapiens.GRCh38.dna.alt.fa.gz".format(workspace_path))
		os.system("mv {0}Homo_sapiens.GRCh38.dna.alt.fa.gz {0}Homo_sapiens.GRCh38.fa".format(workspace_path))
		os.system("bwa index -p {0}Homo_sapiens.GRCh38.fa {0}Homo_sapiens.GRCh38.fa".format(workspace_path))

	#Making sure that Wuhan reference SARS genome BWA index is downloaded
	if not os.path.exists("{0}MN908947_3.fa.bwt".format(workspace_path)):
		os.system("wget --quiet -O {0}MN908947_3.fa https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true".format(workspace_path))
		os.system("bwa index -p {0}MN908947_3.fa {0}MN908947_3.fa".format(workspace_path))
	

	#Making that Artic V4.1 bed file is downloaded
	if not os.path.exists("{0}SARS-CoV-2.primer.bed".format(workspace_path)):
		os.system("wget -O {0}SARS-CoV-2.primer.bed https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V4.1/SARS-CoV-2.primer.bed".format(workspace_path))

	#create file that records samples that have already been analyzed
	#os.system("/bin/rm -f {0}Already_analyzed_samples.txt".format(workspace_path))
	os.system("/bin/touch {0}Already_analyzed_samples.txt".format(workspace_path))
	

	#execute nb_sim_process samples at a time
	liste_index_sample = range(0,len(lst_samples),nb_sim_process)
	for i in liste_index_sample:
		# create a list that will record the samples that have already been analysed
		already_analysed_samples_list=[]
		with open("{0}Already_analyzed_samples.txt".format(workspace_path)) as f:
			already_analyzed_samples_list = f.readlines()
		index_line = 0
		for line in already_analyzed_samples_list:
			already_analyzed_samples_list[index_line] = line.replace("\n", "")
			index_line = index_line + 1

		processes = []
		#define all processes that will be executed parallelly (run only those that have not been executed yet)
		if (((i)<len(lst_samples)) and (not (lst_samples[i] in already_analyzed_samples_list))):
			p1 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p1)
			p1.start()
		if (((i+1)<len(lst_samples)) and (not (lst_samples[i+1] in already_analyzed_samples_list))):
			p2 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+1], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p2)
			p2.start()
		if (((i+2)<len(lst_samples)) and (not (lst_samples[i+2] in already_analyzed_samples_list))):
			p3 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+2], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p3)
			p3.start()
		if (((i+3)<len(lst_samples)) and (not (lst_samples[i+3] in already_analyzed_samples_list))):
			p4 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+3], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p4)
			p4.start()
		if (((i+4)<len(lst_samples)) and (not (lst_samples[i+4] in already_analyzed_samples_list))):
			p5 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+4], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p5)
			p5.start()
		if (((i+5)<len(lst_samples)) and (not (lst_samples[i+5] in already_analyzed_samples_list))):
			p6 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+5], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p6)
			p6.start()
		if (((i+6)<len(lst_samples)) and (not (lst_samples[i+6] in already_analyzed_samples_list))):
			p7 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+6], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p7)
			p7.start()
		if (((i+7)<len(lst_samples)) and (not (lst_samples[i+7] in already_analyzed_samples_list))):
			p8 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+7], nb_threads,sample_path,usherbarcodes, output_folder,))
			processes.append(p8)
			p8.start()
		#if (((i+8)<len(lst_samples)) and (not (lst_samples[i+8] in already_analyzed_samples_list))):
		 	#p9 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+8], nb_threads,))
		 	#processes.append(p9)
		 	#p9.start()
		#if (((i+9)<len(lst_samples)) and (not (lst_samples[i+9] in already_analyzed_samples_list))):
		 	#p10 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+9], nb_threads,))
		 	#processes.append(p10)
		 	#p10.start()
		#if (((i+10)<len(lst_samples)) and (not (lst_samples[i+10] in already_analyzed_samples_list))):
		 	#p11 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+10], nb_threads,))
		 	#processes.append(p11)
		 	#p11.start()
		#if (((i+11)<len(lst_samples)) and (not (lst_samples[i+11] in already_analyzed_samples_list))):
		 	#p12 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+11], nb_threads,))
		 	#processes.append(p12)
			#p12.start()
		#if (((i+12)<len(lst_samples)) and (not (lst_samples[i+12] in already_analyzed_samples_list))):
		 	#p13 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+12], nb_threads,))
		 	#processes.append(p13)
			#p13.start()
		#if (((i+13)<len(lst_samples)) and (not (lst_samples[i+13] in already_analyzed_samples_list))):
		 	#p14 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+13], nb_threads,))
		 	#processes.append(p14)
		 	#p14.start()
		#if (((i+14)<len(lst_samples)) and (not (lst_samples[i+14] in already_analyzed_samples_list))):
		 	#p15 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+14], nb_threads,))
		 	#processes.append(p15)
		 	#p15.start()
		#if (((i+15)<len(lst_samples)) and (not (lst_samples[i+15] in already_analyzed_samples_list))):
		 	#p16 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+15], nb_threads,))
		 	#processes.append(p16)
		 	#p16.start()
		# if (((i+16)<len(lst_samples)) and (not (lst_samples[i+16] in already_analyzed_samples_list))):
		# 	p17 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+16], nb_threads,))
		# 	processes.append(p17)
		# 	p17.start()
		# if (((i+17)<len(lst_samples)) and (not (lst_samples[i+17] in already_analyzed_samples_list))):
		# 	p18 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+17], nb_threads,))
		# 	processes.append(p18)
		# 	p18.start()
		# if (((i+18)<len(lst_samples)) and (not (lst_samples[i+18] in already_analyzed_samples_list))):
		# 	p19 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+18], nb_threads,))
		# 	processes.append(p19)
		# 	p19.start()
		# if (((i+19)<len(lst_samples)) and (not (lst_samples[i+19] in already_analyzed_samples_list))):
		# 	p20 = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[i+19], nb_threads,))
		# 	processes.append(p20)
		# 	p20.start()
		#wait for all processes to finish
		if (len(processes) != 0):
			for p in processes:
				p.join()
			for p in processes:
				while (p.is_alive()):
					time.sleep(5)
