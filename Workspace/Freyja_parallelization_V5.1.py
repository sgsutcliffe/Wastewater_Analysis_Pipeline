#!/bin/python3.6

#Script for analysis of FASTQ files via Freyja
import sys
import time
import multiprocessing as mp
import os
from datetime import date
from datetime import timedelta
import argparse
from argparse import ArgumentParser

#Updates in Version 4
#Automatically remove all *sam files as they become redundant after converting to bam
#Adds in the '-update' or '-fupdate' to decide whether the user would like to update the barcodes

#Updates in Version 5
# Convert iVar output to a snpEff annotated VCF
# Changed the way samples are handled in parallel, can now be modified with option -n (default 8)
# Cleaned up parameters with argparse
# Integrated File Checking


#Time feedback
start_time = time.time()

#A function that execute the recurring step of this pipeline 
def run_core_iPMVC(wp_path,current_sample, nb_t_profiling, smpl_path, usherbarcodes, output):
	# Quality Control: trimming, automatic adapters removal and low complexity reads removal
	os.system("fastp -i {2}{0}_R1_001.fastq.gz -I {2}{0}_R2_001.fastq.gz -o {3}{0}_R1_001_trimmed_1.fastq.gz -O {3}{0}_R2_001_trimmed_2.fastq.gz -l 70 -x --cut_tail --cut_tail_mean_quality 20 --detect_adapter_for_pe --thread {1} --json {3}{0}.fastp.json".format(current_sample, nb_t_profiling, smpl_path, output))
	
	#Decontaminate human reads 
	os.system("bwa mem {3}Homo_sapiens.GRCh38.fa {0}{1}_R1_001_trimmed_1.fastq.gz {0}{1}_R2_001_trimmed_2.fastq.gz -t {2} > {0}{1}.sam".format(output, current_sample, nb_t_profiling, wp_path))
	#Replace with kraken
	os.system("samtools view -bS {0}{1}.sam > {0}{1}.bam && rm {0}{1}.sam ".format(output, current_sample))
	os.system("samtools view -b -f 12 -F 256 {0}{1}.bam > {0}{1}_human_decon.bam ".format(output, current_sample))
	os.system("samtools sort -n -m 5G -@ {2} {0}{1}_human_decon.bam -o {0}{1}_human_decon_sorted.bam ".format(output, current_sample, nb_t_profiling))
	os.system("samtools fastq -@ {2} {0}{1}_human_decon_sorted.bam -1 {0}{1}_R1_001_decon_1.fastq.gz -2 {0}{1}_R2_001_decon_2.fastq.gz".format(output, current_sample, nb_t_profiling))
		
	#Align reads to reference genome 
	os.system("bwa mem {0}MN908947_3.fa {3}{1}_R1_001_decon_1.fastq.gz {3}{1}_R2_001_decon_2.fastq.gz -t {2} > {3}{1}_preprocessed.sam".format(wp_path, current_sample, nb_t_profiling, output))
	
	#Compress, sort and filter alignment
	
	os.system("samtools view -bS {0}{1}_preprocessed.sam > {0}{1}_preprocessed.bam && rm {0}{1}_preprocessed.sam".format(output, current_sample))
	os.system("samtools sort {0}{1}_preprocessed.bam -o {0}{1}_preprocessed_sorted.bam".format(output, current_sample))
	
	os.system("ivar trim -i {2}{1}_preprocessed_sorted.bam -b {0}SARS-CoV-2.primer.bed -p {2}{1}_ivartrim".format(wp_path, current_sample,output))
	# Freyja commnad
	os.system("samtools sort -o {0}{1}_ivartrim_sorted.bam {0}{1}_ivartrim.bam".format(output, current_sample))
	os.system("samtools index {0}{1}_ivartrim_sorted.bam".format(output, current_sample))
	os.system("freyja variants {2}{1}_ivartrim_sorted.bam --variants {2}{1}_variantout --depths {2}{1}_depthout --ref {0}MN908947_3.fa".format(wp_path, current_sample, output))
	os.system("freyja demix {3}{1}_variantout.tsv {3}{1}_depthout --barcodes {0}{2} --output {3}results/{1}_output".format(wp_path, current_sample, usherbarcodes, output))
	# Annotate the iVar variant file
	#Convert to VCF
	os.system("python3 {0}ivar_variants_to_vcf.py --pass_only {1}{2}_variantout.tsv {1}{2}_variantout.vcf".format(wp_path, output, current_sample))
	#Annotate with snpEff
	os.system("java -Xmx8g -jar /opt/snpEff/snpEff.jar MN908947.3 {0}{1}_variantout.vcf > {0}{1}_variantout.snpEff.annotated.vcf".format(output, current_sample))
	#Make a text readible output file
	os.system("python3 {0}convert_recurrent_snpeff_data.py -f {1}{2}_variantout.snpEff.annotated.vcf -s {2} -o {1}{2}_snpEff_converted".format(wp_path, output, current_sample))

	if os.path.exists("{1}results/{0}_output".format(current_sample,output)):
		#glob.glob("date | date '+%F' >> {1}results/{0}_output".format(current_sample, output))
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
	
	#Initializing parameters and options
	parser = argparse.ArgumentParser(prog= 'Wastewater bioinformatics pipeline',
				  description= "This is a script that takes FASTQ files from tiled amplicon-sequencing of wastewater and then runs Freyja to assign lineage/relative abundance measurements, as well as identify mutations in the sample")
	parser.add_argument("workspace_path",
		     help= "The path to directory where: reference files, databases, sample list, usherbarcodes are.") #Required
	parser.add_argument("samples_list_file",
		     help = "Text file with samples to be run. Used to identify input and output files. The pattern is <sample_name>_R[1-2]_001.fastq.gz") #Required
	parser.add_argument("-s", "--sample_path", default = '.',
		     help = "Default: Current working directory.") #Optional
	parser.add_argument("-t", "--threads", default = 4,
		     help= "Default: 4 Threads available should be > threads multiplied by samples run in parallel") #Optional
	parser.add_argument("-b", "--barcode",
		     help = "Default: Will download/build most recent. Format: usher_barcodes_yyyy-mm-dd.csv") #Optional
	parser.add_argument("-u", "--update", action="store_true",
		     help = "Will update barcodes the latest barcodes") #Optional
	parser.add_argument("-o", "--output", default= '.',
		     help="Default: Current working directory. Will make directory specified if not found") #Optional
	parser.add_argument("-n", "--parallel", default = 4,
		     help = "Default: 4 Number of samples run in parallel.") #Optional
	parser.add_argument("-f", "--file_check", action="store_true",
		     help = "Option generates a file checking file/figure to confirm if everything was created") #Optional
	parser.add_argument("-d", "--date",
		     help = "If you want to download a specific usher barcode date. Put date in yyyy-mm-dd") #Optional
	args = parser.parse_args()

	#### ARGUMENTS && MISC FILE CHECKING ####

	#Workspace (Path) [required]
	#Workspace is the directory where all files should be located (references, databases, scripts, and barcodes)
	#This is the default location of all the files unless stated
	workspace_path = args.workspace_path #commandline argument
	#Ensure "/" at the end of the workspace absolute path
	if ((workspace_path[len(workspace_path)-1]) != "/"):
		workspace_path = "{0}/".format(workspace_path)
	if not os.path.isdir(workspace_path):
		print("Workspace directory does not exist")
		raise SystemExit(1)
	print("Chosen Workspace is : '{0}'".format(workspace_path))
	
	#FASTQ file directory of samples [optional]: Default current directory
	#If sequence data is stored in a different location
	sample_path = args.sample_path 
	#Ensure "/" at the end of the workspace absolute path, and check that path specified exists
	if ((sample_path[len(sample_path)-1]) != "/"):
		sample_path = "{0}/".format(sample_path)
	if not os.path.isdir(sample_path):
		print("Sample directory does not exist")
		raise SystemExit(1)
	print("Sample Path is : '{0}'".format(sample_path))

	#Sample list file should be a simple text file with each sample on a new line and located in the workspace path
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
	#Will check if barcode is up to date and inform user wants to update the barcodes either with 1) freyja update or 2) no-update
	usherbarcodes = args.barcode
	today = date.today()
	down_date = today - timedelta(days = 1) #ushebarcode dates are always one day old
	down_date =  "usher_barcodes_{0}.csv".format(down_date)
	if not(args.barcode) or args.update or args.date:
		if args.update:
			print("Using usher barcode update pulling the most recent usher barcodes" + str(down_date))
			usherbarcodes = down_date
			os.system("python3 {0}/usher_update_V4.py {0}".format(workspace_path))
		elif args.date:
			down_date2 =  "usher_barcodes_{0}.csv".format(args.date)
			print("Downloading the specified barcodes" + str(down_date2) + "Reminder: Date should be formatted yyyy-mm-dd")
			usherbarcodes = down_date2
			os.system("python3 {0}/usher_update_V4.py {0} -d {1}".format(workspace_path, args.date))
		elif not(args.barcode):
			print("No barcode information provided Downloading most recent barcodes" + str(down_date))
			usherbarcodes = down_date
			os.system("python3 {0}/usher_update_V4.py {0}".format(workspace_path))
	if usherbarcodes != down_date:
			print("Warning: Expected to use usher barcode " + down_date + " To run using the most recent barcode run using --update")

	print("Barcodes used in analysis: {0}".format(usherbarcodes))
	
	#Path for directory output of pipeline
	output = args.output
	if not(args.output):
		output = os.getcwd()
	#Making sure to have "/" at the end of the workspace absolute path
	if ((output[len(output)-1]) != "/"):
		output = "{0}/".format(output)
	print("Output Path is : '{0}'".format(output))
	
	#Make directory output if does not exist
	if not os.path.exists(output + "results"):
		os.system("mkdir -p {0}results/".format(output))

	#Making sure that Human refrence genome BWA index is downloaded
	if not os.path.exists("{0}Homo_sapiens.GRCh38.fa.bwt".format(workspace_path)):
		os.system("wget --quiet -O {0}Homo_sapiens.GRCh38.dna.alt.fa.gz http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz".format(workspace_path))
		os.system("gzip -d {0}Homo_sapiens.GRCh38.dna.alt.fa.gz".format(workspace_path))
		os.system("mv {0}Homo_sapiens.GRCh38.dna.alt.fa {0}Homo_sapiens.GRCh38.fa".format(workspace_path))
		os.system("bwa index -p {0}Homo_sapiens.GRCh38.fa {0}Homo_sapiens.GRCh38.fa".format(workspace_path))

	#Making sure that Wuhan reference SARS genome BWA index is downloaded
	if not os.path.exists("{0}MN908947_3.fa.bwt".format(workspace_path)):
		os.system("wget --quiet -O {0}MN908947_3.fa https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true".format(workspace_path))
		os.system("bwa index -p {0}MN908947_3.fa {0}MN908947_3.fa".format(workspace_path))
	
	#Making that Artic V4.1 bed file is downloaded
	if not os.path.exists("{0}SARS-CoV-2.primer.bed".format(workspace_path)):
		os.system("wget -O {0}SARS-CoV-2.primer.bed https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V4.1/SARS-CoV-2.primer.bed".format(workspace_path))

	#Create file that records samples that have already been analyzed
	os.system("/bin/touch {0}Already_analyzed_samples.txt".format(workspace_path))
	
	#Number of samples being analysed in parallel (Default 4 samples at a time)
	nb_sim_process = int(args.parallel)
	print('Running ' + str(nb_sim_process) + ' of samples simultaneously')
	
    #Number of threads to use for multithread tasks (fastp and bwa): Default 4
	threads = int(args.threads)
	print('Using ' + str(threads) + 'threads per ' + str(nb_sim_process) + ' samples analyzed simultaneously')
	
	### STARTING ANALYSIS ###
	print("File checking complete, running analysis")
	#Analyze samples in parrallel (4 by default or whichever number selected by --parallel)
	liste_index_sample = range(0, num_samples, nb_sim_process)
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
		if (i + nb_sim_process) < num_samples:
			sub_list_index = range(i, (i+nb_sim_process), 1)
			for t in sub_list_index:
				if (not(lst_samples[t] in already_analyzed_samples_list)):
					print("Starting analysis on sample : " + lst_samples[t])
					globals()["p"+str(t)] = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[t], threads,sample_path,usherbarcodes, output,))
					processes.append(globals()["p"+str(t)])
					globals()["p"+str(t)].start()
				else:
					print(lst_samples[t] + "has already been analyzed, if re-analyzing, remove from Already_analyzed_samples.txt")
		else:
			sub_list_index = range(i,num_samples,1)
			for t in sub_list_index:
				if (not(lst_samples[t] in already_analyzed_samples_list)):
					print("Starting analysis on sample : " + lst_samples[t])
					globals()["p"+str(t)] = mp.Process(target=run_core_iPMVC, args=(workspace_path, lst_samples[t], threads,sample_path,usherbarcodes, output,))
					processes.append(globals()["p"+str(t)])
					globals()["p"+str(t)].start()
				else:
					print(lst_samples[t] + "has already been analyzed, if re-analyzing, remove from Already_analyzed_samples.txt")
		if (len(processes) != 0):
			for p in processes:
				p.join()
			for p in processes:
				while (p.is_alive()):
					time.sleep(5)

	#Run file checking
	if (args.file_check):
		os.system("python3 {0}file_checking_V2.py {0}{1} {2} {3}".format(workspace_path, samples_list_file, output, sample_path))
		os.system("Rscript {1}file_check_visualize_V1.R {0}File_Check_Output.txt {0}".format(output, workspace_path))