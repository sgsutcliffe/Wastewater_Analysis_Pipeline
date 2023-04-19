#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:05:01 2023

@author: Sutcliffe
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 13:11:11 2022
#Version 2: Integrates into Freyja parallelization pipeline
@author: Sutcliffe
"""

import os
import pandas as pd
import glob
import subprocess
import sys


#Using the the sample list file 
path_to_file = sys.argv[1]
path_to_output = sys.argv[2]
path_to_samples = sys.argv[3]
filename = path_to_output + 'File_Check_Output.txt'
data = pd.read_csv(path_to_file, header = None, sep='\t')

fastq_output = []
fastp_output = []
bwa_decon = []
first_align = []
second_align = []
freyja_variants = []
freyja_demix = []
freya_criteria_fail = []

# with open(path_to_file) as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter='\t')
#     csv_reader = csv_reader.loc[csv_reader[1] == run_number]
    
for index, row in data.iterrows():
    sample = row[0]
    path_to_results = path_to_output
    #1 Check the FASTQ files (Forward & Reverse)
    if glob.glob(path_to_samples + sample + "_R2_001.fastq.gz"):
        if glob.glob(path_to_samples + sample + "_R1_001.fastq.gz"):
            fastq_output.append("True")
        else:
            fastq_output.append('False')
    elif glob.glob(path_to_samples + sample + "_R1_001.fastq.gz"):
        if glob.glob(path_to_samples + sample + "_R2_001.fastq.gz"):
            fastq_output.append("True")
        else:
            fastq_output.append('False')
    else:
        fastq_output.append('False')
    #2 Check for fastp output trimmed files (forward & Reverse)
    if glob.glob(path_to_results + sample + "_R1_001_trimmed_1.fastq") or glob.glob(path_to_results + sample + "_R1_001_trimmed_1.fastq.gz"):
        if glob.glob(path_to_results + sample + "_R2_001_trimmed_2.fastq*") or glob.glob(path_to_results + sample + "_R1_001_trimmed_1.fastq.gz"):
            fastp_output.append("True")
        else:
            fastp_output.append("False")
    else:
        fastp_output.append('False')
    #3 Check for decontaminated human reads
    if glob.glob(path_to_results + sample + "*_R1_001_decon_1.fastq") or glob.glob(path_to_results + sample + "*R1_001_decon_1.fastq.gz"):
        if glob.glob(path_to_results + sample + "*_R2_001_decon_2.fastq*") or glob.glob(path_to_results + sample + "*_R2_001_decon_2.fastq.gz"):
            bwa_decon.append("True")
        else:
            bwa_decon.append("False")
    else:
        bwa_decon.append('False')
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
            first_align.append("True")
        else:
            first_align.append('False')
    else:
        first_align.append('False')
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
            second_align.append("True")

        else:
            second_align.append('False')

    else:
        second_align.append('False')
    #6 Check Freyja variant output
    depthout = "ls " + path_to_results + sample + "*_depthout"
    depthout = subprocess.getoutput(depthout)
    variantout = "ls " + path_to_results + sample + "*_variantout.tsv"
    variantout = subprocess.getoutput(variantout)
    if glob.glob(depthout) and glob.glob(variantout):
        if (os.path.getsize(depthout) == 0) or (os.path.getsize(variantout) == 0):
            freyja_variants.append("False")
        else:
            freyja_variants.append('True')
    else:
        freyja_variants.append('False')

    #7 Check Freyja demix output exists
    if glob.glob(path_to_results + "results/" + sample + "*_output"):
        outputfile = "ls " + path_to_results + "results/" + sample + "*_output"
        outputfile = subprocess.getoutput(outputfile)
        if (os.path.getsize(outputfile) == 0):
            freyja_demix.append("False")
            freya_criteria_fail.append("False")
        else: 
            freyja_demix.append("True")
            freyja_out = pd.read_csv(outputfile, header = None, sep='\t', skiprows=1)
            print(freyja_out)
            cov = float(freyja_out.iloc[4][1])
            if cov <= 0:
                freya_criteria_fail.append("False")
            else:
                freya_criteria_fail.append("True")
                
            
        
    else:
        freyja_demix.append("False")
        freya_criteria_fail.append("False")

        
 
data['FASTQ'] = fastq_output
data['fastp'] = fastp_output
data['bwa_decontaminate'] = bwa_decon
data['preiVar_sars_alignmt'] = first_align
data['postiVar_sars_alignmt'] = second_align
data['freyja_variant_outputs'] = freyja_variants
data['freya_final_output_exists'] = freyja_demix
data['freya_output_meets_mincriteria'] = freya_criteria_fail 
data.to_csv(filename, sep='\t')