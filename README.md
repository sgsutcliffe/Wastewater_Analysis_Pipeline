# Running Wastewater Analysis Instructions

A pipeline for SARS-CoV-2 analysis on wastewater samples. Created specifically for my needs of analyzing pair-short-read sequences from tiled-amplicon PCR using Arctic V4.1 primer-set of SARS-CoV-2. The pipeline is under active development as of April 19, 2023, and I will be adding features and updates to make it more versatile. It is largely to python scripts. The first includes sample preprocessing, Freyja relative abundance estimates of lineages, and mutation calling. The second runs a sample quality control pipeline and outputs a MultiQC report and quantative report.

You can install all the software required (below) to run the pipeline, which is a lot, or build the singularity container to run everything.

Dependencies (versions listed have been tested):
- Python 3 (tested on 3.9.12)
   - pyvcf 0.6.8
   - numpy 1.22.3
   - pandas 1.4.3
- R version 4.2.2
- Samtools 1.16.1
- Picard 2.27.4
- iVar 1.3.1
- Freyja 1.3.12
- BWA 0.7.17
- Fastp 0.23.2
- FastQC 0.11.9
- Kraken2 2.1.2
- MultiQC 1.13
- snpEff 4.5 covid19
- bedtools 2.30.0

Optional:
   - Singularity > 3.0 (tested on 3.6.4)

Requirements:
- USHER barcode update requires at least 8GB of memory. (If you get error lineagePaths.txt not found likely a memory issue)

## Get reference files required files

In the workspace directory, or the directory where you store all the files Freyja pipeline needs to run:

Files required:
* sample_list #This you will need to make yourself depending on the samples [See more section below]

Files will be downloaded into the working directory if not present when running Freyja_parallelization.py:    
* Homo_sapiens.GRCh38.fa/.fa.amb/.fa.ann/.fa.bwt/.fa.pac/.fa.sa/ # BWA index for human genome
   - Note: This step takes a long time. I'd suggest downloading a pre-indexed human genome for bwa. For digital research alliance users you can get it at ``` /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/ ```
* MN908947_3.fasta/.amb/.ann/.bwt/.fai/.pac/.fasta.sa/.gb #BWA index of the reference Wuhan strain Accession: MN908947 
* SARS-CoV-2.primer.bed #Artic V4.1 primer BED file for iVar primer trimming

To build yourself:
* Kraken2 database built with contaminants the user wants to look for. Instructions below (needs a minimum SARS-CoV-2 and Human Genome)
* snpEff database [if not using singularity container]

## Building Kraken2 database

```shell
kraken2-build --standard --threads 24  --db /path/to/run/location/Kraken2wViruses
```
Note: This will build the entire Kraken2 database which is quite large and require a lot of memory to build (~100GB). For instructions on customizing the build and speeding up the process using multiple-threads see https://github.com/DerrickWood/kraken2/wiki/Manual I would recommend including all viruses and human-genome in the build at the very least.

I don't know if it is my build or server but I could not get a build to work. I instead used one from here:
https://benlangmead.github.io/aws-indexes/k2
For my analysis it was https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz 

### Recipe for making sample list

The paired reads that come from C3G have the extensions and right now the pipeline expects this:
R1 = <sample_name>_R1_001.fastq.gz
R2 = <sample_name>_R2_001.fastq.gz

Note: I know this is a drag for anyone using reads with a different naming scheme. I will fix this if indeed anyone else uses this pipeline! So ping me.

The 'sample_list' file is just a text file with all the <sample_name> involved in seperate line. I use a quick bash one liner to make the list

```shell
ls <path>/<to>/<reads>/*.fastq.gz | sed "s/<path>/<to>/<reads>\//" | sed "s/_R[0-9]_001.fastq.gz//" | sort | uniq > <path>/<to>/<workspace>/sample_list
```

## Running the pipeline

There are two scripts Freyja_parallelization.py and QC_parallelization.py

The first Freyja_parallelization.py includes the sample preprocessing, Freyja, and snpEff/iVar variant calling.

```shell
python3 Freyja_parallelization.py Workspace Sample-list
```
By default it will use the most up-to-date Usherbarcodes, and run 4 samples at a time with 4 threads per sample.
This assumes everything not in Workspace directory is in the current directory.

### positional arguments:
  workspace_path        The path to directory where: reference files, databases, sample list, usherbarcodes are.
  samples_list_file     Text file with samples to be run. Used to identify input and output files. The pattern is <sample_name>_R[1-2]_001.fastq.gz

### optional arguments:
  - -h, --help            show this help message and exit
  - -s SAMPLE_PATH, --sample_path SAMPLE_PATH
                        Default: Current working directory.
  - -t THREADS, --threads THREADS
                        Default: 4 Threads available should be > threads multiplied by samples run in parallel
  - -b BARCODE, --barcode BARCODE
                        Default: Will download/build most recent. Format: usher_barcodes_yyyy-mm-dd.csv
  - -u, --update          Will update barcodes the latest barcodes
  - -o OUTPUT, --output OUTPUT
                        Default: Current working directory. Will make directory specified if not found
  - -n PARALLEL, --parallel PARALLEL
                        Default: 4 Number of samples run in parallel.
  - -f, --file_check      Option generates a file checking file/figure to confirm if everything was created
  - -d DATE, --date DATE  If you want to download a specific usher barcode date. Put date in yyyy-mm-dd

After you've run this step you can run the QC_parallelization.py

```shell
python3 QC_parallelization.py Workspace Sample-list
```

### positional arguments:
  workspace_path        The path to directory where: reference files, databases, sample list, usherbarcodes are.
  samples_list_file     Text file with samples to be run. Used to identify input and output files. The pattern is <sample_name>_R[1-2]_001.fastq.gz

### optional arguments:
  - -h, --help            show this help message and exit
  - -s SAMPLE_PATH, --sample_path SAMPLE_PATH
                        Default: Current working directory. The location of raw FASTQ files
  - -i INPUT_PATH, --input_path INPUT_PATH
                        Default: Current working directory. The ouput directory location of previous Freyja pipeline
  - -t THREADS, --threads THREADS
                        Default: 4 Threads available should be > threads multiplied by samples run in parallel
  - -o OUTPUT, --output OUTPUT
                        Default: Current working directory. Will make directory specified if not found
  - -n PARALLEL, --parallel PARALLEL
                        Default: 4 Number of samples run in parallel.
  - -k KRAKEN, --kraken KRAKEN
                        Default: Kraken2wViruses, Directory of Kraken database. If running on Digital Alliance servers see option -mugqic
  - -x, --single          When you have too many files for one MultiQC report run each seperately


## Running Freyja using the singularity container

The first step is to build the container/

## Building singularity container

You'll need Singularity Version 3 or greater installed and root permissions

The list of dependencies will be included in Freyja

```shell

sudo singularity build Freyja.sif Singularity-Definition-Files/Freyja_Definition_File

```
Note: This creates a compressed read-only squashfs file system suitable for production. If you want to have a writable file-system (writable ext3 or (ch)root directory) try the options --writable or --sandbox. Respectively. You can also name the container whatever you like!

## Running singularity container

Since everything will be installed within the container you just need execute the container before running the script.

The simplest way is
```shell
singularity exec Freyja.sif python3 ...
```

Some suggestions for a more complex way to run it is
```shell
singularity exec -C -B ${PWD}:${HOME}
```
This keeps other install issues from conflicting. This binds your current working directory to the containers home directory.

One last suggestion is that some tools like MultiQC want to write to the temp directory but the space here is not what you have available on your computer but within the build.
My suggesttion is to make some temp directory on your file system for the container to use.

In this example I made directories /var/tmp and /tmp which are bound to the containers /var/tmp and /tmp respectively then run 

```shell
singularity exec -C -B ${PWD}:${HOME},/your/path/var/tmp:/var/tmp\,/your/path/tmp:/tmp python ...
```

## Overview

![](Misc/Analysis%20Workflow.png)

