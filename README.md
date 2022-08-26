# Running Freyja Instructions

Dependencies:
- Singularity > 3.0 (tested on 3.6.4)
- Python 3 (tested on 3.7.7)

## Required Files for Running

In the workspace directory, or the directort where you store all the files Freyja pipeline needs to run:

Files required:

* Freyja_parallelization_V3.1.py #Python script for parallelization execution of QC, Freyja update, and Freyja. Runs eight samples simultaneously
* sample_list #This you will need to make yourself depending on the samples [See more information below]
* usher_update_V2.py #This is a modified version of Freyja's update function to fix some bugs that came up running it on some servers

Files if not downloaded will be downloaded:    
   Note: The downloading of human genome & building bwa index files take a while. I'd suggest running this first without any samples before analyzing samples.
* Homo_sapiens.GRCh38.fa/.fa.amb/.fa.ann/.fa.bwt/.fa.pac/.fa.sa/ # BWA index for human genome
* MN908947_3.fasta/.amb/.ann/.bwt/.fai/.pac/.fasta.sa/.gb #BWA index of the reference Wuhan strain Accession: MN908947 
* SARS-CoV-2.primer.bed #Artic V4.1 primer BED file for iVar primer trimming

You'll also need the container of singularity for Freyja

* Freyja_V5-1.sif #You'll need to build with sudo permissions
* Freyja_Definition_File_V5.1 #This is the definition file for building the container

### Building singularity container

You'll need Singularity Version 3 or greater installed and root permissions

```shell

sudo singularity build Freyja_V5-1.sif Freyja_Definition_File_V5.1

```

### Recipe for making sample list

The paired reads that come from C3G have the extensions and right now the pipeline expects this:
R1 = <sample_name>_R1_001.fastq.gz
R2 = <sample_name>_R2_001.fastq.gz

Note: Sometimes samples are run in two different lanes and need to be merged. 

The 'sample_list' file is just a text file with all the <sample_name> involved in seperate line. I use a quick bash one liner to make the list

```shell
ls <path>/<to>/<reads>/*.fastq.gz | sed "s/<path>/<to>/<reads>\///" | sed "s/_R[0-9]_001.fastq.gz//" | sort | uniq > <path>/<to>/<workspace>/sample_list
```

## Running Freyja using the singularity container
Note: The pipeline is parallel in the sense it runs 8 samples at a time, regardless of how many blocks-of-memory/threads are available. This will be updated in the snakemake workflow. But keep this in mind when selecting threads per sample and how much memory you allocate per sample.

Right now the order of variables and location of where it is run is very important.

Where you execute the command is also important. You'll need to execute the singularity command in the parent directory of all the folders because you will need to bind this folder and all its sub-folders into the containers /mnt directory.

![](directory_example.jpeg)

directory_from_which_you_excute_the_singularity
    |
    | - mounted_directory (The parent directory to sub-directories involved in analysis. The /mnt directory in the container)
        |
        | - /path/to/container
        | - /path/to/working/directory
        | - /path/to/sequences/
        | - /etc/

For more information on mounting directories see: https://docs.sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html
But when creating variables to run the pipeline do not change the '/mnt' in the path

If you don't have any usherbarcodes already downloaded you can just put:
'usherbarcode.csv' and it will download the most recent version. Right now the pipeline is built to check to make sure the most recent barcode is used. Currently my pipeline does not have a way to bypass this. It only downloads and sets up if up-to-date barcodes are not available in expected naming scheme: 


```shell
mounted_directory=<mounted_directory_name>
container=/path/to/singularity_container/Freyja_V5-1.sif
run_location=/mnt/path/to/workspace/directory/
sample_location=/mnt/path/to/reads/
sample_list=/path/to/sample/list/<sample_list_file>
threads_per_sample=<interger>
barcode=<usher_barcode_file>
output_directory=/mnt/<output_folder_name>

singularity exec --no-home -B ${mounted_directory}:/mnt \
${container} python3 \
${run_location}Freyja_parallelization_V3.1.py ${run_location} ${sample_location} ${sample_list} ${threads_per_sample} ${barcode} ${output_directory}

```
