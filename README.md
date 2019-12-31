# LSF ChIPseq pipeline

This pipeline is a wrapper around the Kundaje lab chipseq_pipeline2 for use with Pegasus at SCCC at the University of Miami.
https://github.com/ENCODE-DCC/chip-seq-pipeline2.  (Optimized for commit 83789e5. Can't guarantee updated versions.)

## Instructions for installation:

1. Install Conda:
	- https://conda.io/miniconda.html
2. Install the ENCODE3 chip-seq-pipeline2:
	- git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
	- cd chip-seq-pipeline2
	- Optional: git checkout 83789e5
	- bash conda/uninstall_dependencies.sh
	- bash conda/install_dependencies.sh
3. Build relevant genomic indicies:
	- bash conda/build_genome_data.sh [hg38/mm10/hg19] /path/to/built/genomes/
4. Install additional dependencies:
	- module rm python share-rpms65
	- conda env create -f chrome_chip_env.yml
	- (soon: conda install -n chrome_chip -c diderote chrome_chip)
5. OPTIONAL SETUP:
	- To allow for contamination screen (fastq screen):
		- Generate the relevant bowtie2 indices and add the path were indicated in options_files/fastq_screen.conf. (many premade by illumina at https://support.illumina.com/sequencing/sequencing_software/igenome.html)
		- In the fastq_screen.conf file included with chrome_chip, unhash relvant lines before each desired DATABASE
		- cp fastq_screen.conf `~/miniconda3/envs/chrome_chip/share/fastq-screen*/`

## Usage

1. Copy the Experimental_Samples.xlsx and chrome_chip_config.yml to a new folder and modify contents for the experiment. 
2. To run (with python3):
	> chrome_chip -f chrome_chip_config.yml -s -p [your_lsf_project]
	or 
	> python chrome_chip -f chrome_chip_config.yml -s -p [your_lsf_project]

	Options: 
	- Add '-t /path/to/ChIPseq.ipynb' if you running as a jupyter notebook
	- Add '-o /path/to/output_ChIPseq.ipynb' to name your output notebook something other than your experimental file name.
	- For all options: 'chrome_chip --help'

Config File Details:
* Restart: (yes/no) Yes restarts the analysis from the beginning.  No picks up from last completed step.
* LSF_project: Project name for LSF.  If none is needed, leave this blank.
* Scratch_folder: locaiton of temporary folder for analysis execution.  If black, will use tmp in current folder.
* Spike_index: Bowtie2 index for aligning spike-in chromatin.
* Genome_tsv: paths to ENCODE3 installed genome tsv files.  If mulitple species in one experiment, separate by commas.

This pipline handles processing and analyses for ChIPseq data on the University of Miami's Pegasus Computer Cluster using and LSF resource manager.  (Can be readily adpated to other resource managers: SGE, SLURM or local)  When fully implemented from start to finish it performs the following tasks:

1. Screening for contamination of other genomes
2. Fastq quality measurements
3. Adapter and quality fastq trimming
4. Spike-in alignment and analysis (if applicable)
5. ENCODE3 pipeline:
	- Genome alignment using BWA (newer versions use bowtie2)
	- Signal to Noise measurements
	- PCR Bottleneck measurements
	- Cross-correlation analysis
	- deduplication
	- peak calling using MACS2
	- signal file generation
	- IDR analysis
6. UMI decovolution and deduplication
7. QC: alignment metrics, gc content, library complexity with preseq, and other qc metrics
8. PCA analysis
9. Overlap analysis
9. Peak annotation
10. Enrichment analysis of annotated peaks and overlapped peaks (enrichr: KEGG, GO Biological Process, ChIP-X, ChEA, OMIM Disease)
Future tasks:
11. Differential binding analysis (using spike-in or not)
12. Lowess normlization of bigwig signal by spike-in

The pipeline handles multiple entry/exit points.  In case of error, the pipeline restarts from the last completed step. Progress is tracked in a .log file in the output directory or a jupyter notebook using papermill.

