# RNASeq-Pipeline
 
 
 RNASeq-Pipeline is a bioinformatics pipeline that can be used for analyzing RNA Sequencing data. 
 RNASeq-Pipeline steps;
 
 - Quality Control with FastQC
 - Trimming with Trim-Galore
 - Mapping with STAR
 - Quantification with Salmon
 - Differential Gene Expression Analysis with DESeq2
 - Enrichment Analysis with pathfindR
 
## Installation

 All required packages can be installed using setup.sh script. setup.sh script creates a conda environment called rnaseq-pipeline and downloads all packages in it. 

'''
bash setup.sh
'''
 
 ## Usage

All arguments can be seen by using help (-h) command.

'''
conda activate rnaseq-pipeline

python rnaseq-pipeline.py -h

usage: rnaseq-pipeline.py [-h] -s  [...] -c  [...] [-r] [-m  [...]] [-t ] [-o ] [-f ] -d  -w  -G  -x  [-T ] [-q ] [-l ] [-i ] [-g ] [-M ] [-n ] [-B ] [-y ] [-p ] [-N ] [-F ] [-L ] [-P ] [-I ] [-C] [-D ]

Run RNASeq Pipeline

optional arguments:
  -h, --help                            show this help message and exit
  -s  [ ...], --samples  [ ...]         Sample names
  -c  [ ...], --controls  [ ...]        Control names
  -r , --replicate                      Technical replicate condition = YES or NO
  -m  [ ...], --matched_samples  [ ...]
                                        if --replicate == YES; Comma seperated matched sample names; sample1_sample2, sample3_sample4
  -t [], --type []                      Sample type; paired (default), single
  -o [], --output []                    Final output; pathway_enrichment (default), trimming, mapping, quantification, diff_exp
  -f [], --is_fastqc []                 Run FastQC or not; YES (default), NO
  -d , --data_dir                       Data folder direction
  -w , --index_whole_genome             Index genome path
  -G , --star_gtf                       GTF file path
  -x , --star_index                     Folder path to create genome index (if already exist, automatically use the indexed genome in the path)
  -T [], --threads []                   Number of thread (default = 16)
  -q [], --quality []                   Trim galore quality option (default = 20)
  -l [], --length []                    Trim galore length option (default = 35)
  -i [], --genomeSAsparseD []           genomeSAsparseD argument for STAR (default = 2)
  -g [], --runMode []                   runMode argument for STAR (default = genomeGenerate)
  -M [], --limitGenomeGenerateRAM []    Available maximum RAM size in bytes (default = 164282414959)
  -n [], --genomeSAindexNbases []       genomeSAindexNbases argument for STAR (default = 14)
  -B [], --genomeChrBinNbits []         genomeChrBinNbits argument for STAR (default = 18)
  -y [], --outSAMtype []                Output SAM file type (default = BAM Unsorted)
  -p [], --twopassMode []               twopassMode argument for STAR (default = Basic)
  -N [], --twopass1readsN []            twopass1readsN argument for STAR (default = -1)
  -F [], --readFilesCommand []          readFilesCommand argument for STAR (default = zcat)
  -L [], --library_type []              Library type for Salmon (default = A)
  -P [], --p_val_threshold []           p-value threshold for pathfindR (default = 0.05)
  -I [], --iterations []                iteration number for pathfindR (default = 25)
  -C , --gene_name_converter            'mart_export.txt' file path
  -D [], --create_DAG []                Create Directed Acyclic Graph (DAG) of commands (default = NO)

'''
