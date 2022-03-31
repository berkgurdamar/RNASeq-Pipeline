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

 All required packages can be installed using `setup.sh` script. setup.sh script creates a conda environment called `rnaseq-pipeline` and downloads all required packages in it. 

```
bash setup.sh
```
 
 ## Usage
 
 For using the pipeline, conda environment should be activated.
 
 ```
conda activate rnaseq-pipeline
```


All arguments can be seen by using `--help` command.

```
python rnaseq-pipeline.py --help
```


```
usage: rnaseq-pipeline.py [-h] -s  [...] -c  [...] [-r] [-m  [...]] [-t ] [-o ] [-f ] -d  -R  -G  -x  [-T ] [-q ] [-l ] [-i ] [-g ] [-M ] [-n ] [-B ] [-y ] [-p ] [-N ] [-F ] [-L ] [-P ] [-I ] [-D ]

RNASeq-Pipeline

optional arguments:
  -h, --help                            show this help message and exit
  -s  [ ...], --samples  [ ...]         Sample names
  -c  [ ...], --controls  [ ...]        Control names
  -r , --replicate                      Technical replicate condition = YES or NO (default)
  -m  [ ...], --replicate_samples  [ ...]
                                        if --replicate == YES; Comma seperated matched sample names; sample1_sample2, sample3_sample4
  -t [], --type []                      Sample type; paired (default), single
  -o [], --output []                    Final output; pathway_enrichment (default), trimming, mapping, quantification, diff_exp
  -f [], --is_fastqc []                 Run FastQC or not; YES (default), NO
  -d , --data_dir                       Data folder direction
  -R , --reference_genome               Reference genome path
  -G , --gtf_file                       GTF file path
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
  -D [], --create_DAG []                Create Directed Acyclic Graph (DAG) of commands (default = NO)


```

RNASeq-Pipeline works with gzipped .fastq files. Input names should be;

- For pair-end samples: `{sample_name}_{1 or 2}.fastq.gz`
- For single-end samples: `{sample_name}.fastq.gz`

## Example Command

If reference genome indexed before, RNASeq-Pipeline automatically uses the path in the `--star_index` argument.If not, index genome will be created into the `--star_index` path.

```
python rnaseq-pipeline.py --samples sample1 sample2 sample3 --controls control1 control2 control3 \
--data_dir /path/to/data/ --reference_genome /path/to/reference_genome --gtf_file /path/to/gtf_file \
--star_index /path/to/target/folder
```

Directed Acyclic Graph (DAG) of the analysis can be created by using `--create_DAG` argument.

```
python rnaseq-pipeline.py --samples sample1 sample2 sample3 --controls control1 control2 control3 \
--data_dir /path/to/data/ --reference_genome /path/to/reference_genome --gtf_file /path/to/gtf_file \
--star_index /path/to/target/folder --create_DAG YES
```

<img src="https://github.com/berkgurdamar/RNASeq-Pipeline/blob/main/workflow/pipeline_dag.pdf?raw=true" style="max-width:100%;" />

## Technical Replicates

If there are technical replicates in the samples, `--replicate` argument should be set as `YES` and technical replicates should be written with using `--replicate_samples` argument. RNASeq-Pipeline automatically combines technical replicates when running DESeq2.

### Example Command
```
python rnaseq-pipeline.py --samples sample1 sample2 sample3 sample4 --controls control1 control2 control3 control4 \
--replicate YES --replicate_samples sample1_sample2 sample3_sample4 control1_control2 control3_control4 \
--data_dir /path/to/data/ --reference_genome /path/to/reference_genome --gtf_file /path/to/gtf_file \
--star_index /path/to/target/folder
```


## Output

All the outputs will be written in seperate folders in the data folder for each step. 



