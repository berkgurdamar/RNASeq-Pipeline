##################################################
## Project: RNASeq-Pipeline
## Purpose: Main function
## Date: September 2021
## Author: Berk Gürdamar
##################################################
import subprocess
import argparse
import yaml
import os

subprocess.call("run_env.sh", shell=True)

os.environ["COLUMNS"] = "1000"

parser = argparse.ArgumentParser(description = "Run RNASeq Pipeline", formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=52))
parser.add_argument("-s", "--samples", type = list, metavar = '', required = True, help = "Sample names")
parser.add_argument("-c", "--controls", type = list, metavar = '', required = True, help = "Control names")
parser.add_argument("-r", "--replicate", type = str, metavar = '', required = True, help = "Replicate condition = YES or NO")
parser.add_argument("-m", "--matched_samples", type = list, metavar = '', required = False, help = "Comma seperated combined matched sample names; sample1_sample2, sample3_sample4")
parser.add_argument("-t", "--type", type = str, metavar = '', required = False, help = "Sample type; paired (default), single")
parser.add_argument("-o", "--output", type = str, metavar = '', required = False, help = "Final output; pathway_enrichment (default), trimming, mapping, quantification, diff_exp")
parser.add_argument("-f", "--is_fastqc", type = str, metavar = '', required = False, help = "Run FastQC or not; YES (default), NO")
parser.add_argument("-d", "--data_dir", type = str, metavar = '', required = True, help = "Data folder direction")
parser.add_argument("-w", "--index_whole_genome", type = str, metavar = '', required = True, help = "Index genome path")
parser.add_argument("-G", "--star_gtf", type = str, metavar = '', required = True, help = "GTF file path")
parser.add_argument("-x", "--star_index", type = str, metavar = '', required = True, help = "Folder path to create genome index")
parser.add_argument("-T", "--threads", type = int, metavar = '', required = False, help = "Number of thread (default = 16)")
parser.add_argument("-q", "--quality", type = int, metavar = '', required = False, help = "Trim galore quality option (default = 20)")
parser.add_argument("-l", "--length", type = int, metavar = '', required = False, help = "Trim galore length option (default = 35)")
parser.add_argument("-i", "--genomeSAsparseD", type = int, metavar = '', required = False, help = "genomeSAsparseD argument for STAR (default = 2)")
parser.add_argument("-g", "--genomeGenerate", type = str, metavar = '', required = False, help = "genomeGenerate argument for STAR")
parser.add_argument("-M", "--limitGenomeGenerateRAM", type = int, metavar = '', required = False, help = "Available maximum RAM size in bytes (default = 164282414959)")
parser.add_argument("-n", "--genomeSAindexNbases", type = int, metavar = '', required = False, help = "genomeSAindexNbases argument for STAR (default = 14)")
parser.add_argument("-B", "--genomeChrBinNbits", type = int, metavar = '', required = False, help = "genomeChrBinNbits argument for STAR (default = 18)")
parser.add_argument("-y", "--outSAMtype", type = str, metavar = '', required = False, help = "Output SAM file type (default = BAM Unsorted)")
parser.add_argument("-p", "--twopassMode", type = str, metavar = '', required = False, help = "twopassMode argument for STAR (default = Basic)")
parser.add_argument("-N", "--twopass1readsN", type = int, metavar = '', required = False, help = "twopass1readsN argument for STAR (default = -1)")
parser.add_argument("-F", "--readFilesCommand", type = str, metavar = '', required = False, help = "readFilesCommand argument for STAR (default = zcat)")
parser.add_argument("-L", "--library_type", type = str, metavar = '', required = False, help = "Library type for Salmon (default = A)")
parser.add_argument("-P", "--p_val_threshold", type = int, metavar = '', required = False, help = "p-value threshold for pathfindR (default = 0.05)")
parser.add_argument("-I", "--iterations", type = int, metavar = '', required = False, help = "iterations number for pathfindR (default = 10)")
args = parser.parse_args()

def run_pipeline(samples, controls, replicate, matched_samples, data_dir, index_whole_genome, star_gtf, star_index,
                 type = "paired", output = "pathway_enrichment", is_fastqc = 'YES', threads = 16, quality = 20,
                 length = 35, genomeSAsparseD = 2, genomeGenerate = "genomeGenerate", limitGenomeGenerateRAM = 164282414959,
                 genomeSAindexNbases = 14, genomeChrBinNbits = 18, outSAMtype = "BAM Unsorted", twopassMode = "Basic",
                 twopass1readsN = -1, readFilesCommand = "zcat", library_type = "A", p_val_threshold = 0.05, iterations = 10):

    analysis = dict(
        analysis = dict(
            sample = ", ".join(samples),
            control = ", ".join(controls),
            replicate = replicate,
            matched_samples = ", ".join(matched_samples),
            type = type,
            output = output,
            is_fastqc = is_fastqc,
            data_dir = data_dir
        ),
        required = dict(
            threads = threads,
            trimgalore = dict(
                quality = quality,
                length = length
            ),
            star = dict(
                index_whole_genome = index_whole_genome,
                star_gtf = star_gtf,
                star_index = star_index,
                genomeSAsparseD = genomeSAsparseD,
                genomeGenerate = genomeGenerate,
                limitGenomeGenerateRAM = limitGenomeGenerateRAM,
                genomeSAindexNbases = genomeSAindexNbases,
                genomeChrBinNbits = genomeChrBinNbits,
                outSAMtype = outSAMtype,
                twopassMode = twopassMode,
                twopass1readsN = twopass1readsN,
                readFilesCommand = readFilesCommand
            ),
            salmon = dict(
                library_type = library_type
            ),
            pathfindr = dict(
                p_val_threshold = p_val_threshold,
                iterations = iterations,
                gene_name_converter = gene_name_converter
            )
        )
    )

    with open('config.yaml', 'w') as outfile:
        yaml.dump(analysis, outfile, default_flow_style=False)

    subprocess.call("run_snakemake.sh", shell=True)



if __name__ == '__main__':
    print(run_pipeline(args.samples, args.controls, args.replicate, args.matched_samples, args.type, args.output,
                       args.is_fastqc, args.data_dir, args.threads, args.quality, args.length, args.index_whole_genome,
                       args.star_gtf, args.star_index, args.genomeSAsparseD, args.genomeGenerate, args.limitGenomeGenerateRAM,
                       args.genomeSAindexNbases, args.genomeChrBinNbits, args.outSAMtype, args.twopassMode, args.twopass1readsN,
                       args.readFilesCommand, args.library_type, args.p_val_threshold, args.iterations))
