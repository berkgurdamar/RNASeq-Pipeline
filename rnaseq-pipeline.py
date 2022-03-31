#!/usr/bin/env python
##################################################
## Project: RNASeq-Pipeline
## Purpose: Main function for pipeline run
## Date: March 2022
## Author: Berk Gürdamar
##################################################

import os
import subprocess
import argparse
import yaml

os.environ["COLUMNS"] = "1000"

parser = argparse.ArgumentParser(description = "RNASeq-Pipeline", formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=40))
parser.add_argument("-s", "--samples", type = str, metavar = '', required = True, help = "Sample names", nargs='+')
parser.add_argument("-c", "--controls", type = str, metavar = '', required = True, help = "Control names", nargs='+')
parser.add_argument("-r", "--replicate", type = str, metavar = '', required = False, default = "NO", help = "Technical replicate condition = YES or NO (default)")
parser.add_argument("-m", "--replicate_samples", type = str, metavar = '', required = False, help = "if --replicate == YES; Comma seperated matched sample names; sample1_sample2, sample3_sample4", nargs='+')
parser.add_argument("-t", "--type", type = str, metavar = '', required = False, nargs='?', default = "paired", help = "Sample type; paired (default), single")
parser.add_argument("-o", "--output", type = str, metavar = '', required = False, nargs='?', default = "pathway_enrichment", help = "Final output; pathway_enrichment (default), trimming, mapping, quantification, diff_exp")
parser.add_argument("-f", "--is_fastqc", type = str, metavar = '', required = False, nargs='?', default = "YES", help = "Run FastQC or not; YES (default), NO")
parser.add_argument("-d", "--data_dir", type = str, metavar = '', required = True, help = "Data folder direction")
parser.add_argument("-R", "--reference_genome", type = str, metavar = '', required = True, help = "Reference genome path")
parser.add_argument("-G", "--gtf_file", type = str, metavar = '', required = True, help = "GTF file path")
parser.add_argument("-x", "--star_index", type = str, metavar = '', required = True, help = "Folder path to create genome index (if already exist, automatically use the indexed genome in the path)")
parser.add_argument("-T", "--threads", type = int, metavar = '', required = False, nargs='?', default = 16, help = "Number of thread (default = 16)")
parser.add_argument("-q", "--quality", type = int, metavar = '', required = False, nargs='?', default = 20, help = "Trim galore quality option (default = 20)")
parser.add_argument("-l", "--length", type = int, metavar = '', required = False, nargs='?', default = 35, help = "Trim galore length option (default = 35)")
parser.add_argument("-i", "--genomeSAsparseD", type = int, metavar = '', required = False, nargs='?', default = 2, help = "genomeSAsparseD argument for STAR (default = 2)")
parser.add_argument("-g", "--runMode", type = str, metavar = '', required = False, nargs='?', default = "genomeGenerate", help = "runMode argument for STAR (default = genomeGenerate)")
parser.add_argument("-M", "--limitGenomeGenerateRAM", type = int, metavar = '', required = False, nargs='?', default = 164282414959, help = "Available maximum RAM size in bytes (default = 164282414959)")
parser.add_argument("-n", "--genomeSAindexNbases", type = int, metavar = '', required = False, nargs='?', default = 14, help = "genomeSAindexNbases argument for STAR (default = 14)")
parser.add_argument("-B", "--genomeChrBinNbits", type = int, metavar = '', required = False, nargs='?', default = 18, help = "genomeChrBinNbits argument for STAR (default = 18)")
parser.add_argument("-y", "--outSAMtype", type = str, metavar = '', required = False, nargs='?', default = "BAM Unsorted", help = "Output SAM file type (default = BAM Unsorted)")
parser.add_argument("-p", "--twopassMode", type = str, metavar = '', required = False, nargs='?', default = "Basic", help = "twopassMode argument for STAR (default = Basic)")
parser.add_argument("-N", "--twopass1readsN", type = int, metavar = '', required = False, nargs='?', default = -1, help = "twopass1readsN argument for STAR (default = -1)")
parser.add_argument("-F", "--readFilesCommand", type = str, metavar = '', required = False, nargs='?', default = "zcat", help = "readFilesCommand argument for STAR (default = zcat)")
parser.add_argument("-L", "--library_type", type = str, metavar = '', required = False, nargs='?', default = "A", help = "Library type for Salmon (default = A)")
parser.add_argument("-P", "--p_val_threshold", type = int, metavar = '', required = False, nargs='?', default = 0.05, help = "p-value threshold for pathfindR (default = 0.05)")
parser.add_argument("-I", "--iterations", type = int, metavar = '', required = False, nargs='?', default = 25, help = "iteration number for pathfindR (default = 25)")
parser.add_argument("-D", "--create_DAG", type = str, metavar = '', required = False, nargs='?', default = "NO", help = "Create Directed Acyclic Graph (DAG) of commands (default = NO)")
args = parser.parse_args()


def run_pipeline(samples, controls, replicate, replicate_samples, data_dir, reference_genome, gtf_file, star_index,
                 type, output, is_fastqc, threads, quality,
                 length, genomeSAsparseD, runMode, limitGenomeGenerateRAM,
                 genomeSAindexNbases, genomeChrBinNbits, outSAMtype, twopassMode,
                 twopass1readsN, readFilesCommand, library_type, p_val_threshold, iterations,
                 create_DAG):


    analysis = dict(
        analysis = dict(
            sample = samples, # ", ".join(samples)
            control = controls, # ", ".join(controls)
            replicate = replicate,
            replicate_samples = replicate_samples,
            data_dir = data_dir,
            type = type,
            output = output,
            is_fastqc = is_fastqc
        ),
        required = dict(
            threads = threads,
            trimgalore = dict(
                quality = quality,
                length = length
            ),
            star = dict(
                reference_genome = reference_genome,
                star_gtf = gtf_file,
                star_index = star_index,
                genomeSAsparseD = genomeSAsparseD,
                genomeGenerate = runMode,
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
                iterations = iterations
            )
        )
    )

    with open('config.yaml', 'w') as outfile:
        yaml.dump(analysis, outfile, default_flow_style=False)

    if create_DAG == "YES":
        subprocess.run(["""snakemake --configfile config.yaml --dag | dot -Tpdf > pipeline_dag.pdf"""], shell=True)
    else:
        subprocess.run(["""snakemake --configfile config.yaml --cores {threads}""".format(threads = threads)], shell=True)



if __name__ == '__main__':
    print(run_pipeline(args.samples, args.controls, args.replicate, args.replicate_samples, args.data_dir, args.reference_genome,
                       args.gtf_file, args.star_index,args.type, args.output, args.is_fastqc,
                       args.threads, args.quality,args.length, args.genomeSAsparseD, args.runMode, args.limitGenomeGenerateRAM,
                       args.genomeSAindexNbases, args.genomeChrBinNbits, args.outSAMtype, args.twopassMode, args.twopass1readsN,
                       args.readFilesCommand, args.library_type, args.p_val_threshold, args.iterations, args.create_DAG))
