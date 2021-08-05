##################################################
## Project: RNASeq
## Purpose: Main Snakefile
## Date: August 2021
## Author: Berk Gürdamar
##################################################



SAMPLES = config["analysis"]["sample"] + config["analysis"]["control"]
PATH_SAMPLE = config["analysis"]["data_dir"]
PATH_QC = os.path.join(PATH_SAMPLE, "qc")
PATH_TRIM = os.path.join(PATH_SAMPLE, "trimming")
PATH_TRANSCRIPT_FA = os.path.join(PATH_SAMPLE, "gffread_file_path")
PATH_STAR_IND = config["required"]["star"]["star_index"]
PATH_STAR_ALIGN = os.path.join(PATH_SAMPLE, "star_alignment")
PATH_SALMON = os.path.join(PATH_SAMPLE, "salmon")
PATH_DESEQ = os.path.join(PATH_SAMPLE, "deseq")
PATH_PATHFINDR = os.path.join(PATH_SAMPLE, "pathfindr")
MATCHED_SAMPLES = config["analysis"]["matched_samples"]



if config["analysis"]["output"] == "trimming":
    if config["analysis"]["type"] == "paired":
        RULE_ALL = expand(os.path.join(PATH_TRIM, "{sample}_1_val_1.fq.gz"), sample = SAMPLES)
        RULE_ALL += expand(os.path.join(PATH_TRIM, "{sample}_2_val_2.fq.gz"), sample = SAMPLES)
    else:
        RULE_ALL = expand(os.path.join(PATH_TRIM, "{sample}_trimmed.fq.gz"), sample = SAMPLES)
elif config["analysis"]["output"] == "mapping":
    RULE_ALL = expand(os.path.join(PATH_STAR_ALIGN, "{sample}", "Aligned.toTranscriptome.out.bam"), sample = SAMPLES)
elif config["analysis"]["output"] == "quantification":
    RULE_ALL = expand(os.path.join(PATH_SALMON, "quant_from_bam", "{sample}", "quant.sf"), sample = SAMPLES)
elif config["analysis"]["output"] == "diff_exp":
    RULE_ALL = os.path.join(PATH_DESEQ, config["analysis"]["sample"][0] + "_" + config["analysis"]["control"][0], "deseq_output.csv")
else:
    RULE_ALL = [PATH_PATHFINDR]

if config["analysis"]["is_fastqc"] == "YES":
    if config["analysis"]["type"] == "paired":
       RULE_ALL += expand(os.path.join(PATH_QC, "{sample}_{FR}_fastqc.html"), sample = SAMPLES, FR = [1, 2])
    else:
       RULE_ALL += expand(os.path.join(PATH_QC, "{sample}_fastqc.html"), sample = SAMPLES)
# RULE_ALL += os.path.join(PATH_MULTIQC, "multiqc_report.html")

rule all:
    input : RULE_ALL





rule run_fastqc:
    input:
        input_fq = os.path.join(PATH_SAMPLE, "{fq_name}.fastq.gz")
    output:
        out_html = os.path.join(PATH_QC, "{fq_name}_fastqc.html")
    params:
        threads = round(config["required"]["threads"]/(len(SAMPLES)*2)),
        outdir = PATH_QC
    shell:
        "fastqc --outdir {params.outdir} --threads {params.threads} --quiet --noextract {input.input_fq}"






# rule run_multiqc:
#     input:
# 		expand(os.path.join(PATH_QC, "{fq_name}_{R12}{paired}_fastqc.html"), fq_name = SAMPLES, R12 = ["R1", "R2"], paired = ["", "_paired"])
# 	output:
# 		multiqc_report = os.path.join(PATH_MULTIQC, "multiqc_report.html")
# 	params:
# 		qc_path = PATH_QC,
#         multiqc_path = PATH_MULTIQC
# 	shell:
#       """
# 		export LC_ALL=en_US.utf-8
# 		export LANG=en_US.utf-8
# 		multiqc {params.qc_path} --outdir {params.multiqc_path} --no-data-dir
# 	    """




if config["analysis"]["type"] == "paired":
    rule trimgalore:
        input:
            sample1 = os.path.join(PATH_SAMPLE, "{fq_name}_1.fastq.gz"),
            sample2 = os.path.join(PATH_SAMPLE, "{fq_name}_2.fastq.gz")
        output:
            out_fq1 = os.path.join(PATH_TRIM, "{fq_name}_1_val_1.fq.gz"),
            out_fq2 = os.path.join(PATH_TRIM, "{fq_name}_2_val_2.fq.gz")
        params:
            cores = round(config["required"]["threads"]/(len(SAMPLES)*2)),
            outdir = PATH_TRIM,
            quality = config["required"]["trimgalore"]["quality"],
            length = config["required"]["trimgalore"]["length"]
        shell:
            """
            trim_galore --phred33 --quality {params.quality} --gzip --length {params.length} \
            --trim-n --output_dir {params.outdir} --retain_unpaired --cores {params.cores} --paired {input.sample1} {input.sample2}
            """
else:
    rule trimgalore:
        input:
            sample1 = os.path.join(PATH_SAMPLE, "{fq_name}.fastq.gz")
        output:
            out_fq1 = os.path.join(PATH_TRIM, "{fq_name}_trimmed.fq.gz")
        params:
            cores = round(config["required"]["threads"]/(len(SAMPLES)*2)),
            outdir = PATH_TRIM,
            quality = config["required"]["trimgalore"]["quality"],
            length = config["required"]["trimgalore"]["length"]
        shell:
            """
            trim_galore --phred33 --quality {params.quality} --gzip --length {params.length} \
            --trim-n --output_dir {params.outdir} --cores {params.cores} {input.sample1}
            """



rule run_star_index:
    input:
        whole_genome = config["required"]["star"]["index_whole_genome"]
    output:
        index_path = directory(PATH_STAR_IND)
    params:
        threads = config["required"]["threads"],
        gtf_file = config["required"]["star"]["star_gtf"],
        genomeSAsparseD = config["required"]["star"]["genomeSAsparseD"],
        genomeSAindexNbases = config["required"]["star"]["genomeSAindexNbases"],
        genomeChrBinNbits = config["required"]["star"]["genomeChrBinNbits"]
    shell:
        """
        STAR --runThreadN {params.threads} --runMode genomeGenerate \
        --genomeDir {output.index_path} --genomeFastaFiles {input.whole_genome} \
        --sjdbGTFfile {params.gtf_file} --sjdbOverhang 100 --limitGenomeGenerateRAM 164282414959
        """




if config["analysis"]["type"] == "paired":
    rule run_star_align_paired:
        input:
            input_fasta1 = os.path.join(PATH_TRIM, "{fq_name}_1_val_1.fq.gz"),
            input_fasta2 = os.path.join(PATH_TRIM, "{fq_name}_2_val_2.fq.gz"),
            genome_idx = PATH_STAR_IND
        output:
            out_name = os.path.join(PATH_STAR_ALIGN, "{fq_name}", "Aligned.toTranscriptome.out.bam")
        params:
            threads = round(config["required"]["threads"]/len(SAMPLES)),
            readFilesCommand = config["required"]["star"]["readFilesCommand"],
            outSAMtype = config["required"]["star"]["outSAMtype"],
            twopassMode = config["required"]["star"]["twopassMode"],
            twopass1readsN = config["required"]["star"]["twopass1readsN"],
            out_file = directory(os.path.join(PATH_STAR_ALIGN, "{fq_name}"))

        shell:
            """
            STAR --runMode alignReads --genomeDir {input.genome_idx} --quantMode TranscriptomeSAM --runThreadN {params.threads} \
            --readFilesIn  {input.input_fasta1}  {input.input_fasta2}  \
            --readFilesCommand {params.readFilesCommand}  \
            --outSAMtype {params.outSAMtype} \
            --twopassMode {params.twopassMode} \
            --outFileNamePrefix {params.out_file}/
            """
else:
    rule run_star_align_single:
        input:
            input_fasta = os.path.join(PATH_TRIM, "{fq_name}_trimmed.fq.gz"),
            genome_idx = PATH_STAR_IND
        output:
            out_name = os.path.join(PATH_STAR_ALIGN, "{fq_name}", "Aligned.toTranscriptome.out.bam")
        params:
            threads = round(config["required"]["threads"]/len(SAMPLES)),
            readFilesCommand = config["required"]["star"]["readFilesCommand"],
            outSAMtype = config["required"]["star"]["outSAMtype"],
            twopassMode = config["required"]["star"]["twopassMode"],
            twopass1readsN = config["required"]["star"]["twopass1readsN"],
            out_file = directory(os.path.join(PATH_STAR_ALIGN, "{fq_name}"))

        shell:
            """
            STAR --runMode alignReads --genomeDir {input.genome_idx} --quantMode TranscriptomeSAM --runThreadN {params.threads} \
            --readFilesIn  {input.input_fasta} \
            --readFilesCommand {params.readFilesCommand}  \
            --outSAMtype {params.outSAMtype} \
            --twopassMode {params.twopassMode} \
            --outFileNamePrefix {params.out_file}/
            """


# rule run_salmon_index:
#         input:
#             cdna_fasta = config["required"]["salmon"]["salmon_cdna"]
#         output:
#             salmon_index_file = os.path.join(PATH_SALMON, 'salmon_index', "sa.bin")
#         params:
#             threads = config["required"]["threads"],
#             salmon_index_dir = os.path.join(PATH_SALMON, 'salmon_index')
#         shell:
#             """
#             salmon index -t {input.cdna_fasta} -i {params.salmon_index_dir} -p {params.threads}
#             """
#
#
#
#
# rule run_salmon_quant:
#         input:
#             index_file = os.path.join(PATH_SALMON, 'salmon_index', "sa.bin"),
#             R1_paired = os.path.join(PATH_TRIM, "{fq_name}_R1_paired.fastq.gz"),
#             R2_paired = os.path.join(PATH_TRIM, "{fq_name}_R2_paired.fastq.gz"),
#         output:
#             out_name = os.path.join(PATH_SALMON, "{fq_name}", "quant.sf"),
#             out_name_gene = os.path.join(PATH_SALMON, "{fq_name}", "quant.genes.sf")
#         params:
#             threads = config["required"]["threads"],
#             index_dir = os.path.join(PATH_SALMON, 'salmon_index'),
#             gtf_file = config["required"]["star"]["star_gtf"],
#             outfolder = directory(os.path.join(PATH_SALMON, "{fq_name}"))
#         shell:
#             """
#             salmon quant -i {params.index_dir} -l A -p {params.threads} -1 {input.R1_paired} -2 {input.R2_paired} --validateMappings -o {params.outfolder} --seqBias --gcBias -g {params.gtf_file}
#
#             """



rule generate_transcript_file_gffread:
    input:
        whole_genome = config["required"]["star"]["index_whole_genome"]
    output:
        transcript_file = os.path.join(PATH_TRANSCRIPT_FA, "all_transcripts.fa")
    params:
        gtf_file = config["required"]["star"]["star_gtf"]
    shell:
        """
        gffread -w {output.transcript_file} -g {input.whole_genome}  {params.gtf_file}
        """


rule run_salmon:
    input:
        transcript_file = os.path.join(PATH_TRANSCRIPT_FA, "all_transcripts.fa"),
        salmon_bam = os.path.join(PATH_STAR_ALIGN, "{fq_name}", "Aligned.toTranscriptome.out.bam")
    output:
        out_name = os.path.join(PATH_SALMON, "quant_from_bam", "{fq_name}", "quant.sf")
    params:
        threads = round(config["required"]["threads"]/len(SAMPLES)),
        library_type = config["required"]["salmon"]["library_type"],
        out_dir = directory(os.path.join(PATH_SALMON, "quant_from_bam", "{fq_name}"))
    shell:
        """
        salmon quant -p {params.threads} -t {input.transcript_file} -l {params.library_type} -a {input.salmon_bam} -o {params.out_dir}
        """



rule run_deseq:
    input:
        deseq_input = expand(os.path.join(PATH_SALMON, "quant_from_bam", "{fq_name}", "quant.sf"), fq_name = SAMPLES)
    output:
        deseq_output = os.path.join(PATH_DESEQ, config["analysis"]["sample"][0] + "_" + config["analysis"]["control"][0], "deseq_output.csv")
    params:
        control_ids = config["analysis"]["control"],
        sample_ids = config["analysis"]["sample"],
        count_table = os.path.join(PATH_DESEQ, config["analysis"]["sample"][0] + "_" + config["analysis"]["control"][0], "count_table.csv"),
        info_table = os.path.join(PATH_DESEQ, config["analysis"]["sample"][0] + "_" + config["analysis"]["control"][0], "info_table.csv"),
        replicates = config["analysis"]["replicate"],
        matched_samples = MATCHED_SAMPLES
    script:
        "scripts/run_deseq.R"



rule run_pathfindr:
    input:
        pathfindr_input = os.path.join(PATH_DESEQ, config["analysis"]["sample"][0] + "_" + config["analysis"]["control"][0], "deseq_output.csv")
    output:
        pathfindr_output = directory(PATH_PATHFINDR)
    params:
        threads = config["required"]["threads"],
        p_val_threshold = config["required"]["pathfindr"]["p_val_threshold"],
        iterations = config["required"]["pathfindr"]["iterations"],
        gene_name_converter = config["required"]["pathfindr"]["gene_name_converter"]
    script:
        "scripts/run_pathfindr.R"
