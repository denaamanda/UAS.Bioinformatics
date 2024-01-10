#!/bin/env bash

# Project: UAS Bioinformatics
# Date: 5 Jan 2024
# Last updated: 5 Jan 2024
# Author: Dena Amanda Ashar - 21611135
# Description: This is a workflow to analyze FASTA file

# Global variables
WORKDIR=$PWD/UAS_Workflow
REF_GENOME_DIR="${WORKDIR}/Reference_Genome"
FASTQ_DIR="${WORKDIR}/Fastq_Data"
TRIMMED_DIR="${WORKDIR}/Trimmed_Data"
RESULTS_DIR="${WORKDIR}/Results"
QC_DIR="${RESULTS_DIR}/QC"
REF_GENOME="RWA2_reference_genome_Dnoxia_1.0-1_20150113.fa"
REF_GENOME_ARCHIVE="RWA2_reference_genome_Current_Release.tgz"
REF_GENOME_URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_001186385.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT"
 
FASTA_R1_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6032457/suppl/GSM6032457_Trinity.fasta.gz"

# Analysis steps
# Step 1: Download Data
# Step 2: Preprocessing data
# Step 3: Trimming
# Step 4: Squence alignment
# Step 5: Calculate statistics about he alignment
# Step 6: Variant calling


# Creating directories
mkdir -p "${WORKDIR}" "${REF_GENOME_DIR}" "${FASTQ_DIR}" "${TRIMMED_DIR}" "${RESULTS_DIR}" "${QC_DIR}"
cd "${WORKDIR}"

# Step 1: Download data
# Downloading the reference genome if it does not exist or is empty
if [ ! -s "${REF_GENOME_DIR}/${REF_GENOME_ARCHIVE}" ]; then
    wget -P "${REF_GENOME_DIR}" $REF_GENOME_URL
    tar -xzvf "${REF_GENOME_DIR}/${REF_GENOME_ARCHIVE}" -C "${REF_GENOME_DIR}"
else
    echo "Reference genome archive already exists and is not empty."
fi

# Downloading FASTQ files if they do not exist or are empty
if [ ! -s "${FASTA_DIR}/$(basename ${FASTA_R1_URL})" ]; then
    wget -P "${FASTA_DIR}" "${FASTA_R1_URL}"
else
    echo "FASTQ R1 file already exists and is not empty."
fi

# Step 2: Preprocessing data using Quality Control with FastQC
fastqc "${FASTA_DIR}"/*.fastq.gz -o "${QC_DIR}"

# Step 3: Read Trimming using Trimmomatic
trimmomatic PE -phred33 \
"${FASTA_DIR}/$(basename ${FASTA_R1_URL})"\
"${TRIMMED_DIR}/trimmed_R1.fasta.gz" "${TRIMMED_DIR}/unpaired_R1.fasta.gz" \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 4: Squence alignment (Read Alignment)
bwa index "${REF_GENOME_DIR}/${REF_GENOME}"
bwa mem "${REF_GENOME_DIR}/${REF_GENOME}" "${TRIMMED_DIR}/trimmed_R1.fastq.gz" "${TRIMMED_DIR}/trimmed_R2.fastq.gz" > "${RESULTS_DIR}/aligned_reads.sam"

# Step 5: Calculate statistics about he alignment (Post-Alignment Processing)
samtools view -bS "${RESULTS_DIR}/aligned_reads.sam" | samtools sort -o "${RESULTS_DIR}/sorted_aligned_reads.bam"
picard MarkDuplicates I="${RESULTS_DIR}/sorted_aligned_reads.bam" O="${RESULTS_DIR}/marked_duplicates.bam" M="${RESULTS_DIR}/marked_dup_metrics.txt"

# Step 6: Variant Calling
freebayes -f "${REF_GENOME_DIR}/${REF_GENOME}" "${RESULTS_DIR}/marked_duplicates.bam" > "${RESULTS_DIR}/variants.vcf"



