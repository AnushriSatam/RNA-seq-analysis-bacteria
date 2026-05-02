#!/bin/bash
set -euo pipefail

WORKDIR="anushri/rna_seq_analysis"
SORTMERNA="softwares/sortmerna_4.3.7/bin/sortmerna"
SORTMERNA_DB_DIR="rna-seq/software/silva"
GFF="anushri/cyno_genome/GCF_000009725.1_ASM972v1_genomic.gff"
GENOME="anushri/cyno_genome/GCF_000009725.1_ASM972v1_genomic.fna"

mkdir -p ${WORKDIR}/results/sortmerna
mkdir -p ${WORKDIR}/results/trimmed_fastqc
mkdir -p ${WORKDIR}/results/counts

# FastQC (Quality check)
for f in ${WORKDIR}/*.fastq.gz
do
    fastqc "$f" -o ${WORKDIR}/raw_fastqc -t 32
done

# fastp (trimming)
mkdir -p ${WORKDIR}/results/fastp

for r1 in ${WORKDIR}/*_R1_001.fastq.gz
do
    BASENAME=$(basename $r1 _R1_001.fastq.gz)

    r2=${WORKDIR}/${BASENAME}_R2_001.fastq.gz

    fastp \
        -i $r1 \
        -I $r2 \
        -o ${WORKDIR}/results/fastp/${BASENAME}_R1_trimmed.fastq.gz \
        -O ${WORKDIR}/results/fastp/${BASENAME}_R2_trimmed.fastq.gz \
        -h ${WORKDIR}/results/fastp/${BASENAME}_fastp.html \
        -j ${WORKDIR}/results/fastp/${BASENAME}_fastp.json \
        -q 30 \
        -l 50 \
        --detect_adapter_for_pe \
        -w 16
done

# Quality check after trimming

for f in ${WORKDIR}/results/fastp/*.fastq.gz
do
    fastqc "$f" -o ${WORKDIR}/results/trimmed_fastqc -t 32
done

# kraken-contamination screening

kraken_db="path_to/Kraken_PlusPFP_16_2025"
kraken_outdir="${WORKDIR}/results/kraken"
mkdir -p ${kraken_outdir}

for r1 in ${WORKDIR}/results/fastp/*_R1_trimmed.fastq.gz
do
    BASENAME=$(basename $r1 _R1_trimmed.fastq.gz)
    r2=${WORKDIR}/results/fastp/${BASENAME}_R2_trimmed.fastq.gz

    kraken2 \
        --db $kraken_db \
        --paired \
        --gzip-compressed \
        --threads 32 \
        --output ${kraken_outdir}/${BASENAME}_kraken2_output.txt \
        --report ${kraken_outdir}/${BASENAME}_kraken2_report.txt \
        --report-minimizer-data \
        $r1 $r2

    echo "Kraken2 done: ${BASENAME}"
done

# SortMeRNA — rRNA removal

for r1 in ${WORKDIR}/results/fastp/*_R1_trimmed.fastq.gz
do
    BASENAME=$(basename $r1 _R1_trimmed.fastq.gz)
    r2=${WORKDIR}/results/fastp/${BASENAME}_R2_trimmed.fastq.gz

    rm -rf /path_to/sortmerna_index/idx/kvdb/* 2>/dev/null || true
    rm -rf /path_to/sortmerna_index/idx/readb/* 2>/dev/null || true

    $SORTMERNA \
        --ref ${SORTMERNA_DB_DIR}/rfam-5s-database-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/rfam-5.8s-database-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-arc-16s-id95.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-arc-23s-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-bac-16s-id90.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-bac-23s-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-euk-18s-id95.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-euk-28s-id98.fasta \
        --reads $r1 \
        --reads $r2 \
        --aligned ${WORKDIR}/results/sortmerna/${BASENAME}_rRNA \
        --other ${WORKDIR}/results/sortmerna/${BASENAME}_mRNA \
        --fastx \
        --paired-in \
        --out2 \
        --threads 32 \
        --index 0 \
        --workdir /path_to/sortmerna_index/idx

    mv ${WORKDIR}/results/sortmerna/${BASENAME}_mRNA_fwd.fq.gz \
       ${WORKDIR}/results/sortmerna/${BASENAME}_R1_mRNA.fastq.gz

    mv ${WORKDIR}/results/sortmerna/${BASENAME}_mRNA_rev.fq.gz \
       ${WORKDIR}/results/sortmerna/${BASENAME}_R2_mRNA.fastq.gz

    echo "SortMeRNA done: ${BASENAME}"
done

# BOWTIE2 indexing and alignment

# indexing was done prior:
# bowtie2-build $GENOME syneno_index

mkdir -p ${WORKDIR}/results/bowtie_alignment
BOWTIE_IDX="${WORKDIR}/results/bowtie_index/syneno_index"
BOWTIE_OUTDIR="${WORKDIR}/results/bowtie_alignment"

for r1 in ${WORKDIR}/results/sortmerna/*_R1_mRNA.fastq.gz
do
    BASENAME=$(basename $r1 _R1_mRNA.fastq.gz)
    r2=${WORKDIR}/results/sortmerna/${BASENAME}_R2_mRNA.fastq.gz

    bowtie2 \
        -x $BOWTIE_IDX \
        -1 $r1 \
        -2 $r2 \
        -p 32 \
        --very-sensitive \
        2> ${BOWTIE_OUTDIR}/${BASENAME}_bowtie2.log \
        | samtools view -bS - \
        | samtools sort -@ 16 -o ${BOWTIE_OUTDIR}/${BASENAME}.sorted.bam

    samtools index ${BOWTIE_OUTDIR}/${BASENAME}.sorted.bam

    echo "Finished: ${BASENAME}"
done

# infer strandedness
BED="anushri/cyno_genome/synechocystis.bed"
OUTFILE="${WORKDIR}/results/strandedness_check.txt"

> ${OUTFILE}

echo "Starting strandedness inference for all samples..."
echo ""

for bam in ${WORKDIR}/results/bowtie_alignment/*.sorted.bam
do
    sample=$(basename $bam .sorted.bam)

    echo "=== ${sample} ===" | tee -a ${OUTFILE}

    infer_experiment.py \
        -r ${BED} \
        -i ${bam} \
        -s 200000 | tee -a ${OUTFILE}

    echo "" | tee -a ${OUTFILE}
done

echo "Done! Results saved to ${OUTFILE}"

# featureCounts — Quantification

BAM_FILES=${WORKDIR}/results/bowtie_alignment/*.sorted.bam

featureCounts \
    -T 16 \
    -a ${GFF} \
    -F GFF \
    -o ${WORKDIR}/results/counts/all_samples_counts.txt \
    -p \
    -s 2 \
    -t gene \
    -g locus_tag \
    -M \
    --fraction \
    --minOverlap 10 \
    ${BAM_FILES}

echo "featureCounts done!"
echo "Count matrix: ${WORKDIR}/results/counts/all_samples_counts.txt"
