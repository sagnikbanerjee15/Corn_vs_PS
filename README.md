# Overview of the project

Maize lines H95 (susceptible to Puccinia sorghi) and Rp1-I (resistant to Puccinia sorghi) were inoculated with Puccinia sorghi (Iowa-2016 isolate) and studied over 7 timepoints and 4 biological replicates. For each sample, there were both control and infected plants. Hence, the total number of samples were 112 = (2 Cultivars \* 7 Timepoints \* 4 Replicates \* 2 Conditions). For S2, the samples were sequenced on a NovaSeq platform with 150 paired ended reads on 2 lanes.

# Location of raw data 

All data (along with relevant metadata) has been upload to Google drive. Please contact mgelmore@iastate.edu if you wish to obtain access to it. Genomes of both the species have been uploaded. (@manju - pls update with version number)

# Data description

## Dateset 1

Metadata file is in https://docs.google.com/spreadsheets/d/1nBV_qJTsBRhHHX6Sze0aFm41hZcUGzdV/edit?rtpof=true

Paired ended (150+150) - S2 (NovaSeq 6000) - Sequenced on 2 lanes (technical replicates)

Number of replicates - 4
Number of timepoints - 7
Number of genotypes - 2 (Susceptible & Resistance)
Number of treatmets - 2 

Total number of samples = 4 * 7 * 2 * 2 = 112 (RNASeq id 1 - 112)

Resting Spores (5 Biological Relicates)
Germinating Spores (5 Biological Relicates)
Samples 113 to 122

## Dataset 2

Metadata file is in https://docs.google.com/spreadsheets/d/1Wv3Na_XafV6C0DZPnSPUNrHGAc5NJuta/edit?usp=drive_web&ouid=113465502642554295866&rtpof=true

Paired ended (250+250) - SP (NovaSeq 6000) - Sequenced on 2 lanes (technical replicates)

Pooled all timepoints and replicates for H95 & Control - Sample no 1
Pooled all timepoints and replicates for Rp1-I & Control - Sample no 2
Pooled all timepoints and replicates for H95 & Infected - Sample no 3 
Pooled all timepoints and replicates for Rp1-I & Infected - Sample no 4
Resting spores-CR1-pooled all 5 replicates - Sample no 5
Germinated spores-CR1-pooled all 5 replicates - Sample no 6

RNA-Seq of 3 fungal genomes from the other project

Sample no. 10
Sample no. 11
Sample no. 12

Soybean Rust

Sample no. 7
Sample no. 8
Sample no. 9

## Genomic References

Maize - AGPv4 B73
Puccinia sorghi (IA16) - Assembled by Katerina Holan 
Puccinia sorghi (Argetina isolate)

# Software requirements

* Singularity
* Python

All other software installations will be provided through singularity. Nothing else needs to be installed.

# Plan of action

## 1. a. Designing a pipeline to detect differentially expressed genes from RNA-Seq data. For this run, we used a combined reference of B73 AGPv4 and Argentina PS isolate. Once FINDER returns the gene models for the IA16 isolate, we will repeat this process all over again with combined reference comprising of B73 AGPv4 and IA16 PS isolate.

Deliver an end-to-end pipeline (written in python and deployed with singularity/docker). The pipeline will perform the following tasks:

* Download reads from NCBI-SRA or provided link or use pre-downloaded samples [Done]
* Preparing reference indices [Done]
* STAR/bowtie2 to reference [Done] - used bowtie2
* Salmon execution to generate the counts [Done]
* Accumulate counts from all the salmon executions [Done]
* Apply the DESeq2 model

## 1b. Perform differential gene analysis between the resting spores and the germinating spores using 5 biological replicates from S2 against IA16.

## 2. Verify the fungal reference genome [Done]

Verification was perfomed by mapping RNA-Seq data (samples 113-122) from S2 and Sample 5 & 6 from SP to the genome assembled by K Holan. We found that salsa_phase_hic was the best.

## 3. Annotate the reference genomes with FINDER [Deliver by End July]

* Preparing reference indices - B73 AGPv4 and IA16 [Done]
* Align reads to reference using STAR [Done]
* Use Spades to denovo assemble unmapped reads [Done]
* Use CD-HIT to cluster all the assembled transcripts [Done for thresholds - 0.95 & 0.90, 0.85 & 0.80]
* Rerun CD-HIT on representative transcripts collected from each cluster
* Run BLAST with NT and discard transcripts that have no hits
* Run PsiClass on the mapped data and obtain gene models
* Prepare final annotation files (GTF and FASTA) after removing Plant reads. Fasta file will contain transcripts from both mapped and unmapped reads.

## 4. Perform final clean up and deliver the results


# Workflow

1. Execute the program `mergeDataFromTwoLanesAndRenameRawFiles` to merge data from 2 lanes. Please note that this program will work only for this project. The program assumes that the same sample will have the same name in both the lanes and that the samples are paired ended. It uses the filenames to merge them. Then it will rename the file to something short. Please note that the renaming convention was chosen to be meaningful for this project only. Please review the code and make meaningful changes to it if you wish to reuse it for other projects. DO NOT execute as is.

2. The reference consisted of transcripts from B73 and Argentina PS isolate. Short reads were aligned to the merged reference using bowtie2 and salmon was used to estimate gene counts. Custom script was written up to accumulate the gene counts from salmon output. DESeq2 will be used to determine the genes that are differentially expressed. No trimmining of reads was performed

# PS IA16 Genome Annotation

## Adapter trimming

Trim adapters using Trimmomatic. No quality trimming was performed. Trimmomatic was executed with a clipping parameter of 2:30:10, that prompted Trimmomatic to allow for 2 mismatches. In case of paired end reads, trimmomatic will look for a score of 30, which is about 50 bases.

## Preparation of Reference 

Genomic reference for Zea mays version 4 was used (AGPv4 B73)
Fungal reference genome was assembled using nanopore reads (ref reqd.)

## Alignment to reference

STAR version 2.7.9a 
Both genomic references were merged. 
Adapter trimmed reads from all 128 samples were mapped. Each execution produced a mapped file in bam format and unmapped reads were output also.

## De novo assembly

Spades version 3.15.4 was used with kmer length set to 101. 
Unmapped reads from all 128 samples were merged together. No other de novo assemblers were explored to reduce complexity. Final assembled transcripts are in /90daydata/maizegdb/sagnik/CORN_VS_PS/PS_SALSA_annotation/denovo_assembly/all_unmapped_merged/.
Spades includes the gene number and isoform number in the name of the assembled transcript, which were later used to map transcript to genes

## Verification of de novo assembled transcripts

Trimmed short reads from all 128 samples were remapped to the de novo assembly using bowtie2 version 0.7.17 (123,792 transripts)
Number of reads mapping to each de novo assembled transcript was calculated and normalized over the length of the transcript, generating the coverage.
Transcripts were retained if they had a coverage of at least 1 and had an amino acid of at least 50 peptides in one of the 6 ORFs. (27098 gene and 37,811 transcripts) [Out of 123K transcripts, only 37,967 transcripts had a coverage > 1]

## Genome guided assembly

Extract read alignments that map to the PS genome using samtools and custom python script.
Generate genome guided assembly using PsiCLASS and also with Stringtie2
Merged annotations from PsiCLASS and Stringtie2 (13043 genes and 32836 transcripts) using gffcompare
These are all PS IA16 genes and transcripts

## Counts generation (Salmon)

Transcripts from PS_IA16 and B73v4 were merged together (106384 genes and 228888 transcripts)
Read counts were generated using salmon


