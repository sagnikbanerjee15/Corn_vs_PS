# Overview of the project

Maize lines H95 (susceptible to Puccinia sorghi) and Rp1-I (resistant to Puccinia sorghi) were inoculated with Puccinia sorghi (Iowa-2016 isolate) and studied over 7 timepoints and 4 biological replicates. For each sample, there were both control and infected plants. Hence, the total number of samples were 112 = (2 Cultivars \* 7 Timepoints \* 4 Replicates \* 2 Conditions). For S2, the samples were sequenced on a NovaSeq platform with 150 paired ended reads on 2 lanes.

# Location of raw data 

All data (along with relevant metadata) has been upload to Google drive. Please contact mgelmore@iastate.edu if you wish to obtain access to it. Genomes of both the species have been uploaded. (@manju - pls update with version number)

# Software requirements

* Singularity
* Python

All other software installations will be provided through singularity. Nothing else needs to be installed.

# Plan of action

1. Designing a pipeline to detect differentially expressed genes from RNA-Seq data

Deliver an end-to-end pipeline (written in python and deployed with singularity/docker). The pipeline will perform the following tasks:

* Download reads from NCBI-SRA or provided link or use pre-downloaded samples
* Trimmomatic to remove adapters
* Preparing reference indices
* STAR/bowtie2 to reference
* Salmon execution to generate the counts
* Accumulating counts from all the salmon executions
* Applying the DESeq2 model

2. Verify the fungal reference genome

3. Annotate the reference genomes with FINDER

4. Perform final clean up and deliver the results


# Workflow

1. Execute the program `mergeDataFromTwoLanesAndRenameRawFiles` to merge data from 2 lanes. Please note that this program will work only for this project. The program assumes that the same sample will have the same name in both the lanes and that the sampoles are paired ended. It uses the filenames to merge them. Then it will rename the file to something short. Please note that the renaming convention was chosen to be meaningful for this project only. Please review the code and make meaningful changes to it if you wish to reuse it for other projects. DO NOT execute as is.

2. 