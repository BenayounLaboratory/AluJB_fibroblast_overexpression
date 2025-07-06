#!/bin/bash

# CHANGE TO DIRECTORY WITH FASTQ FILES AND RUN THIS SCRIPT FROM THERE
# FOR OPTIMAL PERFORMANCE, RUN AS SUDO, IN ORDER TO RESET MEMORY AFTER EACH STAR ALIGNMENT





########### DEFINE INPUTS (CHANGE AS NEEDED)

# Define the number of threads/processors
threads=12
processors=6

# Define genome indices
star_indices="/Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/homo_sapiens_gencode_v44/indices/"

# Define the extensions for the F and R paired end files
R1_extension_Trim=".fastq.gz"


########### TRIM WITH FASTP

#Make new directory for the trimmed reads and for the trimming report
mkdir trimmed_reads
mkdir trimmed_reads/fastp_reports

# Loop over sequencing files and run fastp.
for f in *$R1_extension_Trim; do

	# Run fastp
	/Users/juanb/Documents/Bioinformatic_Tools/fastp_0.20.1/fastp --in1 $f --out1 trimmed_reads/fastp_${f} --disable_quality_filtering --trim_front1 15 --trim_tail1 24 --cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 5 --cut_right_mean_quality 20 --length_required 36 --overrepresentation_analysis --thread $threads --compression 5 --html trimmed_reads/fastp_reports/${f}_report.html --json trimmed_reads/fastp_reports/${f}_report.json

done

#Change to directory with trimmed reads.
cd trimmed_reads


########### QC WITH FASTQC

# Go to directory with processed reads.
mkdir fastqc_reports

# Run fastqc on all .gz files (.gz for more flexibility compared to fq.gz or fastq.gz).
/Users/juanb/Documents/Bioinformatic_Tools/FastQC_0.11.9/fastqc *.gz --outdir fastqc_reports --noextract --threads $threads


########### READ ALIGNMENT WITH STAR

#Make directory to hold STAR alignments (BAM format).
mkdir STAR_alignments

# Loop over sequencing files and align with STAR.
for file in *$R1_extension_Trim; do

	# Define the paired file name and the prefix for the output
	output_prefix=${file%$R1_extension_Trim}

	# Run STAR
	/Users/juanb/Documents/Bioinformatic_Tools/STAR-2.7.3a/source/STAR --runThreadN $processors --genomeDir $star_indices --readFilesCommand gunzip -c --readFilesIn $file --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outFilterMismatchNoverLmax 0.04 --outFileNamePrefix STAR_alignments/$output_prefix
	purge


done # Close FOR loop running STAR alignment over all sequencing files.
