# SeqStrings

This Python script searches for a set of patterns in a database of strings.
Set up for fastq.gz files.

Input files format is:
Depend/mutstrings3.csv

Left side 30 bases,Mutation (Alt base),Right Side 30 bases,Mutation information
CAGACGGAAACCGTAGCTGCCCTGGTAGGT,T,TTCTGGGAAGGGACAGAAGATGACAGGGGC,TP53:chr17:7579384T>A
CGATGGTCATGCACCGCGTGGGTCACAAAC,G,AACGGACGTACCCAAAGGAGAGTTTGCACA,ZFHX3:chr16:g.72993158G>A

Depend/fastq_files.csv
sample1.fastq.gz,sample2.fastq.gz
