#!/bin/bash

for file in sam_filtered/*
do

    basename=${file##*/}
    sampleID=${basename%.sam}

    ## Print mapped segments to BAM format
    samtools view -bS -F 4 -o bam/$sampleID.bam $file

    ## Sort alignments by leftmost coordinates
    samtools sort bam/$sampleID.bam bam/$sampleID.sorted

    # Save as SAM
    samtools view -h -o sam_sorted/$sampleID.sam bam/$sampleID.sorted.bam 

    echo $sampleID
    
done
