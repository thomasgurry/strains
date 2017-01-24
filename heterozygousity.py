"""

Python module for analysing heterozygousity of individual SNPs,
i.e. probability of sharing the same allele on two overlapping reads.

Usage:

1.  Compute SNP loci from folder of kpileups for all samples:

        python heterozygousity.py -l True -k /path/kpileups

2.  Run each sample individually from a SORTED sam file.  Note the sam file must be sorted or the routine will not work.:

        python heterozygousity.py -i input_sorted.sam -s sampleID -o output_dir 

"""

import os, glob
import numpy as np
import cPickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input sam file')
parser.add_argument('-l', help='Compute SNP loci flag (True/False)')
parser.add_argument('-k', help='Input kpileups folder')
parser.add_argument('-s', help='Sample ID')
parser.add_argument('-o', help='output directory')
args = parser.parse_args()

def build_snp_positions_folder(kpileups_folder):
    # Parses kpileups output files and extracts SNP loci
    # Returns a dict of genes with SNP loci positions for each gene
    filenames = glob.glob(os.path.join(kpileups_folder, '*.txt'))
    snp_loci = {}
    for fn in filenames:
        print fn
        for line in open(fn):
            line = line.rstrip().split()
            if len(line) == 10 and line[0] != 'Sample':
                sample = line[0]
                position = line[2]
                gene = line[3]
                allele_freq = float(line[9])
                if allele_freq < 100.0:
                    if gene not in snp_loci:
                        snp_loci[gene] = []
                        if position not in snp_loci[gene]:
                            snp_loci[gene].append(int(position))
    return snp_loci

def build_snp_positions_onefile(kpileups_file):
    # Parses a single kpileups output file and extracts SNP loci
    # Returns a dict of genes with SNP loci positions for each gene
    snp_loci = {}
    for line in open(kpileups_file):
        line = line.rstrip().split()
        if len(line) == 10 and line[0] != 'Sample':
            sample = line[0]
            position = line[2]
            gene = line[3]
            if gene not in snp_loci:
                snp_loci[gene] = []
            if position not in snp_loci[gene]:
                snp_loci[gene].append(int(position))
    return snp_loci


def check_overlap(positions_1, positions_2):
    # Function to compare indicies for [start_position_1, finish_position_1] and 
    # [start_position_2, finish_position_2] to check for overlap.
    # Returns indices of overlap.
    if positions_1[1] < positions_2[0]:
        return None
    elif positions_2[1] < positions_1[0]:
        return None
    else:
        if positions_1[0] <= positions_2[0]:
            if positions_1[1] < positions_2[1]:
                overlapping_indices = [positions_2[0], positions_1[1]]
            elif positions_1[1] >= positions_2[1]:
                overlapping_indices = [positions_2[0], positions_2[1]]
        elif positions_2[0] < positions_1[0]:
            if positions_2[1] < positions_1[1]:
                overlapping_indices = [positions_1[0], positions_2[1]]
            elif positions_2[1] >= positions_1[1]:
                overlapping_indices = [positions_1[0], positions_1[1]]
    return overlapping_indices

def parse_sam_file(samfile, snp_loci):
    # Parses a sorted alignments SAM file for overlapping positions
    # -for each read, extracts gene it maps to, alignment position start and finish
    with open(samfile, 'r') as fid:
        all_lines = fid.readlines()

    # Get start lines for each gene
    gene_firstlines = []
    genes_counted = []
    for i in range(len(all_lines)):
        line = all_lines[i]
        line = line.rstrip().split()
        if len(line) > 3:
            geneID = line[2]
            if geneID not in genes_counted:
                genes_counted.append(geneID)
                gene_firstlines.append(i)

    # Look for read indices of genes of interest
    good_indices = []
    genes_of_interest = snp_loci.keys()
    for gene in genes_of_interest:
        if gene in genes_counted:
            good_indices.append(genes_counted.index(gene))
    
    # Extract overlapping read pairs for each gene, and compute total read count per gene
    heterozygousity_dict = {}
    readcount_dict = {}
    counter = 0
    for i in good_indices:
        counter += 1
        print "Computing " + str(counter) + " of " + str(len(good_indices))
        geneID = genes_counted[i]
        if geneID == 'PN:bwa':
            continue
        bact = '_'.join(geneID.split(';')[0].split('_')[:2])
        gene_name = geneID.split(';')[0].split('_')[2]
        firstline = gene_firstlines[i]
        if i < len(gene_firstlines)-1:
            lastline = gene_firstlines[i+1]
        else:
            lastline = len(all_lines)
        # Compute total number of reads for that gene
        if bact not in readcount_dict:
            readcount_dict[bact] = {}
        if gene_name not in readcount_dict[bact]:
            readcount_dict[bact][gene_name] = lastline - firstline + 1
        
        # Get read positions on the gene
        filtered_reads = {}
        for snp_position_of_interest in snp_loci[geneID]:
            read_positions = {}
            for j in range(firstline, lastline):
                line = all_lines[j]
                line = line.rstrip().split()
                if len(line) > 3:
                    readID = line[0]
                    start_position = int(line[3])
                    end_position = start_position + len(line[9]) - 1  # mapped sequence = line[9]; not allocated for speed purposes
                    read_positions[readID] = [start_position, end_position, j]

            # Get all genotypes and only consider this SNP further if there are at least 2 counts of each allele
            readIDs = read_positions.keys()
            snp_genotypes = []
            for m in range(len(readIDs)):
                readID = readIDs[m]
                if read_positions[readID][0] <= snp_position_of_interest and read_positions[readID][1] >= snp_position_of_interest:
                    read_line = all_lines[read_positions[readID][2]].rstrip().split()
                    read = read_line[9]
                    read_start = int(read_line[3])
                    read_snp_index =  snp_position_of_interest - read_start
                    snp_genotypes.append(read[read_snp_index])
            alleles = np.unique(snp_genotypes)
            # If multiple alleles, check that each has at least two counts; if single allele, keep 
            if len(alleles) > 1:
                if snp_genotypes.count(alleles[0]) > 1 and snp_genotypes.count(alleles[1]) > 1:
                    pass
                else:
                    continue
            else:
                continue
                
            # Check for reads with SNP of interest
            same_allele_counts = 0
            overlapping_read_counts = 0
            good_readIDs = []
            for m in range(len(readIDs)):
                readID = readIDs[m]
                if read_positions[readID][0] <= snp_position_of_interest and read_positions[readID][1] >= snp_position_of_interest:
                    good_readIDs.append(readID)

            # Check for overlap
            for m in range(len(good_readIDs)):
                readID = good_readIDs[m]
                read1_line = all_lines[read_positions[readID][2]].rstrip().split()
                read1 = read1_line[9]
                read1_start = int(read1_line[3])
                read1_snp_index =  snp_position_of_interest - read1_start
                for n in range(m+1, len(good_readIDs)):
                    readID_2 = good_readIDs[n]
                    read2_line = all_lines[read_positions[readID_2][2]].rstrip().split()
                    read2 = read2_line[9]
                    read2_start = int(read2_line[3])
                    read2_snp_index = snp_position_of_interest - read2_start
                    if read1[read1_snp_index] == read2[read2_snp_index]:
                        same_allele_counts += 1
                    overlapping_read_counts += 1

            if overlapping_read_counts > 0:
                if bact not in heterozygousity_dict:
                    heterozygousity_dict[bact] = {}
                if gene_name not in heterozygousity_dict[bact]:
                    heterozygousity_dict[bact][gene_name] = {}
                heterozygousity_dict[bact][gene_name][snp_position_of_interest] = 1.0 - same_allele_counts/float(overlapping_read_counts)
            else:
                continue
                
    return heterozygousity_dict, readcount_dict


def parse_sam_file_individual_snps(samfile, kpileups_file):
    # Parses a sorted alignments SAM file for overlapping positions
    # -for each read, extracts gene it maps to, alignment position start and finish
    with open(samfile, 'r') as fid:
        all_lines = fid.readlines()

    # Get start lines for each gene
    gene_firstlines = []
    genes_counted = []
    for i in range(len(all_lines)):
        line = all_lines[i]
        line = line.rstrip().split()
        if len(line) > 3:
            geneID = line[2]
            if geneID not in genes_counted:
                genes_counted.append(geneID)
                gene_firstlines.append(i)

    # Get SNP positions and build counts arrays for each
    snp_dict = build_snp_positions_onefile(kpileups_file)
    sample_snp_consensus = {} # counts for each SNP where overlapping reads agree on allele
    sample_snp_counts = {} # counts for each SNP when reads overlap at that SNP (snp_consensus/snp_counts = heterozygousity)
    for gene in snp_dict:
        sample_snp_consensus[gene] = np.zeros((1, len(snp_dict[gene])))
        sample_snp_counts[gene] = np.zeros((1, len(snp_dict[gene])))

    # Extract overlapping read pairs for each gene, and compute total read count per gene
    heterozygousity_dict = {}
    readcount_dict = {}
    for i in range(len(gene_firstlines)):
        print "Computing " + str(i) + " of " + str(len(gene_firstlines))
        read_positions = {}
        geneID = genes_counted[i]
        if geneID == 'PN:bwa':
            continue
        bact = '_'.join(geneID.split(';')[0].split('_')[:2])
        gene_name = geneID.split(';')[0].split('_')[2]
        firstline = gene_firstlines[i]
        if i < len(gene_firstlines)-1:
            lastline = gene_firstlines[i+1]
        else:
            lastline = len(all_lines)
        # Compute total number of reads for that gene
        if bact not in readcount_dict:
            readcount_dict[bact] = {}
        if gene_name not in readcount_dict[bact]:
            readcount_dict[bact][gene_name] = lastline - firstline + 1
        # Get read positions on the gene    
        for j in range(firstline, lastline):
            line = all_lines[j]
            line = line.rstrip().split()
            if len(line) > 3:
                readID = line[0]
                start_position = int(line[3])
                end_position = start_position + len(line[9]) - 1  # mapped sequence = line[9]; not allocated for speed purposes
                read_positions[readID] = [start_position, end_position, j]
        
        # Check for overlap
        same_allele_counts = 0
        overlapping_read_counts = 0
        readIDs = read_positions.keys()
        for m in range(len(readIDs)):
            readID = readIDs[m]
            for n in range(m+1, len(readIDs)):
                readID_2 = readIDs[n]
                overlapping_inds = check_overlap(read_positions[readID][:2], read_positions[readID_2][:2])
                if overlapping_inds != None:
                    read1_line = all_lines[read_positions[readID][2]].rstrip().split()
                    read2_line = all_lines[read_positions[readID_2][2]].rstrip().split()
                    read1 = read1_line[9]
                    read2 = read2_line[9]
                    read1_start = int(read1_line[3])
                    read2_start = int(read2_line[3])
                    overlapping_inds_read1 = np.subtract(overlapping_inds, read1_start)
                    overlapping_inds_read2 = np.subtract(overlapping_inds, read2_start)
                    read1_fragment = read1[overlapping_inds_read1[0]:overlapping_inds_read1[1]+1]
                    read2_fragment = read2[overlapping_inds_read2[0]:overlapping_inds_read2[1]+1]
                    if read1_fragment == read2_fragment:
                        same_allele_counts += 1
                    overlapping_read_counts += 1

                    # Check each SNP
                    try:
                        snp_inds = np.intersect1d(range(overlapping_inds[0], overlapping_inds[1]+1), snp_dict[geneID])
                        if snp_inds != []:
                            read1_snp_inds = np.subtract(snp_inds, read1_start)
                            read2_snp_inds = np.subtract(snp_inds, read2_start)
                            for p in range(len(read1_snp_inds)):
                                read1_allele = read1[read1_snp_inds[p]]
                                read2_allele = read2[read2_snp_inds[p]]
                                if read1_allele == read2_allele:
                                    sample_snp_consensus[geneID][0,p] += 1
                                sample_snp_counts[geneID][0,p] += 1
                    except:
                        continue
        if overlapping_read_counts > 0:
            if bact not in heterozygousity_dict:
                heterozygousity_dict[bact] = {}
            heterozygousity_dict[bact][gene_name] = same_allele_counts/float(overlapping_read_counts)
        else:
            continue
        
    # Parse individual SNP heterozygousities
    snp_hz_dict = {}
    for gene in sample_snp_consensus:
        if np.sum(sample_snp_counts[gene]) > 0:
            snp_counts = sample_snp_counts[gene]
            snp_consensus = sample_snp_consensus[gene]
            # Filter counts that don't have an overlapping read pair that have the same allele
            hz_vals = np.divide(sample_snp_consensus[gene], sample_snp_counts[gene])
            comparison_vals = np.divide(float(1), sample_snp_counts[gene])
            hz_vals = hz_vals[hz_vals > comparison_vals]
            if ~np.all(np.isnan(hz_vals)):
                snp_hz_dict[gene] = hz_vals 
    return heterozygousity_dict, readcount_dict, snp_hz_dict, snp_dict

# MAIN
#hz_dict, readcount_dict = parse_sam_file(args.i)


# Get polymorphic loci
if args.l == 'True':
    print "Computing SNP loci from kpileups folder " + str(args.k)
    snp_loci = build_snp_positions_folder(args.k)
    cPickle.dump(snp_loci, open('snp_loci.pkl','wb'))
else:
    print "Computing SNP heterozygousities for sample " + str(args.s)
    snp_loci = cPickle.load(open('snp_loci.pkl','rb'))
    hz_dict, readcount_dict = parse_sam_file(args.i, snp_loci)
    cPickle.dump(hz_dict, open(args.o + '/' + args.s + '_hz.pkl','wb'))
    cPickle.dump(readcount_dict, open(args.o + '/' + args.s + '_readcounts.pkl','wb'))
