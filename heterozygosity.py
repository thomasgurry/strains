"""

Python module for analysing heterozygosity of individual SNPs,
i.e. probability of sharing the same allele on two overlapping reads.
Can be used to compute the heterozygosity WITHIN a sample, as well as
BETWEEN two samples.  A few notes about each use case:

WITHIN: computed as the probability of two distinct, overlapping reads agreeing on a given SNP allele.
   Counts each read pair only once, and does not include pairs where a read is compared to itself.
   Hz = 1.0 - number_of_read_pairs_agreeing/total_number_of_read_pairs

BETWEEN: computed as the product of probabilities of same alleles in each sample, e.g. 
   if sample A = [0.4, 0, 0, 0.6] and B = [0.5, 0.5, 0, 0],
   Hz = 1.0 - 0.4*0.5

Usage:

1.  (WITHIN SAMPLE HETEROZYGOSITY) 

        i.  Compute SNP loci from folder of kpileups for all samples:

            python heterozygosity.py -l True -k /path/kpileups

        ii. Run each sample individually from a SORTED sam file.  Note the sam file must be sorted or the routine will not work:

            python heterozygosity.py -i input_sorted.sam -s sampleID -o output_dir 


2.  (BETWEEN SAMPLE HETEROZYGOSITY)

          i.  Compute SNP loci from folder of kpileups for all samples:

              python heterozygosity.py -l True -k /path/kpileups

         ii.  Compute all start lines for each gene in each SAM file for step (iii) speed-up (note SAM files must be sorted as in WITHIN sample heterozygosity routine):

              python heterozygosity.py -bw True -p True -i /path/to/sam/folders -o /path/to/precomputed/startlines

        iii.  Run each pair of samples as SORTED sam files of interest, either using a multiprocessing pool (function: parse_two_sam_files_multiprocessing()), or as a script as follows:

              python heterozygosity.py -bw True -i1 input_sorted_1.sam -i2 input_sorted_2.sam -ip /path/to/precomputed/startlines -s1 sampleID_1 -s2 sampleID_2 -o output_dir

"""

import os, glob
import numpy as np
import cPickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-bw', help='Between/within flag (True/False; True = compute heterozygosity between subjects)', default='False')
parser.add_argument('-i', help='Input sam file (for within sample heterozygosity), or if used with -bw True and -p True, input folder containing all samples for which to precompute start lines')
parser.add_argument('-l', help='Compute SNP loci flag (True/False)')
parser.add_argument('-k', help='Input kpileups folder')
parser.add_argument('-s', help='Sample ID')
parser.add_argument('-sl', help='Path to SNP loci pickled file from step 1 (Default: snp_loci.pkl in current directory)')
parser.add_argument('-i1', help='Input sam file 1 (use with -bw True)')
parser.add_argument('-i2', help='Input sam file 2 (use with -bw True)')
parser.add_argument('-ip', help='Input folder containing precomputed startline dicts (use with -bw True)')
parser.add_argument('-s1', help='Sample ID 1 (use with -bw True)')
parser.add_argument('-s2', help='Sample ID 2 (use with -bw True)')
parser.add_argument('-p', help='Precompute start lines for each gene in each sample (True/False; use with -bw True)', default='False')
parser.add_argument('-o', help='output directory')
args = parser.parse_args()

def build_snp_positions_folder(kpileups_folder):
    # Parses kpileups output files and extracts SNP loci
    # Returns a dict of genes with SNP loci positions for each gene
    filenames = glob.glob(os.path.join(kpileups_folder, '*.txt'))
    print filenames
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
    heterozygosity_dict = {}
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
                if bact not in heterozygosity_dict:
                    heterozygosity_dict[bact] = {}
                if gene_name not in heterozygosity_dict[bact]:
                    heterozygosity_dict[bact][gene_name] = {}
                heterozygosity_dict[bact][gene_name][snp_position_of_interest] = 1.0 - same_allele_counts/float(overlapping_read_counts)
            else:
                continue
                
    return heterozygosity_dict, readcount_dict


def precompute_samfile_startlines(samfile):
    # Precompute start lines for each gene in a sorted SAM file
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
        if i % 100000 == 0:
            print str(i) + ' of ' + str(len(all_lines))

    # Write results to pickle file
    results_dict = {'genes_counted': genes_counted, 'gene_firstlines': gene_firstlines}
    return results_dict
    
    
def parse_two_sam_files(samfile1, samfile2, snp_loci, precomputed_startlines_folder):
    # Parses two sorted alignments SAM files for overlapping positions
    # -for each read, extracts gene it maps to, alignment position start and finish
    with open(samfile1, 'r') as fid:
        all_lines_1 = fid.readlines()
    with open(samfile2, 'r') as fid:
        all_lines_2 = fid.readlines()
    
    # Get start lines for each gene in both files
    samfile1_basename = samfile1.split('/')[-1]
    samfile2_basename = samfile2.split('/')[-1]
    startlines_dict_1 = cPickle.load(open(os.path.join(precomputed_startlines_folder, samfile1_basename + '.pkl'), 'rb'))
    startlines_dict_2 = cPickle.load(open(os.path.join(precomputed_startlines_folder, samfile2_basename + '.pkl'), 'rb'))
    gene_firstlines_1 = startlines_dict_1['gene_firstlines']
    genes_counted_1 = startlines_dict_1['genes_counted']
    gene_firstlines_2 = startlines_dict_2['gene_firstlines']
    genes_counted_2 = startlines_dict_2['genes_counted']

    # Look for read indices of genes of interest
    good_indices_1 = {}
    good_indices_2 = {}
    genes_of_interest = snp_loci.keys()
    for gene in genes_of_interest:
        if gene in genes_counted_1:
            good_indices_1[gene] = genes_counted_1.index(gene)
        if gene in genes_counted_2:
            good_indices_2[gene] = genes_counted_2.index(gene)

    # Extract overlapping read pairs for each gene, and compute total read count per gene
    heterozygosity_dict = {}
    readcount_dict_1 = {}
    readcount_dict_2 = {}
    counter_1 = 0

    # Loop through each gene in SAM file 1 and find equivalent gene in SAM file 2
    for gene in good_indices_1:

        gene_line_index_1 = good_indices_1[gene]

        try:
            gene_line_index_2 = good_indices_2[gene]

            counter_1 += 1

            print "Computing " + str(counter_1) + " of " + str(len(good_indices_1))
            geneID_1 = genes_counted_1[gene_line_index_1]
            geneID_2 = genes_counted_2[gene_line_index_2]

            if geneID_1 == 'PN:bwa':
                continue
            if geneID_2 == 'PN:bwa':
                continue

            bact_1 = '_'.join(geneID_1.split(';')[0].split('_')[:2])
            gene_name_1 = geneID_1.split(';')[0].split('_')[2]
            firstline_1 = gene_firstlines_1[gene_line_index_1]
            if counter_1 < len(gene_firstlines_1)-1:
                lastline_1 = gene_firstlines_1[gene_line_index_1 + 1]
            else:
                lastline_1 = len(all_lines_1)
            # Compute total number of reads for that gene
            if bact_1 not in readcount_dict_1:
                readcount_dict_1[bact_1] = {}
            if gene_name_1 not in readcount_dict_1[bact_1]:
                readcount_dict_1[bact_1][gene_name_1] = lastline_1 - firstline_1 + 1

            bact_2 = '_'.join(geneID_2.split(';')[0].split('_')[:2])
            gene_name_2 = geneID_2.split(';')[0].split('_')[2]
            firstline_2 = gene_firstlines_2[gene_line_index_2]
            if counter_1 < len(gene_firstlines_2)-1:
                lastline_2 = gene_firstlines_2[gene_line_index_2 + 1]
            else:
                lastline_2 = len(all_lines_2)
            # Compute total number of reads for that gene
            if bact_2 not in readcount_dict_2:
                readcount_dict_2[bact_2] = {}
            if gene_name_2 not in readcount_dict_2[bact_2]:
                readcount_dict_2[bact_2][gene_name_2] = lastline_2 - firstline_2 + 1   

            # Check that both bacteria and genes are the same
            if bact_1 != bact_2 or gene_name_1 != gene_name_2:
                raise NameError('Gene IDs ' + str(geneID_1) + ' and ' + str(geneID_2) + ' do not match!')

            # Get read positions on the gene
            filtered_reads = {}
            for snp_position_of_interest in snp_loci[geneID_1]:

                read_positions_1 = {}
                read_positions_2 = {}

                # SAM file 1
                for j in range(firstline_1, lastline_1):
                    line = all_lines_1[j]
                    line = line.rstrip().split()
                    if len(line) > 3:
                        readID = line[0]
                        start_position = int(line[3])
                        end_position = start_position + len(line[9]) - 1  # mapped sequence = line[9]; not allocated for speed purposes
                        read_positions_1[readID] = [start_position, end_position, j]

                # SAM file 2
                for j in range(firstline_2, lastline_2):
                    line = all_lines_2[j]
                    line = line.rstrip().split()
                    if len(line) > 3:
                        readID = line[0]
                        start_position = int(line[3])
                        end_position = start_position + len(line[9]) - 1  # mapped sequence = line[9]; not allocated for speed purposes
                        read_positions_2[readID] = [start_position, end_position, j]

                # Get all genotypes and only consider this SNP further if there are at least 2 counts of each allele
                readIDs_1 = read_positions_1.keys()
                readIDs_2 = read_positions_2.keys()
                snp_genotypes_1 = []
                snp_genotypes_2 = []
                for m in range(len(readIDs_1)):
                    readID = readIDs_1[m]
                    if read_positions_1[readID][0] <= snp_position_of_interest and read_positions_1[readID][1] >= snp_position_of_interest:
                        read_line = all_lines_1[read_positions_1[readID][2]].rstrip().split()
                        read = read_line[9]
                        read_start = int(read_line[3])
                        read_snp_index =  snp_position_of_interest - read_start
                        snp_genotypes_1.append(read[read_snp_index])
                for m in range(len(readIDs_2)):
                    readID = readIDs_2[m]
                    if read_positions_2[readID][0] <= snp_position_of_interest and read_positions_2[readID][1] >= snp_position_of_interest:
                        read_line = all_lines_2[read_positions_2[readID][2]].rstrip().split()
                        read = read_line[9]
                        read_start = int(read_line[3])
                        read_snp_index =  snp_position_of_interest - read_start
                        snp_genotypes_2.append(read[read_snp_index])
                alleles_1 = np.unique(snp_genotypes_1)
                alleles_2 = np.unique(snp_genotypes_2)

                # If more than one allele, check for counts and only include alleles with more than 1 count.
                good_alleles_1 = []
                good_alleles_2 = []
                if len(alleles_1) > 0:
                    for allele in alleles_1:
                        if snp_genotypes_1.count(allele) > 1:
                            good_alleles_1.append(allele)
                if len(alleles_2) > 0:
                    for allele in alleles_2:
                        if snp_genotypes_2.count(allele) > 1:
                            good_alleles_2.append(allele)

                # Case 1: both have only one 'good' allele that passes the filter
                if len(good_alleles_1) == 1 and len(good_alleles_2) == 1:
                    if bact_1 not in heterozygosity_dict:
                        heterozygosity_dict[bact_1] = {}
                    if gene_name_1 not in heterozygosity_dict[bact_1]:
                        heterozygosity_dict[bact_1][gene_name_1] = {}
                    if good_alleles_1[0] == good_alleles_2[0]:
                        heterozygosity_dict[bact_1][gene_name_1][snp_position_of_interest] = 0.0
                    else:
                        heterozygosity_dict[bact_1][gene_name_1][snp_position_of_interest] = 1.0

                # Case 2: either one or both samples have more than one allele that pass the filter
                if len(good_alleles_1) > 1 or len(good_alleles_1) > 1:
                    if bact_1 not in heterozygosity_dict:
                        heterozygosity_dict[bact_1] = {}
                    if gene_name_1 not in heterozygosity_dict[bact_1]:
                        heterozygosity_dict[bact_1][gene_name_1] = {}

                    allele_freqs_1 = {allele: snp_genotypes_1.count(allele)/float(len(snp_genotypes_1)) for allele in good_alleles_1}
                    allele_freqs_2 = {allele: snp_genotypes_2.count(allele)/float(len(snp_genotypes_2)) for allele in good_alleles_2}

                    probability_of_agreement = 0.0
                    for allele in allele_freqs_1:
                        if allele in allele_freqs_2:
                            probability_of_agreement += allele_freqs_1[allele] * allele_freqs_2[allele]
                    heterozygosity_dict[bact_1][gene_name_1][snp_position_of_interest] = 1.0 - probability_of_agreement

                # Case 3: one or both samples have no allele; do not compute heterozygosity
                if len(good_alleles_1) == 0 or len(good_alleles_1) == 0:
                    continue
        except:
            continue
        
    return heterozygosity_dict, readcount_dict_1, readcount_dict_2


def parse_two_sam_files_multiprocessing(samfile1, samfile2, snp_loci, precomputed_startlines_folder, sampleID1, sampleID2, results_dict):
    # Parses two sorted alignments SAM files for overlapping positions
    # Used for multiprocessing calls using pool.apply_async() from a multiprocessing pool.  Requires 'results_dict' to be a manager dict(), i.e. global across the pool.

    # -for each read, extracts gene it maps to, alignment position start and finish
    with open(samfile1, 'r') as fid:
        all_lines_1 = fid.readlines()
    with open(samfile2, 'r') as fid:
        all_lines_2 = fid.readlines()
    
    # Get start lines for each gene in both files
    samfile1_basename = samfile1.split('/')[-1]
    samfile2_basename = samfile2.split('/')[-1]
    startlines_dict_1 = cPickle.load(open(os.path.join(precomputed_startlines_folder, samfile1_basename + '.pkl'), 'rb'))
    startlines_dict_2 = cPickle.load(open(os.path.join(precomputed_startlines_folder, samfile2_basename + '.pkl'), 'rb'))
    gene_firstlines_1 = startlines_dict_1['gene_firstlines']
    genes_counted_1 = startlines_dict_1['genes_counted']
    gene_firstlines_2 = startlines_dict_2['gene_firstlines']
    genes_counted_2 = startlines_dict_2['genes_counted']

    # Look for read indices of genes of interest
    good_indices_1 = {}
    good_indices_2 = {}
    genes_of_interest = snp_loci.keys()
    for gene in genes_of_interest:
        if gene in genes_counted_1:
            good_indices_1[gene] = genes_counted_1.index(gene)
        if gene in genes_counted_2:
            good_indices_2[gene] = genes_counted_2.index(gene)

    # Extract overlapping read pairs for each gene, and compute total read count per gene
    heterozygosity_dict = {}
    readcount_dict_1 = {}
    readcount_dict_2 = {}
    counter_1 = 0

    # Loop through each gene in SAM file 1 and find equivalent gene in SAM file 2
    for gene in good_indices_1:

        gene_line_index_1 = good_indices_1[gene]

        try:
            gene_line_index_2 = good_indices_2[gene]

            counter_1 += 1

            print "Computing " + str(counter_1) + " of " + str(len(good_indices_1))
            geneID_1 = genes_counted_1[gene_line_index_1]
            geneID_2 = genes_counted_2[gene_line_index_2]

            if geneID_1 == 'PN:bwa':
                continue
            if geneID_2 == 'PN:bwa':
                continue

            bact_1 = '_'.join(geneID_1.split(';')[0].split('_')[:2])
            gene_name_1 = geneID_1.split(';')[0].split('_')[2]
            firstline_1 = gene_firstlines_1[gene_line_index_1]
            if counter_1 < len(gene_firstlines_1)-1:
                lastline_1 = gene_firstlines_1[gene_line_index_1 + 1]
            else:
                lastline_1 = len(all_lines_1)
            # Compute total number of reads for that gene
            if bact_1 not in readcount_dict_1:
                readcount_dict_1[bact_1] = {}
            if gene_name_1 not in readcount_dict_1[bact_1]:
                readcount_dict_1[bact_1][gene_name_1] = lastline_1 - firstline_1 + 1

            bact_2 = '_'.join(geneID_2.split(';')[0].split('_')[:2])
            gene_name_2 = geneID_2.split(';')[0].split('_')[2]
            firstline_2 = gene_firstlines_2[gene_line_index_2]
            if counter_1 < len(gene_firstlines_2)-1:
                lastline_2 = gene_firstlines_2[gene_line_index_2 + 1]
            else:
                lastline_2 = len(all_lines_2)
            # Compute total number of reads for that gene
            if bact_2 not in readcount_dict_2:
                readcount_dict_2[bact_2] = {}
            if gene_name_2 not in readcount_dict_2[bact_2]:
                readcount_dict_2[bact_2][gene_name_2] = lastline_2 - firstline_2 + 1   

            # Check that both bacteria and genes are the same
            if bact_1 != bact_2 or gene_name_1 != gene_name_2:
                raise NameError('Gene IDs ' + str(geneID_1) + ' and ' + str(geneID_2) + ' do not match!')

            # Get read positions on the gene
            filtered_reads = {}
            for snp_position_of_interest in snp_loci[geneID_1]:

                read_positions_1 = {}
                read_positions_2 = {}

                # SAM file 1
                for j in range(firstline_1, lastline_1):
                    line = all_lines_1[j]
                    line = line.rstrip().split()
                    if len(line) > 3:
                        readID = line[0]
                        start_position = int(line[3])
                        end_position = start_position + len(line[9]) - 1  # mapped sequence = line[9]; not allocated for speed purposes
                        read_positions_1[readID] = [start_position, end_position, j]

                # SAM file 2
                for j in range(firstline_2, lastline_2):
                    line = all_lines_2[j]
                    line = line.rstrip().split()
                    if len(line) > 3:
                        readID = line[0]
                        start_position = int(line[3])
                        end_position = start_position + len(line[9]) - 1  # mapped sequence = line[9]; not allocated for speed purposes
                        read_positions_2[readID] = [start_position, end_position, j]

                # Get all genotypes and only consider this SNP further if there are at least 2 counts of each allele
                readIDs_1 = read_positions_1.keys()
                readIDs_2 = read_positions_2.keys()
                snp_genotypes_1 = []
                snp_genotypes_2 = []
                for m in range(len(readIDs_1)):
                    readID = readIDs_1[m]
                    if read_positions_1[readID][0] <= snp_position_of_interest and read_positions_1[readID][1] >= snp_position_of_interest:
                        read_line = all_lines_1[read_positions_1[readID][2]].rstrip().split()
                        read = read_line[9]
                        read_start = int(read_line[3])
                        read_snp_index =  snp_position_of_interest - read_start
                        snp_genotypes_1.append(read[read_snp_index])
                for m in range(len(readIDs_2)):
                    readID = readIDs_2[m]
                    if read_positions_2[readID][0] <= snp_position_of_interest and read_positions_2[readID][1] >= snp_position_of_interest:
                        read_line = all_lines_2[read_positions_2[readID][2]].rstrip().split()
                        read = read_line[9]
                        read_start = int(read_line[3])
                        read_snp_index =  snp_position_of_interest - read_start
                        snp_genotypes_2.append(read[read_snp_index])
                alleles_1 = np.unique(snp_genotypes_1)
                alleles_2 = np.unique(snp_genotypes_2)

                # If more than one allele, check for counts and only include alleles with more than 1 count.
                good_alleles_1 = []
                good_alleles_2 = []
                if len(alleles_1) > 0:
                    for allele in alleles_1:
                        if snp_genotypes_1.count(allele) > 1:
                            good_alleles_1.append(allele)
                if len(alleles_2) > 0:
                    for allele in alleles_2:
                        if snp_genotypes_2.count(allele) > 1:
                            good_alleles_2.append(allele)

                # Case 1: both have only one 'good' allele that passes the filter
                if len(good_alleles_1) == 1 and len(good_alleles_2) == 1:
                    if bact_1 not in heterozygosity_dict:
                        heterozygosity_dict[bact_1] = {}
                    if gene_name_1 not in heterozygosity_dict[bact_1]:
                        heterozygosity_dict[bact_1][gene_name_1] = {}
                    if good_alleles_1[0] == good_alleles_2[0]:
                        heterozygosity_dict[bact_1][gene_name_1][snp_position_of_interest] = 0.0
                    else:
                        heterozygosity_dict[bact_1][gene_name_1][snp_position_of_interest] = 1.0

                # Case 2: either one or both samples have more than one allele that pass the filter
                if len(good_alleles_1) > 1 or len(good_alleles_1) > 1:
                    if bact_1 not in heterozygosity_dict:
                        heterozygosity_dict[bact_1] = {}
                    if gene_name_1 not in heterozygosity_dict[bact_1]:
                        heterozygosity_dict[bact_1][gene_name_1] = {}

                    allele_freqs_1 = {allele: snp_genotypes_1.count(allele)/float(len(snp_genotypes_1)) for allele in good_alleles_1}
                    allele_freqs_2 = {allele: snp_genotypes_2.count(allele)/float(len(snp_genotypes_2)) for allele in good_alleles_2}

                    probability_of_agreement = 0.0
                    for allele in allele_freqs_1:
                        if allele in allele_freqs_2:
                            probability_of_agreement += allele_freqs_1[allele] * allele_freqs_2[allele]
                    heterozygosity_dict[bact_1][gene_name_1][snp_position_of_interest] = 1.0 - probability_of_agreement

                # Case 3: one or both samples have no allele; do not compute heterozygosity
                if len(good_alleles_1) == 0 or len(good_alleles_1) == 0:
                    continue
        except:
            continue

    sample_pair = sampleID1 + '-' + sampleID2
    results_dict[sample_pair] = {'Hz': heterozygosity_dict, 'readcount_dict_1': readcount_dict_1, 'readcount_dict_2': readcount_dict_2}
    
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
    sample_snp_counts = {} # counts for each SNP when reads overlap at that SNP (snp_consensus/snp_counts = heterozygosity)
    for gene in snp_dict:
        sample_snp_consensus[gene] = np.zeros((1, len(snp_dict[gene])))
        sample_snp_counts[gene] = np.zeros((1, len(snp_dict[gene])))

    # Extract overlapping read pairs for each gene, and compute total read count per gene
    heterozygosity_dict = {}
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
            if bact not in heterozygosity_dict:
                heterozygosity_dict[bact] = {}
            heterozygosity_dict[bact][gene_name] = same_allele_counts/float(overlapping_read_counts)
        else:
            continue
        
    # Parse individual SNP heterozygosities
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
    return heterozygosity_dict, readcount_dict, snp_hz_dict, snp_dict


# MAIN
if __name__ == "__main__":

    # Get polymorphic loci
    if args.l == 'True':
        print "Computing SNP loci from kpileups folder " + str(args.k)
        snp_loci = build_snp_positions_folder(args.k)
        cPickle.dump(snp_loci, open('snp_loci.pkl','wb'))
    elif args.bw == 'False':
        print "Computing SNP heterozygosities for sample " + str(args.s)
        if args.sl:
            snp_loci = cPickle.load(open(args.sl,'rb'))
        else:
            snp_loci = cPickle.load(open('snp_loci.pkl','rb'))
        hz_dict, readcount_dict = parse_sam_file(args.i, snp_loci)
        if not os.path.exists(args.o):
            os.system('mkdir ' + str(args.o))
        cPickle.dump(hz_dict, open(os.path.join(args.o, args.s + '_hz.pkl'),'wb'))
        cPickle.dump(readcount_dict, open(os.path.join(args.o, args.s + '_readcounts.pkl'),'wb'))
    elif args.bw == 'True':
        if args.p == 'True':
            print "Precomputing gene start lines in all SAM files for between sample heterozygosity calculations."
            samfiles = os.listdir(args.i)
            if not os.path.exists(args.o):
                os.system('mkdir ' + str(args.o))
            startline_counter = 0
            for filename in samfiles:
                results_dict = precompute_samfile_startlines(os.path.join(args.i, filename))
                cPickle.dump(results_dict, open(os.path.join(args.o, filename + '.pkl'), 'wb'))
                startline_counter += 1
                print "Computed start lines for " + str(startline_counter) + " of " + str(len(samfiles)) + " files."

        elif args.p == 'False':
            print "Computing SNP heterozygosities between samples " + str(args.s1) + " and " + str(args.s2)
            if args.sl:
                snp_loci = cPickle.load(open(args.sl,'rb'))
            else:
                snp_loci = cPickle.load(open('snp_loci.pkl','rb'))
            hz_dict, readcount_dict_1, readcount_dict_2 = parse_two_sam_files(args.i1, args.i2, snp_loci, args.ip)
            if not os.path.exists(args.o):
                os.system('mkdir ' + str(args.o))
            cPickle.dump(hz_dict, open(os.path.join(args.o, args.s1 + '_' + args.s2 + '_hz.pkl'),'wb'))
            cPickle.dump(readcount_dict_1, open(os.path.join(args.o, args.s1 + '_readcounts.pkl'),'wb'))
            cPickle.dump(readcount_dict_1, open(os.path.join(args.o, args.s2 + '_readcounts.pkl'),'wb'))
