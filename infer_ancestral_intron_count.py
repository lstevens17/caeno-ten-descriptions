#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from Bio import SeqIO, AlignIO
from os.path import basename, isfile, abspath, join
import sys

# Assumptions:
# Main assumption is that this scripts assumes that the intron information is contained within a correct format.
# The files we have used are variable, column-length TSV files (I apologise).
# We aim to rewrite these scripts to parse BED/GFF3 file, but when analysing individual orthogroups, we've found the following method to be the quickest/easiest to debug.

# Step 1: Parse alignment using AlignIO (BioPython), returning an alignment object
def parse_alignment(alignment_file):
    with open(alignment_file, 'r') as alignment:
        alignment_object = AlignIO.read(alignment, "fasta")
        return alignment_object

# Step 2: Parse intron information file (custom file created for each Orthogroup using Ensembl Perl API), returning a dict
def parse_introns(intron_position_file):
    with open(intron_position_file, 'r') as intron_positions:
        intron_position_dict = {}
        for line in intron_positions:
            ID = line.rstrip("\n").split("\t")[0]
            list = line.rstrip("\n").split("\t")[1:]
            intron_position_dict[ID] = list
    return intron_position_dict

# Step 3: Convert CDS intron positions to alignment intron positions
def convert_to_relative_intron_positions(alignment_object, intron_position_dict):
    intron_alignment_position_dict = {}
    for record in alignment_object: # for every record in alignment
        ID = record.id # retrieve the ID
        aligned_sequence = record.seq # retrieve the aligned sequence
        intron_position_list = intron_position_dict[ID] # retrieve the intron positions
        base_count, alignment_position, intron_count = 0, 0, 0 # initialise counts
        if len(intron_position_list) == 1 and intron_position_list[0] == '': # if the sequence has no introns
	    intron_alignment_position_dict[ID] = ['-']
	else: # if the sequence has > 0 introns
            for position in aligned_sequence: # for every position in the aligned sequence (including dashes)
                if intron_count < len(intron_position_list): # stop once you've accounted for every intron
                    if not position == '-': # if the position is a base
                        base_count += 1 # increase base_count by 1
                        alignment_position += 1 # increase alignment_position by 1
                        if base_count == int(intron_position_list[intron_count]): # if the base count matches an intron position
                            try:
                                intron_alignment_position_dict[ID].append(alignment_position) # add the alignment position to the list
                            except KeyError:
                                intron_alignment_position_dict[ID] = [alignment_position] # or create the list if it doesn't exist yet
                            intron_count += 1
                    else:
                            alignment_position += 1
    return intron_alignment_position_dict

# Step 4: Declare orthologous introns based on conserved position in alignment
def declare_orthologous_introns(intron_alignment_position_dict):
    presence_absence_dict = {}
    # Step 4.1: Create list of all sites in alignment that contain at least one intron, and a sorted list of sequences
    all_intron_site_list = []
    sequenceID_list = []
    for ID, intron_alignment_position_list in intron_alignment_position_dict.iteritems():
        for intron_alignment_position in intron_alignment_position_list:
            if not intron_alignment_position in all_intron_site_list and not intron_alignment_position == '-':
                all_intron_site_list.append(intron_alignment_position)
        sequenceID_list.append(ID)
    # Step 4.2: For every sequence, ask whether it has a intron present at each intron containing site
    for ID in sorted(sequenceID_list):
        presence_absence_list = []
        present_intron_list = intron_alignment_position_dict[ID]
        for intron_site in all_intron_site_list:
            if intron_site in present_intron_list:
                presence_absence_list.append("1")
            else:
                presence_absence_list.append("0")
        #print ID + "\t" + "\t".join(presence_absence_list)
        presence_absence_dict[ID] = presence_absence_list
    return presence_absence_dict

def find_ancestral_introns(presence_absence_dict):
    for ID, presence_absence_list in presence_absence_dict.iteritems():
        length = len(presence_absence_list)
    ancestral_count = 0
    for i in range(0, length): # for every intron position in list
        basal_flag = 0
        ingroup_flag = 0
        for ID, presence_absence_list in presence_absence_dict.iteritems(): # for every species
            if ID.startswith("CMONO") or ID.startswith("DCORO"): # if there is an intron present at the position in C. monodelphis or Diploscapter coronatus
                if presence_absence_list[i] == '1':
                    basal_flag += 1 # increase basal flag by 1
            else:
                if presence_absence_list[i] == '1': # or if there is an intron present at that position in any other species
                    ingroup_flag += 1 # increase the ingroup flag by 1
        if basal_flag > 0 and ingroup_flag > 0: # if any intron is present in either (or both) of the basal taxa and at least one ingroup taxa
            ancestral_count += 1 # increase count of ancestral introns by 1
        elif basal_flag == 2: # or if intron is present in both the outgroups
            ancestral_count += 1 # increase count of ancestral introns by 1
    print ancestral_count # print the number of ancestral introns

if __name__ == "__main__":

    SCRIPT = basename(__file__)

    aln_file,intron_file = '',''
    try:
        alignment_file = sys.argv[1]
        intron_position_file = sys.argv[2]
    except:
        sys.exit("USAGE: ./%s %s %s" % (SCRIPT, "ALIGNMENT", "INTRONS"))

    alignment_object = parse_alignment(alignment_file)
    intron_position_dict = parse_introns(intron_position_file)
    intron_alignment_position_dict = convert_to_relative_intron_positions(alignment_object, intron_position_dict)
    presence_absence_dict = declare_orthologous_introns(intron_alignment_position_dict)
    find_ancestral_introns(presence_absence_dict)
