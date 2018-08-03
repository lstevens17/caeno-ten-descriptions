# Intron loss analysis

Scripts and data files associated with the intron loss analysis

`cds_alignments`: Directory containing 990 CDS alignments 

`intron_positions`: Directory containing intron position information for 990 orthogroups

`Orthogroups_990_0missing_5hetero_paraloguesremoved_withCDSalignments.txt`: Single copy orthogroups.txt file (note: this is after removing redundant sequences arising from uncollapsed 
heterozygosity and removing orthogroups containing paralogues)

`infer_ancestral_intron_count.py`: Python script to infer ancestral intron counts (take a single CDS alignment and intron position file as input, and prints the inferred ancestral count)

