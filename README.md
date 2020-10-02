# generate_CRISPR_libraries

Parse Bacillus methanolicus genome and native plasmids for PAM sequence in upstream half of all features. Add 20 bps sequence upstream of hits to sgRNA list. Drop sgRNAs with off-target hits >0 of query sequence with 0 & 1 mismatches and off-target hits >10 of query sequence with 2 mismatches.Vary query sequence and mutation sequence (2bp mismatches) to generate 3 libraries: Downstream 13 bps query & mutation sequence resulting in library 0, 20 bp query and downstream 13 bps mutation sequence resulting in library 1 and the whole 20 bps query & mutation sequence resulting in library 2.

# DOI
http://doi.org/10.5281/zenodo.4063436