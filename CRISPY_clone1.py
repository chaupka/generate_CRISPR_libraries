#import os
import os.path
#import json
#from decimal import Decimal
from Bio import SeqIO
import csv
#import matplotlib.pyplot as plt
#import numpy as np
import re

dir = "CRISPRi_library1"
if not os.path.exists(dir):
    os.mkdir(dir)

genbank_list = ['MGA3_chromosome.gb','MGA3_pBM19.gb','MGA3_pBM69.gb']
seq_type_list = ['chromosome','pBM19','pBM69']
record_list = []
for subdir in seq_type_list:
    if not os.path.exists(dir+'/'+subdir):
        os.mkdir(dir+'/'+subdir)

newCSV = 'CRISPRi_library1_overview.csv'
feature_dict = {}
feature_type_list = ['CDS', 'rRNA', 'tRNA', 'misc_feature']
fieldnames = ['seq_type', 'locus_tag', 'product', 'strand', 'start', 'end', 'length', 'sgRNAs', 'sgRNAs with 2bp offtargets']
theor_PAM_list = ["AGG", "CGG", "GGG", "TGG"]
PAM_counter = 0
twobp_counter = 0
sgRNA_list = []
PAM_list = []
strand_list = []
twobp_mismatch_list = []
sgRNA_position_list = []

feature_counter = 0

for genbank_file in genbank_list:
    for gb in SeqIO.parse(genbank_file, 'gb'): #parse through gb file
        record_list.append(gb)
for record_number,record in enumerate(record_list):
    if record.features: #if gb file has features (it should)
        for featnum,feature in enumerate(record.features): #iterate through the features
            if featnum == 2: #limit for testing
                break
            if feature.type not in feature_type_list: #only specified features, skip duplicates
                continue
            else: #collect feature info
                seq_type = seq_type_list[record_number] #chromosome, plasmid...
                try:
                    feature.qualifiers['locus_tag']
                    locus = feature.qualifiers['locus_tag'][0]
                except:
                    continue
                try:
                    feature.qualifiers['product']
                    product = feature.qualifiers['product'][0]
                except:
                    product = ''
                try:
                    feature.location.strand
                    if feature.location.strand == 1:
                        strand = '+'
                    else:
                        strand = '-'
                except:
                    strand = ''
                try:
                    feature.location.start
                    start = feature.location.start + 1
                except:
                    continue
                try:
                    feature.location.end
                    end = feature.location.end
                except:
                    continue
                try:
                    len(feature)
                    length = len(feature)
                except:
                    continue
                if strand == '+' or strand == '':
                    NT = record.seq[start-1:end]
                    T = record.seq[start-1:end].reverse_complement()
                else:
                    NT = record.seq[start-1:end].reverse_complement()
                    T = record.seq[start-1:end]
                if len(T) <= 23:
                    continue
                else:
                    print("Computing target: "+locus)
                    for seq_number, seq in enumerate([NT, T]):
                        if seq_number == 0:
                            target_strand = "NT"
                            target_start = 0
                            target_end = len(seq)//2
                            increment = 1
                        else:
                            target_strand = "T"
                            target_start = len(seq)
                            target_end = len(seq)//2
                            increment = -1
                        for nuc_number in range(target_start, target_end, increment): #iterate through half of the nucleotides
                            PAM = seq[nuc_number:nuc_number+3].reverse_complement()
                            if PAM in theor_PAM_list:
                                sgRNA = seq[
                                        nuc_number + 3:nuc_number + 23].reverse_complement()  # sgRNA is upstream of PAM
                                if len(sgRNA) < 20:
                                    continue
                                print("Found PAM at position " + str(abs(nuc_number - target_start + 2)) + "/" + str(
                                    len(seq)) + " targeting strand " + target_strand)
                                offtarget = 0
                                for nucleotide in range(0,len(sgRNA)): #iterate through all nucleotides of sgRNA
                                    for item in record_list: #iterate through all gb-files
                                        fwd = re.findall(str(sgRNA)[:nucleotide]+r'\w'+str(sgRNA)[(nucleotide+1):], str(item.seq))
                                        rev = re.findall(str(sgRNA)[:nucleotide]+r'\w'+str(sgRNA)[(nucleotide+1):], str(item.seq.reverse_complement()))
                                        #find all the 0 and 1 bp mutated offtargets
                                        offtarget += len(fwd) + len(rev)
                                if offtarget >= 21: #if more than 1 offtarget -> check next sgRNA
                                    continue
                                offtarget = 0
                                for nucleotide in range(0, len(sgRNA)-1): #iterate through all nucleotides of sgRNA -1
                                    for nucleotide2 in range(nucleotide+1,len(sgRNA)): #iterate through all nucleotides of sgRNA -2
                                        for item in record_list: #iterate through all gb-files
                                            fwd = re.findall(str(sgRNA)[:nucleotide]+r'\w'+str(sgRNA)[nucleotide+1:nucleotide2]+r'\w'+str(sgRNA)[(nucleotide2+1):],str(item.seq))
                                            rev = re.findall(str(sgRNA)[:nucleotide]+r'\w'+str(sgRNA)[nucleotide+1:nucleotide2]+r'\w'+str(sgRNA)[(nucleotide2+1):],str(item.seq.reverse_complement()))
                                            # find all the 0,1 and 2 bp mutated offtargets
                                            offtarget += len(fwd)+len(rev)
                                if offtarget >= 201: #if more than 10 offtargets -> check next sgRNA
                                    continue
                                else:
                                    PAM_counter += 1
                                    PAM_list.append(str(PAM))
                                    sgRNA_list.append(str(sgRNA))
                                    strand_list.append(target_strand)
                                    twobp_mismatch_list.append(offtarget-190)
                                    if offtarget >= 191:
                                        twobp_counter += 1
                                    sgRNA_position_list.append(str(abs(nuc_number-target_start+5))+"-"+str(abs(nuc_number-target_start+24)))

                    with open(dir+'/'+seq_type_list[record_number]+'/'+locus+'_sgRNAs.csv', 'w', encoding='utf-8') as new: #save sgRNA-list under feature name
                        w = csv.writer(new)
                        w.writerow(["sgRNA", "PAM", "target strand", "2bp mismatches", "location_revcomp"])
                        for PAM_num, PAM_PAM in enumerate(PAM_list):
                            w.writerow(
                                [sgRNA_list[PAM_num], PAM_list[PAM_num], strand_list[PAM_num], twobp_mismatch_list[PAM_num],
                                 sgRNA_position_list[PAM_num]])
                        new.close()
                        print(locus+"_sgRNAs.csv saved.\n")

                feature_dict[feature_counter] = [seq_type, locus, product, strand, start, end, length, PAM_counter, twobp_counter]
                PAM_counter = 0
                twobp_counter = 0
                feature_counter += 1
                sgRNA_list = []
                PAM_list = []
                strand_list = []
                twobp_mismatch_list = []
                sgRNA_position_list = []

with open(dir+'/'+newCSV,'w', encoding = 'utf-8') as new:
    w = csv.writer(new)
    w.writerow(fieldnames)
    for featnum,feature in enumerate(feature_dict.keys()):
        w.writerow(feature_dict[feature])
    new.close()
print(newCSV+" saved.")
