from Bio import pairwise2
from Bio import SeqIO
import sys
import commands
import os
import csv
import math
import linecache
import pickle
import random
# get mean nanopore read length: 5188
#             OPEN: accuracy_compare.csv, FINAL COLUMN IS BEFORE CORRECT READ LENGTH, GET AVERAGE VALUE
# generate csv file of windows
# perform blast search using settings used to blast the corrected reads:
## ref_blast_params = {"reference": ref_file,
##                "cor_query": correct_fa,
##                "ref_blast_out": refblast_out}
##runfail("blastn -db {reference} -query {cor_query} -outfmt \"6 std qlen slen\" -evalue 1e-10 -reward 5 -penalty -4 -gapopen 8 -gapextend 6 -dust no -task blastn -out {ref_blast_out}".format(**ref_blast_params))
# for each blast output, you get index of start & end of read
# so: use modulus (since the window size is half of mean read length) to figure out the window number the read covers
#     it may cover multiple windows, go through csv file made above and +1 to coverage
# use linecache module to jump to window?

# loop through csv file, looking at ref file and make GC% calculation

# plot it in R

def run_ref_blast(partition, before_corr_reads, coverage_hash, only_include={}):
    # run through each read, run blast against reference
    # loop through each blast result, find windows it covers
    # +1 to the coverage for that window

    for key in before_corr_reads:
        if len(only_include) > 0 and key in only_include:
            continue
        # query id, subject id, % identity, alignment length, mismatches, gap opens,
	# q. start, q. end, s. start, s. end, evalue, bit score
        # -db /home/pa354/c_elegans_data/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
        blastn_command = """/home/pa354/ncbi-blast/bin/blastn \
                            -db /home/pa354/c_elegans_data/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
                            -query <(echo \"{query}\") \
                           -outfmt "6 std qlen slen" -evalue 1e-10 -reward 5
                            -penalty -4 -gapopen 8 -gapextend 6 -dust no -task blastn'
                         """.format(query=before_corr_reads[key])
        blastn_command = '/bin/bash -c \"'+blastn_command+'\" | sort -rn -k 4 | head -1'
        blastn_output = commands.getstatusoutput(blastn_command)[1].split("\n")
        if blastn_output[0] == '':
            # continue as no matches
            continue

        windows_already_covered = []
        for alignment in blastn_output:
            alignment = alignment.split("\t")
            if len(alignment) == 1:
                # continue as this alignment probably returned an error
                continue
            # 8 is start, 9 is end
            start_window = int(math.floor(int(alignment[8]) / 2594))
            end_window = int(math.ceil(int(alignment[9]) / 2594))
            
            if end_window > len(coverage_hash[alignment[1]])-1:
                end_window = len(coverage_hash[alignment[1]]) -1
            for window in range(start_window, end_window+1):
                if window not in windows_already_covered:
                    windows_already_covered.append(window)
                    coverage_hash[alignment[1]][window] = coverage_hash[alignment[1]][window] + 1
    
    return coverage_hash


# an array of seqio objects each of which is a read
def get_records(partition, extension):
        handle = open("/home/pa354/c_elegans_data/0001/p"+partition+extension, "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records

def make_reads_dict(records):
        dict = {}
        for record in records:
                record_name = record.id.replace("_consensus", "")
                dict[record_name] = record.seq
        return dict

coverage_hash = {
	"I" : [0] * 5810,
	"II" : [0] * 5890,
	"III" : [0] * 5313,
	"IV" : [0] * 6743,
	"MtDNA" : [0] * 5,
	"V" : [0] * 8066,
	"X" : [0] * 6830

}

f = open('/home/pa354/c_elegans_data/analyse_nanocorr_data/gcbias_pickle/0483.pickle', 'rb')
coverage_hash = pickle.load(f)

only_include = {}

#f = open('/home/pa354/c_elegans_data/analyse_nanocorr_data/low_accuracy.csv', 'rb')
#f.next()
#for row in f:
#    if len(row) < 5:
#    	continue
#    only_include[row.split(",")[1].replace("\"", "").replace('\n', "")] = 1


print(sys.argv)
for partition_number in range(int(sys.argv[1]),int(sys.argv[2])):
    	print(partition_number)
        partition = '{0:04}'.format(partition_number)
        before_corr_records = get_records(partition, "")
        before_corr_reads = make_reads_dict(before_corr_records)

        coverage_hash = run_ref_blast(partition, before_corr_reads, coverage_hash, only_include)
	# save coverage_hash at this point
        f = open('/home/pa354/c_elegans_data/analyse_nanocorr_data/gcbias_pickle/'+partition+'.pickle', 'wb')
    	pickle.dump(coverage_hash, f)

openpath = '/home/pa354/c_elegans_data/analyse_nanocorr_data/gcbias/'

for key in coverage_hash:
	f = open(openpath+key+".txt", 'w+')
	for window in coverage_hash[key]:
            f.write(str(window)+"\n")
    
