from Bio import pairwise2
from Bio import SeqIO
import sys
import commands
import os
import csv
# load after_corr sequence
# load before_corr sequence
# perform alignment
# investigate output, find indexes for start & end
#       the indexes must correspond to the before_corr read
# extract the bit of before_corr read which matches
# save this to file

def process_reads(partition, before_corr_reads, after_corr_reads):
    print(partition)

    csvfile = open("./make_alignment_csv_output/"+partition+".csv", 'w+')
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow("query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score".split(","))
    for read in after_corr_reads:
        blastn_command = """/home/pa354/ncbi-blast/bin/blastn \
	            	 -query <(echo \"{query}\") \
	                 -subject <(echo \"{subject}\") \
	                 -word_size 11 -outfmt 10 | sort -g -k 11 | tail -n1
		      	 """.format(query=after_corr_reads[read], subject=before_corr_reads[read])
        blastn_command = '/bin/bash -c \"'+blastn_command+'\"'
        blastn_csv_output = commands.getstatusoutput(blastn_command)[1]
        blastn_csv_output = blastn_csv_output.replace("Query_1", read+"_consensus").replace("Subject_1", read)
        csvwriter.writerow(blastn_csv_output.split(","))


def make_reads_dict(records):
        dict = {}
        for record in records:
                record_name = record.id.replace("_consensus", "")
                dict[record_name] = record.seq
        return dict
# an array of seqio objects each of which is a read                                                                                                                    
def get_records(partition, extension):
        handle = open("../0001/p"+partition+extension, "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records



print(sys.argv)

for partition_number in range(int(sys.argv[1]),int(sys.argv[2])):
        if partition_number % 20 == 0:
                print(partition_number)
        partition = '{0:04}'.format(partition_number)
        before_corr_records = get_records(partition, "")
        after_corr_records = get_records(partition, ".blast6.r.fa")
        before_corr_reads = make_reads_dict(before_corr_records)
        after_corr_reads = make_reads_dict(after_corr_records)
        process_reads(partition, before_corr_reads, after_corr_reads)
