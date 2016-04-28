import csv
import sys
import commands
from Bio import SeqIO
# arguments for blast
#    reference
#    query (do it as raw sequence)
#    output (remove)    

# an array of seqio objects each of which is a read
def get_records(partition, extension):
        handle = open("../0001/p"+partition+extension, "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records

def make_reads_dict(records):
        dict = {}
        for record in records:
                record_name = record.id.replace("_consensus", "")
                dict[record_name] = record.seq
        return dict

def get_after_correct_identities(partition):
    after_correct_identities = []
    handle = open("./after_corr_acc_compare/"+partition+".csv", "rU")
    for line in handle:
        line = line.split(",")
        after_correct_identities.append(
       	{"read_name":line[0].replace("_consensus", ""), "identity_after":line[2], "alignment_length_after":line[3], "read_length_after":line[12]})
    return after_correct_identities

def get_before_correct_identities_and_save(partition, before_corr_reads, after_corr_identities):
    # third output is identity
    header = ["partition", "read_name", "identity_after", "alignment_length_after",
	      "read_length_after", "identity_before", "alignment_length_before",
	       "read_length_before"]
    csvfile = open("./accuracy_compare/"+partition+".csv", 'w+')
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(header)
    for identity in after_corr_identities:
        curr_line = [partition]
        curr_line.append(identity["read_name"])
        curr_line.append(identity["identity_after"])
        curr_line.append(identity["alignment_length_after"])
	curr_line.append(identity["read_length_after"])
        blastn_command = """/home/pa354/ncbi-blast/bin/blastn \
			    -db /home/pa354/c_elegans_data/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
		       	    -query <(echo \"{query}\") \
			    -outfmt "6 std qlen slen\" -evalue 1e-10 -reward 5
		       	    -penalty -4 -gapopen 8 -gapextend 6 -dust no -task blastn'
		         """.format(query=before_corr_reads[curr_line[1]])

        blastn_command = '/bin/bash -c \"'+blastn_command+'\" | sort -rn -k 4 | head -1'
        blastn_output = commands.getstatusoutput(blastn_command)[1].split("\t")
	
        if len(blastn_output) == 12:
            curr_line.append(blastn_output[2])
            curr_line.append(blastn_output[3])
	    curr_line.append(len(before_corr_reads[curr_line[1]]))
            csvwriter.writerow(curr_line)
    csvfile.close()

def make_after_correct_identities(partition):
	# unix open partition file load by bitsize
	# split output by new line
	# does the read name already exist
	#     if not add it
	#     if does exist skip
	csvfile = open("./after_corr_acc_compare/"+str(partition)+".csv", 'w+')
	csvwriter = csv.writer(csvfile, delimiter=',')
	comm = "cat ../0001/p"+str(partition)+".blast6.r.refblast6.q | sort -rn -k 12"
	comm_output = commands.getstatusoutput(comm)[1].split("\n")
	exists = []
	for line in comm_output:
		line = line.split("\t")
		if line[0] not in exists:
			exists.append(line[0])
			csvwriter.writerow(line)
	csvfile.close()

print(sys.argv)
for partition_number in range(int(sys.argv[1]),int(sys.argv[2])):
        if partition_number % 20 == 0:
                print(partition_number)
        partition = '{0:04}'.format(partition_number)
        before_corr_records = get_records(partition, "")
        before_corr_reads = make_reads_dict(before_corr_records)
        make_after_correct_identities(partition)
        after_corr_identities = get_after_correct_identities(partition)
        get_before_correct_identities_and_save(partition, before_corr_reads, after_corr_identities)
