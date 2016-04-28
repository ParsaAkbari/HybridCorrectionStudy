from Bio import SeqIO
import csv
import sys
# loop through all the partitions
#   what is the number of A T G C before and after error correction

# for each partition loop through all the reads
# 	for each read what is the number of A T G C before and after error correction


# 1. get_records, input: partition number
#    an array of seqio objects each of which is a read
# 2. make_dictionary_atgc_numbers, input: array of seqio objects
#    loops through the reads one by one passing them into dicationary make
#    create a dicationary, key is read name, value is a dictionary containing number of A T G and Cs

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

def process_reads(partition, before_corr_reads, after_corr_reads):
	print("Part "+partition)
	csv_obj = CsvObj("./prev_post_csv_output/5mer_"+str(partition)+".csv")
	# loop through lines of partition csv
	
	blast_f = open("./make_alignment_csv_output/"+str(partition)+".csv", 'r')
       	blast_f.readline()
	# loop through each read
       	for line in blast_f:
		# grab the indexes from the line
		line = line.split(",")
		if len(line) < 4:
			continue
		# get the correct read from the before/after_corr_reads dictionary using line[1] (read header)
		# splice the reads according to blast output
		fivemer_prev = before_corr_reads[line[1]][(int(line[8])-6):(int(line[8])-1)]
		fivemer_post = before_corr_reads[line[1]][int(line[9]):(int(line[9])+5)]
                output_dict = {'partition':partition, 'read_name':line[1], 'fivemer_prev': fivemer_prev, 'fivemer_post' : fivemer_post}
                csv_obj.add_line_to_csv(output_dict)
        csv_obj.close_csv_handle()

def output_line_to_csv(read_name, dict):
	csv_string = ','.join("%s=%r" % (key,val) for (key,val) in dict.iteritems())
	csv_string = read_name + "," + csv_string

class CsvObj:
	header = []
	csvwriter = ''
	csvfile = ''
	def __init__(self, filename):
		self.csvfile = open(filename, 'w+')
		self.csvwriter = csv.writer(self.csvfile, delimiter=',')
				
	def add_line_to_csv(self, dictionary):
		# if we don't have a header array create one using the keys of the dictionary
		if len(self.header) == 0:
			self.create_header_array(dictionary)
		csv_array = []
		for key in self.header:
			csv_array.append(dictionary[key])

		self.write_line_to_csv(csv_array)

	# grab the keys of the dictionary, these will be the header of the csv
	def create_header_array(self, dictionary):
		header = []
		for key in dictionary:
			header.append(key)
		self.write_line_to_csv(header)
		self.header = header

	# write the input line to the csv file
	def write_line_to_csv(self, csv_array):
		self.csvwriter.writerow(csv_array)

	def close_csv_handle(self):
		self.csvfile.close()


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
