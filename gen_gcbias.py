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

def process_reads(partition, before_corr_reads, after_corr_reads):
	csv_obj = CsvObj("./gcbias_illumina/"+str(partition)+".csv")
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

		# if reads aligns reverse line[8] is larger than line[9] so switch them or this messes up indexing
		if int(line[8]) > int(line[9]):
			start = line[9]
			line[9] = line[8]
			line[8] = start
		kept_gc = get_read_gc_percentage(before_corr_reads[line[1]][(int(line[8])-1):(int(line[9])-1)])
		chucked_read = before_corr_reads[line[1]][0:(int(line[8])-1)] + before_corr_reads[line[1]][(int(line[9])-1):-1]
		chucked_gc = get_read_gc_percentage(chucked_read)
                output_dict = {'partition':partition, 'read_name':line[1],'kept_gc':kept_gc,'chucked_gc':chucked_gc}
                csv_obj.add_line_to_csv(output_dict)
        csv_obj.close_csv_handle()

def get_read_gc_percentage(read):
	if len(read) == 0:
		return None
	bases = {'GC': 0, 'others': 0}
	for char in read:
		if char == "G" or char == "C":
			bases['GC'] = bases['GC'] + 1
		else:
			bases['others'] = bases['others'] + 1
	gc_percentage = float(bases['GC'])*100 / (bases['GC']+bases['others'])
	return gc_percentage

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
