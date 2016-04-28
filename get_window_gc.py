from Bio import pairwise2
from Bio import SeqIO
import math
import csv
def make_reads_dict(records):
        dict = {}
        for record in records:
                record_name = record.id.replace("_consensus", "")
                dict[record_name] = record.seq
        return dict

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

gc_hash = {
        "I" : [0] * 5810,
        "II" : [0] * 5890,
        "III" : [0] * 5313,
        "IV" : [0] * 6743,
        "MtDNA" : [0] * 5,
        "V" : [0] * 8066,
        "X" : [0] * 6830
}

handle = open("/home/pa354/c_elegans_data/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa", "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

chromosomes = make_reads_dict(records)
# 2594
for key in chromosomes:
    csvfile = open("/home/pa354/c_elegans_data/analyse_nanocorr_data/gcbias/"+key+"_gc_percentage.csv", 'w+')
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(["chromosome","window","gc_percentage"])
    number_of_windows = int(math.floor(len(chromosomes[key])/2594))
    for window in range(1,number_of_windows +1):
        window_start = (window - 1) * 2594
        window_end = window * 2594
        gc_percentage = get_read_gc_percentage(chromosomes[key][window_start:window_end])
        csvwriter.writerow([key, window, gc_percentage])
    csvfile.close()
