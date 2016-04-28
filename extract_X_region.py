from Bio import pairwise2
from Bio import SeqIO
import pprint
pp = pprint.PrettyPrinter(indent=4)

# an array of seqio objects each of which is a read
def get_records():
        handle = open("/home/pa354/c_elegans_data/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa", "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records

for record in get_records():
	if record.id == 'X':
		X_chromosome = record.seq
		break



f = open('x_region.txt', 'w+')
f.write('>xregion')
f.write(str(X_chromosome[(1045*2594):(1055*2594)]))
f.close()
