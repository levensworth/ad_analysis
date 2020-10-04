import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument(
    "-f", "--file", type=str, help="path to the genbank file", required=True
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    help="[OPTIONAL] path to the fasta output",
    required=False,
)
args = parser.parse_args()

# record = SeqIO.read(args['file'], "genbank")
# reads mRNA sequences in genbanl file
record = SeqIO.parse("./../Downloads/BioData/apoe_2_rna.gb", "genbank")

table = 11
min_pro_len = 100


def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1

    answer.sort()
    return answer


orf_list = []

for s in record:
    seq = s.seq
    seq = seq.transcribe()
    orf_list += find_orfs_with_trans(seq, table, min_pro_len)

sequences = [SeqRecord(Seq(pro), id="APOE") for start, end, strand, pro in orf_list]

try:
    output_file_direction = args.output if args.output else os.path.basename(args.file)
    with open(output_file_direction, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    print(
        "output sequences ins fasta format can be found in: {}".format(
            output_file_direction
        )
    )
    print("DONE")
except Exception as e:
    print("ERROR! {}".format(str(e)))
