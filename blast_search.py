from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import parse
from Bio import SeqIO

records = SeqIO.parse("./apoe.fas", "fasta")

PROGRAM = 'blastp'
DATABASE = 'nr'
for rec in records:
    # query NCBI Blast API
    xml_result  = qblast(PROGRAM, DATABASE, rec.seq)
    # Parse xml result 
    results = parse(xml_result)
    # iterate over each result
    for record in results:

        for alignment in record.alignments:
            print(alignment)