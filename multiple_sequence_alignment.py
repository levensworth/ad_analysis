from Bio import AlignIO

alignments = AlignIO.parse("./../Downloads/BioData/resampled.phy", "phylip")
for alignment in alignments:
    print(alignment)