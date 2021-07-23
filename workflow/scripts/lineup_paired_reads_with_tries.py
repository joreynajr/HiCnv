import pysam
import argparse
import pygtrie

# parsing the command line interface
parser = argparse.ArgumentParser()

parser.add_argument('--r1', required=True, type=str, help='R1 input file.')
parser.add_argument('--r2', required=True, type=str, help='R2 input file.')
parser.add_argument('--r1-out', required=True, type=str, help='R1 output file.')
parser.add_argument('--r2-out', required=True, type=str, help='R2 output file.')
params = parser.parse_args()

# This trie implementation can result in huge memory consumption 
# so I'll be making the elements all equal to the none object below
# and therefore the memory won't increase dramatically
none = None

# load the same files
print("# load the same files")
sam1 = pysam.AlignmentFile(params.r1, "rb")
sam2 = pysam.AlignmentFile(params.r2, "rb")

# getting the read names for the first BAM/SAM
print("# getting the read names for the first BAM/SAM")
rnames1 = pygtrie.Trie()
for i, aln in enumerate(sam1.fetch()):
    rnames1[aln.query_name] = None

# getting the common names while scanning the second BAM/SAM
print("# getting the read names for the second BAM/SAM")
common = pygtrie.Trie()
for i, aln in enumerate(sam2.fetch()):
    if aln.query_name in rnames1:
        common[aln.query_name] = None
    
# write the reads containing the common names for the first BAM/SAM
print("# write the reads containing the common names for the first BAM/SAM")
sam_out1 = pysam.AlignmentFile(params.r1_out, "wb", template=sam1)
for i, aln in enumerate(sam1.fetch()):
    if aln.query_name in common:
        sam_out1.write(aln)
sam_out1.close()

# write the reads containing the common names for the second BAM/SAM
print("# write the reads containing the common names for the second BAM/SAM")
sam_out2 = pysam.AlignmentFile(params.r2_out, "wb", template=sam2)
for i, aln in enumerate(sam2.fetch()):
    if aln.query_name in common:
        sam_out2.write(aln)
sam_out2.close()

# closing the original files as well
print("# closing the original files as well")
sam1.close()
sam2.close()
