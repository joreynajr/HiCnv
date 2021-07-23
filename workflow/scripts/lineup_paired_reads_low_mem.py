import pysam
import argparse

# parsing the command line interface
parser = argparse.ArgumentParser()

parser.add_argument('--r1', required=True, type=str, help='R1 bam input file with query name sort.')
parser.add_argument('--r2', required=True, type=str, help='R2 bam input file with query name sort.')
parser.add_argument('--r1-out', required=True, type=str, help='R1 bam output file.')
parser.add_argument('--r2-out', required=True, type=str, help='R2 bam output file.')
params = parser.parse_args()

sam1 = pysam.AlignmentFile(params.r1, "rb")
sam2 = pysam.AlignmentFile(params.r2, "rb")

# open the files we'll write into
# write the reads containing the common names for the first BAM/SAM
sam_out1 = pysam.AlignmentFile(params.r1_out, "wb", template=sam1)
sam_out2 = pysam.AlignmentFile(params.r2_out, "wb", template=sam2)

# initial iteration 
iter1 = sam1.fetch(until_eof=True)
iter2 = sam2.fetch(until_eof=True)

aln1 = next(iter1)
aln2 = next(iter2)
        
# all other iterations
keep_going = True
while keep_going:  
        
    # update aln1 and aln1 but completely stop if no more values
    # in either one
    # write out the matching pairs
    if aln1.query_name == aln2.query_name:
        sam_out1.write(aln1)
        sam_out2.write(aln2)

        # update aln1 or completely stop if no more values
        try: 
            aln1 = next(iter1)
        except StopIteration:
            keep_going = False
            
        # update aln2 or completely stop if no more values
        try:
            aln2 = next(iter2)
        except StopIteration:
            keep_going = False


    # update aln1 or completely stop if no more values
    elif aln1.query_name < aln2.query_name:
        try: 
            aln1 = next(iter1)
        except StopIteration:
            keep_going = False

    # update aln2 or completely stop if no more values
    elif aln1.query_name > aln2.query_name:
        try: 
            aln2 = next(iter2)
        except StopIteration:
            keep_going = False

# closing the handles
sam1.close()
sam2.close()        
sam_out1.close()
sam_out2.close()




