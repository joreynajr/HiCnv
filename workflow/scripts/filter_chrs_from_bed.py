import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bed-like',
                    help='A bedlike file where the first'
                         + ' column contains chromosomes.',
                    type=str)
parser.add_argument('--include-list',
                    help='File where each line is a chromosome'
                         + ' which will be included.',
                    default=None,
                    type=str)
parser.add_argument('--exclude-list',
                    help='File where each line is a chromosome'
                         + ' which will be included.',
                    default=None,
                    type=str)
parser.add_argument('-o',
                    help='A bedlike file where the first column'
                         + ' contains chromosomes.',
                    type=str)
params = parser.parse_args()

# read the include and exclude files when present
include = []
if params.include_list:
    include = open(params.include_list).readlines()
    include = [x.strip() for x in include]
exclude = []
if params.exclude_list:
    exclude = open(params.exclude_list).readlines()
    exclude = [x.strip() for x in exclude]

# print('include:', include)
# print('exclude:', exclude)

with open(params.bed_like) as fr, open(params.o, 'w') as fw:
    for line in fr:
        info = line.strip().split()
        chrom = info[0]
        if chrom in include and chrom not in exclude:
            fw.write(line)
