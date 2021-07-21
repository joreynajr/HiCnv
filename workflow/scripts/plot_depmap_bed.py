import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import argparse
import pytls

# ################# Parsing the command line arguments ##################
parser = argparse.ArgumentParser()
parser.add_argument('--bed', type=str)
parser.add_argument('--outfn', type=str)
parser.add_argument('--max-cn', type=int, default=11)
params = parser.parse_args()

# Tesster params
# params.bedgraph = '../../results/main/22Rv1/hicnv/run/22Rv1_SRR7760384_hicnv/CNV_Estimation/22Rv1.SRR7760384.cnv.bedGraph'
# params.outfn = 'comprehensive_cnv_bedgraph.png'
# params = parser.parse_known_args()[0]

# ################# Loading and parsing the data ##################
print('Loading and parsing the data')
hicnv = pd.read_table(params.bedgraph, header=None)
hicnv.columns = ['chr', 'start', 'end', 'counts', 'cn', 'cn_cat']
hicnv.loc[:, 'chr'] = ['chr{}'.format(x) for x in
                       hicnv['chr'].apply(pytls.chrName_to_chrNum)]

# if there is a . replace with a -1 (easily spot the missing data when negative)
if '.' in hicnv.loc[:, 'cn'].values:
    hicnv.loc[:, 'cn'] = hicnv.loc[:, 'cn'].replace('.', '-1')
    hicnv.loc[:, 'cn'] = hicnv.loc[:, 'cn'].astype(int)

# group the regionss by their chromosome
hicnv_grps = hicnv.groupby('chr')
chr_dict = [(int(x.replace('chr', '')), x) for x in hicnv_grps.groups.keys()]
chr_dict = sorted(chr_dict, key=lambda x: x[0])

# ################# Plotting all chromosomes on a single image ################
# making the figure + axes and unravelling axes for plotting
fig, axes = plt.subplots(figsize=(8, 11),
                         nrows=6,
                         ncols=4,
                         gridspec_kw={'hspace': 0.8, 'wspace': 0.5})
axes = np.ravel(axes)

# plotting each chromosome onto it's own axis
for chrom_num, chrom in chr_dict:

    print('Plotting: chr', chrom_num)

    # get the ax
    ax = axes[chrom_num - 1]

    # extract the data for the current chrom
    hicnv_grp_df = hicnv_grps.get_group(chrom)

    # extract the x and corresponding y points
    x_vals = []
    y_vals = []
    for (start, end, cn) in hicnv_grp_df[['start', 'end', 'cn']].values:

        # plot each cn segment
        ax.hline(cn, start, end, color='blue')

    # set the axis labels
    ax.set_title('{}'.format(chrom))
    if chrom_num > 19:
        ax.set_xlabel('Genomic\nCoordinate', labelpad=20)

    if chrom_num % 4 == 1:
        ax.set_ylabel('CN')

    # set the ytick labels
    ax.set_ylim(-2, params.max_cn)
    ax.yaxis.set_ticks(range(0, params.max_cn, 2))
    ax.yaxis.set_ticks(range(1, params.max_cn, 2), minor=True)

    # set a grid on the y-axis for easier visualization
    ax.grid(which='both', axis='y')

    # set the xticks at every 10th quartile
    # chrom_len = int(hicnv_grp_df['end'].max()) * 2
    # interval = int(chrom_len / 10)
    # r = range(0, chrom_len + 1, interval)
    # ax.set_xticks(r, minor=False)

axes[23].set_visible(False)
fn = os.path.join(params.outfn)
fig.savefig(fn, dpi=200)
