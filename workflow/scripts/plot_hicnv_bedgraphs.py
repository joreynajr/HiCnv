import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import argparse
import pytls

# ################# Parsing the command line arguments ##################
# Current coord-type is only handling genomic data since the resolution we
# are using is fragment based and dynamic.
# I am thinking of implementing a bin based approach where you just use the
# index associated with a fragment but that is currently not implemented.
parser = argparse.ArgumentParser()
parser.add_argument('--bedgraph', type=str)
parser.add_argument('--outdir', type=str)
parser.add_argument('--coord-type', default='bin', choices=['bin', 'genomic'])
params = parser.parse_args()

# ################# Loading and parsing the data ##################
print('Loading and parsing the data')
hicnv = pd.read_table(params.bedgraph, header=None)
hicnv.columns = ['chr', 'start', 'end', 'counts', 'cn', 'cn_cat']
hicnv.loc[:, 'chr'] = ['chr{}'.format(x) for x in
                       hicnv['chr'].apply(pytls.chrName_to_chrNum)]
# hicnv.loc[:, 'start_bin'] = hicnv['start'] / params.resolution
# hicnv.loc[:, 'end_bin'] = hicnv['end'] / params.resolution

lw = 1

hicnv_grps = hicnv.groupby('chr')
chr_dict = [(int(x.replace('chr', '')), x) for x in hicnv_grps.groups.keys()]
chr_dict = sorted(chr_dict, key=lambda x: x[0])

## ################# Plotting each chromosome individually ##################
#max_cn = hicnv.cn.max() + 1
#for chrom_num, chrom in chr_dict:
#
#    # plotting the individual copy numbers
#    fig, ax = plt.subplots(figsize=(8, 6))
#
#    hicnv_grp_df = hicnv_grps.get_group(chrom)
#
#    if params.coord_type == 'bin':
#        ax.hlines(y='cn',
#                  xmin='start_bin',
#                  xmax='end_bin',
#                  data=hicnv_grp_df,
#                  color='blue')
#
#        i = 1
#        while i < hicnv_grp_df.shape[0]:
#            prev = hicnv_grp_df.iloc[i - 1]
#            curr = hicnv_grp_df.iloc[i]
#            yvals = [prev.cn, curr.cn]
#            min_y = min(yvals)
#            max_y = max(yvals)
#            ax.vlines(curr.start_bin, min_y, max_y, color='blue', linewidth=lw)
#            i += 1
#
#        ax.set_title('{}'.format(chrom))
#        ax.set_xlabel('Bin coordinate')
#        ax.set_ylabel('Copy Number')
#        ax.set_ylim(-0.2, max_cn)
#
#    else:
#        ax.hlines(y='cn',
#                  xmin='start',
#                  xmax='end',
#                  data=hicnv_grp_df,
#                  color='blue')
#
#        i = 1
#        while i < hicnv_grp_df.shape[0]:
#            prev = hicnv_grp_df.iloc[i - 1]
#            curr = hicnv_grp_df.iloc[i]
#            yvals = [prev.cn, curr.cn]
#            min_y = min(yvals)
#            max_y = max(yvals)
#            ax.vlines(curr.start, min_y, max_y, color='blue', linewidth=lw)
#            i += 1
#
#        # set new xaxis ticks
#        new_ticks = []
#        for tick in ax.get_xticks():
#            new_label = int(float(tick) / 1000000)
#            new_label = '{}mb'.format(new_label)
#            new_ticks.append(new_label)
#        ax.set_xticklabels(new_ticks)
#
#        # set titles
#        ax.set_title('{}'.format(chrom))
#        ax.set_xlabel('Genomic coordinate')
#        ax.set_ylabel('Copy Number')
#        ax.set_ylim(-0.2, max_cn)
#
#    fn = os.path.join(params.outdir, '{}_cnv_bedgraph.png'.format(chrom))
#    fig.savefig(fn, dpi=200)

# ################# Plotting all chromosomes on a single image ################
fig, axes = plt.subplots(figsize=(8, 11),
                         nrows=6,
                         ncols=4,
                         gridspec_kw={'hspace': 0.8, 'wspace': 0.5})
axes = np.ravel(axes)

for chrom_num, chrom in chr_dict:

    ax = axes[chrom_num - 1]

    # plotting the individual copy numbers
    hicnv_grp_df = hicnv_grps.get_group(chrom)
    ax.hlines(y='cn',
              xmin='start_bin',
              xmax='end_bin',
              data=hicnv_grp_df,
              color='blue')

    i = 1
    while i < hicnv_grp_df.shape[0]:
        prev = hicnv_grp_df.iloc[i - 1]
        curr = hicnv_grp_df.iloc[i]
        yvals = [prev.cn, curr.cn]
        min_y = min(yvals)
        max_y = max(yvals)
        ax.vlines(curr.start_bin, min_y, max_y, color='blue', linewidth=lw)
        i += 1

    ax.set_title('{}'.format(chrom))
    if chrom_num > 20:
        ax.set_xlabel('Bin')

    if chrom_num % 4 == 1:
        ax.set_ylabel('CN')

    ax.set_ylim(-0.2, max_cn)
    ax.yaxis.set_ticks(range(0, 7, 2))
    ax.yaxis.set_ticks(range(1, 8, 2), minor=True)
    ax.grid(which='both', axis='y')
    ax.grid(which='minor', axis='x')

    r = range(0,
              int(hicnv_grp_df['end_bin'].max()) + params.resolution,
              params.resolution * 2)
    ax.set_xticks(r, minor=True)
    break

axes[23].set_visible(False)
fn = os.path.join(params.outdir, 'comprehensive_cnv_bedgraph.png')
fig.savefig(fn, dpi=200)
