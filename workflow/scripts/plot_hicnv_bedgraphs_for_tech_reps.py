from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import argparse
import seaborn as sns
import pytls

# ################# Parsing the command line arguments ##################
parser = argparse.ArgumentParser()
parser.add_argument('--bedgraphs', type=str, nargs="+")
parser.add_argument('--outfn', type=str)
parser.add_argument('--max-cn', default=6)
params = parser.parse_args()

# Test params
# params = parser.parse_known_args()[0]
# fn = '../../results/main/22Rv1/hicnv/tech_run_auto_bandwidth/22Rv1_*_hicnv/'
# fn += figures/22Rv1.*.bedGraph'
# params.bedgraphs = glob.glob()
# params.outfn = 'comprehensive_cnv_bedgraph.png'

# ################# Plotting all chromosomes on a single image ################
# making the figure + axes and unravelling axes for plotting
fig, axes = plt.subplots(figsize=(8, 11),
                         nrows=6,
                         ncols=4,
                         gridspec_kw={'hspace': 0.8, 'wspace': 0.5})
axes = np.ravel(axes)

colors = sns.color_palette(n_colors=len(params.bedgraphs), as_cmap=True)

for i, tech_rep_fn in enumerate(params.bedgraphs):

    # loading and parsing the data
    print('Loading and parsing the data for {}'.format(tech_rep_fn))
    hicnv = pd.read_table(tech_rep_fn, header=None)
    hicnv.columns = ['chr', 'start', 'end', 'counts', 'cn', 'cn_cat']
    hicnv.loc[:, 'chr'] = ['chr{}'.format(x) for x in
                           hicnv['chr'].apply(pytls.chrName_to_chrNum)]

    # group the regions by their chromosome
    hicnv_grps = hicnv.groupby('chr')
    chr_dict = [(int(x.replace('chr', '')), x) for x in hicnv_grps.groups.keys()]
    chr_dict = sorted(chr_dict, key=lambda x: x[0])

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

            # add the start
            x_vals.append(start)
            y_vals.append(cn)

            # add the end
            x_vals.append(end)
            y_vals.append(cn)

        # plot a line plot
        ax.plot(x_vals, y_vals, color=colors[i])

        # set the axis labels
        ax.set_title('{}'.format(chrom))
        if chrom_num > 19:
            ax.set_xlabel('Genomic Coordinate')

        if chrom_num % 4 == 1:
            ax.set_ylabel('CN')

        # set the ytick labels
        ax.set_ylim(-0.2, params.max_cn)
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
fig.savefig(params.outfn, dpi=200)
