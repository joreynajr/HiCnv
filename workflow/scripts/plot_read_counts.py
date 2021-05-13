import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cov-bed', type=str)
parser.add_argument('--prefix', type=str)
params = parser.parse_args()

counts_data = pd.read_table(params.cov_bed, header=None, sep='\t')
counts_data.columns = ['chr', 'start', 'end', 'count']

mu_cov = counts_data['count'].mean()
print('The mean count is: {}'.format(mu_cov))

colors = sns.color_palette(n_colors=23)
i = 0
for chrom, chrom_df in counts_data.groupby('chr'):
    color = colors[i]
    fig, ax = plt.subplots()

    xdata = chrom_df['end']
    ydata = chrom_df['count']
    sns.scatterplot(x=xdata, y=ydata, color=color)
    ax.set_title('{}'.format(chrom))

    ax.set_xlabel('Coordinate')
    ax.set_ylabel('Coverage (within 40kb bins)')
    i += 1

    ax.hlines(mu_cov, xdata.min(), xdata.max(), color='red')

    fn = '{}.{}.png'.format(params.prefix, chrom)
    fig.savefig(fn, dpi=200)
