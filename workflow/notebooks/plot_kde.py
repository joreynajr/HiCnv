import pandas as pd 
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

fn = '../../results/main/22Rv1/hicnv/tech_run/22Rv1_SRR7760384_hicnv/Kernel_Smoothing/22Rv1.SRR7760384.chr1.kde2d_x.txt'
xdata = pd.read_table(fn, header=None, squeeze=True)

fn = '../../results/main/22Rv1/hicnv/tech_run/22Rv1_SRR7760384_hicnv/Kernel_Smoothing/22Rv1.SRR7760384.chr1.kde2d_y.txt'
ydata = pd.read_table(fn, header=None, squeeze=True)

fn = '../../results/main/22Rv1/hicnv/tech_run/22Rv1_SRR7760384_hicnv/Kernel_Smoothing/22Rv1.SRR7760384.chr1.kde2d_z.txt'
zdata = pd.read_table(fn, header=None, squeeze=True)

# rename the columns and rows
zdata.columns = ydata
zdata.index = xdata 
zdata = zdata.iloc[:, 0:1000]

# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
xi, yi = np.mgrid[zdata.index.min():zdata.index.max() + 1:1, zdata.columns.min():zdata.columns.max():1]
cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["midnightblue", "navy", "blue", 
                                                         "darkgreen", "green", "lime",
                                                         "greenyellow", "yellow"])

# calculate the top quant of nonzero entries
nonzeros = zdata.values.flatten()
nonzeros = nonzeros[nonzeros > 0]
top_quant = np.quantile(nonzeros, 0.95)

# Create a figure
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(21, 5))

fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), ax=ax)
 
# plot a density
ax.pcolormesh(xi, yi, zdata.values.reshape(xi.shape), shading='gouraud', vmin=0, vmax=top_quant, cmap=cmap)

fn = 'test.plot.png'
fig.savefig(fn)
