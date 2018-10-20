## Python helper file to plot and store graphs to represent neuronal activity
# Graphs produced:
# - Per population Inter-spike interval histogram
# - ISI coefficient of variation across populations (http://neuronaldynamics.epfl.ch/online/Ch7.S3.html)

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec

experiment = "test_0.8"
num_processes = 4
sim_time = 0.5

ISI_filename = experiment + "__" + str(num_processes) + "_ISI"
spikes_filename = experiment + "__" + str(num_processes) + "_spikes"
isicv_filename = experiment + "__" + str(num_processes) + "_ISICV"

# general plot settings
image_format = 'pdf'
plt.rcParams['figure.facecolor'] = 'white'

# from VAbenchmark_graphs.py in Brette model
def plot_hist(panel, hist, bins, width, xlabel=None, ylabel=None,
              label=None, xticks=None, xticklabels=None, xmin=None, ymax=None):
    if xlabel: panel.set_xlabel(xlabel)
    if ylabel: panel.set_ylabel(ylabel)
    for t, n in zip(bins[:-1], hist):
        panel.bar(t, n, width=width, color=None)
    if xmin: panel.set_xlim(xmin=xmin)
    if ymax: panel.set_ylim(ymax=ymax)
    if xticks is not None: panel.set_xticks(xticks)
    if xticklabels: panel.set_xticklabels(xticklabels)
    panel.text(0.8, 0.8, label, transform=panel.transAxes)


def plot_isi_hist(panel, data, label, hide_axis_labels=False):
    print("plotting ISI histogram (%s)" % label)
    bin_width = 0.2
    bins_log = np.arange(0, 8, 0.2)
    bins = np.exp(bins_log)
    isihist, bins = np.histogram(data, bins)
    xlabel = "Inter-spike interval (ms)"
    ylabel = "n in bin"
    if hide_axis_labels:
        xlabel = None
        ylabel = None
    plot_hist(panel, isihist, bins_log, bin_width, label=label,
        xlabel=xlabel, xticks=np.log([10, 100, 1000]),
        xticklabels=['10', '100', '1000'], xmin=np.log(2),
        ylabel=ylabel)
    


def plot_cvisi_hist(panel, data, label, hide_axis_labels=False):
    print("plotting CV(ISI) histogram (%s)" % label)

    bin_width = 0.1
    bins = np.arange(0, 2, bin_width)
    cvhist, bins = np.histogram(data, bins)
    xlabel = "CV(ISI)"
    ylabel = "n in bin"
    if hide_axis_labels:
        xlabel = None
        ylabel = None
    plot_hist(panel, cvhist, bins, bin_width, label=label,
              xlabel=xlabel, xticks=np.arange(0, 2, 0.5),
              ylabel=ylabel)
    plt.savefig("ISICV_" + label + "." + image_format,format=image_format,dpi=1000)



# Inter-spike-interval histograms
# Histograms of coefficients of variation of ISI
fig_settings = {  # pass these in a configuration file?
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'small',
        'legend.fontsize': 'small',
        'font.size': 14,
        'savefig.dpi': 200,
}
plt.rcParams.update(fig_settings)
CM = 1 / 2.54
plt.figure(1, figsize=(15 * CM * 2, 20 * CM))
gs = gridspec.GridSpec(2, 2 * 1, hspace=0.25, wspace=0.25)

data = np.loadtxt(ISI_filename + "_Exc")
plot_isi_hist(plt.subplot(gs[0, 0]), data, 'exc', False)
data  = np.loadtxt(isicv_filename + "_Exc")
plot_cvisi_hist(plt.subplot(gs[1, 0]), data, 'exc', False)

data = np.loadtxt(ISI_filename + "_Inh")
plot_isi_hist(plt.subplot(gs[0, 1]), data, 'inh', False)
data  = np.loadtxt(isicv_filename + "_Inh")
plot_cvisi_hist(plt.subplot(gs[1, 1]), data, 'inh', False)
plt.savefig("output." + image_format,format=image_format,dpi=1000)

    


