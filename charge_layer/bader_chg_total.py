#!/usr/bin/env python3
import pandas as pd
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

b_ion_data = pd.read_fwf('b_ion.dat')
geometry = [0.5*(210 - 30 - 20)/25.4, 80/25.4]

with PdfPages('bader_chg_total.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11.0,
        'figure.figsize': geometry,
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}',
        'legend.frameon': False,
        'legend.fancybox': False
        })

    fig, axs = plt.subplots(
            nrows=1,
            ncols=1,
            constrained_layout=True
    )

    axs.scatter(b_ion_data['NNO'],
            b_ion_data['layer'],
            c='black',
            marker='o',
            label='NNO')

    axs.set(
        ylabel='Layers',
        ylim=[5.5, 0.5],
        xlabel=r'$\Delta$Charge (bulk $\%$)'
    )

    axs.xaxis.set_minor_locator(MultipleLocator(1))

    axs.tick_params(
        axis='y',
        which='major',
        length=10.0,
        labelleft=False,
        labelright=False,
        left=False,
        right=False
    )

    output.savefig()
