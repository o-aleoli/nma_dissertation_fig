#!/usr/bin/env python3
import pandas as pd
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

Na_bader_data = pd.read_fwf('na_ion.dat')
O_bader_data = pd.read_fwf('o_ion.dat')

geometry = [0.75*(210 - 30 - 20)/25.4, 70/25.4]

with PdfPages('charge_layer_na_o.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'figure.figsize': geometry,
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}',
        'legend.frameon': False,
        'legend.fancybox': False
        })

    fig, axs = plt.subplots(
            nrows=1,
            ncols=2,
            sharex=True,
            constrained_layout=True
    )

    axs[0].scatter(
            Na_bader_data['NNO'],
            Na_bader_data['layer'],
            c='black',
            marker='o',
            label='NNO'
    )

    axs[1].scatter(
            O_bader_data['NNO'],
            O_bader_data['layer'],
            c='black',
            marker='o',
            label='NNO'
    )

    axs[0].set(
        ylabel='Layers',
        ylim=[5.5, 0],
        yticks=[],
        xlabel=r'$\Delta$Charge (bulk $\%$)',
        xlim=[-0.2, 0.4]
    )

    axs[1].set(
        ylim=[5.5, 0],
        xlim=[-0.25, 0.75],
        xlabel=r'$\Delta$Charge (bulk $\%$)'
    )

    axs[1].tick_params(
        axis='y',
        which='major',
        length=10.0,
        labelleft=False,
        labelright=False,
        left=False,
        right=False
    )

    axs[0].xaxis.set_minor_locator(tck.AutoMinorLocator())
    axs[1].xaxis.set_minor_locator(tck.AutoMinorLocator())

    axs[0].text(0.05, 0.9, 'a', transform=axs[0].transAxes, fontsize='large')
    axs[1].text(0.05, 0.9, 'b', transform=axs[1].transAxes, fontsize='large')

    output.savefig()
