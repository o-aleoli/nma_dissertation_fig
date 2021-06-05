#!/usr/bin/env python3
import pandas as pd
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

pl_data = pd.read_fwf('pl.dat')
baur_index_data = pd.read_fwf('baur.dat')
average_bond_data = pd.read_fwf('bond.dat')

width = 210 - 30 - 20
height = 60

with PdfPages('nn_pl_baur_strain_evolution.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'figure.figsize': [width/25.4, height/25.4],
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}',
        'legend.frameon': False,
        'legend.fancybox': False,
        'font.size': 11
        })

    fig, axs = plt.subplots(
            nrows=1,
            ncols=3,
            sharex=True,
            constrained_layout=True
    )
# Planarity
    axs[0].scatter(
        pl_data['strain'],
        pl_data['nno'],
        c='#4CB276',
        marker='o'
    )
# Baur's DI
    axs[1].scatter(
        baur_index_data['strain'],
        baur_index_data['nno_ss'],
        c='#4CB276',
        marker='o'
    )
    axs[2].scatter(
        baur_index_data['strain'],
        baur_index_data['nno_bl'],
        c='#4CB276',
        marker='o'
    )

    for axis in axs:
        axis.set(
            xlim=(-6, 6),
            # xticks=[-5, 0, 5]
        )

    axs[0].set(ylabel=r'Surface $P_L$ (\%)')
    axs[0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0].yaxis.set_minor_locator(AutoMinorLocator())
    axs[1].set(ylabel=r'Subsurface $\Delta_d$ (\%)', xlabel=r'Biaxial strain (\%)')
    axs[1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[2].set(ylabel=r'Internal $\Delta_d$ (\%)')
    axs[2].yaxis.set_minor_locator(AutoMinorLocator())

    # plt.minorticks_off()

    axs[0].annotate('a', xy=(0.0 + 0.110, 0.875), xycoords='figure fraction')
    axs[1].annotate('b', xy=(1/3 + 0.140, 0.875), xycoords='figure fraction')
    axs[2].annotate('c', xy=(2/3 + 0.140, 0.875), xycoords='figure fraction')

    output.savefig()
