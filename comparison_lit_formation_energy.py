#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pgf import PdfPages

width, height = 210 - 30 - 20, 120

A_vacancy_shigemi = [
    pd.read_fwf('vNa_shigemi_wada_first-principles_2005.dat', sep=' '),
    pd.read_fwf('vK_shigemi_wada_first-principles_2005.dat', sep=' '),
    pd.read_fwf('vLi_shigemi_wada_first-principles_2005.dat', sep=' '),
]

# A_vacancy_choi = pd.read_fwf('vNa_NT_choi_et_al-first-principles-2008.dat', sep=' ')

O_vacancy_shigemi = [
    pd.read_fwf('vO_NN_shigemi_wada_first-principles_2005.dat', sep=' '),
    pd.read_fwf('vO_KN_shigemi_wada_first-principles_2005.dat', sep=' '),
    pd.read_fwf('vO_LN_shigemi_wada_first-principles_2005.dat', sep=' ')
]

# O_vacancy_choi = pd.read_fwf('vO_NT_choi_et_al-first-principles-2008.dat', sep=' ')

Na_vacancy = pd.read_fwf('vNa_NN.dat', set=' ')
O_vacancy = pd.read_fwf('vO_NN.dat', set=' ')

with PdfPages('comparison_formation_energies.pdf') as output:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": "11",
        "figure.figsize": [width/25.4, height/25.4],
        "savefig.dpi": 1200,
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "pdflatex",
        "pgf.preamble": r"\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}\usepackage[separate-uncertainty = true]{siunitx}",
        "legend.frameon": False,
        "legend.fancybox": False,
        "legend.fontsize": "small"
    })

    fig, axs = plt.subplots(
        nrows=2,
        ncols=1,
        constrained_layout=True,
        sharex=True
    )

    axs[0].scatter(
        Na_vacancy['Point'],
        Na_vacancy['Enthalpy'],
        marker='s',
        color='black',
        label=r'v$_{\text{Na}}^0$ (o-NaNbO$_3$)'
    )

    axs[0].scatter(
        A_vacancy_shigemi[0]['Point'],
        A_vacancy_shigemi[0]['Enthalpy'],
        marker='o',
        color='green',
        label=r'v$_{\text{Na}}^0$ (c-NaNbO$_3$)${}^{*}$'
    )

    axs[0].scatter(
        A_vacancy_shigemi[1]['Point'],
        A_vacancy_shigemi[1]['Enthalpy'],
        marker='x',
        color='blue',
        label=r'v$_{\text{K}}^0$ (c-KNbO$_3$)${}^{*}$'
    )

    axs[0].scatter(
        A_vacancy_shigemi[2]['Point'],
        A_vacancy_shigemi[2]['Enthalpy'],
        marker='+',
        color='purple',
        label=r'v$_{\text{Li}}^0$ (t-LiNbO$_3$)${}^{*}$'
    )

    axs[1].scatter(
        O_vacancy['Point'],
        O_vacancy['Enthalpy'],
        marker='s',
        color='black',
        label=r'v$_{\text{O}}^0$ (o-NaNbO$_3$)'
    )

    axs[1].scatter(
        O_vacancy_shigemi[0]['Point'],
        O_vacancy_shigemi[0]['Enthalpy'],
        marker='o',
        color='green',
        label=r'v$_{\text{O}}^0$ (c-NaNbO$_3$)${}^{*}$'
    )

    axs[1].scatter(
        O_vacancy_shigemi[1]['Point'],
        O_vacancy_shigemi[1]['Enthalpy'],
        marker='x',
        color='blue',
        label=r'v$_{\text{O}}^0$ (c-KNbO$_3$)${}^{*}$'
    )

    axs[1].scatter(
        O_vacancy_shigemi[2]['Point'],
        O_vacancy_shigemi[2]['Enthalpy'],
        marker='+',
        color='purple',
        label=r'v$_{\text{O}}^0$ (t-LiNbO$_3$)${}^{*}$'
    )

    axs[0].set(ylabel=r'$\Delta E_f$ (\si{\electronvolt})')
    axs[1].set(ylabel=r'$\Delta E_f$ (\si{\electronvolt})')
    axs[1].set(xlabel='Chemical equilibrium',xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9],xticklabels=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])
    axs[0].minorticks_on()
    axs[0].tick_params(axis='x', which='minor', bottom=False)
    axs[1].minorticks_on()
    axs[1].tick_params(axis='x', which='minor', bottom=False)
    axs[0].legend(ncol=2, fontsize="x-small")
    axs[1].legend(ncol=1, fontsize="x-small")
    axs[0].text(0.025, 0.925, "a", transform=axs[0].transAxes)
    axs[1].text(0.025, 0.925, "b", transform=axs[1].transAxes)

    output.savefig()