#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pgf import PdfPages

output_file = 'potentials'
width = 210 - 30 - 20
height = 297 - 30 - 20 - 30

defect_center_Na = [0.475900773386, 0.250000000000, 0.504444322758]
defect_center_O =  [0.013688917455, 0.250000000000, 0.431125063290]

for supercell in ['111', '222', '333']:
    bulk_potential_x = pd.read_csv(supercell + '_bulk_pot_x.dat', sep=' ', header=None)
    bulk_potential_y = pd.read_csv(supercell + '_bulk_pot_y.dat', sep=' ', header=None)
    bulk_potential_z = pd.read_csv(supercell + '_bulk_pot_z.dat', sep=' ', header=None)

    with PdfPages(supercell + '_' + output_file + '.pdf') as output:
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
            nrows=4,
            ncols=3,
            constrained_layout=True,
            sharey=True
        )
        
        i = 0

        for defect in ['Na0', 'Na1', 'O0', 'O2']:
            defect_potential_x = pd.read_csv(supercell + '_' + defect + '_pot_x.dat', sep=' ', header=None)
            defect_potential_y = pd.read_csv(supercell + '_' + defect + '_pot_y.dat', sep=' ', header=None)
            defect_potential_z = pd.read_csv(supercell + '_' + defect + '_pot_z.dat', sep=' ', header=None)
            v_db_x = pd.read_csv(supercell + '_' + defect + '_diff_pot_x.dat', sep=' ', header=None)
            v_db_y = pd.read_csv(supercell + '_' + defect + '_diff_pot_y.dat', sep=' ', header=None)
            v_db_z = pd.read_csv(supercell + '_' + defect + '_diff_pot_z.dat', sep=' ', header=None)
    
            axs[i, 0].plot(
                bulk_potential_x[0],
                bulk_potential_x[1],
                color='black',
                linestyle='solid'
            )
    
            axs[i, 0].plot(
                defect_potential_x[0],
                defect_potential_x[1],
                color='red',
                linestyle='solid',
                linewidth=1.0
            )
    
            axs[i, 0].plot(
                v_db_x[0],
                v_db_x[1],
                color='blue',
                linestyle='solid'
            )
    
            axs[i, 1].plot(
                bulk_potential_y[0],
                bulk_potential_y[1],
                color='black',
                linestyle='solid'
            )
    
            axs[i, 1].plot(
                defect_potential_y[0],
                defect_potential_y[1],
                color='red',
                linestyle='solid',
                linewidth=1.0
            )
    
            axs[i, 1].plot(
                v_db_y[0],
                v_db_y[1],
                color='blue',
                linestyle='solid'
            )
    
            axs[i, 2].plot(
                bulk_potential_z[0],
                bulk_potential_z[1],
                color='black',
                linestyle='solid'
            )
    
            axs[i, 2].plot(
                defect_potential_z[0],
                defect_potential_z[1],
                color='red',
                linestyle='solid',
                linewidth=1.0
            )
    
            axs[i, 2].plot(
                v_db_z[0],
                v_db_z[1],
                color='blue',
                linestyle='solid'
            )

            i += 1

        for j in range(2):
            for i in range(3):
                axs[j, i].vlines(
                    defect_center_Na[i],
                    -25,
                    5,
                    color='black',
                    linestyle='dashed',
                    linewidth=1.0
                )

                axs[j + 2, i].vlines(
                    defect_center_O[i],
                    -25,
                    5,
                    color='black',
                    linestyle='dashed',
                    linewidth=1.0
                )
        
        for i in range(3):
            for j in range(1, 3):
                axs[i, j].set(
                    xlim=[0, np.ceil(np.max(v_db_x[0]))],
                    xticklabels=[],
                    yticklabels=[]
                )

        for i in range(4):
            axs[i, 0].set(
                xlim=[0, np.ceil(np.max(v_db_x[0]))],
                ylim=[-25, 5],
                ylabel=r"Potential (\si{\electronvolt})",
                yticks=range(-25, 10, 5),
                yticklabels=['-25', '-20', '-15', '-10', '-5', '0', '5']
            )

        for i in range(3):
            axs[i, 0].set(xticklabels=[])

        axs[3, 0].set(
            ylabel=r"Potential (\si{\electronvolt})",
            xlabel=r'Position in $[100]$ (bulk \%)',
            xlim=[0, np.ceil(np.max(v_db_x[0]))]
        )
        axs[3, 1].set(
            xlabel=r'Position in $[010]$ (bulk \%)',
            xlim=[0, np.ceil(np.max(v_db_y[0]))]
        )
        axs[3, 2].set(
            xlabel=r'Position in $[001]$ (bulk \%)',
            xlim=[0, np.ceil(np.max(v_db_z[0]))]
        )

        axs[3, 0].legend(
            handles=[Line2D(
                [0], [0],
                color='black',
                linewidth=1.0,
                linestyle='solid',
                label=r'Pristine'
            )],
            loc='lower center',
            bbox_to_anchor=(0.50, -0.5)
        )

        axs[3, 1].legend(
            handles=[Line2D(
                [0], [0],
                color='red',
                linewidth=1.0,
                linestyle='solid',
                label=r'With vacancy'
            )],
            loc='lower center',
            bbox_to_anchor=(0.50, -0.5)
        )

        axs[3, 2].legend(
            handles=[Line2D(
                [0], [0],
                color='blue',
                linewidth=1.0,
                linestyle='solid',
                label=r'Difference'
            )],
            loc='lower center',
            bbox_to_anchor=(0.50, -0.5)
        )

        axs[0, 0].set_title(r'$[100]$')
        axs[0, 1].set_title(r'$[010]$')
        axs[0, 2].set_title(r'$[001]$')

        for i in range(3):
            axs[0, i].text(0.80, 0.05, r'v$_{\text{Na}}^0$', transform=axs[0, i].transAxes)
            axs[1, i].text(0.80, 0.05, r'v$_{\text{Na}}^{1-}$', transform=axs[1, i].transAxes)
            axs[2, i].text(0.80, 0.05, r'v$_{\text{O}}^0$', transform=axs[2, i].transAxes)
            axs[3, i].text(0.80, 0.05, r'v$_{\text{O}}^{2+}$', transform=axs[3, i].transAxes)

        output.savefig()
        