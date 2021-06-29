#!/usr/bin/env python3
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.optimize as spop
from matplotlib.backends.backend_pgf import PdfPages

def enthalpy (x, a0, a1, a2):
    return a0 + a1*x + a2*x**3

width = 210.0 - 20.0 - 30.0
height = 90.0

# input_files = [
#     "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/Na0/OUTCAR",
#     "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/Na0/OUTCAR",
#     "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/Na0/OUTCAR"
# ]

# bulk_file = "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/1x1x1/OUTCAR"

# acc = []

# with open(bulk_file) as file:
#     for line in file:
#         match_total_energy_bulk = re.search(r"entropy\=\s+(-\d+.\d+)", line)
#         if match_total_energy_bulk:
#             acc.append(float(match_total_energy_bulk.group(1)))

# total_energy_bulk = acc[-1]

output_file = "fit_formation_energy"

# data_Na0 = []

# for file in input_files:
#     volume = []
#     total_energy_defect = []

#     with open(file) as outcar:
#         for line in outcar:
#             match_volume = re.search(r"cell\s\:\s+(\d+.\d+)", line)
#             match_total_energy_defect = re.search(r"entropy\=\s+(-\d+.\d+)", line)
#             if match_volume:
#                 volume.append(float(match_volume.group(1)))
#             if match_total_energy_defect:
#                 total_energy_defect.append(float(match_total_energy_defect.group(1)))

#     data_Na0.append([volume[-1]**(-1/3), total_energy_defect[-1]])
total_energy_bulk = -153.31568141

data_Na0 = np.array([
    [(246.12)**(-1/3) , -147.71953909],
    [(1968.92)**(-1/3), -1220.89734937],
    [(6645.12)**(-1/3), -4133.51574279]
])
data_Na0[..., 1] = data_Na0[..., 1] - [total_energy_bulk, total_energy_bulk*8, total_energy_bulk*27]

data_Na1 = np.array([
   [(246.12)**(-1/3) , -146.83710152],
   [(1968.92)**(-1/3), -1219.65005724],
   [(6645.12)**(-1/3), -4132.23294839]
])
data_Na1[..., 1] = data_Na1[..., 1] - np.array([total_energy_bulk, total_energy_bulk*8, total_energy_bulk*27])
data_Na1[..., 1] = data_Na1[..., 1] + (-1)*np.array([
                                               1.312 + (1.2808 - 1.0188),
                                               1.312 + (1.6886 - 1.3196) + 0.394715,
                                               1.312 + (1.6449 - 1.2959) + 0.413149
                                           ])

data_O0 = np.array([
    [(246.12)**(-1/3) ,  -143.04635301],
    [(1968.92)**(-1/3), -1217.784801],
    [(6645.12)**(-1/3), -4130.66535093]
])
data_O0[..., 1] = data_O0[..., 1] - [total_energy_bulk, total_energy_bulk*8, total_energy_bulk*27]

# data_O2 = np.array([
#     [(246.12)**(-1/3) , -150.456820190],
#     [(1968.92)**(-1/3), -1224.45990861],
#     [(6645.12)**(-1/3), -4137.24936489]
# ])
# data_O2[..., 1] = data_O2[..., 1] - np.array([total_energy_bulk, total_energy_bulk*8, total_energy_bulk*27])
# data_O2[..., 1] = data_O2[..., 1] + (+2)*np.array([
#                                                 1.312 + (1.5831 - 1.3780) + 0.153511,
#                                                 1.312 + (1.6497 - 1.3465) + ,
#                                                 1.312 + (1.6886 - 1.3196) + 0.408538
# ])

data_O2 = np.array([
    [(246.12)**(-1/3) , -150.456820190],
    [(1968.92)**(-1/3), -1224.45990861],
    [(6645.12)**(-1/3), -4137.24936489]
])
data_O2[..., 1] = data_O2[..., 1] - np.array([total_energy_bulk, total_energy_bulk*8, total_energy_bulk*27])
data_O2[..., 1] = data_O2[..., 1] + (+2)*np.array([
                                                1.312 + (1.5831 - 1.3780),
                                                1.312 + (1.6497 - 1.3465),
                                                1.312 + (1.6886 - 1.3196)
])

initial_guess_Na0 = np.array([4.0, -20., 445.])

(final_guess_Na0, pcov) = spop.curve_fit(
    enthalpy,
    data_Na0[..., 0],
    data_Na0[..., 1],
    initial_guess_Na0
)

initial_guess_Na1 = np.array([4.0, -20., 445.])

(final_guess_Na1, pcov) = spop.curve_fit(
   enthalpy,
   data_Na1[..., 0],
   data_Na1[..., 1],
   initial_guess_Na1
)

initial_guess_O0 = np.array([4.0, -20., 445.])

(final_guess_O0, pcov) = spop.curve_fit(
    enthalpy,
    data_O0[..., 0],
    data_O0[..., 1],
    initial_guess_O0
)

initial_guess_O2 = np.array([4.0, -20., 445.])

(final_guess_O2, pcov) = spop.curve_fit(
    enthalpy,
    data_O2[..., 0],
    data_O2[..., 1],
    initial_guess_O2
)

with open(output_file+".txt", "w") as output:
    output.write("Formation enthalpy fitted with ΔHf (L) = A0 + A1*L⁻¹ + A2*L⁻³\n")
    output.write("v Na 0\n")
    output.write("ΔHf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_Na0[0], initial_guess_Na0[0]))
    output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_Na0[1], initial_guess_Na0[1]))
    output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_Na0[2], initial_guess_Na0[2]))
    output.write("v Na 1-\n")
    output.write("ΔHf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_Na1[0], initial_guess_Na1[0]))
    output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_Na1[1], initial_guess_Na1[1]))
    output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_Na1[2], initial_guess_Na1[2]))
    output.write("v O 0\n")
    output.write("ΔHf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_O0[0], initial_guess_O0[0]))
    output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_O0[1], initial_guess_O0[1]))
    output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_O0[2], initial_guess_O0[2]))
    output.write("v O 2+\n")
    output.write("ΔHf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_O2[0], initial_guess_O2[0]))
    output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_O2[1], initial_guess_O2[1]))
    output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_O2[2], initial_guess_O2[2]))


with PdfPages(output_file+".pdf") as output:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size":"11",
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
        ncols=2,
        constrained_layout=True
    )
    axs[0, 0].plot(
        np.linspace(0, data_Na0[0, 0], num=100),
        enthalpy(np.linspace(0, data_Na0[0, 0], num=100), *final_guess_Na0),
        color='black',
        linestyle='dashed',
        linewidth=1.0
    )
    axs[0, 0].plot(
        data_Na0[...,0],
        data_Na0[...,1],
        "ko",
    )
    axs[0, 1].plot(
        np.linspace(0, data_O0[0, 0], num=100),
        enthalpy(np.linspace(0, data_O0[0, 0], num=100), *final_guess_O0),
        color='black',
        linestyle='dashed',
        linewidth=1.0
    )
    axs[0, 1].plot(
        data_O0[...,0],
        data_O0[...,1],
        "ko",
    )
    axs[1, 0].plot(
        np.linspace(0, data_Na1[0, 0], num=100),
        enthalpy(np.linspace(0, data_Na1[0, 0], num=100), *final_guess_Na1),
        color='black',
        linestyle='dashed',
        linewidth=1.0
    )
    axs[1, 0].plot(
        data_Na1[...,0],
        data_Na1[...,1],
        "ko",
    )
    axs[1, 1].plot(
        np.linspace(0, data_O2[0, 0], num=100),
        enthalpy(np.linspace(0, data_O2[0, 0], num=100), *final_guess_O2),
        color='black',
        linestyle='dashed',
        linewidth=1.0
    )
    axs[1, 1].plot(
        data_O2[...,0],
        data_O2[...,1],
        "ko",
    )
    labels = [
        Line2D(
            [0], [0],
            color='black',
            linewidth=1.0,
            linestyle='dashed',
            label=r"Fitted Eq.~1.13"
        )
    ]
    axs[1, 1].legend(
        handles=labels,
        loc='lower center',
        bbox_to_anchor=(0.75, -0.75),
        ncol=2
    )
    axs[0, 0].set(
        ylabel=r"$\Delta E_f$ (\si{\electronvolt})",
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xticklabels=[]
    )
    secondary_00 = axs[0, 0].secondary_xaxis("top")
    secondary_00.set_xlim([0, 0.175])
    secondary_00.set_xticks([0.00, 0.05, 0.08, 0.16])
    secondary_00.set_xticklabels([r"$\infty$", r"$3\times3\times3$", r"$2\times2\times2$", r"$1\times1\times1$"], rotation=90)
    axs[0, 1].set(
        xlim=[0, 0.175],
        ylim=[8.5, 11.25],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xticklabels=[]
    )
    secondary_01 = axs[0, 1].secondary_xaxis("top")
    secondary_01.set_xlim([0, 0.175])
    secondary_01.set_xticks([0.00, 0.05, 0.08, 0.16])
    secondary_01.set_xticklabels([r"$\infty$", r"$3\times3\times3$", r"$2\times2\times2$", r"$1\times1\times1$"], rotation=90)
    axs[1, 0].set(
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xlabel=r"$L^{-1}$ (\si{\per\angstrom})",
        ylabel=r"$\Delta E_f$ (\si{\electronvolt})"
    )
    axs[1, 1].set(
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xlabel=r"$L^{-1}$ (\si{\per\angstrom})",
    )

    axs[0, 0].text(0.90, 0.85, r"v$_{\text{Na}}^x$", transform=axs[0, 0].transAxes)
    axs[0, 1].text(0.90, 0.85, r"v$_{\text{O}}^x$", transform=axs[0, 1].transAxes)
    axs[1, 0].text(0.90, 0.85, r"v$_{\text{Na}}^{1-}$", transform=axs[1, 0].transAxes)
    axs[1, 1].text(0.90, 0.85, r"v$_{\text{O}}^{2+}$", transform=axs[1, 1].transAxes)

    output.savefig()
