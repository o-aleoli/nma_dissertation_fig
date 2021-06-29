#!/usr/bin/env python3
"""
Vacancy formation energy estimation from infinite-size supercell extrapolation

First, vacancy formation energy is obtained from each supercell size, from 1x1x1 to 3x3x3 repetitions of conventional cell.
Secondly, the a0, a1 and a3 parameters for a scale law are fitted for each vacancy formation energy.
Third, chemical potential for Na and O atoms are changed to describe variation on chemical conditions.

The vacancy formation energy is obtained from a system in chemical and electronic equilibria with a mass and electron resevoir, respectively.
Chemical equilibria for NaNbO3-X-Y compounds takes X and Y from the NaNBO3 ternary phase diagram, obtaining mu_O and mu_Na from mu_NaNbO3 = mu_X = mu_Y.
Electronic transfer occurs from VBM to a electron reservoir with chemical potential = Fermi level.
From the dilute limit picture, electron origin is taken from pristine bulk's VBM and its destiny from the vacancy supercell's Fermi level.
The Fermi level in the defected supercell is referenced to system's VBM, preserving the point defect description without finite-sized supercell errors as they would affect VBM and Fermi levels equally.
No alignment from crystal lattice potential is done.

DFT description of electronic defect levels excessevely delocalizes them, underestimating their filling which impacts supercell's total energy.
Correction for this is taken as described in Freysold, 2014 [DOI: 10.1103/RevModPhys.86.253] review.
"""
from pymatgen.io.vasp.outputs import Oszicar, Eigenval, Vasprun
from pymatgen.electronic_structure.core import Spin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.optimize as spop
from matplotlib.backends.backend_pgf import PdfPages

def enthalpy (x, a0, a1, a2):
    return a0 + a1*x + a2*x**3

# Graph inputs
width = 210.0 - 20.0 - 30.0
height = 80.0
output_file = "fit_formation_enthalpy_no_alignment_I_new"

# Bulk reference
bulk_path = "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/1x1x1"
bulk_total_energy = Oszicar(bulk_path + '/OSZICAR').final_energy
bulk_band_gap, bulk_cbm, bulk_vbm, bulk_direct_band_gap = Eigenval(bulk_path + '/EIGENVAL').eigenvalue_band_properties

# Chemical equilibria chosen: O2--Nb2O5--NaNbO3
mu_Na = -4.4891
mu_O = -4.9297
#############################################################################################################################

path_Na0 = [
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/Na0",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/Na0/nkpts",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/Na0"
]

total_energy_Na0 = np.empty(len(path_Na0))
band_filling_correction_vb_Na0 = np.empty(len(path_Na0))
for directory in range(len(path_Na0)):
    total_energy_Na0[directory] = Oszicar(path_Na0[directory] + "/OSZICAR").final_energy
    kpoints_weights = []
    kpoints_weights = Eigenval(path_Na0[directory] + '/EIGENVAL').kpoints_weights
    eigenvalues = Vasprun(path_Na0[directory] + '/vasprun.xml', parse_dos=False, parse_potcar_file=False).eigenvalues[Spin.up]
    vbm_eigenvalues = []
    vbm_eigenvalues = [eigenvalues[kpt][eigenvalues[kpt, ..., 1] != 0.0] for kpt in range(eigenvalues.shape[0])]
    band_filling_vb = []
    for kpt in range(eigenvalues.shape[0]):
        acc = 0
        for band in range(len(vbm_eigenvalues[kpt])):
            acc += kpoints_weights[kpt]*(1 - vbm_eigenvalues[kpt][band, 1])*(vbm_eigenvalues[kpt][band, 0] - bulk_vbm)
        band_filling_vb.append(acc)
    band_filling_correction_vb_Na0[directory] = np.sum(band_filling_vb)

data_Na0 = np.array([
    [(246.12)**(-1/3) , total_energy_Na0[0] - bulk_total_energy*1],
    [(1968.92)**(-1/3), total_energy_Na0[1] - bulk_total_energy*8],
    [(6645.12)**(-1/3), total_energy_Na0[2] - bulk_total_energy*27]
])
data_Na0[..., 1] = data_Na0[..., 1] + band_filling_correction_vb_Na0 + mu_Na
#############################################################################################################################

path_Na1 = [
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/Na1-",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/Na1-/nkpts",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/Na1-"
]

total_energy_Na1 = np.empty(len(path_Na1))
band_gap_Na1 = np.empty(len(path_Na1))
cbm_Na1 = np.empty(len(path_Na1))
vbm_Na1 = np.empty(len(path_Na1))
direct_band_gap_Na1 = np.empty(len(path_Na1))
fermi_level_Na1 = np.empty(len(path_Na1))
for directory in range(len(path_Na1)):
    total_energy_Na1[directory] = Oszicar(path_Na1[directory] + '/OSZICAR').final_energy
    band_gap_Na1[directory], cbm_Na1[directory], vbm_Na1[directory], direct_band_gap_Na1[directory] = Eigenval(path_Na1[directory] + '/EIGENVAL').eigenvalue_band_properties
    fermi_level_Na1[directory] = Vasprun(path_Na1[directory] + '/vasprun.xml', parse_potcar_file=False).efermi

data_Na1 = np.array([
    [(246.12)**(-1/3) , total_energy_Na1[0] - bulk_total_energy*1],
    [(1968.92)**(-1/3), total_energy_Na1[1] - bulk_total_energy*8],
    [(6645.12)**(-1/3), total_energy_Na1[2] - bulk_total_energy*27]
])

data_Na1[..., 1] = data_Na1[..., 1] + (-1)*np.array([
                                                bulk_vbm + (fermi_level_Na1[0] - vbm_Na1[0]),
                                                bulk_vbm + (fermi_level_Na1[1] - vbm_Na1[1]),
                                                bulk_vbm + (fermi_level_Na1[2] - vbm_Na1[2])
]) + mu_Na

#############################################################################################################################

path_O0 = [
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/O0",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/O0/nkpts",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/O0apical"
]

total_energy_O0 = np.empty(len(path_O0))
band_filling_correction_cb_O0 = np.empty(len(path_O0))
for directory in range(len(path_O0)):
    total_energy_O0[directory] = Oszicar(path_O0[directory] + "/OSZICAR").final_energy
    kpoints_weights = []
    kpoints_weights = Eigenval(path_O0[directory] + '/EIGENVAL').kpoints_weights
    eigenvalues = Vasprun(path_O0[directory] + '/vasprun.xml', parse_dos=False, parse_potcar_file=False).eigenvalues[Spin.up]
    cbm_eigenvalues = []
    cbm_eigenvalues = [eigenvalues[kpt][eigenvalues[kpt, ..., 1] < 1.0] for kpt in range(eigenvalues.shape[0])]
    band_filling_cb = []
    for kpt in range(eigenvalues.shape[0]):
        acc = 0
        for band in range(len(cbm_eigenvalues[kpt])):
            acc -= kpoints_weights[kpt]*cbm_eigenvalues[kpt][band, 1]*(cbm_eigenvalues[kpt][band, 0] - bulk_cbm)
        band_filling_cb.append(acc)
    band_filling_correction_cb_O0[directory] = np.sum(band_filling_cb)    

data_O0 = np.array([
    [(246.12)**(-1/3) , total_energy_O0[0] - bulk_total_energy*1],
    [(1968.92)**(-1/3), total_energy_O0[1] - bulk_total_energy*8],
    [(6645.12)**(-1/3), total_energy_O0[2] - bulk_total_energy*27]
])

data_O0[..., 1] = data_O0[..., 1] + band_filling_correction_cb_O0 + mu_O

#############################################################################################################################

path_O2 = [
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/O2+",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/O2+/nkpts",
    "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/O2+apical"
]

total_energy_O2 = np.empty(len(path_O2))
band_gap_O2 = np.empty(len(path_O2))
cbm_O2 = np.empty(len(path_O2))
vbm_O2 = np.empty(len(path_O2))
direct_band_gap_O2 = np.empty(len(path_O2))
fermi_level_O2 = np.empty(len(path_O2))
for directory in range(len(path_O2)):
    total_energy_O2[directory] = Oszicar(path_O2[directory] + '/OSZICAR').final_energy
    band_gap_O2[directory], cbm_O2[directory], vbm_O2[directory], direct_band_gap_O2[directory] = Eigenval(path_O2[directory] + '/EIGENVAL').eigenvalue_band_properties
    fermi_level_O2[directory] = Vasprun(path_O2[directory] + '/vasprun.xml', parse_potcar_file=False).efermi
# Gambiarra feia
fermi_level_O2[2] = Vasprun(path_O2[2] + '/locpot/vasprun.xml', parse_potcar_file=False).efermi

data_O2 = np.array([
    [(246.12)**(-1/3) , total_energy_O2[0] - bulk_total_energy*1],
    [(1968.92)**(-1/3), total_energy_O2[1] - bulk_total_energy*8],
    [(6645.12)**(-1/3), total_energy_O2[2] - bulk_total_energy*27]
])

data_O2[..., 1] = data_O2[..., 1] + (+2)*np.array([
                                                bulk_vbm + (fermi_level_O2[0] - vbm_O2[0]),
                                                bulk_vbm + (fermi_level_O2[1] - vbm_O2[1]),
                                                bulk_vbm + (fermi_level_O2[2] - vbm_O2[2])
]) + mu_O

#############################################################################################################################

initial_guess_Na0 = np.array([3., -18., -400.])
(final_guess_Na0, pcov) = spop.curve_fit(
    enthalpy,
    data_Na0[..., 0],
    data_Na0[..., 1],
    initial_guess_Na0
)

initial_guess_Na1 = np.array([2., -21., -450.])
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
    output.write("formation energy fitted with ΔHf (L) = A0 + A1*L⁻¹ + A2*L⁻³\n")
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
            label=r"Fitted Eq.~1.14"
        )
    ]
    axs[1, 1].legend(
        handles=labels,
        loc='lower center',
        bbox_to_anchor=(0.75, -0.8),
        ncol=1
    )
    axs[0, 0].set(
        ylabel=r"$\Delta E_f$ @ O$_2$ (\si{\electronvolt})",
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xticklabels=[]
    )
    secondary_00 = axs[0, 0].secondary_xaxis("top")
    secondary_00.set_xlim([0, 0.175])
    secondary_00.set_xticks([0.00, 0.05, 0.08, 0.16])
    secondary_00.set_xticklabels([r"$\infty$", r"$3\!\times\!3\!\times\!3$", r"$2\!\times\!2\!\times\!2$", r"$1\!\times\!1\!\times\!1$"], fontsize=9)
    axs[0, 1].set(
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xticklabels=[]
    )
    secondary_01 = axs[0, 1].secondary_xaxis("top")
    secondary_01.set_xlim([0, 0.175])
    secondary_01.set_xticks([0.00, 0.05, 0.08, 0.16])
    secondary_01.set_xticklabels([r"$\infty$", r"$3\!\times\!3\!\times\!3$", r"$2\!\times\!2\!\times\!2$", r"$1\!\times\!1\!\times\!1$"], fontsize=9)
    axs[1, 0].set(
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xlabel=r"$L^{-1}$ (\si{\per\angstrom})",
        ylabel=r"$\Delta E_f$ @ O$_2$ (\si{\electronvolt})"
    )
    axs[1, 1].set(
        xlim=[0, 0.175],
        xticks=(0.00, 0.05, 0.08, 0.16),
        xlabel=r"$L^{-1}$ (\si{\per\angstrom})",
    )

    axs[0, 0].annotate(
        r"\SI{3.41}{\electronvolt}",
        xy=(0.0000, 3.4088),
        xycoords='data',
        xytext=(0.5, 0.75),
        arrowprops=dict(
            facecolor='black',
            width=1,
            headwidth=4,
            shrink=0.05
        ),
        horizontalalignment='left',
        verticalalignment='top',
        textcoords='axes fraction'
    )
    axs[1, 0].annotate(
        r"\SI{2.87}{\electronvolt}",
        xy=(0.0000, 2.8681),
        xycoords='data',
        xytext=(0.5, 0.75),
        arrowprops=dict(
            facecolor='black',
            width=1,
            headwidth=4,
            shrink=0.05
        ),
        horizontalalignment='left',
        verticalalignment='top',
        textcoords='axes fraction'
    )

    axs[0, 1].annotate(
        r"\SI{5.32}{\electronvolt}",
        xy=(0.0000, 5.3180),
        xycoords='data',
        xytext=(0.5, 0.75),
        arrowprops=dict(
            facecolor='black',
            width=1,
            headwidth=4,
            shrink=0.05
        ),
        horizontalalignment='left',
        verticalalignment='top',
        textcoords='axes fraction'
    )
    axs[1, 1].annotate(
        r"\SI{1.76}{\electronvolt}",
        xy=(0.0000, 1.7654),
        xycoords='data',
        xytext=(0.5, 0.75),
        arrowprops=dict(
            facecolor='black',
            width=1,
            headwidth=4,
            shrink=0.05
        ),
        horizontalalignment='left',
        verticalalignment='top',
        textcoords='axes fraction'
    )

    axs[0, 0].text(0.95, 0.85, "a", transform=axs[0, 0].transAxes)
    axs[0, 1].text(0.95, 0.85, "b", transform=axs[0, 1].transAxes)
    axs[1, 0].text(0.95, 0.85, "c", transform=axs[1, 0].transAxes)
    axs[1, 1].text(0.95, 0.85, "d", transform=axs[1, 1].transAxes)

    output.savefig()
