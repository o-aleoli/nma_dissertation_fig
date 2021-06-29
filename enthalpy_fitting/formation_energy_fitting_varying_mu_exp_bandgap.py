#!/usr/bin/env python3
"""
Vacancy formation energy estimation from infinite-size supercell extrapolation

First, vacancy formation energy is obtained from each supercell size, from
1x1x1 to 3x3x3 repetitions of conventional cell.
Secondly, the a0, a1 and a3 parameters for a scale law are fitted for each
vacancy formation energy.
Third, chemical potential for Na and O atoms are changed to describe
variation on chemical conditions.

The vacancy formation energy is obtained from a system in chemical and
electronic equilibria with a mass and electron resevoir, respectively. Chemical
equilibria for NaNbO3-X-Y compounds takes X and Y from the NaNbO3 ternary
diagram, obtaining chemical_environment[mu][1] and chemical_environment[mu][0].
Electronic transfer occurs from VBM to a electron reservoir with chemical
potential = CBM energy level.
Experimental band gap from the pristine orthorhombic NaNbO3 were used instead
of the DFT one.

Aproach is richly described on [DOI: 10.1103/PhysRevB.73.035215].

Alexandre Olivieri Kraus
olivieri.alexandre0@gmail.com
"""
from pymatgen.io.vasp.outputs import Oszicar, Eigenval, Vasprun
from pymatgen.electronic_structure.core import Spin
import numpy as np
import scipy.optimize as spop

def scaling_law (x, a0, a1, a3):
    return a0 + a1*x + a3*x**3

output_file = "fit_formation_energy_exp_bandgap_"

# Bulk reference
bulk_path = "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/1x1x1"
bulk_total_energy = Oszicar(bulk_path + '/OSZICAR').final_energy
bulk_band_gap, bulk_cbm, bulk_vbm, bulk_direct_band_gap = Eigenval(bulk_path + '/EIGENVAL').eigenvalue_band_properties
# fermi_energy = bulk_cbm - bulk_vbm
fermi_energy = 3.42

chemical_environment = [
    [-3.5748, -4.9297, 'A'],
    [-3.3838, -5.0252, 'B'],
    [-2.8167, -5.5923, 'C'],
    [-1.3107, -8.6042, 'D'],
    [-1.3107, -8.9756, 'E'],
    [-1.7146, -8.8409, 'F'],
    [-2.6306, -8.3829, 'G'],
    [-2.8944, -8.1192, 'H'],
    [-4.4891, -4.9297, 'I']
#    µ_Na,   µ_O,     env
]

for mu in range(len(chemical_environment)):
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
        kpoints_weights = Eigenval(path_Na0[directory] + '/EIGENVAL').kpoints_weights
        eigenvalues = Vasprun(path_Na0[directory] + '/vasprun.xml', parse_dos=False, parse_potcar_file=False).eigenvalues[Spin.up]
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
    data_Na0[..., 1] = data_Na0[..., 1] + band_filling_correction_vb_Na0 + chemical_environment[mu][0]
    #############################################################################################################################

    path_Na1 = [
        "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/Na1-",
        "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/Na1-/nkpts",
        "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/Na1-"
    ]

    total_energy_Na1 = [Oszicar(path + "/OSZICAR").final_energy for path in path_Na1]

    data_Na1 = np.array([
        [(246.12)**(-1/3) , total_energy_Na1[0] - bulk_total_energy*1],
        [(1968.92)**(-1/3), total_energy_Na1[1] - bulk_total_energy*8],
        [(6645.12)**(-1/3), total_energy_Na1[2] - bulk_total_energy*27]
    ])

    data_Na1[..., 1] = data_Na1[..., 1] + (-1)*fermi_energy + chemical_environment[mu][0]

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
        kpoints_weights = Eigenval(path_O0[directory] + '/EIGENVAL').kpoints_weights
        eigenvalues = Vasprun(path_O0[directory] + '/vasprun.xml', parse_dos=False, parse_potcar_file=False).eigenvalues[Spin.up]
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

    data_O0[..., 1] = data_O0[..., 1] + band_filling_correction_cb_O0 + chemical_environment[mu][1]

    #############################################################################################################################

    path_O2 = [
        "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/111/O2+",
        "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/222/O2+/nkpts",
        "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/vacancy/defects/333/O2+apical"
    ]

    total_energy_O2 = [Oszicar(path + "/OSZICAR").final_energy for path in path_O2]
    
    data_O2 = np.array([
        [(246.12)**(-1/3) , total_energy_O2[0] - bulk_total_energy*1],
        [(1968.92)**(-1/3), total_energy_O2[1] - bulk_total_energy*8],
        [(6645.12)**(-1/3), total_energy_O2[2] - bulk_total_energy*27]
    ])

    data_O2[..., 1] = data_O2[..., 1] + (+2)*fermi_energy + chemical_environment[mu][1]

    #############################################################################################################################

    initial_guess_Na0 = np.array([3., -18., -400.])
    (final_guess_Na0, pcov) = spop.curve_fit(
        scaling_law,
        data_Na0[..., 0],
        data_Na0[..., 1],
        initial_guess_Na0
    )

    initial_guess_Na1 = np.array([2., -21., -450.])
    (final_guess_Na1, pcov) = spop.curve_fit(
        scaling_law,
        data_Na1[..., 0],
        data_Na1[..., 1],
        initial_guess_Na1
    )

    initial_guess_O0 = np.array([4.0, -20., 445.])
    (final_guess_O0, pcov) = spop.curve_fit(
        scaling_law,
        data_O0[..., 0],
        data_O0[..., 1],
        initial_guess_O0
    )

    initial_guess_O2 = np.array([4.0, -20., 445.])
    (final_guess_O2, pcov) = spop.curve_fit(
        scaling_law,
        data_O2[..., 0],
        data_O2[..., 1],
        initial_guess_O2
    )

    with open(output_file + chemical_environment[mu][2] + ".txt", "w") as output:
        output.write("formation energy fitted with ΔE_f (L) = A0 + A1*L⁻¹ + A3*L⁻³\n")
        output.write("v Na 0\n")
        output.write("A0 = ΔE_f (∞) = %.4f (%.4f) eV\n" % (final_guess_Na0[0], initial_guess_Na0[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_Na0[1], initial_guess_Na0[1]))
        output.write("A3  =  %.4f (%.4f) eV Å³\n" % (final_guess_Na0[2], initial_guess_Na0[2]))
        output.write("v Na 1-\n")
        output.write("A0 = ΔE_f (∞) = %.4f (%.4f) eV\n" % (final_guess_Na1[0], initial_guess_Na1[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_Na1[1], initial_guess_Na1[1]))
        output.write("A3  =  %.4f (%.4f) eV Å³\n" % (final_guess_Na1[2], initial_guess_Na1[2]))
        output.write("v O 0\n")
        output.write("A0 = ΔE_f (∞) = %.4f (%.4f) eV\n" % (final_guess_O0[0], initial_guess_O0[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_O0[1], initial_guess_O0[1]))
        output.write("A3  =  %.4f (%.4f) eV Å³\n" % (final_guess_O0[2], initial_guess_O0[2]))
        output.write("v O 2+\n")
        output.write("A0 = ΔE_f (∞) = %.4f (%.4f) eV\n" % (final_guess_O2[0], initial_guess_O2[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_O2[1], initial_guess_O2[1]))
        output.write("A3  =  %.4f (%.4f) eV Å³\n" % (final_guess_O2[2], initial_guess_O2[2]))
