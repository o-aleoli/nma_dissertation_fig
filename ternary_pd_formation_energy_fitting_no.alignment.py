#!/usr/bin/env python3
from os import environ
from pymatgen.io.vasp.outputs import Oszicar, Eigenval, Vasprun
from pymatgen.electronic_structure.core import Spin
import numpy as np
import scipy.optimize as spop

def enthalpy (x, a0, a1, a2):
    return a0 + a1*x + a2*x**3

# Graph inputs
output_file = [
    "fit_formation_enthalpy_no_alignment_A",
    "fit_formation_enthalpy_no_alignment_B",
    "fit_formation_enthalpy_no_alignment_C",
    "fit_formation_enthalpy_no_alignment_D",
    "fit_formation_enthalpy_no_alignment_E",
    "fit_formation_enthalpy_no_alignment_F",
    "fit_formation_enthalpy_no_alignment_G",
    "fit_formation_enthalpy_no_alignment_H",
    "fit_formation_enthalpy_no_alignment_I"
]

# Bulk reference
bulk_path = "/home/olivieri/Documents/Academic/Mest/calculos/o-NaNbO3/bulk/1x1x1"
bulk_total_energy = Oszicar(bulk_path + '/OSZICAR').final_energy
bulk_band_gap, bulk_cbm, bulk_vbm, bulk_direct_band_gap = Eigenval(bulk_path + '/EIGENVAL').eigenvalue_band_properties

# Chemical equilibria across the ternary phase diagram
#                A          B          C          D          E          F          G          H          I
mu_Na = [-3.574830, -3.383777, -2.816660, -1.310727, -1.310727, -1.714551, -2.630598, -2.894356, -4.489096]
mu_O =  [-4.929685, -5.025211, -5.592328, -8.604195, -8.975554, -8.840946, -8.382922, -8.119165, -4.929685]

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
                                            ])

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
])

#############################################################################################################################

initial_guess_Na0 = np.array([3., -18., -400.])
initial_guess_Na1 = np.array([3., -18., -400.])
initial_guess_O0 =  np.array([3., -18., -400.])
initial_guess_O2 =  np.array([3., -18., -400.])
for environment in range(len(mu_Na)):
    data_Na0_ternary_pd = data_Na0[..., 1] + band_filling_correction_vb_Na0 + mu_Na[environment]
    (final_guess_Na0_ternary_pd, pcov) = spop.curve_fit(
        enthalpy,
        data_Na0[..., 0],
        data_Na0_ternary_pd,
        initial_guess_Na0
    )
    data_Na1_ternary_pd = data_Na1[..., 1] + mu_Na[environment]
    (final_guess_Na1_ternary_pd, pcov) = spop.curve_fit(
        enthalpy,
        data_Na1[..., 0],
        data_Na1_ternary_pd,
        initial_guess_Na1
    )
    data_O0_ternary_pd = data_O0[..., 1] + band_filling_correction_cb_O0 + mu_O[environment]
    (final_guess_O0_ternary_pd , pcov) = spop.curve_fit(
        enthalpy,
        data_O0[..., 0],
        data_O0_ternary_pd,
        initial_guess_O0
    )
    data_O2_ternary_pd = data_O2[..., 1] + mu_O[environment]
    (final_guess_O2_ternary_pd , pcov) = spop.curve_fit(
        enthalpy,
        data_O2[..., 0],
        data_O2_ternary_pd,
        initial_guess_O2
    )
    with open(output_file[environment] + ".txt", "w") as output:
        output.write("formation energy fitted with ΔEf (L) = A0 + A1*L⁻¹ + A2*L⁻³\n")
        output.write("v Na 0\n")
        output.write("ΔEf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_Na0_ternary_pd[0], initial_guess_Na0[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_Na0_ternary_pd[1], initial_guess_Na0[1]))
        output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_Na0_ternary_pd[2], initial_guess_Na0[2]))
        output.write("v Na 1-\n")
        output.write("ΔEf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_Na1_ternary_pd[0], initial_guess_Na1[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_Na1_ternary_pd[1], initial_guess_Na1[1]))
        output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_Na1_ternary_pd[2], initial_guess_Na1[2]))
        output.write("v O 0\n")
        output.write("ΔEf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_O0_ternary_pd[0], initial_guess_O0[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_O0_ternary_pd[1], initial_guess_O0[1]))
        output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_O0_ternary_pd[2], initial_guess_O0[2]))
        output.write("v O 2+\n")
        output.write("ΔEf (∞) = A0 = %.4f (%.4f) eV\n" % (final_guess_O2_ternary_pd[0], initial_guess_O2[0]))
        output.write("A1  =  %.4f (%.4f) eV Å\n" % (final_guess_O2_ternary_pd[1], initial_guess_O2[1]))
        output.write("A2  =  %.4f (%.4f) eV Å³\n" % (final_guess_O2_ternary_pd[2], initial_guess_O2[2]))

