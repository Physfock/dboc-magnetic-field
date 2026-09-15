#!/usr/bin/env python
# DBOC in a magnetic field parallel to the molecular axis: UHF vs CISD
# UHF branch : complex ChronusQ integrals, determinant overlap
# CISD branch: real part of the same integrals -> real UHF -> UCISD,
#              overlap via pyscf.ci.ucisd.overlap()
# Laplacian from phase-corrected overlaps |S|, Eqs. (22)-(23) of the paper.

import os
import h5py
import numpy as np
import subprocess
from functools import reduce

from pyscf import gto, scf, ci
from pyscf.gto.basis import parse_gaussian
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
MOLECULE = 'H2'                          # 'H2' or 'LiH'
R_min, R_max, R_step = 0.5, 5.5, 0.1     # Bohr
dR = 0.001                               # finite-difference displacement (Bohr)
B_FIELD = -0.2                           # field along z (a.u.), passed to ChronusQ
BASIS_FILE = 'unc-cc-pvtz.gbs'
CONV_TOL = 1e-14

# atoms (atom 1 at the origin, atom 2 at (0, 0, R)) and nuclear masses (m_e)
MOLECULES = {
    'H2':  (('H', 'H'),  (1836.152673426, 1836.152673426)),
    'LiH': (('Li', 'H'), (12652.669, 1836.152673426)),
}
(A, B), (M_1, M_2) = MOLECULES[MOLECULE]
factor_1 = 1.0 / M_1
factor_2 = 1.0 / M_2
BASIS_PATH = os.path.abspath(BASIS_FILE)

R_values = []
DBOC_UHF_values, PEC_UHF_values = [], []
DBOC_CISD_values, PEC_CISD_values = [], []


def run_chronusq(mol):
    bohr_to_angstrom = 0.529177210903
    atom_string = "\n".join([
        f"{mol.atom_symbol(i)} {c[0]*bohr_to_angstrom:.10f} "
        f"{c[1]*bohr_to_angstrom:.10f} {c[2]*bohr_to_angstrom:.10f}"
        for i, c in enumerate(mol.atom_coords())
    ])

    input_content = f"""
[Molecule]
charge = {mol.charge}
mult = {mol.spin + 1}
geom:
{atom_string}

[QM]
reference = COMPLEX UHF
job = SCF

[BASIS]
basis = {BASIS_PATH}
basisType = GIAO

[SCF]
Field = Magnetic 0.0 0.0 {B_FIELD}
ENETOL = 1e-12
DENTOL = 1.e-12

[MISC]
nsmp = 1
mem = 1000 MB

[INTS]
alg = incore
    """

    input_file = 'chronusq_input.inp'
    with open(input_file, 'w') as f:
        f.write(input_content)

    try:
        subprocess.run(["chronusq", input_file],
                       capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error calling ChronusQ: {e}\n{e.stdout}\n{e.stderr}")
        raise

    with h5py.File('chronusq_input.bin', 'r') as f:
        h_core         = f['/INTS/CORE_HAMILTONIAN_SCALAR'][:]
        overlap_matrix = f['/INTS/OVERLAP'][:]
        eri_tensor     = f['/ERI'][:]
        chronus_scf    = f['/SCF/TOTAL_ENERGY'][:][0]

    snorm           = 1 / overlap_matrix[0, 0].real
    h_core          = h_core * (snorm / 2)
    overlap_matrix  = overlap_matrix * snorm
    eri_tensor      = eri_tensor * (snorm * snorm)

    return h_core, overlap_matrix, eri_tensor, chronus_scf


def make_mol(atom_str):
    return gto.M(
        atom=atom_str,
        basis={sym: parse_gaussian.load(BASIS_FILE, sym) for sym in (A, B)},
        verbose=0, unit='Bohr', charge=0, spin=0
    )


def build_uhf(mol, h_core, ovlp, eri):
    mf = mol.apply(scf.UHF)
    mf.get_hcore = lambda *args: h_core
    mf.get_ovlp  = lambda *args: ovlp
    mf._eri      = eri
    mf.conv_tol      = CONV_TOL
    mf.conv_tol_grad = CONV_TOL
    mf.kernel()
    return mf


def uhf_det_overlap(mol_ref, mol_pert, myhf_ref, myhf_pert):
    nalpha, nbeta = myhf_ref.mol.nelec
    s12 = gto.intor_cross('cint1e_ovlp_sph', mol_ref, mol_pert)

    m0a = myhf_ref.mo_coeff[0][:, :nalpha]
    m0b = myhf_ref.mo_coeff[1][:, :nbeta]
    mpa = myhf_pert.mo_coeff[0][:, :nalpha]
    mpb = myhf_pert.mo_coeff[1][:, :nbeta]

    det_a = np.linalg.det(reduce(np.dot, (m0a.T, s12, mpa)))
    det_b = np.linalg.det(reduce(np.dot, (m0b.T, s12, mpb)))
    return det_a * det_b


def cisd_overlap(mol_ref, mol_pert, myhf_ref, ci_ref, myhf_pert, ci_pert):
    s12 = gto.intor_cross('cint1e_ovlp_sph', mol_ref, mol_pert)
    mo1a, mo1b = myhf_ref.mo_coeff
    mo2a, mo2b = myhf_pert.mo_coeff

    s12_mo = (reduce(np.dot, (mo1a.T, s12, mo2a)),
              reduce(np.dot, (mo1b.T, s12, mo2b)))

    nmo  = (myhf_pert.mo_energy[0].size, myhf_pert.mo_energy[1].size)
    nocc = myhf_pert.nelec

    return ci.ucisd.overlap(ci_ref.ci, ci_pert.ci, nmo, nocc, s12_mo)


def compute_laplace_both(mol_neg, mol_pos, mol_zero,
                         myhf0_c, myhf0_r, ci0):
    h_neg, S_neg, eri_neg, _ = run_chronusq(mol_neg)
    h_pos, S_pos, eri_pos, _ = run_chronusq(mol_pos)

    myhf1_c = build_uhf(mol_neg, h_neg, S_neg, eri_neg)
    myhf2_c = build_uhf(mol_pos, h_pos, S_pos, eri_pos)

    Sm_uhf = uhf_det_overlap(mol_zero, mol_neg, myhf0_c, myhf1_c)
    Sp_uhf = uhf_det_overlap(mol_zero, mol_pos, myhf0_c, myhf2_c)
    laplace_uhf = (abs(Sp_uhf) + abs(Sm_uhf) - 2.0) / (dR ** 2)       # Eq. (23)

    h_neg_r, S_neg_r, eri_neg_r = h_neg.real, S_neg.real, eri_neg.real
    h_pos_r, S_pos_r, eri_pos_r = h_pos.real, S_pos.real, eri_pos.real

    myhf1_r = build_uhf(mol_neg, h_neg_r, S_neg_r, eri_neg_r)
    myhf2_r = build_uhf(mol_pos, h_pos_r, S_pos_r, eri_pos_r)

    ci1 = ci.UCISD(myhf1_r); ci1.conv_tol = CONV_TOL; ci1.kernel()
    ci2 = ci.UCISD(myhf2_r); ci2.conv_tol = CONV_TOL; ci2.kernel()

    Sm_cisd = cisd_overlap(mol_zero, mol_neg, myhf0_r, ci0, myhf1_r, ci1)
    Sp_cisd = cisd_overlap(mol_zero, mol_pos, myhf0_r, ci0, myhf2_r, ci2)
    laplace_cisd = (abs(Sp_cisd) + abs(Sm_cisd) - 2.0) / (dR ** 2)    # Eq. (23)

    print(f"    UHF : Sp={Sp_uhf}  Sm={Sm_uhf}  Laplace={laplace_uhf}")
    print(f"    CISD: Sp={Sp_cisd:.10f}  Sm={Sm_cisd:.10f}  Laplace={laplace_cisd:.6f}")

    return laplace_uhf, laplace_cisd


with open("DBOC_UHF.txt", "w") as f_dboc_uhf, \
     open("DBOC_CISD.txt", "w") as f_dboc_cisd, \
     open("PEC_UHF.txt", "w") as f_pec_uhf, \
     open("PEC_CISD.txt", "w") as f_pec_cisd:

    for R in np.arange(R_min, R_max + R_step, R_step):
        print(f"\n=== R = {R:.3f} Bohr ===")

        mol_zero = make_mol(f"{A} 0 0 0; {B} 0 0 {R}")

        h0, S0, eri0, chronus_scf = run_chronusq(mol_zero)

        myhf0_c = build_uhf(mol_zero, h0, S0, eri0)

        h0_r, S0_r, eri0_r = h0.real, S0.real, eri0.real
        myhf0_r = build_uhf(mol_zero, h0_r, S0_r, eri0_r)
        ci0 = ci.UCISD(myhf0_r); ci0.conv_tol = CONV_TOL; ci0.kernel()

        energy_uhf  = myhf0_c.e_tot.real
        energy_cisd = ci0.e_tot

        print(f"  ChronusQ UHF : {chronus_scf:.10f}")
        print(f"  PySCF UHF    : {energy_uhf:.10f}")
        print(f"  PySCF CISD   : {energy_cisd:.10f}")

        # -------------------- atom 1 --------------------
        lap_xx_uhf, lap_xx_cisd = compute_laplace_both(
            make_mol(f"{A} {-dR} 0 0; {B} 0 0 {R}"),
            make_mol(f"{A} {+dR} 0 0; {B} 0 0 {R}"),
            mol_zero, myhf0_c, myhf0_r, ci0)
        lap_yy_uhf, lap_yy_cisd = compute_laplace_both(
            make_mol(f"{A} 0 {-dR} 0; {B} 0 0 {R}"),
            make_mol(f"{A} 0 {+dR} 0; {B} 0 0 {R}"),
            mol_zero, myhf0_c, myhf0_r, ci0)
        lap_zz_uhf, lap_zz_cisd = compute_laplace_both(
            make_mol(f"{A} 0 0 {-dR}; {B} 0 0 {R}"),
            make_mol(f"{A} 0 0 {+dR}; {B} 0 0 {R}"),
            mol_zero, myhf0_c, myhf0_r, ci0)

        laplace_atom1_uhf  = lap_xx_uhf  + lap_yy_uhf  + lap_zz_uhf
        laplace_atom1_cisd = lap_xx_cisd + lap_yy_cisd + lap_zz_cisd

        # -------------------- atom 2 --------------------
        lap_xx_uhf, lap_xx_cisd = compute_laplace_both(
            make_mol(f"{A} 0 0 0; {B} {-dR} 0 {R}"),
            make_mol(f"{A} 0 0 0; {B} {+dR} 0 {R}"),
            mol_zero, myhf0_c, myhf0_r, ci0)
        lap_yy_uhf, lap_yy_cisd = compute_laplace_both(
            make_mol(f"{A} 0 0 0; {B} 0 {-dR} {R}"),
            make_mol(f"{A} 0 0 0; {B} 0 {+dR} {R}"),
            mol_zero, myhf0_c, myhf0_r, ci0)
        lap_zz_uhf, lap_zz_cisd = compute_laplace_both(
            make_mol(f"{A} 0 0 0; {B} 0 0 {R - dR}"),
            make_mol(f"{A} 0 0 0; {B} 0 0 {R + dR}"),
            mol_zero, myhf0_c, myhf0_r, ci0)

        laplace_atom2_uhf  = lap_xx_uhf  + lap_yy_uhf  + lap_zz_uhf
        laplace_atom2_cisd = lap_xx_cisd + lap_yy_cisd + lap_zz_cisd

        # DBOC = -1/2 sum_I <Laplacian_I> / M_I, Eq. (16)
        dboc_uhf  = -0.5 * (factor_1 * laplace_atom1_uhf  + factor_2 * laplace_atom2_uhf)
        dboc_cisd = -0.5 * (factor_1 * laplace_atom1_cisd + factor_2 * laplace_atom2_cisd)

        print(f"  DBOC(UHF)  = {dboc_uhf:.10f} a.u.")
        print(f"  DBOC(CISD) = {dboc_cisd:.10f} a.u.")

        R_values.append(R)
        DBOC_UHF_values.append(dboc_uhf)
        DBOC_CISD_values.append(dboc_cisd)
        PEC_UHF_values.append(energy_uhf)
        PEC_CISD_values.append(energy_cisd)

        f_dboc_uhf.write(f"{R:.2f}\t{dboc_uhf:.10f}\n")
        f_dboc_cisd.write(f"{R:.2f}\t{dboc_cisd:.10f}\n")
        f_pec_uhf.write(f"{R:.2f}\t{energy_uhf:.10f}\n")
        f_pec_cisd.write(f"{R:.2f}\t{energy_cisd:.10f}\n")


plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 10,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

R_fine = np.linspace(R_min, R_max, 500)
spl_uhf  = make_interp_spline(R_values, DBOC_UHF_values, k=3)
spl_cisd = make_interp_spline(R_values, DBOC_CISD_values, k=3)

fig, ax = plt.subplots(figsize=(3.4, 2.8))

ax.plot(R_fine, spl_uhf(R_fine), '-', color='tab:blue', lw=1.4, label='UHF')
ax.plot(R_values, DBOC_UHF_values, 'o', ms=3.5, mfc='white',
        mec='tab:blue', mew=1.0)

ax.plot(R_fine, spl_cisd(R_fine), '--', color='tab:red', lw=1.4, label='CISD')
ax.plot(R_values, DBOC_CISD_values, 's', ms=3.5, mfc='white',
        mec='tab:red', mew=1.0)

ax.set_xlabel(r"$R$ (Bohr)")
ax.set_ylabel(r"DBOC (a.u.)")
ax.text(0.90, 0.90, {'H2': r"H$_2$", 'LiH': r"LiH"}[MOLECULE], transform=ax.transAxes,
        fontsize=10, va='top', ha='right')
ax.legend(frameon=False, fontsize=9, loc='best')
ax.tick_params(which='both', labelsize=9)

fig.tight_layout(pad=0.4)
fig.savefig(f"DBOC_UHF_vs_CISD_{MOLECULE}_PRA.png", dpi=600)
plt.show()
