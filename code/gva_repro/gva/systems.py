"""Molecular systems, reference calculations and active-space Hamiltonians.

All electronic-structure steps use PySCF. Small molecules are treated in the
full STO-3G orbital space. Cisplatin is reduced to a CAS(10e,10o) Hamiltonian
built from symmetry-adapted MP2 natural orbitals (def2-SVP, def2 ECP on Pt).
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from pyscf import gto, scf, mp, ao2mo, symm, cc, ci, fci

ANG = 1.0

# Optional PySCF logging (added for the v2 runs). When PYSCF_LOG["dir"] is set,
# every Mole built here writes verbose-4 output to <dir>/<stem>__mol<k>.log.
PYSCF_LOG = {"dir": None, "stem": "run", "k": 0}


def set_pyscf_log(directory, stem):
    import os
    os.makedirs(directory, exist_ok=True)
    PYSCF_LOG.update({"dir": directory, "stem": stem, "k": 0})


def _mole(**kw):
    if PYSCF_LOG["dir"] is not None:
        import os
        PYSCF_LOG["k"] += 1
        kw["verbose"] = 4
        kw["output"] = os.path.join(PYSCF_LOG["dir"],
                                    f"{PYSCF_LOG['stem']}__mol{PYSCF_LOG['k']}.log")
    return gto.M(**kw)


# ---------------------------------------------------------------- geometries
def small_molecule_atoms(name: str, scale: float = 1.0):
    """Linear test molecules; equilibrium values follow the original draft."""
    if name == "H2":
        r = 0.7414 * scale
        return [("H", (0, 0, -r / 2)), ("H", (0, 0, r / 2))], "D2h"
    if name == "H4":
        r = 1.2300 * scale
        return [("H", (0, 0, (k - 1.5) * r)) for k in range(4)], "D2h"
    if name == "LiH":
        r = 1.5950 * scale
        return [("Li", (0, 0, 0)), ("H", (0, 0, r))], "C2v"
    if name == "BeH2":
        r = 1.3264 * scale
        return [("H", (0, 0, -r)), ("Be", (0, 0, 0)), ("H", (0, 0, r))], "D2h"
    raise ValueError(name)


def cisplatin_atoms(r_ptcl1: float = 2.32, r_ptcl2: float = 2.32,
                    r_ptn: float = 2.05, r_nh: float = 1.02,
                    ang_ptnh: float = 109.5):
    """Idealised square-planar cis-[Pt(NH3)2Cl2].

    Pt at the origin, molecular plane = xy, C2 axis along x. One N-H bond of each
    ammine lies in the molecular plane pointing away from the chlorides; the other
    two are placed symmetrically above and below the plane. The second NH3 is the
    mirror image (y -> -y) of the first, so the structure is exactly C2v when
    r_ptcl1 == r_ptcl2 and Cs (molecular plane retained) otherwise.
    """
    s = np.sqrt(0.5)
    cl1 = np.array([r_ptcl1 * s, r_ptcl1 * s, 0.0])
    cl2 = np.array([r_ptcl2 * s, -r_ptcl2 * s, 0.0])
    n1 = np.array([-r_ptn * s, r_ptn * s, 0.0])
    u = n1 / np.linalg.norm(n1)                    # Pt -> N direction
    perp = np.array([-u[1], u[0], 0.0])            # in-plane, perpendicular to u
    if perp @ np.array([-1.0, 1.0, 0.0]) < 0:      # point away from chlorides side
        perp = -perp
    ez = np.array([0.0, 0.0, 1.0])
    th = np.radians(180.0 - ang_ptnh)              # angle between N-H and u
    hs = []
    for phi in (0.0, 120.0, 240.0):
        p = np.radians(phi)
        d = np.cos(th) * u + np.sin(th) * (np.cos(p) * perp + np.sin(p) * ez)
        hs.append(n1 + r_nh * d)
    mirror = np.array([1.0, -1.0, 1.0])
    atoms = [("Pt", (0.0, 0.0, 0.0)), ("Cl", tuple(cl1)), ("Cl", tuple(cl2)),
             ("N", tuple(n1)), ("N", tuple(n1 * mirror))]
    atoms += [("H", tuple(h)) for h in hs]
    atoms += [("H", tuple(h * mirror)) for h in hs]
    return atoms


def build_cisplatin_mol(r_ptcl1=2.32, r_ptcl2=2.32, symmetry=True):
    atoms = cisplatin_atoms(r_ptcl1, r_ptcl2)
    mol = _mole(atom=atoms, basis={"Pt": "def2-svp", "Cl": "def2-svp",
                                   "N": "def2-svp", "H": "def2-svp"},
                ecp={"Pt": "def2-svp"}, charge=0, spin=0,
                symmetry=symmetry, verbose=0)
    return mol


def build_small_mol(name, scale=1.0):
    atoms, sub = small_molecule_atoms(name, scale)
    mol = _mole(atom=atoms, basis="sto-3g", symmetry=True,
                symmetry_subgroup=sub, verbose=0)
    return mol


# ---------------------------------------------------------- data container
@dataclass
class ActiveHamiltonian:
    label: str
    norb: int
    nelec: tuple
    ecore: float
    h1: np.ndarray            # (norb, norb) in reference canonical orbitals
    eri: np.ndarray           # (norb,)*4 chemist notation
    orbsym: np.ndarray        # abelian irrep ids (XOR-able)
    irrep_names: list
    e_rhf_total: float        # full-molecule RHF energy
    extra: dict = field(default_factory=dict)

    @property
    def nocc(self):
        return self.nelec[0]


# ------------------------------------------------- active-space SCF helpers
def _active_rhf_symmetric(h1, eri, nocc, orbsym, max_cycle=200, tol=1e-11):
    """Closed-shell SCF restricted to rotations within each irrep block."""
    n = h1.shape[0]
    C = np.eye(n)
    occ = np.arange(nocc)
    e_old = 0.0
    for it in range(max_cycle):
        D = 2.0 * C[:, occ] @ C[:, occ].T
        J = np.einsum("pqrs,rs->pq", eri, D)
        K = np.einsum("prqs,rs->pq", eri, D)
        F = h1 + J - 0.5 * K
        e = 0.5 * np.sum(D * (h1 + F))
        # diagonalise F inside each irrep block
        eps = np.zeros(n)
        Cn = np.zeros((n, n))
        sym_of = np.zeros(n, dtype=int)
        col = 0
        for g in np.unique(orbsym):
            idx = np.where(orbsym == g)[0]
            w, v = np.linalg.eigh(F[np.ix_(idx, idx)])
            for k in range(len(idx)):
                Cn[idx, col] = v[:, k]
                eps[col] = w[k]
                sym_of[col] = g
                col += 1
        order = np.argsort(eps, kind="stable")
        C, eps, sym_of = Cn[:, order], eps[order], sym_of[order]
        if abs(e - e_old) < tol and it > 2:
            break
        e_old = e
    return C, eps, sym_of, e


def _transform(h1, eri, C):
    h = C.T @ h1 @ C
    g = np.einsum("pqrs,pi,qj,rk,sl->ijkl", eri, C, C, C, C, optimize=True)
    return h, g


# --------------------------------------------------------- small molecules
def small_molecule_hamiltonian(name, scale=1.0):
    mol = build_small_mol(name, scale)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    C = mf.mo_coeff
    orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, C))
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), C.shape[1])
    n = C.shape[1]
    nocc = mol.nelectron // 2
    names = {i: nm for i, nm in zip(mol.irrep_id, mol.irrep_name)}
    return ActiveHamiltonian(label=f"{name}_s{scale:.2f}", norb=n,
                             nelec=(nocc, nocc), ecore=mol.energy_nuc(),
                             h1=h1, eri=eri, orbsym=orbsym,
                             irrep_names=[names[s] for s in orbsym],
                             e_rhf_total=mf.e_tot,
                             extra={"mol": mol, "mf": mf, "scale": scale})


# --------------------------------------------------------------- cisplatin
def mp2_natural_orbitals(mol, mf):
    """Symmetry-adapted MP2 natural orbitals (all electrons outside the ECP)."""
    pt = mp.MP2(mf).run()
    dm_mo = pt.make_rdm1()                     # MO basis, block diagonal by irrep
    C = mf.mo_coeff
    orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, C))
    nos, occs, syms = [], [], []
    for g in np.unique(orbsym):
        idx = np.where(orbsym == g)[0]
        w, v = np.linalg.eigh(dm_mo[np.ix_(idx, idx)])
        for k in range(len(idx)):
            nos.append(C[:, idx] @ v[:, k])
            occs.append(w[k])
            syms.append(g)
    nos = np.array(nos).T
    occs = np.array(occs)
    syms = np.array(syms)
    order = np.argsort(-occs, kind="stable")
    return nos[:, order], occs[order], syms[order], pt.e_corr


def select_active_by_occupation(occs, n_occ_act, n_vir_act):
    strong = np.where(occs > 1.0)[0]
    weak = np.where(occs <= 1.0)[0]
    act_o = strong[-n_occ_act:]                 # least occupied of the strong set
    act_v = weak[:n_vir_act]                    # most occupied of the weak set
    return np.concatenate([act_o, act_v])


def select_active_by_overlap(mol_prev, C_act_prev, mol, nos, occs, n_occ_act,
                             n_vir_act):
    """Continuity along a scan: pick NOs that best span the previous space."""
    S = gto.intor_cross("int1e_ovlp", mol_prev, mol)
    O = C_act_prev.T @ S @ nos
    w = np.sum(O ** 2, axis=0)
    strong = np.where(occs > 1.0)[0]
    weak = np.where(occs <= 1.0)[0]
    act_o = strong[np.argsort(-w[strong])[:n_occ_act]]
    act_v = weak[np.argsort(-w[weak])[:n_vir_act]]
    return np.sort(np.concatenate([act_o, act_v]))


def cas_hamiltonian_from_orbitals(mol, mf, nos, occs, syms, act, label,
                                  extra=None):
    """Frozen-core CAS Hamiltonian; core = strongly occupied NOs not in act."""
    n_act = len(act)
    core = np.array([k for k in np.where(occs > 1.0)[0] if k not in set(act)],
                    dtype=int)
    C_core = nos[:, core]
    C_act = nos[:, act]
    hcore = mf.get_hcore()
    e_nuc = mol.energy_nuc()
    if len(core):
        dm_core = 2.0 * C_core @ C_core.T
        vhf = mf.get_veff(mol, dm_core)
        e_core = e_nuc + np.einsum("ij,ji", dm_core, hcore) + \
            0.5 * np.einsum("ij,ji", dm_core, vhf)
        h_eff = C_act.T @ (hcore + vhf) @ C_act
    else:
        e_core, h_eff = e_nuc, C_act.T @ hcore @ C_act
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C_act), n_act)
    n_el_act = mol.nelectron - 2 * len(core)
    nocc = n_el_act // 2
    orbsym = syms[act]
    # reorder active NOs: strongly occupied first (NO order is already by occ)
    C, eps, sym_of, e_ref_el = _active_rhf_symmetric(h_eff, eri, nocc, orbsym)
    h1, g = _transform(h_eff, eri, C)
    names = {i: nm for i, nm in zip(mol.irrep_id, mol.irrep_name)}
    ex = {"mol": mol, "mf": mf, "C_act_ao": C_act @ C, "act_no_index": act,
          "act_no_occ": occs[act], "n_core": len(core),
          "e_ref_total": e_core + e_ref_el, "active_eps": eps}
    if extra:
        ex.update(extra)
    return ActiveHamiltonian(label=label, norb=n_act, nelec=(nocc, nocc),
                             ecore=e_core, h1=h1, eri=g, orbsym=sym_of,
                             irrep_names=[names[s] for s in sym_of],
                             e_rhf_total=mf.e_tot, extra=ex)


def run_rhf(mol, dm0=None):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10
    mf.max_cycle = 300
    mf.kernel(dm0)
    if not mf.converged:
        mf = scf.newton(mf)
        mf.kernel(mf.make_rdm1())
    # stability check: follow internal instabilities once if present
    mo_new = mf.stability()[0]
    if not np.allclose(mo_new, mf.mo_coeff, atol=1e-6):
        dm = mf.make_rdm1(mo_new, mf.mo_occ)
        mf.kernel(dm)
    return mf


def cisplatin_cas_hamiltonian(n_occ_act=5, n_vir_act=5, r1=2.32, r2=2.32,
                              prev=None, label=None, symmetry=True):
    mol = build_cisplatin_mol(r1, r2, symmetry=symmetry)
    dm0 = None
    if prev is not None:
        dm0 = prev["mf"].make_rdm1() if prev["mol"].natm == mol.natm else None
    mf = run_rhf(mol, dm0)
    nos, occs, syms, e_mp2 = mp2_natural_orbitals(mol, mf)
    if prev is None:
        act = select_active_by_occupation(occs, n_occ_act, n_vir_act)
    else:
        act = select_active_by_overlap(prev["mol"], prev["C_act_ao"], mol, nos,
                                       occs, n_occ_act, n_vir_act)
    label = label or f"cisplatin_CAS({n_occ_act*2}e,{n_occ_act+n_vir_act}o)"
    return cas_hamiltonian_from_orbitals(
        mol, mf, nos, occs, syms, act, label,
        extra={"e_mp2_corr": e_mp2, "r_ptcl1": r1, "r_ptcl2": r2,
               "point_group": mol.groupname})


# --------------------------------------------------- classical references
def _fake_scf(H: ActiveHamiltonian):
    n = H.norb
    mol = _mole(verbose=0)
    mol.nelectron = sum(H.nelec)
    mol.incore_anyway = True
    mol.nao_nr = lambda *a, **k: n
    mf = scf.RHF(mol)
    mf.get_hcore = lambda *a: H.h1
    mf.get_ovlp = lambda *a: np.eye(n)
    mf._eri = ao2mo.restore(8, H.eri, n)
    mf.energy_nuc = lambda *a: H.ecore
    mf.mo_coeff = np.eye(n)
    occ = np.zeros(n)
    occ[:H.nocc] = 2
    mf.mo_occ = occ
    dm = mf.make_rdm1()
    f = mf.get_fock(dm=dm)
    mf.mo_energy = np.diag(f).copy()
    mf.e_tot = mf.energy_tot(dm=dm)
    mf.converged = True
    return mf


def classical_references(H: ActiveHamiltonian):
    """HF (reference determinant), CISD, CCSD, CCSD(T), CASCI=FCI(active)."""
    mf = _fake_scf(H)
    out = {"E_ref": mf.e_tot}
    myci = ci.CISD(mf)
    myci.conv_tol = 1e-10
    myci.kernel()
    out["E_CISD"] = myci.e_tot
    mycc = cc.CCSD(mf)
    mycc.conv_tol = 1e-11
    mycc.max_cycle = 500
    mycc.kernel()
    out["E_CCSD"] = mycc.e_tot
    out["CCSD_converged"] = bool(mycc.converged)
    out["E_CCSD(T)"] = mycc.e_tot + mycc.ccsd_t()
    e_fci, civec = fci.direct_spin1.kernel(H.h1, H.eri, H.norb, H.nelec,
                                           ecore=H.ecore, conv_tol=1e-12,
                                           max_cycle=500, nroots=1)
    out["E_FCI"] = e_fci
    out["civec"] = civec
    out["t1_diag"] = float(np.sqrt(np.sum(mycc.t1 ** 2) / sum(H.nelec)))
    c0 = abs(civec[0, 0])
    out["c0_ref_weight"] = float(c0 ** 2)
    return out
