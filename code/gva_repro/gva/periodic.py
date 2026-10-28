"""Periodic systems: Bloch-orbital Hamiltonians, crystal-momentum screening,
time-reversal and spin-flip tying, and a complex determinant-space simulator.

Model systems are chains of H2 units (two H atoms per cell, intra-unit
distance r, inter-unit distance d, lattice constant a = r + d) in the STO-3G
basis. The cell is three-dimensionally periodic with 10 A of vacuum between
parallel chains and is sampled on a Gamma-centred N_k x 1 x 1 mesh along the
chain. The k-point Hamiltonian describes the Born-von Karman supercell of N_k
cells. Exchange is treated without a divergence correction (exxdiv=None) so
that the integrals, the KRHF energy and the supercell energy are consistent.

Orbital ordering. Spatial orbitals are ordered occupied first, then virtual;
within each block by k-point index and then band. Spin-orbital 2p+s (s = 0
alpha, 1 beta) as in the molecular code.

Gauge. Orbitals at k and -k are chosen complex conjugate, phi_{n,-k} =
phi_{n,k}^*, and orbitals at time-reversal-invariant k-points are real. The
real supercell basis used for the UCCSD reference ansatz and for the classical
reference calculations is obtained from the Bloch basis by the unitary
c = (phi_k + phi_-k)/sqrt2, s = (phi_k - phi_-k)/(i sqrt2), which leaves the
Hartree-Fock determinant unchanged.

Integral convention: eri[p,q,r,s] = (pq|rs) = int phi_p^*(1) phi_q(1)
phi_r^*(2) phi_s(2) / r12 over the supercell, h1[p,q] = <p|h|q>.
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize

from pyscf.fci import cistring, direct_nosym, spin_op

from .ansatz import Excitation, uccsd_pool, cnot_count_ladder

PBC_LOG = {"dir": None, "stem": "pbc"}


def set_pbc_log(directory, stem):
    os.makedirs(directory, exist_ok=True)
    PBC_LOG.update({"dir": directory, "stem": stem})


# ================================================================ cell, SCF
def make_cell(r, d, basis="sto-3g", vacuum=10.0, verbose=0, output=None):
    from pyscf.pbc import gto as pgto
    cell = pgto.Cell()
    a = r + d
    cell.atom = f"H 0 0 0; H {r:.10f} 0 0"
    cell.a = np.array([[a, 0, 0], [0, vacuum, 0], [0, 0, vacuum]])
    cell.dimension = 3          # chains 10 A apart; k-sampling along the chain only
    cell.unit = "Angstrom"
    cell.basis = basis
    cell.precision = 1e-12
    cell.verbose = verbose
    if output is not None:
        cell.output = output
    cell.build()
    return cell


@dataclass
class PeriodicHamiltonian:
    label: str
    r: float
    d: float
    nk: int
    norb: int
    nelec: tuple
    ecore: float
    h1: np.ndarray                  # complex, Bloch basis
    eri: np.ndarray                 # complex, Bloch basis, (pq|rs)
    h1_real: np.ndarray             # real supercell basis
    eri_real: np.ndarray
    kidx: np.ndarray                # k-point index of each spatial orbital
    band: np.ndarray
    partner: np.ndarray             # index of the time-reversed orbital
    mo_energy: np.ndarray
    e_krhf_cell: float
    U: np.ndarray                   # Bloch -> real basis
    extra: dict = field(default_factory=dict)

    @property
    def nocc(self):
        return self.nelec[0]


def _gauge_fixed_orbitals(mf, nk):
    """Solve F C = S C e at each k with C(-k) = C(k)^* and real C at TRIM k."""
    fock = mf.get_fock()
    ovlp = mf.get_ovlp()
    C, E = [None] * nk, [None] * nk
    for j in range(nk):
        jb = (-j) % nk
        if C[j] is not None:
            continue
        if jb == j:
            f, s = np.asarray(fock[j]), np.asarray(ovlp[j])
            assert np.abs(f.imag).max() < 1e-9 and np.abs(s.imag).max() < 1e-9
            e, c = sla.eigh(f.real, s.real)
            C[j], E[j] = c.astype(complex), e
        else:
            e, c = sla.eigh(fock[j], ovlp[j])
            C[j], E[j] = c, e
            C[jb], E[jb] = c.conj(), e.copy()
    return C, E


def periodic_hamiltonian(r, d, nk, basis="sto-3g", conv_tol=1e-12, log_stem=None):
    from pyscf.pbc import scf as pscf
    out = None
    if PBC_LOG["dir"] is not None:
        out = os.path.join(PBC_LOG["dir"], f"{log_stem or PBC_LOG['stem']}__krhf.log")
    cell = make_cell(r, d, basis, verbose=4 if out else 0, output=out)
    kpts = cell.make_kpts([nk, 1, 1], wrap_around=False)
    mf = pscf.KRHF(cell, kpts, exxdiv=None).density_fit()
    mf.conv_tol = conv_tol
    mf.max_cycle = 200
    mf.kernel()
    assert mf.converged, "KRHF did not converge"
    C, E = _gauge_fixed_orbitals(mf, nk)
    nmo = C[0].shape[1]
    nocc_cell = cell.nelectron // 2
    hcore = mf.get_hcore()
    # orbital list: occupied (k-major, band-minor), then virtual
    orbs = [(j, n) for j in range(nk) for n in range(nocc_cell)] + \
           [(j, n) for j in range(nk) for n in range(nocc_cell, nmo)]
    norb = len(orbs)
    kidx = np.array([o[0] for o in orbs])
    band = np.array([o[1] for o in orbs])
    pos = {o: p for p, o in enumerate(orbs)}
    partner = np.array([pos[((-j) % nk, n)] for j, n in orbs])
    mo_energy = np.array([E[j][n] for j, n in orbs])
    # check aufbau: highest occupied below lowest virtual
    ehomo = mo_energy[:nk * nocc_cell].max()
    elumo = mo_energy[nk * nocc_cell:].min()
    # one-electron part
    h1 = np.zeros((norb, norb), complex)
    for j in range(nk):
        hm = C[j].conj().T @ hcore[j] @ C[j]
        idx = [pos[(j, n)] for n in range(nmo)]
        h1[np.ix_(idx, idx)] = hm
    # two-electron part, momentum-conserving blocks only
    eri = np.zeros((norb,) * 4, complex)
    for k1 in range(nk):
        for k2 in range(nk):
            for k3 in range(nk):
                k4 = (k1 - k2 + k3) % nk
                blk = mf.with_df.ao2mo((C[k1], C[k2], C[k3], C[k4]),
                                       (kpts[k1], kpts[k2], kpts[k3], kpts[k4]),
                                       compact=False).reshape((nmo,) * 4) / nk
                i1 = [pos[(k1, n)] for n in range(nmo)]
                i2 = [pos[(k2, n)] for n in range(nmo)]
                i3 = [pos[(k3, n)] for n in range(nmo)]
                i4 = [pos[(k4, n)] for n in range(nmo)]
                eri[np.ix_(i1, i2, i3, i4)] = blk
    ecore = nk * cell.energy_nuc()
    # Bloch -> real supercell basis
    U = np.zeros((norb, norb), complex)
    done = set()
    for p in range(norb):
        if p in done:
            continue
        q = partner[p]
        if q == p:
            U[p, p] = 1.0
            done.add(p)
        else:
            lo, hi = min(p, q), max(p, q)
            U[lo, lo] = U[hi, lo] = 1 / np.sqrt(2)           # cos-like at lo
            U[lo, hi] = -1j / np.sqrt(2)                       # sin-like at hi
            U[hi, hi] = 1j / np.sqrt(2)
            done.update((p, q))
    h1r = U.conj().T @ h1 @ U
    erir = np.einsum("ap,bq,cr,ds,abcd->pqrs", U.conj(), U, U.conj(), U, eri,
                     optimize=True)
    imag_max = max(np.abs(h1r.imag).max(), np.abs(erir.imag).max())
    assert imag_max < 1e-8, f"real-basis integrals not real ({imag_max:.1e})"
    h1r, erir = h1r.real.copy(), erir.real.copy()
    nel = nk * nocc_cell
    H = PeriodicHamiltonian(
        label=f"H2chain r={r:.2f} d={d:.2f} Nk={nk}", r=r, d=d, nk=nk, norb=norb,
        nelec=(nel, nel), ecore=ecore, h1=h1, eri=eri, h1_real=h1r, eri_real=erir,
        kidx=kidx, band=band, partner=partner, mo_energy=mo_energy,
        e_krhf_cell=float(mf.e_tot), U=U,
        extra={"homo": float(ehomo), "lumo": float(elumo), "a": r + d,
               "kpts_frac": [j / nk for j in range(nk)],
               "imag_max_real_basis": float(imag_max), "n_ao_cell": nmo,
               "krhf_cycles": None})
    H.extra["E_HF_from_integrals"] = hf_energy(H.h1, H.eri, H.nocc, H.ecore)
    H.extra["E_HF_from_real_integrals"] = hf_energy(h1r, erir, H.nocc, H.ecore)
    H.extra["E_KRHF_supercell"] = float(mf.e_tot) * nk
    H.extra["momentum_violation_max"] = momentum_violation(H)
    H.extra["hermiticity_max"] = float(max(
        np.abs(eri - eri.transpose(1, 0, 3, 2).conj()).max(),
        np.abs(eri - eri.transpose(2, 3, 0, 1)).max(),
        np.abs(h1 - h1.conj().T).max()))
    H.extra["time_reversal_max"] = float(max(
        np.abs(h1[np.ix_(partner, partner)] - h1.conj()).max(),
        np.abs(eri[np.ix_(partner, partner, partner, partner)] - eri.conj()).max()))
    return H


def hf_energy(h1, eri, nocc, ecore):
    o = slice(0, nocc)
    e = 2 * np.trace(h1[o, o]).real
    J = np.einsum("iijj->", eri[o, o, o, o]).real
    K = np.einsum("ijji->", eri[o, o, o, o]).real
    return float(e + 2 * J - K + ecore)


def momentum_violation(H):
    k = H.kidx
    nk = H.nk
    dk = (k[:, None, None, None] - k[None, :, None, None] + k[None, None, :, None]
          - k[None, None, None, :]) % nk
    mask = dk != 0
    v1 = np.abs(H.h1[(k[:, None] - k[None, :]) % nk != 0]).max(initial=0.0)
    return float(max(v1, np.abs(H.eri[mask]).max(initial=0.0)))


def supercell_gamma_check(r, d, nk, basis="sto-3g", log_stem=None):
    """Independent route: Gamma-point RHF + FCI on an explicit N_k-cell supercell."""
    from pyscf.pbc import scf as pscf, tools
    from pyscf import fci, ao2mo
    out = None
    if PBC_LOG["dir"] is not None:
        out = os.path.join(PBC_LOG["dir"], f"{log_stem or PBC_LOG['stem']}__supercell.log")
    cell = make_cell(r, d, basis, verbose=0)
    sc = tools.super_cell(cell, [nk, 1, 1])
    sc.verbose = 4 if out else 0
    if out:
        sc.output = out
        sc.build()
    mf = pscf.RHF(sc, exxdiv=None).density_fit()
    mf.conv_tol = 1e-12
    mf.kernel()
    C = mf.mo_coeff
    nmo = C.shape[1]
    h1 = C.T @ mf.get_hcore() @ C
    eri = mf.with_df.ao2mo(C, compact=False).reshape((nmo,) * 4)
    nel = sc.nelectron // 2
    e, _ = fci.direct_spin1.kernel(h1, eri, nmo, (nel, nel), ecore=sc.energy_nuc(),
                                   conv_tol=1e-12, max_cycle=500)
    return {"E_RHF_supercell": float(mf.e_tot), "E_FCI_supercell": float(e)}


# ======================================================== excitation symmetry
def exc_momentum(e: Excitation, kidx, nk):
    return (sum(int(kidx[p]) for p in e.vir) - sum(int(kidx[p]) for p in e.occ)) % nk


def _canon(kind, spin, occ, vir):
    """Canonical Excitation and the sign relating it to the given index order."""
    sign = 1
    if kind == "Dss":
        if occ[0] > occ[1]:
            occ, sign = (occ[1], occ[0]), -sign
        if vir[0] > vir[1]:
            vir, sign = (vir[1], vir[0]), -sign
        if occ[0] == occ[1] or vir[0] == vir[1]:
            return None, 0
    return Excitation(kind, spin, tuple(occ), tuple(vir)), sign


def spin_partner(e: Excitation):
    return e.spin_flip(), 1


def tr_partner(e: Excitation, partner):
    occ = tuple(int(partner[p]) for p in e.occ)
    vir = tuple(int(partner[p]) for p in e.vir)
    return _canon(e.kind, e.spin, occ, vir)


@dataclass
class Orbit:
    members: list            # [(Excitation, sign, conj)]
    mode: str                # 'complex' | 'real' | 'imag'

    @property
    def n_params(self):
        return 2 if self.mode == "complex" else 1


def build_orbits(pool, kidx, nk, partner, use_k=True, use_spin=True, use_tr=True,
                 allowed=None):
    """Group allowed excitations into orbits of the tying group.

    `allowed(e)` overrides the one-dimensional momentum test (used for 2D/3D
    meshes in the count checks).

    Member relation t_m = sign * (t or conj(t)). An orbit whose relations force
    t = conj(t) (up to sign) carries one real parameter; 'zero' orbits are
    dropped and counted.
    """
    if allowed is None:
        allowed = lambda e: exc_momentum(e, kidx, nk) == 0
    kept = [e for e in pool if (not use_k) or allowed(e)]
    keyset = set(kept)
    seen = set()
    orbits, n_zero = [], 0
    for e0 in kept:
        if e0 in seen:
            continue
        rel = {e0: (1, False)}
        stack = [e0]
        constraint = None          # None, 'real', 'imag', 'zero'
        while stack:
            e = stack.pop()
            s, c = rel[e]
            nbrs = []
            if use_spin:
                f, sf = spin_partner(e)
                nbrs.append((f, s * sf, c))
            if use_tr:
                f, sf = tr_partner(e, partner)
                if f is not None:
                    nbrs.append((f, s * sf, not c))
            for f, sf, cf in nbrs:
                if f not in keyset:
                    raise RuntimeError("tying partner outside the screened pool")
                if f in rel:
                    s2, c2 = rel[f]
                    if (s2, c2) == (sf, cf):
                        continue
                    if c2 == cf:              # s2 X = sf X with s2 != sf -> t = 0
                        new = "zero"
                    else:                     # s2 t = sf t^*  ->  real or imaginary t
                        new = "real" if s2 == sf else "imag"
                    if constraint is None or constraint == new:
                        constraint = new
                    else:
                        constraint = "zero"
                    continue
                rel[f] = (sf, cf)
                stack.append(f)
        seen.update(rel)
        if constraint == "zero":
            n_zero += 1
            continue
        members = [(e, s, c) for e, (s, c) in rel.items()]
        orbits.append(Orbit(members, constraint or "complex"))
    return orbits, n_zero


def periodic_pool(H: PeriodicHamiltonian):
    return uccsd_pool(H.nocc, H.norb)


VARIANTS_K = {"k only": (True, False, False), "k + spin": (True, True, False),
              "k + TR": (True, False, True), "GVA-k": (True, True, True)}


def variant_orbits(H, name):
    return build_orbits(periodic_pool(H), H.kidx, H.nk, H.partner, *VARIANTS_K[name])


def orbit_counts(orbits):
    return {"n_params": int(sum(o.n_params for o in orbits)),
            "n_generators": int(sum(len(o.members) for o in orbits)),
            "cnot_1L": int(sum(cnot_count_ladder(m[0]) for o in orbits for m in o.members)),
            "n_orbits": len(orbits),
            "n_real_orbits": int(sum(o.mode != "complex" for o in orbits))}


# ======================================================= complex simulator
class CDeterminantSpace:
    """Determinant space with a complex Hermitian Hamiltonian (PySCF direct_nosym)."""

    def __init__(self, norb, nelec, h1, eri, ecore):
        self.norb = norb
        self.na, self.nb = nelec
        self.nelec = nelec
        self.ecore = ecore
        self.sa = cistring.num_strings(norb, self.na)
        self.sb = cistring.num_strings(norb, self.nb)
        self.dim = self.sa * self.sb
        self.h2e = direct_nosym.absorb_h1e(np.asarray(h1, complex), np.asarray(eri, complex),
                                           norb, nelec, 0.5)
        self.link = (cistring.gen_linkstr_index(range(norb), self.na),
                     cistring.gen_linkstr_index(range(norb), self.nb))
        # single-excitation matrices (real signs), shared with the real engine
        import scipy.sparse as sp
        self._E = {}
        for spin, n in ((0, self.na), (1, self.nb)):
            link = cistring.gen_linkstr_index(range(norb), n)
            ns = link.shape[0]
            rows, cols, vals, keys = [], [], [], []
            for s0 in range(ns):
                for a, i, s1, sg in link[s0]:
                    keys.append((a, i))
                    rows.append(s1)
                    cols.append(s0)
                    vals.append(sg)
            keys = np.array(keys)
            rows, cols, vals = map(np.array, (rows, cols, vals))
            for a in range(norb):
                for i in range(norb):
                    if a == i:
                        continue
                    m = (keys[:, 0] == a) & (keys[:, 1] == i)
                    self._E[(spin, a, i)] = sp.csr_matrix(
                        (vals[m].astype(float), (rows[m], cols[m])), shape=(ns, ns))

    from .engine import DeterminantSpace as _DS
    pairs = _DS.pairs

    def hamiltonian(self, v):
        c = direct_nosym.contract_2e(self.h2e, np.asarray(v, complex).reshape(self.sa, self.sb),
                                     self.norb, self.nelec, self.link)
        return np.asarray(c).ravel()

    def energy(self, v):
        return float(np.vdot(v, self.hamiltonian(v)).real) + self.ecore

    def reference(self):
        v = np.zeros(self.dim, complex)
        v[0] = 1.0
        return v

    def s2(self, v):
        v = v.reshape(self.sa, self.sb)
        sr = spin_op.contract_ss(np.ascontiguousarray(v.real), self.norb, self.nelec)
        si = spin_op.contract_ss(np.ascontiguousarray(v.imag), self.norb, self.nelec)
        return float(np.vdot(v.ravel(), (sr + 1j * si).ravel()).real)

    def ground_state(self, k=1, tol=1e-12):
        from scipy.sparse.linalg import LinearOperator, eigsh
        op = LinearOperator((self.dim, self.dim), matvec=self.hamiltonian, dtype=complex)
        w, V = eigsh(op, k=k, which="SA", tol=tol, v0=self.reference() + 1e-3)
        o = np.argsort(w)
        return w[o] + self.ecore, V[:, o]


def _sinc_terms(r):
    if r < 1e-6:
        return 1 - r * r / 6, -1 / 3 + r * r / 30
    s, c = np.sin(r), np.cos(r)
    return s / r, (r * c - s) / r ** 3


class CAnsatz:
    """Product of complex-amplitude excitation exponentials, orbit-tied, `layers` times.

    For one excitation with amplitude t: exp(t tau - t^* tau^dagger). On each pair
    (up, dn) with tau|up> = sg|dn> it acts as
        [v_up, v_dn] -> [cos|t| v_up - sg t^* sinc v_dn,  cos|t| v_dn + sg t sinc v_up].
    """

    def __init__(self, space: CDeterminantSpace, orbits, layers=1):
        self.space, self.orbits, self.layers = space, orbits, layers
        tab = {}
        for o in orbits:
            for e, _, _ in o.members:
                if e not in tab:
                    tab[e] = space.pairs(e)
        self.ops = []
        npar = 0
        self.param_map = []
        for L in range(layers):
            for oi, o in enumerate(orbits):
                if o.mode == "complex":
                    px, py = npar, npar + 1
                    npar += 2
                elif o.mode == "real":
                    px, py = npar, -1
                    npar += 1
                else:
                    px, py = -1, npar
                    npar += 1
                self.param_map.append((L + 1, oi, px, py))
                for e, s, c in o.members:
                    up, dn, sg = tab[e]
                    self.ops.append((px, py, float(s), -1.0 if c else 1.0, up, dn, sg))
        self.n_params = npar

    @staticmethod
    def _amp(theta, px, py, s, cy):
        x = theta[px] if px >= 0 else 0.0
        y = theta[py] if py >= 0 else 0.0
        return s * x, s * cy * y

    @staticmethod
    def _apply(v, up, dn, sg, u, w):
        t = u + 1j * w
        r = abs(t)
        if r == 0.0:
            return
        f1, _ = _sinc_terms(r)
        c = np.cos(r)
        vu, vd = v[up], v[dn]
        v[up] = c * vu - sg * (np.conj(t) * f1) * vd
        v[dn] = c * vd + sg * (t * f1) * vu

    def state(self, theta):
        v = self.space.reference()
        for px, py, s, cy, up, dn, sg in self.ops:
            u, w = self._amp(theta, px, py, s, cy)
            self._apply(v, up, dn, sg, u, w)
        return v

    def energy_and_grad(self, theta):
        psi = self.state(theta)
        hpsi = self.space.hamiltonian(psi)
        e = float(np.vdot(psi, hpsi).real) + self.space.ecore
        grad = np.zeros(self.n_params)
        lam = hpsi.copy()
        for px, py, s, cy, up, dn, sg in reversed(self.ops):
            u, w = self._amp(theta, px, py, s, cy)
            self._apply(psi, up, dn, sg, -u, -w)          # psi <- U^-1 psi
            r = np.hypot(u, w)
            f1, f2 = _sinc_terms(r)
            t = u + 1j * w
            pu, pd = psi[up], psi[dn]
            lu, ld = lam[up], lam[dn]
            a13 = np.vdot(lu, pu) + np.vdot(ld, pd)       # <lam|psi> on the pair
            a2 = np.vdot(lu, sg * pd)
            a4 = np.vdot(ld, sg * pu)
            # dM/du_c: diag dc, (up,dn) entry -sg*dtc, (dn,up) entry +sg*dt
            for comp, der in ((0, u), (1, w)):
                dc = -f1 * der
                dt = (1.0 if comp == 0 else 1j) * f1 + t * f2 * der
                dtc = (1.0 if comp == 0 else -1j) * f1 + np.conj(t) * f2 * der
                g = 2.0 * (dc * a13 - dtc * a2 + dt * a4).real
                if comp == 0 and px >= 0:
                    grad[px] += s * g
                if comp == 1 and py >= 0:
                    grad[py] += s * cy * g
            self._apply(lam, up, dn, sg, -u, -w)
        return e, grad


def _minimize_cached(fun_grad, theta0, maxiter, gtol, ftol):
    """L-BFGS-B with the energy history taken from cached evaluations."""
    cache = {"x": None, "f": None}
    hist, nfev = [], [0]

    def fun(t):
        nfev[0] += 1
        f, g = fun_grad(t)
        cache["x"], cache["f"] = np.array(t, copy=True), f
        return f, g

    def cb(xk):
        if cache["x"] is not None and np.array_equal(xk, cache["x"]):
            hist.append(cache["f"])
        else:
            hist.append(fun_grad(xk)[0])

    e0 = fun_grad(theta0)[0]
    res = minimize(fun, theta0, jac=True, method="L-BFGS-B", callback=cb,
                   options={"maxiter": maxiter, "gtol": gtol, "ftol": ftol, "maxcor": 30})
    return res, np.array([e0] + hist), nfev[0]


def run_vqe_c(ansatz, theta0=None, e_exact=None, gtol=1e-7, maxiter=3000,
              chem_acc=1.6e-3, ftol=1e-12):
    """VQE driver for CAnsatz (complex) or the real engine's Ansatz.

    Same optimiser settings as engine.run_vqe; the energy history is read from
    the evaluations L-BFGS-B already made instead of being recomputed.
    """
    if theta0 is None:
        theta0 = np.zeros(ansatz.n_params)
    t0 = time.time()
    res, hist, nfev = _minimize_cached(ansatz.energy_and_grad, theta0, maxiter, gtol, ftol)
    wall = time.time() - t0
    out = {"E": float(res.fun), "theta": res.x, "nit": int(res.nit), "nfev": nfev,
           "wall_s": wall, "history": hist, "success": bool(res.success),
           "message": str(res.message), "grad_norm": float(np.linalg.norm(res.jac))}
    if e_exact is not None:
        out["dE"] = out["E"] - e_exact
        hit = np.where(hist - e_exact <= chem_acc)[0]
        out["iter_to_chem_acc"] = int(hit[0]) if len(hit) else None
    psi = ansatz.state(res.x)
    out["S2"] = float(ansatz.space.s2(psi))
    out["psi"] = psi
    return out


def pool_gradients(space: CDeterminantSpace, psi, excitations):
    """|dE/dt| at t = 0 for excitations appended after psi (2 real components)."""
    hpsi = space.hamiltonian(psi)
    out = []
    for e in excitations:
        up, dn, sg = space.pairs(e)
        # d/du: up -> -sg psi_dn ; dn -> +sg psi_up ; d/dw: up -> +i sg psi_dn, dn -> i sg psi_up
        gu = 2 * (np.vdot(hpsi[up], -sg * psi[dn]) + np.vdot(hpsi[dn], sg * psi[up])).real
        gw = 2 * (np.vdot(hpsi[up], 1j * sg * psi[dn]) + np.vdot(hpsi[dn], 1j * sg * psi[up])).real
        out.append(np.hypot(gu, gw))
    return np.array(out)


def adapt_vqe_c(space: CDeterminantSpace, orbits, e_exact=None, grad_tol=1e-3,
                max_ops=150, e_stall=1e-7, stall_window=5, chem_acc=1.6e-3,
                max_wall=None):
    """ADAPT-VQE over tied orbits (complex amplitudes); full re-optimisation.

    `max_wall` (s): stop before adding a further operator once this wall time is
    exceeded; the state reached so far is returned with stop_reason 'time limit'.
    """
    chosen, theta, trace = [], np.zeros(0), []
    t0 = time.time()
    while True:
        A = CAnsatz(space, chosen, 1) if chosen else None
        psi = A.state(theta) if chosen else space.reference()
        e = space.energy(psi)
        g = np.zeros(len(orbits))
        hpsi = space.hamiltonian(psi)
        for n, o in enumerate(orbits):
            gu = gw = 0.0
            for ex, s, c in o.members:
                up, dn, sg = space.pairs(ex)
                cy = -1.0 if c else 1.0
                du = 2 * (np.vdot(hpsi[up], -sg * psi[dn]) + np.vdot(hpsi[dn], sg * psi[up])).real
                dw = 2 * (np.vdot(hpsi[up], 1j * sg * psi[dn]) + np.vdot(hpsi[dn], 1j * sg * psi[up])).real
                gu += s * du
                gw += s * cy * dw
            if o.mode == "real":
                gw = 0.0
            elif o.mode == "imag":
                gu = 0.0
            g[n] = np.hypot(gu, gw)
        gnorm = float(np.linalg.norm(g))
        trace.append({"n_ops": len(chosen), "E": e, "grad_norm": gnorm,
                      "wall_s": time.time() - t0})
        stalled = len(trace) > stall_window and trace[-stall_window - 1]["E"] - e < e_stall
        timed_out = max_wall is not None and time.time() - t0 > max_wall
        if gnorm < grad_tol or len(chosen) >= max_ops or stalled or timed_out:
            break
        best = int(np.argmax(g))
        chosen = chosen + [orbits[best]]
        A = CAnsatz(space, chosen, 1)
        th0 = np.concatenate([theta, np.zeros(A.n_params - len(theta))])
        res = minimize(A.energy_and_grad, th0, jac=True, method="L-BFGS-B",
                       options={"maxiter": 1500, "gtol": 1e-7, "ftol": 1e-12})
        theta = res.x
    stop = ("gradient" if gnorm < grad_tol else "max_ops" if len(chosen) >= max_ops
            else "stall" if stalled else "time limit")
    out = {"E": e, "n_ops": len(chosen), "trace": trace, "chosen": chosen,
           "stop_reason": stop, "theta": theta, "wall_s": time.time() - t0}
    if e_exact is not None:
        out["dE"] = e - e_exact
    return out


# ============================================================ resource counts
def mesh_points(mesh):
    g = np.array(np.meshgrid(*[np.arange(n) for n in mesh], indexing="ij")).reshape(len(mesh), -1).T
    return g


def periodic_counts(mesh, o, v):
    """Exact parameter and CNOT counts for a crystal with o occupied and v
    virtual bands per cell on a Gamma-centred mesh (no band crossing).

    Returns counts for the real supercell UCCSD and for k-screened variants.
    Orbital order: occupied (k-major, band-minor) then virtual; interleaved spins.
    CNOT model as in the molecular code: single 4(q-p); double 16(q-p+s-r+1).
    """
    mesh = tuple(mesh)
    K = mesh_points(mesh)
    nk = len(K)
    radix = np.array([int(np.prod(mesh[i + 1:])) for i in range(len(mesh))])
    kid = (K * radix).sum(1)

    def add(a, b):
        return ((K[a] + K[b]) % mesh * radix).sum(-1)

    def neg(a):
        return ((-K[a]) % mesh * radix).sum(-1)

    no, nv = nk * o, nk * v
    ko = np.repeat(np.arange(nk), o)       # k index of occupied spatial orbital
    kv = np.repeat(np.arange(nk), v)
    po = np.arange(no)
    pv = no + np.arange(nv)

    out = {"mesh": list(mesh), "n_k": nk, "o": o, "v": v, "n_qubits": 2 * (no + nv)}
    # ---------------- real supercell UCCSD (all excitations)
    from math import comb
    out["UCCSD_params"] = 2 * no * nv + 2 * comb(no, 2) * comb(nv, 2) + (no * nv) ** 2

    def pair_tables(pos, kk, same_spin):
        if same_spin:
            i, j = np.triu_indices(len(pos), 1)
            span = 2 * (pos[j] - pos[i])
        else:   # (alpha at x, beta at y), all ordered pairs
            i, j = np.meshgrid(np.arange(len(pos)), np.arange(len(pos)), indexing="ij")
            i, j = i.ravel(), j.ravel()
            span = np.abs((2 * pos[j] + 1) - 2 * pos[i])
        Ksum = add(kk[i], kk[j])
        return Ksum, span, i, j

    # singles
    ia_k = (ko[:, None] == kv[None, :])
    span_s = 2 * (pv[None, :] - po[:, None])
    cnot_s_all = 2 * 4 * span_s.sum()
    cnot_s_k = 2 * 4 * span_s[ia_k].sum()
    S = int(ia_k.sum())
    tabs = {}
    for name, ss in (("ss", True), ("ab", False)):
        Ko, so, _, _ = pair_tables(po, ko, ss)
        Kv, sv, _, _ = pair_tables(pv, kv, ss)
        no_K = np.bincount(Ko, minlength=nk).astype(float)
        nv_K = np.bincount(Kv, minlength=nk).astype(float)
        So_K = np.bincount(Ko, weights=so, minlength=nk)
        Sv_K = np.bincount(Kv, weights=sv, minlength=nk)
        tot_o, tot_v = no_K.sum(), nv_K.sum()
        all_cnot = 16 * (tot_v * So_K.sum() + tot_o * Sv_K.sum() + tot_o * tot_v)
        k_cnot = 16 * (nv_K * So_K + no_K * Sv_K + no_K * nv_K).sum()
        tabs[name] = (no_K, nv_K, all_cnot, k_cnot)
    mult_ss = 2
    out["UCCSD_cnot_1L"] = int(cnot_s_all + mult_ss * tabs["ss"][2] + tabs["ab"][2])
    D_ss = int((tabs["ss"][0] * tabs["ss"][1]).sum())
    D_ab = int((tabs["ab"][0] * tabs["ab"][1]).sum())
    # self spin-partners: i=j, a=b with 2(k_a - k_i) = 0
    two_zero = np.array([np.all((2 * K[q]) % mesh == 0) for q in range(nk)])
    dq = ((K[kv][None, :, :] - K[ko][:, None, :]) % mesh * radix).sum(-1)
    n_self = int(two_zero[dq].sum())
    gen_k = 2 * S + 2 * D_ss + D_ab
    out.update({"S_per_spin": S, "D_ss_per_spin": D_ss, "D_ab": D_ab, "n_self": n_self,
                "k_generators": gen_k,
                "k_only_params": 2 * gen_k,
                "k_spin_params": 2 * (S + D_ss + (D_ab + n_self) // 2),
                "k_TR_params": gen_k,
                "GVA_k_params": S + D_ss + (D_ab + n_self) // 2,
                "k_cnot_1L": int(cnot_s_k + mult_ss * tabs["ss"][3] + tabs["ab"][3])})
    return out


# ============================================================ qubit mapping
def qubit_hamiltonian_c(h1, eri, ecore, prune=None, exact=True):
    """Jordan-Wigner image of a (possibly complex) Hermitian Hamiltonian.

    OpenFermion removes terms whose coefficient falls below its EQ_TOLERANCE
    (1e-8) while operators are accumulated, and InteractionOperator inputs built
    with spinorb_from_spatial are truncated at the same level. Here that
    pruning is switched off during construction, so the operator is exact;
    `prune` (e.g. 1e-8) removes small Pauli terms explicitly afterwards, which
    is what the Pauli-term counts use. exact=False keeps OpenFermion's default
    pruning during construction (much faster for 20 qubits; used for the counts,
    as for the molecular Hamiltonians).
    Returns the dictionary format of gva.qubit (symplectic masks, real coefficients).
    """
    from openfermion import InteractionOperator, jordan_wigner
    import openfermion.ops.operators.symbolic_operator as so
    n = h1.shape[0]
    two = eri.transpose(0, 2, 3, 1)                   # (ps|qr) -> two[p,q,r,s]
    h1s = np.zeros((2 * n, 2 * n), complex)
    h2s = np.zeros((2 * n,) * 4, complex)
    for s_ in (0, 1):
        h1s[s_::2, s_::2] = h1
    for a in (0, 1):
        for b in (0, 1):
            h2s[a::2, b::2, b::2, a::2] = two
    orig = so.SymbolicOperator.__dict__["_issmall"]
    if exact:
        so.SymbolicOperator._issmall = staticmethod(
            lambda val, tol=0.0, _f=orig.__func__: _f(val, tol))
    try:
        qop = jordan_wigner(InteractionOperator(ecore, h1s, 0.5 * h2s))
    finally:
        so.SymbolicOperator._issmall = orig
    if prune:
        qop.compress(prune)
    nq = 2 * n
    const, coefs, xm, zm, ny, imag_max = 0.0, [], [], [], [], 0.0
    for term, c in qop.terms.items():
        imag_max = max(imag_max, abs(c.imag))
        if len(term) == 0:
            const += c.real
            continue
        x = z = y = 0
        for q, p in term:
            bit = 1 << (nq - 1 - q)
            if p == "X":
                x |= bit
            elif p == "Z":
                z |= bit
            else:
                x |= bit
                z |= bit
                y += 1
        coefs.append(c.real)
        xm.append(x)
        zm.append(z)
        ny.append(y)
    assert imag_max < 1e-10, imag_max
    return {"n_qubits": nq, "const": const, "coef": np.array(coefs),
            "x": np.array(xm, dtype=np.int64), "z": np.array(zm, dtype=np.int64),
            "ny": np.array(ny), "imag_max": imag_max}


def civec_to_statevector_c(civec, norb, nelec):
    """Complex version of gva.qubit.civec_to_statevector (same sign rule)."""
    na, nb = nelec
    n = 2 * norb
    civec = np.asarray(civec).reshape(cistring.num_strings(norb, na), -1)
    sa = cistring.make_strings(range(norb), na)
    sb = cistring.make_strings(range(norb), nb)
    psi = np.zeros(2 ** n, complex)
    for ia, a in enumerate(sa):
        occ_a = [p for p in range(norb) if (a >> p) & 1]
        idx_a = sum(1 << (n - 1 - 2 * p) for p in occ_a)
        for ib, b in enumerate(sb):
            c = civec[ia, ib]
            if c == 0:
                continue
            occ_b = [p for p in range(norb) if (b >> p) & 1]
            inv = sum(1 for j in occ_b for i in occ_a if i > j)
            idx = idx_a + sum(1 << (n - 1 - (2 * p + 1)) for p in occ_b)
            psi[idx] = c * (-1) ** inv
    return psi


def momentum_sector_weights(H, civec):
    """Weight of a CI vector (Bloch basis) in each total-momentum sector."""
    sa = cistring.make_strings(range(H.norb), H.nelec[0])
    sb = cistring.make_strings(range(H.norb), H.nelec[1])
    ka = np.array([sum(int(H.kidx[p]) for p in range(H.norb) if (s >> p) & 1) for s in sa])
    kb = np.array([sum(int(H.kidx[p]) for p in range(H.norb) if (s >> p) & 1) for s in sb])
    Ktot = (ka[:, None] + kb[None, :]) % H.nk
    w = np.abs(np.asarray(civec).reshape(len(sa), len(sb))) ** 2
    return np.bincount(Ktot.ravel(), weights=w.ravel(), minlength=H.nk), int(Ktot[0, 0])


def synthetic_mesh_structure(mesh, o, v):
    """k labels, time-reversal partners and a momentum test for an o/v-band
    insulator on `mesh` (used to check the closed-form counts by enumeration)."""
    K = mesh_points(tuple(mesh))
    nk = len(K)
    mesh = np.array(mesh)
    kv_of = {tuple(k): i for i, k in enumerate(K)}
    orbs = [(j, n) for j in range(nk) for n in range(o)] + \
           [(j, n) for j in range(nk) for n in range(o, o + v)]
    pos = {x: p for p, x in enumerate(orbs)}
    kidx = np.array([x[0] for x in orbs])
    partner = np.array([pos[(kv_of[tuple((-K[j]) % mesh)], n)] for j, n in orbs])
    kvec = K[kidx]

    def allowed(e):
        tot = sum(kvec[p] for p in e.vir) - sum(kvec[p] for p in e.occ)
        return bool(np.all(tot % mesh == 0))
    return kidx, partner, allowed, nk * o, nk * (o + v)
