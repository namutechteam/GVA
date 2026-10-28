"""Exact (noiseless) simulation of Trotterised UCC-type ansaetze in the
fixed-(N_alpha, N_beta) determinant space, with adjoint gradients.

Each fermionic generator A = T - T^dagger (T = a^dag a^dag a a or a^dag a with
disjoint index sets) satisfies A^3 = -A, so exp(theta A) is a set of Givens
rotations on determinant pairs (|up>, |down>) with T|up> = s|down>:
    up'   = cos(theta) up - s sin(theta) down
    down' = cos(theta) down + s sin(theta) up
Sign conventions are taken from PySCF's string link tables, so the Hamiltonian
action (pyscf.fci.direct_spin1) and the generators are mutually consistent.
"""
from __future__ import annotations

import time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
from pyscf.fci import cistring, direct_spin1, spin_op

from .ansatz import Excitation


class DeterminantSpace:
    def __init__(self, norb, nelec, h1, eri, ecore):
        self.norb = norb
        self.na, self.nb = nelec
        self.nelec = nelec
        self.ecore = ecore
        self.h1, self.eri = h1, eri
        self.sa = cistring.num_strings(norb, self.na)
        self.sb = cistring.num_strings(norb, self.nb)
        self.dim = self.sa * self.sb
        self.h2e = direct_spin1.absorb_h1e(h1, eri, norb, nelec, 0.5)
        self._E = {}                                   # single-excitation mats
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
                        (vals[m].astype(float), (rows[m], cols[m])),
                        shape=(ns, ns))

    # ---------------------------------------------------------------- H
    def hamiltonian(self, v):
        c = direct_spin1.contract_2e(self.h2e, v.reshape(self.sa, self.sb),
                                     self.norb, self.nelec)
        return c.ravel()

    def energy(self, v):
        return float(v @ self.hamiltonian(v)) + self.ecore

    def reference(self):
        v = np.zeros(self.dim)
        v[0] = 1.0                                     # lowest strings = HF
        return v

    def s2(self, v):
        return spin_op.spin_square(v.reshape(self.sa, self.sb), self.norb,
                                   self.nelec)[0]

    # --------------------------------------------------- generator tables
    def pairs(self, exc: Excitation):
        """Index arrays (up, down, sign) with T|up> = sign|down>."""
        sa, sb = self.sa, self.sb
        if exc.kind == "S":
            (i,), (a,) = exc.occ, exc.vir
            M = self._E[(exc.spin, a, i)].tocoo()
            r, c, s = M.row, M.col, M.data
            if exc.spin == 0:
                ob = np.arange(sb)
                up = (c[:, None] * sb + ob[None, :]).ravel()
                dn = (r[:, None] * sb + ob[None, :]).ravel()
                sg = np.repeat(s, sb)
            else:
                oa = np.arange(sa)
                up = (oa[:, None] * sb + c[None, :]).ravel()
                dn = (oa[:, None] * sb + r[None, :]).ravel()
                sg = np.tile(s, sa)
        elif exc.kind == "Dss":
            (i, j), (a, b) = exc.occ, exc.vir
            E1 = self._E[(exc.spin, a, i)]
            E2 = self._E[(exc.spin, b, j)]
            M = (E1 @ E2).tocoo()          # a^dag_a a^dag_b a_j a_i = E_ai E_bj
            M.eliminate_zeros()
            r, c, s = M.row, M.col, M.data
            if exc.spin == 0:
                ob = np.arange(sb)
                up = (c[:, None] * sb + ob[None, :]).ravel()
                dn = (r[:, None] * sb + ob[None, :]).ravel()
                sg = np.repeat(s, sb)
            else:
                oa = np.arange(sa)
                up = (oa[:, None] * sb + c[None, :]).ravel()
                dn = (oa[:, None] * sb + r[None, :]).ravel()
                sg = np.tile(s, sa)
        else:                               # alpha-beta: E^a_ai (x) E^b_bj
            (i, j), (a, b) = exc.occ, exc.vir
            Ma = self._E[(0, a, i)].tocoo()
            Mb = self._E[(1, b, j)].tocoo()
            up = (Ma.col[:, None] * sb + Mb.col[None, :]).ravel()
            dn = (Ma.row[:, None] * sb + Mb.row[None, :]).ravel()
            sg = (Ma.data[:, None] * Mb.data[None, :]).ravel()
        return up.astype(np.int64), dn.astype(np.int64), sg.astype(float)


class Ansatz:
    """Product of tied generator groups, repeated `layers` times."""

    def __init__(self, space: DeterminantSpace, groups, layers=1):
        self.space = space
        self.groups = groups
        self.layers = layers
        tab = {}
        for g in groups:
            for e in g:
                if e not in tab:
                    tab[e] = space.pairs(e)
        self.ops = []                    # sequence of (param_index, up, dn, sg)
        for L in range(layers):
            for k, g in enumerate(groups):
                for e in g:
                    up, dn, sg = tab[e]
                    self.ops.append((L * len(groups) + k, up, dn, sg))
        self.n_params = layers * len(groups)

    @staticmethod
    def _rot(v, up, dn, sg, c, s):
        vu = v[up]
        vd = v[dn]
        v[up] = c * vu - s * sg * vd
        v[dn] = c * vd + s * sg * vu

    def state(self, theta):
        v = self.space.reference()
        for k, up, dn, sg in self.ops:
            t = theta[k]
            if t != 0.0:
                self._rot(v, up, dn, sg, np.cos(t), np.sin(t))
        return v

    def energy_and_grad(self, theta):
        psi = self.state(theta)
        hpsi = self.space.hamiltonian(psi)
        e = float(psi @ hpsi) + self.space.ecore
        grad = np.zeros(self.n_params)
        lam = hpsi.copy()
        for k, up, dn, sg in reversed(self.ops):
            # dE/dtheta = 2 <lam| A psi_k>, (A psi)[dn] += sg psi[up],
            #                                (A psi)[up] -= sg psi[dn]
            grad[k] += 2.0 * (lam[dn] @ (sg * psi[up]) - lam[up] @ (sg * psi[dn]))
            t = theta[k]
            if t != 0.0:
                c, s = np.cos(t), np.sin(t)
                self._rot(psi, up, dn, sg, c, -s)      # undo U_k
                self._rot(lam, up, dn, sg, c, -s)
        return e, grad


def run_vqe(ansatz: Ansatz, theta0=None, e_exact=None, gtol=1e-7,
            maxiter=3000, chem_acc=1.6e-3, record=True, ftol=1e-12):
    if theta0 is None:
        theta0 = np.zeros(ansatz.n_params)
    hist = []
    nfev = [0]

    def fun(t):
        nfev[0] += 1
        return ansatz.energy_and_grad(t)

    def cb(xk):
        if record:
            hist.append(ansatz.energy_and_grad(xk)[0])

    t0 = time.time()
    e0 = ansatz.energy_and_grad(theta0)[0]
    res = minimize(fun, theta0, jac=True, method="L-BFGS-B", callback=cb,
                   options={"maxiter": maxiter, "gtol": gtol, "ftol": ftol,
                            "maxcor": 30})
    wall = time.time() - t0
    hist = [e0] + hist
    out = {"E": float(res.fun), "theta": res.x, "nit": int(res.nit),
           "nfev": nfev[0], "wall_s": wall, "history": np.array(hist),
           "success": bool(res.success), "message": str(res.message),
           "grad_norm": float(np.linalg.norm(res.jac))}
    if e_exact is not None:
        err = out["history"] - e_exact
        out["dE"] = out["E"] - e_exact
        hit = np.where(err <= chem_acc)[0]
        out["iter_to_chem_acc"] = int(hit[0]) if len(hit) else None
    psi = ansatz.state(res.x)
    out["S2"] = float(ansatz.space.s2(psi))
    return out


def adapt_vqe(space: DeterminantSpace, pool_groups, e_exact=None,
              grad_tol=1e-3, max_ops=150, chem_acc=1.6e-3, verbose=False,
              e_stall=1e-7, stall_window=5):
    """Fermionic ADAPT-VQE with a (tied) operator pool; full re-optimisation.

    Stops when the pool-gradient norm drops below grad_tol, when the energy
    gained over the last `stall_window` additions is below e_stall, or at
    max_ops operators.
    """
    tab = {e: space.pairs(e) for g in pool_groups for e in g}
    chosen = []
    theta = np.zeros(0)
    trace = []
    t0 = time.time()
    while True:
        ans = Ansatz.__new__(Ansatz)
        ans.space, ans.groups, ans.layers = space, chosen, 1
        ans.ops = [(k, *tab[e]) for k, g in enumerate(chosen) for e in g]
        ans.n_params = len(chosen)
        psi = ans.state(theta) if chosen else space.reference()
        hpsi = space.hamiltonian(psi)
        e = float(psi @ hpsi) + space.ecore
        # pool gradients at theta_new = 0 appended at the end: 2<Hpsi|A psi>
        g = np.zeros(len(pool_groups))
        for n, grp in enumerate(pool_groups):
            for ex in grp:
                up, dn, sg = tab[ex]
                g[n] += 2.0 * (hpsi[dn] @ (sg * psi[up]) - hpsi[up] @ (sg * psi[dn]))
        gnorm = float(np.linalg.norm(g))
        trace.append({"n_ops": len(chosen), "E": e, "grad_norm": gnorm,
                      "wall_s": time.time() - t0})
        if verbose:
            print(len(chosen), e - (e_exact or 0), gnorm, flush=True)
        stalled = (len(trace) > stall_window and
                   trace[-stall_window - 1]["E"] - e < e_stall)
        if gnorm < grad_tol or len(chosen) >= max_ops or stalled:
            break
        best = int(np.argmax(np.abs(g)))
        chosen = chosen + [pool_groups[best]]
        theta = np.concatenate([theta, [0.0]])
        ans = Ansatz.__new__(Ansatz)
        ans.space, ans.groups, ans.layers = space, chosen, 1
        ans.ops = [(k, *tab[e]) for k, gg in enumerate(chosen) for e in gg]
        ans.n_params = len(chosen)
        res = minimize(ans.energy_and_grad, theta, jac=True, method="L-BFGS-B",
                       options={"maxiter": 1500, "gtol": 1e-7, "ftol": 1e-12})
        theta = res.x
    stop = ("gradient" if gnorm < grad_tol else
            "max_ops" if len(chosen) >= max_ops else "stall")
    out = {"E": e, "n_ops": len(chosen), "trace": trace, "chosen": chosen,
           "stop_reason": stop,
           "theta": theta, "wall_s": time.time() - t0}
    if e_exact is not None:
        out["dE"] = e - e_exact
    return out
