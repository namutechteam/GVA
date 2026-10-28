"""Jordan-Wigner qubit Hamiltonians, Pauli grouping and shot estimates.

Pauli strings are stored in symplectic form (x-mask, z-mask) with qubit q on
bit (n-1-q) of the computational-basis index (OpenFermion ordering).
"""
from __future__ import annotations

import numpy as np
from openfermion import InteractionOperator, jordan_wigner
from openfermion.chem.molecular_data import spinorb_from_spatial
from pyscf.fci import cistring


def qubit_hamiltonian(H):
    h1s, h2s = spinorb_from_spatial(H.h1, H.eri.transpose(0, 2, 3, 1))
    op = InteractionOperator(H.ecore, h1s, 0.5 * h2s)
    qop = jordan_wigner(op)
    qop.compress(1e-12)
    n = 2 * H.norb
    const = 0.0
    coefs, xm, zm, ny = [], [], [], []
    for term, c in qop.terms.items():
        if len(term) == 0:
            const += c.real
            continue
        x = z = 0
        y = 0
        for q, p in term:
            bit = 1 << (n - 1 - q)
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
    return {"n_qubits": n, "const": const, "coef": np.array(coefs),
            "x": np.array(xm, dtype=np.int64), "z": np.array(zm, dtype=np.int64),
            "ny": np.array(ny)}


def civec_to_statevector(civec, norb, nelec):
    """Map a PySCF CI vector (alpha-string x beta-string) to the JW basis.

    PySCF orders a determinant as (alpha creators ascending)(beta creators
    ascending)|0>; JW/OpenFermion uses ascending interleaved spin-orbitals.
    The reordering sign is (-1)^{#pairs (i in alpha, j in beta) with i > j}.
    """
    na, nb = nelec
    n = 2 * norb
    sa = cistring.make_strings(range(norb), na)
    sb = cistring.make_strings(range(norb), nb)
    psi = np.zeros(2 ** n)
    for ia, a in enumerate(sa):
        occ_a = [p for p in range(norb) if (a >> p) & 1]
        idx_a = sum(1 << (n - 1 - 2 * p) for p in occ_a)
        for ib, b in enumerate(sb):
            c = civec[ia, ib]
            if c == 0.0:
                continue
            occ_b = [p for p in range(norb) if (b >> p) & 1]
            inv = sum(1 for j in occ_b for i in occ_a if i > j)
            idx = idx_a + sum(1 << (n - 1 - (2 * p + 1)) for p in occ_b)
            psi[idx] = c * (-1) ** inv
    return psi


def _popcount(a):
    return np.bitwise_count(a)


def apply_pauli(psi, x, z, ny, idx=None):
    """P|b> = i^ny (-1)^{popcount(z&b)} |b^x>, returned as complex array."""
    if idx is None:
        idx = np.arange(psi.shape[0], dtype=np.int64)
    phase = 1.0 - 2.0 * (_popcount(idx & z) & 1)
    out = np.empty_like(psi, dtype=complex)
    out[idx ^ x] = (1j ** ny) * phase * psi
    return out


def expectation(qh, psi):
    idx = np.arange(psi.shape[0], dtype=np.int64)
    e = qh["const"]
    vals = np.zeros(len(qh["coef"]))
    for k in range(len(qh["coef"])):
        v = np.vdot(psi, apply_pauli(psi, qh["x"][k], qh["z"][k],
                                     qh["ny"][k], idx)).real
        vals[k] = v
        e += qh["coef"][k] * v
    return e, vals


# ------------------------------------------------------------- grouping
def conflict_matrix(qh, kind):
    x, z = qh["x"], qh["z"]
    m = len(x)
    C = np.zeros((m, m), dtype=bool)
    for i in range(m):
        if kind == "GC":           # anticommute -> conflict
            C[i] = (_popcount((x[i] & z) ^ (z[i] & x)) & 1).astype(bool)
        elif kind == "QWC":
            s_i = x[i] | z[i]
            s = x | z
            C[i] = (((x[i] ^ x) | (z[i] ^ z)) & s_i & s) != 0
        else:
            raise ValueError(kind)
    np.fill_diagonal(C, False)
    return C


def greedy_coloring(C, order):
    m = C.shape[0]
    color = -np.ones(m, dtype=int)
    for v in order:
        used = set(color[C[v] & (color >= 0)])
        c = 0
        while c in used:
            c += 1
        color[v] = c
    return color


def dsatur_coloring(C):
    m = C.shape[0]
    deg = C.sum(1)
    color = -np.ones(m, dtype=int)
    nb_colors = [set() for _ in range(m)]
    sat = np.zeros(m, dtype=int)
    for _ in range(m):
        cand = np.where(color < 0)[0]
        key = sat[cand] * (m + 1) + deg[cand]
        v = cand[np.argmax(key)]
        c = 0
        while c in nb_colors[v]:
            c += 1
        color[v] = c
        for u in np.where(C[v] & (color < 0))[0]:
            if c not in nb_colors[u]:
                nb_colors[u].add(c)
                sat[u] += 1
    return color


def sorted_insertion(qh, C):
    """Crawford et al. (Quantum 5, 385): insert terms by |coefficient|."""
    order = np.argsort(-np.abs(qh["coef"]))
    groups = []
    color = -np.ones(len(order), dtype=int)
    for v in order:
        placed = False
        for g, members in enumerate(groups):
            if not C[v, members].any():
                members.append(v)
                color[v] = g
                placed = True
                break
        if not placed:
            groups.append([v])
            color[v] = len(groups) - 1
    return color


def group_terms(qh, kind):
    C = conflict_matrix(qh, kind)
    deg = C.sum(1)
    cols = {"largest_first": greedy_coloring(C, np.argsort(-deg, kind="stable")),
            "dsatur": dsatur_coloring(C),
            "sorted_insertion": sorted_insertion(qh, C)}
    return cols


def group_sigmas(qh, psi, color):
    """sigma_g = sqrt(Var(H_g)) for each group under state psi."""
    idx = np.arange(psi.shape[0], dtype=np.int64)
    ng = color.max() + 1
    sig = np.zeros(ng)
    for g in range(ng):
        members = np.where(color == g)[0]
        hv = np.zeros_like(psi, dtype=complex)
        for k in members:
            hv += qh["coef"][k] * apply_pauli(psi, qh["x"][k], qh["z"][k],
                                              qh["ny"][k], idx)
        mean = np.vdot(psi, hv).real
        sq = np.vdot(hv, hv).real
        sig[g] = np.sqrt(max(sq - mean ** 2, 0.0))
    return sig


def shots_for_precision(sigmas, eps):
    """Optimal allocation n_g ~ sigma_g: N = (sum_g sigma_g)^2 / eps^2."""
    return float(np.sum(sigmas) ** 2 / eps ** 2)
