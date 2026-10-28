"""Excitation pools, GVA construction (point-group screening + spin tying),
and Jordan-Wigner CNOT counting.

Spin-orbital convention: spatial orbital p, spin s (0=alpha, 1=beta) maps to
qubit 2p+s (interleaved ordering, as in OpenFermion).
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb


@dataclass(frozen=True)
class Excitation:
    kind: str          # 'S' single, 'Dss' same-spin double, 'Dab' alpha-beta double
    spin: int          # for S/Dss: 0/1 ; for Dab: -1
    occ: tuple         # spatial occupied indices (i,) or (i, j)
    vir: tuple         # spatial virtual indices (a,) or (a, b)

    def spin_orbitals(self):
        """Return (annihilated spin-orbitals, created spin-orbitals)."""
        if self.kind == "S":
            (i,), (a,) = self.occ, self.vir
            return (2 * i + self.spin,), (2 * a + self.spin,)
        if self.kind == "Dss":
            (i, j), (a, b) = self.occ, self.vir
            s = self.spin
            return (2 * i + s, 2 * j + s), (2 * a + s, 2 * b + s)
        (i, j), (a, b) = self.occ, self.vir     # a,i alpha ; b,j beta
        return (2 * i, 2 * j + 1), (2 * a, 2 * b + 1)

    def irrep_product(self, orbsym):
        g = 0
        for p in self.occ + self.vir:
            g ^= int(orbsym[p])
        return g

    def spin_flip(self):
        if self.kind in ("S", "Dss"):
            return Excitation(self.kind, 1 - self.spin, self.occ, self.vir)
        (i, j), (a, b) = self.occ, self.vir
        return Excitation("Dab", -1, (j, i), (b, a))


def uccsd_pool(nocc: int, norb: int):
    """Spin-conserving UCCSD generators, ordered singles then doubles."""
    occ = range(nocc)
    vir = range(nocc, norb)
    pool = []
    for s in (0, 1):
        for i in occ:
            for a in vir:
                pool.append(Excitation("S", s, (i,), (a,)))
    for s in (0, 1):
        for i, j in combinations(occ, 2):
            for a, b in combinations(vir, 2):
                pool.append(Excitation("Dss", s, (i, j), (a, b)))
    for i in occ:
        for j in occ:
            for a in vir:
                for b in vir:
                    pool.append(Excitation("Dab", -1, (i, j), (a, b)))
    return pool


def uccsd_parameter_count(no: int, nv: int) -> int:
    """N = 2*no*nv + 2*C(no,2)*C(nv,2) + (no*nv)^2  (spin-conserving UCCSD)."""
    return 2 * no * nv + 2 * comb(no, 2) * comb(nv, 2) + (no * nv) ** 2


def gva_groups(pool, orbsym, use_pg=True, use_spin=True):
    """Return list of tied groups (each a tuple of Excitation objects).

    1. Point-group screening: keep generators whose irrep product (XOR of the
       abelian irrep ids of the four or two orbitals) is totally symmetric.
    2. Spin tying: a generator and its alpha<->beta partner share one parameter.
       Self-partnered alpha-beta doubles (i=j and a=b) keep a single parameter.
    Partners are placed next to each other in the product.
    """
    kept = [e for e in pool if (not use_pg) or e.irrep_product(orbsym) == 0]
    if not use_spin:
        return [(e,) for e in kept]
    keyset = set(kept)
    seen = set()
    groups = []
    for e in kept:
        if e in seen:
            continue
        f = e.spin_flip()
        if f == e or f not in keyset:
            groups.append((e,))
            seen.add(e)
        else:
            groups.append((e, f))
            seen.update((e, f))
    return groups


# ------------------------------------------------------------ CNOT counting
def jw_pauli_weights(exc: Excitation):
    """Weights of the Pauli strings in the JW image of (T - T^dagger).

    Single p<q: 2 strings of weight q-p+1.
    Double with sorted spin-orbitals p<q<r<s: 8 strings of weight
    (q-p+1)+(s-r+1) (the Z chain between q and r cancels).
    """
    ann, cre = exc.spin_orbitals()
    idx = sorted(ann + cre)
    if len(idx) == 2:
        w = idx[1] - idx[0] + 1
        return [w] * 2
    p, q, r, s = idx
    w = (q - p + 1) + (s - r + 1)
    return [w] * 8


def cnot_count_ladder(exc: Excitation) -> int:
    """Standard CNOT-ladder compilation: 2(w-1) CNOTs per Pauli exponential."""
    return sum(2 * (w - 1) for w in jw_pauli_weights(exc))


def cnot_count_groups(groups, layers=1) -> int:
    return layers * sum(cnot_count_ladder(e) for g in groups for e in g)


def counts_summary(nocc, norb, orbsym):
    pool = uccsd_pool(nocc, norb)
    out = {}
    for name, pg, sp in (("UCCSD", False, False), ("PG only", True, False),
                         ("spin only", False, True), ("GVA", True, True)):
        g = gva_groups(pool, orbsym, pg, sp)
        out[name] = {"n_params": len(g),
                     "n_generators": sum(len(x) for x in g),
                     "cnot_1L": cnot_count_groups(g)}
    out["formula_UCCSD"] = uccsd_parameter_count(nocc, norb - nocc)
    return out


# ------------------------------------------------------------ labels
def excitation_label(e: Excitation) -> str:
    """Spin-orbital label, e.g. 'D 2a,3b->7a,9b' (orbital index + spin)."""
    ann, cre = e.spin_orbitals()
    f = lambda q: f"{q // 2}{'ab'[q % 2]}"
    return f"{'S' if e.kind == 'S' else 'D'} {','.join(map(f, ann))}->{','.join(map(f, cre))}"


def parameter_map(groups, layers=1):
    """Label of every variational parameter, in the order used by Ansatz."""
    out = []
    for L in range(layers):
        for g in groups:
            out.append({"layer": L + 1, "generators": [excitation_label(e) for e in g]})
    return out
