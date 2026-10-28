"""Closed-form parameter and CNOT counts for UCCSD and GVA, checked against
explicit enumeration.

Notation: n_o occupied and n_v virtual spatial orbitals (closed-shell
reference, N_alpha = N_beta = n_o). For an abelian point group, o_G and v_G are
the numbers of occupied and virtual orbitals in irrep G.

(1) Spin-conserving UCCSD
    N_UCCSD = 2 n_o n_v + 2 C(n_o,2) C(n_v,2) + (n_o n_v)^2
              singles    same-spin doubles     alpha-beta doubles
    For large n_o, n_v:  N_UCCSD ~ (3/2) (n_o n_v)^2  = O(M^4).

(2) Spin-flip tying only (alpha<->beta partners share one parameter)
    N_spin = n_o n_v + C(n_o,2) C(n_v,2) + n_o n_v (n_o n_v + 1) / 2

(3) Point-group screening only (generator kept if the irrep product is A1)
    N_PG = 2 S + 2 D_ss + D_ab, with
      S    = sum_G o_G v_G
      D_ab = sum_G [ X(G) ]^2,      X(G) = sum_{G1 x G2 = G} o_G1 v_G2
      D_ss = sum_G P_o(G) P_v(G),   P_o(G) = #{i<j : G_i x G_j = G}
                                    P_v(G) = #{a<b : G_a x G_b = G}

(4) GVA (screening + tying)
    N_GVA = S + D_ss + (D_ab + S_self) / 2,
    where S_self = sum_G' ... = number of self-partnered alpha-beta doubles
    (i = j, a = b), all of which are totally symmetric:  S_self = n_o n_v.

(5) CNOT count, Jordan-Wigner with interleaved spin orbitals and the standard
    CNOT-ladder circuit for each Pauli exponential (2(w-1) CNOTs for weight w):
      single  (spin orbitals p < q):            4 (q - p)
      double  (sorted spin orbitals p<q<r<s):   16 (q - p + s - r + 1)
"""
import sys
from collections import Counter
from itertools import combinations
from math import comb

from gva.ansatz import uccsd_pool, gva_groups, cnot_count_groups


def n_uccsd(no, nv):
    return 2 * no * nv + 2 * comb(no, 2) * comb(nv, 2) + (no * nv) ** 2


def n_spin(no, nv):
    return no * nv + comb(no, 2) * comb(nv, 2) + no * nv * (no * nv + 1) // 2


def n_pg_terms(occ_irreps, vir_irreps):
    oc, vc = Counter(occ_irreps), Counter(vir_irreps)
    irreps = set(oc) | set(vc) | {g1 ^ g2 for g1 in list(oc) + list(vc)
                                  for g2 in list(oc) + list(vc)}
    S = sum(oc[g] * vc[g] for g in irreps)
    X = Counter()
    for g1, n1 in oc.items():
        for g2, n2 in vc.items():
            X[g1 ^ g2] += n1 * n2
    D_ab = sum(v * v for v in X.values())
    Po = Counter(a ^ b for a, b in combinations(occ_irreps, 2))
    Pv = Counter(a ^ b for a, b in combinations(vir_irreps, 2))
    D_ss = sum(Po[g] * Pv[g] for g in Po)
    return S, D_ss, D_ab


def n_pg(occ_irreps, vir_irreps):
    S, D_ss, D_ab = n_pg_terms(occ_irreps, vir_irreps)
    return 2 * S + 2 * D_ss + D_ab


def n_gva(occ_irreps, vir_irreps):
    S, D_ss, D_ab = n_pg_terms(occ_irreps, vir_irreps)
    no, nv = len(occ_irreps), len(vir_irreps)
    return S + D_ss + (D_ab + no * nv) // 2


def check(orbsym, nocc, label):
    norb = len(orbsym)
    pool = uccsd_pool(nocc, norb)
    occ, vir = list(orbsym[:nocc]), list(orbsym[nocc:])
    enum = {k: len(gva_groups(pool, orbsym, pg, sp)) for k, (pg, sp) in
            {"UCCSD": (0, 0), "spin": (0, 1), "PG": (1, 0), "GVA": (1, 1)}.items()}
    form = {"UCCSD": n_uccsd(nocc, norb - nocc), "spin": n_spin(nocc, norb - nocc),
            "PG": n_pg(occ, vir), "GVA": n_gva(occ, vir)}
    cn = {k: cnot_count_groups(gva_groups(pool, orbsym, pg, sp)) for k, (pg, sp) in
          {"UCCSD": (0, 0), "GVA": (1, 1)}.items()}
    ok = enum == form
    print(f"{label:28s} enum={enum} formula={form} CNOT/layer={cn} {'OK' if ok else 'MISMATCH'}")
    return ok


if __name__ == "__main__":
    import numpy as np
    from gva.systems import small_molecule_hamiltonian
    allok = True
    for m in ("H2", "H4", "LiH", "BeH2"):
        H = small_molecule_hamiltonian(m, 1.0)
        allok &= check(np.array(H.orbsym), H.nocc, f"{m} STO-3G")
    try:
        d = np.load("results/cisplatin_cas10_10.npz")
        allok &= check(np.array(d["orbsym"]), 5, "Cisplatin CAS(10e,10o) C2v")
    except FileNotFoundError:
        pass
    print("table of N_UCCSD(n_o, n_v):")
    for k in (2, 3, 4, 5, 6, 8, 10):
        print(f"  CAS({2*k}e,{2*k}o): n_o=n_v={k}: N_UCCSD={n_uccsd(k, k)}, N_spin={n_spin(k, k)}")
    sys.exit(0 if allok else 1)
