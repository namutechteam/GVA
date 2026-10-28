"""Main-text section on crystals (Results, Table 3, Figs 9-10) and its Methods.

Imported by ms_build.py. Every number is read from results/ via facts.periodic().
"""
import numpy as np

from facts import periodic
from refs import cite

P = periodic()
S = P["per"]["systems"]
stat, tag, sp = P["stat"], P["tag"], P["scal"]
ONE = {x["n_k"]: x for x in sp["one_d"] if x["o"] == 1}
THREE = {x["mesh"][0]: x for x in sp["three_d"] if x["o"] == 2}
MC = {x["mesh"][0]: x for x in sp["molecular_crystal"]}
T5 = [tag(0.74, 1.50, 5), tag(1.0, 1.0, 5), tag(1.5, 1.5, 5), tag(2.0, 2.0, 5)]
SYS = P["sys"]
MEAS = {(x["system"], x["basis"]): x for x in (P["measurement_periodic"] or [])}


def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}×10^{{{int(e)}}}"


def mh(x, d=3):
    return f"{x:.{d}f}"


def results(ms, F):
    cisc = F["cis"]["counts"]
    g5 = SYS[T5[0]]["counts"]
    ms.heading("Translational symmetry in crystals", 2)
    ms.para(
        "In a crystal the same construction uses the translation group. A mesh of "
        "*N*_{*k*} *k*-points describes a Born–von Kármán supercell whose translation "
        "group has order *N*_{*k*}, and in a basis of Bloch orbitals an excitation can "
        "connect the reference to the ground state only if it conserves crystal momentum: "
        "**k**_{*a*} + **k**_{*b*} − **k**_{*i*} − **k**_{*j*} ≡ 0 for a double and "
        "**k**_{*a*} = **k**_{*i*} for a single excitation. The rule is exact in the same "
        "sense as point-group screening (Supplementary Note 11), and momentum-conserving "
        f"unitary coupled cluster has been used for periodic systems before{cite(50, 51)}. "
        "Bloch orbitals are complex, so each kept excitation carries a complex amplitude "
        "*t* and enters as exp(*t*T̂ − *t*^{*}T̂^{†}); the phase of *t* is a single-qubit "
        "rotation and costs no CNOT gates. Two tyings reduce the parameters further. "
        "Spin-flip partners share *t*, as for molecules, and the time-reversal partner of "
        "an excitation, obtained by replacing every **k** by −**k**, carries ±*t*^{*} when the "
        "orbitals at **k** and −**k** are chosen complex conjugate, which is possible for any "
        "spin-free Hamiltonian. With both tyings each spin-flip group carries one real "
        "parameter. We call the result GVA-k; its counts follow closed forms that extend "
        "equations (2)–(4) to any abelian group (Methods).",
        align="justify")
    cis10 = [r for r in F["scal"] if r["n_active_orbitals"] == 20][0]
    ms.para(
        "Measured against UCCSD applied to the supercell in a real orbital basis, the "
        "parameter count of GVA-k approaches 1/(2*N*_{*k*}) and the CNOT count 1/*N*_{*k*} "
        "(Fig. 9a,b). The molecular results follow the same two laws with *N*_{*k*} "
        "replaced by the order *h* of the point group; for cisplatin active spaces in "
        f"*C*_{{2v}} (*h* = 4) the ratios reach {cis10['GVA_params'] / cis10['UCCSD_params']:.2f} "
        f"and {cis10['GVA_cnot_1L'] / cis10['UCCSD_cnot_1L']:.2f} at 20 active orbitals. The "
        "difference lies in *h*. No abelian point group has an order above eight, whereas "
        "the translation group grows with the *k*-point mesh. On a 4×4×4 mesh with two "
        f"occupied and two virtual bands per cell ({THREE[4]['n_qubits']} qubits), one UCCSD "
        f"layer has {sci(THREE[4]['UCCSD_params'])} parameters and "
        f"{sci(THREE[4]['UCCSD_cnot_1L'])} CNOT gates, and GVA-k has "
        f"{sci(THREE[4]['GVA_k_params'])} and {sci(THREE[4]['k_cnot_1L'])} (Fig. 9c). "
        "Explicit enumeration of the pool and of the tied orbits reproduced the closed "
        "forms on 12 meshes in one, two and three dimensions (Supplementary Note 13).",
        align="justify")
    ms.figure("figures/Fig9.png", "Fig. 9 | Savings from symmetry against the order of the "
              "symmetry group.",
              "**a**, Parameters per layer of GVA divided by those of UCCSD. Molecules: "
              "point groups of order *h* = 4 (*C*_{2v}) or 8 (*D*_{2h}); cisplatin active "
              "spaces with 10–20 orbitals. Crystals: GVA-k against UCCSD on the supercell, "
              "with *h* = *N*_{*k*}, for one-dimensional meshes (one occupied and one virtual "
              "band per cell) and cubic meshes (two and two, or eight and eight bands). The "
              "small molecules lie above the line because their few orbitals occupy only two "
              "or three irreps. Dashed line: 1/(2*h*). **b**, The same for CNOT gates; dashed "
              "line: 1/*h*. **c**, Parameters (solid) and CNOT gates (dashed) per layer for "
              "cubic meshes from 1^{3} to 5^{3} with two occupied and two virtual bands per "
              "cell.", width_mm=170)

    c0 = [SYS[t]["c0_sq"] for t in T5]
    d1 = P["gva_minus_uccsd_1L_mHa"]
    ms.para(
        "We tested GVA-k on chains of H_{2} units in the STO-3G basis, with one unit per "
        "cell and *N*_{*k*} = 2–5 *k*-points along the chain (8–20 qubits). A dimerized "
        "chain, with H–H distances of 0.74 Å inside a unit and 1.50 Å between units, is a "
        "minimal molecular crystal. Uniform chains with spacings of 1.0, 1.5 and 2.0 Å "
        f"take the hydrogen chain{cite(52)} from weak to strong correlation: at "
        "*N*_{*k*} = 5 the weight of the reference determinant in the FCI state is "
        f"{c0[0]:.2f} for the dimerized chain and {c0[1]:.2f}, {c0[2]:.2f} and {c0[3]:.2f} "
        "for the uniform chains. The *k*-point integrals reproduce an independent "
        "calculation on an explicit supercell; the FCI energies of the two routes agree to "
        f"better than {sci(P['max_fci_supercell_diff'], 0)} Ha (Supplementary Note 12). At "
        f"20 qubits GVA-k uses {g5['GVA-k']['n_params']} parameters and "
        f"{g5['GVA-k']['cnot_1L']:,} CNOT gates per layer against "
        f"{g5['UCCSD']['n_params']} and {g5['UCCSD']['cnot_1L']:,} for supercell UCCSD, "
        "which is the size of the cisplatin UCCSD circuit; GVA needed "
        f"{cisc['GVA']['n_params']} parameters and {cisc['GVA']['cnot_1L']:,} CNOT gates "
        "there.",
        align="justify")
    _bench_paragraph(ms)
    _table3(ms)
    ms.figure("figures/Fig10.png", "Fig. 10 | VQE on periodic hydrogen chains.",
              "**a**, Error relative to FCI for the dimerized chain (*r* = 0.74 Å, *d* = "
              "1.50 Å) against the number of *k*-points. Single-layer ansätze and ADAPT-VQE: "
              "Hartree–Fock start. Two layers: median (symbol) and range (bar) over three "
              "starts (two-layer UCCSD for the 20-qubit 1.0 and 2.0 Å chains: one start). "
              "Lines: CCSD and CCSD(T), plotted as absolute values (both lie below FCI at 1.5 "
              "and 2.0 Å). Dotted: 1.6 mHa. **b**, Uniform chains at three spacings with "
              "*N*_{*k*} = 5 (solid) and 3 (faint). *c*_{0}^{2}: reference weight at "
              "*N*_{*k*} = 5. **c**, Error against CNOT gates for the four 20-qubit chains "
              "(two layers: best start; shading from dark to light: dimerized chain, 1.0, 1.5 "
              "and 2.0 Å). UCCSD acts on the real supercell basis; GVA-k and ADAPT-VQE on the "
              "Bloch "              "basis.", width_mm=170)


def _bench_paragraph(ms):
    tD, t10, t15, t20 = T5
    sD, s10, s15, s20 = (SYS[t] for t in T5)
    ccsd = lambda s: (s["E_CCSD"] - s["E_FCI"]) * 1e3
    ccsdt = lambda s: (s["E_CCSD(T)"] - s["E_FCI"]) * 1e3
    u2D, g2D = stat(tD, "UCCSD", 2), stat(tD, "GVA-k", 2)
    u210, g210 = stat(t10, "UCCSD", 2), stat(t10, "GVA-k", 2)
    adD = stat(tD, "ADAPT-VQE", None)
    nweak = len(P["weak"])
    ms.para(
        f"In the {['', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten'][nweak]} weakly correlated systems (the dimerized chain at every *N*_{{*k*}} "
        f"and the uniform 1.0 Å chain, *c*_{{0}}^{{2}} ≥ {P['min_c0_weak']:.2f}), single-layer "
        "GVA-k started from the Hartree–Fock state reproduced single-layer supercell UCCSD "
        f"to within {P['d1_weak_uHa']:.1f} µHa, and two layers of either ansatz reached "
        "chemical accuracy. At *N*_{*k*} = 5 the best two-layer errors were "
        f"{mh(u2D['min'])} mHa for UCCSD and {mh(g2D['min'])} mHa for GVA-k on the dimerized "
        f"chain, and {mh(u210['min'], 2)} and {mh(g210['min'], 2)} mHa on the 1.0 Å chain, "
        f"where CCSD is already off by {mh(ccsd(s10), 2)} mHa (Fig. 10a,b and Table 3). "
        "ADAPT-VQE with the GVA-k orbits as pool stopped at "
        f"{mh(adD['hf'])} mHa on the dimerized chain with {adD['params']} real parameters "
        f"and {adD['cnot']:,} CNOT gates, about the cost of one GVA-k layer. The Bloch basis "
        "is more expensive to measure than the real one: for the 20-qubit dimerized chain "
        "general-commutativity grouping needs "
        f"{MEAS[(tD, 'Bloch')]['GC']['shots'] / MEAS[(tD, 'real')]['GC']['shots']:.1f} times "
        "as many shots per energy evaluation (Supplementary Note 15), against a factor of "
        "five saved in CNOT gates.",
        align="justify")
    u115, g115 = stat(t15, "UCCSD", 1)["hf"], stat(t15, "GVA-k", 1)["hf"]
    u120, g120 = stat(t20, "UCCSD", 1)["hf"], stat(t20, "GVA-k", 1)["hf"]
    u215, g215 = stat(t15, "UCCSD", 2), stat(t15, "GVA-k", 2)
    u220, g220 = stat(t20, "UCCSD", 2), stat(t20, "GVA-k", 2)
    tr20, sp20 = P["tr_cost_mHa"][t20], P["spin_cost_mHa"][t20]
    sp_all = max(abs(v) for v in P["spin_cost_mHa"].values())
    u1r = [x for x in stat(t20, "UCCSD", 1)["runs"] if x["S2"] > 0.5]
    ms.para(
        "Stretching the uniform chain makes the problem strongly correlated: at "
        f"*N*_{{*k*}} = 5 the reference weight falls to {s15['c0_sq']:.2f} at 1.5 Å and "
        f"{s20['c0_sq']:.2f} at 2.0 Å, and CCSD and CCSD(T) become non-variational "
        f"({mh(ccsd(s15), 0)} and {mh(ccsdt(s15), 0)} mHa at 1.5 Å). At *N*_{{*k*}} = 5 no "
        "ansatz reaches chemical accuracy with one or two layers. Time-reversal tying is then no "
        "longer free. At one layer GVA-k lies "
        f"{mh(g115 - u115, 2)} mHa (1.5 Å) and {mh(g120 - u120, 1)} mHa (2.0 Å) above "
        "supercell UCCSD. The cost comes from the time-reversal relation: "
        f"at 2.0 Å adding it to the untied ansatz raises the energy by {mh(tr20, 1)} mHa, "
        "whereas spin-flip tying alone changes it by at most "
        f"{mh(sp_all, 1)} mHa in any system (Supplementary Fig. 2a). With a second layer "
        "the gap shrinks or reverses: the best two-layer errors were "
        f"{mh(u215['min'], 1)} (UCCSD) and {mh(g215['min'], 1)} mHa (GVA-k) at 1.5 Å and "
        f"{mh(u220['min'], 1)} (UCCSD, Hartree–Fock start only) and {mh(g220['min'], 1)} mHa "
        "(GVA-k) at 2.0 Å, with five times fewer CNOT gates for GVA-k. At *N*_{*k*} = 3 the "
        "comparison goes either way ("
        f"{mh(stat(tag(1.5, 1.5, 3), 'GVA-k', 2)['min'], 2)} against "
        f"{mh(stat(tag(1.5, 1.5, 3), 'UCCSD', 2)['min'], 2)} mHa at 1.5 Å, "
        f"{mh(stat(tag(2.0, 2.0, 3), 'GVA-k', 2)['min'], 1)} against "
        f"{mh(stat(tag(2.0, 2.0, 3), 'UCCSD', 2)['min'], 1)} mHa at 2.0 Å, where the UCCSD "
        f"states have ⟨Ŝ^{{2}}⟩ up to {stat(tag(2.0, 2.0, 3), 'UCCSD', 2)['S2max']:.2f}). "
        "In the same regime, "
        "randomly started supercell UCCSD and the untied momentum-screened ansatz converged "
        f"to states with ⟨Ŝ^{{2}}⟩ up to {P['S2_max']['UCCSD']:.1f}, lower in energy than the "
        "spin-pure solution"
        + (f" ({mh(u1r[0]['dE'] * 1e3, 1)} against {mh(u120, 1)} mHa for single-layer "
           "UCCSD at 2.0 Å)" if u1r else "")
        + f"; the largest ⟨Ŝ^{{2}}⟩ among GVA-k runs was {P['S2_max']['GVA-k']:.2f}.",
        align="justify")


def _table3(ms):
    rows = []
    for t in T5:
        s = SYS[t]
        u1, g1 = stat(t, "UCCSD", 1), stat(t, "GVA-k", 1)
        u2, g2 = stat(t, "UCCSD", 2), stat(t, "GVA-k", 2)
        ad = stat(t, "ADAPT-VQE", None)
        rows.append([f"{s['r']:.2f}/{s['d']:.2f}", f"{s['c0_sq']:.3f}",
                     f"{s['E_corr'] * 1e3:.1f}",
                     f"{(s['E_CCSD'] - s['E_FCI']) * 1e3:.3f}" + ("" if s["CCSD_converged"] else "*"),
                     f"{(s['E_CCSD(T)'] - s['E_FCI']) * 1e3:.3f}",
                     f"{u1['hf']:.3f}", f"{g1['hf']:.3f}", f"{u2['min']:.3f}", f"{g2['min']:.3f}",
                     f"{ad['hf']:.3f}" + ("†" if ad["runs"][0].get("stop_reason") == "time limit" else "")])
    g5 = SYS[T5[0]]["counts"]
    ms.table("**Table 3 | Periodic H_{2} chains with *N*_{*k*} = 5 (20 qubits).** *r*/*d*: "
             "distances inside and between H_{2} units (Å). *c*_{0}^{2}: reference weight in "
             "the FCI state. *E*_{corr}: FCI correlation energy of the five-cell supercell "
             "(mHa). Other columns: errors relative to FCI (mHa); 1L from the Hartree–Fock "
             "start, 2L best of three starts (UCCSD at 1.0 and 2.0 Å: Hartree–Fock start "
             "only). Per layer, supercell UCCSD has "
             f"{g5['UCCSD']['n_params']} parameters and {g5['UCCSD']['cnot_1L']:,} CNOT gates "
             f"and GVA-k {g5['GVA-k']['n_params']} and {g5['GVA-k']['cnot_1L']:,}."
             + (" *CCSD not converged." if any(not SYS[t]["CCSD_converged"] for t in T5) else "")
             + " †ADAPT-VQE stopped at the 30-min wall-time limit.",
             ["*r*/*d*", "*c*_{0}^{2}", "*E*_{corr}", "CCSD", "CCSD(T)", "UCCSD 1L",
              "GVA-k 1L", "UCCSD 2L", "GVA-k 2L", "ADAPT"], rows, size=7.5,
             col_widths=[17, 12, 15, 19, 18, 16, 16, 16, 16, 17])


def methods(ms):
    ms.heading("Periodic systems", 2)
    ms.para(
        "Chains of H_{2} units were built with two H atoms per cell on the *x* axis (H–H "
        "distance *r* inside a unit, *d* between units, lattice constant *r* + *d*), 10 Å "
        "between parallel chains, the STO-3G basis and Γ-centred meshes of *N*_{*k*} × 1 × 1 "
        f"*k*-points. Restricted Hartree–Fock with *k*-points{cite(25, 53)} used Gaussian "
        f"density fitting{cite(54)} without a correction for the exchange divergence, so "
        "that the integrals, the KRHF energy and the supercell energy belong to one "
        "Hamiltonian. Orbitals at **k** and −**k** were made complex conjugate and those at "
        "time-reversal-invariant points real; the supercell integrals were assembled from "
        "momentum-conserving blocks of the *k*-point integrals divided by *N*_{*k*}. The "
        "real supercell basis used for UCCSD and for the CISD, CCSD and CCSD(T) references "
        "was obtained by the unitary φ_{c} = (φ_{**k**} + φ_{−**k**})/√2, φ_{s} = (φ_{**k**} − "
        "φ_{−**k**})/(*i*√2), which leaves the Hartree–Fock determinant unchanged. Uniform "
        "chains were run only with odd *N*_{*k*}, for which the closed-shell reference is "
        "unique. The complex Hamiltonian was applied with the PySCF kernel for complex "
        "integrals; each complex excitation acts on a pair of determinants as a rotation "
        "by |*t*| with phase arg *t*, and gradients with respect to Re *t* and Im *t* "
        "follow from the same adjoint sweep as for real amplitudes (Supplementary Note 11). "
        "Optimizer settings were those of the molecular runs; random starts drew all real "
        "parameters from N(0, 0.05^{2}) with NumPy seeds 3000 + *k*. Each system was run "
        "with supercell UCCSD (one and two layers, three starts), the untied "
        "momentum-screened ansatz (one layer, three starts), each tying alone (one layer, "
        "Hartree–Fock start), GVA-k (one and two layers, three starts) and ADAPT-VQE over "
        "the GVA-k orbits. Two exceptions were made for cost: for the 20-qubit uniform "
        "chains at 1.0 and 2.0 Å two-layer UCCSD was run from the Hartree–Fock start only, "
        "and ADAPT-VQE stopped adding operators after 30 min of wall time, which ended the "
        "runs for the 20-qubit chains at 1.5 and 2.0 Å.",
        align="justify")
    ms.para(
        "For an abelian group *G* with orbital labels *g*_{*p*}, let *X*(*q*) be the number "
        "of occupied–virtual pairs with *g*_{*a*}*g*_{*i*}^{−1} = *q*, *P*_{o}(*Q*) and "
        "*P*_{v}(*Q*) the numbers of occupied and virtual pairs of distinct orbitals with "
        "product *Q*, and *n*_{self} the number of pairs (*i*, *a*) with (*g*_{*a*}"
        "*g*_{*i*}^{−1})^{2} = 1. Then *S* = *X*(1), *D*_{ss} = Σ_{*Q*} *P*_{o}(*Q*)"
        "*P*_{v}(*Q*), *D*_{αβ} = Σ_{*q*} *X*(*q*)*X*(*q*^{−1}) and",
        align="justify")
    ms.equation("*N*_{GVA-k} = *S* + *D*_{ss} + (*D*_{αβ} + *n*_{self})/2,", 6)
    ms.para(
        "which reduces to equation (4) for point groups, where every element is its own "
        "inverse and *n*_{self} = *n*_{o}*n*_{v}. For the translation group *q*^{−1} "
        "corresponds to −**q**. Without tying, the momentum-screened ansatz has 2(2*S* + "
        "2*D*_{ss} + *D*_{αβ}) real parameters; time-reversal tying alone gives 2*S* + "
        "2*D*_{ss} + *D*_{αβ}, and spin-flip tying alone twice equation (6). For an "
        "insulator with *o* occupied and *v* virtual bands per cell, *D*_{αβ} = "
        "*N*_{*k*}^{3}*o*^{2}*v*^{2} and *N*_{GVA-k} ≈ ¾*N*_{*k*}^{3}*o*^{2}*v*^{2}, against "
        "1.5*N*_{*k*}^{4}*o*^{2}*v*^{2} for UCCSD on the supercell. CNOT counts were summed "
        "with equation (5) over the momentum-conserving generators, grouping occupied and "
        "virtual pairs by total momentum.",
        align="justify")
