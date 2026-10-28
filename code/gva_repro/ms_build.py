"""Build the main manuscript (docx). All numbers come from facts.compute()."""
import json

import numpy as np

from docbuilder import Manuscript
from facts import compute, HA2KCAL, orbital_sensitivity
from refs import REFS, cite
import ms_periodic as MP

F = compute()
cis = F["cis"]
CM = F["cis_med"]
red = F["red"]
two = F["two"]
meas = {m["system"]: m for m in F["meas"]}
pes = F["pes"]
noise = F["noise"]


def z3(x):
    s = f"{x:.3f}"
    return "0.000" if s in ("-0.000", "0.000") else s


def f2(x):
    return f"{x:.2f}"


def f3(x):
    return f"{x:.3f}"


def pct(a, b):
    return f"{100 * (1 - b / a):.0f}"


def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}×10^{{{int(e)}}}"


def best(rec, kind, key):
    return min(rec[kind].values(), key=lambda v: v[key])[key]


def nb(r, key):
    return [x for x in noise if x["system"].startswith(r[0]) and x["ansatz"] == r[1]
            and x["layers"] == r[2]][0][key]


ms = Manuscript("Symmetry-screened unitary coupled-cluster circuits", author="")
import numbering
ms.remap = numbering.remap
MOVED = {}          # sections placed in the SI; called by si_build with (doc, cite)
import re as _re
_jw = open("results/logs/check_jw_cisplatin.log").read()
_v = max(float(x) for x in _re.findall(r"= ([0-9.]+e[+-][0-9]+) Ha", _jw))
_e = int(np.floor(np.log10(_v)))
JWC = f"{int(np.ceil(_v / 10 ** _e))}×10^{{{_e}}}"

# ------------------------------------------------------------------ front
ms.title("Symmetry-screened unitary coupled-cluster circuits for molecules and crystals")
ms.para("[Author names]^{1}", space_after=2)
ms.para("^{1}[Affiliation, address]", size=10, space_after=2)
ms.para("Correspondence: [e-mail address]", size=10, space_after=10)

ms.heading("Abstract")
gva_c = cis["counts"]
g2r = F["cis_runs"][("GVA", 2)]
g5 = MP.SYS[MP.T5[0]]["counts"]
ab = (
    "Unitary coupled-cluster singles and doubles (UCCSD) circuits for variational quantum "
    "eigensolvers grow as the fourth power of the active-space size. Our geometric "
    "variational ansatz (GVA) removes generators forbidden by an abelian symmetry group and "
    "ties symmetry-related generators to one parameter, which reduces parameters to about "
    "1/(2*h*) and CNOT gates to 1/*h* of UCCSD for a group of order *h*. Molecular point "
    "groups have *h* ≤ 8: in a 20-qubit cisplatin active space GVA cut each layer from "
    f"{gva_c['UCCSD']['n_params']} to {gva_c['GVA']['n_params']} parameters and from "
    f"{gva_c['UCCSD']['cnot_1L']:,} to {gva_c['GVA']['cnot_1L']:,} CNOT gates without loss of "
    "accuracy. In crystals *h* is the number of *k*-points. With crystal-momentum screening "
    "and time-reversal tying, 20-qubit hydrogen chains needed "
    f"{g5['GVA-k']['n_params']} instead of {g5['UCCSD']['n_params']} parameters and "
    f"{g5['GVA-k']['cnot_1L']:,} instead of {g5['UCCSD']['cnot_1L']:,} CNOT gates per layer. "
    f"One layer matched supercell UCCSD to within {MP.P['d1_weak_uHa']:.1f} µHa at weak "
    "correlation, whereas at strong correlation time-reversal tying cost accuracy. Spin-flip "
    "tying excluded spin-contaminated minima."
)
assert len(ab.replace("_{", "").replace("*", "").split()) <= 150, len(ab.split())
ms.para(ab, align="justify")

# ------------------------------------------------------------------ main
ms.heading("Introduction")
ms.para(
    "Quantum computers may one day prepare molecular ground states whose classical "
    f"description grows exponentially with system size{cite(3, 4)}. On current hardware the "
    "most studied route is the variational quantum eigensolver (VQE), in which a "
    "parameterized circuit prepares a trial state and a classical optimizer lowers its "
    f"measured energy{cite(1, 2)}. The choice of circuit family, or ansatz, fixes how many "
    "parameters must be optimized, how many two-qubit gates each energy evaluation costs, "
    f"and whether the optimizer can reach the ground state at all{cite(7, 39)}.",
    align="justify")
no, nv = 5, 5
ms.para(
    "Unitary coupled cluster with singles and doubles (UCCSD) is the usual chemically "
    f"motivated ansatz{cite(5)}. Its generators are fermionic excitations, so the trial "
    "state keeps the particle number of the reference, and products of such exponentials "
    f"can represent the exact state when enough of them are applied{cite(33)}. The "
    "drawback is size. For *n*_{o} occupied and *n*_{v} virtual spatial orbitals, one "
    "spin-conserving UCCSD layer has 2*n*_{o}*n*_{v} + 2C(*n*_{o},2)C(*n*_{v},2) + "
    "(*n*_{o}*n*_{v})^{2} parameters, about 1.5(*n*_{o}*n*_{v})^{2} for large spaces, and "
    f"its Jordan–Wigner circuit{cite(40)} for a 20-qubit active space already contains "
    f"{gva_c['UCCSD']['cnot_1L']:,} CNOT gates. Hardware-efficient circuits{cite(8)} are "
    "far shallower but break particle-number and spin conservation and are prone to "
    f"vanishing gradients{cite(7)}. Particle-conserving gate fabrics{cite(9, 10)} and "
    f"adaptive operator selection{cite(6)} lie between these two extremes.",
    align="justify")
ms.para(
    "Symmetry reduces this cost. An excitation that changes the irreducible "
    "representation of the reference cannot contribute to a totally symmetric ground "
    "state, and the α and β amplitudes of a closed-shell singlet are linked by spin-flip "
    "symmetry. Both facts are textbook coupled-cluster theory and have been used in "
    f"quantum algorithms to taper qubits{cite(12, 13)}, to prune excitation "
    f"lists{cite(11)} and to build spin-adapted operator pools{cite(6, 14)}. In a crystal "
    "the translation group plays the role of the point group, and crystal momentum "
    f"restricts the excitations in the same way{cite(50, 51)}. Less is known about how much "
    "each constraint saves on its own, what tying costs once the product of exponentials "
    f"is Trotterized{cite(34)}, how the saving depends on the size of the symmetry group, "
    "and how the resulting circuits compare with the error rates of present devices.",
    align="justify")
ms.para(
    "Here we combine symmetry screening and parameter tying in a geometric variational "
    "ansatz (GVA) and test it in exact state-vector simulations, from H_{2} to a 20-qubit "
    "CAS(10e,10o) model of the anticancer drug "
    f"cisplatin{cite(30)} and to periodic hydrogen chains with up to 20 qubits, for which "
    "crystal-momentum screening is combined with spin-flip and time-reversal tying. The "
    "savings follow from the order of the symmetry group, which is at most eight for a "
    "molecule and equals the number of *k*-points for a crystal. Every VQE energy is "
    "compared with full configuration "
    "interaction (FCI) or complete-active-space configuration interaction (CASCI) in the "
    "same orbital space, every parameter and gate count follows from closed-form "
    "expressions that we check against explicit enumeration, and all scripts and source "
    "data are released with the paper. We also estimate the two-qubit error rates such "
    "circuits would need under a simple depolarizing model.",
    align="justify")

# ------------------------------------------------------------------ results
ms.heading("Results")
ms.heading("Construction of the ansatz", 2)
ms.para(
    "GVA starts from the spin-conserving UCCSD generators Ĝ_{k} = T̂_{k} − T̂_{k}^{†}, "
    "where T̂_{k} is a single excitation *a*^{†}_{*a*}*a*_{*i*} or a double excitation "
    "*a*^{†}_{*a*}*a*^{†}_{*b*}*a*_{*j*}*a*_{*i*} from occupied (*i*, *j*) to virtual "
    "(*a*, *b*) spin orbitals of the reference determinant |Φ_{0}⟩ (Fig. 1a). Two "
    "filters are then applied. For an abelian point group with orbital irreducible "
    "representations Γ_{*p*}, a double excitation is kept only if Γ_{*a*}⊗Γ_{*b*}⊗"
    "Γ_{*i*}⊗Γ_{*j*} is totally symmetric (Γ_{*a*}⊗Γ_{*i*} for singles). Each surviving "
    "generator then shares one parameter with its spin-flipped partner: αα with ββ, "
    "*a*_{α}←*i*_{α} with *a*_{β}←*i*_{β}, and *a*_{α}*b*_{β}←*i*_{α}*j*_{β} with "
    "*b*_{α}*a*_{β}←*j*_{α}*i*_{β}. The trial state is",
    align="justify")
ms.equation("|Ψ(**θ**)⟩ = ∏_{*l*=1}^{*L*} ∏_{*k*} exp(θ_{*lk*} Ĝ_{*k*}) |Φ_{0}⟩,", 1)
ms.para(
    "with *L* layers of independent parameters and tied partners applied next to each "
    "other. Because each Ĝ_{*k*} satisfies Ĝ_{*k*}^{3} = −Ĝ_{*k*}, every exponential "
    "is a set of Givens rotations between pairs of determinants, which makes exact "
    "simulation in the fixed-particle-number space straightforward (Methods).",
    align="justify")
ms.para(
    "Screening removes generators without removing reachable states. If |Ψ⟩ is totally "
    "symmetric and Ĝ transforms as a non-trivial irrep, ⟨Ψ|[Ĥ, Ĝ]|Ψ⟩ changes sign under "
    "some group operation and therefore vanishes (Supplementary Note 2). A forbidden "
    "generator thus has zero gradient along any trajectory generated by allowed "
    "generators, and an optimizer started at **θ** = 0 never moves it. Tying is an "
    "approximation of a different kind. A tied pair maps onto itself under the α↔β "
    "exchange, so the ansatz keeps the exchange symmetry *c*(*I*_{α}, *I*_{β}) = "
    "*c*(*I*_{β}, *I*_{α}) of a closed-shell singlet, apart from ordering effects between "
    "non-commuting αβ partners. That symmetry excludes triplet components, but a tied "
    "pair is not a complete spin-adapted double excitation, so the energy cost of tying "
    "has to be measured. For the "
    "cisplatin active space in *C*_{2v}, screening removes "
    f"{gva_c['UCCSD']['n_params'] - gva_c['PG only']['n_params']} of the "
    f"{gva_c['UCCSD']['n_params']} UCCSD generators and tying reduces the remaining "
    f"{gva_c['PG only']['n_params']} to {gva_c['GVA']['n_params']} parameters; the CNOT "
    f"count falls from {gva_c['UCCSD']['cnot_1L']:,} to {gva_c['GVA']['cnot_1L']:,} per "
    "layer (Fig. 1b). Closed-form expressions for all four counts (UCCSD, screening only, "
    "tying only and GVA) are given in Methods and reproduce the enumerated values for "
    "every system in this work (Supplementary Note 4).",
    align="justify")
ms.figure("figures/Fig1.png", "Fig. 1 | Construction of GVA.",
          "**a**, Workflow. Molecular orbitals and their irreducible representations "
          "define the UCCSD generator pool; generators whose irrep product is not "
          "totally symmetric are removed, spin-flip partners share one parameter, and "
          "the remaining exponentials are applied in *L* layers. **b**, Parameters and "
          "CNOT gates per layer for cisplatin CAS(10e,10o) in *C*_{2v} after each step. "
          "Spin tying leaves the gate count unchanged because both members of a pair "
          "are still applied.", width_mm=165)

# ---------------- small molecules
ms.heading("Accuracy and optimization on small molecules", 2)
sd = F["sysd"]
ms.para(
    "We first treated H_{2}, linear H_{4}, LiH and BeH_{2} in the STO-3G basis at "
    "reference bond lengths and with all bonds stretched by factors of 1.5 and 2.0 "
    "(Table 1, Fig. 2). Stretching increases static correlation: the weight of the "
    "Hartree–Fock determinant in the FCI wavefunction falls from "
    f"{sd[('H4', 1.0)]['c0_sq']:.2f} to {sd[('H4', 2.0)]['c0_sq']:.2f} for H_{{4}} and "
    f"from {sd[('BeH2', 1.0)]['c0_sq']:.2f} to {sd[('BeH2', 2.0)]['c0_sq']:.2f} for "
    "BeH_{2}, and the CCSD error reaches "
    f"{(sd[('H4', 2.0)]['E_CCSD'] - sd[('H4', 2.0)]['E_FCI']) * 1e3:.1f} mHa and "
    f"{(sd[('BeH2', 2.0)]['E_CCSD'] - sd[('BeH2', 2.0)]['E_FCI']) * 1e3:.1f} mHa, "
    "respectively. Per layer, GVA used "
    f"{red['H2'][1]}, {red['H4'][1]}, {red['LiH'][1]} and {red['BeH2'][1]} parameters "
    f"against {red['H2'][0]}, {red['H4'][0]}, {red['LiH'][0]} and {red['BeH2'][0]} for "
    f"UCCSD, and {pct(red['H2'][2], red['H2'][3])}–{pct(red['BeH2'][2], red['BeH2'][3])}% "
    "fewer CNOT gates.",
    align="justify")
d1 = sorted(F["gva_minus_uccsd_1L_seed0"], key=lambda x: -x[0])[0]
u_lih = [r for r in F["small_groups"][("LiH", 2.0, "UCCSD", 1)] if r["seed"] == 0][0]["dE"] * 1e3
g_lih = [r for r in F["small_groups"][("LiH", 2.0, "GVA", 1)] if r["seed"] == 0][0]["dE"] * 1e3
worst1 = max(abs(r["dE"]) for r in F["small_groups"][("BeH2", 2.0, "GVA", 1)]) * 1e3
max2_u = max(two[(m, s, "UCCSD")][2] for m in ("H2", "H4", "LiH", "BeH2") for s in (1.0, 1.5, 2.0))
max2_g = max(two[(m, s, "GVA")][2] for m in ("H2", "H4", "LiH", "BeH2") for s in (1.0, 1.5, 2.0))
assert all(two[("BeH2", s, "GVA")][0] < two[("BeH2", s, "UCCSD")][0] for s in (1.0, 1.5, 2.0))
assert all(r["nit"] >= 3000 for s in (1.0, 1.5, 2.0) for r in F["small_groups"][("BeH2", s, "UCCSD", 2)])
ms.para(
    "Started from the Hartree–Fock state, single-layer GVA reproduced the single-layer "
    f"UCCSD energy to within {F['max_gva_penalty_1L_mHa']:.3f} mHa at all twelve "
    f"geometries; the largest difference, for LiH ×2.0 ({g_lih:.3f} versus "
    f"{u_lih:.3f} mHa), comes from tying, since screening alone gives the UCCSD value "
    "(Supplementary Table S4). Both single-layer ansätze missed the 1.6 mHa threshold at "
    "the four most strongly correlated geometries (H_{4} and BeH_{2} at ×1.5 and ×2.0), "
    f"with errors up to {worst1:.1f} mHa. A second layer closed this gap for both: over "
    f"five starts per geometry the largest error was {max2_u:.3f} mHa for UCCSD and "
    f"{max2_g:.3f} mHa for GVA. For BeH_{{2}}, two-layer UCCSD stopped at the 3,000-"
    "iteration limit in all fifteen runs, whereas two-layer GVA converged within the "
    "limit in "
    f"{sum(1 for s in (1.0, 1.5, 2.0) for r in F['small_groups'][('BeH2', s, 'GVA', 2)] if r['nit'] < 3000)}"
    " of fifteen and reached lower median errors at every BeH_{2} geometry.",
    align="justify")
h4u = [r for r in F["small_groups"][("H4", 2.0, "UCCSD", 1)] if r["seed"] == 1][0]
h4u0 = [r for r in F["small_groups"][("H4", 2.0, "UCCSD", 1)] if r["seed"] == 0][0]
beu = [r for r in F["small_groups"][("BeH2", 2.0, "UCCSD", 1)] if r["seed"] == 1][0]
beu0 = [r for r in F["small_groups"][("BeH2", 2.0, "UCCSD", 1)] if r["seed"] == 0][0]
ncont = len(F["spin_contaminated"])
n_main_runs = sum(len(v) for k, v in F["small_groups"].items() if k[2] in ("UCCSD", "GVA"))
assert all(an == "UCCSD" and sd_ > 0 for _, _, an, _, sd_, _, _ in F["spin_contaminated"])
ms.para(
    "Random starting parameters revealed a difference that the Hartree–Fock start "
    "hides. At H_{4} ×2.0 and BeH_{2} ×2.0, all four random starts drove single-layer "
    f"UCCSD to energies {(h4u0['dE'] - h4u['dE']) * 1e3:.1f} and "
    f"{(beu0['dE'] - beu['dE']) * 1e3:.1f} mHa below its symmetric optimum, but the "
    f"states found had ⟨Ŝ^{{2}}⟩ = {h4u['S2']:.2f} and {beu['S2']:.2f}, that is, they "
    "mixed higher spin states into the singlet (Fig. 2a and Supplementary Fig. 1). The "
    "truncated ansatz lowers its energy by breaking spin symmetry, much as unrestricted "
    "Hartree–Fock does, and the lower energy belongs to a physically wrong state. In "
    "the stretched H_{4} chain the lowest triplet lies only "
    f"{F['h4_triplet_gap_mHa']:.1f} mHa above the singlet, which makes such mixing easy. "
    "The CI vectors of singlet and quintet states are symmetric under α↔β exchange and "
    "those of triplets antisymmetric, so the exchange symmetry kept by tying removes the "
    "triplet channel. GVA converged to the same spin-pure solution from every start "
    "(largest ⟨Ŝ^{2}⟩ "
    f"{sci(F['gva_max_S2'])}). In total {ncont} of the {n_main_runs} UCCSD and GVA runs ended "
    "with ⟨Ŝ^{2}⟩ > 0.01, all of them UCCSD from random starts.",
    align="justify")
ad_h415 = F["small_groups"][("H4", 1.5, "ADAPT-VQE", None)][0]
ad_h42 = F["small_groups"][("H4", 2.0, "ADAPT-VQE", None)][0]
ad_be2 = F["small_groups"][("BeH2", 2.0, "ADAPT-VQE", None)][0]
g2_be2 = [r for r in F["small_groups"][("BeH2", 2.0, "GVA", 2)] if r["seed"] == 0][0]
ms.para(
    f"ADAPT-VQE{cite(6)} drawing operators from the same screened and tied pool "
    f"reached {sci(abs(ad_h415['dE']) * 1e3)} mHa for H_{{4}} ×1.5 with {ad_h415['n_params']} "
    f"operators, but stopped at {ad_h42['dE'] * 1e3:.2f} mHa for H_{{4}} ×2.0 when the "
    "largest pool gradient fell below 10^{−3} Ha. For BeH_{2} ×2.0 it reached "
    f"{ad_be2['dE'] * 1e3:.3f} mHa with {ad_be2['n_params']} operators and "
    f"{ad_be2['cnot']:,} CNOT gates, against {g2_be2['dE'] * 1e3:.4f} mHa and "
    f"{g2_be2['cnot']:,} CNOT gates for two-layer GVA (Table 1).",
    align="justify")

# Table 1
rows = []
for m, lab in (("H2", "H_{2}"), ("H4", "H_{4}"), ("LiH", "LiH"), ("BeH2", "BeH_{2}")):
    for sc in (1.0, 1.5, 2.0):
        s = sd[(m, sc)]
        g = F["small_groups"]
        u1 = [r for r in g[(m, sc, "UCCSD", 1)] if r["seed"] == 0][0]["dE"] * 1e3
        g1 = [r for r in g[(m, sc, "GVA", 1)] if r["seed"] == 0][0]["dE"] * 1e3
        u2 = two[(m, sc, "UCCSD")][0]
        g2 = two[(m, sc, "GVA")][0]
        ad = g[(m, sc, "ADAPT-VQE", None)][0]

        def e(x):
            return "<10^{−6}" if abs(x) < 1e-6 else (f"{x:.3f}" if abs(x) >= 0.01 else sci(x))
        rows.append([f"{lab} ×{sc:.1f}" if sc != 1.0 else f"{lab} ×1.0",
                     f"{s['counts']['UCCSD']['n_params']} / {s['counts']['GVA']['n_params']}",
                     f"{s['c0_sq']:.2f}",
                     z3((s['E_CCSD'] - s['E_FCI']) * 1e3),
                     e(u1), e(g1), e(u2), e(g2), f"{e(ad['dE'] * 1e3)} ({ad['n_params']})"])
ms.table("**Table 1 | Small-molecule benchmark (STO-3G).** Errors relative to FCI in mHa. "
         "Single-layer (1L) values are for the Hartree–Fock start; two-layer (2L) values are "
         "medians over five starts. Parameters are per layer.",
         ["System", "Params UCCSD / GVA", "*c*_{0}^{2}", "CCSD", "UCCSD 1L", "GVA 1L",
          "UCCSD 2L", "GVA 2L", "ADAPT (ops)"], rows,
         footnote="*c*_{0}^{2}, weight of the Hartree–Fock determinant in the FCI state. "
                  "CCSD errors below zero are non-variational. Full per-run data, including "
                  "minima, maxima, iteration counts and ⟨Ŝ^{2}⟩, are in Supplementary "
                  "Table S4 and Supplementary Data 1.")
ms.figure("figures/Fig2.png", "Fig. 2 | Small-molecule benchmark.",
          "**a**, Energy error relative to FCI for UCCSD and GVA with one and two layers "
          "and for ADAPT-VQE with the GVA pool. Symbols are medians over five starts "
          "(Hartree–Fock start and four random starts); vertical bars span the minimum "
          "and maximum. Values below 10^{−6} mHa are drawn at that floor. The dotted line "
          "marks 1.6 mHa. Where UCCSD 1L error bars extend below GVA (H_{4} and BeH_{2} "
          "×2.0) the lower-energy UCCSD states are spin-contaminated (main text). "
          "**b**,**c**, Parameters (**b**) and CNOT gates (**c**) per layer at the "
          "reference geometries; percentages give the CNOT reduction. **d**, Energy "
          "error during L-BFGS optimization for BeH_{2} with bonds stretched ×2.0, "
          "Hartree–Fock start.")

# ---------------- measurement (Supplementary Note 7)
def _mv_measurement(ms, cite):
    mc = meas["Cisplatin CAS(10e,10o)"]
    ms.para(
        "Energy estimation adds a second cost that does not depend on the ansatz. We grouped "
        "the Jordan–Wigner Pauli terms of each Hamiltonian either by qubit-wise commutativity "
        "(QWC), which needs only single-qubit basis changes, or by general commutativity (GC), "
        f"which needs a Clifford circuit per group{cite(16, 17, 18, 21)}, using three "
        "colouring heuristics for each (Methods). For the cisplatin Hamiltonian, with "
        f"{mc['n_pauli_terms']:,} Pauli terms including the identity, the fewest groups found "
        f"were {best(mc, 'QWC', 'n_groups'):,} (QWC) and {best(mc, 'GC', 'n_groups')} (GC) "
        "(Fig. 3a). Group counts overstate the gain. Using the exact ground state to evaluate "
        "the variance of every group and allocating shots optimally"
        f"{cite(19, 20)}, a standard error of 1.6 mHa on one energy evaluation required "
        f"{sci(best(mc, 'QWC', 'shots'))} shots with QWC and "
        f"{sci(best(mc, 'GC', 'shots'))} with GC, a factor of "
        f"{best(mc, 'QWC', 'shots') / best(mc, 'GC', 'shots'):.1f}, whereas the group counts "
        f"differ by a factor of {best(mc, 'QWC', 'n_groups') / best(mc, 'GC', 'n_groups'):.0f} "
        "(Fig. 3b). For H_{2} the two schemes needed the same number of shots. Apart from "
        "H_{2}, where all partitions coincide, sorted insertion"
        f"{cite(19)} gave the lowest shot count for both partition types without ever giving "
        "the fewest groups (Supplementary Table S6).",
        align="justify")
    ms.figure("figures/Fig3.png", "Fig. 3 | Measurement groups and shot counts.",
              "**a**, Fewest measurement groups found for qubit-wise commuting (QWC) and "
              "general commuting (GC) partitions of the Jordan–Wigner Hamiltonians. "
              "**b**, Shots per energy evaluation for a standard error of 1.6 mHa under "
              "optimal allocation, evaluated with the exact ground state; labels give the "
              "GC/QWC ratio. The Clifford circuits required by GC groups are not included.",
              width_mm=140)



MOVED['measurement'] = _mv_measurement

# ---------------- cisplatin
ms.heading("Cisplatin CAS(10e,10o)", 2)
occ = cis["active_NO_occupations"]
nc = cis["n_core_orbitals"]
orb_irr = [o["irrep"] for o in json.load(open("results/cisplatin_orbitals.json"))["orbitals"]]


def irr_txt(irreps):
    from collections import Counter
    c = Counter(irreps)
    parts = []
    for g in sorted(c):
        lab = g[0].lower() + "_{" + g[1:] + "}"
        parts.append(f"{c[g]} {lab}" if c[g] > 1 else lab)
    return " + ".join(parts)
ms.para(
    f"Cisplatin, *cis*-[Pt(NH_{{3}})_{{2}}Cl_{{2}}], was modelled at an idealized "
    "*C*_{2v} geometry (Pt–Cl 2.32 Å, Pt–N 2.05 Å; Supplementary Table S1) with the "
    f"def2-SVP basis and a 60-electron effective core potential on Pt{cite(27, 28)}, "
    f"leaving {cis['n_electrons_total']} explicit electrons in {cis['n_ao']} basis "
    f"functions. The active space consists of symmetry-adapted MP2 natural "
    f"orbitals{cite(29)}: of the {cis['n_ao']} natural orbitals ranked by occupation, "
    f"numbers {nc + 1}–{nc + 5} (occupations "
    f"{min(o for o in occ if o > 1):.3f}–{max(o for o in occ if o > 1):.3f}; "
    f"{irr_txt(orb_irr[:5])}) and {nc + 6}–{nc + 10} (occupations "
    f"{min(o for o in occ if o < 1):.3f}–{max(o for o in occ if o < 1):.3f}; "
    f"{irr_txt(orb_irr[5:])}), that is, the ten whose occupations deviate most from "
    "2 or 0 (Fig. 4a). "
    "The strongly occupied set consists of ammine and chloride σ-donor combinations and "
    "a Pt d-dominated orbital; the weakly occupied set contains the Pt d–ligand σ^{*} "
    f"orbital (occupation {max(o for o in occ if o < 1):.3f}), "
    "Pt s- and p-type virtuals and an N–H antibonding combination (Supplementary "
    "Table S9). The "
    f"remaining {cis['n_core_orbitals']} doubly occupied orbitals were frozen. In this "
    f"space the correlation energy is {F['cis_corr_mHa']:.2f} mHa, the Hartree–Fock "
    f"determinant carries {cis['c0_sq'] * 100:.1f}% of the CASCI wavefunction, and "
    f"CISD, CCSD and CCSD(T) deviate from CASCI by {F['cis_E_CISD_err_mHa']:.3f}, "
    f"{F['cis_E_CCSD_err_mHa']:.3f} and {F['cis_E_CCSD(T)_err_mHa']:.3f} mHa (Table 2). "
    "At equilibrium the active space is therefore weakly correlated and classically easy; "
    "the benchmark tests circuit cost and optimization and does not address quantum "
    "advantage.",
    align="justify")
OS = {r["choice"][0]: r for r in orbital_sensitivity()}
assert abs(OS["A"]["E_CASCI"] - cis["E_CASCI"]) < 1e-8
ms.para(
    "The choice of active orbitals changes the problem substantially, so we repeated the "
    "equilibrium calculation with two other sets of ten orbitals (Supplementary Note 10, "
    "Supplementary Table S10). With the canonical RHF orbitals HOMO−4 to LUMO+4, the "
    f"correlation energy in the active space falls to {OS['B']['E_corr_active_mHa']:.2f} mHa: "
    "the five lowest canonical virtuals are dominated by H s, Pt s and Pt p functions "
    "and correlate the occupied orbitals only weakly. Every method then looks accurate; "
    f"single-layer GVA deviates from CASCI by {OS['B']['GVA_1L']['dE_mHa']:.3f} mHa and "
    f"CCSD by {OS['B']['CCSD_err_mHa']:.3f} mHa. With CASSCF(10e,10o) orbitals optimized from "
    f"the MP2 guess, the active-space correlation energy is "
    f"{OS['C']['E_corr_active_mHa']:.1f} mHa, the single-layer GVA error rises to "
    f"{OS['C']['GVA_1L']['dE_mHa']:.2f} mHa and the CCSD error to "
    f"{OS['C']['CCSD_err_mHa']:.2f} mHa, while two-layer GVA reaches "
    f"{OS['C']['GVA_2L']['dE_mHa']:.3f} mHa. GVA and UCCSD agreed to within "
    f"{max(abs(OS[k]['GVA_1L']['dE_mHa'] - OS[k]['UCCSD_1L']['dE_mHa']) for k in 'ABC') * 1e3:.1f} µHa "
    "with one layer for all three orbital sets. An energy error in millihartree is "
    "therefore informative only together with the orbital definition, and a fixed 1.6 mHa "
    "threshold is easily met in an active space that holds little correlation.",
    align="justify")


def cm(k):
    return CM[k]


same1 = abs(cm(('GVA', 1))[0] - cm(('UCCSD', 1))[0]) < 5e-5 and \
    cm(('GVA', 1))[2] - cm(('GVA', 1))[1] < 5e-5 and cm(('UCCSD', 1))[2] - cm(('UCCSD', 1))[1] < 5e-5
s1 = (f"the same energy, {cm(('GVA', 1))[0]:.4f} mHa above CASCI, from all three starts each"
      if same1 else
      f"{cm(('GVA', 1))[0]:.4f} and {cm(('UCCSD', 1))[0]:.4f} mHa above CASCI (medians over "
      "three starts)")
assert 5.5 < cm(('UCCSD', 1))[0] / cm(('GVA', 2))[0] < 6.5
ms.para(
    f"Single-layer GVA and single-layer UCCSD reached {s1} "
    "(Table 2, Fig. 4b). With two layers the median errors were "
    f"{cm(('GVA', 2))[0]:.4f} mHa for GVA ({cm(('GVA', 2))[3]} parameters, "
    f"{cm(('GVA', 2))[4]:,} CNOT gates) and {cm(('UCCSD', 2))[0]:.4f} mHa for UCCSD "
    f"({cm(('UCCSD', 2))[3]:,} parameters, {cm(('UCCSD', 2))[4]:,} CNOT gates), both below "
    f"the CCSD error of {F['cis_E_CCSD_err_mHa']:.3f} mHa. Three layers gave "
    f"{cm(('GVA', 3))[0]:.4f} and {cm(('UCCSD', 3))[0]:.4f} mHa. At equal depth GVA thus "
    f"stayed within {abs(cm(('GVA', 2))[0] - cm(('UCCSD', 2))[0]):.3f} mHa of UCCSD (two "
    f"layers, medians) and within {abs(cm(('GVA', 3))[0] - cm(('UCCSD', 3))[0]):.4f} mHa "
    "(three layers) while using about a fifth of the parameters and a third of the CNOT "
    f"gates. The two-layer GVA circuit has fewer CNOT gates ({cm(('GVA', 2))[4]:,}) than "
    f"single-layer UCCSD ({cm(('UCCSD', 1))[4]:,}) and an error six times smaller (Fig. 5). "
    "Running only the screening step or only the tying step gave "
    f"{cm(('PG only', 1))[0]:.4f} mHa with one layer in both cases, so neither constraint "
    "costs accuracy here.",
    align="justify")
ad = cis["runs"][[r["ansatz"] for r in cis["runs"]].index("ADAPT-VQE")]
stop_txt = {"max_ops": "at the limit of 120 operators", "gradient":
            "when the pool gradient fell below 10^{−3} Ha", "stall":
            "when five additions gained less than 10^{−7} Ha"}[ad["stop_reason"]]
g1 = CM[("GVA", 1)]
ms.para(
    f"ADAPT-VQE with the GVA pool stopped {stop_txt}, at {ad['dE'] * 1e3:.3f} mHa with "
    f"{ad['cnot']:,} CNOT gates (Fig. 4c), close to single-layer GVA ({g1[0]:.3f} mHa, "
    f"{g1[4]:,} CNOT gates) with {100 * (1 - ad['cnot'] / g1[4]):.0f}% fewer gates but "
    f"{ad['n_params']} sequential re-optimizations. Wall times on a two-core virtual machine were "
    f"{cm(('GVA', 2))[7] / 60:.0f} min for a two-layer GVA optimization and "
    f"{cm(('UCCSD', 2))[7] / 60:.0f} min for two-layer UCCSD (Supplementary Table S3); "
    "these times refer to the classical simulator.",
    align="justify")

# Table 2
rows = [["RHF reference", "–", "–", f"{(cis['E_ref'] - cis['E_CASCI']) * 1e3:.3f}", "–"],
        ["CISD", "–", "–", f"{F['cis_E_CISD_err_mHa']:.4f}", "–"],
        ["CCSD", "–", "–", f"{F['cis_E_CCSD_err_mHa']:.4f}", "–"],
        ["CCSD(T)", "–", "–", f"{F['cis_E_CCSD(T)_err_mHa']:.4f}", "–"]]
for an in ("UCCSD", "GVA"):
    for L in (1, 2, 3):
        v = CM[(an, L)]
        rng = f"{v[0]:.4f}" if v[5] == 1 else f"{v[0]:.4f} ({v[1]:.4f}–{v[2]:.4f})"
        rows.append([f"{an} {L}L", f"{v[3]:,}", f"{v[4]:,}", rng, f"{v[5]}"])
for an, lab in (("PG only", "Screening only 1L"), ("spin only", "Tying only 1L")):
    v = CM[(an, 1)]
    rows.append([lab, f"{v[3]:,}", f"{v[4]:,}", f"{v[0]:.4f}", "1"])
rows.append(["ADAPT-VQE (GVA pool)", f"{ad['n_params']}", f"{ad['cnot']:,}",
             f"{ad['dE'] * 1e3:.4f}", "1"])
ms.table("**Table 2 | Cisplatin CAS(10e,10o), *C*_{2v} geometry.** Energies relative to "
         f"CASCI (E_{{CASCI}} = {cis['E_CASCI']:.6f} Ha). For VQE rows the median over "
         "starts is given with the range in parentheses.",
         ["Method", "Parameters", "CNOT", "*E* − *E*_{CASCI} (mHa)", "Starts"], rows,
         footnote="All methods use the same frozen-core CAS(10e,10o) Hamiltonian. CNOT "
                  "counts use the ladder compilation described in Methods. Individual runs "
                  "are listed in Supplementary Table S3.")
ms.figure("figures/Fig4.png", "Fig. 4 | Cisplatin CAS(10e,10o).",
          "**a**, Deviation of the active MP2 natural-orbital occupations from 2 "
          "(strongly occupied, orange) or 0 (weakly occupied, blue), labelled by *C*_{2v} "
          "irrep. **b**, Energy error relative to CASCI during optimization from the "
          "Hartree–Fock start for UCCSD (dashed) and GVA (solid) with one to three layers. "
          "Dash-dotted lines give the CISD, CCSD and CCSD(T) errors; the dotted line marks "
          "1.6 mHa. **c**, ADAPT-VQE with the GVA pool: error after each operator is added "
          "and the parameters are re-optimized.")
def _mv_resources(ms, cite):
    ms.figure("figures/Fig5.png", "Fig. 5 | Accuracy against circuit resources for cisplatin.",
              "Energy error relative to CASCI against CNOT count (**a**) and parameter count "
              "(**b**). Symbol size grows with the number of layers. Error bars span all "
              "starts. Diamonds show single-layer circuits with only screening or only tying; "
              "the triangle is ADAPT-VQE. Dash-dotted lines mark the CCSD and CCSD(T) errors "
              "and the dotted line 1.6 mHa.", width_mm=140)



MOVED['resources'] = _mv_resources

# ---------------- condensed summary of the sections moved to the SI
cis_rows = {(r["ansatz"], r["layers"]): r for r in noise if r["system"].startswith("Cisplatin")}
gl2 = cis_rows[("GVA", 2)]
gl1 = cis_rows[("GVA", 1)]
ul2 = cis_rows[("UCCSD", 2)]
mc_ = meas["Cisplatin CAS(10e,10o)"]
pes_c0 = min(r["c0_sq"] for r in pes)
pes_g2 = max(abs(r["GVA_2L"]["dE"]) for r in pes) * 1e3
pes_ccsd = max(abs(r["E_CCSD"] - r["E_CASCI"]) for r in pes) * 1e3
ms.heading("Measurement, scaling, noise and bond stretching", 2)
ms.para(
    "Four further analyses of the molecular systems are reported in the Supplementary "
    "Information. When the Pauli terms of the cisplatin Hamiltonian were grouped by general "
    f"commutativity instead of qubit-wise commutativity{cite(16, 17, 18, 19, 20)}, the number "
    f"of groups fell by a factor of {best(mc_, 'QWC', 'n_groups') / best(mc_, 'GC', 'n_groups'):.0f} "
    "but the number of shots for a 1.6 mHa standard error only by a factor of "
    f"{best(mc_, 'QWC', 'shots') / best(mc_, 'GC', 'shots'):.1f} (Supplementary Note 7 and "
    f"Fig. 3). For cisplatin active spaces of {F['scal'][0]['n_active_orbitals']} to "
    f"{F['scal'][-1]['n_active_orbitals']} orbitals the GVA/UCCSD ratios fall towards the "
    "values expected for a group of order four and reach "
    f"{F['scal'][-1]['GVA_params'] / F['scal'][-1]['UCCSD_params']:.2f} for parameters and "
    f"{F['scal'][-1]['GVA_cnot_1L'] / F['scal'][-1]['UCCSD_cnot_1L']:.2f} for CNOT gates at "
    f"{F['scal'][-1]['n_active_orbitals']} orbitals (Supplementary Note 4, Fig. 6). Under a global depolarizing model, two GVA "
    f"layers tolerate a two-qubit error probability of {sci(gl2['p_max_unmitigated'])} "
    f"without mitigation and {sci(gl2['p_max_rescaling_overhead'])} with idealized "
    f"mitigation{cite(22, 23, 24)}, far below the infidelity of 1.8×10^{{−3}} reported for "
    f"a trapped-ion processor{cite(44)} (Supplementary Note 8, Fig. 7). Stretching one "
    f"Pt–Cl bond to {max(r['R_PtCl'] for r in pes):.1f} Å lowers the symmetry to *C*_{{s}}, "
    f"which leaves {pes[0]['counts']['GVA']['n_params']} GVA parameters per layer instead "
    f"of {cis['counts']['GVA']['n_params']}; in this active space the reference weight "
    f"stays above {pes_c0:.2f}, CCSD stays within {pes_ccsd:.2f} mHa of CASCI and "
    f"two-layer GVA within {pes_g2:.3f} mHa (Supplementary Note 9, Fig. 8).",
    align="justify")

# ---------------- scaling (Supplementary Note 4)
def _mv_scaling(ms, cite):
    sc = F["scal"]
    last = sc[-1]
    ms.para(
        "To see how the savings develop, we counted parameters and CNOT gates for cisplatin "
        f"active spaces from CAS({sc[0]['n_active_orbitals']}e,{sc[0]['n_active_orbitals']}o) "
        f"to CAS({last['n_active_orbitals']}e,{last['n_active_orbitals']}o), all drawn from "
        "the same set of MP2 natural orbitals (Fig. 6). Fitted over the six largest spaces, "
        "UCCSD parameters grow as *N*^{"
        f"{F['slope_UCCSD_params']:.2f}" "} and GVA parameters as *N*^{"
        f"{F['slope_GVA_params']:.2f}" "} with the number of active orbitals *N*; the "
        "corresponding CNOT exponents are "
        f"{F['slope_UCCSD_cnot_1L']:.2f} and {F['slope_GVA_cnot_1L']:.2f}. The GVA/UCCSD "
        f"parameter ratio settles near {np.mean(F['ratio_params'][-4:]):.2f} and the CNOT ratio "
        f"near {np.mean(F['ratio_cnot'][-4:]):.2f}, close to the 1/(2|G|) and 1/|G| values "
        "expected for a group of order |G| = 4 when orbitals are spread evenly over irreps. "
        "The smaller fitted exponent for GVA reflects a ratio that is still falling over this "
        "range; once the ratio has settled, symmetry lowers the prefactor and leaves the "
        "asymptotic growth unchanged. "
        f"At {last['n_qubits']} qubits a single GVA layer still has "
        f"{last['GVA_params']:,} parameters and {last['GVA_cnot_1L']:,} CNOT gates.",
        align="justify")
    ms.figure("figures/Fig6.png", "Fig. 6 | Resource scaling for cisplatin active spaces.",
              "Parameters (**a**) and CNOT gates (**b**) per layer for UCCSD and GVA as the "
              "number of active MP2 natural orbitals grows (half occupied, half virtual), on "
              "logarithmic axes. Slopes are fitted to the six largest spaces. Grey squares "
              "(right axis) give the GVA/UCCSD ratio.", width_mm=130)



MOVED['scaling'] = _mv_scaling

# ---------------- noise (Supplementary Note 8)
cis_rows = {(r["ansatz"], r["layers"]): r for r in noise if r["system"].startswith("Cisplatin")}
gl2 = cis_rows[("GVA", 2)]
gl1 = cis_rows[("GVA", 1)]
ul2 = cis_rows[("UCCSD", 2)]


def _mv_noise(ms, cite):
    ms.para(
        "Lower gate counts matter only if they bring a circuit within reach of real hardware. "
        "We estimated the tolerable two-qubit error with a global depolarizing approximation "
        "in which each CNOT carries depolarizing probability *p*, so that the energy relaxes "
        "towards the mean eigenvalue Tr(Ĥ)/2^{*n*} as the circuit fidelity "
        "*F* = (1 − 16*p*/15)^{*N*} decays, *N* being the number of CNOT gates (Methods). "
        "This model ignores the "
        "structure of real noise and serves only as an order-of-magnitude budget. For "
        "cisplatin, Tr(Ĥ)/2^{*n*} lies "
        f"{gl2['gap_Ha']:.2f} Ha above the ground state, so even a small loss of fidelity "
        "produces a large bias (Fig. 7a). Keeping the unmitigated bias below 1.6 mHa requires "
        f"*p* ≤ {sci(gl1['p_max_unmitigated'])} for one GVA layer and "
        f"*p* ≤ {sci(gl2['p_max_unmitigated'])} for two, against {sci(ul2['p_max_unmitigated'])} "
        "for two UCCSD layers. Richardson extrapolation over noise scale factors 1, 3 and 5"
        f"{cite(22, 23, 24)} relaxes this to {sci(gl2['p_max_richardson'])}, and exact "
        "rescaling by a known fidelity with a sampling overhead of at most 100 tolerates "
        f"*p* ≤ {sci(gl2['p_max_rescaling_overhead'])}. For comparison, an average two-qubit "
        f"gate infidelity of 1.8×10^{{−3}} has been reported for a trapped-ion processor{cite(44)}, "
        "above all of these thresholds and, without mitigation, by more than "
        f"{int(np.floor(np.log10(1.84e-3 / gl2['p_max_unmitigated'])))} orders of magnitude (Fig. 7b). The threefold gate "
        "reduction from symmetry raises the tolerable error rate by the same factor, which "
        "is useful for resource estimates but does not by itself move cisplatin-sized "
        f"circuits into the near-term regime{cite(37)}.",
        align="justify")
    ms.figure("figures/Fig7.png", "Fig. 7 | Error budget under a global depolarizing model.",
              "**a**, Unmitigated energy bias for cisplatin circuits as a function of the "
              "two-qubit depolarizing probability *p*. The dashed line marks the average "
              f"two-qubit infidelity of 1.8×10^{{−3}} reported for a trapped-ion processor{cite(44)}; the dotted "
              "line marks 1.6 mHa. **b**, "
              "Largest *p* that keeps the bias below 1.6 mHa without mitigation, with "
              "Richardson extrapolation (λ = 1, 3, 5), or with exact fidelity rescaling at a "
              "sampling overhead of at most 100, for the cheapest GVA circuit that reaches "
              "1.6 mHa in noiseless simulation (layers and CNOT count below each label). The "
              "dashed line repeats the reported infidelity from **a**.",
              width_mm=140)



MOVED['noise'] = _mv_noise

# ---------------- PES (Supplementary Note 9)
def _mv_pes(ms, cite):
    R = [r["R_PtCl"] for r in pes]
    c0min = min(pes, key=lambda r: r["c0_sq"])
    ccsd_err = [(r["E_CCSD"] - r["E_CASCI"]) * 1e3 for r in pes]
    ccsdt_err = [(r["E_CCSD(T)"] - r["E_CASCI"]) * 1e3 for r in pes]
    g2err = [r["GVA_2L"]["dE"] * 1e3 for r in pes]
    g1err = [r["GVA_1L"]["dE"] * 1e3 for r in pes]
    u1err = [r["UCCSD_1L"]["dE"] * 1e3 for r in pes]
    iw = int(np.argmax(np.abs(ccsd_err)))
    i_eq = [k for k, r in enumerate(pes) if abs(r["R_PtCl"] - 2.32) < 1e-6][0]
    rise = (pes[-1]["E_CASCI"] - pes[i_eq]["E_CASCI"]) * HA2KCAL
    ecorr = [(r["E_CASCI"] - r["E_ref"]) * 1e3 for r in pes]
    ms.para(
        f"We stretched one Pt–Cl bond from {min(R):.1f} to {max(R):.1f} Å with the "
        "other ligands fixed, which lowers the symmetry to *C*_{s}. The scan describes "
        "gas-phase bond dissociation; aquation, the first step in the activation of "
        "cisplatin, would need an incoming water molecule and solvent and lies outside this "
        "model. The active orbitals were selected at 2.32 Å and followed by maximum overlap "
        f"along the scan (Methods). Between 2.32 and {max(R):.1f} Å the CASCI energy rose by "
        f"{rise:.0f} kcal mol^{{−1}} and the active-space correlation energy changed from "
        f"{ecorr[i_eq]:.1f} to {ecorr[-1]:.1f} mHa, but the reference weight stayed between "
        f"{min(r['c0_sq'] for r in pes):.3f} and {max(r['c0_sq'] for r in pes):.3f} (Fig. 8c). "
        "Within this frozen-core space built on a closed-shell reference, the stretched bond "
        "therefore does not develop strong static correlation: CCSD stayed within "
        f"{max(abs(x) for x in ccsd_err):.2f} mHa of CASCI and CCSD(T) within "
        f"{max(abs(x) for x in ccsdt_err):.2f} mHa. Two-layer GVA stayed within "
        f"{max(abs(x) for x in g2err):.3f} mHa at every point, single-layer GVA within "
        f"{max(abs(x) for x in g1err):.3f} mHa and single-layer UCCSD within "
        f"{max(abs(x) for x in u1err):.3f} mHa (Fig. 8b). A description of homolytic "
        "cleavage would need the Pt–Cl σ and σ^{*} orbitals to be optimized at each geometry, "
        "for example by CASSCF, which we did not attempt. In *C*_{s} the screening is "
        f"weaker: GVA used {pes[0]['counts']['GVA']['n_params']} parameters per layer "
        f"against {pes[0]['counts']['UCCSD']['n_params']} for UCCSD.",
        align="justify")
    ms.figure("figures/Fig8.png", "Fig. 8 | Pt–Cl bond stretching in cisplatin, CAS(10e,10o).",
              "**a**, Energies relative to the CASCI minimum along the Pt–Cl distance "
              "(*C*_{s}). **b**, Errors relative to CASCI for CISD, CCSD, CCSD(T), "
              "single-layer UCCSD and GVA, and two-layer GVA (Hartree–Fock start); the dotted "
              "line marks 1.6 mHa. The single-layer UCCSD and GVA curves coincide on this "
              "scale. **c**, Weight of the reference determinant in the CASCI state (left axis) "
              "and active-space correlation energy (right axis).")



MOVED['pes'] = _mv_pes

# ------------------------------------------------------------------ crystals
MP.results(ms, F)

# ------------------------------------------------------------------ discussion
ms.heading("Discussion")
ms.para(
    "Point-group screening is exact in the sense that matters for optimization: a "
    "forbidden generator has zero gradient along every symmetric trajectory, so removing "
    "it changes the parameter count and the gate count but not the set of states the "
    "optimizer can reach from the Hartree–Fock start. Crystal-momentum screening is exact "
    "for the same reason. Spin-flip tying is an "
    "approximation, but its cost was below "
    f"{F['max_gva_penalty_1L_mHa'] * 1e3:.0f} µHa in every case we examined, and it "
    "removed a failure mode that we did not anticipate: from random starts, unconstrained "
    "UCCSD found spin-contaminated states with lower energy than the spin-pure optimum. "
    "Monitoring ⟨Ŝ^{2}⟩ or restricting the ansatz is therefore necessary when random "
    "initialization is used; the symmetry constraints remove spurious minima as well as "
    "redundant parameters.",
    align="justify")
ms.para(
    "None of the constraints is new on its own. Point-group restrictions of excitation "
    f"lists{cite(11)}, symmetry-based qubit tapering{cite(12, 13)}, spin-adapted "
    f"pools{cite(6, 14)} and momentum-conserving unitary coupled cluster for "
    f"crystals{cite(50, 51)} are established, and tapering could be applied on top of GVA "
    "to remove further qubits. Our contribution is a quantitative account of what each "
    "constraint saves and costs, a closed form for the counts that holds for any abelian "
    "group and shows that the saving is set by the group order, time-reversal tying of "
    "complex amplitudes, and a benchmark in which every number can be regenerated from "
    "the released scripts.",
    align="justify")
ms.para(
    "Several limitations apply. We use abelian point groups only; degenerate orbitals in "
    "non-abelian groups would allow further tying with fixed relative signs, which we "
    "have not implemented. All circuits were simulated exactly, without shot noise or "
    "gate errors, and the CNOT counts use a simple ladder compilation; optimized "
    f"fermionic-excitation circuits{cite(15)} and better orbital orderings would lower "
    "them without changing the comparison between ansätze. We did not benchmark "
    "particle-conserving hardware-efficient fabrics at the cisplatin scale. Absolute "
    "errors depend on how the active orbitals are chosen (Supplementary Note 10), although "
    "the comparison between GVA and UCCSD did not. The cisplatin "
    "active space is weakly correlated at equilibrium, where CCSD(T) is already accurate "
    "to "
    f"{abs(F['cis_E_CCSD(T)_err_mHa']):.2f} mHa, so these results say nothing about a "
    f"quantum advantage for this molecule{cite(31)}. The periodic tests use "
    "one-dimensional chains in a minimal basis with at most five cells, do not correct "
    "finite-size errors, and use UCCSD on the supercell as the only reference ansatz.",
    align="justify")
cc = cis["counts"]
s_big = F["scal"][-1]
gva1 = F["cis_med"][("GVA", 1)]
mc4 = MP.MC[4]
ms.para(
    "The size of the saving is set by the order of the symmetry group, and this decides "
    "where GVA is useful. For molecules the order is at most eight. Screening is then "
    "most effective for symmetric metal complexes, macrocycles and idealized cluster "
    "models of metalloenzyme active sites, whereas most drug-like molecules and "
    f"active-site fragments cut from protein structures{cite(45)} have no spatial "
    "symmetry. In *C*_{1} only spin-flip tying remains: treating the cisplatin active "
    f"space as *C*_{{1}}, tying alone reduced the parameter count from {cc['UCCSD']['n_params']} "
    f"to {cc['spin only']['n_params']} and left the CNOT count at "
    f"{cc['spin only']['cnot_1L']:,}, because each tied pair is still applied as two "
    "exponentials. For such fragments, which in practice would be defined by QM/MM "
    f"partitioning or an embedding scheme{cite(46, 47)}, GVA has to be combined with "
    "methods that shorten the operator list by other criteria, such as ADAPT-VQE (here "
    f"{ad['dE'] * 1e3:.3f} mHa with {ad['n_params']} operators and {ad['cnot']:,} CNOT "
    f"gates) or the removal of excitations between distant localized orbitals{cite(48, 49)}. "
    "Approximate local symmetry should not be used for screening, because a generator that "
    "is forbidden only approximately has a small but non-zero gradient. Crystals are the "
    "opposite case. The translation group is exact and its order grows with the "
    "*k*-point mesh, so the reduction in parameters and CNOT gates grows linearly with "
    "*N*_{*k*} and does not depend on any spatial symmetry of the cell. For crystals of "
    "pharmaceutical interest, polymorphs of molecular solids whose lattice energies often "
    f"differ by only a few kJ mol^{{−1}}{cite(55)}, each cell holds several molecules and "
    "would need a large active space. On a 4×4×4 mesh with eight occupied and eight "
    f"virtual bands per cell, GVA-k still has {MP.sci(mc4['GVA_k_params'])} parameters per "
    "layer, so the factor *N*_{*k*} has to be combined with a reduction of the per-cell "
    "space before such crystals become accessible.",
    align="justify")
ms.para(
    "The error budget sets the practical context. A two-layer GVA circuit for a 20-qubit "
    f"active space needs {gl2['cnot']:,} CNOT gates in the ladder compilation, and under "
    f"the depolarizing model it tolerates *p* ≈ {sci(gl2['p_max_unmitigated'])} without "
    f"mitigation and {sci(gl2['p_max_rescaling_overhead'])} with idealized mitigation. "
    "For molecules, symmetry screening is therefore best viewed as a constant-factor "
    f"saving for early fault-tolerant resource estimates{cite(32)} and for classical "
    "simulation of variational circuits. For crystals the factor grows with the number "
    "of *k*-points, although the scaling with the size of the per-cell active space is "
    "unchanged. The suppression of triplet contamination is a practical benefit at any "
    "scale.",
    align="justify")

# ------------------------------------------------------------------ methods
ms.heading("Methods")
ms.heading("Molecular systems and integrals", 2)
ms.para(
    f"Integrals and reference calculations used PySCF 2.14{cite(25)}. "
    "H_{2} (*R* = 0.7414 Å), linear H_{4} (equal spacing 1.2300 Å), LiH (*R* = 1.5950 Å) "
    "and BeH_{2} (*R* = 1.3264 Å) were treated in the STO-3G basis with all orbitals "
    "active (4, 8, 12 and 14 qubits); stretched geometries scale every bond by 1.5 or "
    "2.0. Orbital symmetries were assigned in *D*_{2h} (H_{2}, H_{4}, BeH_{2}) or "
    "*C*_{2v} (LiH). Cisplatin was built with Pt at the origin, a square-planar PtN_{2}Cl_{2} "
    "core (Cl–Pt–Cl and N–Pt–N 90°), N–H 1.02 Å and Pt–N–H 109.5°, with one N–H bond of "
    "each ammine in the molecular plane (Supplementary Table S1). All atoms used def2-SVP "
    "and Pt the 60-electron def2 effective core potential, giving 72 explicit electrons "
    f"and 126 basis functions; the restricted Hartree–Fock (RHF) energy is "
    f"{cis['E_RHF']:.6f} Ha.",
    align="justify")
ms.heading("Active space", 2)
ms.para(
    "The MP2 one-particle density matrix was diagonalized within each irrep block to give "
    "symmetry-adapted natural orbitals, which were ranked by occupation. The active space "
    "contains the five strongly occupied natural orbitals with the lowest occupations and "
    f"the five weakly occupied ones with the highest occupations, numbers {nc + 1}–{nc + 10} "
    f"of {cis['n_ao']} ({irr_txt(orb_irr)} in *C*_{{2v}}); their occupations and atomic "
    "composition are listed in Supplementary Table S9. The other "
    f"{nc} strongly occupied orbitals form a frozen core whose mean field enters the "
    "one-electron active-space Hamiltonian. "
    "Inside the active space a closed-shell self-consistent field calculation, with "
    "rotations restricted to each irrep block, gives the canonical active orbitals and "
    "the reference determinant; its energy reproduces the full RHF energy to 10^{−9} Ha. "
    "For the Pt–Cl scan (*C*_{s}, 2.0–3.8 Å) the space was chosen in the same way at "
    "2.32 Å and then propagated outward and inward by selecting, at each new geometry, "
    "the strongly and weakly occupied natural orbitals with the largest projection onto "
    "the previous active space, using the overlap matrix between the two basis sets. "
    "For comparison at the equilibrium geometry, two further active spaces of the same "
    f"size were built: canonical RHF orbitals {OS['B']['selected'].split(' ')[2].replace('-', '–')} "
    "(HOMO−4 to LUMO+4), and CASSCF(10e,10o) natural orbitals obtained by full "
    "orbital optimization of the ground state starting from the MP2 natural orbitals "
    f"(*E*_{{CASSCF}} = {OS['C']['E_CASSCF']:.6f} Ha). The coefficients of the active "
    "orbitals used in the main calculations are distributed with the code.",
    align="justify")
ms.heading("Ansätze and closed-form counts", 2)
ms.para(
    "The UCCSD pool contains all spin-conserving singles and doubles relative to the "
    "reference, ordered as singles, same-spin doubles and αβ doubles. For *n*_{o} "
    "occupied and *n*_{v} virtual spatial orbitals the numbers of parameters are",
    align="justify")
ms.equation("*N*_{UCCSD} = 2*n*_{o}*n*_{v} + 2C(*n*_{o},2)C(*n*_{v},2) + (*n*_{o}*n*_{v})^{2},", 2)
ms.equation("*N*_{spin} = *n*_{o}*n*_{v} + C(*n*_{o},2)C(*n*_{v},2) + *n*_{o}*n*_{v}(*n*_{o}*n*_{v} + 1)/2,", 3)
ms.equation("*N*_{GVA} = *S* + *D*_{ss} + (*D*_{αβ} + *n*_{o}*n*_{v})/2,", 4)
ms.para(
    "where *S* = Σ_{Γ} *o*_{Γ}*v*_{Γ} counts allowed singles per spin, *D*_{αβ} = "
    "Σ_{Γ} [Σ_{Γ'⊗Γ''=Γ} *o*_{Γ'}*v*_{Γ''}]^{2} allowed αβ doubles and *D*_{ss} = "
    "Σ_{Γ} *P*_{o}(Γ)*P*_{v}(Γ) allowed doubles per spin, with *o*_{Γ} and *v*_{Γ} the "
    "numbers of occupied and virtual orbitals in irrep Γ and *P*_{o}(Γ), *P*_{v}(Γ) the "
    "numbers of occupied and virtual orbital pairs whose product is Γ. The term "
    "*n*_{o}*n*_{v} in equation (4) counts the αβ doubles *i*_{α}*i*_{β}→*a*_{α}*a*_{β}, "
    "which are their own spin-flip partners. Screening alone gives 2*S* + 2*D*_{ss} + "
    "*D*_{αβ}. Layers repeat the full product with new parameters. The reference "
    "determinant is the initial state in all cases.",
    align="justify")
ms.heading("Simulation and optimization", 2)
ms.para(
    "States were represented exactly in the space of determinants with fixed numbers of "
    "α and β electrons (63,504 determinants for cisplatin, equivalent to the "
    "particle-number- and *S*_{z}-conserving sector of the 20-qubit register). For a "
    "generator with T̂|*u*⟩ = *s*|*d*⟩, exp(θĜ) maps |*u*⟩ → cos θ|*u*⟩ + *s* sin θ|*d*⟩ "
    "and |*d*⟩ → cos θ|*d*⟩ − *s* sin θ|*u*⟩ and leaves other determinants unchanged; "
    "signs were taken from the PySCF string tables so that generators and Hamiltonian "
    "share one convention. The Hamiltonian was applied with the PySCF direct CI kernel. "
    "Mapping the CASCI vector to the Jordan–Wigner basis reproduced the CASCI energy "
    "with the OpenFermion qubit Hamiltonian"
    f"{cite(26)} to within 10^{{−13}} Ha for the small molecules and {JWC} Ha for "
    "cisplatin. Gradients were computed exactly by reverse-mode "
    "(adjoint) differentiation through the product of rotations, which equals the "
    f"infinite-shot parameter-shift gradient{cite(36)}, and checked against central "
    "finite differences. Parameters were optimized with L-BFGS-B"
    f"{cite(35)} (SciPy 1.18) until the projected gradient fell below 10^{{−7}} or the "
    "relative energy change below 10^{−12}, with at most 3,000 (small molecules) or "
    "2,500 (cisplatin) iterations. Seed 0 starts from **θ** = 0; other seeds draw "
    "θ ~ N(0, 0.05^{2}) from NumPy generators seeded with 1000 + *k* (small molecules) or "
    "2000 + *k* (cisplatin). ADAPT-VQE added at each step the pool member with the "
    "largest energy gradient and re-optimized all parameters, stopping when the gradient "
    "norm fell below 10^{−3} Ha, when five additions lowered the energy by less than "
    "10^{−7} Ha, or at 120 operators. CISD, CCSD and CCSD(T) were run with PySCF in the "
    "same active space, and FCI or CASCI was converged to 10^{−12} Ha.",
    align="justify")
ms.heading("Circuit counts", 2)
ms.para(
    "With interleaved spin orbitals (α of orbital *p* on qubit 2*p*, β on 2*p* + 1), "
    "the Jordan–Wigner image of a single excitation between spin orbitals *p* < *q* "
    "contains two Pauli strings of weight *q* − *p* + 1, and a double excitation with "
    "sorted indices *p* < *q* < *r* < *s* contains eight strings of weight "
    "(*q* − *p* + 1) + (*s* − *r* + 1). Compiling each Pauli exponential with a CNOT "
    "ladder costs 2(*w* − 1) CNOT gates for weight *w*, giving",
    align="justify")
ms.equation("*N*_{CNOT}(single) = 4(*q* − *p*),   *N*_{CNOT}(double) = 16(*q* − *p* + *s* − *r* + 1)", 5)
ms.para(
    "per generator, independent of tying. All-to-all connectivity is assumed. A complex "
    "amplitude adds two single-qubit rotations and no CNOT gates (Supplementary Note 11).",
    align="justify")
MP.methods(ms)
ms.heading("Measurement grouping and shots", 2)
ms.para(
    "Pauli terms were stored as binary symplectic vectors. Two terms conflict under QWC "
    "if they act with different non-identity Paulis on a common qubit and under GC if "
    "they anticommute. Groups were formed by largest-degree-first greedy colouring, by "
    f"DSATUR{cite(42)} and by sorted insertion{cite(19)} of the conflict graph. For each "
    "grouping, the standard deviation σ_{*g*} of every group operator was evaluated "
    "exactly in the FCI or CASCI ground state, and the shots for standard error ε on one "
    f"energy evaluation follow from optimal allocation{cite(20)} as "
    "*N* = (Σ_{*g*} σ_{*g*})^{2}/ε^{2}, with ε = 1.6 mHa. Covariances between groups do "
    "not enter because groups are measured in separate circuits.",
    align="justify")
ms.heading("Depolarizing error budget", 2)
ms.para(
    "Each CNOT is followed by two-qubit depolarizing noise that applies one of the 15 "
    "non-identity Pauli pairs with total probability *p*, contracting the Bloch "
    "components by *f* = 1 − 16*p*/15. Approximating the accumulated noise as global, "
    "the state becomes *F*ρ + (1 − *F*)𝟙/2^{*n*} with *F* = *f*^{*N*} for *N* CNOT "
    "gates, and the "
    "measured energy is *E*(*p*) = *F**E* + (1 − *F*)*c*_{0}, where *c*_{0} = "
    "Tr(Ĥ)/2^{*n*} is the identity coefficient of the qubit Hamiltonian. The "
    "unmitigated bias is (1 − *F*)(*c*_{0} − *E*). For zero-noise extrapolation, unitary "
    f"folding{cite(24)} with scale factors λ = 1, 3, 5 and Richardson weights "
    "(15/8, −5/4, 3/8) leaves a bias |1 − Σ_{λ}γ_{λ}*F*^{λ}|(*c*_{0} − *E*). If *F* were "
    "known exactly, rescaling would remove the bias at a sampling cost of 1/*F*^{2}. "
    "Thresholds were found by bisection on log *p*.",
    align="justify")
ms.heading("Software and hardware", 2)
ms.para(
    "Calculations used Python 3.13, NumPy 2.5, SciPy 1.18, PySCF 2.14 and OpenFermion "
    "1.8.1 on a cloud virtual machine with two CPU cores and 8 GB memory. No quantum "
    "hardware was used.",
    align="justify")

# ------------------------------------------------------------------ end matter
ms.heading("Code availability", 2)
ms.para(
    "The scripts that generate all results and figures (gva_repro, including uccsd_counts.py "
    "for the closed-form counts) are provided as Supplementary Software and will be "
    "deposited in a public repository with a persistent identifier. [Repository URL, "
    "release tag and licence to be added.]",
    align="justify")
ms.heading("Data availability")
ms.para(
    "Source data for every figure and table, all optimization runs with their energy "
    "histories and final variational parameters, the cisplatin active-space integrals, "
    "the Bloch-basis and real-basis integrals of every periodic system and the Cartesian "
    "coordinates are provided as Supplementary Data 1 (GVA_Source_Data.xlsx) and in the "
    "code repository, together with full PySCF output of the molecular and periodic "
    "reference calculations and the run logs. Every final state can be rebuilt from the "
    "stored parameters with the scripts verify_theta.py and verify_periodic.py, which "
    "reproduce all stored energies. [Repository URL and DOI to be added on deposition.]",
    align="justify")
ms.heading("References")
ms.references(REFS.formatted())
ms.heading("Acknowledgements")
ms.para("[To be completed by the authors.]")
ms.heading("Author contributions")
ms.para("[To be completed by the authors.]")
ms.heading("Competing interests")
ms.para("[To be completed by the authors.]")
ms.save("GVA_npjQI_Manuscript.docx")
print("manuscript written")
