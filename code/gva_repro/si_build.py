"""Build the Supplementary Information (docx) from results/."""
import json

import numpy as np

from docbuilder import Manuscript
from facts import compute
from refs import MASTER, REFS, SREFS, citeS
import numbering
import ms_build  # noqa: F401  builds the manuscript first so that main-text reference numbers are fixed

F = compute()
cis = F["cis"]
g = F["small_groups"]
sd = F["sysd"]
orb = json.load(open("results/cisplatin_orbitals.json"))


def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}×10^{{{int(e)}}}"


def esci(x, d=1, floor=1e-9):
    """Energy error in mHa; values below the numerical floor are shown as a bound."""
    return "<10^{−9}" if abs(x) < floor else sci(x, d)


MT = {"H2": "H_{2}", "H4": "H_{4}", "LiH": "LiH", "BeH2": "BeH_{2}"}
AN = {"PG only": "Screening only", "spin only": "Tying only"}


def mt(name):
    return MT.get(name, name)

si = Manuscript("Supplementary Information: symmetry-screened UCC circuits", author="",
                line_numbers=False, spacing=1.15, size=10.5)
si.remap = numbering.remap
si.title("Supplementary Information", size=15)
si.para("**Symmetry-screened unitary coupled-cluster circuits for molecules and crystals**",
        size=12, space_after=4)
si.para("[Author names]", space_after=10)

si.heading("Contents")
contents = [
    "Supplementary Note 1. Generators, sign conventions and the rotation identity",
    "Supplementary Note 2. Zero gradient of symmetry-forbidden generators",
    "Supplementary Note 3. Spin-flip tying and spin contamination",
    "Supplementary Note 4. Closed-form parameter counts",
    "Supplementary Note 5. CNOT counts under the Jordan–Wigner mapping",
    "Supplementary Note 6. Cisplatin geometry and active space",
    "Supplementary Note 7. Measurement grouping and shot estimates",
    "Supplementary Note 8. Depolarizing error budget and zero-noise extrapolation",
    "Supplementary Note 9. Pt–Cl bond-stretching scan",
    "Supplementary Note 10. Dependence on the choice of active orbitals",
    "Supplementary Note 11. Crystal momentum and time reversal",
    "Supplementary Note 12. Periodic integrals and their validation",
    "Supplementary Note 13. Counts for crystals",
    "Supplementary Note 14. Benchmark details for the chains",
    "Supplementary Note 15. Measurement in the Bloch and real bases",
    "Supplementary Figs. 1–7",
    "Supplementary Tables S1–S15",
    "Supplementary References",
    "Supplementary Data 1 and Supplementary Software",
]
for c in contents:
    si.para(c, space_after=1)

# ------------------------------------------------------------------ Note 1
si.heading("Supplementary Note 1. Generators, sign conventions and the rotation identity")
si.para(
    "Spin orbitals are labelled by spatial index *p* and spin σ ∈ {α, β}. Determinants "
    "are ordered as in PySCF: |*I*_{α}*I*_{β}⟩ = (∏ *a*^{†}_{*p*α})(∏ *a*^{†}_{*q*β})|0⟩ "
    "with creation operators in ascending order inside each spin block. Single "
    "excitations *E*^{σ}_{*ai*} = *a*^{†}_{*a*σ}*a*_{*i*σ} act on strings with the "
    "signs of the PySCF link tables. For distinct indices, the same-spin double "
    "excitation satisfies *a*^{†}_{*a*}*a*^{†}_{*b*}*a*_{*j*}*a*_{*i*} = "
    "*E*_{*ai*}*E*_{*bj*}, and the αβ double satisfies *a*^{†}_{*a*α}*a*^{†}_{*b*β}"
    "*a*_{*j*β}*a*_{*i*α} = *E*^{α}_{*ai*}*E*^{β}_{*bj*}; both identities follow from "
    "two anticommutations and fix every sign without separate bookkeeping.",
    align="justify")
si.para(
    "Because occupied and virtual index sets are disjoint, T̂^{2} = 0 for every "
    "generator. Writing Ĝ = T̂ − T̂^{†}, one has Ĝ^{2} = −(T̂T̂^{†} + T̂^{†}T̂) = −P̂, where "
    "P̂ projects onto the determinants connected by T̂, and Ĝ^{3} = −Ĝ. Hence",
    align="justify")
si.equation("exp(θĜ) = 1 + sin θ Ĝ + (1 − cos θ) Ĝ^{2},", "S1")
si.para(
    "which rotates each connected pair (|*u*⟩, |*d*⟩) with T̂|*u*⟩ = *s*|*d*⟩, *s* = ±1, "
    "by the angle θ and leaves all other determinants unchanged. The energy gradient "
    "with respect to θ_{*k*} is 2⟨λ_{*k*}|Ĝ_{*k*}|ψ_{*k*}⟩, where |ψ_{*k*}⟩ is the state "
    "after the *k*-th factor and ⟨λ_{*k*}| = ⟨Ψ|Ĥ*U*_{>*k*}; both are obtained in one "
    "backward sweep. Analytic gradients agreed with central finite differences "
    "(step 10^{−6}) to eight significant figures for all four small molecules. The "
    "determinant-space vector was mapped to the Jordan–Wigner basis with the reordering "
    "sign (−1)^{*n*}, where *n* counts pairs (*i* ∈ α, *j* ∈ β) with *i* > *j*; the qubit "
    "Hamiltonian then reproduced the FCI energies of the small molecules to better than "
    "10^{−13} Ha and the cisplatin CASCI energy within the 10^{−8} Ha tolerance asserted "
    "by run_measure.py.",
    align="justify")

# ------------------------------------------------------------------ Note 2
si.heading("Supplementary Note 2. Zero gradient of symmetry-forbidden generators")
si.para(
    "Let *G* be an abelian point group with unitary representation *R*(*g*) on Fock "
    "space, [Ĥ, *R*(*g*)] = 0, and let |Φ_{0}⟩ be totally symmetric, *R*(*g*)|Φ_{0}⟩ = "
    "|Φ_{0}⟩. A generator built from orbitals with irreps Γ_{*a*}, Γ_{*b*}, Γ_{*i*}, "
    "Γ_{*j*} transforms as *R*(*g*)Ĝ*R*(*g*)^{†} = χ(*g*)Ĝ with χ = χ_{Γa}χ_{Γb}χ_{Γi}"
    "χ_{Γj} = ±1. Suppose that every generator except Ĝ_{*k*} is totally symmetric and "
    "θ_{*k*} = 0. Then |ψ_{*k*}⟩ and the dressed Hamiltonian H̃_{*k*} = *U*_{>*k*}^{†}"
    "Ĥ*U*_{>*k*} are invariant under *R*(*g*), and",
    align="justify")
si.equation("∂*E*/∂θ_{*k*} = ⟨ψ_{*k*}|[H̃_{*k*}, Ĝ_{*k*}]|ψ_{*k*}⟩ = χ(*g*)⟨ψ_{*k*}|[H̃_{*k*}, Ĝ_{*k*}]|ψ_{*k*}⟩.", "S2")
si.para(
    "If χ(*g*) = −1 for some *g*, the derivative vanishes. A gradient-based optimizer "
    "started with all forbidden angles at zero therefore never changes them, and the "
    "energy along its trajectory is identical with or without the forbidden factors. "
    "Removing them is exact for optimizations started at **θ** = 0. For random starts the "
    "forbidden angles are non-zero and the state can leave the totally symmetric sector. "
    "We did not analyse the spatial symmetry of the random-start UCCSD solutions "
    "separately; Supplementary Note 3 shows that leaving the spin-symmetric sector does "
    "lower the energy of a truncated ansatz in two stretched geometries.",
    align="justify")

# ------------------------------------------------------------------ Note 3
cont = F["spin_contaminated"]
si.heading("Supplementary Note 3. Spin-flip tying and spin contamination")
si.para(
    "For a closed-shell singlet the exact CI vector satisfies *c*(*I*_{α}, *I*_{β}) = "
    "*c*(*I*_{β}, *I*_{α}); in the same convention triplet vectors are antisymmetric "
    "and quintet vectors symmetric, which we verified numerically for the lowest six "
    "roots of H_{4} ×2.0. Let *P̂* exchange α and β strings. A tied pair satisfies "
    "*P̂*Ĝ_{*k*}*P̂* = Ĝ_{*k'*}, so *P̂* maps the GVA product onto the same product with "
    "the members of each pair interchanged. An α single and its β partner act on "
    "different strings and commute, and so do an αα double and its ββ partner, so for "
    "these the interchange changes nothing; αβ partners need not commute, and only their "
    "ordering can break the exchange symmetry. The ansatz therefore suppresses triplet admixture. It does not "
    "enforce *Ŝ*^{2} exactly, because a tied αβ, αα and ββ set is not a complete "
    "spin-adapted double excitation, and the energy cost of tying is the difference "
    "between the GVA and screening-only rows of Supplementary Table S4 (at most "
    f"{F['max_gva_penalty_1L_mHa'] * 1e3:.0f} µHa over all geometries). The largest "
    f"⟨Ŝ^{{2}}⟩ over all {sum(len(v) for k, v in g.items() if k[2] == 'GVA')} GVA runs was "
    f"{sci(F['gva_max_S2'])}. For H_{{4}} ×2.0 the lowest triplet lies "
    f"{F['h4_triplet_gap_mHa']:.1f} mHa above the singlet ground state.",
    align="justify")
txt = "; ".join(sorted({f"{mt(m)} ×{sc:.1f}, {an} {L}L: ⟨Ŝ^{{2}}⟩ = {s2:.3f}"
                        for m, sc, an, L, sd_, s2, de in cont}))
si.para(
    f"Runs with ⟨Ŝ^{{2}}⟩ > 0.01 ({len(cont)} in total, all UCCSD from random starts): "
    f"{txt}. In every case the Hartree–Fock start of the same ansatz converged to a "
    "spin-pure state with a higher energy, and GVA converged to that spin-pure state from "
    "all five starts (Supplementary Fig. 1).",
    align="justify")
si.figure("figures/FigS1.png", "Supplementary Fig. 1 | Spin contamination.",
          "⟨Ŝ^{2}⟩ of the final state against energy error relative to FCI for every "
          "single- and two-layer UCCSD and GVA run on the small molecules "
          f"({sum(len(v) for k, v in g.items() if k[2] in ('UCCSD', 'GVA'))} runs). Open "
          "symbols, UCCSD; filled symbols, GVA. Values of ⟨Ŝ^{2}⟩ below 10^{−12} are drawn "
          "at that floor.", width_mm=110)

# ------------------------------------------------------------------ Note 4
si.heading("Supplementary Note 4. Closed-form parameter counts")
ms_build.MOVED["scaling"](si, citeS)
si.para(
    "With *n*_{o} occupied and *n*_{v} virtual spatial orbitals, the spin-conserving UCCSD "
    "pool has 2*n*_{o}*n*_{v} singles (α and β), 2C(*n*_{o},2)C(*n*_{v},2) same-spin "
    "doubles and (*n*_{o}*n*_{v})^{2} αβ doubles, which gives equation (2) of the main "
    "text. For large spaces *N*_{UCCSD} → (3/2)(*n*_{o}*n*_{v})^{2}. Spin-flip tying "
    "pairs α and β singles and αα with ββ doubles; among the (*n*_{o}*n*_{v})^{2} αβ "
    "doubles, the *n*_{o}*n*_{v} with *i* = *j* and *a* = *b* are self-partnered and the "
    "rest form pairs, which gives equation (3). Point-group screening keeps a single when "
    "Γ_{*a*} = Γ_{*i*}, so *S* = Σ_{Γ}*o*_{Γ}*v*_{Γ} per spin. An αβ double (*i*α, *j*β → "
    "*a*α, *b*β) is allowed when Γ_{*i*}⊗Γ_{*a*} = Γ_{*j*}⊗Γ_{*b*}; with *X*(Γ) the "
    "number of occupied–virtual pairs of product Γ, *D*_{αβ} = Σ_{Γ}*X*(Γ)^{2}. A "
    "same-spin double is allowed when the occupied pair and the virtual pair have equal "
    "products, *D*_{ss} = Σ_{Γ}*P*_{o}(Γ)*P*_{v}(Γ). Self-partnered αβ doubles are always "
    "allowed, giving equation (4). The script uccsd_counts.py evaluates these "
    "expressions and compares them with explicit enumeration; Supplementary Table S5 "
    "lists the results.",
    align="justify")

# ------------------------------------------------------------------ Note 5
si.heading("Supplementary Note 5. CNOT counts under the Jordan–Wigner mapping")
si.para(
    "With interleaved ordering (qubit 2*p* for *p*α, 2*p* + 1 for *p*β), the single "
    "excitation *a*^{†}_{*q*}*a*_{*p*} − h.c. (*p* < *q*) maps to (*i*/2)(*X*_{*p*}"
    "*Z*…*Z* *Y*_{*q*} − *Y*_{*p*}*Z*…*Z* *X*_{*q*}), two commuting strings of weight "
    "*q* − *p* + 1. A double excitation with sorted indices *p* < *q* < *r* < *s* maps "
    "to eight commuting strings, each with *Z* chains between *p* and *q* and between "
    "*r* and *s*; the chain between *q* and *r* cancels, so the weight is (*q* − *p* + 1) "
    "+ (*s* − *r* + 1). Exponentiating a weight-*w* string with a CNOT ladder costs "
    "2(*w* − 1) CNOT gates, which gives equation (5). Tied generators are still applied "
    "individually, so tying does not change the count. The counts assume all-to-all "
    "connectivity and no cancellation between neighbouring ladders; optimized circuits for "
    f"fermionic excitations{citeS(15)} reduce them and would lower UCCSD and GVA "
    "by similar factors.",
    align="justify")

# ------------------------------------------------------------------ Note 6
si.heading("Supplementary Note 6. Cisplatin geometry and active space")
ms_build.MOVED["resources"](si, citeS)
si.para(
    "The geometry is given in Supplementary Table S1. Symmetry detection in PySCF finds "
    f"*C*_{{2v}} at the reference geometry. The RHF energy is {cis['E_RHF']:.6f} Ha and the "
    f"MP2 correlation energy (all non-ECP electrons) is {cis['E_MP2_corr']:.4f} Ha. The "
    "ten active natural orbitals, their occupations and their dominant atomic-orbital "
    "contributions are listed in Supplementary Table S9. In the active space, "
    f"*E*_{{ref}} = {cis['E_ref']:.6f} Ha, *E*_{{CASCI}} = {cis['E_CASCI']:.6f} Ha, "
    f"the weight of the reference determinant is {cis['c0_sq']:.4f} and the CCSD *T*_{{1}} "
    f"diagnostic is {cis['T1_diag']:.4f}. The active-space Hamiltonian (core energy, "
    "one- and two-electron integrals in the canonical active orbitals, orbital irreps) "
    "and the CASCI vector are distributed as cisplatin_cas10_10.npz.",
    align="justify")

# ------------------------------------------------------------------ Note 7
si.heading("Supplementary Note 7. Measurement grouping and shot estimates")
ms_build.MOVED["measurement"](si, citeS)
si.para(
    "For a group *g* with operator Ĥ_{*g*} = Σ_{*j*∈*g*}*h*_{*j*}*P̂*_{*j*}, the variance "
    "of a single-shot estimate is σ_{*g*}^{2} = ⟨Ĥ_{*g*}^{2}⟩ − ⟨Ĥ_{*g*}⟩^{2}. Allocating "
    "*n*_{*g*} ∝ σ_{*g*} shots minimizes the total for a fixed standard error ε, giving "
    "*N* = (Σ_{*g*}σ_{*g*})^{2}/ε^{2}. Without grouping each term is its own group with "
    "σ_{*j*} = |*h*_{*j*}|(1 − ⟨*P̂*_{*j*}⟩^{2})^{1/2}. We evaluated ⟨Ĥ_{*g*}^{2}⟩ "
    "exactly as ‖Ĥ_{*g*}|Ψ⟩‖^{2} on the 2^{*n*}-dimensional state vector. Supplementary "
    "Table S6 lists group and shot counts for all three colouring heuristics. The circuit "
    "that diagonalizes a GC group is a Clifford circuit with *O*(*n*^{2}/log *n*) "
    f"gates{citeS(21)}; its cost and the noise it adds are not included.",
    align="justify")

# ------------------------------------------------------------------ Note 8
si.heading("Supplementary Note 8. Depolarizing error budget and zero-noise extrapolation")
ms_build.MOVED["noise"](si, citeS)
si.para(
    "The two-qubit channel ℰ(ρ) = (1 − *p*)ρ + (*p*/15)Σ_{*P*≠*II*}*P*ρ*P* multiplies every "
    "non-identity Pauli expectation value by *f* = 1 − 16*p*/15. If the accumulated "
    "noise is replaced by a global channel with the same contraction, ρ → *F*ρ + "
    "(1 − *F*)𝟙/2^{*n*} with *F* = *f*^{*N*} for *N* CNOT gates, and ⟨Ĥ⟩ → *F**E* + "
    "(1 − *F*)*c*_{0} with *c*_{0} = Tr(Ĥ)/2^{*n*}. The bias is (1 − *F*)(*c*_{0} − *E*), "
    "which for small *p* is ≈ (16/15)*pN*(*c*_{0} − *E*). Under unitary folding the effective "
    "fidelity at scale λ is *F*^{λ}; Richardson extrapolation from λ = 1, 3, 5 uses "
    "weights γ = (15/8, −5/4, 3/8), cancels the first and second orders in ln *F* and "
    "leaves a bias |1 − Σγ_{λ}*F*^{λ}|(*c*_{0} − *E*) ≈ (5/2)|ln *F*|^{3}(*c*_{0} − *E*). "
    "Exact rescaling (*E*(1) − *c*_{0})/*F* + *c*_{0} removes the bias but multiplies the "
    "variance by 1/*F*^{2}. The global approximation is optimistic for structured noise "
    "and pessimistic where errors cancel; it is used only to place circuits on an "
    "order-of-magnitude scale. Supplementary Table S8 lists all thresholds.",
    align="justify")

# ------------------------------------------------------------------ Note 9
pes = F["pes"]
si.heading("Supplementary Note 9. Pt–Cl bond-stretching scan")
ms_build.MOVED["pes"](si, citeS)
si.para(
    "One Pt–Cl distance was varied with all other coordinates fixed. All points, "
    "including 2.32 Å, were computed in *C*_{s} (molecular plane retained) so that the "
    "ansatz definition is the same along the scan; at 2.32 Å the *C*_{s} and *C*_{2v} "
    "calculations give the same CASCI energy "
    f"({sci(abs([r for r in pes if abs(r['R_PtCl'] - 2.32) < 1e-6][0]['E_CASCI'] - cis['E_CASCI']))} Ha "
    "difference). The RHF solution at each point was started from the converged density "
    "of the neighbouring point and checked for internal instabilities. Active orbitals "
    "were propagated by maximum overlap (Methods). Supplementary Table S7 gives the "
    "energies, errors and diagnostics.",
    align="justify")

# ------------------------------------------------------------------ Note 10
from facts import orbital_sensitivity
OSL = orbital_sensitivity()
SEL = {"A": "(A) natural orbitals 32–41 of 126 ranked by occupation",
       "B": "(B) canonical MOs 32–41, HOMO = 36",
       "C": "(C) CASSCF converged, *E*_{CASSCF} = "
            + f"{[r for r in OSL if r['choice'][0] == 'C'][0]['E_CASSCF']:.6f} Ha"}
assert [r for r in OSL if r["choice"][0] == "C"][0]["selected"].startswith("CASSCF converged=True")
assert "32-41" in [r for r in OSL if r["choice"][0] == "A"][0]["selected"]
assert "32-41 (HOMO = 36)" in [r for r in OSL if r["choice"][0] == "B"][0]["selected"]
si.heading("Supplementary Note 10. Dependence on the choice of active orbitals")
si.para(
    "All three active spaces use the same *C*_{2v} geometry, basis, ECP and number of "
    "orbitals (ten) and electrons (ten). (A) Symmetry-adapted MP2 natural orbitals, as in "
    "the main text. (B) Canonical RHF orbitals HOMO−4 to LUMO+4. (C) CASSCF(10e,10o) "
    "orbitals, optimized for the ground state with PySCF starting from (A) and expressed as "
    "active natural orbitals. In each space the reference determinant is the closed-shell "
    "determinant of the five most strongly occupied active orbitals, re-optimized inside "
    "the active space (Methods), and VQE starts from it. Supplementary Table S10 lists the "
    "selected orbitals, their character and the results. The canonical virtual orbitals "
    "in (B) are dominated by H s, Pt s and Pt p functions and the canonical occupied "
    "orbitals by Cl p and Pt d functions; the space contains little dynamic or static "
    "correlation. CASSCF (C) concentrates correlation in the active space and gives the "
    "largest correlation energy and the largest single-layer errors.",
    align="justify")

si.page_break()
rows = []
for r in OSL:
    rows.append([r["choice"].split(":")[0], r["choice"].split(": ")[1],
                 f"{abs(r['E_corr_active_mHa']):.2f}", f"{r['c0_sq']:.4f}",
                 f"{r['CCSD_err_mHa']:.3f}", f"{r['CCSD(T)_err_mHa']:.3f}",
                 f"{r['counts']['GVA']['n_params']} / {r['counts']['GVA']['cnot_1L']:,}",
                 f"{r['UCCSD_1L']['dE_mHa']:.4f}", f"{r['GVA_1L']['dE_mHa']:.4f}",
                 f"{r['GVA_2L']['dE_mHa']:.4f}"])
si.table("**Supplementary Table S10 | Cisplatin CAS(10e,10o) at the *C*_{2v} geometry with "
         "three choices of active orbitals.** Correlation energy *E*_{CASCI} − *E*_{ref} "
         "in the active space (mHa, sign omitted); errors relative to CASCI in mHa; VQE from "
         "the reference determinant. UCCSD has 875 parameters and 108,800 CNOT gates per "
         "layer in all three cases; GVA counts depend on the orbital irreps.",
         ["", "Active orbitals", "−*E*_{corr}", "*c*_{0}^{2}", "CCSD", "CCSD(T)",
          "GVA params / CNOT", "UCCSD 1L", "GVA 1L", "GVA 2L"], rows, size=7.5,
         col_widths=[6, 32, 13, 13, 12, 15, 24, 15, 13, 13])
rows = []
for r in OSL:
    for k, (ir, comp) in enumerate(zip(r["active_irreps"], r["active_composition"])):
        rows.append([r["choice"].split(":")[0] if k == 0 else "", k + 1, ir, comp])
si.table("**Supplementary Table S10 (continued) | Character of the active orbitals in each "
         "choice.** Canonical active orbitals after the closed-shell optimization inside the "
         "active space, in order of orbital energy; Mulliken gross populations summed by "
         "element and angular momentum (three largest). Selection details: "
         + "; ".join(SEL[r["choice"][0]] for r in OSL) + ".",
         ["", "Orbital", "Irrep", "Main contributions"], rows, size=7.5,
         col_widths=[8, 16, 16, 120])

# ------------------------------------------------------------------ crystals
import si_periodic
si_periodic.add(si)

# ------------------------------------------------------------------ Tables
si.page_break()
si.heading("Supplementary Tables")
rows = []
for mol, atoms in F["geo"].items():
    name = {"cisplatin_eq": "Cisplatin (*C*_{2v})", "H2": "H_{2}", "H4": "H_{4}",
            "LiH": "LiH", "BeH2": "BeH_{2}"}[mol]
    for k, (a, c) in enumerate(atoms):
        rows.append([name if k == 0 else "", a, f"{c[0]:.6f}", f"{c[1]:.6f}", f"{c[2]:.6f}"])
si.table("**Supplementary Table S1 | Cartesian coordinates (Å) of the reference "
         "geometries.** Stretched small-molecule geometries multiply all coordinates by "
         "1.5 or 2.0. In the Pt–Cl scan the first Cl atom moves along its Pt–Cl axis.",
         ["Molecule", "Atom", "*x*", "*y*", "*z*"], rows, col_widths=[40, 15, 35, 35, 35])

rows = [
    ["Optimizer", "L-BFGS-B (SciPy 1.18), history size 30"],
    ["Gradient", "Exact adjoint (reverse-mode) gradient of the state-vector energy"],
    ["Convergence", "Projected gradient < 10^{−7} or relative energy change < 10^{−12}"],
    ["Iteration limit", "3,000 (small molecules, chains); 2,500 (cisplatin)"],
    ["Starts", "Seed 0: **θ** = 0; seeds *k* ≥ 1: θ ~ N(0, 0.05^{2}), NumPy default_rng(1000 + *k*), (2000 + *k*) or (3000 + *k*) (chains; real and imaginary parts drawn separately)"],
    ["Starts per setting", "Small molecules: 5 (UCCSD, GVA), 1 (ablations); cisplatin: 3 (1L, GVA 2L), 2 (UCCSD 2L), 1 (3L); chains: 3 (UCCSD, *k* only, GVA-k), 1 (*k* + spin, *k* + TR, ADAPT-VQE, and two-layer UCCSD for the 1.0 and 2.0 Å chains with *N*_{*k*} = 5)"],
    ["Generator order", "Singles, same-spin doubles, αβ doubles; tied partners adjacent"],
    ["ADAPT-VQE pool", "GVA pool (screened, tied); full re-optimization after each addition"],
    ["ADAPT stop", "Pool-gradient norm < 10^{−3} Ha, or < 10^{−7} Ha gain over 5 additions, or 120 operators"],
    ["Reference solvers", "PySCF CISD, CCSD, CCSD(T); FCI/CASCI converged to 10^{−12} Ha"],
    ["Chains", "KRHF with Gaussian density fitting, exxdiv = None, conv_tol 10^{−12}; 10 Å vacuum; complex Hamiltonian applied with fci.direct_nosym; two processes, OMP_NUM_THREADS = 1"],
    ["Hardware", "Cloud virtual machine, 2 CPU cores, 8 GB memory; no quantum hardware"],
]
si.table("**Supplementary Table S2 | Optimization and run settings.**", ["Item", "Setting"], rows, col_widths=[40, 120])

rows = []
for r in cis["runs"]:
    rows.append([AN.get(r["ansatz"], r["ansatz"]), r["layers"] if r["layers"] else "–", r["seed"], f"{r['n_params']:,}",
                 f"{r['cnot']:,}", f"{r['E']:.8f}", f"{r['dE'] * 1e3:.5f}",
                 r.get("nit", "–"), sci(r['S2']) if "S2" in r else "–",
                 f"{r['wall_s']:.0f}"])
si.table("**Supplementary Table S3 | All cisplatin CAS(10e,10o) runs at the *C*_{2v} "
         "geometry.** Δ*E* relative to CASCI. For ADAPT-VQE the parameter column gives "
         "the number of operators.",
         ["Ansatz", "Layers", "Seed", "Params", "CNOT", "*E* (Ha)", "Δ*E* (mHa)", "Iter.",
          "⟨Ŝ^{2}⟩", "Wall (s)"], rows,
         col_widths=[22, 12, 10, 14, 17, 28, 17, 12, 17, 11])


def stat(key):
    rr = g[key]
    e = np.array([r["dE"] for r in rr]) * 1e3
    if len(rr) == 1 or np.abs(e).max() < 1e-9:
        return esci(e[0], 2)
    return f"{esci(np.median(e), 2)} [{esci(e.min())}, {esci(e.max())}]"


rows = []
for m in ("H2", "H4", "LiH", "BeH2"):
    for sc in (1.0, 1.5, 2.0):
        for an, L in (("UCCSD", 1), ("GVA", 1), ("PG only", 1), ("spin only", 1),
                      ("UCCSD", 2), ("GVA", 2)):
            rr = g[(m, sc, an, L)]
            rows.append([f"{mt(m)} ×{sc:.1f}", f"{AN.get(an, an)} {L}L", rr[0]["n_params"], f"{rr[0]['cnot']:,}",
                         stat((m, sc, an, L)),
                         f"{int(np.median([r['nit'] for r in rr]))}",
                         "<10^{−12}" if max(r['S2'] for r in rr) < 1e-12 else sci(max(r['S2'] for r in rr))])
        ad = g[(m, sc, "ADAPT-VQE", None)][0]
        rows.append([f"{mt(m)} ×{sc:.1f}", "ADAPT-VQE", ad["n_params"], f"{ad['cnot']:,}",
                     esci(ad['dE'] * 1e3, 2), "–", "–"])
si.table("**Supplementary Table S4 | Small-molecule results (STO-3G).** Δ*E* relative to "
         "FCI in mHa: median [minimum, maximum] over starts. Iterations are medians; "
         "⟨Ŝ^{2}⟩ is the maximum over starts.",
         ["System", "Ansatz", "Params", "CNOT", "Δ*E* (mHa)", "Iter.", "max ⟨Ŝ^{2}⟩"], rows,
         size=7.5, col_widths=[20, 28, 14, 16, 50, 13, 19])

rows = []
for s in F["small_systems"]:
    if s["scale"] != 1.0:
        continue
    c = s["counts"]
    rows.append([mt(s["system"]), s["n_qubits"], c["formula_UCCSD"], c["UCCSD"]["n_params"],
                 c["spin only"]["n_params"], c["PG only"]["n_params"], c["GVA"]["n_params"],
                 f"{c['UCCSD']['cnot_1L']:,}", f"{c['GVA']['cnot_1L']:,}"])
c = cis["counts"]
for r in F["scal"]:
    k = r["n_active_orbitals"]
    if k == 10:
        assert r["GVA_params"] == c["GVA"]["n_params"] and r["GVA_cnot_1L"] == c["GVA"]["cnot_1L"]
    rows.append([f"Cisplatin CAS({k}e,{k}o)", r["n_qubits"], f"{r['UCCSD_formula']:,}",
                 f"{r['UCCSD_params']:,}", f"{r['spinonly_params']:,}",
                 f"{r['PGonly_params']:,}", f"{r['GVA_params']:,}", f"{r['UCCSD_cnot_1L']:,}",
                 f"{r['GVA_cnot_1L']:,}"])
si.table("**Supplementary Table S5 | Parameter and CNOT counts per layer.** 'Formula' is "
         "equation (2); all other parameter counts were obtained by enumeration and agree "
         "with equations (3) and (4) and the screening-only expression. Tying alone does not "
         "change the CNOT count, which equals the UCCSD value; this is also the GVA count for "
         "a molecule without spatial symmetry (*C*_{1}).",
         ["System", "Qubits", "UCCSD formula", "UCCSD", "Tying only", "Screening only", "GVA",
          "CNOT UCCSD", "CNOT GVA"], rows, size=7.5,
         col_widths=[36, 12, 16, 15, 15, 17, 13, 18, 18])

rows = []
for m in F["meas"]:
    for kind in ("QWC", "GC"):
        for meth, lab in (("largest_first", "largest first"), ("dsatur", "DSATUR"),
                          ("sorted_insertion", "sorted insertion")):
            rows.append([mt(m["system"]), f"{m['n_pauli_terms']:,}", kind, lab,
                         f"{m[kind][meth]['n_groups']:,}", sci(m[kind][meth]['shots'], 2)])
    rows.append([mt(m["system"]), f"{m['n_pauli_terms']:,}", "none", "–", f"{m['n_pauli_terms'] - 1:,}",
                 sci(m['shots_ungrouped'], 2)])
si.table("**Supplementary Table S6 | Measurement groups and shots per energy evaluation "
         "(standard error 1.6 mHa, optimal allocation, exact ground state).**",
         ["System", "Pauli terms", "Partition", "Heuristic", "Groups", "Shots"], rows, size=7.5,
         col_widths=[40, 20, 18, 32, 18, 32])

rows = []
for r in pes:
    ef = r["E_CASCI"]
    rows.append([f"{r['R_PtCl']:.2f}", f"{ef:.6f}", f"{(r['E_ref'] - ef) * 1e3:.2f}",
                 f"{(r['E_CISD'] - ef) * 1e3:.3f}", f"{(r['E_CCSD'] - ef) * 1e3:.3f}",
                 f"{(r['E_CCSD(T)'] - ef) * 1e3:.3f}", f"{r['UCCSD_1L']['dE'] * 1e3:.4f}",
                 f"{r['GVA_1L']['dE'] * 1e3:.4f}", f"{r['GVA_2L']['dE'] * 1e3:.4f}",
                 f"{r['c0_sq']:.3f}", f"{r['T1_diag']:.3f}"])
si.table("**Supplementary Table S7 | Pt–Cl scan (*C*_{s}, CAS(10e,10o)).** Errors in mHa "
         "relative to CASCI. In *C*_{s} one layer has "
         f"{pes[0]['counts']['UCCSD']['n_params']} UCCSD and "
         f"{pes[0]['counts']['GVA']['n_params']} GVA parameters.",
         ["*R* (Å)", "*E*_{CASCI} (Ha)", "RHF", "CISD", "CCSD", "CCSD(T)", "UCCSD 1L",
          "GVA 1L", "GVA 2L", "*c*_{0}^{2}", "*T*_{1}"], rows, size=7.5,
         col_widths=[12, 24, 12, 12, 12, 14, 15, 14, 14, 13, 12])

rows = []
for r in F["noise"]:
    rows.append([mt(r["system"]), f"{r['ansatz']} {r['layers']}L", f"{r['cnot']:,}",
                 esci(r['dE_noiseless_mHa']), f"{r['gap_Ha']:.3f}",
                 sci(r['p_max_unmitigated']), sci(r['p_max_richardson']),
                 sci(r['p_max_rescaling_overhead'])])
si.table("**Supplementary Table S8 | Largest tolerable two-qubit depolarizing probability "
         "under the global model (bias ≤ 1.6 mHa, or sampling overhead ≤ 100 for "
         "rescaling).**",
         ["System", "Circuit", "CNOT", "Noiseless Δ*E* (mHa)", "*c*_{0} − *E* (Ha)",
          "Unmitigated", "Richardson", "Rescaling"], rows, size=7.5,
         col_widths=[36, 20, 18, 20, 16, 17, 17, 16])

rows = []
for k, o in enumerate(orb["orbitals"], 1):
    rows.append([k, o["irrep"], f"{o['occupation']:.4f}", o["composition"]])
si.table("**Supplementary Table S9 | Active MP2 natural orbitals of cisplatin (*C*_{2v}).** "
         "Composition: Mulliken gross populations of the normalized orbital summed by "
         "element and angular momentum over symmetry-equivalent atoms (four largest).",
         ["NO", "Irrep", "Occupation", "Main contributions"], rows, size=8,
         col_widths=[12, 16, 22, 110])

si.heading("Supplementary Data 1 and Supplementary Software")
si.para(
    "Supplementary Data 1 (GVA_Source_Data.xlsx) contains the source data for every "
    "figure and table, every optimization run with its final parameters, the geometries, "
    "and live spreadsheet formulas for equations (2), (3) and (5). Supplementary Software "
    "(gva_repro) contains the Python package and the scripts run_small.py, "
    "run_cisplatin.py, run_measure.py, run_scaling.py, run_pes.py, run_periodic.py, "
    "run_scaling_periodic.py, run_measure_periodic.py, uccsd_counts.py, make_figures.py, "
    "make_figures_periodic.py and make_workbook.py, which regenerate all results from "
    "scratch, together with the run logs and full PySCF output.",
    align="justify")
si.heading("Supplementary References")
si.references(SREFS.formatted())
si.save("GVA_npjQI_Supplementary_Information.docx")
print("SI written")
