"""Supplementary Notes 11-15 and Supplementary Tables S11-S15 (crystals).

Imported by si_build.py; adds content to the SI document object passed in.
All numbers are read from results/ through facts.periodic().
"""
import numpy as np

from facts import periodic
from refs import REFS, citeS


def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}×10^{{{int(e)}}}"


def bound(x):
    """Upper bound in powers of ten for small check values."""
    if x <= 0:
        return "10^{−15}"
    e = int(np.ceil(np.log10(x)))
    return f"10^{{{e}}}"


def add(si):
    P = periodic()
    S = P["per"]["systems"]
    stat, tag = P["stat"], P["tag"]
    sp = P["scal"]
    r = lambda k: citeS(k)

    # ------------------------------------------------------------ Note 11
    si.heading("Supplementary Note 11. Crystal momentum and time reversal")
    si.para(
        "Bloch orbitals φ_{*n*​**k**}(**r** + **T**) = e^{*i*​**k**·**T**}φ_{*n*​**k**}(**r**) "
        "sampled on a Γ-centred mesh of *N*_{*k*} points describe a Born–von Kármán "
        "supercell of *N*_{*k*} cells. A lattice translation multiplies *a*^{†}_{*n*​**k**} by "
        "e^{*i*​**k**·**T**}, so a determinant of Bloch orbitals has the crystal momentum "
        "**K** = Σ_{occ} **k** (modulo a reciprocal lattice vector). The Hamiltonian commutes "
        "with every translation; only integrals with **k**_{*p*} − **k**_{*q*} + **k**_{*r*} "
        "− **k**_{*s*} ≡ 0 are non-zero, and the ground state of each system studied here "
        "lies in the sector of the closed-shell reference, **K** = 0 (the FCI weight outside "
        f"it is below {bound(P['max_outside_K'])}). An excitation changes **K** by Δ**k** = "
        "Σ**k**_{created} − Σ**k**_{annihilated}. For Δ**k** ≢ 0 it maps the **K** = 0 sector "
        "onto another sector, so ⟨Ψ|[Ĥ, Ĝ]|Ψ⟩ = 0 whenever |Ψ⟩ was built from the reference "
        "with momentum-conserving generators only. The argument of Supplementary Note 2 "
        "therefore carries over with the translation group in place of the point group: "
        "momentum-violating generators have zero gradient along every trajectory that "
        "starts at **θ** = 0, and removing them does not change the reachable set. In the "
        "benchmark the largest gradient of a momentum-violating excitation at a converged "
        f"two-layer state was {P['max_forbidden_grad']:.1e} Ha.",
        align="justify")
    si.para(
        "The translation group is abelian with *N*_{*k*} one-dimensional irreps labelled by "
        "**k**, so the counting of Supplementary Note 4 applies with one change. In the "
        "abelian point groups used for molecules every irrep is its own inverse; in the "
        "translation group the inverse of **k** is −**k**. Writing *X*(**q**) for the number "
        "of occupied–virtual pairs (*i*, *a*) with **k**_{*a*} − **k**_{*i*} = **q**, the "
        "αβ doubles that conserve momentum number *D*_{αβ} = Σ_{**q**} *X*(**q**)*X*(−**q**), "
        "which reduces to the molecular expression Σ *X*(Γ)^{2} when every element is its "
        "own inverse.",
        align="justify")
    si.para(
        "Time reversal. For a spin-free Hamiltonian, complex conjugation *K̂* in a real "
        "orbital basis commutes with Ĥ. With orbitals chosen so that φ_{*n*,−**k**} = "
        "φ^{*}_{*n*​**k**} and real at time-reversal-invariant **k**, *K̂* maps "
        "*a*^{†}_{*n*​**k**} to *a*^{†}_{*n*,−**k**}. Applied to exp(*t*T̂ − *t*^{*}T̂^{†}) it "
        "gives exp(*t*^{*}T̂′ − *t*T̂′^{†}), where T̂′ is T̂ with every **k** replaced by −**k**; "
        "rewriting T̂′ in the canonical index order of the pool introduces a sign *s* = ±1. "
        "A non-degenerate ground state can be chosen invariant under *K̂*, and an exact "
        "(untrotterized) unitary coupled-cluster state is invariant when *t*(T̂′) = "
        "*s* *t*(T̂)^{*}. GVA-k imposes this relation together with the spin-flip relation "
        "*t*(flip(T̂)) = *t*(T̂). The two relations generate orbits of up to four "
        "excitations. Following the relations around an orbit can force *t* = ±*t*^{*}, in "
        "which case the orbit carries one real parameter; no orbit was forced to zero in "
        "any system or mesh examined. Each spin-flip group then contributes exactly one real "
        "parameter, which gives *N*_{GVA-k} = *S* + *D*_{ss} + (*D*_{αβ} + *n*_{self})/2 "
        "(Methods). As with spin-flip tying, the relations hold exactly only without "
        "Trotterization; their accuracy in the untied optimum is reported in Supplementary "
        "Note 14.",
        align="justify")
    si.para(
        "Complex amplitudes on hardware. For *t* = |*t*|e^{*i*φ} and an excitation that "
        "creates spin orbital *a* (and does not annihilate it), exp(*i*φ*n̂*_{*a*})T̂"
        "exp(−*i*φ*n̂*_{*a*}) = e^{*i*φ}T̂, so exp(*t*T̂ − *t*^{*}T̂^{†}) = exp(*i*φ*n̂*_{*a*})"
        "exp(|*t*|(T̂ − T̂^{†}))exp(−*i*φ*n̂*_{*a*}). Under the Jordan–Wigner mapping "
        "*n̂*_{*a*} = (1 − *Ẑ*_{*a*})/2, so the phase costs two single-qubit *Z* rotations "
        "and no CNOT gates, and the CNOT counts of equation (5) apply unchanged. In the "
        "simulator, each excitation acts on the pair (|*u*⟩, |*d*⟩) with T̂|*u*⟩ = *σ*|*d*⟩ "
        "as |*u*⟩ → cos|*t*||*u*⟩ + *σt* sinc|*t*||*d*⟩, |*d*⟩ → cos|*t*||*d*⟩ − "
        "*σt*^{*} sinc|*t*||*u*⟩ (sinc *x* = sin *x*/*x*), and gradients with respect to Re *t* "
        "and Im *t* are obtained by the same adjoint sweep as for real amplitudes.",
        align="justify")

    # ------------------------------------------------------------ Note 12
    si.heading("Supplementary Note 12. Periodic integrals and their validation")
    si.para(
        "Each cell contains two H atoms on the *x* axis, separated by *r*; neighbouring "
        "units are separated by *d*, so the lattice constant is *a* = *r* + *d*. Parallel "
        "chains are 10 Å apart in *y* and *z*, and the Brillouin zone is sampled only along "
        "the chain. The STO-3G basis gives two bands per cell, one occupied and one "
        f"virtual. Restricted Hartree–Fock with *k*-points{r(25)} used Gaussian density "
        "fitting with the default auxiliary basis and no correction for the exchange "
        "divergence (exxdiv = None), so that the integrals, the KRHF energy and the "
        "supercell energy refer to the same Hamiltonian; finite-size errors of the exchange "
        "energy are therefore not removed, which does not affect comparisons made within one "
        "Hamiltonian. Orbitals at **k** and −**k** were made complex conjugate by solving the "
        "Fock eigenproblem at one point of each pair and conjugating, and orbitals at "
        "time-reversal-invariant points were taken real. The two-electron integrals of the "
        "supercell were assembled from momentum-conserving blocks of the *k*-point "
        "integrals divided by *N*_{*k*}.",
        align="justify")
    si.para(
        "Uniform chains (*r* = *d*) were used only with odd *N*_{*k*}. With two atoms per "
        "cell the bands of a uniform chain are folded and touch at the zone boundary; for "
        "even *N*_{*k*} that point is part of the mesh and the closed-shell reference is "
        "ambiguous.",
        align="justify")
    si.para(
        "For every system we checked (Supplementary Table S11): (i) the Hartree–Fock energy "
        "evaluated from the integrals equals *N*_{*k*} times the KRHF energy per cell, to "
        f"{sci(P['max_hf_diff'])} Ha; (ii) integrals that violate momentum conservation "
        + ("are exactly zero" if P['max_mom_viol'] == 0 else f"vanish (largest {P['max_mom_viol']:.1e})")
        + " and the time-reversal relations hold to "
        + f"{sci(P['max_tr_viol'])}; (iii) three FCI energies agree, one from the Bloch basis "
        "with complex integrals (PySCF direct_nosym, Lanczos), one from the real supercell "
        "basis (direct_spin1) and one from an explicit Γ-point calculation on the "
        "*N*_{*k*}-cell supercell built independently with pyscf.pbc.tools.super_cell; the "
        f"largest differences were {sci(P['max_fci_bloch_diff'])} and "
        f"{sci(P['max_fci_supercell_diff'])} Ha. Adjoint gradients agreed with central "
        "finite differences to 10^{−9} for all four Bloch-basis variants with one and two "
        "layers (tests_periodic.py).",
        align="justify")
    v = P["verify_periodic"]
    if v is not None:
        jw = [x["jw_abs_diff"] for x in v["runs"] if "jw_abs_diff" in x]
        si.para(
            f"All {len(v['runs'])} final states were rebuilt from the stored parameters "
            "(verify_periodic.py); "
            + ("energies and ⟨Ŝ^{2}⟩ agreed with the stored values exactly"
               if max(x['abs_dE'] for x in v['runs']) == 0 and max(x['abs_dS2'] for x in v['runs']) == 0
               else f"energies and ⟨Ŝ^{{2}}⟩ agreed with the stored values to "
                    f"{sci(max(x['abs_dE'] for x in v['runs']))} Ha and "
                    f"{sci(max(x['abs_dS2'] for x in v['runs']))}")
            + ". For the Hartree–Fock-start "
            "single-layer UCCSD and GVA-k states of the systems with up to 16 qubits, the "
            "energy evaluated with the exact Jordan–Wigner qubit Hamiltonian agreed to "
            f"{sci(max(jw))} Ha. Building that operator requires switching off the removal "
            "of terms below 10^{−8} that OpenFermion applies while operators are "
            "accumulated; with that default the uniform chains, which have many integrals "
            "below 10^{−8}, showed differences of about 10^{−9} Ha in testing. For the "
            "20-qubit systems the check made during the runs, with the default operator, "
            f"agreed to {sci(max(r['jw_abs_diff'] for r in P['per']['runs'] if r.get('jw_abs_diff') is not None and r['n_k'] == 5))} Ha.",
            align="justify")
    rows = []
    for s in S:
        rows.append([f"{s['r']:.2f}/{s['d']:.2f}", s["n_k"], s["n_qubits"],
                     f"{s['E_FCI']:.8f}", f"{s['E_corr'] * 1e3:.2f}",
                     f"{s['E_corr'] * 1e3 / s['n_k']:.2f}", f"{s['c0_sq']:.3f}",
                     f"{(s['E_CCSD'] - s['E_FCI']) * 1e3:.3f}" + ("" if s["CCSD_converged"] else "*"),
                     f"{(s['E_CCSD(T)'] - s['E_FCI']) * 1e3:.3f}",
                     sci(abs(s["E_FCI_bloch"] - s["E_FCI"]), 0),
                     sci(abs(s["E_FCI_supercell_gamma"] - s["E_FCI"]), 0)])
    si.table("**Supplementary Table S11 | Periodic H_{2} chains.** *r*/*d* in Å; energies of "
             "the *N*_{*k*}-cell supercell. *E*_{corr} = *E*_{FCI} − *E*_{RHF}. The last two "
             "columns give |*E*_{FCI}(Bloch) − *E*_{FCI}(real)| and |*E*_{FCI}(explicit "
             "supercell) − *E*_{FCI}(real)|. *CCSD not converged in 500 iterations.",
             ["*r*/*d*", "*N*_{*k*}", "Qubits", "*E*_{FCI} (Ha)", "*E*_{corr} (mHa)",
              "per cell", "*c*_{0}^{2}", "CCSD (mHa)", "CCSD(T) (mHa)", "Bloch", "Supercell"],
             rows, size=7.5, col_widths=[16, 9, 11, 25, 17, 13, 12, 17, 18, 14, 16])

    # ------------------------------------------------------------ Note 13
    si.heading("Supplementary Note 13. Counts for crystals")
    one = {x["n_k"]: x for x in sp["one_d"] if x["o"] == 1}
    three = {x["mesh"][0]: x for x in sp["three_d"] if x["o"] == 2}
    mc = {x["mesh"][0]: x for x in sp["molecular_crystal"]}
    si.para(
        "For an insulator with *o* occupied and *v* virtual bands at every **k**, the "
        "quantities in equation (6) are *S* = *N*_{*k*}*ov*, *D*_{αβ} = "
        "*N*_{*k*}^{3}*o*^{2}*v*^{2}, *D*_{ss} = ¼Σ_{**K**}(*N*_{*k*}*o*^{2} − *o*σ_{**K**})"
        "(*N*_{*k*}*v*^{2} − *v*σ_{**K**}) with σ_{**K**} the number of **k** with 2**k** ≡ "
        "**K**, and *n*_{self} = *N*_{*k*}*ov* × 2^{*m*}, where *m* is the number of mesh "
        "directions of even length. For large meshes *N*_{GVA-k} ≈ ¾*N*_{*k*}^{3}*o*^{2}"
        "*v*^{2}, whereas UCCSD on the supercell has *N*_{UCCSD} ≈ "
        "1.5*N*_{*k*}^{4}*o*^{2}*v*^{2}, so the ratio tends to 1/(2*N*_{*k*}). The CNOT "
        "count of the momentum-conserving generators was obtained by grouping occupied and "
        "virtual pairs by total momentum and summing the ladder lengths of equation (5) "
        "within each group; for odd meshes the ratio to supercell UCCSD equals 1/*N*_{*k*} "
        "to within rounding because the index separations are distributed in the same way "
        "in every momentum sector. Explicit enumeration of the pool and of the tied orbits "
        "reproduced all closed-form counts for 12 meshes in one, two and three dimensions, "
        "including meshes with time-reversal-invariant points other than Γ and with *o* ≠ "
        "*v* (Supplementary Table S12). For orientation, a 4×4×4 mesh with eight occupied "
        f"and eight virtual bands per cell ({mc[4]['n_qubits']:,} qubits) gives "
        f"{sci(mc[4]['UCCSD_params'])} UCCSD parameters and {sci(mc[4]['UCCSD_cnot_1L'])} "
        f"CNOT gates per layer, against {sci(mc[4]['GVA_k_params'])} and "
        f"{sci(mc[4]['k_cnot_1L'])} for GVA-k (Supplementary Table S13).",
        align="justify")
    rows = []
    for x in sp["enumeration_check"]:
        rows.append(["×".join(map(str, x["mesh"])), f"{x['o']}+{x['v']}", f"{x['pool']:,}",
                     f"{x['k only']['n_params']:,}", f"{x['k + spin']['n_params']:,}",
                     f"{x['k + TR']['n_params']:,}", f"{x['GVA-k']['n_params']:,}",
                     f"{x['GVA-k']['cnot_1L']:,}", x["GVA-k"]["n_real_orbits"],
                     "yes" if x["agree"] else "NO"])
    si.table("**Supplementary Table S12 | Closed-form counts against enumeration.** Bands "
             "*o*+*v* per cell; pool = supercell UCCSD generators. Parameter columns are real "
             "parameters per layer; 'real' = orbits with one real parameter. 'Agree' means "
             "every enumerated count equals the closed form, including the UCCSD parameter and "
             "CNOT counts.",
             ["Mesh", "Bands", "Pool", "*k* only", "*k* + spin", "*k* + TR", "GVA-k",
              "GVA-k CNOT", "Real", "Agree"], rows, size=7.5,
             col_widths=[16, 12, 17, 15, 16, 16, 15, 19, 11, 12])
    rows = []
    for nk in (2, 3, 4, 5, 8, 16, 32):
        x = one[nk]
        rows.append([f"1D, {nk}", x["n_qubits"], f"{x['UCCSD_params']:,}", f"{x['GVA_k_params']:,}",
                     f"{x['GVA_k_params'] / x['UCCSD_params']:.4f}", f"{x['UCCSD_cnot_1L']:,}",
                     f"{x['k_cnot_1L']:,}", f"{x['k_cnot_1L'] / x['UCCSD_cnot_1L']:.4f}"])
    for n in (2, 3, 4, 5):
        x = three[n]
        rows.append([f"3D, {n}^{{3}}", x["n_qubits"], sci(x["UCCSD_params"]), sci(x["GVA_k_params"]),
                     f"{x['GVA_k_params'] / x['UCCSD_params']:.4f}", sci(x["UCCSD_cnot_1L"]),
                     sci(x["k_cnot_1L"]), f"{x['k_cnot_1L'] / x['UCCSD_cnot_1L']:.4f}"])
    for n in (2, 3, 4):
        x = mc[n]
        rows.append([f"3D, {n}^{{3}}, 8+8", x["n_qubits"], sci(x["UCCSD_params"]),
                     sci(x["GVA_k_params"]), f"{x['GVA_k_params'] / x['UCCSD_params']:.4f}",
                     sci(x["UCCSD_cnot_1L"]), sci(x["k_cnot_1L"]),
                     f"{x['k_cnot_1L'] / x['UCCSD_cnot_1L']:.4f}"])
    si.table("**Supplementary Table S13 | Counts per layer for crystals.** 1D: one occupied "
             "and one virtual band per cell; 3D: two and two unless stated. The full tables "
             "(*N*_{*k*} = 1–32 in 1D, *n* = 1–5 in 3D, all variants) are in Supplementary "
             "Data 1.",
             ["Mesh", "Qubits", "UCCSD params", "GVA-k params", "Ratio", "UCCSD CNOT",
              "GVA-k CNOT", "Ratio"], rows, size=7.5,
             col_widths=[24, 13, 22, 22, 15, 24, 22, 15])

    # ------------------------------------------------------------ Note 14
    si.heading("Supplementary Note 14. Benchmark details for the chains")
    pen = P["tying_penalty_1L_mHa"]
    d1 = P["gva_minus_uccsd_1L_mHa"]
    tc = P["tying_check_periodic"]
    si.para(
        "Every system was run with UCCSD on the supercell in the real orbital basis (one and "
        "two layers), the untied momentum-screened ansatz (*k* only, one layer), the two "
        "partial tyings (*k* + spin, *k* + TR, one layer, Hartree–Fock start), GVA-k (one "
        "and two layers) and ADAPT-VQE with the GVA-k orbits as pool. Seed 0 starts from "
        "zero parameters; seeds 1 and 2 draw every parameter (for complex amplitudes, the "
        "real and imaginary parts) from N(0, 0.05^{2}) with NumPy seeds 3001 and 3002. "
        "For the uniform 1.0 and 2.0 Å chains with *N*_{*k*} = 5, two-layer supercell UCCSD "
        "was run from the Hartree–Fock start only: from random starts at 20 qubits it took "
        "25–33 min per run and stopped at the 3,000-iteration limit. ADAPT-VQE added at "
        "most 120 orbits and, for the systems run last, stopped adding operators after "
        "30 min of wall time; runs that ended this way are marked 'time limit' in "
        "Supplementary Data 1. "
        "Optimizer settings are those of the molecular runs (Supplementary Table S2); the "
        "energy history is taken from evaluations the optimizer already made. "
        f"In total {P['n_runs']} runs were made (Supplementary Table S14 and Supplementary "
        "Data 1).",
        align="justify")
    trc, spc = P["tr_cost_mHa"], P["spin_cost_mHa"]
    lab = lambda x: f"{SYSD[x]['r']:.1f} Å, *N*_{{*k*}} = {SYSD[x]['n_k']}"
    SYSD = P["sys"]
    si.para(
        "Tying cost. In the weakly correlated systems (reference weight above 0.85), "
        "single-layer GVA-k from the Hartree–Fock start differed from single-layer "
        f"supercell UCCSD by at most {P['d1_weak_uHa']:.1f} µHa and from the untied "
        f"*k*-only ansatz by at most {max(abs(trc[x]) + abs(spc[x]) for x in P['weak']) * 1e3:.1f} "
        "µHa. In the strongly correlated systems the time-reversal relation costs energy "
        "at one layer: adding it to the untied ansatz raised the energy by "
        + ", ".join(f"{trc[x]:.2f} mHa ({lab(x)})" for x in P["strong"])
        + ", whereas spin-flip tying alone changed it by "
        + ", ".join(f"{spc[x]:+.2f}" for x in P["strong"]) + " mHa."
        + (" In the untied optimum, the amplitudes of time-reversal partners deviated from "
           "the imposed relation by at most "
           f"{100 * max(x['rel_dev_time_reversal'] for x in tc if x['system'] in P['weak']):.1f}% "
           "of the largest amplitude in the weakly correlated systems and by up to "
           f"{100 * max(x['rel_dev_time_reversal'] for x in tc if x['system'] in P['strong']):.0f}% "
           "in the strongly correlated ones; for spin-flip partners the largest deviations were "
           f"{100 * max(x['rel_dev_spin_flip'] for x in tc if x['system'] in P['weak']):.1f}% and "
           f"{100 * max(x['rel_dev_spin_flip'] for x in tc if x['system'] in P['strong']):.0f}%. "
           "Both relations hold exactly only for the untrotterized exponential, and the "
           "ordering error grows with the size of the amplitudes." if tc else ""),
        align="justify")
    si.para(
        "A small floor. For the smallest system (*N*_{*k*} = 2), every ansatz stopped "
        f"{stat(tag(0.74, 1.50, 2), 'UCCSD', 1)['hf'] * 1e3:.1f} µHa above FCI. From the "
        "Hartree–Fock start the optimizer stays in a parameter subspace fixed by the "
        "symmetry of the starting point and converges to a stationary point there; started "
        "from random parameters with full BFGS, the same ansatz reaches FCI to 10^{−14} Ha. "
        "The floor is three orders of magnitude below chemical accuracy; the optimizer "
        "settings were kept and the values are reported as obtained.",
        align="justify")
    cont, tot, gross = P["contaminated"], P["n_by_ansatz"], P["S2_gross"]
    order = [k for k in ("UCCSD", "k only", "k + spin", "k + TR", "GVA-k", "ADAPT-VQE") if k in tot]
    si.para(
        "Spin contamination. Runs with ⟨Ŝ^{2}⟩ > 0.5 numbered "
        + ", ".join(f"{gross.get(k, 0)} of {tot[k]} ({k})" for k in order)
        + "; the largest values were "
        + ", ".join(f"{P['S2_max'][k]:.2f} ({k})" for k in order)
        + ". Smaller contamination (⟨Ŝ^{2}⟩ between 0.01 and 0.5) appeared only in the "
        "strongly correlated chains, for "
        + ", ".join(k for k in order if 0.01 < P["S2_max"][k] <= 0.5 or
                    any(0.01 < r["S2"] <= 0.5 for r in P["per"]["runs"] if r["ansatz"] == k))
        + ". For GVA-k it arises because spin-flip tying keeps the α↔β exchange symmetry only "
        "up to ordering effects between non-commuting partners (Supplementary Note 3), and "
        "those effects grow with the amplitudes "
        "(Supplementary Fig. 2b). One two-layer GVA-k run from the zero start (1.0 Å, "
        "*N*_{*k*} = 3) stopped at the single-layer energy, at the symmetric stationary "
        "point described below; the two random starts did not.",
        align="justify")
    si.figure("figures/FigS2.png", "Supplementary Fig. 2 | Tying ablation for the chains.",
              "**a**, Energy of single-layer ansätze started from the Hartree–Fock state, "
              "relative to the untied momentum-screened ansatz (*k* only); symmetric "
              "logarithmic axis, linear between ±10^{−3} mHa. UCCSD acts on the real "
              "supercell basis; the other ansätze act on the Bloch basis with momentum "
              "screening and the indicated tying. "
              "**b**, Largest ⟨Ŝ^{2}⟩ over the starts made for each single-layer ansatz.",
              width_mm=165)
    rows = []
    for s in S:
        t = s["system"]
        u1, g1 = stat(t, "UCCSD", 1), stat(t, "GVA-k", 1)
        u2, g2 = stat(t, "UCCSD", 2), stat(t, "GVA-k", 2)
        ad = stat(t, "ADAPT-VQE", None)
        rows.append([f"{s['r']:.2f}/{s['d']:.2f}", s["n_k"],
                     f"{u1['params']}/{g1['params']}", f"{u1['cnot']:,}/{g1['cnot']:,}",
                     f"{u1['hf']:.3f}", f"{g1['hf']:.3f}",
                     f"{u2['min']:.3f}–{u2['max']:.3f}", f"{g2['min']:.3f}–{g2['max']:.3f}",
                     f"{ad['hf']:.3f} ({ad['params']}, {ad['cnot']:,})"])
    si.table("**Supplementary Table S14 | VQE on the chains.** Errors relative to FCI in "
             "mHa. 1L: Hartree–Fock start. 2L: range over three starts (one for UCCSD at 1.0/1.0 "
             "and 2.0/2.0 Å with *N*_{*k*} = 5). Parameters and CNOT "
             "gates per layer, UCCSD/GVA-k. ADAPT: error (real parameters, CNOT gates). "
             "Every run, with ⟨Ŝ^{2}⟩, overlap with the FCI state, iterations and final "
             "parameters, is listed in Supplementary Data 1.",
             ["*r*/*d*", "*N*_{*k*}", "Params", "CNOT", "UCCSD 1L", "GVA-k 1L", "UCCSD 2L",
              "GVA-k 2L", "ADAPT"], rows, size=7.0,
             col_widths=[14, 9, 16, 26, 14, 14, 23, 23, 30])

    # ------------------------------------------------------------ Note 15
    m = P["measurement_periodic"]
    if m:
        si.heading("Supplementary Note 15. Measurement in the Bloch and real bases")
        mm = {(x["system"], x["basis"]): x for x in m}
        t1, t2 = tag(0.74, 1.50, 5), tag(1.0, 1.0, 5)
        si.para(
            "The qubit Hamiltonian depends on the orbital basis. For the 20-qubit chains we "
            "built it in the Bloch basis used by GVA-k and in the real supercell basis used by "
            "UCCSD, removed Pauli terms below 10^{−8}, and grouped the terms by sorted "
            "insertion under qubit-wise (QWC) and general (GC) commutativity "
            "(Supplementary Table S15). The Bloch basis is the more expensive one to measure. "
            "For the dimerized chain it has "
            f"{mm[(t1, 'Bloch')]['n_pauli_terms']:,} terms against "
            f"{mm[(t1, 'real')]['n_pauli_terms']:,} in the real basis and a 1-norm larger by a "
            f"factor of {mm[(t1, 'Bloch')]['one_norm'] / mm[(t1, 'real')]['one_norm']:.2f}, "
            "and GC grouping needs "
            f"{sci(mm[(t1, 'Bloch')]['GC']['shots'])} against "
            f"{sci(mm[(t1, 'real')]['GC']['shots'])} shots per energy evaluation at 1.6 mHa "
            f"(factor {mm[(t1, 'Bloch')]['GC']['shots'] / mm[(t1, 'real')]['GC']['shots']:.2f}; "
            f"{mm[(t2, 'Bloch')]['GC']['shots'] / mm[(t2, 'real')]['GC']['shots']:.2f} for the "
            "1.0 Å chain). With QWC grouping the factor is "
            f"{mm[(t1, 'Bloch')]['QWC']['shots'] / mm[(t1, 'real')]['QWC']['shots']:.1f} and "
            f"{mm[(t2, 'Bloch')]['QWC']['shots'] / mm[(t2, 'real')]['QWC']['shots']:.1f}. Each "
            "complex integral contributes Pauli strings for its real and imaginary parts, which "
            "accounts for the larger number of terms. The penalty is smaller than the factor of "
            "five saved in CNOT gates per layer at *N*_{*k*} = 5, but it does not vanish, and "
            "how it grows with the mesh was not examined here.",
            align="justify")
        rows = []
        for x in m:
            rows.append([f"{x['system'].split('_r')[1].replace('_d', '/').replace('_Nk', ', N_{k}=')}",
                         x["basis"], f"{x['n_pauli_terms']:,}", f"{x['one_norm']:.2f}",
                         f"{x['QWC']['n_groups']:,}", sci(x["QWC"]["shots"]),
                         f"{x['GC']['n_groups']:,}", sci(x["GC"]["shots"])])
        si.table("**Supplementary Table S15 | Measurement cost of the 20-qubit chains.** Pauli "
                 "terms after removing coefficients below 10^{−8} (including the identity); "
                 "1-norm of the non-identity coefficients in Ha; groups and shots for a 1.6 mHa "
                 "standard error with sorted insertion and optimal allocation, evaluated with "
                 "the exact ground state.",
                 ["System", "Basis", "Terms", "1-norm", "QWC groups", "QWC shots", "GC groups",
                  "GC shots"], rows, size=7.5, col_widths=[34, 14, 15, 14, 17, 18, 17, 18])
    return P
