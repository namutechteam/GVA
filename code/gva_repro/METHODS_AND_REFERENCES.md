# Methods, equations and literature behind each part of the code

For every module and function this file states what it computes, the
equation it implements (with the equation numbers used in the manuscript,
"eq. n", and the Supplementary Information, "eq. Sn" / "SI Note n"), the
numerical settings, and the literature the method rests on. Methods that are
our own construction are marked **[this work]**. Literature numbers refer to
the list in Section 9; all entries were checked against Crossref, OpenAlex,
arXiv or the publisher page on 23 October 2026.

Notation: *n*_o occupied and *n*_v virtual spatial orbitals; spin orbital
(*p*, σ) is qubit 2*p* + σ (σ = 0 for α, 1 for β); Ĝ = T̂ − T̂†.

---

## 1. `gva/systems.py` — molecules, reference calculations, active spaces

| Function | What it does | Equations / settings | Basis |
|---|---|---|---|
| `small_molecule_atoms` | H2 0.7414 Å, linear H4 (spacing 1.2300 Å), LiH 1.5950 Å, BeH2 1.3264 Å; all bonds scaled by 1.0/1.5/2.0 | geometry values from the earlier draft | — |
| `cisplatin_atoms` | Idealized *cis*-[Pt(NH3)2Cl2]: Pt–Cl 2.32 Å, Pt–N 2.05 Å, N–H 1.02 Å, Pt–N–H 109.5°, Cl–Pt–Cl = N–Pt–N = 90°; second NH3 is the mirror image of the first | exact C2v (Cs when one Pt–Cl is changed); coordinates in SI Table S1 | **[this work]** |
| `build_small_mol`, `build_cisplatin_mol` | PySCF molecule; STO-3G (small), def2-SVP with the 60-electron def2 ECP on Pt (cisplatin); abelian point group (D2h, C2v, Cs) | 72 explicit electrons, 126 basis functions for cisplatin | [1] PySCF; [2] STO-3G; [3] def2-SVP; [4] Pt ECP |
| `run_rhf` (cisplatin), RHF inside `small_molecule_hamiltonian` | Restricted Hartree–Fock; cisplatin `conv_tol = 1e-10`, second-order solver if not converged, one internal-stability check and re-convergence; small molecules `conv_tol = 1e-12` | standard Roothaan–Hall SCF | [1] |
| `mp2_natural_orbitals` | MP2 one-particle density matrix (all non-ECP electrons) in the MO basis, diagonalized inside each irrep block → symmetry-adapted natural orbitals (NOs) | γ^MP2 = C n C^T per irrep | [5] MP2 NOs as active-space generator |
| `select_active_by_occupation` | the *k* strongly occupied NOs with lowest occupation and *k* weakly occupied NOs with highest occupation | CAS(10e,10o): *k* = 5 | **[this work]**, after [5] |
| `select_active_by_overlap` | along the Pt–Cl scan: NOs with the largest projection Σ\|⟨φ_prev\|φ_new⟩\|² onto the previous active space, using the cross-geometry AO overlap | maximum-overlap tracking | **[this work]** |
| `set_pyscf_log`, `_mole` | optional PySCF output (verbose 4) for every molecule built; used by the v2 runs | — | [1] |
| `cas_hamiltonian_from_orbitals` | frozen-core active-space Hamiltonian | *E*_core = *E*_nuc + Tr[*D*_c(*h* + ½*V*_c)], *h*_eff = *C*_act^T(*h* + *V*_c)*C*_act, (pq\|rs) by AO→MO transformation | standard CASCI construction [1] |
| `_active_rhf_symmetric` | closed-shell SCF inside the active space with rotations restricted to irrep blocks → canonical active orbitals and the reference determinant | *F* = *h* + *J* − ½*K*; tolerance 1e-11 Ha | **[this work]**; reproduces the full RHF energy (difference 5e-12 Ha) |
| `_fake_scf`, `classical_references` | CISD, CCSD, CCSD(T) and FCI/CASCI in the same active space; *c*_0², *T*_1 diagnostic | CISD `conv_tol = 1e-10`; CCSD `conv_tol = 1e-11`, `max_cycle = 500`; FCI Davidson to 1e-12 Ha | [1]; CCSD [6]; (T) [7]; determinant FCI [8]; Davidson [9] |

## 2. `gva/ansatz.py` — generator pool, screening, tying, circuit counts

| Function | What it does | Equations | Basis |
|---|---|---|---|
| `uccsd_pool` | all spin-conserving singles and doubles, ordered singles → same-spin doubles → αβ doubles | eq. 2: *N*_UCCSD = 2*n*_o*n*_v + 2C(*n*_o,2)C(*n*_v,2) + (*n*_o*n*_v)² | UCCSD for VQE [10]; Trotterized product [11, 12] |
| `Excitation.irrep_product`, `gva_groups(use_pg=True)` | keep Ĝ only if Γ_a⊗Γ_b⊗Γ_i⊗Γ_j (Γ_a⊗Γ_i) is totally symmetric (XOR of abelian irrep ids) | SI Note 2, eq. S2 (zero gradient of forbidden generators) | selection rule known in CC theory; prior use in VQE [13]; related qubit tapering [14, 15] |
| `Excitation.spin_flip`, `gva_groups(use_spin=True)` | tie each generator to its α↔β partner; αβ doubles with *i* = *j*, *a* = *b* are self-partnered; partners adjacent in the product | eq. 3 (tying only), eq. 4 (GVA); SI Note 3 (exchange symmetry *c*(*I*_α,*I*_β) = *c*(*I*_β,*I*_α)) | **[this work]**; related spin-complemented pools [16] and spin projection [17] |
| `uccsd_parameter_count`, `counts_summary` | closed-form and enumerated counts for UCCSD, screening only, tying only, GVA | eqs 2–4; screening only 2*S* + 2*D*_ss + *D*_αβ | **[this work]**; checked by `uccsd_counts.py` |
| `excitation_label`, `parameter_map` | human-readable label of each generator (e.g. `D 0a,1b->2a,3b`: annihilate spin orbitals 0α, 1β, create 2α, 3β) and the generator list behind every parameter, in the order used by `Ansatz` | — | **[this work]** (v2) |
| `jw_pauli_weights`, `cnot_count_ladder`, `cnot_count_groups` | Jordan–Wigner Pauli-string weights; CNOT ladder of 2(*w* − 1) per Pauli exponential | eq. 5: single 4(*q* − *p*); double 16(*q* − *p* + *s* − *r* + 1); SI Note 5 | JW [18]; ladder compilation as in [10]; cheaper circuits exist [19] |

## 3. `gva/engine.py` — exact simulator and optimizers

| Function | What it does | Equations / settings | Basis |
|---|---|---|---|
| `DeterminantSpace.__init__` | α and β string lists and single-excitation matrices *E*_ai from PySCF link tables | sign convention of PySCF; *a*†_a*a*†_b*a*_j*a*_i = *E*_ai*E*_bj (same spin), *E*^α_ai*E*^β_bj (αβ); SI Note 1 | determinant-based CI strings [8]; [1] |
| `DeterminantSpace.hamiltonian` | *H*\|ψ⟩ with `pyscf.fci.direct_spin1.contract_2e` (one-electron part absorbed) | — | [8]; [1] |
| `DeterminantSpace.pairs` | determinant pairs (*u*, *d*) with T̂\|*u*⟩ = *s*\|*d*⟩ | — | **[this work]** |
| `Ansatz._rot`, `Ansatz.state` | exp(θĜ) = 1 + sin θ Ĝ + (1 − cos θ)Ĝ², i.e. a Givens rotation on each pair | eq. S1; follows from Ĝ³ = −Ĝ | identity for fermionic excitations, cf. [19] |
| `Ansatz.energy_and_grad` | energy and exact gradient by one backward (adjoint) sweep | ∂*E*/∂θ_k = 2⟨λ_k\|Ĝ_k\|ψ_k⟩; equals the infinite-shot parameter-shift gradient | adjoint differentiation in state-vector simulation [20]; parameter shift [21] |
| `run_vqe` | L-BFGS-B (SciPy), history 30, `gtol = 1e-7`, `ftol = 1e-12`, `maxiter` 3000 (small) / 2500 (cisplatin); records the energy after every iteration; final ⟨Ŝ²⟩ | starts: θ = 0 (seed 0) or θ ~ N(0, 0.05²) | VQE [22, 23]; L-BFGS [24, 25]; SciPy [26] |
| `adapt_vqe` | ADAPT-VQE: add the pool member with the largest \|∂*E*/∂θ\|, re-optimize all parameters (L-BFGS-B, `maxiter = 1500`, `ftol = 1e-12`); stop at gradient norm < 1e-3 Ha, gain < 1e-7 Ha over five additions, or 120 operators | pool gradient 2⟨*H*ψ\|Ĝψ⟩ | ADAPT-VQE [16] |
| `DeterminantSpace.s2` | ⟨Ŝ²⟩ via `pyscf.fci.spin_op.spin_square` | — | [1] |

## 4. `gva/qubit.py` — qubit Hamiltonian, grouping, shots

| Function | What it does | Equations / settings | Basis |
|---|---|---|---|
| `qubit_hamiltonian` | spin-orbital integrals with `spinorb_from_spatial`, `InteractionOperator` with ½·two-body term, Jordan–Wigner transform; Pauli strings stored as symplectic bit masks | qubit *q* on bit *n* − 1 − *q* | OpenFermion [27]; JW [18]; symplectic representation [28] |
| `civec_to_statevector` | PySCF CI vector → JW basis with reordering sign (−1)^{#(i∈α, j∈β, i>j)} | SI Note 1; checked: qubit energy = FCI to < 1e-13 Ha (small), CASCI within the asserted 1e-8 Ha (cisplatin) | **[this work]** |
| `apply_pauli`, `expectation` | *P*\|*b*⟩ = *i*^{n_Y}(−1)^{popcount(z&b)}\|*b*⊕*x*⟩ | — | [28] |
| `conflict_matrix` | QWC conflict (different non-identity Paulis on a shared qubit) and GC conflict (anticommutation) | anticommute ⇔ popcount((x₁&z₂)⊕(z₁&x₂)) odd | QWC/clique cover [29]; general commuting families [30, 31] |
| `greedy_coloring` (largest degree first), `dsatur_coloring`, `sorted_insertion` | three partition heuristics | — | Welsh–Powell [32]; DSATUR [33]; sorted insertion [34] |
| `group_sigmas`, `shots_for_precision` | σ_g² = ⟨Ĥ_g²⟩ − ⟨Ĥ_g⟩² in the exact ground state; optimal allocation *N* = (Σ_g σ_g)²/ε², ε = 1.6 mHa | SI Note 7 | shot allocation [22, 35]; [34] |

The Clifford circuits that diagonalize GC groups (O(*n*²/log *n*) gates [28])
are not counted.

## 5. `gva/noise.py` — depolarizing error budget

| Function | What it does | Equations | Basis |
|---|---|---|---|
| `fidelity` | per-CNOT two-qubit depolarizing channel with probability *p*; contraction *f* = 1 − 16*p*/15; global approximation *F* = *f*^{N_CNOT·λ} | SI Note 8 | depolarizing channel [36]; global-depolarizing approximation [37] |
| `bias_unmitigated` | (1 − *F*)(*c*_0 − *E*), *c*_0 = Tr(Ĥ)/2^n = identity coefficient | Methods | [37] |
| `bias_richardson` | Richardson extrapolation from λ = 1, 3, 5 with weights (15/8, −5/4, 3/8): \|1 − Σγ_λ*F*^λ\|(*c*_0 − *E*) | SI Note 8 | ZNE [38, 39]; unitary folding [40] |
| `log_rescaling_overhead` | ln(1/*F*²), sampling cost of exact fidelity rescaling | — | [37] |
| `p_max`, `budget` | largest *p* with bias ≤ 1.6 mHa (or overhead ≤ 100), Brent root search on log₁₀ *p* | — | SciPy [26]; context [41]; hardware infidelity 1.84×10⁻³ [42] |

## 5a. `gva/periodic.py` — crystals (added in v3)

| Function | What it computes | Equation / method | Literature |
|---|---|---|---|
| `make_cell`, `periodic_hamiltonian` | Chain of H2 units (two H per cell, r, d, a = r + d), STO-3G, 10 Å vacuum between chains, Γ-centred N_k×1×1 mesh; KRHF with Gaussian density fitting and `exxdiv=None`; integrals in a Bloch-orbital basis for the Born–von Kármán supercell | h_pq = ⟨φ_p|ĥ|φ_q⟩ (block-diagonal in k); (pq\|rs) non-zero only if k_p − k_q + k_r − k_s ≡ 0 (mod reciprocal lattice); PySCF k-point ERIs divided by N_k | [44, 45, 1] |
| `_gauge_fixed_orbitals` | Orbitals with φ_{n,−k} = φ*_{n,k}, real at time-reversal-invariant k | generalized eigenproblem F(k)C = S(k)Cε solved for one k of each ±k pair; C(−k) = C(k)* | standard |
| real supercell basis (`U`) | Real orbitals c = (φ_k + φ_−k)/√2, s = (φ_k − φ_−k)/(i√2); same HF determinant | unitary transformation of h and (pq\|rs); imaginary parts checked < 1e-8 | — |
| `supercell_gamma_check` | Independent route: explicit N_k-cell supercell (`pbc.tools.super_cell`), Γ-point RHF and FCI | same Hamiltonian; agreement ~1e-13 Ha | [1] |
| `exc_momentum`, `build_orbits` | Crystal-momentum screening (keep excitations with Σk_created − Σk_annihilated ≡ 0) and tying of spin-flip partners (t equal) and time-reversal partners (t → s t*, s = reordering sign); orbits forced to t = t* or t = −t* carry one real parameter | exact selection rule of the translation group; time reversal = complex conjugation for a spin-free Hamiltonian | [46, 47, 2] |
| `CDeterminantSpace` | Complex Hermitian Hamiltonian on the determinant space (PySCF `fci.direct_nosym`, real and imaginary parts contracted separately) | — | [1] |
| `CAnsatz` | exp(tT̂ − t*T̂†) with complex t as a complex Givens rotation on each determinant pair; adjoint gradients with respect to Re t and Im t | block [[cos\|t\|, −s t* sinc\|t\|],[s t sinc\|t\|, cos\|t\|]]; on hardware the phase of t is a single-qubit Z rotation around a real Givens rotation, so the CNOT count is that of a real excitation | [20] |
| `run_vqe_c`, `adapt_vqe_c` | L-BFGS-B (same settings as molecular runs; energy history taken from cached evaluations) and ADAPT-VQE over tied orbits; `max_wall` stops adding operators after a wall-time limit (1,800 s in the runs; reached for the 20-qubit 1.5 and 2.0 Å chains) | — | [24, 25, 16] |
| `periodic_counts` | Closed forms for o occupied and v virtual bands per cell on any mesh: S = Σ_k o v; D_αβ = Σ_q X(q)X(−q); D_ss = Σ_K P_o(K)P_v(K); n_self = #{(i,a): 2(k_a − k_i) ≡ 0}; N_GVA-k = S + D_ss + (D_αβ + n_self)/2; CNOT sums grouped by pair momentum | generalizes eqs. (3)–(4) from Z₂ⁿ point groups (every element its own inverse) to any abelian group | this work |
| `qubit_hamiltonian_c` | JW image of a complex Hamiltonian. `exact=True`: OpenFermion's internal 1e-8 pruning is disabled during construction (with it, ~1e-9 Ha discrepancies appear for systems with many small integrals); `exact=False`: default pruning, used for the Pauli-term counts | — | [27, 18] |
| `momentum_sector_weights` | Weight of a CI vector in each total-momentum sector | — | — |

Scripts: `run_periodic.py` (benchmark, resumable, two parallel processes merged by
`merge_periodic.py`), `run_scaling_periodic.py` (closed forms + enumeration check on
12 meshes), `run_measure_periodic.py` (Bloch vs real-basis Pauli grouping),
`verify_periodic.py` (rebuild every state from θ, exact JW energies),
`analyze_tying_periodic.py` (untied optimum vs imposed relations),
`tests_periodic.py`, `make_figures_periodic.py` (Figs. 9, 10, Supplementary Fig. 2).

## 6. Run scripts

| Script | Calls | Settings |
|---|---|---|
| `run_small.py` | sections 1–3 for H2, H4, LiH, BeH2 × {1.0, 1.5, 2.0} | UCCSD and GVA with 1 and 2 layers × 5 starts; screening-only and tying-only 1 layer × 1 start; ADAPT-VQE |
| `run_cisplatin.py` | sections 1–3, cisplatin C2v | GVA 1L/2L × 3 starts, 3L × 1; UCCSD 1L × 3, 2L × 2, 3L × 1; ablations; ADAPT-VQE; resumable |
| `run_measure.py` | section 4 | exact ground states; three heuristics × {QWC, GC} |
| `run_scaling.py` | sections 1–2 | CAS(*k*e,*k*o), *k* = 4–20 from the same MP2 NOs |
| `run_pes.py` | sections 1–3, Cs | 9 Pt–Cl distances; tracking; GVA 1L/2L and UCCSD 1L from the Hartree–Fock start |
| `run_orbitals.py` | section 1 | NO irreps, occupations, gross populations by element and angular momentum |
| `log_references.py` | section 1 with PySCF `verbose = 4` | re-runs all reference steps and compares with stored values (tolerance 1e-7 Ha) |
| `uccsd_counts.py` | section 2 | closed forms vs enumeration; non-zero exit on mismatch |
| `verify_theta.py` | sections 1–3 | rebuilds every final state from the stored θ (and, for ADAPT-VQE, the stored operator labels); energy must agree to 1e-9 Ha and ⟨Ŝ²⟩ to 1e-8; non-zero exit on mismatch |
| `compare_v1_v2.py` | — | compares the v2 rerun with the original v1 results run by run |
| `tests_engine.py`, `tests_qubit.py` | sections 3–4 | gradient vs central finite difference (step 1e-6); 2-layer convergence to FCI; JW energy check |

## 7. Figures, workbook, documents

`fig1_schematic.py`, `make_figures.py` (matplotlib; colour-blind-safe palette
[43]), `make_workbook.py` (openpyxl), `ms_build.py`, `si_build.py` and
`docbuilder.py` (python-docx). Every number quoted in the manuscript and SI is
read from `results/` through `facts.py`; several statements are guarded by
`assert` so that the documents cannot be built if the data contradict the text.

## 8. Choices that are ours and could reasonably be made differently

* Idealized cisplatin geometry (no optimization) and def2-SVP basis.
* Active space from MP2 natural orbitals (five strongly and five weakly occupied).
* Generator order singles → same-spin doubles → αβ doubles; Trotterized product
  results depend on the order [12].
* CNOT ladder compilation with all-to-all connectivity.
* Global depolarizing noise model, which is an order-of-magnitude estimate only.
* Optimal shot allocation evaluated with the exact ground state.
* Periodic model: chains of H2 units in STO-3G with k-sampling along the chain only
  and 10 Å vacuum; `exxdiv=None` so that integrals, KRHF and supercell energies are
  consistent (finite-size errors of the exchange energy are not corrected).
* Uniform chains only with odd N_k (closed-shell reference; for even N_k the folded
  bands touch at the zone edge).
* Reference ansatz for crystals: UCCSD on the supercell in a real orbital basis.

## 9. References

1. Sun, Q. et al. Recent developments in the PySCF program package. *J. Chem. Phys.* **153**, 024109 (2020). https://doi.org/10.1063/5.0006074
2. Hehre, W. J., Stewart, R. F. & Pople, J. A. Self-consistent molecular-orbital methods. I. Use of Gaussian expansions of Slater-type atomic orbitals. *J. Chem. Phys.* **51**, 2657–2664 (1969). https://doi.org/10.1063/1.1672392
3. Weigend, F. & Ahlrichs, R. Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for H to Rn: design and assessment of accuracy. *Phys. Chem. Chem. Phys.* **7**, 3297–3305 (2005). https://doi.org/10.1039/b508541a
4. Andrae, D., Häußermann, U., Dolg, M., Stoll, H. & Preuß, H. Energy-adjusted ab initio pseudopotentials for the second and third row transition elements. *Theor. Chim. Acta* **77**, 123–141 (1990). https://doi.org/10.1007/BF01114537
5. Jensen, H. J. Aa., Jørgensen, P., Ågren, H. & Olsen, J. Second-order Møller–Plesset perturbation theory as a configuration and orbital generator in multiconfiguration self-consistent field calculations. *J. Chem. Phys.* **88**, 3834–3839 (1988). https://doi.org/10.1063/1.453884
6. Purvis, G. D. III & Bartlett, R. J. A full coupled-cluster singles and doubles model: the inclusion of disconnected triples. *J. Chem. Phys.* **76**, 1910–1918 (1982). https://doi.org/10.1063/1.443164
7. Raghavachari, K., Trucks, G. W., Pople, J. A. & Head-Gordon, M. A fifth-order perturbation comparison of electron correlation theories. *Chem. Phys. Lett.* **157**, 479–483 (1989). https://doi.org/10.1016/S0009-2614(89)87395-6
8. Knowles, P. J. & Handy, N. C. A new determinant-based full configuration interaction method. *Chem. Phys. Lett.* **111**, 315–321 (1984). https://doi.org/10.1016/0009-2614(84)85513-X
9. Davidson, E. R. The iterative calculation of a few of the lowest eigenvalues and corresponding eigenvectors of large real-symmetric matrices. *J. Comput. Phys.* **17**, 87–94 (1975). https://doi.org/10.1016/0021-9991(75)90065-0
10. Romero, J. et al. Strategies for quantum computing molecular energies using the unitary coupled cluster ansatz. *Quantum Sci. Technol.* **4**, 014008 (2019). https://doi.org/10.1088/2058-9565/aad3e4
11. Evangelista, F. A., Chan, G. K.-L. & Scuseria, G. E. Exact parameterization of fermionic wave functions via unitary coupled cluster theory. *J. Chem. Phys.* **151**, 244112 (2019). https://doi.org/10.1063/1.5133059
12. Grimsley, H. R., Claudino, D., Economou, S. E., Barnes, E. & Mayhall, N. J. Is the Trotterized UCCSD ansatz chemically well-defined? *J. Chem. Theory Comput.* **16**, 1–6 (2020). https://doi.org/10.1021/acs.jctc.9b01083
13. Cao, C. et al. Progress toward larger molecular simulation on a quantum computer: simulating a system with up to 28 qubits accelerated by point-group symmetry. *Phys. Rev. A* **105**, 062452 (2022). https://doi.org/10.1103/PhysRevA.105.062452
14. Setia, K. et al. Reducing qubit requirements for quantum simulations using molecular point group symmetries. *J. Chem. Theory Comput.* **16**, 6091–6097 (2020). https://doi.org/10.1021/acs.jctc.0c00113
15. Bravyi, S., Gambetta, J. M., Mezzacapo, A. & Temme, K. Tapering off qubits to simulate fermionic Hamiltonians. Preprint at https://arxiv.org/abs/1701.08213 (2017).
16. Grimsley, H. R., Economou, S. E., Barnes, E. & Mayhall, N. J. An adaptive variational algorithm for exact molecular simulations on a quantum computer. *Nat. Commun.* **10**, 3007 (2019). https://doi.org/10.1038/s41467-019-10988-2
17. Tsuchimochi, T., Mori, Y. & Ten-no, S. L. Spin-projection for quantum computation: a low-depth approach to strong correlation. *Phys. Rev. Res.* **2**, 043142 (2020). https://doi.org/10.1103/PhysRevResearch.2.043142
18. Jordan, P. & Wigner, E. Über das Paulische Äquivalenzverbot. *Z. Phys.* **47**, 631–651 (1928). https://doi.org/10.1007/BF01331938
19. Yordanov, Y. S., Arvidsson-Shukur, D. R. M. & Barnes, C. H. W. Efficient quantum circuits for quantum computational chemistry. *Phys. Rev. A* **102**, 062612 (2020). https://doi.org/10.1103/PhysRevA.102.062612
20. Jones, T. & Gacon, J. Efficient calculation of gradients in classical simulations of variational quantum algorithms. Preprint at https://arxiv.org/abs/2009.02823 (2020).
21. Wierichs, D., Izaac, J., Wang, C. & Lin, C. Y.-Y. General parameter-shift rules for quantum gradients. *Quantum* **6**, 677 (2022). https://doi.org/10.22331/q-2022-03-30-677
22. McClean, J. R., Romero, J., Babbush, R. & Aspuru-Guzik, A. The theory of variational hybrid quantum-classical algorithms. *New J. Phys.* **18**, 023023 (2016). https://doi.org/10.1088/1367-2630/18/2/023023
23. Peruzzo, A. et al. A variational eigenvalue solver on a photonic quantum processor. *Nat. Commun.* **5**, 4213 (2014). https://doi.org/10.1038/ncomms5213
24. Liu, D. C. & Nocedal, J. On the limited memory BFGS method for large scale optimization. *Math. Program.* **45**, 503–528 (1989). https://doi.org/10.1007/BF01589116
25. Byrd, R. H., Lu, P., Nocedal, J. & Zhu, C. A limited memory algorithm for bound constrained optimization. *SIAM J. Sci. Comput.* **16**, 1190–1208 (1995). https://doi.org/10.1137/0916069
26. Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* **17**, 261–272 (2020). https://doi.org/10.1038/s41592-019-0686-2
27. McClean, J. R. et al. OpenFermion: the electronic structure package for quantum computers. *Quantum Sci. Technol.* **5**, 034014 (2020). https://doi.org/10.1088/2058-9565/ab8ebc
28. Aaronson, S. & Gottesman, D. Improved simulation of stabilizer circuits. *Phys. Rev. A* **70**, 052328 (2004). https://doi.org/10.1103/PhysRevA.70.052328
29. Verteletskyi, V., Yen, T.-C. & Izmaylov, A. F. Measurement optimization in the variational quantum eigensolver using a minimum clique cover. *J. Chem. Phys.* **152**, 124114 (2020). https://doi.org/10.1063/1.5141458
30. Gokhale, P. et al. Minimizing state preparations in variational quantum eigensolver by partitioning into commuting families. Preprint at https://arxiv.org/abs/1907.13623 (2019).
31. Yen, T.-C., Verteletskyi, V. & Izmaylov, A. F. Measuring all compatible operators in one series of single-qubit measurements using unitary transformations. *J. Chem. Theory Comput.* **16**, 2400–2409 (2020). https://doi.org/10.1021/acs.jctc.0c00008
32. Welsh, D. J. A. & Powell, M. B. An upper bound for the chromatic number of a graph and its application to timetabling problems. *Comput. J.* **10**, 85–86 (1967). https://doi.org/10.1093/comjnl/10.1.85
33. Brélaz, D. New methods to color the vertices of a graph. *Commun. ACM* **22**, 251–256 (1979). https://doi.org/10.1145/359094.359101
34. Crawford, O. et al. Efficient quantum measurement of Pauli operators in the presence of finite sampling error. *Quantum* **5**, 385 (2021). https://doi.org/10.22331/q-2021-01-20-385
35. Wecker, D., Hastings, M. B. & Troyer, M. Progress towards practical quantum variational algorithms. *Phys. Rev. A* **92**, 042303 (2015). https://doi.org/10.1103/PhysRevA.92.042303
36. Nielsen, M. A. & Chuang, I. L. *Quantum Computation and Quantum Information: 10th Anniversary Edition* (Cambridge Univ. Press, 2010). https://doi.org/10.1017/CBO9780511976667
37. Vovrosh, J. et al. Simple mitigation of global depolarizing errors in quantum simulations. *Phys. Rev. E* **104**, 035309 (2021). https://doi.org/10.1103/PhysRevE.104.035309
38. Temme, K., Bravyi, S. & Gambetta, J. M. Error mitigation for short-depth quantum circuits. *Phys. Rev. Lett.* **119**, 180509 (2017). https://doi.org/10.1103/PhysRevLett.119.180509
39. Li, Y. & Benjamin, S. C. Efficient variational quantum simulator incorporating active error minimization. *Phys. Rev. X* **7**, 021050 (2017). https://doi.org/10.1103/PhysRevX.7.021050
40. Giurgica-Tiron, T., Hindy, Y., LaRose, R., Mari, A. & Zeng, W. J. Digital zero noise extrapolation for quantum error mitigation. In *2020 IEEE International Conference on Quantum Computing and Engineering (QCE)* 306–316 (IEEE, 2020). https://doi.org/10.1109/QCE49297.2020.00045
41. Stilck França, D. & García-Patrón, R. Limitations of optimization algorithms on noisy quantum devices. *Nat. Phys.* **17**, 1221–1227 (2021). https://doi.org/10.1038/s41567-021-01356-3
42. Moses, S. A. et al. A race-track trapped-ion quantum processor. *Phys. Rev. X* **13**, 041052 (2023). https://doi.org/10.1103/PhysRevX.13.041052
43. Okabe, M. & Ito, K. Color Universal Design (CUD): how to make figures and presentations that are friendly to colorblind people. *J\*FLY* https://jfly.uni-koeln.de/color/ (2002, updated 2008).
44. McClain, J., Sun, Q., Chan, G. K.-L. & Berkelbach, T. C. Gaussian-based coupled-cluster theory for the ground-state and band structure of solids. *J. Chem. Theory Comput.* **13**, 1209–1218 (2017). https://doi.org/10.1021/acs.jctc.7b00049
45. Sun, Q., Berkelbach, T. C., McClain, J. D. & Chan, G. K.-L. Gaussian and plane-wave mixed density fitting for periodic systems. *J. Chem. Phys.* **147**, 164119 (2017). https://doi.org/10.1063/1.4998644
46. Manrique, D. Z., Khan, I. T., Yamamoto, K., Wichitwechkarn, V. & Muñoz Ramo, D. Momentum-space unitary coupled cluster and translational quantum subspace expansion for periodic systems on quantum computers. Preprint at https://arxiv.org/abs/2008.08694 (2020).
47. Yoshioka, N., Sato, T., Nakagawa, Y. O., Ohnishi, Y.-y. & Mizukami, W. Variational quantum simulation for periodic materials. *Phys. Rev. Res.* **4**, 013052 (2022). https://doi.org/10.1103/PhysRevResearch.4.013052
48. Motta, M. et al. Towards the solution of the many-electron problem in real materials: equation of state of the hydrogen chain with state-of-the-art many-body methods. *Phys. Rev. X* **7**, 031059 (2017). https://doi.org/10.1103/PhysRevX.7.031059
49. Beran, G. J. O. Modeling polymorphic molecular crystals with electronic structure theory. *Chem. Rev.* **116**, 5567–5613 (2016). https://doi.org/10.1021/acs.chemrev.5b00648

