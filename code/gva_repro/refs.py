"""Reference list (verified against Crossref/OpenAlex/arXiv, 23 Oct 2026).

Keys are fixed IDs; numbers are assigned in order of first citation.
"""

MASTER = {
    1: "Peruzzo, A. et al. A variational eigenvalue solver on a photonic quantum processor. *Nat. Commun.* **5**, 4213 (2014).",
    2: "McClean, J. R., Romero, J., Babbush, R. & Aspuru-Guzik, A. The theory of variational hybrid quantum-classical algorithms. *New J. Phys.* **18**, 023023 (2016).",
    3: "Cao, Y. et al. Quantum chemistry in the age of quantum computing. *Chem. Rev.* **119**, 10856–10915 (2019).",
    4: "Bauer, B., Bravyi, S., Motta, M. & Chan, G. K.-L. Quantum algorithms for quantum chemistry and quantum materials science. *Chem. Rev.* **120**, 12685–12717 (2020).",
    5: "Romero, J. et al. Strategies for quantum computing molecular energies using the unitary coupled cluster ansatz. *Quantum Sci. Technol.* **4**, 014008 (2019).",
    6: "Grimsley, H. R., Economou, S. E., Barnes, E. & Mayhall, N. J. An adaptive variational algorithm for exact molecular simulations on a quantum computer. *Nat. Commun.* **10**, 3007 (2019).",
    7: "McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R. & Neven, H. Barren plateaus in quantum neural network training landscapes. *Nat. Commun.* **9**, 4812 (2018).",
    8: "Kandala, A. et al. Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets. *Nature* **549**, 242–246 (2017).",
    9: "Anselmetti, G.-L. R., Wierichs, D., Gogolin, C. & Parrish, R. M. Local, expressive, quantum-number-preserving VQE ansätze for fermionic systems. *New J. Phys.* **23**, 113010 (2021).",
    10: "Gard, B. T. et al. Efficient symmetry-preserving state preparation circuits for the variational quantum eigensolver algorithm. *npj Quantum Inf.* **6**, 10 (2020).",
    11: "Cao, C. et al. Progress toward larger molecular simulation on a quantum computer: simulating a system with up to 28 qubits accelerated by point-group symmetry. *Phys. Rev. A* **105**, 062452 (2022).",
    12: "Setia, K. et al. Reducing qubit requirements for quantum simulations using molecular point group symmetries. *J. Chem. Theory Comput.* **16**, 6091–6097 (2020).",
    13: "Bravyi, S., Gambetta, J. M., Mezzacapo, A. & Temme, K. Tapering off qubits to simulate fermionic Hamiltonians. Preprint at https://arxiv.org/abs/1701.08213 (2017).",
    14: "Tsuchimochi, T., Mori, Y. & Ten-no, S. L. Spin-projection for quantum computation: a low-depth approach to strong correlation. *Phys. Rev. Res.* **2**, 043142 (2020).",
    15: "Yordanov, Y. S., Arvidsson-Shukur, D. R. M. & Barnes, C. H. W. Efficient quantum circuits for quantum computational chemistry. *Phys. Rev. A* **102**, 062612 (2020).",
    16: "Verteletskyi, V., Yen, T.-C. & Izmaylov, A. F. Measurement optimization in the variational quantum eigensolver using a minimum clique cover. *J. Chem. Phys.* **152**, 124114 (2020).",
    17: "Gokhale, P. et al. Minimizing state preparations in variational quantum eigensolver by partitioning into commuting families. Preprint at https://arxiv.org/abs/1907.13623 (2019).",
    18: "Yen, T.-C., Verteletskyi, V. & Izmaylov, A. F. Measuring all compatible operators in one series of single-qubit measurements using unitary transformations. *J. Chem. Theory Comput.* **16**, 2400–2409 (2020).",
    19: "Crawford, O. et al. Efficient quantum measurement of Pauli operators in the presence of finite sampling error. *Quantum* **5**, 385 (2021).",
    20: "Wecker, D., Hastings, M. B. & Troyer, M. Progress towards practical quantum variational algorithms. *Phys. Rev. A* **92**, 042303 (2015).",
    21: "Aaronson, S. & Gottesman, D. Improved simulation of stabilizer circuits. *Phys. Rev. A* **70**, 052328 (2004).",
    22: "Temme, K., Bravyi, S. & Gambetta, J. M. Error mitigation for short-depth quantum circuits. *Phys. Rev. Lett.* **119**, 180509 (2017).",
    23: "Li, Y. & Benjamin, S. C. Efficient variational quantum simulator incorporating active error minimization. *Phys. Rev. X* **7**, 021050 (2017).",
    24: "Giurgica-Tiron, T., Hindy, Y., LaRose, R., Mari, A. & Zeng, W. J. Digital zero noise extrapolation for quantum error mitigation. In *2020 IEEE International Conference on Quantum Computing and Engineering (QCE)* 306–316 (IEEE, 2020).",
    25: "Sun, Q. et al. Recent developments in the PySCF program package. *J. Chem. Phys.* **153**, 024109 (2020).",
    26: "McClean, J. R. et al. OpenFermion: the electronic structure package for quantum computers. *Quantum Sci. Technol.* **5**, 034014 (2020).",
    27: "Weigend, F. & Ahlrichs, R. Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for H to Rn: design and assessment of accuracy. *Phys. Chem. Chem. Phys.* **7**, 3297–3305 (2005).",
    28: "Andrae, D., Häußermann, U., Dolg, M., Stoll, H. & Preuß, H. Energy-adjusted ab initio pseudopotentials for the second and third row transition elements. *Theor. Chim. Acta* **77**, 123–141 (1990).",
    29: "Jensen, H. J. Aa., Jørgensen, P., Ågren, H. & Olsen, J. Second-order Møller–Plesset perturbation theory as a configuration and orbital generator in multiconfiguration self-consistent field calculations. *J. Chem. Phys.* **88**, 3834–3839 (1988).",
    30: "Kelland, L. The resurgence of platinum-based cancer chemotherapy. *Nat. Rev. Cancer* **7**, 573–584 (2007).",
    31: "Lee, S. et al. Evaluating the evidence for exponential quantum advantage in ground-state quantum chemistry. *Nat. Commun.* **14**, 1952 (2023).",
    32: "Reiher, M., Wiebe, N., Svore, K. M., Wecker, D. & Troyer, M. Elucidating reaction mechanisms on quantum computers. *Proc. Natl Acad. Sci. USA* **114**, 7555–7560 (2017).",
    33: "Evangelista, F. A., Chan, G. K.-L. & Scuseria, G. E. Exact parameterization of fermionic wave functions via unitary coupled cluster theory. *J. Chem. Phys.* **151**, 244112 (2019).",
    34: "Grimsley, H. R., Claudino, D., Economou, S. E., Barnes, E. & Mayhall, N. J. Is the Trotterized UCCSD ansatz chemically well-defined? *J. Chem. Theory Comput.* **16**, 1–6 (2020).",
    35: "Liu, D. C. & Nocedal, J. On the limited memory BFGS method for large scale optimization. *Math. Program.* **45**, 503–528 (1989).",
    36: "Wierichs, D., Izaac, J., Wang, C. & Lin, C. Y.-Y. General parameter-shift rules for quantum gradients. *Quantum* **6**, 677 (2022).",
    37: "Stilck França, D. & García-Patrón, R. Limitations of optimization algorithms on noisy quantum devices. *Nat. Phys.* **17**, 1221–1227 (2021).",
    38: "Lee, J., Huggins, W. J., Head-Gordon, M. & Whaley, K. B. Generalized unitary coupled cluster wave functions for quantum computation. *J. Chem. Theory Comput.* **15**, 311–324 (2019).",
    39: "Bittel, L. & Kliesch, M. Training variational quantum algorithms is NP-hard. *Phys. Rev. Lett.* **127**, 120502 (2021).",
    40: "Jordan, P. & Wigner, E. Über das Paulische Äquivalenzverbot. *Z. Phys.* **47**, 631–651 (1928).",
    44: "Moses, S. A. et al. A race-track trapped-ion quantum processor. *Phys. Rev. X* **13**, 041052 (2023).",
    42: "Brélaz, D. New methods to color the vertices of a graph. *Commun. ACM* **22**, 251–256 (1979).",
    # added 2026-10-24 (Discussion, scope for C1 fragments); verified via Crossref
    45: "Goings, J. J. et al. Reliably assessing the electronic structure of cytochrome P450 on today's classical computers and tomorrow's quantum computers. *Proc. Natl Acad. Sci. USA* **119**, e2203533119 (2022).",
    46: "Knizia, G. & Chan, G. K.-L. Density matrix embedding: a simple alternative to dynamical mean-field theory. *Phys. Rev. Lett.* **109**, 186404 (2012).",
    47: "Manby, F. R., Stella, M., Goodpaster, J. D. & Miller, T. F. III. A simple, exact density-functional-theory embedding scheme. *J. Chem. Theory Comput.* **8**, 2564–2568 (2012).",
    48: "Riplinger, C. & Neese, F. An efficient and near linear scaling pair natural orbital based local coupled cluster method. *J. Chem. Phys.* **138**, 034106 (2013).",
    49: "Kottmann, J. S., Schleich, P., Tamayo-Mendoza, T. & Aspuru-Guzik, A. Reducing qubit requirements while maintaining numerical precision for the variational quantum eigensolver: a basis-set-free approach. *J. Phys. Chem. Lett.* **12**, 663–673 (2021).",
    # added 2026-10-24 (crystals); verified via Crossref/OpenAlex/arXiv
    50: "Manrique, D. Z., Khan, I. T., Yamamoto, K., Wichitwechkarn, V. & Muñoz Ramo, D. Momentum-space unitary coupled cluster and translational quantum subspace expansion for periodic systems on quantum computers. Preprint at https://arxiv.org/abs/2008.08694 (2020).",
    51: "Yoshioka, N., Sato, T., Nakagawa, Y. O., Ohnishi, Y.-y. & Mizukami, W. Variational quantum simulation for periodic materials. *Phys. Rev. Res.* **4**, 013052 (2022).",
    52: "Motta, M. et al. Towards the solution of the many-electron problem in real materials: equation of state of the hydrogen chain with state-of-the-art many-body methods. *Phys. Rev. X* **7**, 031059 (2017).",
    53: "McClain, J., Sun, Q., Chan, G. K.-L. & Berkelbach, T. C. Gaussian-based coupled-cluster theory for the ground-state and band structure of solids. *J. Chem. Theory Comput.* **13**, 1209–1218 (2017).",
    54: "Sun, Q., Berkelbach, T. C., McClain, J. D. & Chan, G. K.-L. Gaussian and plane-wave mixed density fitting for periodic systems. *J. Chem. Phys.* **147**, 164119 (2017).",
    55: "Beran, G. J. O. Modeling polymorphic molecular crystals with electronic structure theory. *Chem. Rev.* **116**, 5567–5613 (2016).",
}


class Registry:
    def __init__(self):
        self.order = []

    def num(self, key):
        if key not in self.order:
            self.order.append(key)
        return self.order.index(key) + 1

    def index_of(self, key):
        return self.num(key)

    def formatted(self):
        return [MASTER[k] for k in self.order]


REFS = Registry()


def _compress(nums):
    nums = sorted(set(nums))
    out, start, prev = [], nums[0], nums[0]
    for n in nums[1:] + [None]:
        if n is not None and n == prev + 1:
            prev = n
            continue
        if prev - start >= 2:
            out.append(f"{start}–{prev}")
        elif prev - start == 1:
            out += [str(start), str(prev)]
        else:
            out.append(str(start))
        if n is not None:
            start = prev = n
    return ",".join(out)


def cite(*keys):
    return "^{" + _compress([REFS.num(k) for k in keys]) + "}"


# Supplementary Information keeps its own reference list (Nature/npj style), numbered
# in order of first citation in the SI.
SREFS = Registry()


def citeS(*keys):
    return "^{" + _compress([SREFS.num(k) for k in keys]) + "}"
