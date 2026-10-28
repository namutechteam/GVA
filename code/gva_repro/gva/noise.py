"""Closed-form error budget under a global depolarizing approximation.

Each CNOT carries two-qubit depolarizing noise with probability p (one of the
15 non-identity Pauli pairs). Its polarisation factor is f = 1 - 16p/15, and
after N_CNOT gates the state is approximated by
    rho = F rho_ideal + (1 - F) I / 2^n,   F = f^N_CNOT.
The energy then reads E(p) = F E_var + (1 - F) c0 with c0 = Tr(H)/2^n, which is
the identity coefficient of the qubit Hamiltonian. The model ignores error
structure, coherent errors and single-qubit gate errors; it is used here only
as an order-of-magnitude budget.
"""
import numpy as np
from scipy.optimize import brentq

GAMMA_135 = np.array([15 / 8, -5 / 4, 3 / 8])      # Richardson weights, lambda=1,3,5
LAMBDAS = np.array([1.0, 3.0, 5.0])


def fidelity(p, n_cnot, lam=1.0):
    return (1.0 - 16.0 * p / 15.0) ** (n_cnot * lam)


def bias_unmitigated(p, n_cnot, gap):
    """gap = c0 - E_var (> 0)."""
    return (1.0 - fidelity(p, n_cnot)) * gap


def bias_richardson(p, n_cnot, gap):
    F = np.array([fidelity(p, n_cnot, l) for l in LAMBDAS])
    return abs(1.0 - GAMMA_135 @ F) * gap


def log_rescaling_overhead(p, n_cnot):
    """ln of the sampling overhead 1/F^2 when the global fidelity is known."""
    return -2.0 * n_cnot * np.log1p(-16.0 * p / 15.0)


def p_max(fun, target, lo=1e-14, hi=0.05):
    g = lambda lp: fun(10 ** lp) - target
    a, b = np.log10(lo), np.log10(hi)
    if g(a) > 0:
        return lo
    if g(b) < 0:
        return hi
    return 10 ** brentq(g, a, b, xtol=1e-6)


def budget(n_cnot, gap, eps=1.6e-3, overhead=100.0):
    return {
        "p_max_unmitigated": p_max(lambda p: bias_unmitigated(p, n_cnot, gap), eps),
        "p_max_richardson": p_max(lambda p: bias_richardson(p, n_cnot, gap), eps),
        "p_max_rescaling_overhead": p_max(lambda p: log_rescaling_overhead(p, n_cnot),
                                          np.log(overhead)),
    }
