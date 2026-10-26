"""Gabriel (1996) four-term Cole-Cole tissue model.

eps_c(w) = eps_inf + sum_n d_eps_n / (1 + (j w tau_n)^(1 - alpha_n)) + sigma_i / (j w eps0)
         = eps_r - j sigma / (w eps0)

Parameters from S. Gabriel, R. W. Lau, C. Gabriel, Phys. Med. Biol. 41 (1996) 2271-2293
(Part III, Table 1). Values at the course frequencies match the Week 1 table and IT'IS.
"""
import numpy as np

EPS0 = 8.854187817e-12
MU0 = 4e-7 * np.pi
C0 = 1.0 / np.sqrt(MU0 * EPS0)
ETA0 = np.sqrt(MU0 / EPS0)

# name: (eps_inf, sigma_ionic [S/m], [(d_eps, tau [s], alpha) x 4])
GABRIEL = {
    "skin": (4.0, 0.0002, [(32.0, 7.234e-12, 0.0), (1100.0, 32.481e-9, 0.20),
                           (0.0, 159.155e-6, 0.20), (0.0, 15.915e-3, 0.20)]),
    "fat": (2.5, 0.01, [(3.0, 7.958e-12, 0.20), (15.0, 15.915e-9, 0.10),
                        (3.3e4, 159.155e-6, 0.05), (1.0e7, 7.958e-3, 0.01)]),
    "muscle": (4.0, 0.2, [(50.0, 7.234e-12, 0.10), (7000.0, 353.678e-9, 0.10),
                          (1.2e6, 318.31e-6, 0.10), (2.5e7, 2.274e-3, 0.0)]),
    "bone": (2.5, 0.02, [(10.0, 13.263e-12, 0.20), (180.0, 79.577e-9, 0.20),
                         (5000.0, 159.155e-6, 0.20), (1.0e5, 15.915e-3, 0.0)]),
}
GABRIEL["skin_dry"] = GABRIEL["skin"]
GABRIEL["fat_not_infiltrated"] = GABRIEL["fat"]
GABRIEL["bone_cortical"] = GABRIEL["bone"]

# Mass density [kg/m^3] for SAR (muscle: course value; others: IT'IS)
DENSITY = {"skin": 1109.0, "fat": 911.0, "muscle": 1090.0, "bone": 1908.0}


def gabriel_eps(tissue, f_hz):
    """Complex relative permittivity eps_r - j sigma/(w eps0) of a Gabriel tissue."""
    eps_inf, sig_i, terms = GABRIEL[tissue.lower()]
    w = 2 * np.pi * f_hz
    eps = eps_inf + sig_i / (1j * w * EPS0)
    for d_eps, tau, alpha in terms:
        eps = eps + d_eps / (1 + (1j * w * tau) ** (1 - alpha))
    return eps


def gabriel_props(tissue, f_hz):
    """(eps_r, sigma [S/m]) of a Gabriel tissue at f_hz."""
    eps = gabriel_eps(tissue, f_hz)
    return eps.real, -eps.imag * 2 * np.pi * f_hz * EPS0


def complex_eps(eps_r, sigma, f_hz):
    """Complex relative permittivity from eps_r and effective conductivity."""
    return eps_r - 1j * sigma / (2 * np.pi * f_hz * EPS0)


def wave_params(eps_c, f_hz):
    """k [1/m] (k = beta - j alpha), eta [ohm], delta [m] (field 1/e), lambda [m]."""
    n = np.sqrt(eps_c + 0j)
    if n.imag > 0:  # pick the root with a decaying wave for e^{jwt}
        n = -n
    k = 2 * np.pi * f_hz / C0 * n
    eta = ETA0 / n
    alpha = -k.imag
    delta = np.inf if alpha == 0 else 1 / alpha
    return k, eta, delta, 2 * np.pi / k.real
