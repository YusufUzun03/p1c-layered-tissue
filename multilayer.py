"""Plane wave at normal incidence on N planar layers: exact transmission-line solution.

Geometry follows the Feko set-up: air for z > 0, the plane wave comes from +z,
tissue occupies z < 0. Internally we use depth x = -z (x = 0 at the skin surface).
Convention: e^{jwt}, peak phasors, so time-averaged quantities carry a factor 1/2.

Each layer is a line section with wave impedance eta_i and propagation constant
k_i = beta_i - j alpha_i. In layer i (local depth u from its top):
    E(u) = a_i e^{-j k u} + a_i Gamma_i e^{-j k (2 d - u)}
    H(u) = [a_i e^{-j k u} - a_i Gamma_i e^{-j k (2 d - u)}] / eta_i
where Gamma_i is the reflection coefficient at the layer's lower interface.
A finite last layer is backed by an air half-space (as in a Feko slab);
a layer with d = None is a half-space.
"""
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from tissue import EPS0, ETA0, C0, DENSITY, complex_eps, gabriel_eps, wave_params


@dataclass
class Layer:
    name: str
    d: Optional[float]            # thickness [m]; None = half-space
    eps_c: complex                # complex relative permittivity
    f: float                      # frequency [Hz] the eps_c belongs to
    density: Optional[float] = None
    # filled in by solve()
    x_top: float = 0.0
    k: complex = 0j
    eta: complex = 0j
    gamma_bottom: complex = 0j    # reflection at lower interface
    a: complex = 0j               # forward amplitude at top of layer
    P_abs: float = 0.0            # absorbed power per area [W/m^2] (Poynting drop)
    P_abs_int: float = 0.0        # same, from integrating 1/2 sigma |E|^2

    @property
    def eps_r(self):
        return self.eps_c.real

    @property
    def sigma(self):
        return -self.eps_c.imag * 2 * np.pi * self.f * EPS0

    @property
    def delta(self):
        return wave_params(self.eps_c, self.f)[2]

    @property
    def wavelength(self):
        return wave_params(self.eps_c, self.f)[3]

    @property
    def x_bottom(self):
        return np.inf if self.d is None else self.x_top + self.d


def make_layer(spec, f_hz):
    """Build a Layer from a dict: {name, d_mm, and either tissue or eps_r + sigma}.

    d_mm = null (None) means half-space. 'tissue' uses the Gabriel model at f_hz;
    explicit eps_r / sigma override it (cite your source in the report).
    """
    name = spec.get("name", spec.get("tissue", "layer"))
    d = None if spec.get("d_mm") is None else spec["d_mm"] * 1e-3
    if "eps_r" in spec:
        eps_c = complex_eps(spec["eps_r"], spec.get("sigma", 0.0), f_hz)
    else:
        eps_c = gabriel_eps(spec["tissue"], f_hz)
    tissue = spec.get("tissue", name).lower()
    density = spec.get("density", DENSITY.get(tissue))
    return Layer(name=name, d=d, eps_c=complex(eps_c), f=f_hz, density=density)


@dataclass
class Result:
    f: float
    S_inc: float
    E0: float
    Gamma: complex
    layers: List[Layer]
    backing: Optional[Layer]
    P_trans: float                # power per area leaving through the bottom [W/m^2]
    notes: List[str] = field(default_factory=list)

    @property
    def P_in(self):
        """Power per area entering the tissue, S_inc (1 - |Gamma|^2)."""
        return self.S_inc * (1 - abs(self.Gamma) ** 2)

    @property
    def conservation_error(self):
        """Relative error of 1 - |Gamma|^2 = sum(P_abs)/S_inc + P_trans/S_inc."""
        total = sum(l.P_abs for l in self.layers) + self.P_trans
        return abs(total - self.P_in) / self.P_in

    @property
    def integral_error(self):
        """Max relative error between Poynting drop and integral of 1/2 sigma |E|^2."""
        errs = [abs(l.P_abs_int - l.P_abs) / l.P_abs for l in self.layers if l.P_abs > 1e-15]
        return max(errs) if errs else 0.0

    # ---- fields ---------------------------------------------------------
    def _all_media(self):
        return self.layers + ([self.backing] if self.backing is not None else [])

    def fields(self, z):
        """Total E and H phasors (peak) and layer index at positions z [m] (z < 0 in tissue).

        Index -1 = air above. At an interface the point belongs to the layer below.
        """
        z = np.atleast_1d(np.asarray(z, dtype=float))
        E = np.zeros(z.shape, complex)
        H = np.zeros(z.shape, complex)
        idx = np.full(z.shape, -1, int)
        x = -z
        k0 = 2 * np.pi * self.f / C0
        above = x < 0
        E[above] = self.E0 * (np.exp(-1j * k0 * x[above]) + self.Gamma * np.exp(1j * k0 * x[above]))
        H[above] = self.E0 * (np.exp(-1j * k0 * x[above]) - self.Gamma * np.exp(1j * k0 * x[above])) / ETA0
        for i, m in enumerate(self._all_media()):
            sel = (x >= m.x_top) & (x < m.x_bottom)
            u = x[sel] - m.x_top
            fwd = m.a * np.exp(-1j * m.k * u)
            if m.d is None:
                bwd = 0
            else:
                bwd = m.a * m.gamma_bottom * np.exp(-1j * m.k * (2 * m.d - u))
            E[sel] = fwd + bwd
            H[sel] = (fwd - bwd) / m.eta
            idx[sel] = i
        return E, H, idx

    def poynting(self, z):
        """Time-averaged power flux density into the tissue (-z direction) [W/m^2]."""
        E, H, _ = self.fields(z)
        return 0.5 * np.real(E * np.conj(H))

    def loss_density(self, z):
        """Absorbed power density 1/2 sigma |E|^2 [W/m^3] (0 in air)."""
        E, _, idx = self.fields(z)
        media = self._all_media()
        sig = np.array([media[i].sigma if i >= 0 else 0.0 for i in idx])
        return 0.5 * sig * np.abs(E) ** 2

    def sar(self, z):
        """Local SAR = sigma |E|^2 / (2 rho) [W/kg] (nan where no density is known)."""
        _, _, idx = self.fields(z)
        media = self._all_media()
        rho = np.array([(media[i].density or np.nan) if i >= 0 else np.nan for i in idx], float)
        return self.loss_density(z) / rho

    def summary(self):
        rows = []
        for l in self.layers:
            rows.append({
                "layer": l.name,
                "d_mm": None if l.d is None else l.d * 1e3,
                "eps_r": l.eps_r, "sigma_Spm": l.sigma,
                "delta_mm": l.delta * 1e3, "lambda_mm": l.wavelength * 1e3,
                "P_abs_Wpm2": l.P_abs,
                "frac_of_incident": l.P_abs / self.S_inc,
                "frac_of_entering": l.P_abs / self.P_in,
            })
        return {
            "f_MHz": self.f / 1e6, "S_inc_Wpm2": self.S_inc, "E0_Vpm_peak": self.E0,
            "Gamma_re": self.Gamma.real, "Gamma_im": self.Gamma.imag, "abs_Gamma": abs(self.Gamma),
            "P_in_Wpm2": self.P_in, "P_trans_Wpm2": self.P_trans,
            "conservation_error": self.conservation_error,
            "integral_error": self.integral_error,
            "layers": rows, "notes": self.notes,
        }


def solve(layers: List[Layer], f_hz: float, S_inc: float = 10.0) -> Result:
    """Exact normal-incidence solution for the given stack (top layer first)."""
    if not layers:
        raise ValueError("need at least one layer")
    for l in layers[:-1]:
        if l.d is None:
            raise ValueError(f"only the last layer may be a half-space ({l.name})")
    layers = list(layers)
    backing = None
    if layers[-1].d is not None:
        backing = Layer("air_below", None, 1 + 0j, f_hz)
    media = layers + ([backing] if backing else [])

    x = 0.0
    for m in media:
        m.k, m.eta, _, _ = wave_params(m.eps_c, f_hz)
        m.x_top = x
        x = x + (m.d or 0.0)

    # Backward sweep: input impedance at the top of each medium
    Zin = [0j] * len(media)
    Zin[-1] = media[-1].eta
    media[-1].gamma_bottom = 0j
    for i in range(len(media) - 2, -1, -1):
        m = media[i]
        m.gamma_bottom = (Zin[i + 1] - m.eta) / (Zin[i + 1] + m.eta)
        g_top = m.gamma_bottom * np.exp(-2j * m.k * m.d)
        Zin[i] = m.eta * (1 + g_top) / (1 - g_top)
    Gamma = (Zin[0] - ETA0) / (Zin[0] + ETA0)

    # Forward sweep: amplitudes from continuity of the total E field
    E0 = np.sqrt(2 * ETA0 * S_inc)
    V = E0 * (1 + Gamma)
    for m in media:
        if m.d is None:
            m.a = V
        else:
            g_top = m.gamma_bottom * np.exp(-2j * m.k * m.d)
            m.a = V / (1 + g_top)
            V = m.a * np.exp(-1j * m.k * m.d) * (1 + m.gamma_bottom)

    res = Result(f=f_hz, S_inc=S_inc, E0=E0, Gamma=complex(Gamma),
                 layers=layers, backing=backing, P_trans=0.0)

    # Absorbed power per layer: drop in Poynting flux, cross-checked by integration
    for l in layers:
        s_top = res.poynting(-(l.x_top + 1e-12))[0]
        if l.d is None:
            s_bot = 0.0
            if not np.isfinite(l.delta):  # lossless half-space: nothing absorbed
                l.P_abs, l.P_abs_int = s_top - s_bot, 0.0
                continue
            u = np.linspace(0, 30 * l.delta, 30001)  # e^-60 of the power is left
        else:
            s_bot = res.poynting(-(l.x_bottom - 1e-12))[0]
            u = np.linspace(0, l.d, max(2001, int(l.d / (l.delta / 200)) + 1))
        l.P_abs = s_top - s_bot
        fwd = l.a * np.exp(-1j * l.k * u)
        bwd = 0 if l.d is None else l.a * l.gamma_bottom * np.exp(-1j * l.k * (2 * l.d - u))
        l.P_abs_int = np.trapezoid(0.5 * l.sigma * np.abs(fwd + bwd) ** 2, u)
    if backing is not None:
        res.P_trans = res.poynting(-(backing.x_top + 1e-12))[0]

    for l in layers:
        if l.d is not None and l is layers[-1] and l.d < 2.5 * l.delta:
            res.notes.append(f"{l.name}: {l.d*1e3:.1f} mm < 2.5 delta = {2.5*l.delta*1e3:.1f} mm "
                             f"(back-face reflection affects the profile)")
    return res


def fit_delta(z, absE, z_top, z_bottom, skip=5e-3, skip_bottom=None):
    """Penetration depth from a log-linear fit of |E| between z_top and z_bottom [m].

    Week 1 Demo 3 procedure: skip `skip` metres at each face, fit ln|E| against depth,
    slope = -alpha, delta = 1/alpha. Works on script or Feko data.

    In a finite slab backed by air the wave reflected at the back face (|Gamma| ~ 0.76
    for muscle/air) adds a ripple that biases the fit by about +-2 % at 2.45 GHz with
    a 5 mm skip. Pass skip_bottom (e.g. 1.5 delta) to fit only the part of the slab
    where the back-face wave has decayed.
    """
    z = np.asarray(z, float)
    absE = np.asarray(absE, float)
    hi, lo = max(z_top, z_bottom), min(z_top, z_bottom)
    sel = (z <= hi - skip) & (z >= lo + (skip if skip_bottom is None else skip_bottom))
    if sel.sum() < 3:
        raise ValueError("not enough points in the fit window")
    depth = -z[sel]
    slope, _ = np.polyfit(depth, np.log(absE[sel]), 1)
    return -1.0 / slope
