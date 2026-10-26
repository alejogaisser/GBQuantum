"""quantum1d: motor de referencia de Schrödinger 1D dependiente del tiempo.

Unidades naturales hbar = m = 1. Resuelve i dpsi/dt = [-(1/2) d^2/dx^2 + V(x)] psi con
split-step de Strang sobre una grilla periódica. Ojo: `quantum1d.evolve` es la FUNCIÓN
(el módulo se alcanza con `from quantum1d.evolve import evolve`).
"""

from . import diagnostics, potentials
from .diagnostics import energy, mean_p, mean_x, norm, std_x
from .evolve import evolve
from .grid import Grid, check_resolution
from .potentials import barrier, free, harmonic
from .wavepacket import gaussian_packet

__all__ = [
    "Grid",
    "check_resolution",
    "gaussian_packet",
    "evolve",
    "potentials",
    "diagnostics",
    "free",
    "barrier",
    "harmonic",
    "norm",
    "mean_x",
    "mean_p",
    "std_x",
    "energy",
]
