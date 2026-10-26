"""Potenciales V(x) reales sobre una Grid. Unidades hbar = m = 1 (energía = longitud^-2)."""

import math

import numpy as np


def free(grid):
    """Partícula libre: V = 0 en todos los nodos. Devuelve ndarray float64."""
    return np.zeros(grid.n, dtype=np.float64)


def barrier(grid, center, width, height):
    """Barrera rectangular: V = height en los nodos con |x - center| <= width/2 (cerrado), 0 afuera.

    Discretización: cada nodo representa una celda de ancho dx, así que el ancho efectivo es
    count_nonzero(V)*dx. Si un borde center +- width/2 cae justo sobre un nodo, su inclusión
    la decide el redondeo y el ancho cambia en +-dx (y T ~ exp(-2*kappa*a) cambia mucho):
    elegir center para que ningún borde caiga en un nodo. Lanza ValueError si width <= 0.
    """
    width = float(width)
    if not (math.isfinite(width) and width > 0.0):
        raise ValueError(f"width debe ser finito y positivo, se recibió {width}")
    return np.where(np.abs(grid.x - float(center)) <= width / 2, float(height), 0.0).astype(np.float64)


def harmonic(grid, omega, center=0.0):
    """Oscilador armónico de masa 1: V = (1/2) omega^2 (x - center)^2. Lanza ValueError si omega <= 0.

    En el anillo V tiene un quiebre en +-L/2; es inofensivo si psi ~ 0 ahí.
    """
    omega = float(omega)
    if not (math.isfinite(omega) and omega > 0.0):
        raise ValueError(f"omega debe ser finito y positivo, se recibió {omega}")
    return 0.5 * omega**2 * (grid.x - float(center)) ** 2
