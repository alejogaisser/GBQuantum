"""Paquete de ondas gaussiano. Unidades hbar = m = 1."""

import numpy as np

from .grid import check_resolution


def gaussian_packet(grid, x0, k0, sigma):
    """Paquete gaussiano mínimo: psi0 = (2 pi sigma^2)^(-1/4) exp(-(x-x0)^2/(4 sigma^2)) exp(i k0 x).

    x0 es <x>, k0 es <p> (= velocidad) y sigma la desviación estándar de |psi0|^2; en momento
    sigma_k = 1/(2 sigma). Se evalúa la fórmula exacta en los nodos, sin renormalizar: la norma
    discreta es 1 salvo la cola fuera del dominio. Llama a check_resolution (ValueError si la
    grilla no resuelve el paquete o sigma <= 0). Devuelve ndarray complex128.
    """
    check_resolution(grid, k0, sigma)
    x0, k0, sigma = float(x0), float(k0), float(sigma)
    envelope = (2 * np.pi * sigma**2) ** (-0.25) * np.exp(-((grid.x - x0) ** 2) / (4 * sigma**2))
    return (envelope * np.exp(1j * k0 * grid.x)).astype(np.complex128)
