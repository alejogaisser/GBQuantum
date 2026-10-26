"""Diagnósticos de psi sobre una Grid. Unidades hbar = m = 1.

N = dx*sum|psi|^2 es la norma (integral, sin raíz). Los valores esperados se dividen por N,
así que no exigen psi normalizado; si N = 0 lanzan ValueError.
"""

import numpy as np


def norm(psi, grid):
    """Norma N = dx * sum_j |psi_j|^2 (la integral, sin raíz). Invariante: vale 1 en toda la evolución."""
    return float(grid.dx * np.sum(np.abs(psi) ** 2))


def _checked_norm(psi, grid):
    """Devuelve N y lanza ValueError si es 0 (no se puede normalizar un estado nulo)."""
    n = norm(psi, grid)
    if n == 0.0:
        raise ValueError("psi tiene norma 0: los valores esperados no están definidos")
    return n


def mean_x(psi, grid):
    """Posición media <x> = dx * sum_j x_j |psi_j|^2 / N (coordenada de la grilla en [-L/2, L/2))."""
    n = _checked_norm(psi, grid)
    return float(grid.dx * np.sum(grid.x * np.abs(psi) ** 2) / n)


def std_x(psi, grid):
    """Desviación estándar de x, en dos pasadas: sqrt(dx * sum (x_j - <x>)^2 |psi_j|^2 / N)."""
    n = _checked_norm(psi, grid)
    mu = mean_x(psi, grid)
    return float(np.sqrt(grid.dx * np.sum((grid.x - mu) ** 2 * np.abs(psi) ** 2) / n))


def mean_p(psi, grid):
    """Momento medio <p> = sum_m k_m |FFT(psi)_m|^2 / sum_m |FFT(psi)_m|^2 (p = k con hbar = 1)."""
    _checked_norm(psi, grid)
    weights = np.abs(np.fft.fft(psi)) ** 2
    return float(np.sum(grid.k * weights) / np.sum(weights))


def energy(psi, V, grid):
    """Energía <H> = <k^2/2> (cinética, vía FFT) + <V> (potencial en x), ambas divididas por la norma."""
    n = _checked_norm(psi, grid)
    weights = np.abs(np.fft.fft(psi)) ** 2
    kinetic = np.sum(0.5 * grid.k**2 * weights) / np.sum(weights)
    potential = grid.dx * np.sum(np.asarray(V) * np.abs(psi) ** 2) / n
    return float(kinetic + potential)
