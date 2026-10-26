"""Grilla periódica 1D y guarda de resolución.

Unidades naturales: hbar = m = 1. El dominio es el anillo [-L/2, L/2): lo que sale por la
derecha entra por la izquierda, y la FFT asume esa periodicidad.
"""

import math
import operator

import numpy as np


class Grid:
    """Discretización uniforme de [-L/2, L/2) con n nodos y sus números de onda de FFT.

    Atributos:
        n: cantidad de nodos (int, >= 2).
        length: largo L del dominio (float, finito y positivo).
        dx: paso espacial L/n.
        x: coordenadas x_j = -L/2 + j*dx (el punto +L/2 no está: equivale a -L/2).
        k: números de onda 2*pi*fftfreq(n, d=dx), en el orden de la FFT; el máximo
            representable (Nyquist) es pi/dx.
    """

    def __init__(self, n, length):
        """Construye la grilla; lanza ValueError si n < 2 o length no es finito y positivo."""
        try:
            n = operator.index(n)
        except TypeError as exc:
            raise ValueError(f"n debe ser un entero, se recibió {n!r}") from exc
        if n < 2:
            raise ValueError(f"n debe ser >= 2, se recibió {n}")
        length = float(length)
        if not (math.isfinite(length) and length > 0.0):
            raise ValueError(f"length debe ser finito y positivo, se recibió {length}")
        self.n = n
        self.length = length
        self.dx = length / n
        self.x = -length / 2 + self.dx * np.arange(n)
        self.k = 2 * np.pi * np.fft.fftfreq(n, d=self.dx)


def check_resolution(grid, k0, sigma):
    """Verifica que un paquete gaussiano (k0, sigma) esté resuelto por la grilla.

    Si el paquete tiene contenido más allá de k_max = pi/dx, ese contenido se pliega a otro k
    (aliasing) y el resultado es incorrecto sin aviso. Condición exigida:
    pi/dx >= 4*(|k0| + 3/sigma). Devuelve None si pasa; ValueError si no pasa o sigma <= 0.
    """
    k0 = float(k0)
    sigma = float(sigma)
    if not (math.isfinite(k0) and math.isfinite(sigma) and sigma > 0.0):
        raise ValueError(f"k0 debe ser finito y sigma finito y positivo (k0={k0}, sigma={sigma})")
    k_max = math.pi / grid.dx
    # Factor 4 y 3/sigma: ver physics-spec sección 8 (|phi|^2 < e^-18 fuera de |k0|+3/sigma
    # y al menos 8 puntos por la longitud de onda más corta relevante).
    required = 4.0 * (abs(k0) + 3.0 / sigma)
    if k_max < required:
        raise ValueError(
            f"grilla insuficiente: pi/dx = {k_max:.6g} < {required:.6g} requerido "
            f"(4*(|k0| + 3/sigma)); aumente n o reduzca length"
        )
    return None
