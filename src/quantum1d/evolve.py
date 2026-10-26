"""Evolución temporal por split-step de Strang. Unidades hbar = m = 1."""

import math
import operator

import numpy as np


def evolve(psi, V, grid, dt, steps):
    """Avanza psi 'steps' pasos de Strang de tamaño dt bajo H = -(1/2) d^2/dx^2 + V(x).

    Un paso: psi <- exp(-i V dt/2) . IFFT[ exp(-i k^2 dt/2) . FFT[ exp(-i V dt/2) psi ] ].
    Cada factor es una fase (unitario) y la FFT conserva la norma por Parseval, así que la
    norma se conserva para cualquier V real y dt, salvo redondeo. El método es simétrico: error
    global O(dt^2), y dt negativo da la inversa exacta. Equivale a evolucionar exactamente un
    Hamiltoniano modificado, por lo que <H> oscila O(dt^2) sin deriva.

    No modifica psi ni V, no renormaliza, no fusiona medios pasos (evolve(s1) luego evolve(s2)
    es idéntico a evolve(s1+s2)) y no verifica la resolución espectral de psi. Devuelve un
    ndarray complex128 nuevo. ValueError si las formas no son (grid.n,), V no es real, dt no
    es un real finito o steps no es un entero >= 0.
    """
    if isinstance(steps, (bool, np.bool_)):
        raise ValueError("steps debe ser un entero >= 0, no un booleano")
    try:
        steps = operator.index(steps)
    except TypeError as exc:
        raise ValueError(f"steps debe ser un entero >= 0, se recibió {steps!r}") from exc
    if steps < 0:
        raise ValueError(f"steps debe ser >= 0, se recibió {steps}")
    try:
        dt = float(dt)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"dt debe ser un real finito, se recibió {dt!r}") from exc
    if not math.isfinite(dt):
        raise ValueError(f"dt debe ser finito, se recibió {dt}")
    psi = np.asarray(psi)
    V = np.asarray(V)
    if psi.shape != (grid.n,):
        raise ValueError(f"psi debe tener forma ({grid.n},), tiene {psi.shape}")
    if V.shape != (grid.n,):
        raise ValueError(f"V debe tener forma ({grid.n},), tiene {V.shape}")
    if not np.isrealobj(V):
        raise ValueError("V debe ser real")

    psi = np.array(psi, dtype=np.complex128, copy=True)
    half_v = np.exp(-0.5j * V * dt)  # medio paso de potencial (fase unitaria)
    kinetic = np.exp(-0.5j * grid.k**2 * dt)  # paso completo de cinética en k (fase unitaria)
    for _ in range(steps):
        psi = half_v * psi
        psi = np.fft.ifft(kinetic * np.fft.fft(psi))
        psi = half_v * psi
    return psi
