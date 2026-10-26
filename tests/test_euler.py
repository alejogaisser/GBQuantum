"""Euler explícito vive solo en tests/: demostración negativa de unitariedad."""

import numpy as np


def explicit_euler_eigenmode(psi, energy, dt, steps):
    """Euler sobre Hψ = Eψ: ψ siguiente = ψ − i·dt·Eψ.

    Un modo plano libre es un autoestado exacto. Usar su acción analítica
    evita que modos parásitos FFT de alta frecuencia oculten el crecimiento
    que se quiere medir. No hay renormalización ni acceso al motor.
    """
    state = np.array(psi, dtype=complex, copy=True)
    for _ in range(steps):
        state = state - 1j * dt * energy * state
    return state


def test_euler_explicito_aumenta_la_probabilidad_como_predice_su_factor():
    n, length, mode = 512, 64.0, 10
    dx = length / n
    x = -length / 2 + np.arange(n) * dx
    k = 2 * np.pi * mode / length
    energy = k**2 / 2
    psi = np.exp(1j * k * x) / np.sqrt(length)
    dt, steps = 0.05 / energy, 100
    evolved = explicit_euler_eigenmode(psi, energy, dt, steps)
    probability = np.sum(np.abs(evolved) ** 2) * dx
    expected = (1 + (dt * energy) ** 2) ** steps
    # La predicción es algebraica: ~100 operaciones dejan error relativo <2e−12.
    assert abs(probability / expected - 1) < 2e-12
    # (1+0.05²)^100 ~1.284: crecimiento >20%, lejos del error de redondeo.
    assert probability > 1.2


def test_split_step_conserva_la_probabilidad_en_el_mismo_caso_que_euler(q):
    grid = q.Grid(512, 64.0)
    k = 2 * np.pi * 10 / grid.length
    energy = k**2 / 2
    psi = np.exp(1j * k * grid.x) / np.sqrt(grid.length)
    evolved = q.evolve(psi, np.zeros(grid.n), grid, dt=0.05 / energy, steps=100)
    # El modo libre es exacto con Strang; 1e−10 cubre ampliamente 100 FFT.
    assert abs(np.sum(np.abs(evolved) ** 2) * grid.dx - 1) < 1e-10
