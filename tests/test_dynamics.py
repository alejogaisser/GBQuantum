"""Norma, dispersión y oscilador con predicciones exactas independientes."""

import numpy as np
import pytest

from analytic import free_center, free_width, gaussian_density, harmonic_center


@pytest.mark.parametrize("potential_kind", ["free", "barrier", "harmonic", "smooth_real"])
def test_norma_se_conserva_tras_dos_mil_pasos(q, potential_kind):
    grid = q.Grid(2048, 128.0)
    psi = np.sqrt(gaussian_density(grid.x, -3, 1.4)) * np.exp(1.2j * grid.x)
    if potential_kind == "free":
        potential = q.potentials.free(grid)
    elif potential_kind == "barrier":
        potential = q.potentials.barrier(grid, center=0, width=1, height=3)
    elif potential_kind == "harmonic":
        potential = q.potentials.harmonic(grid, omega=1)
    else:
        potential = 0.3 * np.sin(0.4 * grid.x) + 0.2 * np.cos(0.7 * grid.x)
    evolved = q.evolve(psi, potential, grid, dt=0.004, steps=2000)
    # Los factores de Strang son unitarios; 1e−10 cubre redondeo acumulado FFT.
    assert abs(np.sum(np.abs(evolved) ** 2) * grid.dx - 1) < 1e-10
    assert abs(q.diagnostics.norm(evolved, grid) - 1) < 1e-10


@pytest.mark.parametrize("time", [0.5, 2.0, 6.0])
def test_centro_del_paquete_libre_sigue_la_velocidad_de_grupo(q, time):
    grid = q.Grid(2048, 160.0)
    x0, k0, sigma = -8.0, 1.25, 1.3
    psi = np.sqrt(gaussian_density(grid.x, x0, sigma)) * np.exp(1j * k0 * grid.x)
    evolved = q.evolve(psi, q.potentials.free(grid), grid, dt=0.01, steps=round(time / 0.01))
    # Sin potencial, separación exacta; colas y redondeo quedan muy bajo 2e−8.
    assert abs(q.diagnostics.mean_x(evolved, grid) - free_center(x0, k0, time)) < 2e-8


@pytest.mark.parametrize("time", [0.5, 2.0, 6.0])
def test_ancho_libre_sigue_la_dispersion_exacta(q, time):
    grid = q.Grid(2048, 160.0)
    sigma = 1.3
    psi = np.sqrt(gaussian_density(grid.x, -8, sigma)) * np.exp(1.25j * grid.x)
    evolved = q.evolve(psi, q.potentials.free(grid), grid, dt=0.01, steps=round(time / 0.01))
    # El paquete ocupa <4 unidades de ancho frente a L = 160; no hay envoltura.
    assert abs(q.diagnostics.std_x(evolved, grid) - free_width(sigma, time)) < 2e-8


def test_estado_base_armonico_mantiene_la_densidad(q):
    grid = q.Grid(1024, 32.0)
    sigma = np.sqrt(1 / 2)
    density = gaussian_density(grid.x, 0, sigma)
    psi = np.sqrt(density).astype(complex)
    evolved = q.evolve(psi, q.potentials.harmonic(grid, omega=1), grid, dt=0.002, steps=1000)
    # Strang introduce respiración O(dt²); 2e−6 admite ese error físico numérico.
    assert np.max(np.abs(np.abs(evolved) ** 2 - density)) < 2e-6


def test_centro_armonico_desplazado_sigue_el_coseno_en_varios_instantes(q):
    grid = q.Grid(1024, 32.0)
    x0, omega = 2.0, 1.0
    psi = np.sqrt(gaussian_density(grid.x, x0, np.sqrt(1 / (2 * omega)))).astype(complex)
    potential = q.potentials.harmonic(grid, omega=omega)
    previous_time = 0.0
    for time in [0.5, 1.5, 3.0, 4.5, 6.0]:
        psi = q.evolve(psi, potential, grid, dt=0.002, steps=round((time - previous_time) / 0.002))
        # El desfase de Verlet/Strang crece como t·dt²; cota prefijada a t <= 6.
        assert abs(q.diagnostics.mean_x(psi, grid) - harmonic_center(x0, omega, time)) < 2e-5
        previous_time = time
