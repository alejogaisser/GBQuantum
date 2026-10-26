"""Colisión completa con una barrera: se compara un paquete, no una onda plana."""

import numpy as np
import pytest
from scipy.special import ndtr

from analytic import barrier_transmission, free_center, free_width, gaussian_density, packet_transmission


@pytest.fixture(scope="module")
def collision(q):
    """Una sola trayectoria compartida; sin usar integradores del auditor."""
    grid = q.Grid(8192, 256.0)
    x0, k0, sigma = -35.0, 2.0, 2.0
    center, width, height = grid.dx / 2, 1.0, 3.0
    psi = np.sqrt(gaussian_density(grid.x, x0, sigma)) * np.exp(1j * k0 * grid.x)
    initial = psi.copy()
    potential = q.potentials.barrier(grid, center=center, width=width, height=height)
    reference = packet_transmission(initial, grid.dx, height, width)
    outside = np.abs(grid.x) > 0.4 * grid.length
    edge_masses = [float(np.sum(np.abs(psi[outside]) ** 2) * grid.dx)]
    # dt acota la fase cinética por paso en el Nyquist: dt*(pi/dx)^2/2 = 2.53 rad <= pi.
    # Con dt=0.0025 era 12.6 rad (> 2pi): la barrera discontinua alimenta modos de k alto
    # que resuenan con el paquete, <H> sube ~32% y llega masa al borde (4.5e-5).
    # Mismo tiempo total: 16 controles de 5000 pasos = 2.5 unidades de tiempo cada uno.
    for _ in range(16):
        psi = q.evolve(psi, potential, grid, dt=0.0005, steps=5000)
        edge_masses.append(float(np.sum(np.abs(psi[outside]) ** 2) * grid.dx))
    probability = np.abs(psi) ** 2 * grid.dx
    # El límite es el borde derecho de la barrera; no se integra su interior.
    measured = float(np.sum(probability[grid.x > center + width / 2]))
    return {
        "measured": measured,
        "reference": reference,
        "monochromatic": barrier_transmission(k0**2 / 2, height, width),
        "central_mass": float(np.sum(probability[np.abs(grid.x - center) < 8])),
        "max_edge_mass": max(edge_masses),
        "final_probability": float(np.sum(probability)),
    }


def test_transmision_medida_coincide_con_el_promedio_espectral_del_paquete(collision):
    # 1% absoluto presupuestado antes de correr: barrera espacial discreta y
    # separación temporal con potencial discontinuo; no se exige exactitud espectral.
    discrepancy = abs(collision["measured"] - collision["reference"])
    assert discrepancy < 0.01, (
        f"T medida={collision['measured']:.9g}; T paquete={collision['reference']:.9g}; "
        f"T(E0)={collision['monochromatic']:.9g}; diferencia={discrepancy:.9g}"
    )


def test_medicion_del_tunel_ocurre_con_la_colision_terminada(collision):
    # <0.2% entre x = centro±8 evita medir peso que aún interactúa con la barrera.
    assert collision["central_mass"] < 0.002, collision


def test_colision_no_acumula_masa_apreciable_cerca_del_borde_periodico(collision):
    # <1e−6 en la franja exterior del 20% en controles cada 2.5 unidades de tiempo;
    # muy por debajo del presupuesto de 1% absoluto de transmisión.
    assert collision["max_edge_mass"] < 1e-6, collision


def test_probabilidad_total_se_conserva_en_la_colision_larga(collision):
    # Dieciséis mil pasos acumulan más redondeo que el test de 2000: margen 1e−9.
    assert abs(collision["final_probability"] - 1) < 1e-9, collision


def test_parametros_del_tunel_dejan_las_colas_libres_lejos_del_borde():
    times = np.linspace(0, 40, 401)
    centers = free_center(-35, 2, times)
    sigmas = free_width(2, times)
    safe_boundary = 0.4 * 256
    tails = ndtr((-safe_boundary - centers) / sigmas) + ndtr((centers - safe_boundary) / sigmas)
    # Es una comprobación previa de escala, no una prueba de la dinámica dispersada.
    # El corte exterior y el presupuesto 1e−6 coinciden con el control de trayectoria.
    assert np.max(tails) < 1e-6
