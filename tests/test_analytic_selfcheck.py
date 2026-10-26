"""Autoconsistencia ejecutable aun sin pyproject.toml ni quantum1d."""

import numpy as np
import pytest
from scipy.integrate import quad

from analytic import (
    barrier_transmission,
    free_center,
    free_width,
    gaussian_density,
    harmonic_center,
    packet_transmission,
)


def test_ancho_libre_inicial_es_sigma_de_la_densidad():
    assert free_width(1.7, 0) == 1.7


def test_ancho_libre_es_par_en_el_tiempo():
    times = np.array([0.0, 0.5, 3.0, 12.0])
    np.testing.assert_array_equal(free_width(1.2, times), free_width(1.2, -times))


def test_densidad_gaussiana_tiene_masa_y_varianza_correctas():
    sigma = 1.7
    mass = quad(lambda x: gaussian_density(x, 0, sigma), -np.inf, np.inf)[0]
    variance = quad(lambda x: x**2 * gaussian_density(x, 0, sigma), -np.inf, np.inf)[0]
    # La cuadratura adaptativa de estas funciones suaves alcanza ~1e−9 o mejor.
    assert abs(mass - 1) < 2e-9
    assert abs(variance - sigma**2) < 2e-9


def test_centro_libre_respeta_la_velocidad_de_grupo():
    np.testing.assert_array_equal(free_center(-2.0, 3.0, [0, 1, 2]), [-2, 1, 4])


def test_centro_armonico_cierra_un_periodo_y_cambia_de_signo():
    # cos de múltiplos de π tiene redondeo de unas pocas ulps.
    np.testing.assert_allclose(
        harmonic_center(2.0, 1.0, [0, np.pi / 2, np.pi, 2 * np.pi]),
        [2, 0, -2, 2], rtol=0, atol=2e-14,
    )


def test_transmision_de_barrera_nula_es_uno():
    energies = np.array([0.0, 0.1, 3.0, 100.0])
    np.testing.assert_array_equal(barrier_transmission(energies, 0, 2), np.ones(4))


def test_transmision_de_ancho_nulo_es_uno():
    np.testing.assert_array_equal(barrier_transmission([0, 1, 3, 9], 3, 0), np.ones(4))


def test_transmision_a_energia_cero_es_cero():
    assert barrier_transmission(0, 3, 1) == 0


def test_transmision_a_energia_grande_tiende_a_uno():
    # La cota 1−T <= V₀²/[4E(E−V₀)] es < 3e−16 para E = 1e8, V₀ = 3.
    assert abs(barrier_transmission(1e8, 3, 2) - 1) < 2e-15


def test_transmision_es_continua_en_el_umbral():
    height, width = 3.0, 1.2
    expected = 1 / (1 + height * width**2 / 2)
    # Un acercamiento relativo 1e−8 deja un cambio físico O(1e−8), no singular.
    energies = height * np.array([1 - 1e-8, 1, 1 + 1e-8])
    np.testing.assert_allclose(
        barrier_transmission(energies, height, width), expected, rtol=0, atol=2e-8,
    )


def test_transmision_sobre_barrera_tiene_resonancias_unitarias():
    height, width = 3.0, 1.2
    resonant_energy = height + (np.pi / width) ** 2 / 2
    # sin(π) deja solo error de redondeo; el límite físico es T = 1.
    assert abs(barrier_transmission(resonant_energy, height, width) - 1) < 2e-14


def test_transmision_subbarrera_disminuye_con_el_ancho():
    values = [barrier_transmission(1.0, 3.0, width) for width in [0.1, 1, 10]]
    assert 1 > values[0] > values[1] > values[2] >= 0


def test_transmision_es_finita_y_probabilistica_incluso_en_barrera_opaca():
    values = barrier_transmission([0, 1e-300, 0.5, 1, 2, 1e8], 1, 1e4)
    assert np.all(np.isfinite(values))
    assert np.all((values >= 0) & (values <= 1))
    assert values[2] == 0  # Exp(−20000) subdesborda: opacidad física a precisión doble.


def test_transmision_vectorizada_coincide_con_escalar():
    energies = np.array([[0, 0.3, 1], [1.1, 3, 100]])
    expected = np.array([[barrier_transmission(e, 1, 2) for e in row] for row in energies])
    np.testing.assert_array_equal(barrier_transmission(energies, 1, 2), expected)


@pytest.mark.parametrize("energy", [0.4, 1.5, 2.8, 3.2, 4.5, 10.0])
def test_formula_de_barrera_coincide_con_continuidad_de_onda_y_derivada(energy):
    height, width = 3.0, 1.2
    k = np.sqrt(2 * energy)
    internal_k = np.sqrt(2 * (energy - height) + 0j)
    forward = np.exp(1j * internal_k * width)
    backward = np.exp(-1j * internal_k * width)
    transmitted_phase = np.exp(1j * k * width)
    # Incógnitas (r, A, B, t): continuidad de ψ y ψ' en x = 0 y x = a.
    # Es una referencia algebraica independiente de las fórmulas con sinh/sin.
    matching = np.array([
        [1, -1, -1, 0],
        [-k, -internal_k, internal_k, 0],
        [0, forward, backward, -transmitted_phase],
        [0, internal_k * forward, -internal_k * backward, -k * transmitted_phase],
    ], dtype=complex)
    reflection, _, _, transmission = np.linalg.solve(matching, [-1, -k, 0, 0])
    # La matriz es moderadamente condicionada en estos casos; 2e−12 cubre solve.
    assert abs(abs(reflection) ** 2 + abs(transmission) ** 2 - 1) < 2e-12
    assert abs(barrier_transmission(energy, height, width) - abs(transmission) ** 2) < 2e-12


@pytest.mark.parametrize("energy,height,width", [(-1, 1, 1), (1, -1, 1), (1, 1, -1), (np.nan, 1, 1)])
def test_formula_rechaza_parametros_fuera_de_su_dominio(energy, height, width):
    with pytest.raises(ValueError):
        barrier_transmission(energy, height, width)


def test_promedio_fft_coincide_con_integral_continua_independiente():
    n, length = 8192, 256.0
    dx = length / n
    x = -length / 2 + dx * np.arange(n)
    sigma, k0 = 2.0, 2.0
    psi = np.sqrt(gaussian_density(x, -35, sigma)) * np.exp(1j * k0 * x)
    spectral_density = lambda k: np.exp(-2 * sigma**2 * (k - k0) ** 2)
    # epsabs = 1e−11 deja el error de cuadratura muy por debajo de la comparación.
    numerator = quad(
        lambda k: spectral_density(k) * barrier_transmission(k**2 / 2, 3, 1),
        0, np.inf, epsabs=1e-11,
    )[0]
    denominator = quad(spectral_density, 0, np.inf, epsabs=1e-11)[0]
    # Dominio y resolución amplios: discrepancia FFT/cuadratura prevista < 2e−9.
    assert abs(packet_transmission(psi, dx, 3, 1) - numerator / denominator) < 2e-9


def test_promedio_de_un_modo_plano_es_su_transmision():
    n, length, mode = 1024, 80.0, 20
    dx = length / n
    x = -length / 2 + dx * np.arange(n)
    k = 2 * np.pi * mode / length
    psi = np.exp(1j * k * x) / np.sqrt(length)
    # Un único modo FFT reduce el promedio a T(E); margen para redondeo.
    assert abs(packet_transmission(psi, dx, 3, 1) - barrier_transmission(k**2 / 2, 3, 1)) < 2e-13


def test_promedio_rechaza_ausencia_de_peso_espectral_positivo():
    with pytest.raises(ValueError):
        packet_transmission(np.zeros(64, dtype=complex), 0.1, 1, 1)
