"""Contrato público y factores físicos, sin importar módulos internos."""

import numpy as np
import pytest

from analytic import gaussian_density


def test_grilla_es_periodica_y_semiabierta(q):
    grid = q.Grid(256, 32.0)
    assert grid.n == 256
    assert grid.length == 32.0
    assert grid.dx == 32.0 / 256
    expected = -16.0 + np.arange(256) * 0.125
    # Muestras representables exactamente; margen mínimo para aritmética equivalente.
    np.testing.assert_allclose(grid.x, expected, rtol=0, atol=2e-13)
    assert grid.x[-1] < grid.length / 2


def test_numeros_de_onda_siguen_el_orden_fft(q):
    grid = q.Grid(256, 32.0)
    expected = 2 * np.pi * np.fft.fftfreq(grid.n, d=grid.dx)
    # Solo intervienen multiplicaciones y divisiones de precisión doble.
    np.testing.assert_allclose(grid.k, expected, rtol=0, atol=2e-13)


def test_potencial_libre_es_real_y_nulo(q):
    grid = q.Grid(256, 32.0)
    potential = q.potentials.free(grid)
    assert np.isrealobj(potential)
    np.testing.assert_array_equal(potential, np.zeros(grid.n))


def test_barrera_tiene_altura_ancho_y_centro_solicitados(q):
    grid = q.Grid(256, 32.0)
    center, width, height = grid.dx / 2, 2.0, 3.0
    expected = np.where(np.abs(grid.x - center) < width / 2, height, 0)
    potential = q.potentials.barrier(grid, center=center, width=width, height=height)
    assert np.isrealobj(potential)
    # Los bordes caen entre muestras: no se presupone < frente a <=.
    np.testing.assert_array_equal(potential, expected)


def test_potencial_armonico_respeta_masa_y_frecuencia(q):
    grid = q.Grid(256, 32.0)
    omega, center = 0.7, 1.3
    potential = q.potentials.harmonic(grid, omega=omega, center=center)
    assert np.isrealobj(potential)
    # V = mω²(x−c)²/2, m = 1; margen para asociatividad de productos.
    np.testing.assert_allclose(
        potential, 0.5 * omega**2 * (grid.x - center) ** 2, rtol=0, atol=2e-13,
    )


def test_paquete_gaussiano_respeta_la_formula_y_es_complejo(q):
    grid = q.Grid(2048, 80.0)
    x0, k0, sigma = -5.0, 1.75, 1.2
    psi = q.gaussian_packet(grid, x0=x0, k0=k0, sigma=sigma)
    expected = np.sqrt(gaussian_density(grid.x, x0, sigma)) * np.exp(1j * k0 * grid.x)
    assert psi.shape == (grid.n,)
    assert np.iscomplexobj(psi)
    # Colas a >25σ del borde: normalización discreta y continua coinciden.
    np.testing.assert_allclose(psi, expected, rtol=0, atol=2e-12)


def test_paquete_inicial_tiene_probabilidad_unidad(q):
    grid = q.Grid(2048, 80.0)
    psi = q.gaussian_packet(grid, x0=-5, k0=1.75, sigma=1.2)
    # La cuadratura espectral de esta gaussiana resuelta deja solo redondeo.
    assert abs(np.sum(np.abs(psi) ** 2) * grid.dx - 1) < 2e-12


def test_posicion_media_inicial_coincide_con_el_centro_gaussiano(q):
    grid = q.Grid(2048, 80.0)
    psi = np.sqrt(gaussian_density(grid.x, -5, 1.2)).astype(complex)
    # Colas despreciables y cuadratura espectral: margen de redondeo 2e−10.
    assert abs(q.diagnostics.mean_x(psi, grid) + 5) < 2e-10


@pytest.mark.parametrize("k0", [-2.0, 0.0, 2.0])
def test_momento_medio_inicial_respeta_el_signo_de_la_fase(q, k0):
    grid = q.Grid(2048, 80.0)
    x0, sigma = -5.0, 1.2
    psi = np.sqrt(gaussian_density(grid.x, x0, sigma)) * np.exp(1j * k0 * grid.x)
    # Colas despreciables y kmax ~80: error espectral de momentos < 2e−10.
    assert abs(q.diagnostics.mean_p(psi, grid) - k0) < 2e-10


def test_desviacion_inicial_es_sigma_de_la_densidad(q):
    grid = q.Grid(2048, 80.0)
    sigma = 1.2
    psi = np.sqrt(gaussian_density(grid.x, -5, sigma)).astype(complex)
    # Las colas a >25σ y dx ~0.04 hacen despreciable el error de muestreo.
    assert abs(q.diagnostics.std_x(psi, grid) - sigma) < 2e-10


def test_evolucion_no_modifica_la_entrada_y_devuelve_otro_array(q):
    grid = q.Grid(1024, 64.0)
    psi = np.sqrt(gaussian_density(grid.x, -4, 1.2)).astype(complex)
    original = psi.copy()
    evolved = q.evolve(psi, q.potentials.free(grid), grid, dt=0.01, steps=10)
    np.testing.assert_array_equal(psi, original)
    assert evolved.shape == psi.shape
    assert np.iscomplexobj(evolved)
    assert not np.shares_memory(evolved, psi)


def test_modo_plano_adquiere_la_fase_cinetica_y_potencial_exactas(q):
    grid = q.Grid(1024, 64.0)
    k = 2 * np.pi * 7 / grid.length
    potential = np.full(grid.n, 0.7)
    psi = np.exp(1j * k * grid.x) / np.sqrt(grid.length)
    dt, steps = 0.01, 100
    expected = psi * np.exp(-1j * (k**2 / 2 + 0.7) * dt * steps)
    evolved = q.evolve(psi, potential, grid, dt=dt, steps=steps)
    # V constante conmuta con T: no hay error temporal, solo redondeo FFT.
    np.testing.assert_allclose(evolved, expected, rtol=0, atol=2e-12)


def test_evolucion_es_reversible_por_conjugacion_con_potencial_real(q):
    grid = q.Grid(1024, 64.0)
    psi = np.sqrt(gaussian_density(grid.x, -2, 1.2)) * np.exp(0.7j * grid.x)
    potential = q.potentials.harmonic(grid, omega=0.7)
    forward = q.evolve(psi.copy(), potential, grid, dt=0.01, steps=100)
    # Con H real, conjugar permite invertir el tiempo usando solo dt positivo.
    recovered = np.conj(q.evolve(np.conj(forward), potential, grid, dt=0.01, steps=100))
    # Ida y vuelta de un separador simétrico: solo redondeo, sin error O(dt²).
    np.testing.assert_allclose(recovered, psi, rtol=0, atol=2e-11)


def test_guardas_rechazan_una_resolucion_insuficiente(q):
    k0, sigma = 2.0, 1.0
    check = getattr(q, "check_resolution", None)
    if callable(check):
        # kmax = π/(64/128) ~6.28 < 4(2+3) = 20.
        with pytest.raises(ValueError):
            check(q.Grid(128, 64.0), k0, sigma)
    else:
        # Se construye el estado fuera de gaussian_packet para probar evolve/Grid.
        with pytest.raises(ValueError):
            grid = q.Grid(128, 64.0)
            psi = np.sqrt(gaussian_density(grid.x, 0, sigma)) * np.exp(1j * k0 * grid.x)
            q.evolve(psi, np.zeros(grid.n), grid, dt=0.001, steps=1)


def test_helper_de_resolucion_acepta_el_umbral_y_una_grilla_suficiente(q):
    check = getattr(q, "check_resolution", None)
    if callable(check):
        # dx = π/16 da kmax = 16; 4(2 + 3/1.5) = 16 sin cancelaciones.
        grid = q.Grid(256, 16 * np.pi)
        check(grid, k0=2.0, sigma=1.5)
        check(q.Grid(2048, 64.0), k0=2.0, sigma=1.0)
    else:
        # La alternativa contractual solo ofrece aceptación vía evolve.
        grid = q.Grid(2048, 64.0)
        psi = np.sqrt(gaussian_density(grid.x, 0, 1.0)) * np.exp(2j * grid.x)
        evolved = q.evolve(psi, np.zeros(grid.n), grid, dt=0.001, steps=1)
        assert np.all(np.isfinite(evolved))
