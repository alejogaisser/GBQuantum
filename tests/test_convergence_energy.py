"""Orden temporal cuantitativo y energía, sin pedir conservación exacta de H."""

import numpy as np
import pytest

from analytic import gaussian_density


def test_strang_converge_con_pendiente_dos_en_norma_compleja(q):
    grid = q.Grid(1024, 32.0)
    psi = np.sqrt(gaussian_density(grid.x, 1.4, np.sqrt(0.5))) * np.exp(0.7j * grid.x)
    potential = q.potentials.harmonic(grid, omega=1)
    final_time = 1.28
    reference_dt = 0.000625
    reference = q.evolve(psi.copy(), potential, grid, dt=reference_dt, steps=2048)
    dts = np.array([0.08, 0.04, 0.02, 0.01])
    errors = []
    for dt in dts:
        evolved = q.evolve(psi.copy(), potential, grid, dt=float(dt), steps=round(final_time / dt))
        errors.append(np.sqrt(grid.dx * np.sum(np.abs(evolved - reference) ** 2)))
    errors = np.array(errors)
    assert np.all(np.isfinite(errors))
    # 1e−9 separa señal de error temporal del piso de redondeo de la referencia.
    assert np.all(errors > 1e-9), f"Errores sin señal temporal resoluble: {errors}"
    assert np.all(np.diff(errors) < 0), f"No hay disminución monótona: {errors}"
    slope, intercept = np.polyfit(np.log(dts), np.log(errors), deg=1)
    fitted = slope * np.log(dts) + intercept
    observed = np.log(errors)
    r_squared = 1 - np.sum((observed - fitted) ** 2) / np.sum((observed - observed.mean()) ** 2)
    # El dt más fino supera 16 veces el de referencia: contaminación <0.4%.
    # [1.8, 2.2] admite términos superiores y excluye orden 1; R² exige ley de potencia.
    assert 1.8 < slope < 2.2, f"Pendiente={slope:.6g}; errores={errors}"
    assert r_squared > 0.995, f"R²={r_squared:.6g}; errores={errors}"


def test_energia_libre_inicial_incluye_la_dispersion_de_momento(q):
    grid = q.Grid(2048, 80.0)
    k0, sigma = 1.7, 1.2
    psi = np.sqrt(gaussian_density(grid.x, -4, sigma)) * np.exp(1j * k0 * grid.x)
    expected = k0**2 / 2 + 1 / (8 * sigma**2)
    # <p²>/2 = k₀²/2 + Δp²/2, Δp = 1/(2σ); FFT resuelta a precisión doble.
    assert abs(q.diagnostics.energy(psi, q.potentials.free(grid), grid) - expected) < 2e-10


@pytest.mark.parametrize("x0", [0.0, 2.0])
def test_energia_armonica_inicial_coincide_con_el_valor_exacto(q, x0):
    grid = q.Grid(1024, 32.0)
    omega = 1.0
    psi = np.sqrt(gaussian_density(grid.x, x0, np.sqrt(1 / (2 * omega)))).astype(complex)
    expected = omega / 2 + omega**2 * x0**2 / 2
    # La base aporta ω/2; el desplazamiento agrega solo energía potencial clásica.
    assert abs(q.diagnostics.energy(psi, q.potentials.harmonic(grid, omega=omega), grid) - expected) < 2e-10


def test_deriva_de_energia_armonica_es_pequena_en_toda_la_trayectoria(q):
    grid = q.Grid(1024, 32.0)
    psi = np.sqrt(gaussian_density(grid.x, 2, np.sqrt(0.5))).astype(complex)
    potential = q.potentials.harmonic(grid, omega=1)
    initial_energy = q.diagnostics.energy(psi, potential, grid)
    deviations = []
    for _ in range(40):
        psi = q.evolve(psi, potential, grid, dt=0.005, steps=50)
        deviations.append(abs(q.diagnostics.energy(psi, potential, grid) - initial_energy) / initial_energy)
    # dt² = 2.5e−5: la oscilación de H para este coherente queda bajo 2e−5.
    assert max(deviations) < 2e-5, f"Máxima deriva relativa: {max(deviations):.6g}"
