"""Experimento 02: Euler explícito vs split-step de Strang (conservación de la norma).

Unidades hbar = m = 1. Euler explícito, psi_{n+1} = psi_n - i dt H psi_n, cumple
||psi_{n+1}||^2 = ||psi_n||^2 + dt^2 ||H psi_n||^2: la norma CRECE en cada paso, para cualquier
dt (incondicionalmente inestable). El split-step, en cambio, es producto de fases y FFT
unitarias: conserva la norma hasta redondeo. Este es el ÚNICO lugar con Euler (decisión 0005).
Genera outputs/02_euler_vs_splitstep.png e imprime tres líneas de resumen.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # backend sin ventana: solo escribe archivos
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from quantum1d import Grid, evolve, gaussian_packet, harmonic, norm  # noqa: E402

# --- Parámetros (physics-spec 11.2) ---
N = 128  # nodos; dx=0.1406, k_max=22.3 cubre la guarda 4(1+4.243)=20.97 a lo largo de la trayectoria
L = 18.0  # largo del anillo; el estado coherente en x0=1 con sigma=0.71 queda lejos de los bordes
OMEGA = 1.0  # frecuencia del oscilador
X0 = 1.0  # desplazamiento inicial del estado coherente
K0 = 0.0  # momento medio inicial
SIGMA = float(np.sqrt(0.5))  # sigma del estado base del oscilador omega=1 (estado coherente)
DT = 2.5e-4  # paso chico: con dt 4 veces mayor los modos altos de Euler explotarían
STEPS = 8000  # t_final = 2
RECORD_EVERY = 20  # pasos entre puntos de la curva graficada
LOG_FLOOR = 1e-17  # piso para log10 de errores nulos (bajo el redondeo ~1e-16); solo para graficar

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "outputs"


def euler_step(psi, V, grid, dt):
    """Un paso de Euler explícito: psi - i dt H psi, con T psi calculado espectralmente (FFT)."""
    h_psi = np.fft.ifft(0.5 * grid.k**2 * np.fft.fft(psi)) + V * psi
    return psi - 1j * dt * h_psi


def main():
    """Evoluciona el mismo estado con Euler y con split-step, grafica |N-1| e imprime resumen."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    grid = Grid(N, L)
    V = harmonic(grid, OMEGA)
    psi0 = gaussian_packet(grid, X0, K0, SIGMA)

    psi_e = psi0
    euler_norms = [norm(psi0, grid)]  # norma de Euler en CADA paso (para verificar monotonía)
    for _ in range(STEPS):
        psi_e = euler_step(psi_e, V, grid, DT)
        euler_norms.append(norm(psi_e, grid))
    euler_norms = np.array(euler_norms)

    psi_s = psi0
    split_norms = [norm(psi0, grid)]
    for _ in range(STEPS // RECORD_EVERY):
        psi_s = evolve(psi_s, V, grid, DT, RECORD_EVERY)
        split_norms.append(norm(psi_s, grid))
    split_norms = np.array(split_norms)

    euler_err = float(abs(euler_norms[-1] - 1.0))
    split_err = float(abs(split_norms[-1] - 1.0))
    # Estrictamente creciente en cada paso: lo exige ||psi_{n+1}||^2 = ||psi_n||^2 + dt^2||H psi_n||^2
    monotonic = bool(np.all(np.diff(euler_norms) > 0.0))

    t_euler = DT * np.arange(STEPS + 1)[::RECORD_EVERY]
    t_split = DT * RECORD_EVERY * np.arange(len(split_norms))
    with plt.style.context("dark_background"):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.semilogy(t_euler, np.maximum(np.abs(euler_norms[::RECORD_EVERY] - 1.0), LOG_FLOOR),
                    color="tomato", label="Euler explícito")
        ax.semilogy(t_split, np.maximum(np.abs(split_norms - 1.0), LOG_FLOOR),
                    color="cyan", label="split-step (Strang)")
        ax.set_xlabel("t")
        ax.set_ylabel(r"$|N-1|$")
        ax.set_title("Error de norma: Euler crece, split-step conserva")
        ax.legend()
        fig.savefig(OUTPUT_DIR / "02_euler_vs_splitstep.png", dpi=150)  # dpi 150: legible en pantalla
        plt.close(fig)

    print(f"euler_norm_error_final={euler_err:.6e}")
    print(f"splitstep_norm_error_final={split_err:.6e}")
    print(f"euler_norm_monotonic={monotonic}")


if __name__ == "__main__":
    main()
