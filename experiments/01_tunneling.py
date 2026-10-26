"""Experimento 01: efecto túnel de un paquete gaussiano contra una barrera rectangular.

Unidades hbar = m = 1. Un paquete (x0=-40, k0=1.3, sigma=4) con energía media E0 = k0^2/2 =
0.845 incide sobre una barrera de altura V0 = 1 y ancho 2. Se mide qué fracción de la
probabilidad queda a la derecha (T medida) y se compara con dos predicciones:
  * T(E0): coeficiente de onda plana evaluado en la energía media (ingenuo);
  * T_paquete: promedio de T(E(k)) pesado por |phi(k)|^2 (correcto por linealidad).
Las fórmulas analíticas viven SOLO en este archivo (decisión 0005). Genera
outputs/01_tunneling.png y outputs/report.md.
"""

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # backend sin ventana: solo escribe archivos
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from quantum1d import (  # noqa: E402
    Grid,
    barrier,
    energy,
    evolve,
    gaussian_packet,
    norm,
)

# --- Parámetros (physics-spec 11.1) ---
N = 4096  # nodos; con L=409.6 da dx=0.1 (20 celdas en la barrera)
L = 409.6  # largo del anillo; deja ~205 de margen a cada lado del paquete
DT = 0.005  # paso temporal; fase cinética por paso 2.47 rad <= pi (spec 6.4); O(dt^2) se verifica con CONV_DT
T_FINAL = 100.0  # tiempo final: la colisión ocurre en t~31 y ya terminó
X0 = -40.0  # posición inicial: 9.75 sigma del borde de la barrera, sin solapamiento
K0 = 1.3  # momento medio (= velocidad); E0 = 0.845 < V0
SIGMA = 4.0  # ancho; sigma_k = 0.125 hace |T_paquete - T(E0)| ~ 0.02 (distinguible)
V0 = 1.0  # altura de la barrera
A = 2.0  # ancho nominal de la barrera
BARRIER_CENTER = 0.01  # desplazado para que ningún borde (-0.99, 1.01) caiga sobre un nodo
SNAPSHOT_TIMES = (0.0, 30.0, 50.0, 100.0)  # instantes de la figura
RECORD_EVERY = 100  # pasos entre registros de norma/energía (0.5 de tiempo con DT=0.005)
EDGE_MARGIN = 20.0  # franja junto a +-L/2 donde se mide la probabilidad que "da la vuelta"
CONV_N = (8192, 0.00125)  # (n, dt) de la corrida de convergencia en dx: dx=0.05 y misma fase por paso (2.47 rad)
CONV_DT = (4096, 0.0025)  # (n, dt) de la corrida de convergencia en dt: mismo dx, dt a la mitad
# Fase cinética máxima por paso dt*kmax^2/2. La resonancia numérica entre modos de la grilla y el
# paquete aparece cuando (k^2/2 - E)*dt llega a 2*pi; pi deja un margen x2 (physics-spec 6.4, decisión 0007).
MAX_PHASE_PER_STEP = math.pi
LOG_FLOOR = 1e-17  # piso para log10 de errores nulos (bajo el redondeo ~1e-16); solo para graficar
VIEW_HALF_WIDTH = 150.0  # los paneles de |psi|^2 muestran x en [-150, 150]

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "outputs"


def transmission_plane_wave(E, V0, a):
    """Coeficiente de transmisión T(E) de una onda plana contra una barrera rectangular.

    E < V0: T = [1 + V0^2 sinh^2(kappa a) / (4 E (V0-E))]^-1, kappa = sqrt(2(V0-E)) (túnel).
    E > V0: T = [1 + V0^2 sin^2(q a) / (4 E (E-V0))]^-1, q = sqrt(2(E-V0)) (por encima).
    E = V0: límite común T = [1 + V0 a^2 / 2]^-1. Acepta escalares o arrays (vectorizada).
    """
    E = np.asarray(E, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):  # ramas descartadas por np.where
        kappa = np.sqrt(2.0 * (V0 - E))
        q = np.sqrt(2.0 * (E - V0))
        t_below = 1.0 / (1.0 + V0**2 * np.sinh(kappa * a) ** 2 / (4.0 * E * (V0 - E)))
        t_above = 1.0 / (1.0 + V0**2 * np.sin(q * a) ** 2 / (4.0 * E * (E - V0)))
    t_equal = 1.0 / (1.0 + V0 * a**2 / 2.0)
    return np.where(E < V0, t_below, np.where(E > V0, t_above, t_equal))


def transmission_packet(psi0, grid, V0, a):
    """Transmisión predicha para un paquete: sum_{k>0} |phi|^2 T(k^2/2) / sum_{k>0} |phi|^2.

    Cada componente de onda plana se transmite con su propia T (linealidad); phi = FFT(psi0).
    """
    weights = np.abs(np.fft.fft(psi0)) ** 2
    positive = grid.k > 0
    t_k = transmission_plane_wave(grid.k[positive] ** 2 / 2.0, V0, a)
    return float(np.sum(weights[positive] * t_k) / np.sum(weights[positive]))


def phase_per_step(grid, dt):
    """Fase cinética máxima acumulada por paso, dt*kmax^2/2 (rad), con kmax = pi/dx (hbar = m = 1)."""
    return float(dt * (np.pi / grid.dx) ** 2 / 2.0)


def run(n, dt, record):
    """Evoluciona el paquete hasta T_FINAL en una grilla de n nodos y paso dt.

    Si record es True, guarda cada RECORD_EVERY pasos el tiempo, la norma y la energía, y
    copias de psi en SNAPSHOT_TIMES. Devuelve (grid, V, psi0, psi_final, history).
    """
    grid = Grid(n, L)
    phase = phase_per_step(grid, dt)
    if phase > MAX_PHASE_PER_STEP:
        raise ValueError(
            f"fase por paso {phase:.4g} rad > {MAX_PHASE_PER_STEP:.4g}: riesgo de resonancia; reduzca dt o n"
        )
    V = barrier(grid, BARRIER_CENTER, A, V0)
    psi0 = gaussian_packet(grid, X0, K0, SIGMA)
    total_steps = int(round(T_FINAL / dt))
    history = {"t": [], "norm": [], "energy": [], "snapshots": {}}
    if not record:
        return grid, V, psi0, evolve(psi0, V, grid, dt, total_steps), history
    psi = psi0
    for step in range(0, total_steps + 1, RECORD_EVERY):
        t = step * dt
        history["t"].append(t)
        history["norm"].append(norm(psi, grid))
        history["energy"].append(energy(psi, V, grid))
        for ts in SNAPSHOT_TIMES:
            if abs(t - ts) < 0.5 * dt:  # el instante cae en un paso exacto de registro
                history["snapshots"][ts] = psi.copy()
        if step < total_steps:
            psi = evolve(psi, V, grid, dt, RECORD_EVERY)
    return grid, V, psi0, psi, history


def measure(psi, grid, a_eff):
    """Devuelve (T, R, P_dentro) con T=dx*sum_{x>c+a/2}|psi|^2, R análogo a la izquierda."""
    density = np.abs(psi) ** 2
    right = grid.x > BARRIER_CENTER + a_eff / 2
    left = grid.x < BARRIER_CENTER - a_eff / 2
    inside = ~(right | left)
    return (
        float(grid.dx * np.sum(density[right])),
        float(grid.dx * np.sum(density[left])),
        float(grid.dx * np.sum(density[inside])),
    )


def write_report(path, rows):
    """Escribe el reporte Markdown con tablas '| clave | valor |' (floats con repr exacto)."""
    lines = ["# Reporte del experimento 01: efecto túnel", ""]
    for title, items in rows:
        lines += [f"## {title}", "", "| clave | valor |", "|---|---|"]
        for key, value in items:
            # repr(float) es exacto (round-trip): el redondeo de impresión es 0 y AC18 mide el error real
            text = repr(float(value)) if isinstance(value, float) else str(value)
            lines.append(f"| {key} | {text} |")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _phase_curve(ax, x, density, phase, norm_phase, cmap):
    """Dibuja |psi|^2 como curva coloreada por la fase, con relleno suave bajo la curva."""
    points = np.column_stack([x, density])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    mid_phase = 0.5 * (phase[:-1] + phase[1:])
    quads = np.stack(
        [
            np.column_stack([x[:-1], np.zeros(len(x) - 1)]),
            segments[:, 0],
            segments[:, 1],
            np.column_stack([x[1:], np.zeros(len(x) - 1)]),
        ],
        axis=1,
    )
    fill = PolyCollection(quads, cmap=cmap, norm=norm_phase, alpha=0.35, linewidths=0)  # alpha: relleno suave
    fill.set_array(mid_phase)
    line = LineCollection(segments, cmap=cmap, norm=norm_phase, linewidths=1.6)
    line.set_array(mid_phase)
    ax.add_collection(fill)
    ax.add_collection(line)
    return line


def make_figure(path, grid, V, history, values):
    """Figura 3x2: 4 instantes de |psi|^2 (color = fase), errores log y barras de transmisión."""
    cmap = plt.get_cmap("hsv")  # cíclico y brillante sobre fondo oscuro (twilight es oscuro cerca de fase 0)
    norm_phase = Normalize(-np.pi, np.pi)
    with plt.style.context("dark_background"):
        fig = plt.figure(figsize=(13, 12))
        gs = GridSpec(3, 2, figure=fig, hspace=0.38, wspace=0.22)
        window = np.abs(grid.x) <= VIEW_HALF_WIDTH
        x = grid.x[window]
        line = None
        for idx, ts in enumerate(SNAPSHOT_TIMES):
            ax = fig.add_subplot(gs[idx // 2, idx % 2])
            psi = history["snapshots"][ts][window]
            density = np.abs(psi) ** 2
            line = _phase_curve(ax, x, density, np.angle(psi), norm_phase, cmap)
            ax.axvspan(
                BARRIER_CENTER - values["a_effective"] / 2,
                BARRIER_CENTER + values["a_effective"] / 2,
                color="0.6",
                alpha=0.6,
                label="barrera V(x)",
            )
            ax.set_xlim(-VIEW_HALF_WIDTH, VIEW_HALF_WIDTH)
            ax.set_ylim(0, 1.15 * float(np.max(density)))
            ax.set_title(f"t = {ts:g}")
            ax.set_xlabel("x")
            ax.set_ylabel(r"$|\psi|^2$")
            if idx == 0:
                ax.legend(loc="upper right", fontsize=8)
        fig.colorbar(line, ax=fig.axes[:4], label=r"fase de $\psi$ (rad)", ticks=[-np.pi, 0, np.pi],
                     fraction=0.025, pad=0.02).ax.set_yticklabels(["-π", "0", "π"])

        ax_err = fig.add_subplot(gs[2, 0])
        t = np.array(history["t"])
        n_err = np.maximum(np.abs(np.array(history["norm"]) - 1.0), LOG_FLOOR)
        e_err = np.maximum(
            np.abs(np.array(history["energy"]) - values["H0"]) / abs(values["H0"]), LOG_FLOOR
        )
        ax_err.semilogy(t, n_err, color="cyan", label=r"$|N-1|$")
        ax_err.semilogy(t, e_err, color="orange", label=r"$|E-E_0|/|E_0|$")
        ax_err.set_xlabel("t")
        ax_err.set_ylabel("error")
        ax_err.set_title("Conservación: norma y energía")
        ax_err.legend(fontsize=8)

        ax_bar = fig.add_subplot(gs[2, 1])
        labels = ["T medida", "T(E0)", "T paquete"]
        heights = [values["T_measured"], values["T_E0"], values["T_packet"]]
        bars = ax_bar.bar(labels, heights, color=["#4fc3f7", "#ffb74d", "#aed581"])
        for bar, h in zip(bars, heights):
            ax_bar.text(bar.get_x() + bar.get_width() / 2, h, f"{h:.4f}", ha="center", va="bottom")
        ax_bar.set_ylim(0, 1.2 * max(heights))
        ax_bar.set_ylabel("transmisión")
        ax_bar.set_title("Transmisión: medida vs predicciones")
        fig.suptitle("Efecto túnel de un paquete gaussiano (color = fase de ψ)", fontsize=14)
        fig.savefig(path, dpi=150)  # dpi 150: legible en pantalla sin pesar demasiado
        plt.close(fig)


def main():
    """Corre la simulación base y las dos de convergencia; escribe figura y reporte."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    grid, V, psi0, psi_f, hist = run(N, DT, record=True)
    a_eff = float(np.count_nonzero(V) * grid.dx)
    t_meas, r_meas, p_inside = measure(psi_f, grid, a_eff)
    density = np.abs(psi_f) ** 2
    p_edges = float(grid.dx * np.sum(density[np.abs(grid.x) > L / 2 - EDGE_MARGIN]))

    e0 = K0**2 / 2.0
    h0 = hist["energy"][0]
    t_e0 = float(transmission_plane_wave(e0, V0, a_eff))
    t_packet = transmission_packet(psi0, grid, V0, a_eff)

    grid_n, V_n, _, psi_n, _ = run(*CONV_N, record=False)
    t_n = measure(psi_n, grid_n, float(np.count_nonzero(V_n) * grid_n.dx))[0]
    grid_d, V_d, _, psi_d, _ = run(*CONV_DT, record=False)
    t_d = measure(psi_d, grid_d, float(np.count_nonzero(V_d) * grid_d.dx))[0]

    norm_err = np.abs(np.array(hist["norm"]) - 1.0)
    e_dev = np.abs(np.array(hist["energy"]) - h0) / abs(h0)

    values = {
        "a_effective": a_eff,
        "H0": h0,
        "T_measured": t_meas,
        "T_E0": t_e0,
        "T_packet": t_packet,
    }
    rows = [
        ("Parámetros", [
            ("N", N), ("L", L), ("dx", grid.dx), ("dt", DT), ("steps", int(round(T_FINAL / DT))),
            ("t_final", T_FINAL), ("x0", X0), ("k0", K0), ("sigma", SIGMA), ("V0", V0), ("a", A),
            ("barrier_center", BARRIER_CENTER), ("n_barrier_cells", int(np.count_nonzero(V))),
            ("a_effective", a_eff), ("kmax", float(np.pi / grid.dx)),
            ("guard_rhs", 4.0 * (abs(K0) + 3.0 / SIGMA)),
        ]),
        ("Energías", [("E0", e0), ("H0", h0)]),
        ("Probabilidades al tiempo final", [
            ("T_measured", t_meas), ("R_measured", r_meas),
            ("P_inside_barrier", p_inside), ("P_near_edges", p_edges),
        ]),
        ("Predicciones", [
            ("T_E0", t_e0), ("T_packet", t_packet),
            ("abs_diff_measured_packet", abs(t_meas - t_packet)),
            ("abs_diff_measured_E0", abs(t_meas - t_e0)),
        ]),
        ("Convergencia", [
            ("T_measured_convN", t_n), ("T_measured_convdt", t_d),
            ("convN_N", CONV_N[0]), ("convN_dt", CONV_N[1]),
            ("convdt_N", CONV_DT[0]), ("convdt_dt", CONV_DT[1]),
            ("conv_diff_N", abs(t_n - t_meas)), ("conv_diff_dt", abs(t_d - t_meas)),
        ]),
        ("Conservación", [
            ("norm_error_max", float(np.max(norm_err))),
            ("energy_rel_dev_max", float(np.max(e_dev))),
            ("energy_rel_dev_final", float(e_dev[-1])),
            ("phase_per_step_base", phase_per_step(grid, DT)),
            ("phase_per_step_convN", phase_per_step(grid_n, CONV_N[1])),
            ("phase_per_step_convdt", phase_per_step(grid_d, CONV_DT[1])),
        ]),
    ]
    make_figure(OUTPUT_DIR / "01_tunneling.png", grid, V, hist, values)
    write_report(OUTPUT_DIR / "report.md", rows)
    print(f"T_measured={t_meas:.6f} T_E0={t_e0:.6f} T_packet={t_packet:.6f}")
    print(f"conv_diff_N={abs(t_n - t_meas):.3e} conv_diff_dt={abs(t_d - t_meas):.3e}")
    print(f"norm_error_max={np.max(norm_err):.3e} energy_rel_dev_max={np.max(e_dev):.3e}")


if __name__ == "__main__":
    main()
