"""Referencias analíticas independientes en unidades naturales, ℏ = m = 1.

No importa quantum1d. σ siempre es la desviación de la densidad |ψ|²,
no la anchura del exponente de la amplitud. La energía es E = k²/2.
"""

import numpy as np


def free_center(x0, k0, time):
    """Centro libre: la velocidad de grupo es p/m = k₀ en estas unidades."""
    return x0 + k0 * np.asarray(time)


def free_width(sigma0, time):
    """Dispersión exacta de una gaussiana libre; σ(0) = σ₀."""
    return sigma0 * np.sqrt(1.0 + (np.asarray(time) / (2 * sigma0**2)) ** 2)


def gaussian_density(x, center, sigma):
    """Densidad normalizada sobre la recta; sus unidades son longitud⁻¹."""
    return np.exp(-0.5 * ((np.asarray(x) - center) / sigma) ** 2) / (
        np.sqrt(2 * np.pi) * sigma
    )


def harmonic_center(x0, omega, time):
    """Centro de un estado coherente inicialmente en reposo: x₀ cos(ωt)."""
    return x0 * np.cos(omega * np.asarray(time))


def barrier_transmission(energy, height, width):
    """Transmisión exacta de ondas planas por una barrera rectangular.

    Se toma incidencia desde la izquierda, E >= 0, V₀ >= 0 y ancho >= 0.
    En E = V₀ se usa el límite continuo 1/(1 + V₀ a²/2).
    Bajo barrera se evalúa log(sinh) para evitar overflow: una barrera
    muy opaca debe dar transmisión cero a precisión máquina, nunca NaN.
    En la barrera nula se adopta T = 1 incluso en el límite E = 0.
    """
    values = np.asarray(energy, dtype=float)
    if (
        not np.isfinite(height)
        or not np.isfinite(width)
        or height < 0
        or width < 0
        or np.any(~np.isfinite(values))
        or np.any(values < 0)
    ):
        raise ValueError("E, altura y ancho deben ser finitos y no negativos")
    flat = values.reshape(-1)
    result = np.zeros_like(flat)
    if height == 0 or width == 0:
        result.fill(1.0)
    else:
        threshold = flat == height
        result[threshold] = 1.0 / (1.0 + height * width**2 / 2.0)
        below = (flat > 0) & (flat < height)
        above = flat > height
        with np.errstate(divide="ignore", under="ignore"):
            delta = height - flat[below]
            z = np.sqrt(2.0) * np.sqrt(delta) * width
            log_sinh = z + np.log(-np.expm1(-2 * z)) - np.log(2.0)
            log_ratio = (
                2 * np.log(height) + 2 * log_sinh - np.log(4.0)
                - np.log(flat[below]) - np.log(delta)
            )
            result[below] = np.exp(-np.logaddexp(0.0, log_ratio))
            delta = flat[above] - height
            phase = np.sqrt(2.0) * np.sqrt(delta) * width
            log_ratio = (
                2 * np.log(height) + 2 * np.log(np.abs(np.sin(phase)))
                - np.log(4.0) - np.log(flat[above]) - np.log(delta)
            )
            result[above] = np.exp(-np.logaddexp(0.0, log_ratio))
    result = result.reshape(values.shape)
    return result.item() if result.ndim == 0 else result


def packet_transmission(psi, dx, height, width):
    """Promedio de T(E) con el espectro FFT inicial, restringido a k > 0.

    Para una grilla uniforme, Δk y los factores de normalización de la
    transformada se cancelan entre numerador y denominador. No se usa
    ningún diagnóstico ni número de onda provisto por el motor auditado.
    """
    psi = np.asarray(psi)
    k = 2 * np.pi * np.fft.fftfreq(psi.size, d=dx)
    positive = k > 0
    weights = np.abs(np.fft.fft(psi)[positive]) ** 2
    total = np.sum(weights)
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Se necesita peso espectral positivo finito")
    transmission = barrier_transmission(k[positive] ** 2 / 2, height, width)
    return float(np.sum(weights * transmission) / total)
