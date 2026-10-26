<!-- dev-team plan -->
# Plan Fase 0, ROL A: motor de referencia `quantum1d`

## Handoff

### Goal
Un paquete Python `quantum1d` (src layout) que resuelve Schrödinger 1D dependiente del tiempo con split-step de Strang según la API del contrato, más dos experimentos que generan figuras y `outputs/report.md`, de modo que el ROL B pueda auditarlo a ciegas y el usuario entienda la física leyendo docstrings en español.

### Context
- `C:/dev/quantum-A/docs/physics-spec.md`: física, convenciones exactas de la API (sección 12), parámetros y cuentas de los experimentos (sección 11). Es la fuente de verdad para el coder; solo lectura.
- `C:/dev/quantum-A/docs/decisions/0001`...`0006`: decisiones aceptadas (checks con PYTHONPATH, guarda, diagnostics, barrera, ubicación de fórmulas, semántica de evolve).
- `C:/dev/quantum-A/README.md`: único archivo previo; no se toca.
- Contrato del usuario (mensaje del orquestador): API pública, física a-g, reglas permanentes 1-6, reparto A/B.

### Constraints
- La API del contrato es fija: nombres, firmas y semántica de la sección 12 de physics-spec.md. Nada de parámetros ni funciones públicas extra.
- Unidades ħ = m = 1; `float64`/`complex128`; FFT con `numpy.fft`; código compatible con Python 3.12 (nada exclusivo de 3.13+), aunque corre en el 3.14 local.
- No instalar nada. Todo se ejecuta desde `C:/dev/quantum-A` como `PYTHONPATH=src python ...` (en PowerShell: `$env:PYTHONPATH='src'; python ...`).
- Identificadores en inglés; docstrings y comentarios en español que expliquen qué representa, unidades e invariante.
- Cada tolerancia numérica lleva un comentario que la justifica. Si un check falla, no se relaja la tolerancia: se reporta la causa.
- src/ sin fórmulas analíticas de verificación ni Euler; las fórmulas de T solo en `experiments/01_tunneling.py`; Euler solo en `experiments/02_euler_vs_splitstep.py`.
- No tocar `tests/` ni `docs/verification.md` (ROL B). No hacer merge, rebase ni push. No commitear `__pycache__`.
- Los scripts imprimen solo ASCII en stdout; los archivos se escriben con `encoding="utf-8"`. Backend de matplotlib `Agg`. Las rutas de salida se resuelven con `Path(__file__).resolve().parent.parent / "outputs"`.
- Sanity checks del tester: scripts en su scratchpad, fuera del repo, sin commitear (decisión 0005).

### Files to touch
- `pyproject.toml` (create): metadatos y stack para uv (paso 1).
- `AGENTS.md` (create): las 6 reglas permanentes textuales, reparto A/B y cómo correr (paso 2).
- `CLAUDE.md` (create): solo remite a AGENTS.md (paso 2).
- `src/quantum1d/__init__.py` (create): reexporta la API pública (paso 8).
- `src/quantum1d/grid.py` (create): `Grid` y `check_resolution` (paso 3).
- `src/quantum1d/potentials.py` (create): `free`, `barrier`, `harmonic` (paso 4).
- `src/quantum1d/wavepacket.py` (create): `gaussian_packet` (paso 5).
- `src/quantum1d/evolve.py` (create): `evolve` (paso 6).
- `src/quantum1d/diagnostics.py` (create): `norm`, `mean_x`, `mean_p`, `std_x`, `energy` (paso 7).
- `experiments/01_tunneling.py` (create): experimento de túnel, figura y reporte (paso 9).
- `experiments/02_euler_vs_splitstep.py` (create): Euler vs split-step (paso 10).
- `outputs/01_tunneling.png`, `outputs/02_euler_vs_splitstep.png`, `outputs/report.md` (create, generados): resultado de correr los experimentos (paso 11).
- `STATUS.md` (create): qué está hecho, qué falla, "Lo que NO verifiqué", mejoras no implementadas (paso 12).

### Out of scope
- Web, UI, 2D, spin, qubits, GPU, port a TypeScript.
- `tests/`, `tests/analytic.py`, `docs/verification.md` y cualquier test pytest (ROL B).
- Editar `docs/physics-spec.md`, `docs/plan.md`, `docs/decisions/` (arquitecto). Si el coder encuentra un error en la spec, lo reporta.
- Guardas de resolución dentro de `evolve`, renormalización, métodos de orden superior, potenciales extra, `.gitignore`, CLI, logging.

### Acceptance criteria
Comandos en bash desde `C:/dev/quantum-A`; `P=` abrevia `PYTHONPATH=src python`. Los checks marcados "script" los escribe el tester en su scratchpad con exactamente el setup indicado.
- AC1: La API es importable desde `quantum1d` y `quantum1d.evolve` es la función. Check: `P -c "import quantum1d as q; n=['Grid','check_resolution','gaussian_packet','evolve','potentials','diagnostics','free','barrier','harmonic','norm','mean_x','mean_p','std_x','energy']; print([x for x in n if not hasattr(q,x)], callable(q.evolve)); from quantum1d.potentials import free,barrier,harmonic; from quantum1d.diagnostics import norm,mean_x,mean_p,std_x,energy; from quantum1d.grid import Grid,check_resolution"`. Pass: imprime `[] True` y sale con código 0.
- AC2: `Grid(8, 4.0)` cumple `array_equal(x, -2.0+0.5*arange(8))`, `dx == 0.5`, `array_equal(k, 2*pi*fftfreq(8, d=0.5))`, `n == 8`, `length == 4.0`. Check: `P -c` con esas 5 comparaciones impresas. Pass: 5 veces `True`.
- AC3: `Grid(1, 1.0)`, `Grid(16, 0.0)`, `Grid(16, -1.0)` y `Grid(16, float('inf'))` lanzan `ValueError`. Check: `P -c` con un `try/except ValueError` por caso. Pass: los 4 lanzan `ValueError`.
- AC4: Con `g = Grid(256, 20.0)` (π/dx = 40.21): `check_resolution(g, 5.0, 1.0)` devuelve `None` (32 ≤ 40.21); `check_resolution(g, 8.0, 1.0)`, `check_resolution(g, -8.0, 1.0)` (44 > 40.21) y `gaussian_packet(g, 0.0, 8.0, 1.0)` lanzan `ValueError`. Check: `P -c` por caso. Pass: lo indicado en cada caso.
- AC5: Con `g = Grid(1024, 100.0)` y `psi = gaussian_packet(g, -5.0, 1.3, 2.0)`: dtype `complex128`; `max|psi - fórmula| <= 1e-14`; `|norm-1|`, `|mean_x+5|`, `|std_x-2|`, `|mean_p-1.3|` y `|energy(psi, free(g), g) - 0.87625|` (= 1.3²/2 + 1/32) todos `<= 1e-10`. Check: script. Pass: todas las condiciones `True` (la suma discreta de una gaussiana bien muestreada es exacta salvo redondeo, ~1e-14).
- AC6: Con `g = Grid(8192, 409.6)`: `barrier(g, 0.01, 2.0, 1.0)` es `float64`, valores en {0.0, 1.0}, `count_nonzero == 40`, todo nodo con V=1 cumple `|x-0.01| <= 1.0`; `harmonic(g, 2.0, center=1.0)` coincide con `2.0*(x-1)**2` (rtol 1e-14); `free(g)` es `zeros(8192)` `float64`. Check: script. Pass: todo `True`.
- AC7: `evolve` no muta y devuelve un array nuevo. Con `g = Grid(256, 20.0)`, `psi = gaussian_packet(g, 1.0, 0.5, 0.7071)`, `V = harmonic(g, 1.0)`: después de `out = evolve(psi, V, g, 0.01, 10)`, `psi` y `V` son iguales a sus copias previas, `out is not psi`, `not shares_memory(out, psi)`, `out.dtype == complex128`; `evolve(psi, V, g, 0.01, 0)` es igual a `psi` sin compartir memoria; `max|evolve(evolve(psi,V,g,0.01,3),V,g,0.01,4) - evolve(psi,V,g,0.01,7)| <= 1e-14`. Check: script. Pass: todo `True`.
- AC8: `evolve` lanza `ValueError` con: `steps=-1`, `steps=1.5`, `psi[:-1]`, `V[:-1]`, `V+0j`, `dt=nan`, `dt=inf` (mismo setup que AC7). Check: script con `try/except ValueError` por caso. Pass: los 7 lanzan `ValueError`.
- AC9: Con `z = zeros(256, complex)` en `Grid(256, 20.0)`: `norm(z, g) == 0.0`; `mean_x`, `std_x`, `mean_p` y `energy(z, free(g), g)` lanzan `ValueError`; `barrier(g, 0, 0.0, 1)`, `barrier(g, 0, -1.0, 1)`, `harmonic(g, 0.0)`, `harmonic(g, -1.0)`, `gaussian_packet(g, 0, 0, 0.0)` y `check_resolution(g, 0, -1.0)` lanzan `ValueError`. Check: script. Pass: todo según lo indicado.
- AC10: Reversibilidad. `g = Grid(2048, 200.0)`, `V = barrier(g, 0.01, 2.0, 1.0)`, `psi = gaussian_packet(g, -10.0, 1.3, 4.0)`, `back = evolve(evolve(psi,V,g,0.02,500),V,g,-0.02,500)`. Check: script. Pass: `max|back-psi| < 1e-10` (redondeo esperado ~1e-13).
- AC11 (a): Norma tras 2000 pasos, `|norm-1| < 1e-10` en 3 casos: (i) `Grid(2048, 200.0)`, `barrier(g,0.01,2.0,1.0)`, `gaussian_packet(g,-30,1.3,4)`, dt=0.02; (ii) `Grid(256, 20.0)`, `harmonic(g,1.0)`, `gaussian_packet(g,2,0,sqrt(0.5))`, dt=0.01; (iii) `Grid(512, 50.0)`, `V = default_rng(12345).uniform(0, 2, 512)`, `gaussian_packet(g,0,1,2)`, dt=0.01. Check: script. Pass: los 3 casos `< 1e-10`.
- AC12 (b): Libre: `Grid(1024, 200.0)`, `gaussian_packet(g,-20,2,2)`, dt=0.01, diagnósticos a t = 5, 10, 15, 20 (500 pasos por tramo). Check: script. Pass: en los 4 tiempos `|mean_x-(x0+k0 t)| < 1e-8`, `|std_x/σ(t)-1| < 1e-8` con σ(t) = 2·sqrt(1+(t/8)²), `|mean_p-2| < 1e-10` (con V=0 el método es exacto: solo redondeo ~1e-12).
- AC13 (c): Estado base: `Grid(256, 20.0)`, `harmonic(g,1.0)`, `gaussian_packet(g,0,0,sqrt(0.5))`, dt=0.01, 2000 pasos, registro cada 10. Check: script. Pass: `max|std_x/σ0-1| < 1e-4` (predicción ≤ dt²/8 = 1.25e-5) y `max|mean_x| < 1e-10` (simetría).
- AC14 (c): Desplazado: igual que AC13 pero `gaussian_packet(g,2,0,sqrt(0.5))`. Check: script. Pass: `max_t |mean_x - 2cos t| < 5e-4` (predicción x0·t·dt²/24 = 1.7e-4 en t=20).
- AC15 (d): Convergencia: `Grid(256, 20.0)`, `harmonic(g,1.0)`, `psi0 = gaussian_packet(g,1.0,0.5,sqrt(0.5))`, t=2, dt ∈ {0.04, 0.02, 0.01, 0.005}, referencia dt=0.000625; error `sqrt(dx·Σ|psi-psi_ref|²)`. Check: script con `polyfit(log dt, log err, 1)`. Pass: pendiente en [1.9, 2.1] y cada cociente err(dt)/err(dt/2) en [3.6, 4.4].
- AC16 (e): Energía: setup de AC14 con dt=0.01 (2000 pasos) y dt=0.005 (4000 pasos), `energy` cada 0.1 de tiempo. Check: script. Pass: pico a pico con dt=0.01 en [4.5e-5, 5.5e-5] (predicción dt²x0²/8 = 5.0e-5, spec 6.3); cociente pico-a-pico(0.01)/pico-a-pico(0.005) en [3.6, 4.4]; `max|E-E0|` en t ∈ [10,20] ≤ 1.1 × el de t ∈ [0,10] (sin deriva).
- AC17: `P experiments/01_tunneling.py` sale con código 0 y deja `outputs/01_tunneling.png` (> 0 bytes) y `outputs/report.md` con todas las claves del paso 9.6. Check: correr el script; `P -c` que parsea el reporte con la regex `^\| (\w+) \| ([^|]+) \|$`. Pass: código 0, archivos presentes, ninguna clave faltante.
- AC18 (f): Valores de `outputs/report.md`: `abs_diff_measured_packet < 0.01`; `conv_diff_N < 1e-3`; `conv_diff_dt < 1e-3`; `P_near_edges < 1e-10`; `P_inside_barrier < 1e-6`; `|T_measured+R_measured+P_inside_barrier-1| < 1e-10`; `|T_E0-0.221138| < 1e-5`; `0.236 < T_packet < 0.246`; `abs_diff_measured_E0 > 0.01`; `norm_error_max < 1e-10`; `energy_rel_dev_max < 1e-3`; `n_barrier_cells == 40`; `a_effective == 2.0`; `kmax >= guard_rhs`. Check: `P -c` sobre el reporte parseado. Pass: las 14 condiciones `True`.
- AC19: Fórmula de T del experimento, cargada con `importlib.util.spec_from_file_location("tun", "experiments/01_tunneling.py")` (sin ejecutar `main`): `transmission_plane_wave(0.845,1,2)` = 0.221138 ± 1e-5; `(1.0,1,2)` = 1/3 ± 1e-12; `(1.125,1,2)` = 0.442713 ± 1e-5; `(0.5,1,2)` = 0.070651 ± 1e-5; `|T(1±1e-9)-1/3| < 1e-6`; con `E = array([0.5,1.0,1.125])` devuelve forma `(3,)`. Check: script. Pass: todo `True` y la importación no genera archivos.
- AC20: La figura 01 muestra fondo oscuro, 4 paneles de |psi|² en t = 0, 30, 50, 100 con la curva y un relleno coloreados por la fase (colormap cíclico y barra de color de -π a π), la barrera sombreada, un panel log de error de norma y de energía vs t, y un panel de 3 barras rotuladas (T medida, T(E0), T paquete) con su valor. Check: abrir `outputs/01_tunneling.png` y revisarla. Pass: los 6 elementos presentes y legibles.
- AC21 (g): `P experiments/02_euler_vs_splitstep.py` sale con código 0, crea `outputs/02_euler_vs_splitstep.png` e imprime `euler_norm_error_final=<v>` con v en [3.75e-4, 1.5e-3] (predicción 7.5e-4, spec 11.2), `splitstep_norm_error_final=<v>` con v < 1e-10, y `euler_norm_monotonic=True`. Check: correr y leer stdout. Pass: las 3 líneas cumplen y el PNG existe.
- AC22: src/ no contiene Euler ni fórmulas de T. Check: `grep -rniE "def [a-z_]*euler|sinh|transmission" src/`. Pass: sin salida.
- AC23: El ROL A no tocó lo del ROL B. Check: `git diff --name-only main -- tests docs/verification.md` y `git status --porcelain -- tests docs/verification.md`. Pass: ambas salidas vacías.
- AC24: `CLAUDE.md` tiene ≤ 5 líneas no vacías, menciona `AGENTS.md` y no repite reglas; `AGENTS.md` contiene, textuales, `Unidades naturales: hbar = m = 1`, `Nunca declares algo "correcto" sin un test que pase`, `Cada tolerancia numérica lleva un comentario`, `Sin features extra`, `Lo que NO verifiqué`, `No gastes más del 10% del esfuerzo en documentación`. Check: `P -c` leyendo ambos en UTF-8. Pass: todo `True`.
- AC25: `pyproject.toml` (leído con `tomllib`) tiene `requires-python == ">=3.12"`, dependencias que empiezan con numpy, scipy y matplotlib, `dependency-groups.dev` con pytest, `build-backend == "hatchling.build"` y `tool.pytest.ini_options.pythonpath == ["src"]`. Check: `P -c`. Pass: todo `True`.
- AC26: `STATUS.md` tiene los encabezados `## Hecho`, `## Falla`, `## Lo que NO verifiqué`, `## Mejoras no implementadas` y cita T_measured y T_packet del reporte. Check: `grep -c` por encabezado. Pass: cada uno aparece 1 vez.
- AC27: Toda clase y función pública de `src/quantum1d/*.py` y toda función de `experiments/*.py` tiene docstring (no vacío, en español), y cada literal de tolerancia (`1e-`, `rtol`, `atol`, umbrales) en src/ y experiments/ tiene un comentario de justificación en su línea o la anterior. Check: script con `ast` para docstrings; `grep -rnE "1e-|tol" src experiments` revisado a mano. Pass: 0 sin docstring, 0 tolerancias sin justificar.

## Diseño y trade-offs
- Módulos uno a uno con el contrato; `__init__` reexporta funciones y módulos (`from .evolve import evolve` deja `quantum1d.evolve` como función, documentado en spec 12). Alternativa rechazada: exponer solo módulos; B podría importar `from quantum1d import norm`.
- Guarda en `check_resolution` llamada por `gaussian_packet` (decisión 0002). Rechazado: estimar el espectro dentro de `evolve` (umbral arbitrario, falsos positivos, API implícita).
- `evolve` sin fusionar medios pasos: un poco más lento, pero cada llamada es exactamente N pasos de Strang y dividir en tramos es idéntico a una sola llamada; además es más legible para el usuario. Rechazado: fusionar `e^{-iVdt}` entre pasos.
- `dt` negativo permitido: Strang simétrico da la inversa exacta (AC10). Validaciones mínimas (formas, finitud, enteros) para que nada falle en silencio; no cambian firmas.
- Barrera con intervalo cerrado y ancho efectivo (decisión 0004); el experimento desplaza el centro 0.01.
- Fórmulas de T dentro del experimento 01 como funciones importables; sanity checks del tester fuera del repo (decisión 0005). Rechazado: un `experiments/theory.py` (archivo fuera del layout pedido).
- Parámetros del experimento 01 (spec 11.1): σ = 4 da |T_paquete − T(E0)| ≈ 0.02, el doble de la tolerancia, y todavía permite un dominio de 409.6. σ = 10 da una diferencia de 0.003 (no distinguible); σ = 3 da más contraste, pero dispersa más y deja ~25 % del paquete por encima de V0.

## Riesgos y casos borde
- Bordes de barrera sobre nodos: mitigado con center=0.01 (válido para dx = 0.1, 0.05, 0.025, 0.0125).
- Energía con barrera: los conmutadores de 6.3 son grandes en el salto de V; se estima ~1e-4 relativo, con un umbral de 1e-3. Si falla, se reporta (no se relaja).
- Python 3.14 local vs 3.12 declarado: no usar sintaxis nueva (PEP 695, `type X = ...`).
- Consola de Windows: imprimir solo ASCII.
- `log(0)` en los gráficos de error: piso de 1e-17 con comentario.
- El ROL B podría esperar un `ValueError` de `evolve` por resolución: queda documentado en spec 8 y 12 (pregunta abierta 1).
- Tiempo de cómputo del exp01 (3 corridas, ~40k pasos con FFT de 8k-16k): estimado < 1 min.

## Pasos ordenados
1. `pyproject.toml`: `[build-system] requires=["hatchling"]`, `build-backend="hatchling.build"`; `[project] name="quantum1d"`, `version="0.1.0"`, `readme="README.md"`, `requires-python=">=3.12"`, `dependencies=["numpy>=1.26","scipy>=1.11","matplotlib>=3.8"]`; `[dependency-groups] dev=["pytest>=8"]`; `[tool.hatch.build.targets.wheel] packages=["src/quantum1d"]`; `[tool.pytest.ini_options] testpaths=["tests"]`, `pythonpath=["src"]`.
2. `AGENTS.md`: las 6 reglas copiadas textuales del contrato, el reparto de dueños A/B, el comando `PYTHONPATH=src python ...` y punteros a `docs/physics-spec.md` y `docs/decisions/`. `CLAUDE.md`: una línea `@AGENTS.md` y una frase "Las reglas del proyecto están en AGENTS.md".
3. `grid.py`: `class Grid` con `__init__(self, n, length)`: `n = operator.index(n)` (un `TypeError` se convierte en `ValueError`), exigir `n >= 2`; `length = float(length)` finito y `> 0`; `dx = length/n`; `x = -length/2 + dx*np.arange(n)`; `k = 2*np.pi*np.fft.fftfreq(n, d=dx)`. `check_resolution(grid, k0, sigma) -> None` según spec 8; el mensaje de error incluye π/dx y el valor requerido.
4. `potentials.py`: `free(grid)`, `barrier(grid, center, width, height)` con `np.where(np.abs(grid.x-center) <= width/2, float(height), 0.0)`, `harmonic(grid, omega, center=0.0) = 0.5*omega**2*(grid.x-center)**2`; validaciones de spec 12. Docstrings: efecto de discretización (spec 5).
5. `wavepacket.py`: `gaussian_packet(grid, x0, k0, sigma)`: llama `check_resolution`, evalúa la fórmula de spec 4 en `grid.x` y devuelve `complex128` sin renormalizar.
6. `evolve.py`: `evolve(psi, V, grid, dt, steps)`: validar (spec 12; `steps` con `operator.index` sin aceptar `bool`; `dt` con `float()` y `math.isfinite`; `np.isrealobj(V)`); `psi = np.array(psi, dtype=np.complex128, copy=True)`; precalcular `half_v = np.exp(-0.5j*V*dt)` y `kinetic = np.exp(-0.5j*grid.k**2*dt)`; por paso: `psi = half_v*psi`; `psi = np.fft.ifft(kinetic*np.fft.fft(psi))`; `psi = half_v*psi`. Docstring: unitariedad, orden 2, H modificado (spec 6).
7. `diagnostics.py`: las 5 funciones con las definiciones exactas de spec 7; helper privado `_checked_norm` que lanza `ValueError` si la norma es 0.
8. `__init__.py`: reexporta todo lo de spec 12 y define `__all__`; docstring del paquete con unidades.
9. `experiments/01_tunneling.py`:
   1. Constantes con comentario: `N=8192`, `L=409.6`, `DT=0.01`, `T_FINAL=100.0`, `X0=-40.0`, `K0=1.3`, `SIGMA=4.0`, `V0=1.0`, `A=2.0`, `BARRIER_CENTER=0.01`, `SNAPSHOT_TIMES=(0,30,50,100)`, `RECORD_EVERY=50`, `EDGE_MARGIN=20.0`.
   2. `transmission_plane_wave(E, V0, a)`: vectorizada con las 3 ramas de spec 10 (`E == V0` exacto usa el límite) y `np.errstate` para las ramas descartadas.
   3. `transmission_packet(psi0, grid, V0, a) -> float`: `phi = fft(psi0)`, suma sobre `grid.k > 0` según spec 10.
   4. `run(n, dt, record)`: arma la grilla, la barrera y el paquete; evoluciona por tramos de `RECORD_EVERY` pasos; si `record`, guarda `t`, `norm`, `energy` y snapshots; devuelve el psi final y el historial.
   5. Regiones: `right = x > c+a/2`, `left = x < c-a/2`, `inside` el resto; `near_edges = |x| > L/2 - EDGE_MARGIN`.
   6. `write_report(path, rows)`: Markdown con tablas `| clave | valor |` (floats con `f"{v:.10g}"`), con exactamente estas claves: N, L, dx, dt, steps, t_final, x0, k0, sigma, V0, a, barrier_center, n_barrier_cells, a_effective, kmax, guard_rhs, E0, H0, T_measured, R_measured, P_inside_barrier, P_near_edges, T_E0, T_packet, abs_diff_measured_packet, abs_diff_measured_E0, T_measured_N16384, T_measured_dt0005, conv_diff_N, conv_diff_dt, norm_error_max, energy_rel_dev_max, energy_rel_dev_final. `a` en las fórmulas es `a_effective`.
   7. `make_figure(...)`: `dark_background`, GridSpec 3×2. Cuatro snapshots en x ∈ [-150, 150] con `LineCollection` y `PolyCollection` (quads bajo la curva, alpha ≈ 0.35) coloreados por `np.angle(psi)` con `twilight` y `Normalize(-π, π)`; `axvspan` gris para la barrera; colorbar de fase. Panel log con `|N-1|` y `|E-E0|/|E0|` y piso de 1e-17. Panel de 3 barras con su valor anotado. `dpi=150`.
   8. `main()`: corrida base con registro; corridas de convergencia (`n=16384`; `dt=0.005`) sin registro; escribe la figura y el reporte; imprime un resumen ASCII. Guarda `if __name__ == "__main__": main()`.
10. `experiments/02_euler_vs_splitstep.py`: parámetros de spec 11.2 (`Grid(128, 18.0)`, `harmonic(g,1.0)`, `gaussian_packet(g,1.0,0.0,sqrt(0.5))`, `dt=2.5e-4`, 8000 pasos, registro cada 20); `euler_step(psi, V, grid, dt)` local (Hψ con cinética vía FFT más Vψ); split-step con `evolve` en tramos de 20; figura semilog de `|N-1|` para ambos métodos (piso de 1e-17) en `outputs/02_euler_vs_splitstep.png`; imprime las 3 líneas de AC21.
11. Correr ambos experimentos y commitear `outputs/`. Comprobar a mano AC17, AC18 y AC21.
12. `STATUS.md`: Hecho / Falla / Lo que NO verifiqué / Mejoras no implementadas, con los números del reporte. Lo actualiza el coder al final del ciclo con el resultado del tester.
13. Tester: scripts de sanity en su scratchpad para AC2-AC16, AC19 y AC27; corre AC1, AC17-AC18 y AC20-AC26. No escribe en `tests/`.

## Preguntas abiertas y supuestos
1. (No bloqueante) Supuesto: el ROL B lee `docs/physics-spec.md`. Si no lo lee, conviene pasarle la sección 12, sobre todo que `evolve` no lanza `ValueError` por resolución (sí lo hacen `gaussian_packet` y `check_resolution`) y que `norm` es la integral sin raíz.
2. (No bloqueante) Las validaciones extra (`Grid` con n<2, `width<=0`, `omega<=0`, `dt` no finito, `V` complejo) no cambian firmas, pero no están en el contrato textual.
3. Supuesto: `T_packet ≈ 0.241` sale de una cuadratura a mano con paso 0.05 en k; el valor que manda es el del script (AC18 acepta [0.236, 0.246]).
4. Supuesto: `main` existe como rama local para AC23.

## Memory
- Decisiones escritas: `docs/decisions/0001-stack-uv-y-checks-locales.md`, `0002-guarda-de-resolucion.md`, `0003-convenciones-de-diagnostics.md`, `0004-barrera-discretizada.md`, `0005-ubicacion-de-formulas-y-sanity-checks.md`, `0006-semantica-de-evolve.md`.
- Lecciones: 3 agregadas a `docs/lessons-learned.md` (convergencia con V=0, bordes de barrera sobre nodos, sin uv en el entorno).
