# Estado de la Fase 0 (ROL A)

## Hecho

- Paquete `quantum1d` (`Grid`, `check_resolution`, `gaussian_packet`, `free`, `barrier`, `harmonic`, `evolve`, diagnósticos) según `docs/physics-spec.md` sección 12.
- `experiments/01_tunneling.py` y `experiments/02_euler_vs_splitstep.py`, con `outputs/01_tunneling.png`, `outputs/02_euler_vs_splitstep.png` y `outputs/report.md`.
- Verificación rápida del coder con scripts temporales fuera del repo (no son tests): AC1-AC17, AC19, AC21, AC22, AC24, AC25, AC27 pasan. Números: norma `|N-1|` < 5e-13 en AC11; reversibilidad 3.3e-14; convergencia en dt con pendiente 2.007; pico a pico de energía del oscilador 4.9999e-5.
- Experimento 01 (`outputs/report.md`): T_measured = 0.2412214, T_packet = 0.2412562 (diferencia 3.5e-5), T_E0 = 0.2211382 (diferencia con la medida 0.0201); conv_diff_N = 7.0e-5, conv_diff_dt = 6.0e-5; norm_error_max = 3.3e-12.
- Experimento 02: `euler_norm_error_final` = 7.509e-4, `splitstep_norm_error_final` = 2.46e-12, `euler_norm_monotonic=True`.

## Falla

AC18 falla en 3 de 14 condiciones (no se tocaron parámetros ni tolerancias):

| Clave | Medido | Umbral |
|---|---|---|
| energy_rel_dev_max | 0.1505 | < 1e-3 |
| P_near_edges | 5.43e-6 | < 1e-10 |
| P_inside_barrier | 1.136e-6 | < 1e-6 |

Hipótesis: con dx = 0.05 y dt = 0.01, el número de onda de Nyquist k_max = 62.8 da k_max^2 dt / 2 = 19.7 rad por paso, muy por encima de pi. La discontinuidad de V excita componentes de alto k (P(|k|>5) ~ 1.4e-4 a t >= 50) que el split-step propaga con fase mal resuelta (resonancia de Strang). Esa radiación de alto k pesa mucho en la energía cinética (k^2/2 ~ 2000), recorre el anillo (v = k) y llega a los bordes. La transmisión no se ve afectada (T_measured coincide con T_packet a 3.5e-5). Evidencia: con dt = 0.0025 (mismo n) la desviación de energía baja a ~1e-7 relativo y P(|k|>5) a 8e-17 en t = 100; con n = 4096 (dx = 0.1) y dt = 0.01 la energía se conserva a ~5e-7. Las condiciones de AC18 sobre energía y bordes probablemente suponen que ese efecto no ocurre (physics-spec 6.3 advierte que la serie de Hamiltoniano modificado es formal si dt k_max^2/2 >> 1 y lo estimaba en ~1e-4). Decisión pendiente del arquitecto: cambiar parámetros (dt menor o dx mayor) o criterios.

## Lo que NO verifiqué

- AC20 (figura 01) y AC23 sólo por inspección visual / `git`; no hay test automático.
- No corrí pytest (no hay tests en esta rama; son del ROL B).
- Compatibilidad con Python 3.12: sólo corrí en 3.14.
- El efecto de la resonancia de Strang sólo lo probé con 3 corridas de diagnóstico (n y dt variados), no con un barrido sistemático.
- No instalé el paquete con hatchling/uv; sólo se usó `PYTHONPATH=src`.

## Mejoras no implementadas

- Guarda de estabilidad espectral en `evolve` (aviso si dt * k_max^2 / 2 >> pi).
- Filtro de alto k o paso adaptativo para barreras de borde abrupto.
- Parámetros de CLI para los experimentos.
