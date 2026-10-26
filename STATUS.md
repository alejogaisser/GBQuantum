# Estado de la Fase 0 (ROL A)

## Hecho

- Paquete `quantum1d` (`Grid`, `check_resolution`, `gaussian_packet`, `free`, `barrier`, `harmonic`, `evolve`, diagnósticos) según `docs/physics-spec.md` sección 12.
- `experiments/01_tunneling.py` y `experiments/02_euler_vs_splitstep.py`, con `outputs/01_tunneling.png`, `outputs/02_euler_vs_splitstep.png` y `outputs/report.md`.
- Verificación rápida del coder con scripts temporales fuera del repo (no son tests): AC1-AC17, AC19, AC21, AC22, AC24, AC25, AC27 pasan. Números: norma `|N-1|` < 5e-13 en AC11; reversibilidad 3.3e-14; convergencia en dt con pendiente 2.007; pico a pico de energía del oscilador 4.9999e-5.
- Experimento 01 (revisión 1: N=4096, dx=0.1, dt=0.005; `outputs/report.md`): T_measured = 0.2407910, T_packet = 0.2412562 (diferencia 4.7e-4), T_E0 = 0.2211382 (diferencia con la medida 0.0197); conv_diff_N = 3.5e-4 (N=8192, dt=0.00125), conv_diff_dt = 4.7e-6 (dt=0.0025); norm_error_max = 5.7e-12; energy_rel_dev_max = 3.2e-6; P_near_edges = 1.8e-13; P_inside_barrier = 4.0e-13; fase por paso 2.467 rad (base y convN), 1.234 rad (convdt).
- Experimento 02: `euler_norm_error_final` = 7.509e-4, `splitstep_norm_error_final` = 2.46e-12, `euler_norm_monotonic=True`.

## Falla

Nada conocido. AC18 pasa las 15 condiciones (el reporte imprime floats con `repr`, sin redondeo de impresión; `|T+R+P-1|` medido sobre el reporte es del orden de `norm_error_max` = 5.7e-12).

Resuelto en la revisión 1: las condiciones de energía y bordes de AC18 fallaban (0.15, 5e-6) por fase cinética por paso de 19.7 rad (N=8192, dt=0.01), que resonaba con modos de k alto. Con dt*kmax^2/2 <= pi (spec 6.4, decisión 0007) quedan en 3.2e-6 y 1.8e-13. Resuelto en la revisión 2: el formato `.10g` hacía fallar la condición 6 por redondeo de impresión.

## Lo que NO verifiqué

- AC20 (figura 01) y AC23 sólo por inspección visual / `git`; no hay test automático.
- No corrí pytest (no hay tests en esta rama; son del ROL B).
- Compatibilidad con Python 3.12: sólo corrí en 3.14.
- El efecto de la resonancia de Strang sólo lo probé con 3 corridas de diagnóstico (n y dt variados), no con un barrido sistemático.
- No instalé el paquete con hatchling/uv; sólo se usó `PYTHONPATH=src`.

## Mejoras no implementadas

- Guarda o advertencia de fase por paso en `evolve` (aviso si dt * k_max^2 / 2 > pi).
- Filtro de alto k o paso adaptativo para barreras de borde abrupto.
- Parámetros de CLI para los experimentos.
