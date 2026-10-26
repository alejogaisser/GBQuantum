# Reporte del experimento 01: efecto túnel

## Parámetros

| clave | valor |
|---|---|
| N | 4096 |
| L | 409.6 |
| dx | 0.1 |
| dt | 0.005 |
| steps | 20000 |
| t_final | 100 |
| x0 | -40 |
| k0 | 1.3 |
| sigma | 4 |
| V0 | 1 |
| a | 2 |
| barrier_center | 0.01 |
| n_barrier_cells | 20 |
| a_effective | 2 |
| kmax | 31.41592654 |
| guard_rhs | 8.2 |

## Energías

| clave | valor |
|---|---|
| E0 | 0.845 |
| H0 | 0.8528125 |

## Probabilidades al tiempo final

| clave | valor |
|---|---|
| T_measured | 0.2407909679 |
| R_measured | 0.7592090322 |
| P_inside_barrier | 4.031401776e-13 |
| P_near_edges | 1.798952666e-13 |

## Predicciones

| clave | valor |
|---|---|
| T_E0 | 0.2211381545 |
| T_packet | 0.2412562067 |
| abs_diff_measured_packet | 0.0004652388107 |
| abs_diff_measured_E0 | 0.0196528133 |

## Convergencia

| clave | valor |
|---|---|
| T_measured_convN | 0.241141476 |
| T_measured_convdt | 0.2407956763 |
| convN_N | 8192 |
| convN_dt | 0.00125 |
| convdt_N | 4096 |
| convdt_dt | 0.0025 |
| conv_diff_N | 0.0003505081936 |
| conv_diff_dt | 4.708466218e-06 |

## Conservación

| clave | valor |
|---|---|
| norm_error_max | 5.705880213e-12 |
| energy_rel_dev_max | 3.187376838e-06 |
| energy_rel_dev_final | 3.397795054e-14 |
| phase_per_step_base | 2.4674011 |
| phase_per_step_convN | 2.4674011 |
| phase_per_step_convdt | 1.23370055 |
