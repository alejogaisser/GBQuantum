# Reporte del experimento 01: efecto túnel

## Parámetros

| clave | valor |
|---|---|
| N | 8192 |
| L | 409.6 |
| dx | 0.05 |
| dt | 0.01 |
| steps | 10000 |
| t_final | 100 |
| x0 | -40 |
| k0 | 1.3 |
| sigma | 4 |
| V0 | 1 |
| a | 2 |
| barrier_center | 0.01 |
| n_barrier_cells | 40 |
| a_effective | 2 |
| kmax | 62.83185307 |
| guard_rhs | 8.2 |

## Energías

| clave | valor |
|---|---|
| E0 | 0.845 |
| H0 | 0.8528125 |

## Probabilidades al tiempo final

| clave | valor |
|---|---|
| T_measured | 0.2412213765 |
| R_measured | 0.7587774871 |
| P_inside_barrier | 1.136417414e-06 |
| P_near_edges | 5.431160617e-06 |

## Predicciones

| clave | valor |
|---|---|
| T_E0 | 0.2211381545 |
| T_packet | 0.2412562067 |
| abs_diff_measured_packet | 3.48302061e-05 |
| abs_diff_measured_E0 | 0.02008322191 |

## Convergencia

| clave | valor |
|---|---|
| T_measured_N16384 | 0.2412912035 |
| T_measured_dt0005 | 0.241161228 |
| conv_diff_N | 6.982707537e-05 |
| conv_diff_dt | 6.014850376e-05 |

## Conservación

| clave | valor |
|---|---|
| norm_error_max | 3.269162718e-12 |
| energy_rel_dev_max | 0.1504681634 |
| energy_rel_dev_final | 0.1500554153 |
