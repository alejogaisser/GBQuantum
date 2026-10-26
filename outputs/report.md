# Reporte del experimento 01: efecto túnel

## Parámetros

| clave | valor |
|---|---|
| N | 4096 |
| L | 409.6 |
| dx | 0.1 |
| dt | 0.005 |
| steps | 20000 |
| t_final | 100.0 |
| x0 | -40.0 |
| k0 | 1.3 |
| sigma | 4.0 |
| V0 | 1.0 |
| a | 2.0 |
| barrier_center | 0.01 |
| n_barrier_cells | 20 |
| a_effective | 2.0 |
| kmax | 31.41592653589793 |
| guard_rhs | 8.2 |

## Energías

| clave | valor |
|---|---|
| E0 | 0.8450000000000001 |
| H0 | 0.8528125 |

## Probabilidades al tiempo final

| clave | valor |
|---|---|
| T_measured | 0.24079096785063403 |
| R_measured | 0.7592090321546687 |
| P_inside_barrier | 4.031401775726209e-13 |
| P_near_edges | 1.7989526657358672e-13 |

## Predicciones

| clave | valor |
|---|---|
| T_E0 | 0.22113815454967434 |
| T_packet | 0.24125620666135408 |
| abs_diff_measured_packet | 0.00046523881072005335 |
| abs_diff_measured_E0 | 0.019652813300959693 |

## Convergencia

| clave | valor |
|---|---|
| T_measured_convN | 0.24114147604423963 |
| T_measured_convdt | 0.24079567631685217 |
| convN_N | 8192 |
| convN_dt | 0.00125 |
| convdt_N | 4096 |
| convdt_dt | 0.0025 |
| conv_diff_N | 0.00035050819360560315 |
| conv_diff_dt | 4.708466218139584e-06 |

## Conservación

| clave | valor |
|---|---|
| norm_error_max | 5.7058802127585295e-12 |
| energy_rel_dev_max | 3.187376837641679e-06 |
| energy_rel_dev_final | 3.397795053744708e-14 |
| phase_per_step_base | 2.4674011002723395 |
| phase_per_step_convN | 2.4674011002723395 |
| phase_per_step_convdt | 1.2337005501361697 |
