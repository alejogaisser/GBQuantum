# 0007. Resolución temporal de fase: dt*(pi/dx)**2/2 <= pi en el experimento de túnel
- Date: 2026-10-26
- Status: accepted
- Decision: toda corrida del experimento 01 cumple dt*(pi/dx)**2/2 <= pi y lo registra en outputs/report.md (claves phase_per_step_*). La base pasa a N=4096 (dx=0.1) con dt=0.005; la convergencia usa N=8192 con dt=0.00125 y N=4096 con dt=0.0025. evolve no verifica esta condición; una guarda en evolve queda como mejora no implementada.
- Reason: si algún modo de la grilla tiene fase cinética por paso >= 2π, la barrera lo excita en resonancia numérica. La norma sigue exacta, pero <H> subió 15 % y llegó radiación a los bordes (con N=8192, dt=0.01, fase 19.7). El barrido medido pasa con fase 5.9 y falla con 6.9, consistente con el umbral 2π; π da margen ×2.
