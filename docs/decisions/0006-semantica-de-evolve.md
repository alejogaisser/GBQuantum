# 0006. Semántica de evolve y gaussian_packet
- Date: 2026-10-26
- Status: accepted
- Decision: evolve copia psi a complex128, nunca modifica psi ni V, no renormaliza y aplica en cada paso V/2, T, V/2 sin fusionar medios pasos. Acepta dt real finito (negativo = hacia atrás) y steps entero >= 0 (0 devuelve una copia); lanza ValueError si las formas no son (grid.n,), si V no es real, si dt no es finito o si steps es inválido. gaussian_packet devuelve la fórmula exacta, sin renormalizar.
- Reason: renormalizar escondería errores de norma. Strang es simétrico, así que dt negativo da la inversa exacta, útil para probar reversibilidad. Con la fórmula exacta, B puede comparar psi0 sin ambigüedad.
