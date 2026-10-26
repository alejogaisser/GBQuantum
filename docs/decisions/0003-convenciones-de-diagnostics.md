# 0003. Firmas y definiciones de diagnostics
- Date: 2026-10-26
- Status: accepted
- Decision: norm(psi, grid), mean_x(psi, grid), mean_p(psi, grid), std_x(psi, grid) y energy(psi, V, grid) devuelven float. norm = dx*sum(|psi|^2), sin raíz. Los valores esperados se dividen por la norma; <p> y la parte cinética de <H> se calculan con |fft(psi)|^2. std_x usa dos pasadas. Si la norma es 0, las funciones que dividen lanzan ValueError. Definiciones exactas en docs/physics-spec.md sección 7.
- Reason: el ROL B programa a ciegas contra la API; cualquier variante (raíz, sin normalizar) cambiaría resultados en tests con tolerancias de 1e-10.
