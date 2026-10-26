# Verificación independiente — Fase 0, rol B

## Predicciones registradas antes de importar el motor

Este registro se escribe desde la especificación pública, sin leer `src/`, los
experimentos, la especificación del constructor ni su rama. Un test pendiente no
es evidencia de que el motor pase. Se usan unidades naturales: ℏ = m = 1.

| Chequeo | Predicción / referencia independiente | Tolerancia prefijada y motivo | Fallo que detectaría |
| --- | --- | --- | --- |
| Grilla y potenciales | x_j = −L/2+jL/N; k = 2π·fftfreq; V libre = 0, barrera = altura en su interior, V armónico = ω²(x−c)²/2 | 2e−13 absoluto para aritmética de grilla/potenciales; pruebas sin puntos en los extremos de la barrera | Dominio cerrado, orden FFT incorrecto, factor físico equivocado |
| Paquete inicial | Densidad gaussiana con desviación σ, norma 1, ⟨x⟩ = x₀, ⟨p⟩ = k₀ | 2e−12 para norma y forma; 2e−10 para momentos: dominio amplio y resolución espectral hacen despreciables las colas | Confundir σ de amplitud con σ de densidad, signo de fase |
| Norma (a) | ∫\|ψ\|² dx = 1 tras 2000 pasos; V libre, barrera y armónico | Error < 1e−10; acumulación de redondeo FFT, no error de truncamiento temporal de un operador unitario | Normalización FFT, exponentes reales, pérdidas o ganancias de probabilidad |
| Movimiento libre (b) | ⟨x⟩ = x₀+k₀t, σ(t) = σ₀√(1+[t/(2σ₀²)]²) | 2e−8 absoluto; V = 0 elimina el error de separación y las colas quedan lejos del borde | Factor 1/2 cinético, signo de evolución, dispersión incorrecta |
| Base armónica (c) | Para ω = 1 y σ = √(1/2), la densidad permanece fija | Máximo error puntual < 2e−6 con dt = 0.002, t = 2; Strang aproxima H y produce una oscilación de tamaño O(dt²) | Anchura base o potencial incorrectos; no se exige estacionariedad exacta al separador |
| Centro armónico (c) | Paquete desplazado: ⟨x⟩ = x₀ cos(ωt) | 2e−5 absoluto en varios tiempos hasta t = 6 con dt = 0.002; margen sobre el desfase O(t·dt²) | Frecuencia, signo de fuerza o masa incorrectos |
| Convergencia (d) | Error L² complejo frente a una referencia temporal fina ∝ dt² | Pendiente log-log en [1.8, 2.2], R² > 0.995; dt = 0.08, 0.04, 0.02, 0.01 frente a dt = 0.000625 a t = 1.28 | Lie de orden 1, secuencia de Strang incorrecta, ausencia de convergencia |
| Energía (e) | ⟨H⟩ inicial libre = k₀²/2+1/(8σ²); base = ω/2; coherente desplazado = ω/2+ω²x₀²/2 | 2e−10 absoluto para valores iniciales espectrales | Confundir p con energía cinética o perder el factor 1/2 |
| Deriva de energía (e) | Energía casi constante en evolución armónica desplazada | Máxima deriva relativa < 2e−5 con dt = 0.005 hasta t = 10; margen sobre el error O(dt²) del Hamiltoniano modificado | Integración disipativa, fases equivocadas o diagnósticos inconsistentes |
| Túnel (f) | T medida ≈ Σ_{k>0}\|FFT(ψ₀)\|²T(k²/2) / Σ_{k>0}\|FFT(ψ₀)\|² | 0.01 absoluto; dx = 0.03125, barrera de ancho 1 con bordes entre muestras; error espacial y separación temporal en un potencial discontinuo | Comparar solo con T(E₀), barrera mal escalada o probabilidades incorrectas |
| Ventana de colisión (f) | A t = 40, reflexión y transmisión separadas; masa despreciable cerca de bordes | Masa central (±8) < 0.002 y exterior (\|x\| > 0.4L) < 1e−6 en todos los controles; presupuesto pequeño respecto del 1% de T | Medir durante colisión o después de envoltura periódica |
| Norma en colisión larga | Probabilidad total 1 tras los 16000 pasos del túnel | 1e−9 absoluto; presupuesto de redondeo mayor que en 2000 pasos, sin renormalizar | Pérdida acumulada de probabilidad durante dispersión |
| Euler (g) | Un modo plano libre crece como (1+[E·dt]²)^steps en norma cuadrada | 2e−12 relativo para la predicción y crecimiento > 20%; dt·E = 0.05, 100 pasos | Mostrar una comparación que renormaliza Euler y esconde su inestabilidad |
| Aliasing | kmax < 4(k₀+3/σ) debe rechazar el paquete | ValueError en resolución insuficiente; aceptación justo en igualdad y con margen cuando existe helper | Aliasing silencioso o desigualdad invertida |
| Contrato de evolución | Sin mutación del estado de entrada; salida nueva; un modo plano adquiere exp[−i(k²/2+V)t] | Identidad exacta de entrada; 2e−12 absoluto de fase; 2e−11 para reversibilidad | Cambios in place, fase global omitida, pasos o FFT mal aplicados |

Parámetros del túnel fijados: N = 8192, L = 256, x₀ = −35, k₀ = 2,
σ = 2, V₀ = 3, ancho = 1, centro = dx/2, dt = 0.0025, t final = 40.
El centro evita depender de la convención de inclusión de puntos en los bordes.
Se controla la masa exterior cada 2.5 unidades de tiempo y también al comienzo.
La separación prevista es aproximadamente ±45, muy lejos de ±128; además se
comprueban colas analíticas del paquete libre en toda la ventana.
Esto acota el riesgo de envoltura; no prueba ausencia matemática de
amplitudes infinitesimales ni observa cada paso de un potencial discontinuo.

## Fórmulas y decisiones

`tests/analytic.py` no importa el motor. Para E = V₀ se usa el límite
T = [1+V₀a²/2]⁻¹, para E = 0 con barrera positiva T = 0, y para
V₀ = 0 o a = 0 se toma T = 1. La fórmula subbarrera se evalúa en
logaritmos para evitar overflow en barreras opacas. El promedio espectral usa
únicamente k > 0; Δk y la convención de normalización FFT se cancelan en el cociente.

La energía no se conserva exactamente porque los operadores cinético y potencial
no conmutan: Strang conserva la norma por sus factores unitarios, pero aproxima la
evolución del Hamiltoniano con error global O(dt²). El caso libre debe ser exacto
salvo redondeo. La prueba de orden compara funciones complejas completas, sin
descartar la fase global, sobre una grilla compartida bien resuelta.

Ambigüedades explícitas: el contrato no fija inclusión de los extremos de la
barrera, ni si `norm` devuelve ∫\|ψ\|² o su raíz. Las pruebas con norma 1 sirven
para ambas convenciones; en el resto se integra la probabilidad directamente.
`potentials` y `diagnostics` se acceden desde el paquete público `quantum1d`.
La guarda se prueba mediante `quantum1d.check_resolution` si existe; en caso
contrario `evolve` debe rechazar un estado construido fuera de resolución.
La guarda especificada usa k₀ con signo; los tests del umbral usan incidencia
positiva y no presuponen que el constructor sustituya k₀ por |k₀|.
No se presupone validación adicional de argumentos fuera del contrato.

## Ejecución y resultados

2026-10-26: uv 0.13.5, Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1,
pytest 9.1.1, en entorno temporal externo al repositorio. Las predicciones
se registraron primero en el commit `98284fa`; no se ajustó ninguna tolerancia.
El ejecutable usado fue `C:\Users\ASUS\AppData\Local\Temp\gbquantum-auditor-uv\bin\uv.exe`;
los comandos siguientes lo abrevian como `uv` (no se instaló globalmente).

```powershell
uv run --no-project --python 3.12 --with numpy --with scipy --with pytest pytest -p no:cacheprovider tests/test_analytic_selfcheck.py
uv run --no-project --python 3.12 --with numpy --with scipy --with pytest pytest -p no:cacheprovider tests -q --tb=line
```

Resultado final de la suite: **29 pasaron, 39 errores de preparación**, todos
`ModuleNotFoundError: No module named 'quantum1d'`. Ningún test del motor alcanzó
sus aserciones. Los 29 aprobados son 27 autochecks analíticos, la demostración
negativa de Euler y el margen espacial analítico del túnel. No se usa `skip` para
ocultar el paquete ausente. La suite completa recolecta 68 casos.

Predicciones numéricas independientes: T(E₀) = 0.19184490127131365,
T_paquete = 0.20723017884137554; la diferencia supera un punto porcentual y
justifica usar el promedio. Euler deja probabilidad 1.283624888738461.
La transmisión medida queda pendiente. También se validó la fórmula de barrera
contra un sistema lineal de continuidad de ψ y ψ′, por debajo y encima del umbral.

Tras integrar e instalar el paquete del constructor: `uv run pytest tests -q`.
Un fallo debe investigarse y reportarse; esta ejecución no autoriza relajar las
tolerancias. No se copiará ni se leerá implementación para hacerlas pasar.

## Lo que NO verifiqué

- El motor del rol A, su documentación, sus experimentos y sus salidas.
- El port futuro a TypeScript, la interfaz y cualquier física fuera de 1D sin spin.
- Valores de tolerancia fuera de los parámetros concretos declarados arriba.
