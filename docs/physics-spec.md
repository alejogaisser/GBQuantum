# GBQuantum Fase 0: especificación física del motor `quantum1d`

Estado: final (ROL A, arquitecto). Fecha: 2026-10-26.

Este documento es la referencia común del ROL A (implementa `src/quantum1d/`) y del ROL B
(verifica a ciegas en `tests/`). La sección 12 resume las convenciones de la API y es vinculante:
si algo de este documento contradice a la API del contrato del usuario, manda el contrato y este
documento tiene un error.

---

## 1. Unidades naturales: $\hbar = m = 1$

Todo el motor trabaja en unidades donde la constante de Planck reducida y la masa de la partícula
valen 1. Consecuencias directas:

| Magnitud | Fórmula general | Con $\hbar = m = 1$ |
|---|---|---|
| Momento | $p = \hbar k$ | $p = k$ |
| Energía cinética de una onda plana | $E = \hbar^2 k^2 / 2m$ | $E = k^2/2$ |
| Velocidad de grupo | $v = \hbar k / m$ | $v = k$ |
| Longitud de onda de de Broglie | $\lambda = 2\pi/k$ | igual |

Una sola unidad de longitud $\ell$ fija las demás: tiempo $m\ell^2/\hbar$, energía $\hbar^2/(m\ell^2)$.
Ejemplo: un electrón con $\ell = 1\,\text{nm}$ tiene unidad de energía $\approx 76.2\,\text{meV}$ y
unidad de tiempo $\approx 8.64\,\text{fs}$. Todas las cantidades son `float64` / `complex128`.

## 2. La ecuación que se resuelve

$$
i\,\frac{\partial \psi(x,t)}{\partial t} = \hat H\,\psi(x,t),
\qquad
\hat H = \hat T + \hat V = -\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x).
$$

- $\psi(x,t)$: amplitud de probabilidad compleja (unidades de $\text{longitud}^{-1/2}$).
- $|\psi(x,t)|^2$: densidad de probabilidad de encontrar la partícula en $x$.
- $V(x)$: energía potencial real e independiente del tiempo.

Invariantes del problema continuo (los que el motor debe respetar o aproximar):
1. **Norma** $\int |\psi|^2\,dx = 1$ en todo tiempo (la evolución es unitaria).
2. **Energía** $\langle H\rangle$ constante (porque $V$ no depende de $t$).

## 3. Grilla periódica y FFT (`Grid`)

`Grid(n, length)` discretiza el dominio $[-L/2,\,L/2)$ con $n$ puntos:

$$
\Delta x = \frac{L}{n},\qquad x_j = -\frac{L}{2} + j\,\Delta x,\quad j = 0,\dots,n-1
\quad(\text{el punto } x = +L/2 \text{ no está: es el mismo que } -L/2).
$$

Números de onda en el **orden de la FFT** (`k = 2*pi*numpy.fft.fftfreq(n, d=dx)`):
$k = [0, \Delta k, \dots, (n/2-1)\Delta k, -\tfrac{n}{2}\Delta k, \dots, -\Delta k]$ para $n$ par, con
$\Delta k = 2\pi/L$. El número de onda máximo representable (Nyquist) es $k_{\max} = \pi/\Delta x$.

**El dominio es un anillo**: lo que sale por la derecha entra por la izquierda. Toda simulación
debe mantener $|\psi|^2$ despreciable cerca de $\pm L/2$.

Convención de FFT (la de NumPy): $\tilde\psi_m = \sum_j \psi_j\, e^{-2\pi i\, jm/n}$. Como la grilla
empieza en $-L/2$ y no en 0, $\tilde\psi_m$ difiere de la transformada continua por una fase
$e^{i k_m L/2}$ y un factor $1/\Delta x$; ninguno importa para $|\tilde\psi_m|^2$ ni para multiplicar
por una función de $k$. Teorema de Parseval discreto: $\sum_m |\tilde\psi_m|^2 = n \sum_j |\psi_j|^2$.
La derivada es espectral: $\partial_x \leftrightarrow$ multiplicar por $i k_m$.

## 4. Paquete gaussiano (`gaussian_packet`)

$$
\psi_0(x) = (2\pi\sigma^2)^{-1/4}\, \exp\!\left[-\frac{(x-x_0)^2}{4\sigma^2}\right] e^{i k_0 x}
$$

- $x_0$: posición media $\langle x\rangle$. $k_0$: momento medio $\langle p\rangle$ (= velocidad).
- $\sigma$: desviación estándar de $|\psi_0|^2 = (2\pi\sigma^2)^{-1/2} e^{-(x-x_0)^2/(2\sigma^2)}$.
- En momento: $|\phi(k)|^2 \propto \exp[-(k-k_0)^2/(2\sigma_k^2)]$ con $\sigma_k = 1/(2\sigma)$.
  El producto $\sigma\,\sigma_k = 1/2$ es el mínimo que permite Heisenberg.
- Energía de un paquete libre: $\langle H\rangle = k_0^2/2 + 1/(8\sigma^2)$
  (porque $\langle p^2\rangle = k_0^2 + \sigma_k^2$).

Implementación: se evalúa **exactamente esa fórmula** en los $x_j$, sin renormalizar
numéricamente y sin periodizar. La norma discreta vale 1 salvo la masa de la cola que cae fuera
de $[-L/2, L/2)$, que es exponencialmente chica si $x_0$ está a muchas $\sigma$ de los bordes.
Antes de construir $\psi$, `gaussian_packet` llama a `check_resolution(grid, k0, sigma)` (sección 8).

## 5. Potenciales (`potentials`)

- `free(grid)`: $V = 0$.
- `harmonic(grid, omega, center=0.0)`: $V = \tfrac12 \omega^2 (x-c)^2$ (oscilador de masa 1).
  En el anillo, $V$ tiene un quiebre en $\pm L/2$; es inofensivo si $\psi \approx 0$ ahí.
- `barrier(grid, center, width, height)`: $V = V_0$ (= `height`) en los nodos con
  $|x_j - c| \le w/2$ (intervalo **cerrado**) y $V = 0$ en el resto.

**Efecto de discretización de la barrera.** En la grilla la barrera ocupa
$n_{\text{celdas}}$ nodos y su ancho efectivo es $a_{\text{ef}} = n_{\text{celdas}}\,\Delta x$
(cada nodo representa una celda de ancho $\Delta x$). Si $w/\Delta x$ es entero y ningún borde
$c \pm w/2$ cae sobre un nodo, entonces $n_{\text{celdas}} = w/\Delta x$ exactamente y
$a_{\text{ef}} = w$. Si un borde cae justo sobre un nodo (por ejemplo $c=0$, $w=2$, $\Delta x=0.05$),
la inclusión de ese nodo la decide el redondeo de punto flotante y el ancho puede variar en
$\pm\Delta x$; como $T \sim e^{-2\kappa a}$, eso cambia $T$ en $\sim 2\kappa\Delta x$ relativo
($\approx 6\,\%$ con $\kappa = 0.56$, $\Delta x = 0.05$, o sea $\approx 0.013$ absoluto: más que la tolerancia).
Regla: (1) elegir $c$ y $w$ para que ningún borde caiga en un nodo; (2) al comparar con la $T$
analítica, usar $a = a_{\text{ef}} = $ `count_nonzero(V) * dx`. Además, un escalón muestreado solo
está resuelto con precisión $O(\Delta x)$; el error residual se mide con convergencia en $n$.

## 6. Evolución: split-step de Strang (`evolve`)

Un paso de tamaño $\Delta t$:

$$
\psi \;\leftarrow\; e^{-iV\Delta t/2}\;\mathcal F^{-1}\!\left[e^{-i k^2 \Delta t/2}\;\mathcal F\!\left[e^{-iV\Delta t/2}\,\psi\right]\right]
$$

1. medio paso de potencial en $x$ (multiplicar por la fase $e^{-iV_j\Delta t/2}$);
2. paso completo de energía cinética en $k$ (FFT, multiplicar por $e^{-ik_m^2\Delta t/2}$, FFT inversa);
3. otro medio paso de potencial en $x$.

Contrato de comportamiento: devuelve un array nuevo `complex128`; **no modifica** `psi` ni `V`;
**no renormaliza** nunca (renormalizar escondería errores); `steps = 0` devuelve una copia;
`dt` negativo evoluciona hacia atrás. Llamar `evolve` dos veces con `s1` y `s2` pasos da lo mismo
que una vez con `s1 + s2` (los medios pasos no se fusionan entre pasos).

### 6.1 Por qué conserva la norma exactamente

Cada factor es unitario: multiplicar por $e^{-i\theta_j}$ con $\theta_j$ real no cambia $|\psi_j|$, y la
FFT es unitaria salvo el factor $\sqrt n$ que la inversa deshace (Parseval). Entonces
$\Delta x\sum_j|\psi_j|^2$ se conserva **exactamente en aritmética exacta, para cualquier $V$ real y
cualquier $\Delta t$**. En punto flotante solo queda redondeo: $\sim 10^{-16}$ por operación, que tras
2000 pasos acumula $\lesssim 10^{-13}$. Por eso la tolerancia $|N-1| < 10^{-10}$ tiene un margen
$\gtrsim 10^3$ y una falla indicaría un bug, no ruido.

### 6.2 Por qué es de orden 2 (Strang + BCH)

Con $X = -iV\Delta t$ e $Y = -iT\Delta t$, la fórmula simétrica de Baker-Campbell-Hausdorff da

$$
e^{X/2}e^{Y}e^{X/2} = \exp\!\Big(X + Y - \tfrac{1}{24}[X,[X,Y]] + \tfrac{1}{12}[Y,[Y,X]] + O(\Delta t^5)\Big).
$$

Por la simetría no aparece término en $\Delta t^2$: el error local es $O(\Delta t^3)$ y, sumado sobre
$t/\Delta t$ pasos, el error global es $O(\Delta t^2)$. Al dividir $\Delta t$ por 2 el error cae $\approx 4$
veces; en log-log, error vs $\Delta t$ tiene pendiente 2.

**Caso especial importante:** si $[T,V]=0$ (con $V=0$ o constante) los conmutadores se anulan y el
método es **exacto** para cualquier $\Delta t$ (salvo redondeo). Por eso una prueba de convergencia
en $\Delta t$ con potencial libre solo mide redondeo; hay que usar un $V$ con $V'\ne 0$ donde vive $\psi$
(por ejemplo el oscilador armónico). Por la misma razón, la evolución libre (criterio b) reproduce
la solución analítica hasta redondeo y truncamiento.

### 6.3 Por qué $\langle H\rangle$ no se conserva exacto (pero no deriva)

La fórmula de 6.2 dice que un paso es la evolución **exacta** de un Hamiltoniano modificado

$$
\tilde H = H + \Delta t^2\left(\tfrac{1}{24}[V,[V,T]] - \tfrac{1}{12}[T,[T,V]]\right) + O(\Delta t^4).
$$

El método conserva $\langle\tilde H\rangle$, no $\langle H\rangle$. Entonces
$\langle H\rangle(t) = \langle \tilde H\rangle - \Delta t^2\langle C\rangle(t)$: oscila con amplitud $O(\Delta t^2)$
**acotada y sin deriva secular** (no crece con $t$). Con $V=0$, $\tilde H = H$. Advertencia honesta:
en una grilla con $\Delta t\,k_{\max}^2/2 \gg 1$ la serie es formal; el argumento vale para la parte del
espectro que el paquete ocupa. La comprobación práctica es medir que no hay deriva.

**Predicción exacta para el oscilador** ($V=\omega^2x^2/2$, todo cuadrático, Ehrenfest es exacto): los
momentos evolucionan como el integrador de Verlet clásico, que conserva exactamente
$J = \tfrac12\langle p^2\rangle + \tfrac12(1-\omega^2\Delta t^2/4)\,\omega^2\langle x^2\rangle$. De ahí:
- $\langle H\rangle(t) = J + \tfrac{\omega^4\Delta t^2}{8}\langle x^2\rangle(t)$. Para un estado desplazado a
  $x_0$ (con $k_0=0$) la oscilación pico a pico de $\langle H\rangle$ es $\omega^4\Delta t^2 x_0^2/8$
  ($\omega=1$, $x_0=2$, $\Delta t=0.01$: $5.0\times10^{-5}$, relativo a $\langle H\rangle=2.5$: $2\times10^{-5}$).
- Frecuencia numérica $\tilde\omega$ con $\cos(\tilde\omega\Delta t) = 1-\omega^2\Delta t^2/2$, o sea
  $\tilde\omega \approx \omega(1+\omega^2\Delta t^2/24)$: $\langle x\rangle(t) = x_0\cos(\tilde\omega t)$ y
  $|\langle x\rangle - x_0\cos\omega t| \le x_0\,\omega^3\Delta t^2\, t/24$ ($x_0=2$, $t=20$, $\Delta t=0.01$: $1.7\times10^{-4}$).
- El estado base analítico ($\sigma_0 = 1/\sqrt{2\omega}$) no es exactamente estacionario para el
  método: `std_x` oscila entre $\sigma_0$ y $\sigma_0(1-\omega^2\Delta t^2/4)^{-1/2} \approx \sigma_0(1+\omega^2\Delta t^2/8)$
  (desvío relativo máximo $1.25\times10^{-5}$ con $\Delta t = 0.01$).

## 7. Diagnósticos (`diagnostics`)

Con $N = \Delta x\sum_j|\psi_j|^2$ y $\tilde\psi = $ `fft(psi)`:

- `norm(psi, grid)`: $N = \Delta x\sum_j |\psi_j|^2$ (la integral, **sin raíz cuadrada**).
- `mean_x(psi, grid)`: $\langle x\rangle = \Delta x\sum_j x_j|\psi_j|^2 / N$.
- `std_x(psi, grid)`: $\sqrt{\Delta x\sum_j (x_j-\langle x\rangle)^2|\psi_j|^2 / N}$ (dos pasadas, sin restar $\langle x\rangle^2$).
- `mean_p(psi, grid)`: $\langle p\rangle = \sum_m k_m|\tilde\psi_m|^2 / \sum_m|\tilde\psi_m|^2$.
- `energy(psi, V, grid)`: $\langle H\rangle = \sum_m \tfrac{k_m^2}{2}|\tilde\psi_m|^2/\sum_m|\tilde\psi_m|^2 \;+\; \Delta x\sum_j V_j|\psi_j|^2/N$.

Todas devuelven `float`. Los valores esperados se dividen por la norma (para un $\psi$ normalizado
es lo mismo); si $N = 0$, `mean_x`, `std_x`, `mean_p` y `energy` lanzan `ValueError`. `norm` nunca lanza.
$\langle x\rangle$ usa la coordenada de la grilla en $[-L/2,L/2)$: solo tiene sentido si el paquete está
lejos de los bordes. El nodo de Nyquist cuenta con su $k$ negativo (irrelevante si está resuelto).

## 8. Guarda de resolución y aliasing (`check_resolution`)

En una grilla con paso $\Delta x$, los números de onda $k$ y $k + 2\pi/\Delta x$ son indistinguibles. Si el
paquete tiene contenido más allá de $k_{\max} = \pi/\Delta x$, ese contenido se "pliega" a otro $k$ con otra
velocidad: el resultado es incorrecto **sin ningún aviso** (aliasing). La guarda es

$$
\text{ValueError si}\quad \frac{\pi}{\Delta x} < 4\left(|k_0| + \frac{3}{\sigma}\right).
$$

- $k_c = |k_0| + 3/\sigma = |k_0| + 6\sigma_k$: fuera de ese rango $|\phi(k)|^2 < e^{-18}\approx 1.5\times10^{-8}$ del pico.
- Factor 4: $\pi/\Delta x \ge 4k_c \iff \lambda_c = 2\pi/k_c \ge 8\Delta x$, al menos 8 puntos por la longitud de onda
  más corta relevante, y margen para momento que agregue el potencial.
- Se usa $|k_0|$: el aliasing depende del módulo; un paquete con $k_0<0$ se mueve a la izquierda y
  exige la misma resolución que con $+|k_0|$. La igualdad pasa. $\sigma \le 0$ lanza `ValueError`.
- Dónde vive: `check_resolution(grid, k0, sigma)` en `quantum1d/grid.py`, importable como
  `from quantum1d import check_resolution`; devuelve `None` si pasa. `gaussian_packet` la llama.
  **`evolve` no tiene esta guarda** (no conoce $k_0$ ni $\sigma$; recibe un $\psi$ arbitrario) y `Grid`
  tampoco (no conoce el paquete).
- Límite: solo cubre el paquete inicial. Si el potencial acelera la partícula, hay que aplicar la
  guarda con el momento máximo de la trayectoria; en el oscilador
  $|k|_{\max} \le \sqrt{k_0^2 + \omega^2(x_0-c)^2}$ (ejemplo: $x_0=2$, $\omega=1$, $\sigma_0=1/\sqrt2$ exige
  $\pi/\Delta x \ge 4(2+4.243) = 24.97$).

## 9. Física a verificar (criterios a-g)

| | Qué | Predicción | Por qué la tolerancia |
|---|---|---|---|
| a | Norma, 2000 pasos, cualquier $V$ | $\lvert N-1\rvert < 10^{-10}$ | 6.1: solo redondeo $\lesssim10^{-13}$ |
| b | Paquete libre | $\sigma(t) = \sigma_0\sqrt{1+(t/(2\sigma_0^2))^2}$, $\langle x\rangle = x_0+k_0t$, $\langle p\rangle = k_0$ | 6.2: con $V=0$ el método es exacto; solo redondeo y truncamiento |
| c | Oscilador $\omega=1$ | estado base estacionario; desplazado: $\langle x\rangle = x_0\cos\omega t$ | 6.3: desvíos $O(\Delta t^2)$ con cotas explícitas |
| d | Convergencia en $\Delta t$ | error vs referencia fina cae $\approx4\times$ al dividir $\Delta t$; pendiente 2 | 6.2; requiere $V'\ne0$ |
| e | Energía | oscilación $O(\Delta t^2)$ acotada, sin deriva | 6.3 |
| f | Túnel | transmisión medida $\approx T_{\text{paquete}}$ dentro de 0.01 absoluto | sección 10 |
| g | Euler explícito | la norma crece | sección 11.2 |

## 10. Efecto túnel

Barrera rectangular de alto $V_0$ y ancho $a$; onda plana de energía $E = k^2/2$ que llega desde la
izquierda. Coeficiente de transmisión (probabilidad de cruzar):

$$
E<V_0:\quad T = \left[1 + \frac{V_0^2\sinh^2(\kappa a)}{4E(V_0-E)}\right]^{-1},\quad \kappa = \sqrt{2(V_0-E)}
$$
$$
E>V_0:\quad T = \left[1 + \frac{V_0^2\sin^2(q a)}{4E(E-V_0)}\right]^{-1},\quad q = \sqrt{2(E-V_0)}
$$
$$
E = V_0 \;(\text{límite de ambas}):\quad T = \left[1 + \frac{V_0 a^2}{2}\right]^{-1}
$$

$\kappa$ es la tasa de decaimiento de la onda evanescente dentro de la barrera; $q$, el número de onda
dentro de ella cuando pasa "por encima". Debajo de la barrera, $T \sim e^{-2\kappa a}$: muy sensible a $a$.

**Un paquete no es una onda plana.** Es una superposición de ondas planas con peso $|\phi(k)|^2$; por
linealidad, cada componente se transmite con su propia $T(E(k))$. La predicción para el paquete es

$$
T_{\text{paquete}} = \frac{\sum_{k_m>0} |\phi_m|^2\, T(k_m^2/2)}{\sum_{k_m>0}|\phi_m|^2},\qquad \phi = \text{FFT}(\psi_0).
$$

**Por qué $T_{\text{paquete}} \ne T(E_0)$** con $E_0 = k_0^2/2$: $T(E)$ no es lineal (es convexa y creciente
cerca de $E_0$), y el promedio de una función convexa supera a la función del promedio (Jensen). Además
la parte del paquete con $E > V_0$ pasa "por encima" con $T$ mucho mayor. Las componentes más rápidas
se transmiten más: el paquete transmitido tiene $\langle k\rangle > k_0$ (la barrera filtra momentos).

**Medición** (al tiempo final $t_f$, terminada la colisión), con $c$ el centro y $a$ el ancho:
$T_{\text{med}} = \Delta x\sum_{x_j > c+a/2}|\psi_j|^2$, $R_{\text{med}} = \Delta x\sum_{x_j < c-a/2}|\psi_j|^2$, y el resto
$P_{\text{dentro}}$ queda en la barrera. Validez: el paquete inicial no toca la barrera, nada llega a los
bordes del anillo, y las componentes lentas que se transmitirían ya cruzaron.

**Tolerancia 0.01 absoluta.** Fuentes de error: bordes de barrera discretizados $O(\Delta x)$ y división
temporal $O(\Delta t^2)$ (ambas se miden por convergencia y deben quedar $<10^{-3}$), componentes que no
terminaron de cruzar ($<10^{-6}$) y probabilidad que llegue a los bordes ($<10^{-10}$). Total esperado
$\lesssim 2\times10^{-3}$: 0.01 deja margen $\times5$ y sigue siendo menor que $|T_{\text{paquete}}-T(E_0)| \approx 0.02$
del experimento 01, así que la prueba distingue entre ambas predicciones.

## 11. Experimentos

### 11.1 `experiments/01_tunneling.py`

Parámetros: $n = 8192$, $L = 409.6$ ($\Delta x = 0.05$), $\Delta t = 0.01$, $t_f = 100$ (10000 pasos),
$x_0 = -40$, $k_0 = 1.3$, $\sigma = 4$, $V_0 = 1$, $a = 2$, centro de barrera $c = 0.01$.

Cuentas:
- Guarda: $\pi/\Delta x = 62.83 \ge 4(1.3 + 3/4) = 8.2$. Se usa $\Delta x$ mucho menor que lo que exige la
  guarda para resolver los bordes de la barrera (40 celdas en el ancho).
- Barrera: bordes en $-0.99$ y $1.01$, a $0.2\Delta x$ del nodo más cercano; $n_{\text{celdas}} = 40$,
  $a_{\text{ef}} = 2.0$. Lo mismo vale con $\Delta x = 0.025$ (80 celdas), usado en la convergencia.
- Energía: $E_0 = 0.845 = 0.845\,V_0$; $\langle H\rangle(0) = 0.845 + 1/128 = 0.8528$; $\sigma_k = 0.125$.
  $T(E_0) = 0.2211$. $T_{\text{paquete}} \approx 0.241$ (cuadratura a mano; el script da el valor exacto).
  Un $18\,\%$ de $|\phi|^2$ tiene $E > V_0$ ($k > \sqrt2$).
- Sin solapamiento inicial: el borde izquierdo de la barrera está a $39\,\text{u} = 9.75\sigma$ de $x_0$
  ($P < 10^{-20}$).
- Colisión: el centro llega a la barrera en $t \approx 40/1.3 = 30.8$. Una componente de momento $k$
  transmitida está en $x \approx x_0 + kt$; a $t_f = 100$ ya cruzó si $k > 0.41$, o sea todo salvo
  $k_0 - 7.1\sigma_k$ (peso $\sim10^{-11}$, y con $T(k) \sim 10^{-3}$).
- Bordes: $\sigma(100) = 4\sqrt{1+(100/32)^2} = 13.1$; centro libre equivalente en $x_0+k_0t_f = 90$.
  Para tocar $\pm204.8$ hace falta $k \ge 2.45 = k_0 + 9.2\sigma_k$ (peso $\sim10^{-18}$). Lo reflejado
  está en $\approx -(x_0 + kt)$, simétrico.
- Convergencia que corre el script: $T_{\text{med}}$ con $n = 16384$ (mismo $L$) y con $\Delta t = 0.005$;
  ambas diferencias contra la corrida base deben ser $<10^{-3}$.

**Qué deberías ver y por qué.**
1. $t=0$: una campana centrada en $-40$; el color (fase) da una vuelta completa cada
   $\lambda_0 = 2\pi/k_0 \approx 4.8$: eso es el momento $k_0$ "pintado".
2. $t = 30$ (en plena colisión): a la izquierda de la barrera aparecen franjas de interferencia
   entre la onda incidente y la reflejada, con período $\pi/k_0 \approx 2.4$; dentro de la barrera,
   $|\psi|^2$ decae (onda evanescente, ritmo $2\kappa$).
3. $t = 50$: dos paquetes que se separan: uno transmitido (chico) y uno reflejado (grande).
4. $t = 100$: transmitido $\approx 24\,\%$ a la derecha, reflejado $\approx 76\,\%$ a la izquierda, ambos más
   anchos (dispersión). El transmitido tiene franjas de color algo más apretadas ($\langle k\rangle > k_0$).
5. Error de norma: plano en $\sim10^{-14}$-$10^{-13}$ (solo redondeo, 6.1).
6. Error de energía: un "bulto" mientras el paquete toca la barrera (los conmutadores de 6.3 son
   grandes donde $V$ salta) de orden $10^{-5}$-$10^{-4}$ relativo, que vuelve casi a cero después.
7. Barras: $T_{\text{med}} \approx T_{\text{paquete}} \approx 0.24$, claramente distinta de $T(E_0) \approx 0.22$.

### 11.2 `experiments/02_euler_vs_splitstep.py`

Euler explícito: $\psi_{n+1} = \psi_n - i\Delta t\,H\psi_n$. Entonces
$\|\psi_{n+1}\|^2 = \|\psi_n\|^2 + \Delta t^2\|H\psi_n\|^2 > \|\psi_n\|^2$: la norma crece en **cada** paso, para
cualquier $\Delta t$. Cada autoestado de energía $E$ se amplifica por $\sqrt{1+E^2\Delta t^2} > 1$: Euler es
incondicionalmente inestable para Schrödinger. Para tiempos cortos $N(t) \approx \exp(\Delta t\,\langle H^2\rangle\, t)$.

Parámetros: $n = 128$, $L = 18$ ($\Delta x = 0.1406$, $k_{\max} = 22.34$), $V$ = oscilador $\omega = 1$,
paquete $x_0 = 1$, $k_0 = 0$, $\sigma = 1/\sqrt2$ (estado coherente), $\Delta t = 2.5\times10^{-4}$,
8000 pasos ($t_f = 2$). Guarda a lo largo de la trayectoria: $4(1 + 4.243) = 20.97 \le 22.34$.
$\langle H\rangle = 1$, $\langle H^2\rangle = 1.5$, así que Euler debe llegar a $|N-1| \approx 7.5\times10^{-4}$ en $t=2$.
Los modos de la grilla con $E \lesssim 290$ se amplifican $e^{42}$ en total, partiendo de redondeo
($\sim10^{-32}$ en $|\cdot|^2$): $\sim10^{-14}$, no contaminan. Con $\Delta t$ 4 veces mayor sí explotarían.

**Qué deberías ver:** en escala log, la curva de Euler sube casi en línea recta (crecimiento
proporcional a $t$) hasta $\sim10^{-3}$; la del split-step queda en $\sim10^{-15}$-$10^{-13}$, once órdenes abajo.

## 12. Convenciones de la API (vinculante para A y B)

Todo se importa desde `quantum1d`: `Grid`, `check_resolution`, `gaussian_packet`, `evolve`, los módulos
`potentials` y `diagnostics`, y también `free`, `barrier`, `harmonic`, `norm`, `mean_x`, `mean_p`, `std_x`,
`energy` directamente. `quantum1d.evolve` es la **función** (el módulo es `quantum1d.evolve` solo vía
`from quantum1d.evolve import evolve`).

- `Grid(n, length)`: `n` entero $\ge 2$ y `length` finito y positivo; si no, `ValueError`. Atributos `n` (int),
  `length` (float), `dx = length/n`, `x = -length/2 + dx*arange(n)`, `k = 2*pi*fftfreq(n, d=dx)`.
- `check_resolution(grid, k0, sigma) -> None`: sección 8.
- `gaussian_packet(grid, x0, k0, sigma) -> ndarray complex128`: sección 4; `ValueError` por guarda o $\sigma \le 0$.
- `free(grid)`, `barrier(grid, center, width, height)`, `harmonic(grid, omega, center=0.0)` -> `ndarray float64`.
  `ValueError` si `width <= 0` o `omega <= 0`.
- `evolve(psi, V, grid, dt, steps) -> ndarray complex128`: sección 6. `ValueError` si `psi` o `V` no tienen
  forma `(grid.n,)`, si `V` no es real, si `dt` no es un real finito o si `steps` no es un entero $\ge 0$.
- `norm`, `mean_x`, `mean_p`, `std_x` con firma `(psi, grid)`; `energy(psi, V, grid)`; devuelven `float` (sección 7).
- Transmisión medida: sección 10. Ancho de barrera para la fórmula: `count_nonzero(V) * dx`.
