# Lecciones aprendidas

- 2026-10-26: Con V = 0 o constante, el split-step es exacto para cualquier dt, así que una prueba de convergencia en dt con potencial libre solo mide redondeo; hay que usar un V con V' != 0 donde vive psi. Source: planificación de la Fase 0 (criterio d).
- 2026-10-26: Una barrera muestreada con un borde justo sobre un nodo queda a merced del redondeo (ancho ±dx, ~0.013 en T); hay que elegir center y width para que width/dx sea entero y ningún borde caiga en un nodo. Source: planificación del experimento 01 de túnel.
- 2026-10-26: En este entorno no hay uv ni se puede instalar nada; los checks se escriben como `PYTHONPATH=src python ...` contra el Python local. Source: restricción del entorno en la Fase 0.
- 2026-10-26: Con potenciales discontinuos, el split-step necesita dt*(pi/dx)**2/2 <= pi. Si no se cumple, los modos de k alto entran en resonancia de fase (mod 2π): la norma sigue exacta pero <H> y la radiación explotan, así que la guarda de aliasing en x no alcanza y hay que acotar también la fase por paso. Source: AC18 falló en la revisión 1 (fase 19.7, energía +15 %, P_near_edges 5e-6).
