# Lecciones aprendidas

- 2026-10-26: Con V = 0 o constante, el split-step es exacto para cualquier dt, así que una prueba de convergencia en dt con potencial libre solo mide redondeo; hay que usar un V con V' != 0 donde vive psi. Source: planificación de la Fase 0 (criterio d).
- 2026-10-26: Una barrera muestreada con un borde justo sobre un nodo queda a merced del redondeo (ancho ±dx, ~0.013 en T); hay que elegir center y width para que width/dx sea entero y ningún borde caiga en un nodo. Source: planificación del experimento 01 de túnel.
- 2026-10-26: En este entorno no hay uv ni se puede instalar nada; los checks se escriben como `PYTHONPATH=src python ...` contra el Python local. Source: restricción del entorno en la Fase 0.
