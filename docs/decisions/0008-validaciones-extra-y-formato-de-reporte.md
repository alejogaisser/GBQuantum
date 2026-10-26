# 0008. Validaciones extra de la API y formato exacto del reporte
- Date: 2026-10-26
- Status: accepted
- Decision: se mantienen las validaciones que no figuran en el contrato textual, aprobadas por el usuario. Lanzan ValueError: Grid con n < 2 o con length no finito o <= 0; barrier con width <= 0; harmonic con omega <= 0; evolve con dt no finito, V complejo o formas inválidas. outputs/report.md escribe los floats con repr(float(v)), y el experimento 01 aborta con ValueError si una corrida supera dt*(pi/dx)**2/2 > pi (MAX_PHASE_PER_STEP).
- Reason: evitan fallas silenciosas sin cambiar firmas. repr se recupera exacto al parsear, así que los AC con umbral de 1e-10 miden el error numérico y no el redondeo de impresión (revisión 2).
