# 0004. Barrera con intervalo cerrado y ancho efectivo
- Date: 2026-10-26
- Status: accepted
- Decision: barrier pone height en los nodos con abs(x - center) <= width/2. El ancho efectivo es count_nonzero(V)*dx y es el que se usa en la fórmula de T. Los experimentos eligen width/dx entero y ningún borde sobre un nodo: el experimento 01 usa center = 0.01, width = 2.0.
- Reason: con un borde sobre un nodo, el redondeo decide si el nodo entra y el ancho varía en dx, lo que cambia T en ~0.013 absoluto, más que la tolerancia de 0.01.
