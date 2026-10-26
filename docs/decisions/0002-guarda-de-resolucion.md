# 0002. Guarda de resolución en check_resolution, llamada por gaussian_packet
- Date: 2026-10-26
- Status: accepted
- Decision: check_resolution(grid, k0, sigma) vive en quantum1d/grid.py, se exporta desde quantum1d, devuelve None y lanza ValueError si sigma <= 0 o si pi/dx < 4*(abs(k0) + 3/sigma) (la igualdad pasa). gaussian_packet la llama antes de construir psi. Ni evolve ni Grid hacen esta guarda.
- Reason: evolve recibe un psi arbitrario y no conoce k0 ni sigma; Grid no conoce el paquete. El contrato prevé el helper para ese caso. Se usa abs(k0) porque el aliasing depende del módulo del número de onda.
