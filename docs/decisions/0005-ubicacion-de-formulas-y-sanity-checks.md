# 0005. Dónde viven las fórmulas analíticas y los sanity checks del ROL A
- Date: 2026-10-26
- Status: accepted
- Decision: las fórmulas de T (onda plana y paquete) viven solo en experiments/01_tunneling.py, como funciones importables y con guarda `if __name__ == "__main__":`. Euler explícito vive solo en experiments/02_euler_vs_splitstep.py. Los sanity checks del tester del ROL A viven en su scratchpad, fuera del repo, y no se commitean. src/ no contiene fórmulas analíticas de verificación y tests/ es exclusivo del ROL B.
- Reason: el contrato prohíbe fórmulas de verificación en src/ y reserva tests/ para el auditor ciego; el layout del usuario solo prevé dos scripts en experiments/.
